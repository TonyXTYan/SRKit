//! WebAssembly bindings for SRKit.
//!
//! A plain C ABI, with no generated glue: numbers go in as `f64` arguments or through a shared input buffer,
//! results come back as the return value or through a shared output buffer, and exact times cross as decimal
//! text. The module imports nothing, so it instantiates synchronously from bytes in a browser, a worker or Node.
//! The hand-written loader `js/srkit.js` wraps these exports in ordinary JavaScript objects.
//!
//! ## Buffers
//!
//! * `srkit_in()`: 64 `f64` slots the caller writes before a call (a segment, mostly).
//! * `srkit_out()`: 192 `f64` slots a call writes its results to.
//! * `srkit_text(slot)`: three 256-byte ASCII buffers for exact times (inputs in slots 0 and 1, results in 2).
//! * `srkit_batch(len)`: a growable `f64` buffer for batch calls.
//!
//! ## A segment in the input buffer (21 slots from an offset)
//!
//! | Slot | rest (0) | inertial (1) | hyperbolic (2) | circular (3) | null (4) | steered (5) | burn (6) |
//! |---|---|---|---|---|---|---|---|
//! | 0 | kind | | | | | | |
//! | 1, 2 | `t0`, `t1` | | | | | | |
//! | 3..6 | `p` | `p0` | `p0` | centre | `p0` | `p0` | `p0` |
//! | 6..9 | | `v` | `n` | `r`, `omega`, - | `n` | `u0` | `n` |
//! | 9.. | | | `alpha`, `phi0` | `e1` (8..11), `e2` (11..14), `a0` (14) | `p1` (9..12) | thrust (9..12), `alpha` (12) | `phi0`, `alpha0`, `w` |
//! | 20 | `tau0` (worldline calls only) | | | | | | |
//!
//! Functions that can fail return a status: 0 ok, 1 no intersection or no answer, 2 invalid input, 3 the result
//! does not fit. Functions that return one number return NaN where there is no answer.
//!
//! There are no `unsafe` blocks. `#[no_mangle]` fixes the export names, which is all this crate's ABI is.

use std::cell::RefCell;

use srkit_core::{
    causal::{self, GridVec},
    frames, kepler,
    kinematics::{self, Event},
    optics,
    orbit::CircularOrbit,
    plan::{self, EntryOptions, Phase, Piece, Pieces, Start},
    reach, rocket,
    segment::{Burn, Null, Segment, Steered},
    time::{Duration, Timestamp},
    trip::{FlipBurn, Trip},
    worldline::{Anchored, WorldlineBuf},
    Error, Vec3,
};

/// Bumped whenever an export's meaning or a buffer layout changes.
pub const ABI_VERSION: u32 = 4;
/// Slots a segment occupies in the input buffer.
pub const SEG_SLOTS: usize = 21;
const IN_SLOTS: usize = 64;
const OUT_SLOTS: usize = 192;
const TEXT_BYTES: usize = 256;

thread_local! {
    static IN: RefCell<[f64; IN_SLOTS]> = const { RefCell::new([0.0; IN_SLOTS]) };
    static OUT: RefCell<[f64; OUT_SLOTS]> = const { RefCell::new([0.0; OUT_SLOTS]) };
    static TEXT: RefCell<[[u8; TEXT_BYTES]; 3]> = const { RefCell::new([[0; TEXT_BYTES]; 3]) };
    static BATCH: RefCell<Vec<f64>> = const { RefCell::new(Vec::new()) };
    static WORLDLINES: RefCell<Arena> = const { RefCell::new(Arena { items: Vec::new(), free: Vec::new(), building: Vec::new() }) };
    static LAST_TRIP: RefCell<Option<([f64; 5], Trip)>> = const { RefCell::new(None) };
}

struct Arena {
    items: Vec<Option<WorldlineBuf>>,
    free: Vec<u32>,
    building: Vec<Anchored>,
}

const OK: u32 = 0;
const NONE: u32 = 1;
const INVALID: u32 = 2;
const OVERFLOW: u32 = 3;

fn status(e: Error) -> u32 {
    match e {
        Error::NoIntersection => NONE,
        Error::Overflow | Error::PrecisionLoss => OVERFLOW,
        _ => INVALID,
    }
}

fn out(values: &[f64]) {
    OUT.with(|o| o.borrow_mut()[..values.len()].copy_from_slice(values));
}

fn v3(a: &[f64], i: usize) -> Vec3 {
    Vec3::new(a[i], a[i + 1], a[i + 2])
}

fn segment_from(a: &[f64]) -> Option<Segment> {
    let (t0, t1) = (a[1], a[2]);
    Some(match a[0] as u32 {
        0 => Segment::rest(t0, t1, v3(a, 3)),
        1 => Segment::inertial(t0, t1, v3(a, 3), v3(a, 6)),
        2 => Segment::hyperbolic(t0, t1, v3(a, 3), v3(a, 6), a[9], a[10]),
        3 => Segment::circular(t0, t1, v3(a, 3), a[6], a[7], v3(a, 8), v3(a, 11), a[14]),
        4 => {
            let mut g = Null::new(t0, t1, v3(a, 3), v3(a, 6));
            g.p1 = v3(a, 9);
            Segment::Null(g)
        }
        5 => Segment::Steered(
            Steered::from_proper_velocity(t0, t1, v3(a, 3), v3(a, 6), v3(a, 9), a[12]).ok()?,
        ),
        6 => Segment::Burn(Burn::new(t0, t1, v3(a, 3), v3(a, 6), a[9], a[10], a[11]).ok()?),
        _ => return None,
    })
}

/// The segment at slot 0 of the input buffer.
fn segment() -> Option<Segment> {
    IN.with(|i| segment_from(&i.borrow()[..SEG_SLOTS]))
}

fn write_segment(seg: &Segment, o: &mut [f64]) {
    o[..SEG_SLOTS].fill(0.0);
    o[1] = seg.t0();
    o[2] = seg.t1();
    let mut put = |i: usize, v: Vec3| {
        o[i] = v.x;
        o[i + 1] = v.y;
        o[i + 2] = v.z;
    };
    match seg {
        Segment::Rest(g) => put(3, g.p),
        Segment::Inertial(g) => {
            put(3, g.p0);
            put(6, g.v);
            o[0] = 1.0;
        }
        Segment::Hyperbolic(g) => {
            put(3, g.p0);
            put(6, g.n);
            o[0] = 2.0;
            o[9] = g.alpha;
            o[10] = g.phi0;
        }
        Segment::Circular(g) => {
            put(3, g.centre);
            put(8, g.e1);
            put(11, g.e2);
            o[0] = 3.0;
            o[6] = g.r;
            o[7] = g.omega;
            o[14] = g.a0;
        }
        Segment::Null(g) => {
            put(3, g.p0);
            put(6, g.n);
            put(9, g.p1);
            o[0] = 4.0;
        }
        Segment::Steered(g) => {
            put(3, g.p0);
            put(6, g.u0);
            put(9, g.n);
            o[0] = 5.0;
            o[12] = g.alpha;
        }
        Segment::Burn(g) => {
            put(3, g.p0);
            put(6, g.n);
            o[0] = 6.0;
            o[9] = g.phi0;
            o[10] = g.sigma * g.w / g.t_burn;
            o[11] = g.w;
        }
    }
}

// ---------------------------------------------------------------------------------------------- buffers

#[no_mangle]
pub extern "C" fn srkit_abi_version() -> u32 {
    ABI_VERSION
}
#[no_mangle]
pub extern "C" fn srkit_in() -> usize {
    IN.with(|v| v.borrow().as_ptr() as usize)
}
#[no_mangle]
pub extern "C" fn srkit_out() -> usize {
    OUT.with(|v| v.borrow().as_ptr() as usize)
}
#[no_mangle]
pub extern "C" fn srkit_text(slot: u32) -> usize {
    TEXT.with(|v| v.borrow()[(slot as usize).min(2)].as_ptr() as usize)
}
/// A buffer of at least `len` numbers for batch calls. The pointer is valid until the next call that grows it.
#[no_mangle]
pub extern "C" fn srkit_batch(len: u32) -> usize {
    BATCH.with(|b| {
        let mut b = b.borrow_mut();
        if b.len() < len as usize {
            b.resize(len as usize, 0.0);
        }
        b.as_ptr() as usize
    })
}
/// Set one input slot (for hosts without typed arrays).
#[no_mangle]
pub extern "C" fn srkit_set(index: u32, value: f64) {
    IN.with(|v| {
        if let Some(slot) = v.borrow_mut().get_mut(index as usize) {
            *slot = value;
        }
    });
}
/// Read one output slot.
#[no_mangle]
pub extern "C" fn srkit_get(index: u32) -> f64 {
    OUT.with(|v| v.borrow().get(index as usize).copied().unwrap_or(f64::NAN))
}

// ---------------------------------------------------------------------------------------------- segments

/// Position (out 0..3) and velocity (out 3..6) of the input segment at coordinate time `t`.
#[no_mangle]
pub extern "C" fn seg_at(t: f64) -> u32 {
    let Some(seg) = segment() else { return INVALID };
    let s = seg.at(t);
    out(&[s.x.x, s.x.y, s.x.z, s.v.x, s.v.y, s.v.z]);
    OK
}
/// Proper time from the segment's start to coordinate time `t`.
#[no_mangle]
pub extern "C" fn seg_dtau(t: f64) -> f64 {
    segment().map_or(f64::NAN, |s| s.dtau(t))
}
/// Coordinate time after proper time `s` from the segment's start.
#[no_mangle]
pub extern "C" fn seg_t_of_tau(s: f64) -> f64 {
    segment().map_or(f64::NAN, |g| g.t_of_tau(s))
}
/// State after proper time `s`: out = t, x (3), v (3), phi, acc, jump.
#[no_mangle]
pub extern "C" fn seg_at_tau(s: f64) -> u32 {
    let Some(seg) = segment() else { return INVALID };
    let e = seg.at_tau(s);
    out(&[
        e.t,
        e.x.x,
        e.x.y,
        e.x.z,
        e.v.x,
        e.v.y,
        e.v.z,
        e.phi,
        e.acc,
        if e.jump { 1.0 } else { 0.0 },
    ]);
    OK
}
/// Largest speed on the segment.
#[no_mangle]
pub extern "C" fn seg_vmax() -> f64 {
    segment().map_or(f64::NAN, |s| s.vmax())
}
/// When light from the event `(te, x, y, z)` reaches the segment, within `[lo, hi]`; NaN if it does not.
#[no_mangle]
pub extern "C" fn seg_recv(te: f64, x: f64, y: f64, z: f64, lo: f64, hi: f64) -> f64 {
    segment()
        .and_then(|s| s.recv(Event::new(te, Vec3::new(x, y, z)), lo, hi))
        .unwrap_or(f64::NAN)
}
/// The time on the segment, within `[lo, hi]`, whose light reaches the event; NaN if there is none.
#[no_mangle]
pub extern "C" fn seg_send(te: f64, x: f64, y: f64, z: f64, lo: f64, hi: f64) -> f64 {
    segment()
        .and_then(|s| s.send(Event::new(te, Vec3::new(x, y, z)), lo, hi))
        .unwrap_or(f64::NAN)
}
/// Margin inside the horizon of a for-ever-accelerating segment (`past != 0`: the past horizon); NaN for kinds
/// with no horizon.
#[no_mangle]
pub extern "C" fn seg_horizon_margin(te: f64, x: f64, y: f64, z: f64, past: u32) -> f64 {
    segment()
        .and_then(|s| s.horizon_margin(Event::new(te, Vec3::new(x, y, z)), past != 0))
        .unwrap_or(f64::NAN)
}
/// Proper time the whole segment takes.
#[no_mangle]
pub extern "C" fn seg_proper_duration() -> f64 {
    segment().map_or(f64::NAN, |s| s.proper_duration())
}

// ---------------------------------------------------------------------------------------------- worldlines

/// Start building a worldline. Follow with `wl_push` for each segment, then `wl_finish`.
#[no_mangle]
pub extern "C" fn wl_begin() {
    WORLDLINES.with(|w| w.borrow_mut().building.clear());
}
/// Append the input segment (its `tau0` in slot 20 is kept as given). Returns the number of segments so far, or 0
/// if the segment is invalid.
#[no_mangle]
pub extern "C" fn wl_push() -> u32 {
    let Some((seg, tau0)) = IN.with(|i| {
        let i = i.borrow();
        segment_from(&i[..SEG_SLOTS]).map(|s| (s, i[20]))
    }) else {
        return 0;
    };
    WORLDLINES.with(|w| {
        let mut w = w.borrow_mut();
        w.building.push(Anchored::new(seg, tau0));
        w.building.len() as u32
    })
}
/// Finish the worldline. With `anchor != 0` every segment's proper-time anchor is recomputed so that proper time
/// is `anchor_tau` at `anchor_t`; otherwise the anchors pushed are used. Returns a handle (1 or more), or 0 if
/// the segments are empty or out of order.
#[no_mangle]
pub extern "C" fn wl_finish(anchor: u32, anchor_t: f64, anchor_tau: f64) -> u32 {
    WORLDLINES.with(|w| {
        let mut w = w.borrow_mut();
        let mut segs = std::mem::take(&mut w.building);
        if anchor != 0 {
            srkit_core::worldline::anchor(&mut segs, anchor_t, anchor_tau);
        }
        let Ok(buf) = WorldlineBuf::from_anchored(segs) else {
            return 0;
        };
        if let Some(slot) = w.free.pop() {
            w.items[slot as usize] = Some(buf);
            slot + 1
        } else {
            w.items.push(Some(buf));
            w.items.len() as u32
        }
    })
}
/// Release a worldline.
#[no_mangle]
pub extern "C" fn wl_free(handle: u32) {
    WORLDLINES.with(|w| {
        let mut w = w.borrow_mut();
        let slot = handle.wrapping_sub(1) as usize;
        if w.items.get(slot).is_some_and(|i| i.is_some()) {
            w.items[slot] = None;
            w.free.push(slot as u32);
        }
    });
}
/// How many worldlines are alive (for leak checks).
#[no_mangle]
pub extern "C" fn wl_live() -> u32 {
    WORLDLINES.with(|w| w.borrow().items.iter().filter(|i| i.is_some()).count() as u32)
}

fn with_wl<R>(handle: u32, none: R, f: impl FnOnce(&WorldlineBuf) -> R) -> R {
    WORLDLINES.with(|w| {
        match w
            .borrow()
            .items
            .get(handle.wrapping_sub(1) as usize)
            .and_then(|i| i.as_ref())
        {
            Some(buf) => f(buf),
            None => none,
        }
    })
}

/// Number of segments.
#[no_mangle]
pub extern "C" fn wl_len(handle: u32) -> u32 {
    with_wl(handle, 0, |w| w.segments().len() as u32)
}
/// Copy segment `index` to the output buffer (the input layout, with `tau0` in slot 20).
#[no_mangle]
pub extern "C" fn wl_segment(handle: u32, index: u32) -> u32 {
    with_wl(handle, INVALID, |w| {
        let Some(a) = w.segments().get(index as usize) else {
            return INVALID;
        };
        OUT.with(|o| {
            let mut o = o.borrow_mut();
            write_segment(&a.seg, &mut o[..]);
            o[20] = a.tau0;
        });
        OK
    })
}
/// Index of the segment in force at coordinate time `t`.
#[no_mangle]
pub extern "C" fn wl_find(handle: u32, t: f64) -> u32 {
    with_wl(handle, 0, |w| w.view().find(t) as u32)
}
/// Position (out 0..3), velocity (3..6) and proper time (6) at coordinate time `t`.
#[no_mangle]
pub extern "C" fn wl_at(handle: u32, t: f64) -> u32 {
    with_wl(handle, INVALID, |w| {
        let v = w.view();
        let s = v.at(t);
        out(&[s.x.x, s.x.y, s.x.z, s.v.x, s.v.y, s.v.z, v.tau(t)]);
        OK
    })
}
/// Proper time at coordinate time `t`.
#[no_mangle]
pub extern "C" fn wl_tau(handle: u32, t: f64) -> f64 {
    with_wl(handle, f64::NAN, |w| w.view().tau(t))
}
/// Coordinate time at proper time `tau`.
#[no_mangle]
pub extern "C" fn wl_t_of_tau(handle: u32, tau: f64) -> f64 {
    with_wl(handle, f64::NAN, |w| w.view().t_of_tau(tau))
}
/// State at proper time `tau`: out = t, x (3), v (3), phi, acc, jump.
#[no_mangle]
pub extern "C" fn wl_at_tau(handle: u32, tau: f64) -> u32 {
    with_wl(handle, INVALID, |w| {
        let e = w.view().at_tau(tau);
        out(&[
            e.t,
            e.x.x,
            e.x.y,
            e.x.z,
            e.v.x,
            e.v.y,
            e.v.z,
            e.phi,
            e.acc,
            if e.jump { 1.0 } else { 0.0 },
        ]);
        OK
    })
}
/// When light from the event reaches the worldline (`+inf` if never).
#[no_mangle]
pub extern "C" fn wl_recv(handle: u32, te: f64, x: f64, y: f64, z: f64) -> f64 {
    with_wl(handle, f64::NAN, |w| {
        w.view().recv(Event::new(te, Vec3::new(x, y, z)))
    })
}
/// The time on the worldline whose light reaches the event (`-inf` if none).
#[no_mangle]
pub extern "C" fn wl_send(handle: u32, te: f64, x: f64, y: f64, z: f64) -> f64 {
    with_wl(handle, f64::NAN, |w| {
        w.view().send(Event::new(te, Vec3::new(x, y, z)))
    })
}
/// Largest speed on the worldline.
#[no_mangle]
pub extern "C" fn wl_vmax(handle: u32) -> f64 {
    with_wl(handle, f64::NAN, |w| w.view().vmax())
}
/// 0 if consecutive segments join within the tolerances, 2 if not.
#[no_mangle]
pub extern "C" fn wl_check(handle: u32, position_tol: f64, velocity_tol: f64) -> u32 {
    with_wl(handle, INVALID, |w| {
        match w.view().check_continuity(position_tol, velocity_tol) {
            Ok(()) => OK,
            Err(_) => INVALID,
        }
    })
}
/// Radar coordinates of the event for this observer: out = tau_send, tau_recv.
#[no_mangle]
pub extern "C" fn wl_radar(handle: u32, te: f64, x: f64, y: f64, z: f64) -> u32 {
    with_wl(handle, INVALID, |w| {
        let r = frames::radar(&w.view(), Event::new(te, Vec3::new(x, y, z)));
        out(&[r.tau_send, r.tau_recv]);
        OK
    })
}
/// Batch evaluation: the batch buffer holds `n` coordinate times; they are replaced by `n` rows of seven numbers
/// (x, y, z, vx, vy, vz, tau). The buffer must hold `7 n` numbers.
#[no_mangle]
pub extern "C" fn wl_at_batch(handle: u32, n: u32) -> u32 {
    with_wl(handle, INVALID, |w| {
        BATCH.with(|b| {
            let mut b = b.borrow_mut();
            let n = n as usize;
            if b.len() < 7 * n {
                return INVALID;
            }
            let v = w.view();
            // back to front, so that a row never overwrites a time not yet read
            for i in (0..n).rev() {
                let t = b[i];
                let s = v.at(t);
                let row = [s.x.x, s.x.y, s.x.z, s.v.x, s.v.y, s.v.z, v.tau(t)];
                b[7 * i..7 * i + 7].copy_from_slice(&row);
            }
            OK
        })
    })
}
/// Batch light cones: the batch buffer holds `n` events as rows of four (t, x, y, z); each row's `t` is replaced
/// by the time on the worldline whose light reaches that event (`send != 0`) or at which the event's light
/// reaches the worldline.
#[no_mangle]
pub extern "C" fn wl_cone_batch(handle: u32, n: u32, send: u32) -> u32 {
    with_wl(handle, INVALID, |w| {
        BATCH.with(|b| {
            let mut b = b.borrow_mut();
            let n = n as usize;
            if b.len() < 4 * n {
                return INVALID;
            }
            let v = w.view();
            for i in 0..n {
                let e = Event::new(
                    b[4 * i],
                    Vec3::new(b[4 * i + 1], b[4 * i + 2], b[4 * i + 3]),
                );
                b[4 * i] = if send != 0 { v.send(e) } else { v.recv(e) };
            }
            OK
        })
    })
}

// ---------------------------------------------------------------------------------------------- rides

fn trip(d: f64, a: f64, cap: f64, stops: f64, rt: f64) -> Result<Trip, Error> {
    if !(stops.is_finite() && (0.0..=10_000.0).contains(&stops) && stops.fract() == 0.0)
        || !(rt == 0.0 || rt == 1.0)
    {
        return Err(Error::InvalidInput);
    }
    let key = [d, a, cap, stops, rt];
    LAST_TRIP.with(|last| {
        let mut last = last.borrow_mut();
        if let Some((k, t)) = last.as_ref() {
            // bitwise comparison: the same inputs, including the same infinity
            if k.iter().zip(&key).all(|(p, q)| p.to_bits() == q.to_bits()) {
                return Ok(*t);
            }
        }
        let t = Trip::new(d, a, cap, stops as u32, rt == 1.0)?;
        *last = Some((key, t));
        Ok(t)
    })
}

/// Plan a ride. out = d, a, n, n_out, round_trip, leg_tau, tau, du, dv, u_end, v_end, t_end, t_out, x_end, lag,
/// gamma_peak, phc, leg.a, leg.phc, leg.ta, leg.tc.
#[no_mangle]
pub extern "C" fn trip_plan(d: f64, a: f64, cap: f64, stops: f64, rt: f64) -> u32 {
    match trip(d, a, cap, stops, rt) {
        Ok(p) => {
            out(&[
                p.d,
                p.a,
                f64::from(p.n),
                f64::from(p.n_out),
                if p.round_trip { 1.0 } else { 0.0 },
                p.leg_tau,
                p.tau,
                p.du,
                p.dv,
                p.u_end,
                p.v_end,
                p.t_end,
                p.t_out,
                p.x_end,
                p.lag,
                p.gamma_peak,
                p.leg.phc,
                p.leg.a,
                p.leg.phc,
                p.leg.ta,
                p.leg.tc,
            ]);
            OK
        }
        Err(e) => status(e),
    }
}
/// State on a ride at proper time `tau` (`clamped != 0`: held within the ride). out = tau, u, v, t, x, phi, dir,
/// gamma, doppler_home, doppler_dest.
#[no_mangle]
pub extern "C" fn trip_state(
    d: f64,
    a: f64,
    cap: f64,
    stops: f64,
    rt: f64,
    tau: f64,
    clamped: u32,
) -> u32 {
    match trip(d, a, cap, stops, rt) {
        Ok(p) => {
            let s = if clamped != 0 {
                p.state_clamped(tau)
            } else {
                p.state(tau)
            };
            out(&[
                s.tau,
                s.u,
                s.v,
                s.t,
                s.x,
                s.phi,
                s.dir,
                s.gamma,
                s.doppler_home,
                s.doppler_dest,
            ]);
            OK
        }
        Err(e) => status(e),
    }
}
/// Proper time at which the ride's `u = t - x` (`which = 0`) or `v = t + x` (`which = 1`) equals `value`.
#[no_mangle]
pub extern "C" fn trip_tau_of(
    d: f64,
    a: f64,
    cap: f64,
    stops: f64,
    rt: f64,
    value: f64,
    which: u32,
) -> f64 {
    match trip(d, a, cap, stops, rt) {
        Ok(p) if which == 0 => p.tau_of_u(value),
        Ok(p) => p.tau_of_v(value),
        Err(_) => f64::NAN,
    }
}
/// Proper time at which the ride meets a light cone of the event `(t, along, across)`.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub extern "C" fn trip_light_cone(
    d: f64,
    a: f64,
    cap: f64,
    stops: f64,
    rt: f64,
    t: f64,
    along: f64,
    across: f64,
    send: u32,
) -> f64 {
    trip(d, a, cap, stops, rt)
        .and_then(|p| p.light_cone_tau(t, along, across, send != 0))
        .unwrap_or(f64::NAN)
}
/// Rest-to-rest flip-and-burn. out = t, tau, phi_peak, gamma_peak, lag.
#[no_mangle]
pub extern "C" fn flip_burn(d: f64, a: f64) -> u32 {
    match FlipBurn::new(d, a) {
        Ok(f) => {
            out(&[f.t, f.tau, f.phi_peak, f.gamma_peak, f.lag]);
            OK
        }
        Err(e) => status(e),
    }
}
/// Where a ship can stop: out = along, across. Status 1 if the budget cannot stop it.
#[no_mangle]
pub extern "C" fn reach_point(phi: f64, a: f64, budget: f64, heading: f64) -> u32 {
    match reach::reach(phi, a, budget, heading) {
        Ok(Some(p)) => {
            out(&p);
            OK
        }
        Ok(None) => NONE,
        Err(e) => status(e),
    }
}

// ---------------------------------------------------------------------------------------------- itineraries

const HEADER: usize = 24;
const PIECE: usize = 3 + SEG_SLOTS;

fn phase_code(p: Phase) -> f64 {
    match p {
        Phase::Rest => 0.0,
        Phase::Boost => 1.0,
        Phase::Coast => 2.0,
        Phase::Brake => 3.0,
        Phase::Orbit => 4.0,
        Phase::Jump => 5.0,
    }
}

fn write_piece(p: &Piece, o: &mut [f64]) {
    o[0] = phase_code(p.phase);
    o[1] = p.tau0;
    o[2] = p.dtau;
    write_segment(&p.seg, &mut o[3..]);
}

/// Header: 0 piece count, 1 end.t, 2 end.tau, 3..6 end.x, 6..9 end.v, 9 end.phi, 10 end.acc, then the caller's
/// extras from slot 11. Pieces follow from slot 24, 24 slots each: phase, tau0, dtau, segment.
fn write_plan(pieces: &Pieces, end: &plan::End, extras: &[f64]) {
    OUT.with(|o| {
        let mut o = o.borrow_mut();
        o[..HEADER].fill(0.0);
        o[0] = pieces.len() as f64;
        o[1] = end.t;
        o[2] = end.tau;
        o[3] = end.x.x;
        o[4] = end.x.y;
        o[5] = end.x.z;
        o[6] = end.v.x;
        o[7] = end.v.y;
        o[8] = end.v.z;
        o[9] = end.phi;
        o[10] = end.acc;
        o[11..11 + extras.len()].copy_from_slice(extras);
        for (k, p) in pieces.iter().enumerate() {
            write_piece(p, &mut o[HEADER + k * PIECE..]);
        }
    });
}

/// A straight run (see `srkit_core::plan::line`). In: slots 0..3 start x, 3..6 target, 6..9 default direction.
/// Extras: 11..14 n, 14 dphi, 15 length, 16 phi_peak, 17 coast_length. Status 1 if the run is too short.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub extern "C" fn plan_line(
    t: f64,
    tau: f64,
    a: f64,
    phi_cap: f64,
    phi0: f64,
    phi1: f64,
    a_brake: f64,
    has_default: u32,
) -> u32 {
    let (x, x1, nd) = IN.with(|i| {
        let i = i.borrow();
        (v3(&i[..], 0), v3(&i[..], 3), v3(&i[..], 6))
    });
    let start = Start { t, tau, x };
    match plan::line(
        start,
        x1,
        a,
        phi_cap,
        phi0,
        phi1,
        a_brake,
        (has_default != 0).then_some(nd),
    ) {
        Some(l) => {
            write_plan(
                &l.pieces,
                &l.end,
                &[
                    l.n.x,
                    l.n.y,
                    l.n.z,
                    l.dphi,
                    l.length,
                    l.phi_peak,
                    l.coast_length,
                ],
            );
            OK
        }
        None => NONE,
    }
}

/// A plan of at most one piece, in the plan layout, with the rapidity spent in slot 14.
fn single(piece: Option<Piece>, end: plan::End, dphi: f64) {
    OUT.with(|o| {
        let mut o = o.borrow_mut();
        o[..HEADER].fill(0.0);
        o[0] = if piece.is_some() { 1.0 } else { 0.0 };
        o[1] = end.t;
        o[2] = end.tau;
        o[3] = end.x.x;
        o[4] = end.x.y;
        o[5] = end.x.z;
        o[6] = end.v.x;
        o[7] = end.v.y;
        o[8] = end.v.z;
        o[9] = end.phi;
        o[10] = end.acc;
        o[14] = dphi;
        if let Some(p) = piece {
            write_piece(&p, &mut o[HEADER..]);
        }
    });
}

/// Brake to rest from rapidity `phi0` along the direction in slots 3..6, starting at the position in 0..3.
#[no_mangle]
pub extern "C" fn plan_brake(t: f64, tau: f64, phi0: f64, a: f64) -> u32 {
    let (x, n) = IN.with(|i| {
        let i = i.borrow();
        (v3(&i[..], 0), v3(&i[..], 3))
    });
    let (piece, end, dphi) = plan::brake(Start { t, tau, x }, n, phi0, a);
    single(piece, end, dphi);
    OK
}
/// Leave a circle: brake along the velocity in slots 3..6 from the position in 0..3.
#[no_mangle]
pub extern "C" fn plan_orbit_exit(t: f64, tau: f64, a: f64) -> u32 {
    let (x, v) = IN.with(|i| {
        let i = i.borrow();
        (v3(&i[..], 0), v3(&i[..], 3))
    });
    let (piece, end, dphi) = plan::orbit_exit(Start { t, tau, x }, v, a);
    single(piece, end, dphi);
    OK
}
/// One piece of the given kind from a start: 0 rest (`p1 = dtau`), 1 hyperbolic (direction in 3..6;
/// `p1 = alpha`, `p2 = phi0`, `p3 = phi1`, `p4 = phase code`), 2 coast (velocity in 3..6; `p1 = dtau`), 3 jump
/// (target in 3..6), 4 circle (centre in 3..6, e1 in 6..9, e2 in 9..12; `p1 = r`, `p2 = omega`, `p3 = a0`,
/// `p4 = dtau`). The piece is written at slot 24 as in a plan.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub extern "C" fn plan_piece(
    kind: u32,
    t: f64,
    tau: f64,
    p1: f64,
    p2: f64,
    p3: f64,
    p4: f64,
) -> u32 {
    let (x, b, c, d) = IN.with(|i| {
        let i = i.borrow();
        (v3(&i[..], 0), v3(&i[..], 3), v3(&i[..], 6), v3(&i[..], 9))
    });
    let start = Start { t, tau, x };
    let piece = match kind {
        0 => plan::rest(start, p1),
        1 => plan::hyp_piece(
            start,
            b,
            p1,
            p2,
            p3,
            match p4 as u32 {
                1 => Phase::Boost,
                3 => Phase::Brake,
                _ => Phase::Coast,
            },
        ),
        2 => plan::lin_piece(start, b, p1),
        3 => plan::null_piece(start, b),
        4 => {
            let orb = plan::OrbitElements {
                centre: b,
                r: p1,
                omega: p2,
                e1: c,
                e2: d,
                a0: p3,
                orbit: CircularOrbit::from_radius_speed(p1, (p1 * p2).abs()),
            };
            plan::circ_piece(start, &orb, p4)
        }
        _ => return INVALID,
    };
    single(Some(piece), piece.end(), 0.0);
    OK
}
/// The plane through a centre (slots 0..3) and a point (3..6), with the default `e1` in 6..9 if `has_default`.
/// out = e1, e2, normal.
#[no_mangle]
pub extern "C" fn plan_plane_through(has_default: u32) -> u32 {
    let (c, x, e) = IN.with(|i| {
        let i = i.borrow();
        (v3(&i[..], 0), v3(&i[..], 3), v3(&i[..], 6))
    });
    let p = plan::plane_through(c, x, (has_default != 0).then_some(e));
    out(&[
        p.e1.x, p.e1.y, p.e1.z, p.e2.x, p.e2.y, p.e2.z, p.normal.x, p.normal.y, p.normal.z,
    ]);
    OK
}
/// Enter a circular orbit (see `srkit_core::plan::orbit_entry`). In: 0..3 start x, 3..6 centre, 6..9 the default
/// `e1` (used if `flags & 1`). `flags & 2`: `a0` is given; `flags & 4`: `a_tangent` is given.
/// Extras: 11 dphi, 12..15 e1, 15..18 e2, 18 a0, 19 omega, 20 r, 21 beta, 22 entered from inside. Status 1 if a
/// run is too short.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub extern "C" fn plan_orbit_entry(
    t: f64,
    tau: f64,
    r: f64,
    beta: f64,
    a: f64,
    phi_cap: f64,
    a0: f64,
    a_tangent: f64,
    flags: u32,
) -> u32 {
    let (x, c, e1) = IN.with(|i| {
        let i = i.borrow();
        (v3(&i[..], 0), v3(&i[..], 3), v3(&i[..], 6))
    });
    let opts = EntryOptions {
        e1: (flags & 1 != 0).then_some(e1),
        a0: (flags & 2 != 0).then_some(a0),
        a_tangent: (flags & 4 != 0).then_some(a_tangent),
    };
    match plan::orbit_entry(Start { t, tau, x }, c, r, beta, a, phi_cap, opts) {
        Some(e) => {
            let o = &e.orb;
            write_plan(
                &e.pieces,
                &e.end,
                &[
                    e.dphi,
                    o.e1.x,
                    o.e1.y,
                    o.e1.z,
                    o.e2.x,
                    o.e2.y,
                    o.e2.z,
                    o.a0,
                    o.omega,
                    o.r,
                    o.orbit.beta,
                    if e.from_inside { 1.0 } else { 0.0 },
                ],
            );
            OK
        }
        None => NONE,
    }
}

// ---------------------------------------------------------------------------------------------- orbits, optics

fn orbit_out(o: &CircularOrbit) {
    out(&[
        o.r,
        o.beta,
        o.gamma,
        o.phi,
        o.omega,
        o.period,
        o.lap_tau,
        o.acc,
        o.sagnac(),
        o.thomas_precession_per_lap(),
    ]);
}
/// A circle from its radius and speed. out = r, beta, gamma, phi, omega, period, lap_tau, acc, sagnac, thomas.
#[no_mangle]
pub extern "C" fn orbit_from_speed(r: f64, beta: f64) -> u32 {
    orbit_out(&CircularOrbit::from_radius_speed(r, beta));
    OK
}
/// A circle from its radius and the proper acceleration held aboard (same output).
#[no_mangle]
pub extern "C" fn orbit_from_acceleration(r: f64, acc: f64) -> u32 {
    if !(r.is_finite() && r > 0.0 && acc.is_finite() && acc > 0.0) {
        return INVALID;
    }
    orbit_out(&CircularOrbit::from_radius_acceleration(r, acc));
    OK
}
/// State on a circle (radius, held acceleration) about the origin in the xy plane after proper time `tau`,
/// starting at angle `a0`. out = t, x (3), v (3), heading (3), angle.
#[no_mangle]
pub extern "C" fn orbit_state(r: f64, acc: f64, a0: f64, tau: f64) -> u32 {
    let s = CircularOrbit::from_radius_acceleration(r, acc).state(a0, tau);
    out(&[
        s.t,
        s.x.x,
        s.x.y,
        s.x.z,
        s.v.x,
        s.v.y,
        s.v.z,
        s.heading.x,
        s.heading.y,
        s.heading.z,
        s.angle,
    ]);
    OK
}
/// Proper time at which that orbiting ship meets a light cone of the event `(t, x, y, z)`.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub extern "C" fn orbit_light_cone(
    r: f64,
    acc: f64,
    a0: f64,
    t: f64,
    x: f64,
    y: f64,
    z: f64,
    send: u32,
) -> f64 {
    CircularOrbit::from_radius_acceleration(r, acc).light_cone_tau(
        a0,
        Event::new(t, Vec3::new(x, y, z)),
        send != 0,
    )
}
/// Aberration and Doppler in a ship-aligned basis. out = doppler, theta, theta_seen, forward, right, up.
#[no_mangle]
pub extern "C" fn observe(phi: f64, a: f64, br: f64, bu: f64, aberration: u32) -> u32 {
    match optics::observe(phi, a, br, bu, aberration != 0) {
        Ok(s) => {
            out(&[s.doppler, s.theta, s.theta_seen, s.forward, s.right, s.up]);
            OK
        }
        Err(e) => status(e),
    }
}
/// Rate at which an observer sees a source's clock run. In: 0..3 observer velocity, 3..6 unit vector to the
/// source, 6..9 source velocity.
#[no_mangle]
pub extern "C" fn playback() -> f64 {
    IN.with(|i| {
        let i = i.borrow();
        optics::playback(v3(&i[..], 0), v3(&i[..], 3), v3(&i[..], 6))
    })
}
/// An observer's axes from its velocity (slots 0..3). out = forward, up, right.
#[no_mangle]
pub extern "C" fn basis() -> u32 {
    let v = IN.with(|i| v3(&i.borrow()[..], 0));
    let (f, u, r) = optics::basis(v);
    out(&[f.x, f.y, f.z, u.x, u.y, u.z, r.x, r.y, r.z]);
    OK
}
/// Boost `(dt, dx)` by a rapidity vector. out = x', y', z', t'.
#[no_mangle]
pub extern "C" fn boost(wx: f64, wy: f64, wz: f64, dt: f64, dx: f64, dy: f64, dz: f64) -> u32 {
    let (t, x) = kinematics::boost_rapidity(Vec3::new(wx, wy, wz), dt, Vec3::new(dx, dy, dz));
    out(&[x.x, x.y, x.z, t]);
    OK
}
/// Relativistic composition of velocities `v` then `u`. out = the resulting velocity. Status 2 if either is at or
/// above `c`.
#[no_mangle]
pub extern "C" fn compose_velocities(vx: f64, vy: f64, vz: f64, ux: f64, uy: f64, uz: f64) -> u32 {
    match kinematics::compose_velocities(Vec3::new(vx, vy, vz), Vec3::new(ux, uy, uz)) {
        Ok(w) => {
            out(&[w.x, w.y, w.z]);
            OK
        }
        Err(e) => status(e),
    }
}

// ---------------------------------------------------------------------------------------------- scalars

#[no_mangle]
pub extern "C" fn rocket_mass_ratio(dphi: f64, w: f64) -> f64 {
    rocket::mass_ratio(dphi, w)
}
#[no_mangle]
pub extern "C" fn rocket_rapidity(mass_ratio: f64, w: f64) -> f64 {
    rocket::rapidity_of_mass_ratio(mass_ratio, w)
}
/// Propellant energy per unit of dry rest energy: `(MR - 1) / efficiency`.
#[no_mangle]
pub extern "C" fn rocket_energy(dphi: f64, w: f64, efficiency: f64) -> f64 {
    rocket::energy_per_dry_rest_energy(dphi, w, efficiency)
}
/// Rapidity of a perfect sail given beam energy over its rest energy.
#[no_mangle]
pub extern "C" fn sail_rapidity(energy_per_rest_energy: f64) -> f64 {
    rocket::sail_rapidity(energy_per_rest_energy)
}
#[no_mangle]
pub extern "C" fn kepler_eccentric_anomaly(mean_anomaly: f64, e: f64) -> f64 {
    kepler::eccentric_anomaly(mean_anomaly, e)
}
#[no_mangle]
pub extern "C" fn rayleigh_resolution(aperture: f64, wavelength: f64, distance: f64) -> f64 {
    optics::rayleigh_resolution(aperture, wavelength, distance)
}
#[no_mangle]
pub extern "C" fn friis_received_power(
    power: f64,
    d_t: f64,
    d_r: f64,
    wavelength: f64,
    distance: f64,
) -> f64 {
    optics::friis_received_power(power, d_t, d_r, wavelength, distance)
}
#[no_mangle]
pub extern "C" fn photon_bit_rate(power_w: f64, wavelength_m: f64, bits_per_photon: f64) -> f64 {
    optics::photon_bit_rate(power_w, wavelength_m, bits_per_photon)
}
#[no_mangle]
pub extern "C" fn black_body_luminosity(radius_m: f64, temperature_k: f64) -> f64 {
    optics::black_body_luminosity(radius_m, temperature_k)
}
/// A constant from `srkit_core::units` by index: 0 c, 1 Julian year, 2 light-year, 3 au, 4 parsec, 5 g0, 6 h,
/// 7 hbar, 8 k_B, 9 G, 10 sigma, 11 L_sun, 12 R_sun, 13 GM_sun, 14 g0 per year, 15 au in ly, 16 parsec in ly.
#[no_mangle]
pub extern "C" fn unit_constant(index: u32) -> f64 {
    use srkit_core::units::*;
    [
        C,
        JULIAN_YEAR,
        LIGHT_YEAR,
        AU,
        PARSEC,
        G0,
        PLANCK,
        HBAR,
        BOLTZMANN,
        G_NEWTON,
        STEFAN_BOLTZMANN,
        L_SUN,
        R_SUN,
        GM_SUN,
        G0_PER_YEAR,
        AU_LY,
        PARSEC_LY,
    ]
    .get(index as usize)
    .copied()
    .unwrap_or(f64::NAN)
}

// ---------------------------------------------------------------------------------------------- exact time

fn text(slot: usize, len: u32) -> Option<String> {
    TEXT.with(|t| {
        let t = t.borrow();
        let bytes = t[slot].get(..len as usize)?;
        core::str::from_utf8(bytes).ok().map(str::to_owned)
    })
}

/// Write `s` to the result slot and return its length, or -3 if it does not fit.
fn put_text(s: &str) -> i32 {
    if s.len() > TEXT_BYTES {
        return -(OVERFLOW as i32);
    }
    TEXT.with(|t| t.borrow_mut()[2][..s.len()].copy_from_slice(s.as_bytes()));
    s.len() as i32
}

fn fail(e: Error) -> i32 {
    -(status(e) as i32)
}

fn grid(dx: f64, dy: f64, dz: f64) -> Result<GridVec, Error> {
    GridVec::from_metres(dx, dy, dz)
}

macro_rules! exact_time {
    ($n:literal, $add:ident, $diff:ident, $cmp:ident, $modulo:ident, $seconds:ident, $ticks:ident, $arrival:ident, $arrived:ident, $classify:ident) => {
        /// Timestamp (text slot 0) plus signed duration (slot 1), as decimal milliseconds in slot 2. Returns the
        /// result's length, or minus a status.
        #[no_mangle]
        pub extern "C" fn $add(len_a: u32, len_b: u32) -> i32 {
            let run = || -> Result<String, Error> {
                let a = Timestamp::<$n>::from_decimal(&text(0, len_a).ok_or(Error::InvalidInput)?)?;
                let b = Duration::<$n>::from_decimal(&text(1, len_b).ok_or(Error::InvalidInput)?)?;
                Ok(a.checked_add(b)?.to_string())
            };
            run().map_or_else(fail, |s| put_text(&s))
        }
        /// Timestamp (slot 0) minus timestamp (slot 1), as a signed duration in slot 2.
        #[no_mangle]
        pub extern "C" fn $diff(len_a: u32, len_b: u32) -> i32 {
            let run = || -> Result<String, Error> {
                let a = Timestamp::<$n>::from_decimal(&text(0, len_a).ok_or(Error::InvalidInput)?)?;
                let b = Timestamp::<$n>::from_decimal(&text(1, len_b).ok_or(Error::InvalidInput)?)?;
                Ok(a.signed_duration_since(b).to_string())
            };
            run().map_or_else(fail, |s| put_text(&s))
        }
        /// Order of two timestamps: -1, 0 or 1; -2 if either is invalid.
        #[no_mangle]
        pub extern "C" fn $cmp(len_a: u32, len_b: u32) -> i32 {
            let run = || -> Result<i32, Error> {
                let a = Timestamp::<$n>::from_decimal(&text(0, len_a).ok_or(Error::InvalidInput)?)?;
                let b = Timestamp::<$n>::from_decimal(&text(1, len_b).ok_or(Error::InvalidInput)?)?;
                Ok(a.cmp(&b) as i32)
            };
            run().unwrap_or(-2)
        }
        /// The timestamp in slot 0 modulo a whole number of milliseconds: an exact phase. NaN if invalid.
        #[no_mangle]
        pub extern "C" fn $modulo(len_a: u32, period_ms: f64) -> f64 {
            let run = || -> Result<f64, Error> {
                if !((1.0..=9_007_199_254_740_992.0).contains(&period_ms)
                    && period_ms.fract() == 0.0)
                {
                    return Err(Error::InvalidInput);
                }
                let a = Timestamp::<$n>::from_decimal(&text(0, len_a).ok_or(Error::InvalidInput)?)?;
                Ok(a.modulo_millis(period_ms as u64)? as f64)
            };
            run().unwrap_or(f64::NAN)
        }
        /// Timestamp (slot 0) minus timestamp (slot 1) in seconds, subtracted exactly before the conversion. NaN
        /// if the difference is beyond 2^53 ms.
        #[no_mangle]
        pub extern "C" fn $seconds(len_a: u32, len_b: u32) -> f64 {
            let run = || -> Result<f64, Error> {
                let a = Timestamp::<$n>::from_decimal(&text(0, len_a).ok_or(Error::InvalidInput)?)?;
                let b = Timestamp::<$n>::from_decimal(&text(1, len_b).ok_or(Error::InvalidInput)?)?;
                a.local_seconds_since(b)
            };
            run().unwrap_or(f64::NAN)
        }
        /// Light-travel time over a displacement in metres, rounded up to whole milliseconds, in slot 2.
        #[no_mangle]
        pub extern "C" fn $ticks(dx: f64, dy: f64, dz: f64) -> i32 {
            grid(dx, dy, dz)
                .and_then(causal::light_ticks::<$n>)
                .map_or_else(fail, |d| put_text(&d.to_string()))
        }
        /// The first tick at which a signal emitted at the timestamp in slot 0 has crossed the displacement.
        #[no_mangle]
        pub extern "C" fn $arrival(len_a: u32, dx: f64, dy: f64, dz: f64) -> i32 {
            let run = || -> Result<String, Error> {
                let a = Timestamp::<$n>::from_decimal(&text(0, len_a).ok_or(Error::InvalidInput)?)?;
                Ok(causal::first_arrival(a, grid(dx, dy, dz)?)?.to_string())
            };
            run().map_or_else(fail, |s| put_text(&s))
        }
        /// Whether a signal emitted at the timestamp in slot 0 has crossed the displacement by the one in slot 1:
        /// 1 yes, 0 no, negative if invalid. Exact.
        #[no_mangle]
        pub extern "C" fn $arrived(len_a: u32, len_b: u32, dx: f64, dy: f64, dz: f64) -> i32 {
            let run = || -> Result<i32, Error> {
                let a = Timestamp::<$n>::from_decimal(&text(0, len_a).ok_or(Error::InvalidInput)?)?;
                let b = Timestamp::<$n>::from_decimal(&text(1, len_b).ok_or(Error::InvalidInput)?)?;
                Ok(causal::arrived(a, b, grid(dx, dy, dz)?) as i32)
            };
            run().unwrap_or_else(fail)
        }
        /// Classify the separation of two events a signed duration (slot 0) and a displacement apart: 1 timelike,
        /// 0 null, -1 spacelike; -2 if invalid. Exact.
        #[no_mangle]
        pub extern "C" fn $classify(len_a: u32, dx: f64, dy: f64, dz: f64) -> i32 {
            let run = || -> Result<i32, Error> {
                let a = Duration::<$n>::from_decimal(&text(0, len_a).ok_or(Error::InvalidInput)?)?;
                Ok(match causal::classify(a, grid(dx, dy, dz)?) {
                    causal::Causal::Timelike => 1,
                    causal::Causal::Null => 0,
                    causal::Causal::Spacelike => -1,
                })
            };
            run().unwrap_or(-2)
        }
    };
}

exact_time!(
    6,
    time384_add,
    time384_diff,
    time384_cmp,
    time384_mod,
    time384_seconds,
    time384_light_ticks,
    time384_first_arrival,
    time384_arrived,
    time384_classify
);
exact_time!(
    8,
    time512_add,
    time512_diff,
    time512_cmp,
    time512_mod,
    time512_seconds,
    time512_light_ticks,
    time512_first_arrival,
    time512_arrived,
    time512_classify
);

// ---------------------------------------------------------------------------------------------- cosmology

/// One function of `srkit-flrw`, for the flat ΛCDM universe with Hubble constant `h0_per_year` (per Julian
/// year) and matter density `omega_m`. Times are in years since the Big Bang, lengths in light-years (comoving
/// lengths are proper lengths today). NaN where there is no answer (a signal that never arrives) or the input
/// is invalid.
///
/// | Code | Function | Arguments |
/// |---|---|---|
/// | 0 | age now | |
/// | 1 | scale factor | t |
/// | 2 | time at a scale factor | a |
/// | 3 | Hubble rate | t |
/// | 4 | conformal time | t |
/// | 5 | conformal time at infinity | |
/// | 6, 7 | event horizon: comoving, proper | t |
/// | 8, 9 | particle horizon: comoving, proper | t |
/// | 10 | reach of light sent at t | t |
/// | 11 | last time to signal a target at comoving distance | chi |
/// | 12 | last time to signal and still get a reply | chi |
/// | 13 | arrival time of light | t emitted, chi |
/// | 14 | emission time of light | t received, chi |
/// | 15 | comoving distance light covers | t from, t to |
/// | 16 | lookback time | z |
/// | 17 | time at a redshift | z |
/// | 18 | redshift between two times | t emitted, t received |
/// | 19 | proper velocity (gamma beta) of a coasting body after Hubble drag | u, t from, t to |
/// | 20 | Lorentz factor at which a rocket's thrust balances Hubble drag | alpha |
/// | 21 | de Sitter Hubble rate | |
/// | 22 | de Sitter horizon | |
/// | 23 | de Sitter temperature, K | |
/// | 24 | background temperature, K | t, temperature today |
/// | 25 | time when the background cools to the de Sitter temperature | temperature today |
/// | 26 | mass within reach, kg | t |
/// | 27 | share of the observable universe still reachable | t |
#[no_mangle]
pub extern "C" fn flrw(code: u32, h0_per_year: f64, omega_m: f64, a: f64, b: f64, c: f64) -> f64 {
    let Some(u) = srkit_flrw::Cosmology::per_year(h0_per_year, omega_m) else {
        return f64::NAN;
    };
    let some = |v: Option<f64>| v.unwrap_or(f64::NAN);
    match code {
        0 => u.age(),
        1 => u.scale_factor(a),
        2 => u.time_at_scale_factor(a),
        3 => u.hubble(a),
        4 => u.conformal_time(a),
        5 => u.conformal_time_limit(),
        6 => u.event_horizon(a),
        7 => u.event_horizon_proper(a),
        8 => u.particle_horizon(a),
        9 => u.particle_horizon_proper(a),
        10 => u.reach(a),
        11 => some(u.last_contact_time(a)),
        12 => some(u.last_reply_time(a)),
        13 => some(u.light_arrival_time(a, b)),
        14 => some(u.light_emission_time(a, b)),
        15 => u.comoving_distance(a, b),
        16 => u.lookback_time(a),
        17 => u.time_at_redshift(a),
        18 => u.redshift_between(a, b),
        19 => u.coasting_proper_velocity(a, b, c),
        20 => u.drag_limited_gamma(a),
        21 => u.h_lambda(),
        22 => u.de_sitter_horizon(),
        23 => u.de_sitter_temperature(),
        24 => u.cmb_temperature(a, b),
        25 => u.cmb_crossover_time(a),
        26 => u.reach_mass(a),
        27 => u.reachable_share_of_observable(a),
        _ => f64::NAN,
    }
}
/// The Planck 2018 Hubble constant per Julian year (`which = 0`) and matter density (`which = 1`).
#[no_mangle]
pub extern "C" fn flrw_planck_2018(which: u32) -> f64 {
    let u = srkit_flrw::Cosmology::PLANCK_2018;
    if which == 0 {
        u.h0()
    } else {
        u.omega_m()
    }
}

// ---------------------------------------------------------------------------------------------- gravity

/// One function of `srkit-gr`. SI throughout (kg, m, s, K, W); a spin is the dimensionless `a*`, and where a
/// function takes only a spin, a negative one means a retrograde orbit. NaN outside a function's domain.
///
/// | Code | Function | Arguments |
/// |---|---|---|
/// | 0 | Kerr ISCO radius, in units of GM/c² | spin |
/// | 1 | dt/dtau on the Kerr ISCO | spin |
/// | 2 | accretion efficiency, 1 - E at the ISCO | spin |
/// | 3 | Hawking temperature | mass |
/// | 4, 14 | Hawking lifetime: photons only (black body over the horizon area); Page's photons and gravitons | mass |
/// | 5 | Hawking power, photons only | mass |
/// | 6 | Schwarzschild radius | mass |
/// | 7 | dtau/dt of a static observer | mass, r |
/// | 8 | dt/dtau on a circular orbit (Schwarzschild) | mass, r |
/// | 9 | Shapiro delay, one way | mass, r emitter, r receiver, impact parameter |
/// | 10 | deflection of light | mass, impact parameter |
/// | 11 | gravitational redshift z between static observers | mass, r emitter, r observer |
/// | 12 | Kerr outer horizon radius, in units of GM/c² | spin |
/// | 13 | dtau/dt of a zero-angular-momentum observer (Kerr) | mass, spin, r, polar angle |
/// | 15 | Bekenstein-Hawking entropy, J/K | mass |
/// | 16 | mass whose photon-only lifetime is a given time | lifetime |
/// | 17 | share of a Kerr hole's mass-energy that can be extracted | spin |
#[no_mangle]
pub extern "C" fn gr(code: u32, a: f64, b: f64, c: f64, d: f64) -> f64 {
    use srkit_gr::{hawking, kerr, schwarzschild, weak_field};
    let orbit = if a < 0.0 {
        kerr::Orbit::Retrograde
    } else {
        kerr::Orbit::Prograde
    };
    // ratios to GM/c² do not depend on the mass; one solar mass keeps the arithmetic at ordinary sizes
    let m = srkit_core::units::M_SUN;
    match code {
        0 => kerr::isco_radius(m, a.abs(), orbit) / srkit_gr::gravitational_radius(m),
        1 => kerr::isco_dt_dtau(a.abs(), orbit),
        2 => kerr::accretion_efficiency(a.abs(), orbit),
        3 => hawking::temperature(a),
        4 => hawking::lifetime(a, hawking::EmissionModel::PhotonBlackbody),
        5 => hawking::power(a, hawking::EmissionModel::PhotonBlackbody),
        6 => schwarzschild::schwarzschild_radius(a),
        7 => schwarzschild::static_dtau_dt(a, b),
        8 => schwarzschild::circular_orbit_dt_dtau(a, b),
        9 => weak_field::shapiro_delay(a, b, c, d),
        10 => weak_field::deflection_angle(a, b),
        11 => schwarzschild::gravitational_redshift(a, b, c),
        12 => kerr::outer_horizon_radius(m, a) / srkit_gr::gravitational_radius(m),
        13 => kerr::zamo_dtau_dt(a, b, c, d),
        14 => hawking::lifetime(a, hawking::EmissionModel::PagePhotonsGravitons),
        15 => hawking::entropy(a),
        16 => hawking::mass_for_lifetime(a, hawking::EmissionModel::PhotonBlackbody),
        17 => kerr::extractable_energy_fraction(a),
        _ => f64::NAN,
    }
}

// ---------------------------------------------------------------------------------------------- benchmark

/// The clock-width benchmark (questions R23.8 and R71.4), run inside the module so that the cost measured is
/// the arithmetic and not the crossing: `n` rounds of advancing a time by a tick, taking a difference,
/// comparing, reducing to a yearly phase and making an exact arrival check, on a 384-, 512- or 1,024-bit clock
/// (`width`) with `ticks_per_second` ticks to the second on a grid of `grid_per_metre` steps to the metre.
/// Returns a checksum, or NaN for a width or resolution it does not have.
#[no_mangle]
pub extern "C" fn srkit_bench_clock(
    width: u32,
    ticks_per_second: f64,
    grid_per_metre: f64,
    n: u32,
) -> f64 {
    use srkit_core::causal::Resolution;
    fn run<const N: usize>(res: Resolution, n: u32) -> f64 {
        // 3e107 ticks: 1e97 years at a millisecond, 1e91 at a nanosecond; within every width
        let Ok(epoch) = Timestamp::<N>::from_decimal("300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000") else {
            return f64::NAN;
        };
        let tick = Duration::<N>::from_unsigned_ticks(1);
        let g = res.grid_per_metre() as i128;
        let dx = GridVec::new(40_175_000_000_000_000 * g, -3_000_000_000 * g, 77 * g);
        let year = 31_557_600 * res.ticks_per_second();
        let Ok(arrival) = causal::first_arrival_at(res, epoch, dx) else {
            return f64::NAN;
        };
        let (mut t, mut sum) = (epoch, 0u64);
        for _ in 0..n {
            let Ok(next) = t.checked_add(tick) else {
                return f64::NAN;
            };
            t = next;
            let since = t.signed_duration_since(epoch).magnitude_millis().digits()[0];
            let phase = t.modulo_millis(year).unwrap_or(0);
            let hit = causal::arrived_at(res, epoch, if t > arrival { t } else { arrival }, dx);
            sum = sum.wrapping_add(since ^ phase).wrapping_add(u64::from(hit));
        }
        (sum % (1 << 52)) as f64
    }
    let Ok(res) = Resolution::new(ticks_per_second as u64, grid_per_metre as u64) else {
        return f64::NAN;
    };
    match width {
        384 => run::<6>(res, n),
        512 => run::<8>(res, n),
        1024 => run::<16>(res, n),
        _ => f64::NAN,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(values: &[f64]) {
        for (i, v) in values.iter().enumerate() {
            srkit_set(i as u32, *v);
        }
    }

    #[test]
    fn a_segment_round_trips_through_the_buffers() {
        // hyperbolic from rest at 1 per unit time along +y
        set(&[
            2.0,
            0.0,
            f64::INFINITY,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            1.0,
            0.0,
        ]);
        assert_eq!(seg_at(1.0), OK);
        assert!((srkit_get(1) - (2.0_f64.sqrt() - 1.0)).abs() < 1e-15);
        assert!((seg_dtau(1.0) - 1.0_f64.asinh()).abs() < 1e-15);
        assert!((seg_t_of_tau(1.0_f64.asinh()) - 1.0).abs() < 1e-15);
        // an invalid kind is refused
        srkit_set(0, 99.0);
        assert_eq!(seg_at(1.0), INVALID);
        assert!(seg_dtau(1.0).is_nan());
    }

    #[test]
    fn worldline_handles_build_query_and_free() {
        wl_begin();
        set(&[0.0, 0.0, 10.0, 0.0, 0.0, 0.0]);
        assert_eq!(wl_push(), 1);
        set(&[1.0, 10.0, f64::INFINITY, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0]);
        assert_eq!(wl_push(), 2);
        let h = wl_finish(1, 0.0, 0.0);
        assert!(h > 0 && wl_len(h) == 2 && wl_live() == 1);
        assert_eq!(wl_at(h, 20.0), OK);
        assert!((srkit_get(0) - 6.0).abs() < 1e-15 && (srkit_get(6) - 18.0).abs() < 1e-14);
        // seen from the origin at t = 26: the ship was at x = 6 when the light left, at t = 20
        assert!((wl_send(h, 26.0, 0.0, 0.0, 0.0) - 20.0).abs() < 1e-12);
        assert_eq!(wl_segment(h, 1), OK);
        assert_eq!(srkit_get(0), 1.0);
        assert_eq!(srkit_get(20), 10.0);
        wl_free(h);
        assert_eq!(wl_live(), 0);
        assert!(wl_tau(h, 0.0).is_nan());
        // the slot is reused
        wl_begin();
        set(&[0.0, 0.0, 10.0, 1.0, 2.0, 3.0]);
        wl_push();
        assert_eq!(wl_finish(1, 0.0, 0.0), h);
        wl_free(h);
    }

    #[test]
    fn rides_and_exact_time() {
        assert_eq!(
            trip_plan(4.2465, 1.0322952409169414, f64::INFINITY, 0.0, 0.0),
            OK
        );
        assert!((srkit_get(6) - 3.54).abs() < 5e-3);
        assert_eq!(trip_plan(1.0, 1.0, f64::INFINITY, 0.5, 0.0), INVALID);
        let put = |slot: usize, s: &str| {
            TEXT.with(|t| t.borrow_mut()[slot][..s.len()].copy_from_slice(s.as_bytes()));
            s.len() as u32
        };
        let read = |len: i32| {
            TEXT.with(|t| String::from_utf8(t.borrow()[2][..len as usize].to_vec()).unwrap())
        };
        let epoch = format!("1{}", "0".repeat(100));
        let (la, lb) = (put(0, &epoch), put(1, "-1"));
        let len = time384_add(la, lb);
        assert_eq!(read(len), "9".repeat(100));
        // one light-second is exactly 1000 ticks
        let len = time512_light_ticks(299_792_458.0, 0.0, 0.0);
        assert_eq!(read(len), "1000");
        assert_eq!(time384_classify(put(0, "1000"), 299_792_458.0, 0.0, 0.0), 0);
        assert_eq!(
            time384_classify(put(0, "1000"), 299_792_459.0, 0.0, 0.0),
            -1
        );
    }
}
