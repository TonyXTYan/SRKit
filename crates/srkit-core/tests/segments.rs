//! Segments, worldlines and light cones: the identities that must hold whatever the numbers.
use proptest::prelude::*;
use srkit_core::{
    frames, kinematics,
    kinematics::Event,
    rindler::RindlerHorizon,
    segment::{Burn, Hyperbolic, Segment, Steered},
    worldline::{anchor, Anchored, Worldline, WorldlineBuf},
    Vec3,
};

fn unit(a: f64, b: f64) -> Vec3 {
    Vec3::new(a.cos() * b.cos(), a.sin() * b.cos(), b.sin())
}

/// Segments of every kind with parameters in the ranges the game uses.
fn any_segment() -> impl Strategy<Value = Segment> {
    let dir = (0.0..6.28f64, -1.5..1.5f64).prop_map(|(a, b)| unit(a, b));
    let point =
        (-50.0..50.0f64, -50.0..50.0f64, -50.0..50.0f64).prop_map(|(x, y, z)| Vec3::new(x, y, z));
    prop_oneof![
        (point.clone(), -100.0..100.0f64).prop_map(|(p, t0)| Segment::rest(t0, t0 + 30.0, p)),
        (point.clone(), dir.clone(), 0.0..0.999f64, -100.0..100.0f64)
            .prop_map(|(p, n, b, t0)| Segment::inertial(t0, t0 + 30.0, p, n * b)),
        (
            point.clone(),
            dir.clone(),
            0.01..3.0f64,
            0.0..4.0f64,
            -100.0..100.0f64
        )
            .prop_map(|(p, n, a, phi0, t0)| Segment::hyperbolic(
                t0,
                t0 + 30.0,
                p,
                n,
                a,
                phi0
            )),
        (
            point.clone(),
            dir.clone(),
            0.01..3.0f64,
            0.5..4.0f64,
            -100.0..100.0f64
        )
            .prop_map(|(p, n, a, phi0, t0)| {
                // braking, ending before the turnaround
                Segment::hyperbolic(t0, t0 + 0.8 * phi0.sinh() / a, p, n, -a, phi0)
            }),
        (
            point.clone(),
            0.1..40.0f64,
            0.01..0.99f64,
            0.0..6.28f64,
            -100.0..100.0f64
        )
            .prop_map(|(c, r, beta, a0, t0)| Segment::circular(
                t0,
                t0 + 30.0,
                c,
                r,
                beta / r,
                Vec3::X,
                Vec3::Y,
                a0
            )),
        (
            point.clone(),
            dir.clone(),
            dir.clone(),
            0.0..0.99f64,
            0.01..3.0f64,
            -100.0..100.0f64
        )
            .prop_map(|(p, n, m, b, a, t0)| {
                Segment::Steered(Steered::new(t0, t0 + 30.0, p, n * b, m, a).unwrap())
            }),
        (
            point,
            dir,
            0.0..2.0f64,
            0.05..2.0f64,
            0.05..1.0f64,
            -100.0..100.0f64
        )
            .prop_map(|(p, n, phi0, a0, w, t0)| {
                let g = Burn::new(t0, f64::NAN, p, n, phi0, a0, w).unwrap();
                let mut seg = Segment::Burn(g);
                let t1 = seg.t_of_tau(0.8 * g.t_burn);
                seg.set_t1(t1);
                seg
            }),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(400))]

    /// Proper time and coordinate time are inverses, and the two ways of evaluating a state agree.
    #[test]
    fn proper_time_inverts(seg in any_segment(), f in 0.0..1.0f64) {
        let t = seg.t0() + f * (seg.t1() - seg.t0());
        let s = seg.dtau(t);
        prop_assert!(s >= 0.0 && s <= (t - seg.t0()) * (1.0 + 1e-12));
        let back = seg.t_of_tau(s);
        prop_assert!((back - t).abs() <= 1e-9 * (1.0 + t.abs()), "t {t} -> tau {s} -> t {back}");
        let (a, b) = (seg.at(t), seg.at_tau(s));
        prop_assert!(a.x.dist(b.x) <= 1e-8 * (1.0 + a.x.norm()));
        prop_assert!(a.v.dist(b.v) <= 1e-8);
    }

    /// Nothing outruns light, and proper time never runs faster than coordinate time.
    #[test]
    fn slower_than_light(seg in any_segment(), f in 0.0..1.0f64, h in 1e-4..1e-2f64) {
        let t = seg.t0() + f * (seg.t1() - seg.t0());
        let (a, b) = (seg.at(t), seg.at(t + h));
        prop_assert!(a.v.norm() < 1.0);
        prop_assert!(a.x.dist(b.x) <= h * (1.0 + 1e-9));
        // the velocity is the derivative of the position
        let numeric = (seg.at(t + 1e-6).x - seg.at(t - 1e-6).x) / 2e-6;
        prop_assert!(numeric.dist(a.v) < 1e-4, "{numeric:?} vs {:?}", a.v);
        // dtau/dt = 1/gamma
        let rate = (seg.dtau(t + 1e-6) - seg.dtau(t - 1e-6)) / 2e-6;
        prop_assert!((rate - (1.0 - a.v.norm_squared()).sqrt()).abs() < 1e-4);
    }

    /// A light-cone root really is on the light cone, and what you see is what was sent.
    #[test]
    fn light_cone_roots(seg in any_segment(), ex in -80.0..80.0f64, ey in -80.0..80.0f64, ez in -80.0..80.0f64, et in -150.0..150.0f64) {
        let event = Event::new(et, Vec3::new(ex, ey, ez));
        if let Some(t) = seg.send(event, seg.t0(), seg.t1()) {
            let r = seg.at(t).x.dist(event.x);
            prop_assert!(t <= event.t && t >= seg.t0() && t <= seg.t1());
            prop_assert!((event.t - t - r).abs() <= 1e-9 * (1.0 + r + t.abs()), "send off the cone by {}", event.t - t - r);
            // light leaving the body then arrives back at the event
            let emitted = Event::new(t, seg.at(t).x);
            prop_assert!(emitted.can_reach(event) || (event.t - t - r).abs() < 1e-9 * (1.0 + r));
        }
        if let Some(t) = seg.recv(event, seg.t0(), seg.t1()) {
            let r = seg.at(t).x.dist(event.x);
            prop_assert!(t >= event.t);
            if t > seg.t0().max(event.t) {
                prop_assert!((t - event.t - r).abs() <= 1e-9 * (1.0 + r + t.abs()), "recv off the cone by {}", t - event.t - r);
            }
        }
    }

    /// The interval between two events does not depend on the frame.
    #[test]
    fn interval_is_invariant(w in (-3.0..3.0f64, -3.0..3.0f64, -3.0..3.0f64), dt in -10.0..10.0f64, dx in (-10.0..10.0f64, -10.0..10.0f64, -10.0..10.0f64)) {
        let (w, dx) = (Vec3::new(w.0, w.1, w.2), Vec3::new(dx.0, dx.1, dx.2));
        let (t2, x2) = kinematics::boost_rapidity(w, dt, dx);
        let (s1, s2) = (dt * dt - dx.norm_squared(), t2 * t2 - x2.norm_squared());
        let scale = (dt * dt + dx.norm_squared()) * (2.0 * w.norm()).cosh();
        prop_assert!((s1 - s2).abs() <= 1e-12 * (1.0 + scale));
        // and the inverse boost undoes it
        let (t3, x3) = kinematics::boost_rapidity(-w, t2, x2);
        prop_assert!((t3 - dt).abs() <= 1e-12 * (1.0 + scale) && x3.dist(dx) <= 1e-12 * (1.0 + scale));
    }

    /// Thrust at an angle agrees with thrust along the motion when the angle is zero or a half turn.
    #[test]
    fn steered_matches_hyperbolic(a in 0.01..3.0f64, phi0 in 0.0..3.0f64, dt in 0.0..20.0f64, toward in any::<bool>()) {
        let n = unit(0.7, 0.3);
        let sign = if toward { 1.0 } else { -1.0 };
        let dt = if toward { dt } else { dt.min(0.8 * phi0.sinh() / a) };
        let straight = Segment::hyperbolic(0.0, 50.0, Vec3::ZERO, n, sign * a, phi0);
        let steered = Segment::Steered(Steered::from_proper_velocity(0.0, 50.0, Vec3::ZERO, n * phi0.sinh(), n * sign, a).unwrap());
        let (p, q) = (straight.at(dt), steered.at(dt));
        prop_assert!(p.x.dist(q.x) <= 1e-9 * (1.0 + p.x.norm()), "{:?} vs {:?}", p.x, q.x);
        prop_assert!(p.v.dist(q.v) <= 1e-9 * phi0.cosh().powi(2));
        prop_assert!((straight.dtau(dt) - steered.dtau(dt)).abs() <= 1e-9 * (1.0 + dt) * phi0.cosh());
    }

    /// A steered burn feels constant proper acceleration: the four-acceleration has magnitude alpha.
    #[test]
    fn steered_holds_its_acceleration(v in (0.0..6.28f64, -1.5..1.5f64, 0.0..0.95f64), m in (0.0..6.28f64, -1.5..1.5f64), a in 0.05..2.0f64, s in 0.01..3.0f64) {
        let seg = Segment::Steered(Steered::new(0.0, 100.0, Vec3::ZERO, unit(v.0, v.1) * v.2, unit(m.0, m.1), a).unwrap());
        let h = 1e-4;
        let four_velocity = |tau: f64| { let st = seg.at_tau(tau); let g = 1.0 / (1.0 - st.v.norm_squared()).sqrt(); (g, st.v * g) };
        let (g1, u1) = four_velocity(s - h);
        let (g2, u2) = four_velocity(s + h);
        let (ag, au) = ((g2 - g1) / (2.0 * h), (u2 - u1) / (2.0 * h));
        let magnitude = (au.norm_squared() - ag * ag).sqrt();
        prop_assert!((magnitude / a - 1.0).abs() < 1e-4, "proper acceleration {magnitude} vs {a}");
    }

    /// A constant-thrust burn obeys the rocket equation, and its clocks agree with a direct integration.
    #[test]
    fn burn_obeys_the_rocket_equation(phi0 in 0.0..1.5f64, a0 in 0.05..2.0f64, w in 0.05..1.0f64, f in 0.05..0.9f64, against in any::<bool>()) {
        let g = Burn::new(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, phi0, if against { -a0 } else { a0 }, w).unwrap();
        let seg = Segment::Burn(g);
        let s = f * g.t_burn;
        let st = seg.at_tau(s);
        let dphi = (st.phi - phi0).abs();
        prop_assert!((srkit_core::rocket::mass_ratio(dphi, w) * g.mass_fraction(s) - 1.0).abs() < 1e-11);
        // Simpson's rule on dt = cosh(phi) ds
        let n = 2000;
        let h = s / n as f64;
        let c = |k: usize| seg.at_tau(k as f64 * h).phi.cosh();
        let sum: f64 = (1..n).map(|k| if k % 2 == 1 { 4.0 * c(k) } else { 2.0 * c(k) }).sum::<f64>() + c(0) + c(n);
        prop_assert!((sum * h / 3.0 / st.t - 1.0).abs() < 1e-6, "{} vs {}", sum * h / 3.0, st.t);
        // the acceleration at the start is the one asked for
        prop_assert!((g.acceleration(0.0) / a0 - 1.0).abs() < 1e-14);
    }
}

fn boost_coast_brake() -> WorldlineBuf {
    // from rest at the origin: 1 unit of proper time boosting at a = 1, coast, brake to rest, then stay
    let a = 1.0;
    let boost = Segment::hyperbolic(0.0, 1.0_f64.sinh() / a, Vec3::ZERO, Vec3::X, a, 0.0);
    let s1 = boost.at_tau(1.0);
    let coast = Segment::inertial(s1.t, s1.t + 5.0, s1.x, s1.v);
    let s2 = coast.at(s1.t + 5.0);
    let brake = Segment::hyperbolic(
        s1.t + 5.0,
        s1.t + 5.0 + 1.0_f64.sinh() / a,
        s2.x,
        Vec3::X,
        -a,
        1.0,
    );
    let s3 = brake.at_tau(1.0);
    let stay = Segment::rest(s3.t, f64::INFINITY, s3.x);
    WorldlineBuf::new(
        [
            Segment::rest(-1000.0, 0.0, Vec3::ZERO),
            boost,
            coast,
            brake,
            stay,
        ],
        0.0,
        0.0,
    )
    .unwrap()
}

#[test]
fn worldline_joins_and_clocks() {
    let buf = boost_coast_brake();
    let w = buf.view();
    w.check_continuity(1e-12, 1e-12).unwrap();
    assert_eq!(w.len(), 5);
    // proper time: 1 boosting + 5 / gamma coasting + 1 braking
    let gamma = 1.0_f64.cosh();
    let end = w.segments()[4].seg.t0();
    assert!((w.tau(end) - (2.0 + 5.0 / gamma)).abs() < 1e-13);
    // before launch the traveller's clock runs with coordinate time
    assert_eq!(w.tau(-3.0), -3.0);
    // by proper time and by coordinate time agree all along
    for k in 0..200 {
        let t = -2.0 + 0.06 * k as f64;
        let s = w.at_tau(w.tau(t));
        assert!((s.t - t).abs() < 1e-12 && s.x.dist(w.at(t).x) < 1e-12);
        assert!((w.t_of_tau(w.tau(t)) - t).abs() < 1e-12);
    }
    assert!((w.vmax() - 1.0_f64.tanh()).abs() < 1e-15);
}

#[test]
fn what_is_seen_and_when_signals_arrive() {
    let buf = boost_coast_brake();
    let w = buf.view();
    let end = w.segments()[4].seg.t0();
    let x_end = w.at(end).x;
    // an observer at the origin watches the ship; a signal from the origin chases it
    for k in 0..300 {
        let t = -1.0 + 0.1 * k as f64;
        let watcher = Event::new(t, Vec3::new(0.0, 3.0, 0.0));
        let te = w.send(watcher);
        assert!(te <= t);
        let r = w.at(te).x.dist(watcher.x);
        assert!(
            (t - te - r).abs() < 1e-11,
            "send off the cone by {}",
            t - te - r
        );
        let tr = w.recv(watcher);
        assert!(tr >= t && (tr - t - w.at(tr).x.dist(watcher.x)).abs() < 1e-11);
        // radar: the two cones bracket the observer's "now" for the event
        let radar = frames::radar(&w, watcher);
        assert!(radar.tau_send <= radar.tau_recv);
    }
    // long after the ride the ship is seen at rest at its destination, at the light delay
    let late = Event::new(end + 1000.0, Vec3::ZERO);
    assert!((w.send(late) - (late.t - x_end.norm())).abs() < 1e-10);
    let (te, state) = w.seen_from(late).unwrap();
    assert!(te > end && state.v == Vec3::ZERO);
}

#[test]
fn anchoring_matches_the_owned_worldline() {
    let buf = boost_coast_brake();
    let mut segs: Vec<Anchored> = buf
        .segments()
        .iter()
        .map(|a| Anchored::new(a.seg, 0.0))
        .collect();
    anchor(&mut segs, 0.0, 0.0);
    assert_eq!(&segs[..], buf.segments());
    // anchoring to a later moment shifts every clock by the same amount
    anchor(&mut segs, 2.0, 100.0);
    let w = Worldline::new(&segs).unwrap();
    let shift = w.tau(2.0) - buf.view().tau(2.0);
    assert!((w.tau(2.0) - 100.0).abs() < 1e-12);
    for t in [-5.0, 0.3, 4.0, 9.0] {
        assert!((w.tau(t) - buf.view().tau(t) - shift).abs() < 1e-12);
    }
    assert!(Worldline::new(&[]).is_err());
}

#[test]
fn a_jump_takes_light_time_and_no_proper_time() {
    let target = Vec3::new(3.0, 4.0, 0.0);
    let jump = Segment::Null(srkit_core::segment::Null::between(10.0, Vec3::ZERO, target));
    assert_eq!(jump.t1(), 15.0);
    assert_eq!(jump.dtau(12.0), 0.0);
    let after = jump.at_tau(0.0);
    assert!(after.jump && after.t == 15.0 && after.x == target);
    let buf = WorldlineBuf::new(
        [
            Segment::rest(0.0, 10.0, Vec3::ZERO),
            jump,
            Segment::rest(15.0, f64::INFINITY, target),
        ],
        10.0,
        0.0,
    )
    .unwrap();
    let w = buf.view();
    assert!(w.has_null());
    w.check_continuity(1e-12, 1e-12).unwrap();
    assert_eq!(w.tau(15.0), 0.0);
    assert_eq!(w.tau(16.0), 1.0);
    // an observer at the target sees the traveller leave at the instant they arrive: the image rides the jump
    let seen = w.send(Event::new(15.0, target));
    assert!((10.0..=15.0).contains(&seen));
    // an observer off to the side sees the jump compressed but in order
    let side = Vec3::new(0.0, 10.0, 0.0);
    let mut last = f64::NEG_INFINITY;
    for k in 0..100 {
        let te = w.send(Event::new(20.0 + 0.2 * k as f64, side));
        assert!(te >= last);
        last = te;
    }
    // light from an event reaches the far side of the jump no earlier than the jump allows
    let tr = w.recv(Event::new(9.0, Vec3::new(-1.0, 0.0, 0.0)));
    assert!((tr - 10.0).abs() < 1e-12);
}

#[test]
fn rindler_horizon_hides_what_falls_behind() {
    let a = 2.0;
    let seg = Hyperbolic::new(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, a, 0.0);
    let h = RindlerHorizon::of(&seg);
    assert!((h.distance - 0.5).abs() < 1e-15);
    assert_eq!(h.pivot.x, Vec3::new(-0.5, 0.0, 0.0));
    let ship = Segment::Hyperbolic(seg);
    // a signal from the launch point sent before t = 1/a arrives; one sent after never does
    let last = h.last_signal_time(Vec3::ZERO);
    assert!((last - 0.5).abs() < 1e-15);
    let before = Event::new(last - 1e-3, Vec3::ZERO);
    let after = Event::new(last + 1e-3, Vec3::ZERO);
    assert!(h.can_signal(before) && !h.can_signal(after));
    let t = ship.recv(before, 0.0, f64::INFINITY).unwrap();
    assert!(t.is_finite() && (t - before.t - ship.at(t).x.dist(before.x)).abs() < 1e-9 * t);
    assert_eq!(ship.recv(after, 0.0, f64::INFINITY), None);
    // the same through a worldline: the signal never arrives
    let buf = WorldlineBuf::new([ship], 0.0, 0.0).unwrap();
    assert_eq!(buf.view().recv(after), f64::INFINITY);
    assert!(frames::radar(&buf.view(), after).tau_recv.is_infinite());
}
