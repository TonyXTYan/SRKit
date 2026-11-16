//! `srkit.Worldline`: ordered segments with a proper-time anchor on each.
//!
//! The class here is the native base. `python/srkit/__init__.py` subclasses it to add the vectorised methods,
//! which turn their arguments into flat buffers of doubles and call the `_..._many` methods below.

use pyo3::{
    exceptions::{PyIndexError, PyTypeError},
    prelude::*,
    types::{PyByteArray, PyIterator, PyTuple},
};
use srkit_core::{frames, worldline::WorldlineBuf};

use crate::{
    convert::{event, finite, fl, invalid, repr_fields, tup, type_name, IntoPy, Tup3, V3},
    segment::{Segment, TauState},
    NoIntersection,
};

/// Radar coordinates of an event for an observer: the observer's proper times when a ping would have to leave
/// to reach the event, and when its echo comes back.
#[pyclass(frozen, get_all, module = "srkit", name = "Radar")]
pub struct Radar {
    /// Proper time at which the observer's light reaches the event if sent then. Minus infinity if no light
    /// from the observer can reach the event.
    tau_send: f64,
    /// Proper time at which the event's light reaches the observer. Infinity if it never does (an event
    /// behind a Rindler horizon has no echo).
    tau_recv: f64,
}

#[pymethods]
impl Radar {
    /// Radar time: the midpoint `(tau_recv + tau_send) / 2`.
    #[getter]
    fn time(&self) -> f64 {
        0.5 * (self.tau_recv + self.tau_send)
    }
    /// Radar distance: half the round trip, `(tau_recv - tau_send) / 2`.
    #[getter]
    fn distance(&self) -> f64 {
        0.5 * (self.tau_recv - self.tau_send)
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        repr_fields!(py, "Radar", self, tau_send, tau_recv)
    }
}

/// Read a buffer of native-endian doubles.
fn doubles(data: &[u8], per_item: usize) -> PyResult<Vec<f64>> {
    if data.len() % (8 * per_item) != 0 {
        return Err(invalid(format!(
            "expected a buffer of {per_item} double(s) per item, got {} bytes",
            data.len()
        )));
    }
    Ok(data
        .chunks_exact(8)
        .map(|c| f64::from_ne_bytes(c.try_into().expect("chunks of eight bytes")))
        .collect())
}

/// Write doubles as a native-endian buffer.
fn buffer<'py>(py: Python<'py>, values: &[f64]) -> Bound<'py, PyByteArray> {
    let mut out = Vec::with_capacity(values.len() * 8);
    for v in values {
        out.extend_from_slice(&v.to_ne_bytes());
    }
    PyByteArray::new(py, &out)
}

/// A worldline: segments in order of their start times, with proper time carried along the chain.
///
/// `Worldline(segments, anchor_t=None, anchor_tau=0.0)` takes `Segment` objects in order of `t0` and sets the
/// body's proper time to `anchor_tau` at coordinate time `anchor_t` (by default the start of the first
/// segment).
///
/// The first segment is taken to run back to the infinite past and the last on to the infinite future, so a
/// worldline answers every query; end it with a `Segment.rest` if the body should stay put. Each segment is in
/// force from its own `t0` to the next segment's `t0`. The first segment needs a finite `t0`, and neither the
/// first nor the last may be a burn, which is not defined outside its own span.
///
/// Natural units, c = 1. A worldline cannot be changed once made.
#[pyclass(frozen, subclass, module = "srkit._srkit", name = "Worldline")]
pub struct Worldline {
    buf: WorldlineBuf,
}

#[pymethods]
impl Worldline {
    #[new]
    #[pyo3(signature = (segments, anchor_t=None, anchor_tau=0.0))]
    fn new(segments: &Bound<'_, PyAny>, anchor_t: Option<f64>, anchor_tau: f64) -> PyResult<Self> {
        let mut segs = Vec::new();
        for item in segments.try_iter()? {
            let item = item?;
            let seg = item.cast::<Segment>().map_err(|_| {
                PyTypeError::new_err(format!(
                    "a worldline is built from Segment objects, got {}",
                    type_name(&item)
                ))
            })?;
            segs.push(seg.get().inner);
        }
        let Some(first) = segs.first() else {
            return Err(invalid("a worldline needs at least one segment"));
        };
        if !first.t0().is_finite() {
            return Err(invalid(
                "the first segment needs a finite t0: proper time is counted from it",
            ));
        }
        if segs.windows(2).any(|w| !(w[0].t0() <= w[1].t0())) {
            return Err(invalid("segments must be in order of t0"));
        }
        // The first and last segments are continued for ever, and a burn has no continuation: before its
        // start the core gives NaN, and after its mass is gone it holds the last state.
        let is_burn = |s: &srkit_core::Segment| matches!(s, srkit_core::Segment::Burn(_));
        if is_burn(first) || segs.last().is_some_and(is_burn) {
            return Err(invalid(
                "a burn cannot be the first or the last segment of a worldline: it is not defined before it \
                 starts or after its mass runs out; put a rest or inertial segment before and after it",
            ));
        }
        let anchor_t = finite("anchor_t", anchor_t.unwrap_or(first.t0()))?;
        finite("anchor_tau", anchor_tau)?;
        Ok(Self {
            buf: WorldlineBuf::new(segs, anchor_t, anchor_tau).py()?,
        })
    }

    fn __len__(&self) -> usize {
        self.buf.segments().len()
    }

    fn __getitem__(&self, index: isize) -> PyResult<Segment> {
        let segs = self.buf.segments();
        let n = segs.len() as isize;
        let i = if index < 0 { index + n } else { index };
        if i < 0 || i >= n {
            return Err(PyIndexError::new_err("segment index out of range"));
        }
        Ok(segs[i as usize].seg.into())
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        self.segments(py)?.as_any().try_iter()
    }

    /// The segments, in order, as a tuple.
    #[getter]
    fn segments<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            self.buf
                .segments()
                .iter()
                .map(|a| Segment::from(a.seg))
                .collect::<Vec<_>>(),
        )
    }

    /// The body's proper time at the start of each segment, as a tuple.
    #[getter]
    fn tau0s<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            self.buf
                .segments()
                .iter()
                .map(|a| a.tau0)
                .collect::<Vec<_>>(),
        )
    }

    /// Whether any segment is light-like (a jump).
    #[getter]
    fn has_null(&self) -> bool {
        self.buf.view().has_null()
    }

    /// Index of the segment in force at coordinate time `t`.
    fn find(&self, t: f64) -> usize {
        self.buf.view().find(t)
    }
    /// Index of the segment in force at proper time `tau`.
    fn find_tau(&self, tau: f64) -> usize {
        self.buf.view().find_tau(tau)
    }

    /// Position and velocity at coordinate time `t`, as `(x, v)`.
    fn at(&self, t: f64) -> (Tup3, Tup3) {
        let s = self.buf.view().at(t);
        (tup(s.x), tup(s.v))
    }
    /// The body's proper time at coordinate time `t`.
    fn tau(&self, t: f64) -> f64 {
        self.buf.view().tau(t)
    }
    /// Coordinate time at proper time `tau`.
    fn t_of_tau(&self, tau: f64) -> f64 {
        self.buf.view().t_of_tau(tau)
    }
    /// The state at proper time `tau`.
    fn at_tau(&self, tau: f64) -> TauState {
        self.buf.view().at_tau(tau).into()
    }
    /// The largest speed on the worldline, as a fraction of c.
    fn vmax(&self) -> f64 {
        self.buf.view().vmax()
    }

    /// The coordinate time at which light from the event `(event_t, event_x)` reaches this worldline, or
    /// `inf` if it never does (a body that outruns the light, such as one accelerating away for ever).
    fn recv(&self, event_t: f64, event_x: V3) -> f64 {
        self.buf.view().recv(event(event_t, event_x))
    }

    /// The coordinate time of the moment on this worldline whose light reaches the event
    /// `(event_t, event_x)`: what an observer there sees of this body. `-inf` if no such moment exists.
    fn send(&self, event_t: f64, event_x: V3) -> f64 {
        self.buf.view().send(event(event_t, event_x))
    }

    /// How this body is seen from the event `(event_t, event_x)`: the time the light left, and the body's
    /// position and velocity then, as `(t, x, v)`. Raises `srkit.NoIntersection` if nothing of it can be seen.
    fn seen_from(&self, event_t: f64, event_x: V3) -> PyResult<(f64, Tup3, Tup3)> {
        let (t, s) = self
            .buf
            .view()
            .seen_from(event(event_t, event_x))
            .ok_or_else(|| {
                NoIntersection::new_err("no light from this worldline reaches the event")
            })?;
        Ok((t, tup(s.x), tup(s.v)))
    }

    /// Check that consecutive segments join: each starts where and when the one before ends, at the same
    /// velocity and proper time. Raises `ValueError` if a join is off by more than `position_tol` (a length)
    /// or `velocity_tol` (a fraction of c). A light-like segment is exempt from the velocity check on both
    /// sides.
    #[pyo3(signature = (position_tol=1e-9, velocity_tol=1e-9))]
    fn check_continuity(&self, position_tol: f64, velocity_tol: f64) -> PyResult<()> {
        self.buf
            .view()
            .check_continuity(position_tol, velocity_tol)
            .py()
    }

    /// As `check_continuity`, returning `True` or `False` instead of raising.
    #[pyo3(signature = (position_tol=1e-9, velocity_tol=1e-9))]
    fn is_continuous(&self, position_tol: f64, velocity_tol: f64) -> bool {
        self.buf
            .view()
            .check_continuity(position_tol, velocity_tol)
            .is_ok()
    }

    /// Radar coordinates of the event `(event_t, event_x)` for the observer on this worldline. A side is
    /// infinite where the light cone never meets the worldline.
    fn radar(&self, event_t: f64, event_x: V3) -> Radar {
        let r = frames::radar(&self.buf.view(), event(event_t, event_x));
        Radar {
            tau_send: r.tau_send,
            tau_recv: r.tau_recv,
        }
    }

    /// Light-cone coordinates of the event `(event_t, event_x)`: the observer's proper time when the event is
    /// seen, and the event's distance in the observer's rest frame at that moment, as `(tau, distance)`.
    /// Raises `srkit.NoIntersection` if the event is never seen.
    fn light_cone_coordinates(&self, event_t: f64, event_x: V3) -> PyResult<(f64, f64)> {
        frames::light_cone_coordinates(&self.buf.view(), event(event_t, event_x)).ok_or_else(|| {
            NoIntersection::new_err("light from the event never reaches this worldline")
        })
    }

    /// Positions and velocities at many times. `ts` is a buffer of native-endian doubles; the result is two
    /// buffers of three doubles per time.
    fn _at_many<'py>(
        &self,
        py: Python<'py>,
        ts: &[u8],
    ) -> PyResult<(Bound<'py, PyByteArray>, Bound<'py, PyByteArray>)> {
        let ts = doubles(ts, 1)?;
        let buf = &self.buf;
        // The batch touches no Python object, so other Python threads may run meanwhile.
        let (xs, vs) = py.detach(|| {
            let view = buf.view();
            let (mut xs, mut vs) = (
                Vec::with_capacity(ts.len() * 3),
                Vec::with_capacity(ts.len() * 3),
            );
            for t in ts {
                let s = view.at(t);
                xs.extend_from_slice(&s.x.to_array());
                vs.extend_from_slice(&s.v.to_array());
            }
            (xs, vs)
        });
        Ok((buffer(py, &xs), buffer(py, &vs)))
    }

    /// Proper times at many coordinate times, buffer to buffer.
    fn _tau_many<'py>(&self, py: Python<'py>, ts: &[u8]) -> PyResult<Bound<'py, PyByteArray>> {
        let ts = doubles(ts, 1)?;
        let buf = &self.buf;
        let out: Vec<f64> = py.detach(|| {
            let view = buf.view();
            ts.into_iter().map(|t| view.tau(t)).collect()
        });
        Ok(buffer(py, &out))
    }

    /// Light-cone times for many events. `events` is a buffer of four doubles per event (t, x, y, z).
    fn _cone_many<'py>(
        &self,
        py: Python<'py>,
        events: &[u8],
        send: bool,
    ) -> PyResult<Bound<'py, PyByteArray>> {
        let events = doubles(events, 4)?;
        let buf = &self.buf;
        let out: Vec<f64> = py.detach(|| {
            let view = buf.view();
            events
                .chunks_exact(4)
                .map(|e| {
                    let ev = srkit_core::Event::new(e[0], srkit_core::Vec3::new(e[1], e[2], e[3]));
                    if send {
                        view.send(ev)
                    } else {
                        view.recv(ev)
                    }
                })
                .collect()
        });
        Ok(buffer(py, &out))
    }

    fn __repr__(&self) -> String {
        let segs = self.buf.segments();
        let kinds: Vec<&str> = segs.iter().map(|a| Segment::from(a.seg).kind()).collect();
        format!(
            "Worldline({} segment{}: {}; t0={}, tau0={})",
            segs.len(),
            if segs.len() == 1 { "" } else { "s" },
            kinds.join(", "),
            fl(segs[0].seg.t0()),
            fl(segs[0].tau0)
        )
    }
}
