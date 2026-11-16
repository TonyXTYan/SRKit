//! `srkit.Segment`: one closed-form piece of a worldline, and the small records its methods return.

use pyo3::{prelude::*, types::PyDict};
use srkit_core::{
    rindler,
    segment::{self as core, Burn, Null, Steered},
    Vec3,
};

use crate::{
    convert::{
        event, finite, finite_vec, fl, invalid, repr_fields, repr_values, subluminal, tup,
        unit_vec, IntoPy, Tup3, V3,
    },
    NoIntersection,
};

/// The state of a body after a given proper time: what `at_tau` returns.
#[pyclass(frozen, get_all, module = "srkit", name = "TauState")]
pub struct TauState {
    /// Coordinate time.
    t: f64,
    /// Position `(x, y, z)`.
    x: Tup3,
    /// Velocity as a fraction of c. After a light-like jump it is zero.
    v: Tup3,
    /// Rapidity: signed along the segment's axis on a hyperbolic or burn segment, the speed's rapidity
    /// otherwise.
    phi: f64,
    /// Proper acceleration (signed along the axis on a hyperbolic or burn segment, the magnitude otherwise).
    acc: f64,
    /// True for the state after a light-like jump, which takes no proper time.
    jump: bool,
}

impl From<core::TauState> for TauState {
    fn from(s: core::TauState) -> Self {
        Self {
            t: s.t,
            x: tup(s.x),
            v: tup(s.v),
            phi: s.phi,
            acc: s.acc,
            jump: s.jump,
        }
    }
}

#[pymethods]
impl TauState {
    fn __repr__(&self, py: Python<'_>) -> String {
        repr_fields!(py, "TauState", self, t, x, v, phi, acc, jump)
    }
}

/// The Rindler horizon of a hyperbolic segment continued for ever.
///
/// A ship that holds its proper acceleration outruns light from part of spacetime: the region beyond a plane
/// that trails it at distance `1 / |alpha|` and closes on it at the speed of light.
#[pyclass(frozen, module = "srkit", name = "RindlerHorizon")]
pub struct RindlerHorizon {
    inner: rindler::RindlerHorizon,
}

#[pymethods]
impl RindlerHorizon {
    /// Coordinate time of the event the hyperbola is centred on.
    #[getter]
    fn pivot_t(&self) -> f64 {
        self.inner.pivot.t
    }
    /// Position of the event the hyperbola is centred on.
    #[getter]
    fn pivot_x(&self) -> Tup3 {
        tup(self.inner.pivot.x)
    }
    /// Unit vector along the proper acceleration (the direction the ship is being pushed).
    #[getter]
    fn n(&self) -> Tup3 {
        tup(self.inner.n)
    }
    /// Distance from the ship to the horizon in the ship's own rest frame: `1 / |alpha|`.
    #[getter]
    fn distance(&self) -> f64 {
        self.inner.distance
    }
    /// How far the event lies inside the region whose light still reaches the ship. Positive: the light
    /// arrives. Zero: on the horizon. Negative: it never arrives while the ship keeps accelerating.
    fn margin(&self, event_t: f64, event_x: V3) -> f64 {
        self.inner.margin(event(event_t, event_x))
    }
    /// Whether light from the event ever reaches the ship.
    fn can_signal(&self, event_t: f64, event_x: V3) -> bool {
        self.inner.can_signal(event(event_t, event_x))
    }
    /// The last coordinate time at which a signal sent from the point `x` still reaches the ship.
    fn last_signal_time(&self, x: V3) -> f64 {
        self.inner.last_signal_time(x.0)
    }
    /// How far the event lies inside the region the ship's own signals can reach. Positive means some signal
    /// from the ship arrives there.
    fn margin_to_be_seen(&self, event_t: f64, event_x: V3) -> f64 {
        self.inner.margin_to_be_seen(event(event_t, event_x))
    }
    fn __repr__(&self) -> String {
        format!(
            "RindlerHorizon(pivot_t={}, pivot_x={}, n={}, distance={})",
            fl(self.inner.pivot.t),
            fmt3(self.inner.pivot.x),
            fmt3(self.inner.n),
            fl(self.inner.distance)
        )
    }
}

fn fmt3(v: Vec3) -> String {
    format!("({}, {}, {})", fl(v.x), fl(v.y), fl(v.z))
}

/// One closed-form piece of a worldline, covering coordinate times `t0` to `t1` (`t1` may be infinite).
///
/// Build one with a constructor: `Segment.rest`, `Segment.inertial`, `Segment.hyperbolic`, `Segment.steered`,
/// `Segment.burn`, `Segment.circular`, `Segment.null` or `Segment.jump`. A segment cannot be changed once made.
///
/// Natural units, c = 1: times and lengths share one unit (years and light-years, say), velocities are
/// fractions of c, rapidities have no unit and proper accelerations are rates of change of rapidity per unit of
/// time. Positions, velocities and directions are `(x, y, z)`; any sequence of three numbers is accepted and
/// tuples of floats are returned.
///
/// Evaluation does not clamp: a segment is its own analytic continuation outside `[t0, t1]`. A `Worldline`
/// decides which segment a time belongs to.
#[pyclass(frozen, eq, module = "srkit", name = "Segment")]
#[derive(Clone, Copy, PartialEq)]
pub struct Segment {
    pub(crate) inner: core::Segment,
}

impl From<core::Segment> for Segment {
    fn from(inner: core::Segment) -> Self {
        Self { inner }
    }
}

/// Check the time span of a new segment.
fn span(t0: f64, t1: f64) -> PyResult<()> {
    if t0.is_nan() || t1.is_nan() {
        return Err(invalid("t0 and t1 must not be NaN"));
    }
    if t1 < t0 {
        return Err(invalid(format!(
            "t1 ({}) must not be before t0 ({})",
            fl(t1),
            fl(t0)
        )));
    }
    Ok(())
}

/// Check the time span of a new segment that moves: its start is where the motion is pinned, so it must be
/// finite. (A rest may reach back for ever.)
fn moving_span(t0: f64, t1: f64) -> PyResult<()> {
    span(finite("t0", t0)?, t1)
}

/// The coordinate times between which a burn is defined: from its start to a hair before its mass is gone.
/// The core gives NaN before the start and holds the last state after the end, so nothing outside this range
/// is physics. `None` for the other kinds, which continue analytically both ways.
fn burn_domain(seg: &core::Segment) -> Option<(f64, f64)> {
    match seg {
        core::Segment::Burn(g) => Some((g.t0, seg.t_of_tau(g.t_burn * (1.0 - f64::EPSILON)))),
        _ => None,
    }
}

/// Reject an end time after a burn's mass has run out.
fn check_burn_end(seg: &core::Segment) -> PyResult<()> {
    if let Some((_, limit)) = burn_domain(seg) {
        if !(seg.t1() <= limit) {
            return Err(invalid(format!(
                "the burn must end before the mass runs out: t1 must not be after {}, got {}",
                fl(limit),
                fl(seg.t1())
            )));
        }
    }
    Ok(())
}

impl Segment {
    /// A search range for a light-cone query, kept within the times the segment is defined for.
    fn range(&self, lo: f64, hi: f64) -> (f64, f64) {
        match burn_domain(&self.inner) {
            Some((start, end)) => (lo.max(start), hi.min(end)),
            None => (lo, hi),
        }
    }
}

#[pymethods]
impl Segment {
    /// At rest at the point `p` from `t0` to `t1`.
    #[staticmethod]
    fn rest(t0: f64, t1: f64, p: V3) -> PyResult<Self> {
        span(t0, t1)?;
        Ok(core::Segment::rest(t0, t1, finite_vec("p", p.0)?).into())
    }

    /// Constant velocity `v` (a fraction of c, below 1) from the point `p0` at `t0`.
    #[staticmethod]
    fn inertial(t0: f64, t1: f64, p0: V3, v: V3) -> PyResult<Self> {
        moving_span(t0, t1)?;
        Ok(core::Segment::inertial(t0, t1, finite_vec("p0", p0.0)?, subluminal("v", v.0)?).into())
    }

    /// Constant proper acceleration `alpha` (signed) along the unit vector `n`, from the point `p0` at `t0`
    /// where the rapidity along `n` is `phi0`.
    ///
    /// After proper time `s` the rapidity is `phi0 + alpha s`; the velocity stays along `n`. Rapidity is stored
    /// directly, so a ride at a Lorentz factor of a million loses no precision. For 1 g in years and
    /// light-years pass `alpha = srkit.units.G0_PER_YEAR`.
    #[staticmethod]
    #[pyo3(signature = (t0, t1, p0, n, alpha, phi0=0.0))]
    fn hyperbolic(t0: f64, t1: f64, p0: V3, n: V3, alpha: f64, phi0: f64) -> PyResult<Self> {
        moving_span(t0, t1)?;
        if !(alpha.is_finite() && alpha != 0.0) {
            return Err(invalid(format!(
                "alpha must be finite and not zero, got {}",
                fl(alpha)
            )));
        }
        Ok(core::Segment::hyperbolic(
            t0,
            t1,
            finite_vec("p0", p0.0)?,
            unit_vec("n", n.0)?,
            alpha,
            finite("phi0", phi0)?,
        )
        .into())
    }

    /// Constant proper acceleration `alpha > 0` in a direction fixed aboard, at any angle to the motion.
    ///
    /// The ship starts at `p0` at `t0` with velocity `v0` (a fraction of c) and thrusts along `thrust`, a
    /// direction in its starting rest frame (it need not be a unit vector). With `proper_velocity=True`, `v0`
    /// is instead the starting proper velocity `gamma v`, which keeps its precision near c.
    #[staticmethod]
    #[pyo3(signature = (t0, t1, p0, v0, thrust, alpha, *, proper_velocity=false))]
    fn steered(
        t0: f64,
        t1: f64,
        p0: V3,
        v0: V3,
        thrust: V3,
        alpha: f64,
        proper_velocity: bool,
    ) -> PyResult<Self> {
        moving_span(t0, t1)?;
        let seg = if proper_velocity {
            Steered::from_proper_velocity(t0, t1, p0.0, v0.0, thrust.0, alpha)
        } else {
            Steered::new(t0, t1, p0.0, v0.0, thrust.0, alpha)
        }
        .py()?;
        Ok(core::Segment::Steered(seg).into())
    }

    /// Constant thrust along (or against) `n` while the ship's mass falls: a real rocket.
    ///
    /// The ship starts at `p0` at `t0` with rapidity `phi0` along `n`. `alpha0` is the proper acceleration at
    /// the start (negative to thrust against `n`) and `w` the exhaust speed as a fraction of c (1 for a photon
    /// rocket). The mass would run out after proper time `w / |alpha0|`, so the segment must end before then:
    /// `t1` must be finite and early enough.
    ///
    /// Unlike the other kinds, a burn means something only from its start to the moment its mass is gone.
    /// Outside that span the formula still gives numbers (a heavier, slower ship before `t0`), but `recv` and
    /// `send` search no further than the span, and a `Worldline` needs another segment before a burn and after
    /// it.
    #[staticmethod]
    #[pyo3(signature = (t0, t1, p0, n, phi0, alpha0, w=1.0))]
    fn burn(t0: f64, t1: f64, p0: V3, n: V3, phi0: f64, alpha0: f64, w: f64) -> PyResult<Self> {
        moving_span(t0, t1)?;
        let burn = Burn::new(t0, t1, finite_vec("p0", p0.0)?, n.0, phi0, alpha0, w).py()?;
        let seg = core::Segment::Burn(burn);
        check_burn_end(&seg)?;
        Ok(seg.into())
    }

    /// A circle of radius `r` about `centre` at angular velocity `omega` (radians per unit of coordinate time).
    ///
    /// The circle lies in the plane of the orthonormal pair `(e1, e2)`, by default the x and y axes, and the
    /// body is at angle `a0 + omega (t - t0)` from `e1` towards `e2`. The speed `r |omega|` must be below 1.
    #[staticmethod]
    #[pyo3(signature = (t0, t1, centre, r, omega, e1=None, e2=None, a0=0.0))]
    #[allow(clippy::too_many_arguments)]
    fn circular(
        t0: f64,
        t1: f64,
        centre: V3,
        r: f64,
        omega: f64,
        e1: Option<V3>,
        e2: Option<V3>,
        a0: f64,
    ) -> PyResult<Self> {
        moving_span(t0, t1)?;
        let e1 = unit_vec("e1", e1.map_or(Vec3::X, |v| v.0))?;
        let e2 = unit_vec("e2", e2.map_or(Vec3::Y, |v| v.0))?;
        if e1.dot(e2).abs() > 1e-9 {
            return Err(invalid("e1 and e2 must be at right angles"));
        }
        if !(r.is_finite() && r > 0.0) {
            return Err(invalid(format!(
                "r must be finite and positive, got {}",
                fl(r)
            )));
        }
        finite("omega", omega)?;
        if !((r * omega).abs() < 1.0) {
            return Err(invalid(format!(
                "the speed r |omega| must be below 1 (the speed of light), got {}",
                fl((r * omega).abs())
            )));
        }
        Ok(core::Segment::circular(
            t0,
            t1,
            finite_vec("centre", centre.0)?,
            r,
            omega,
            e1,
            e2,
            finite("a0", a0)?,
        )
        .into())
    }

    /// A light-like run from `p0` along the unit vector `n`, covering a distance `t1 - t0`. No proper time
    /// passes on it. Use it for a signal, or for a jump that moves at c.
    #[staticmethod]
    fn null(t0: f64, t1: f64, p0: V3, n: V3) -> PyResult<Self> {
        span(finite("t0", t0)?, finite("t1", t1)?)?;
        Ok(core::Segment::null(t0, t1, finite_vec("p0", p0.0)?, unit_vec("n", n.0)?).into())
    }

    /// A light-like jump from `p0` at `t0` to the point `p1`, arriving after the light-travel time between
    /// them. The end is pinned exactly on `p1`.
    #[staticmethod]
    fn jump(t0: f64, p0: V3, p1: V3) -> PyResult<Self> {
        finite("t0", t0)?;
        Ok(core::Segment::Null(Null::between(
            t0,
            finite_vec("p0", p0.0)?,
            finite_vec("p1", p1.0)?,
        ))
        .into())
    }

    /// Start of the segment in coordinate time.
    #[getter]
    fn t0(&self) -> f64 {
        self.inner.t0()
    }
    /// End of the segment in coordinate time (may be infinite).
    #[getter]
    fn t1(&self) -> f64 {
        self.inner.t1()
    }
    /// The kind of motion: 'rest', 'inertial', 'hyperbolic', 'steered', 'burn', 'circular' or 'null'.
    #[getter]
    pub(crate) fn kind(&self) -> &'static str {
        match self.inner {
            core::Segment::Rest(_) => "rest",
            core::Segment::Inertial(_) => "inertial",
            core::Segment::Hyperbolic(_) => "hyperbolic",
            core::Segment::Steered(_) => "steered",
            core::Segment::Burn(_) => "burn",
            core::Segment::Circular(_) => "circular",
            core::Segment::Null(_) => "null",
        }
    }
    /// Whether the segment is light-like.
    #[getter]
    fn is_null(&self) -> bool {
        self.inner.is_null()
    }
    /// The numbers that define the segment besides `t0` and `t1`, as a new dict keyed by name. The keys depend
    /// on the kind: for example `p0`, `n`, `alpha` and `phi0` on a hyperbolic segment.
    #[getter]
    fn params<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        match self.inner {
            core::Segment::Rest(g) => {
                d.set_item("p", tup(g.p))?;
            }
            core::Segment::Inertial(g) => {
                d.set_item("p0", tup(g.p0))?;
                d.set_item("v", tup(g.v))?;
            }
            core::Segment::Hyperbolic(g) => {
                d.set_item("p0", tup(g.p0))?;
                d.set_item("n", tup(g.n))?;
                d.set_item("alpha", g.alpha)?;
                d.set_item("phi0", g.phi0)?;
            }
            core::Segment::Steered(g) => {
                d.set_item("p0", tup(g.p0))?;
                d.set_item("u0", tup(g.u0))?;
                d.set_item("gamma0", g.gamma0)?;
                d.set_item("thrust", tup(g.n))?;
                d.set_item("alpha", g.alpha)?;
            }
            core::Segment::Burn(g) => {
                d.set_item("p0", tup(g.p0))?;
                d.set_item("n", tup(g.n))?;
                d.set_item("phi0", g.phi0)?;
                d.set_item("alpha0", g.sigma * g.w / g.t_burn)?;
                d.set_item("w", g.w)?;
                d.set_item("t_burn", g.t_burn)?;
            }
            core::Segment::Circular(g) => {
                d.set_item("centre", tup(g.centre))?;
                d.set_item("r", g.r)?;
                d.set_item("omega", g.omega)?;
                d.set_item("e1", tup(g.e1))?;
                d.set_item("e2", tup(g.e2))?;
                d.set_item("a0", g.a0)?;
            }
            core::Segment::Null(g) => {
                d.set_item("p0", tup(g.p0))?;
                d.set_item("n", tup(g.n))?;
                d.set_item("p1", tup(g.p1))?;
            }
        }
        Ok(d)
    }

    /// The same segment ending at `t1` instead.
    fn with_t1(&self, t1: f64) -> PyResult<Self> {
        span(self.inner.t0(), t1)?;
        let mut seg = self.inner;
        seg.set_t1(t1);
        check_burn_end(&seg)?;
        Ok(seg.into())
    }

    /// Position and velocity at coordinate time `t`, as `(x, v)`. On a null segment `v` is the unit direction
    /// of travel.
    fn at(&self, t: f64) -> (Tup3, Tup3) {
        let s = self.inner.at(t);
        (tup(s.x), tup(s.v))
    }
    /// Proper time elapsed from `t0` to coordinate time `t`.
    fn dtau(&self, t: f64) -> f64 {
        self.inner.dtau(t)
    }
    /// Coordinate time after proper time `s` from `t0`. On a null segment: its end.
    fn t_of_tau(&self, s: f64) -> f64 {
        self.inner.t_of_tau(s)
    }
    /// The state after proper time `s` from `t0`. On a null segment: the state after the jump, at rest.
    fn at_tau(&self, s: f64) -> TauState {
        self.inner.at_tau(s).into()
    }
    /// Proper time the whole segment takes (infinite if it never ends, zero on a null segment).
    fn proper_duration(&self) -> f64 {
        self.inner.proper_duration()
    }
    /// The largest speed reached between `t0` and `t1`, as a fraction of c. An accelerating segment that never
    /// ends approaches 1.
    fn vmax(&self) -> f64 {
        self.inner.vmax()
    }

    /// The coordinate time, between `lo` and `hi`, at which light from the event `(event_t, event_x)` reaches
    /// this segment.
    ///
    /// The root is unique for anything slower than light: closed form on rest, inertial and null segments, a
    /// bracketed Newton iteration otherwise. If the light has already passed at `lo`, the answer is `lo`.
    /// Raises `srkit.NoIntersection` if the light does not arrive within the range.
    ///
    /// By default the range is all of time, on the segment continued beyond `t0` and `t1`; pass `seg.t0` and
    /// `seg.t1` to keep to the segment's own span. A burn is searched only where it is defined.
    #[pyo3(
        signature = (event_t, event_x, lo=f64::NEG_INFINITY, hi=f64::INFINITY),
        text_signature = "($self, event_t, event_x, lo=-1e999, hi=1e999)"
    )]
    fn recv(&self, event_t: f64, event_x: V3, lo: f64, hi: f64) -> PyResult<f64> {
        let (lo, hi) = self.range(lo, hi);
        self.inner
            .recv(event(event_t, event_x), lo, hi)
            .ok_or_else(|| {
                NoIntersection::new_err(
                    "light from the event does not reach the segment between lo and hi",
                )
            })
    }

    /// The coordinate time on this segment, between `lo` and `hi`, of the moment whose light reaches the event
    /// `(event_t, event_x)`: what an observer at the event sees of this body.
    ///
    /// Raises `srkit.NoIntersection` if that moment is not within the range. The range is as for `recv`: all
    /// of time by default, on the segment continued beyond its span.
    #[pyo3(
        signature = (event_t, event_x, lo=f64::NEG_INFINITY, hi=f64::INFINITY),
        text_signature = "($self, event_t, event_x, lo=-1e999, hi=1e999)"
    )]
    fn send(&self, event_t: f64, event_x: V3, lo: f64, hi: f64) -> PyResult<f64> {
        let (lo, hi) = self.range(lo, hi);
        self.inner
            .send(event(event_t, event_x), lo, hi)
            .ok_or_else(|| {
                NoIntersection::new_err(
                    "no moment on the segment between lo and hi sends light that reaches the event",
                )
            })
    }

    /// For a hyperbolic or steered segment continued for ever: how far the event lies inside the region that
    /// can signal the ship (`past=False`, the future horizon) or that the ship's signals reach (`past=True`).
    /// Positive is inside; zero or negative is on or beyond the horizon. `None` for the other kinds, which
    /// have no horizon.
    #[pyo3(signature = (event_t, event_x, past=false))]
    fn horizon_margin(&self, event_t: f64, event_x: V3, past: bool) -> Option<f64> {
        self.inner.horizon_margin(event(event_t, event_x), past)
    }

    /// The Rindler horizon of a hyperbolic segment continued for ever, or `None` for the other kinds.
    fn rindler_horizon(&self) -> Option<RindlerHorizon> {
        match &self.inner {
            core::Segment::Hyperbolic(g) => Some(RindlerHorizon {
                inner: rindler::RindlerHorizon::of(g),
            }),
            _ => None,
        }
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        let (t0, t1) = (self.inner.t0(), self.inner.t1());
        match self.inner {
            core::Segment::Rest(g) => {
                repr_values!(py, "Segment.rest", t0 = t0, t1 = t1, p = g.p)
            }
            core::Segment::Inertial(g) => {
                repr_values!(py, "Segment.inertial", t0 = t0, t1 = t1, p0 = g.p0, v = g.v)
            }
            core::Segment::Hyperbolic(g) => repr_values!(
                py,
                "Segment.hyperbolic",
                t0 = t0,
                t1 = t1,
                p0 = g.p0,
                n = g.n,
                alpha = g.alpha,
                phi0 = g.phi0
            ),
            // The stored start is the proper velocity, which is what `proper_velocity=True` reads back.
            core::Segment::Steered(g) => repr_values!(
                py,
                "Segment.steered",
                t0 = t0,
                t1 = t1,
                p0 = g.p0,
                v0 = g.u0,
                thrust = g.n,
                alpha = g.alpha,
                proper_velocity = true
            ),
            core::Segment::Burn(g) => repr_values!(
                py,
                "Segment.burn",
                t0 = t0,
                t1 = t1,
                p0 = g.p0,
                n = g.n,
                phi0 = g.phi0,
                alpha0 = g.sigma * g.w / g.t_burn,
                w = g.w
            ),
            core::Segment::Circular(g) => repr_values!(
                py,
                "Segment.circular",
                t0 = t0,
                t1 = t1,
                centre = g.centre,
                r = g.r,
                omega = g.omega,
                e1 = g.e1,
                e2 = g.e2,
                a0 = g.a0
            ),
            core::Segment::Null(g) => {
                repr_values!(py, "Segment.null", t0 = t0, t1 = t1, p0 = g.p0, n = g.n)
            }
        }
    }
}
