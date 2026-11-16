//! `srkit.CircularOrbit`: a circle flown under thrust.

use pyo3::prelude::*;
use srkit_core::{orbit as core, Vec3};

use crate::{
    convert::{
        event, finite, finite_vec, fl, invalid, positive, repr_fields, tup, unit_vec, Tup3, V3,
    },
    segment::Segment,
};

/// A point on a `CircularOrbit`: what `CircularOrbit.state` returns.
#[pyclass(frozen, get_all, module = "srkit", name = "OrbitState")]
pub struct OrbitState {
    /// Coordinate time.
    t: f64,
    /// Position `(x, y, z)`.
    x: Tup3,
    /// Velocity as a fraction of c.
    v: Tup3,
    /// Unit vector along the motion.
    heading: Tup3,
    /// Angle round the circle, in radians.
    angle: f64,
}

#[pymethods]
impl OrbitState {
    fn __repr__(&self, py: Python<'_>) -> String {
        repr_fields!(py, "OrbitState", self, t, x, v, heading, angle)
    }
}

/// A circle of radius `r` flown at constant speed: its period, the clocks' rates and what the pilot feels.
///
/// `CircularOrbit(r, beta)` takes the radius and the speed as a fraction of c;
/// `CircularOrbit.from_radius_acceleration(r, acc)` takes the proper acceleration held aboard instead.
/// Natural units, c = 1.
#[pyclass(frozen, eq, module = "srkit", name = "CircularOrbit")]
#[derive(Clone, Copy, PartialEq)]
pub struct CircularOrbit {
    pub(crate) inner: core::CircularOrbit,
}

impl From<core::CircularOrbit> for CircularOrbit {
    fn from(inner: core::CircularOrbit) -> Self {
        Self { inner }
    }
}

/// Check a radius and a speed, and build the orbit.
pub(crate) fn orbit_of_speed(r: f64, beta: f64) -> PyResult<core::CircularOrbit> {
    positive("r", r)?;
    if !(beta > 0.0 && beta < 1.0) {
        return Err(invalid(format!(
            "beta must be above 0 and below 1 (a fraction of c), got {}",
            fl(beta)
        )));
    }
    Ok(core::CircularOrbit::from_radius_speed(r, beta))
}

#[pymethods]
impl CircularOrbit {
    #[new]
    fn new(r: f64, beta: f64) -> PyResult<Self> {
        Ok(orbit_of_speed(r, beta)?.into())
    }

    /// From the radius and the speed `beta` (a fraction of c, above 0 and below 1).
    #[staticmethod]
    fn from_radius_speed(r: f64, beta: f64) -> PyResult<Self> {
        Ok(orbit_of_speed(r, beta)?.into())
    }

    /// From the radius and the proper acceleration held aboard: `gamma beta = sqrt(acc r)`.
    #[staticmethod]
    fn from_radius_acceleration(r: f64, acc: f64) -> PyResult<Self> {
        positive("r", r)?;
        positive("acc", acc)?;
        Ok(core::CircularOrbit::from_radius_acceleration(r, acc).into())
    }

    /// Radius.
    #[getter]
    fn r(&self) -> f64 {
        self.inner.r
    }
    /// Speed as a fraction of c.
    #[getter]
    fn beta(&self) -> f64 {
        self.inner.beta
    }
    /// Lorentz factor.
    #[getter]
    fn gamma(&self) -> f64 {
        self.inner.gamma
    }
    /// Rapidity of the speed.
    #[getter]
    fn phi(&self) -> f64 {
        self.inner.phi
    }
    /// Angular velocity in coordinate time, `beta / r`, in radians per unit of time.
    #[getter]
    fn omega(&self) -> f64 {
        self.inner.omega
    }
    /// Coordinate time for one lap.
    #[getter]
    fn period(&self) -> f64 {
        self.inner.period
    }
    /// Proper time for one lap.
    #[getter]
    fn lap_tau(&self) -> f64 {
        self.inner.lap_tau
    }
    /// Proper acceleration felt aboard: `gamma^2 beta^2 / r`.
    #[getter]
    fn acc(&self) -> f64 {
        self.inner.acc
    }

    /// The Sagnac gap per lap, `2 pi r beta gamma`: how far clocks synchronised round the loop disagree where
    /// the loop closes.
    fn sagnac(&self) -> f64 {
        self.inner.sagnac()
    }
    /// Thomas precession per lap, in radians: a gyroscope carried round turns by `2 pi (gamma - 1)` against
    /// the stars.
    fn thomas_precession_per_lap(&self) -> f64 {
        self.inner.thomas_precession_per_lap()
    }
    /// The smallest and largest Doppler factors a distant source in the orbit's plane shows over a lap, as
    /// `(gamma (1 - beta), gamma (1 + beta))`.
    fn doppler_range(&self) -> (f64, f64) {
        self.inner.doppler_range()
    }
    /// The state after proper time `tau` on a circle about the origin in the xy plane, starting at angle `a0`
    /// (radians) at `t = 0`.
    fn state(&self, a0: f64, tau: f64) -> OrbitState {
        let s = self.inner.state(a0, tau);
        OrbitState {
            t: s.t,
            x: tup(s.x),
            v: tup(s.v),
            heading: tup(s.heading),
            angle: s.angle,
        }
    }
    /// The orbit as a `Segment` about `centre` in the plane of the orthonormal pair `(e1, e2)` (by default the
    /// x and y axes), starting at angle `a0` at time `t0`.
    #[pyo3(signature = (t0, t1, centre, e1=None, e2=None, a0=0.0))]
    fn segment(
        &self,
        t0: f64,
        t1: f64,
        centre: V3,
        e1: Option<V3>,
        e2: Option<V3>,
        a0: f64,
    ) -> PyResult<Segment> {
        if t0.is_nan() || t1.is_nan() || t1 < t0 {
            return Err(invalid("t0 and t1 must be numbers with t1 not before t0"));
        }
        let e1 = unit_vec("e1", e1.map_or(Vec3::X, |v| v.0))?;
        let e2 = unit_vec("e2", e2.map_or(Vec3::Y, |v| v.0))?;
        if e1.dot(e2).abs() > 1e-9 {
            return Err(invalid("e1 and e2 must be at right angles"));
        }
        Ok(self
            .inner
            .segment(
                t0,
                t1,
                finite_vec("centre", centre.0)?,
                e1,
                e2,
                finite("a0", a0)?,
            )
            .into())
    }
    /// The orbiting ship's proper time when it receives light from the event (`send=False`), or when it emits
    /// the light that reaches the event (`send=True`). The ship circles the origin in the xy plane and is at
    /// angle `a0` at `t = 0`, where its proper time is zero.
    #[pyo3(signature = (a0, event_t, event_x, send=false))]
    fn light_cone_tau(&self, a0: f64, event_t: f64, event_x: V3, send: bool) -> f64 {
        self.inner.light_cone_tau(a0, event(event_t, event_x), send)
    }

    fn __repr__(&self) -> String {
        format!(
            "CircularOrbit(r={}, beta={})",
            fl(self.inner.r),
            fl(self.inner.beta)
        )
    }
}
