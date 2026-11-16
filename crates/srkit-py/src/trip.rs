//! `srkit.Trip`, `srkit.FlipBurn` and `srkit.reach`: rest-to-rest rides in null coordinates.

use pyo3::prelude::*;
use srkit_core::{reach as core_reach, trip as core};

use crate::convert::{fl, invalid, positive, repr_fields, IntoPy};

/// One leg of a `Trip`: boost for proper time `ta` to rapidity `phc`, coast for `tc`, brake for `ta`.
#[pyclass(frozen, get_all, module = "srkit", name = "Leg")]
pub struct Leg {
    /// Proper acceleration.
    a: f64,
    /// Peak (cruise) rapidity.
    phc: f64,
    /// Proper time of each burn.
    ta: f64,
    /// Proper time coasting.
    tc: f64,
}

#[pymethods]
impl Leg {
    fn __repr__(&self, py: Python<'_>) -> String {
        repr_fields!(py, "Leg", self, a, phc, ta, tc)
    }
}

/// The ship's state at a proper time on a `Trip`: what `Trip.state` returns.
#[pyclass(frozen, get_all, module = "srkit", name = "TripState")]
pub struct TripState {
    /// The ship's proper time.
    tau: f64,
    /// Null coordinate `u = t - x`: the launch site's clock as the ship sees it (outbound).
    u: f64,
    /// Null coordinate `v = t + x`.
    v: f64,
    /// Coordinate time.
    t: f64,
    /// Position along the line of flight.
    x: f64,
    /// Rapidity along the present direction of flight (never negative).
    phi: f64,
    /// `+1.0` outbound, `-1.0` on the way home.
    dir: f64,
    /// Lorentz factor.
    gamma: f64,
    /// Rate at which the launch site's clock is seen to run: `exp(-phi)` receding, `exp(phi)` approaching.
    doppler_home: f64,
    /// Rate at which the destination's clock is seen to run.
    doppler_dest: f64,
}

impl From<core::TripState> for TripState {
    fn from(s: core::TripState) -> Self {
        Self {
            tau: s.tau,
            u: s.u,
            v: s.v,
            t: s.t,
            x: s.x,
            phi: s.phi,
            dir: s.dir,
            gamma: s.gamma,
            doppler_home: s.doppler_home,
            doppler_dest: s.doppler_dest,
        }
    }
}

#[pymethods]
impl TripState {
    fn __repr__(&self, py: Python<'_>) -> String {
        repr_fields!(
            py,
            "TripState",
            self,
            tau,
            u,
            v,
            t,
            x,
            phi,
            dir,
            gamma,
            doppler_home,
            doppler_dest
        )
    }
}

/// A ride of identical rest-to-rest legs out to distance `d`, and optionally the same legs home again.
///
/// `Trip(d, a, cap=inf, stops=0, round_trip=False)`: distance `d`, proper acceleration `a`, the Lorentz factor
/// capped at `cap` (the ship coasts once it is reached), `stops` evenly spaced stopovers on the way, and a
/// return if `round_trip`.
///
/// Along the line of flight `u = t - x` and `v = t + x` are the null coordinates; each grows by a pure
/// exponential of the rapidity, so the whole ride stays exact to the last digit at any Lorentz factor.
/// Natural units, c = 1: with `d` in light-years, `a` is per year (1 g is `srkit.units.G0_PER_YEAR`) and times
/// are in years.
#[pyclass(frozen, eq, module = "srkit", name = "Trip")]
#[derive(Clone, Copy, PartialEq)]
pub struct Trip {
    inner: core::Trip,
    cap: f64,
}

#[pymethods]
impl Trip {
    #[new]
    #[pyo3(
        signature = (d, a, cap=f64::INFINITY, stops=0, round_trip=false),
        text_signature = "(d, a, cap=1e999, stops=0, round_trip=False)"
    )]
    fn new(d: f64, a: f64, cap: f64, stops: i64, round_trip: bool) -> PyResult<Self> {
        if stops < 0 {
            return Err(invalid(format!("stops must not be negative, got {stops}")));
        }
        let stops = u32::try_from(stops).map_err(|_| {
            pyo3::exceptions::PyOverflowError::new_err("stops must fit an unsigned 32-bit count")
        })?;
        Ok(Self {
            inner: core::Trip::new(d, a, cap, stops, round_trip).py()?,
            cap,
        })
    }

    /// One-way distance.
    #[getter]
    fn d(&self) -> f64 {
        self.inner.d
    }
    /// Proper acceleration.
    #[getter]
    fn a(&self) -> f64 {
        self.inner.a
    }
    /// The cap on the Lorentz factor (infinite for none).
    #[getter]
    fn cap(&self) -> f64 {
        self.cap
    }
    /// Number of stopovers on the way out.
    #[getter]
    fn stops(&self) -> u32 {
        self.inner.n_out - 1
    }
    /// Number of legs in all.
    #[getter]
    fn n(&self) -> u32 {
        self.inner.n
    }
    /// Number of legs on the way out.
    #[getter]
    fn n_out(&self) -> u32 {
        self.inner.n_out
    }
    /// Whether the ride returns to its start.
    #[getter]
    fn round_trip(&self) -> bool {
        self.inner.round_trip
    }
    /// The burn and coast times of one leg, as a `Leg`.
    #[getter]
    fn leg(&self) -> Leg {
        let l = self.inner.leg;
        Leg {
            a: l.a,
            phc: l.phc,
            ta: l.ta,
            tc: l.tc,
        }
    }
    /// Proper time per leg.
    #[getter]
    fn leg_tau(&self) -> f64 {
        self.inner.leg_tau
    }
    /// Proper time for the whole ride.
    #[getter]
    fn tau(&self) -> f64 {
        self.inner.tau
    }
    /// Increase of `u = t - x` over one outbound leg.
    #[getter]
    fn du(&self) -> f64 {
        self.inner.du
    }
    /// Increase of `v = t + x` over one outbound leg.
    #[getter]
    fn dv(&self) -> f64 {
        self.inner.dv
    }
    /// `u = t - x` at the end of the ride.
    #[getter]
    fn u_end(&self) -> f64 {
        self.inner.u_end
    }
    /// `v = t + x` at the end of the ride.
    #[getter]
    fn v_end(&self) -> f64 {
        self.inner.v_end
    }
    /// Coordinate time at the end of the ride.
    #[getter]
    fn t_end(&self) -> f64 {
        self.inner.t_end
    }
    /// Coordinate time at the far end.
    #[getter]
    fn t_out(&self) -> f64 {
        self.inner.t_out
    }
    /// Position at the end of the ride (`d`, or 0 for a round trip).
    #[getter]
    fn x_end(&self) -> f64 {
        self.inner.x_end
    }
    /// How far the ship arrives behind a light signal that left with it: `t_out - d`. It is also the launch
    /// site's clock reading that the ship sees on arrival.
    #[getter]
    fn lag(&self) -> f64 {
        self.inner.lag
    }
    /// Peak Lorentz factor.
    #[getter]
    fn gamma_peak(&self) -> f64 {
        self.inner.gamma_peak
    }
    /// Peak (cruise) rapidity.
    #[getter]
    fn phi_peak(&self) -> f64 {
        self.inner.leg.phc
    }

    /// The state at proper time `tau`, with the ship waiting at rest before launch (`tau < 0`) and after
    /// arrival. A `tau` within one part in a thousand million of the ride's length counts as arrival.
    fn state(&self, tau: f64) -> TripState {
        self.inner.state(tau).into()
    }
    /// The state at proper time `tau` clamped to the ride: before launch it is the launch, after arrival the
    /// arrival.
    fn state_clamped(&self, tau: f64) -> TripState {
        self.inner.state_clamped(tau).into()
    }
    /// The ship's proper time when its `u = t - x` equals `u`: when it sees the launch site's clock read `u`
    /// (outbound). Closed form.
    fn tau_of_u(&self, u: f64) -> f64 {
        self.inner.tau_of_u(u)
    }
    /// The ship's proper time when its `v = t + x` equals `v`. Closed form.
    fn tau_of_v(&self, v: f64) -> f64 {
        self.inner.tau_of_v(v)
    }
    /// The ship's proper time at which a light cone of the event `(t, along, across)` meets its worldline.
    ///
    /// `along` is the event's coordinate on the line of flight and `across >= 0` its distance from that line.
    /// With `send=False` it is when the ship receives the event's light; with `send=True`, when the ship must
    /// emit for its light to reach the event.
    #[pyo3(signature = (t, along, across=0.0, send=false))]
    fn light_cone_tau(&self, t: f64, along: f64, across: f64, send: bool) -> PyResult<f64> {
        self.inner.light_cone_tau(t, along, across, send).py()
    }
    /// Mass ratio (start over end) of the whole ride for a rocket of exhaust speed `w`, a fraction of c
    /// (1 for a photon rocket).
    #[pyo3(signature = (w=1.0))]
    fn mass_ratio(&self, w: f64) -> PyResult<f64> {
        exhaust(w)?;
        Ok(self.inner.mass_ratio(w))
    }

    fn __repr__(&self) -> String {
        format!(
            "Trip(d={}, a={}, cap={}, stops={}, round_trip={})",
            fl(self.inner.d),
            fl(self.inner.a),
            fl(self.cap),
            self.inner.n_out - 1,
            if self.inner.round_trip {
                "True"
            } else {
                "False"
            }
        )
    }
}

/// Check an exhaust speed.
pub(crate) fn exhaust(w: f64) -> PyResult<f64> {
    if w > 0.0 && w <= 1.0 {
        Ok(w)
    } else {
        Err(invalid(format!(
            "w must be above 0 and at most 1 (a fraction of c), got {}",
            fl(w)
        )))
    }
}

/// The symmetric rest-to-rest flip-and-burn over distance `d` at proper acceleration `a`: accelerate to the
/// midpoint, turn, brake.
///
/// `FlipBurn(d, a)`. The forms used stay accurate from a hop between moons to the far side of the observable
/// universe. Natural units, c = 1.
#[pyclass(frozen, eq, module = "srkit", name = "FlipBurn")]
#[derive(Clone, Copy, PartialEq)]
pub struct FlipBurn {
    inner: core::FlipBurn,
}

#[pymethods]
impl FlipBurn {
    #[new]
    fn new(d: f64, a: f64) -> PyResult<Self> {
        Ok(Self {
            inner: core::FlipBurn::new(d, a).py()?,
        })
    }
    /// Distance.
    #[getter]
    fn d(&self) -> f64 {
        self.inner.d
    }
    /// Proper acceleration.
    #[getter]
    fn a(&self) -> f64 {
        self.inner.a
    }
    /// Coordinate time: `sqrt(d^2 + 4 d / a)`.
    #[getter]
    fn t(&self) -> f64 {
        self.inner.t
    }
    /// Proper time: `(2 / a) acosh(1 + a d / 2)`.
    #[getter]
    fn tau(&self) -> f64 {
        self.inner.tau
    }
    /// Peak rapidity.
    #[getter]
    fn phi_peak(&self) -> f64 {
        self.inner.phi_peak
    }
    /// Peak Lorentz factor: `1 + a d / 2`.
    #[getter]
    fn gamma_peak(&self) -> f64 {
        self.inner.gamma_peak
    }
    /// `t - d`: how far behind light the ship arrives, formed without subtracting two nearly equal times.
    #[getter]
    fn lag(&self) -> f64 {
        self.inner.lag
    }
    /// Mass ratio (start over end) for a rocket of exhaust speed `w`, a fraction of c (1 for a photon
    /// rocket): `exp(2 phi_peak / w)`.
    #[pyo3(signature = (w=1.0))]
    fn mass_ratio(&self, w: f64) -> PyResult<f64> {
        exhaust(w)?;
        Ok(self.inner.mass_ratio(w))
    }
    fn __repr__(&self) -> String {
        format!("FlipBurn(d={}, a={})", fl(self.inner.d), fl(self.inner.a))
    }
}

/// Where a ship can come to rest, given how much longer it may fly.
///
/// The ship moves at rapidity `phi` and can hold proper acceleration `a` for a further proper time `budget`.
/// It burns in a direction fixed aboard at angle `heading` (radians) from its velocity, then turns and brakes
/// to rest. Returns the stopping point as `(along, across)` relative to the ship's position and line of
/// flight, or `None` when the budget cannot even cancel the ship's present speed.
#[pyfunction]
pub(crate) fn reach(phi: f64, a: f64, budget: f64, heading: f64) -> PyResult<Option<(f64, f64)>> {
    Ok(core_reach::reach(phi, a, budget, heading)
        .py()?
        .map(|p| (p[0], p[1])))
}

/// Proper time for a rest-to-rest flip-and-burn over distance `d` at proper acceleration `a`, with no coast:
/// `(2 / a) acosh(1 + a d / 2)`.
#[pyfunction]
pub(crate) fn flip_burn_tau(d: f64, a: f64) -> PyResult<f64> {
    positive("a", a)?;
    Ok(core_reach::flip_burn_tau(d, a))
}

/// The distance a rest-to-rest flip-and-burn covers in proper time `tau` at proper acceleration `a`:
/// `(2 / a)(cosh(a tau / 2) - 1)`.
#[pyfunction]
pub(crate) fn flip_burn_distance(tau: f64, a: f64) -> PyResult<f64> {
    positive("a", a)?;
    Ok(core_reach::flip_burn_distance(tau, a))
}
