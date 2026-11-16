//! Rest-to-rest rides in null coordinates: flip-and-burn legs with an optional speed cap, stopovers and a return.
//!
//! Along the line of flight `u = t - x` and `v = t + x` are the null coordinates: `u` is the launch site's clock
//! as the ship sees it, `v - d` the destination's. Each grows by a pure exponential of the rapidity
//! (`du = e^(-phi) dtau`, `dv = e^(phi) dtau`), so a whole ride is sums of `expm1` terms and stays exact to the
//! last digit at any gamma. Rapidity, not velocity, is the state.
//!
//! Natural units: `c = 1`; `a` is the proper acceleration as a rate of change of rapidity.

use crate::{roots::newton, Error};

/// One leg: boost for proper time `ta` to rapidity `phc`, coast for `tc`, brake for `ta`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Leg {
    pub a: f64,
    /// Peak (cruise) rapidity.
    pub phc: f64,
    /// Proper time of each burn.
    pub ta: f64,
    /// Proper time coasting.
    pub tc: f64,
}

impl Leg {
    /// A rest-to-rest leg of length `length` at proper acceleration `a`, the Lorentz factor capped at `cap`
    /// (`f64::INFINITY` for none).
    pub fn new(length: f64, a: f64, cap: f64) -> Self {
        let x = a * length / 2.0; // peak gamma - 1 with no coast
        if cap < 1.0 + x {
            let phc = libm::acosh(cap);
            let tc = (length - 2.0 * (cap - 1.0) / a).max(0.0) / libm::sinh(phc);
            Self {
                a,
                phc,
                ta: phc / a,
                tc,
            }
        } else {
            let phc = crate::plan::acosh1p(x);
            Self {
                a,
                phc,
                ta: phc / a,
                tc: 0.0,
            }
        }
    }
}

/// A ride of identical rest-to-rest legs out to distance `d`, and optionally the same legs home again.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Trip {
    /// One-way distance.
    pub d: f64,
    /// Proper acceleration.
    pub a: f64,
    /// Number of legs in all.
    pub n: u32,
    /// Number of legs on the way out.
    pub n_out: u32,
    /// Whether the ride returns to its start.
    pub round_trip: bool,
    pub leg: Leg,
    /// Proper time per leg.
    pub leg_tau: f64,
    /// Proper time for the whole ride.
    pub tau: f64,
    /// Increase of `u = t - x` over one outbound leg.
    pub du: f64,
    /// Increase of `v = t + x` over one outbound leg.
    pub dv: f64,
    pub u_end: f64,
    pub v_end: f64,
    /// Coordinate time at the end of the ride.
    pub t_end: f64,
    /// Coordinate time at the far end.
    pub t_out: f64,
    /// Position at the end of the ride (`d`, or 0 for a round trip).
    pub x_end: f64,
    /// How far the ship arrives behind a light signal that left with it: `t_out - d`. This is also the launch
    /// site's clock reading that the ship sees on arrival.
    pub lag: f64,
    /// Peak Lorentz factor.
    pub gamma_peak: f64,
}

/// The ship's state at a proper time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TripState {
    pub tau: f64,
    pub u: f64,
    pub v: f64,
    pub t: f64,
    pub x: f64,
    /// Rapidity along the present direction of flight (never negative).
    pub phi: f64,
    /// `+1` outbound, `-1` on the way home.
    pub dir: f64,
    pub gamma: f64,
    /// Rate at which the launch site's clock is seen to run: `e^(-phi)` receding, `e^(phi)` approaching.
    pub doppler_home: f64,
    /// Rate at which the destination's clock is seen to run.
    pub doppler_dest: f64,
}

impl Trip {
    /// Plan a ride over distance `d` at proper acceleration `a`, with the Lorentz factor capped at `cap`
    /// (`f64::INFINITY` for none), `stops` evenly spaced stopovers on the way, and a return if `round_trip`.
    pub fn new(d: f64, a: f64, cap: f64, stops: u32, round_trip: bool) -> Result<Self, Error> {
        if !(d.is_finite() && d >= 0.0)
            || !(a.is_finite() && a > 0.0)
            || cap.is_nan()
            || cap < 1.0
            || (d > 0.0 && cap == 1.0)
        {
            return Err(Error::InvalidInput);
        }
        let n_out = stops.checked_add(1).ok_or(Error::Overflow)?;
        let n = if round_trip {
            n_out.checked_mul(2).ok_or(Error::Overflow)?
        } else {
            n_out
        };
        let leg = Leg::new(d / f64::from(n_out), a, cap);
        let Leg { phc, tc, .. } = leg;
        let (ep, em) = (libm::exp(phc), libm::exp(-phc));
        let du = -2.0 * libm::expm1(-phc) / a + em * tc;
        let dv = 2.0 * libm::expm1(phc) / a + ep * tc;
        let leg_tau = 2.0 * leg.ta + tc;
        let t_out = f64::from(n_out) * (du + dv) / 2.0;
        let t_end = f64::from(n) * (du + dv) / 2.0;
        let out = Self {
            d,
            a,
            n,
            n_out,
            round_trip,
            leg,
            leg_tau,
            tau: f64::from(n) * leg_tau,
            du,
            dv,
            u_end: if round_trip { t_end } else { f64::from(n) * du },
            v_end: if round_trip { t_end } else { f64::from(n) * dv },
            t_end,
            t_out,
            x_end: if round_trip { 0.0 } else { d },
            lag: f64::from(n_out) * du,
            gamma_peak: libm::cosh(phc),
        };
        if [out.tau, out.du, out.dv, out.t_end, out.gamma_peak]
            .iter()
            .all(|x| x.is_finite())
        {
            Ok(out)
        } else {
            Err(Error::Overflow)
        }
    }

    /// The state at proper time `tau`, with the ship waiting at rest before launch (`tau < 0`) and after arrival.
    ///
    /// A `tau` within 1e-9 of the ride's length counts as arrival, so that a slider a hair short of its end
    /// shows the ship at rest.
    pub fn state(&self, tau: f64) -> TripState {
        self.state_impl(tau, true)
    }

    /// The state at proper time `tau` clamped to the ride: before launch it is the launch, after arrival the
    /// arrival.
    pub fn state_clamped(&self, tau: f64) -> TripState {
        let mut s = self.state_impl(tau.max(0.0).min(self.tau), true);
        if tau >= self.tau {
            // the ride's last state keeps the direction it arrived in
            s.dir = if self.round_trip { -1.0 } else { 1.0 };
        }
        s
    }

    fn rest(&self, tau: f64, u: f64, v: f64, x: f64, dir: f64) -> TripState {
        TripState {
            tau,
            u,
            v,
            t: 0.5 * (u + v),
            x,
            phi: 0.0,
            dir,
            gamma: 1.0,
            doppler_home: 1.0,
            doppler_dest: 1.0,
        }
    }

    fn state_impl(&self, tau: f64, snap: bool) -> TripState {
        if tau < 0.0 {
            return self.rest(tau, tau, tau, 0.0, 1.0);
        }
        let tau = if snap && tau < self.tau && self.tau - tau <= 1e-9 * self.tau {
            self.tau
        } else {
            tau
        };
        if tau >= self.tau {
            let extra = tau - self.tau;
            let mut s = self.rest(
                tau,
                self.u_end + extra,
                self.v_end + extra,
                self.x_end,
                if self.round_trip { -1.0 } else { 1.0 },
            );
            s.t = self.t_end + extra;
            return s;
        }
        // Which leg: a tau that is a whole number of legs to within rounding starts the next leg.
        let quotient = tau / self.leg_tau;
        let boundary = libm::round(quotient);
        let quotient =
            if (quotient - boundary).abs() <= 4.0 * f64::EPSILON * quotient.abs().max(1.0) {
                boundary
            } else {
                quotient
            };
        let i = (libm::floor(quotient) as u32).min(self.n - 1);
        let s = (tau - f64::from(i) * self.leg_tau).max(0.0);
        let Leg { a, phc, ta, tc } = self.leg;
        let (ep, em) = (libm::exp(phc), libm::exp(-phc));
        // this leg's own null-coordinate increments, along its direction of travel
        let (lu, lv, phi) = if s <= ta {
            let phi = a * s;
            (-libm::expm1(-phi) / a, libm::expm1(phi) / a, phi)
        } else if s <= ta + tc {
            (
                -libm::expm1(-phc) / a + em * (s - ta),
                libm::expm1(phc) / a + ep * (s - ta),
                phc,
            )
        } else {
            let s2 = (s - ta - tc).min(ta);
            let phi = (phc - a * s2).max(0.0);
            (
                -libm::expm1(-phc) / a + em * tc + em * libm::expm1(phc - phi) / a,
                libm::expm1(phc) / a + ep * tc - ep * libm::expm1(phi - phc) / a,
                phi,
            )
        };
        // Heading home the roles swap: the launch site's light now comes from ahead.
        let back = self.round_trip && i >= self.n_out;
        let (u, v) = if back {
            let j = f64::from(i - self.n_out);
            (
                f64::from(self.n_out) * self.du + j * self.dv + lv,
                f64::from(self.n_out) * self.dv + j * self.du + lu,
            )
        } else {
            (f64::from(i) * self.du + lu, f64::from(i) * self.dv + lv)
        };
        let (ed, emd) = (libm::exp(phi), libm::exp(-phi));
        TripState {
            tau,
            u,
            v,
            t: (u + v) / 2.0,
            x: (v - u) / 2.0,
            phi,
            dir: if back { -1.0 } else { 1.0 },
            gamma: libm::cosh(phi),
            doppler_home: if back { ed } else { emd },
            doppler_dest: if back { emd } else { ed },
        }
    }

    fn inverse_leg(&self, value: f64, v: bool) -> f64 {
        let Leg { a, phc, ta, tc } = self.leg;
        let (ep, em) = (libm::exp(phc), libm::exp(-phc));
        if !v {
            let first = -libm::expm1(-phc) / a;
            let second = first + em * tc;
            if value <= first {
                -libm::log1p(-(a * value.max(0.0)).min(1.0 - em)) / a
            } else if value <= second {
                ta + (value - first) / em
            } else {
                ta + tc + (phc + libm::log((em + a * (value - second)).min(1.0))) / a
            }
        } else {
            let first = libm::expm1(phc) / a;
            let second = first + ep * tc;
            if value <= first {
                libm::log1p(a * value.max(0.0)) / a
            } else if value <= second {
                ta + (value - first) / ep
            } else {
                ta + tc + (phc - libm::log((ep - a * (value - second)).max(1.0))) / a
            }
        }
    }

    fn inverse_null(&self, value: f64, v: bool) -> f64 {
        if value.is_nan() || value <= 0.0 {
            return value;
        }
        let end = if v { self.v_end } else { self.u_end };
        if value >= end {
            return self.tau + value - end;
        }
        let first_step = if v { self.dv } else { self.du };
        let n = f64::from(self.n_out);
        if !self.round_trip || value < n * first_step {
            let i = libm::floor(value / first_step).min(n - 1.0);
            i * self.leg_tau + self.inverse_leg(value - i * first_step, v)
        } else {
            let rest = value - n * first_step;
            let step = if v { self.du } else { self.dv };
            let j = libm::floor(rest / step).min(n - 1.0);
            (n + j) * self.leg_tau + self.inverse_leg(rest - j * step, !v)
        }
    }

    /// The ship's proper time when its `u = t - x` equals `u`: when it sees the launch site's clock read `u`
    /// (outbound), in closed form.
    pub fn tau_of_u(&self, u: f64) -> f64 {
        self.inverse_null(u, false)
    }
    /// The ship's proper time when its `v = t + x` equals `v`.
    pub fn tau_of_v(&self, v: f64) -> f64 {
        self.inverse_null(v, true)
    }

    /// The ship's proper time at which a light cone of the event `(t, along, across)` meets its worldline.
    ///
    /// `along` is the event's coordinate on the line of flight and `across >= 0` its distance from that line.
    /// With `send = false` it is when the ship receives the event's light; with `send = true`, when the ship
    /// must emit for its light to reach the event. On the line of flight this is the closed-form null
    /// inverse; off it, a bracketed Newton iteration between the on-axis answers.
    pub fn light_cone_tau(
        &self,
        t: f64,
        along: f64,
        across: f64,
        send: bool,
    ) -> Result<f64, Error> {
        if !t.is_finite() || !along.is_finite() || !(across.is_finite() && across >= 0.0) {
            return Err(Error::InvalidInput);
        }
        let line = |tt: f64| {
            let a = self.tau_of_u(tt - along);
            let b = self.tau_of_v(tt + along);
            if send {
                a.min(b)
            } else {
                a.max(b)
            }
        };
        let nearest = line(t);
        if across <= 1e-13 * (1.0 + along.abs()) {
            return Ok(nearest);
        }
        let (lo, hi) = if send {
            (line(t - across), nearest)
        } else {
            (nearest, line(t + across))
        };
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return Err(Error::Overflow);
        }
        let mut failed = false;
        let root = newton(
            |tau| {
                // no arrival snap inside a root query
                let s = self.state_impl(tau, false);
                let r = libm::hypot(along - s.x, across);
                let rapidity = s.dir * s.phi;
                let sign = if send { 1.0 } else { -1.0 };
                let value = s.t + sign * r - t;
                let slope = libm::cosh(rapidity) + sign * (s.x - along) * libm::sinh(rapidity) / r;
                if !value.is_finite() || !slope.is_finite() {
                    failed = true;
                }
                (value, slope)
            },
            lo,
            hi,
            nearest,
        );
        if failed || !root.is_finite() {
            Err(Error::Overflow)
        } else {
            Ok(root)
        }
    }

    /// Mass ratio of the whole ride for a rocket of exhaust speed `w`: every leg boosts and brakes through
    /// `phc`, so the rapidity budget is `2 n phc`.
    pub fn mass_ratio(&self, w: f64) -> f64 {
        crate::rocket::mass_ratio(2.0 * f64::from(self.n) * self.leg.phc, w)
    }
}

/// The symmetric rest-to-rest flip-and-burn over distance `d` at proper acceleration `a`, in the forms that stay
/// accurate from a hop between moons to the far side of the observable universe.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlipBurn {
    pub d: f64,
    pub a: f64,
    /// Coordinate time: `sqrt(d² + 4 d / a)`.
    pub t: f64,
    /// Proper time: `(2 / a) acosh(1 + a d / 2)`.
    pub tau: f64,
    /// Peak rapidity.
    pub phi_peak: f64,
    /// Peak Lorentz factor: `1 + a d / 2`.
    pub gamma_peak: f64,
    /// `t - d`: how far behind light the ship arrives, formed without subtracting two nearly equal times.
    pub lag: f64,
}

impl FlipBurn {
    pub fn new(d: f64, a: f64) -> Result<Self, Error> {
        if !(d.is_finite() && d >= 0.0) || !(a.is_finite() && a > 0.0) {
            return Err(Error::InvalidInput);
        }
        if d == 0.0 {
            return Ok(Self {
                d,
                a,
                t: 0.0,
                tau: 0.0,
                phi_peak: 0.0,
                gamma_peak: 1.0,
                lag: 0.0,
            });
        }
        // t² = d² + 4 d / a. Square roots first, so that neither a d nor 1 / a can overflow or underflow.
        let (root_d, root_a) = (libm::sqrt(d), libm::sqrt(a));
        let newtonian = 2.0 * root_d / root_a;
        let t = libm::hypot(d, newtonian);
        let half_sinh = 0.5 * root_a * root_d; // sinh(phi_peak / 2)
        let asinh = libm::asinh(half_sinh);
        let tau = newtonian
            * if half_sinh == 0.0 {
                1.0
            } else {
                (asinh / half_sinh).min(1.0)
            };
        let lag = newtonian * ((newtonian / t) / (1.0 + d / t));
        let out = Self {
            d,
            a,
            t,
            tau,
            phi_peak: 2.0 * asinh,
            gamma_peak: 1.0 + 2.0 * half_sinh * half_sinh,
            lag,
        };
        if out.t.is_finite() && out.tau.is_finite() && out.gamma_peak.is_finite() {
            Ok(out)
        } else {
            Err(Error::Overflow)
        }
    }
    /// Mass ratio for a rocket of exhaust speed `w` (1 for a photon rocket): `exp(2 phi_peak / w)`.
    pub fn mass_ratio(&self, w: f64) -> f64 {
        crate::rocket::mass_ratio(2.0 * self.phi_peak, w)
    }
}
