//! Relativistic Navigator kernels in geometrized units: Julian years and
//! light-years, with c = 1. This compatibility layer retains rapidity and null
//! coordinates at high gamma; it does not convert a rounded velocity back into
//! rapidity or claim an exact integer causal predicate.
use crate::{
    kinematics::{Vec3, C, G0, JULIAN_YEAR},
    Error,
};

pub use crate::navigator_optics::{boost, observe};

pub const G_LY: f64 = G0 * JULIAN_YEAR / C;

#[derive(Clone, Copy, Debug)]
pub struct Leg {
    pub a: f64,
    pub phc: f64,
    pub ta: f64,
    pub tc: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Trip {
    pub d: f64,
    pub a: f64,
    pub n: u32,
    pub n_out: u32,
    pub rt: bool,
    pub leg: Leg,
    pub leg_tau: f64,
    pub tau: f64,
    pub du: f64,
    pub dv: f64,
    pub u_end: f64,
    pub v_end: f64,
    pub t_end: f64,
    pub t_out: f64,
    pub x_end: f64,
    pub lag: f64,
    pub g_peak: f64,
    pub phc: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct State {
    pub tau: f64,
    pub u: f64,
    pub v: f64,
    pub t: f64,
    pub x: f64,
    /// Nonnegative rapidity along the current direction of flight.
    pub phi: f64,
    pub dir: f64,
    pub gamma: f64,
    pub dh: f64,
    pub dd: f64,
}

fn finite(values: &[f64]) -> Result<(), Error> {
    if values.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(Error::Overflow)
    }
}

impl Trip {
    pub fn new(d: f64, g: f64, cap: f64, stops: u32, roundtrip: bool) -> Result<Self, Error> {
        if !d.is_finite()
            || d < 0.0
            || !g.is_finite()
            || g <= 0.0
            || cap.is_nan()
            || cap < 1.0
            || (d > 0.0 && cap == 1.0)
        {
            return Err(Error::InvalidInput);
        }
        let n_out = stops.checked_add(1).ok_or(Error::Overflow)?;
        let n = if roundtrip {
            n_out.checked_mul(2).ok_or(Error::Overflow)?
        } else {
            n_out
        };
        let a = g * G_LY;
        if !a.is_finite() || a == 0.0 {
            return Err(Error::Overflow);
        }
        let length = d / f64::from(n_out);
        let x = a * length / 2.0;
        let g_pk = 1.0 + x;
        // acosh(1+x) expressed without discarding tiny nonzero distances.
        let (phc, tc) = if cap < g_pk {
            let phi = libm::acosh(cap);
            let coast = (length - 2.0 * (cap - 1.0) / a).max(0.0) / libm::sinh(phi);
            (phi, coast)
        } else {
            (2.0 * libm::asinh(libm::sqrt(x / 2.0)), 0.0)
        };
        let leg = Leg {
            a,
            phc,
            ta: phc / a,
            tc,
        };
        let ep = libm::exp(phc);
        let em = libm::exp(-phc);
        let du = -2.0 * libm::expm1(-phc) / a + em * tc;
        let dv = 2.0 * libm::expm1(phc) / a + ep * tc;
        let leg_tau = 2.0 * leg.ta + tc;
        let tau = f64::from(n) * leg_tau;
        let t_out = f64::from(n_out) * (du + dv) / 2.0;
        let t_end = f64::from(n) * (du + dv) / 2.0;
        let out = Self {
            d,
            a,
            n,
            n_out,
            rt: roundtrip,
            leg,
            leg_tau,
            tau,
            du,
            dv,
            u_end: if roundtrip { t_end } else { f64::from(n) * du },
            v_end: if roundtrip { t_end } else { f64::from(n) * dv },
            t_end,
            t_out,
            x_end: if roundtrip { 0.0 } else { d },
            // This subtraction is a display quantity retained for compatibility;
            // du remains the stable one-way null lag.
            lag: f64::from(n_out) * du,
            g_peak: libm::cosh(phc),
            phc,
        };
        finite(&[
            out.a, out.tau, out.du, out.dv, out.t_end, out.g_peak, out.phc,
        ])?;
        Ok(out)
    }

    /// State at ship time, extending the ship at rest before and after the trip.
    /// The prototype snaps the last 1e-9 fraction of a ride to arrival; retaining
    /// that UI convention makes old slider values land at rest.
    pub fn state(&self, tau: f64) -> Result<State, Error> {
        self.state_impl(tau, true)
    }

    fn state_impl(&self, tau: f64, snap_arrival: bool) -> Result<State, Error> {
        if !tau.is_finite() {
            return Err(Error::InvalidInput);
        }
        if tau < 0.0 {
            return Ok(State {
                tau,
                u: tau,
                v: tau,
                t: tau,
                x: 0.0,
                phi: 0.0,
                dir: 1.0,
                gamma: 1.0,
                dh: 1.0,
                dd: 1.0,
            });
        }
        let tau = if snap_arrival && tau < self.tau && self.tau - tau <= 1e-9 * self.tau {
            self.tau
        } else {
            tau
        };
        if tau >= self.tau {
            let extra = tau - self.tau;
            let state = State {
                tau,
                u: self.u_end + extra,
                v: self.v_end + extra,
                t: self.t_end + extra,
                x: self.x_end,
                phi: 0.0,
                dir: if self.rt { -1.0 } else { 1.0 },
                gamma: 1.0,
                dh: 1.0,
                dd: 1.0,
            };
            finite(&[state.u, state.v, state.t])?;
            return Ok(state);
        }
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
        let ep = libm::exp(phc);
        let em = libm::exp(-phc);
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
        let back = self.rt && i >= self.n_out;
        let (u, v) = if back {
            let j = f64::from(i - self.n_out);
            (
                f64::from(self.n_out) * self.du + j * self.dv + lv,
                f64::from(self.n_out) * self.dv + j * self.du + lu,
            )
        } else {
            (f64::from(i) * self.du + lu, f64::from(i) * self.dv + lv)
        };
        let ep = libm::exp(phi);
        let em = libm::exp(-phi);
        let state = State {
            tau,
            u,
            v,
            t: (u + v) / 2.0,
            x: (v - u) / 2.0,
            phi,
            dir: if back { -1.0 } else { 1.0 },
            gamma: libm::cosh(phi),
            dh: if back { ep } else { em },
            dd: if back { em } else { ep },
        };
        finite(&[
            state.u,
            state.v,
            state.t,
            state.x,
            state.gamma,
            state.dh,
            state.dd,
        ])?;
        Ok(state)
    }

    fn inverse_leg(&self, value: f64, v: bool) -> f64 {
        let Leg { a, phc, ta, tc } = self.leg;
        let ep = libm::exp(phc);
        let em = libm::exp(-phc);
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

    fn inverse_null(&self, value: f64, v: bool) -> Result<f64, Error> {
        if !value.is_finite() {
            return Err(Error::InvalidInput);
        }
        if value <= 0.0 {
            return Ok(value);
        }
        let end = if v { self.v_end } else { self.u_end };
        if value >= end {
            return Ok(self.tau + value - end);
        }
        let first_step = if v { self.dv } else { self.du };
        let n = f64::from(self.n_out);
        let tau = if !self.rt || value < n * first_step {
            let i = libm::floor(value / first_step).min(n - 1.0);
            i * self.leg_tau + self.inverse_leg(value - i * first_step, v)
        } else {
            let rest = value - n * first_step;
            let step = if v { self.du } else { self.dv };
            let j = libm::floor(rest / step).min(n - 1.0);
            (n + j) * self.leg_tau + self.inverse_leg(rest - j * step, !v)
        };
        if tau.is_finite() {
            Ok(tau)
        } else {
            Err(Error::Overflow)
        }
    }

    pub fn tau_u(&self, value: f64) -> Result<f64, Error> {
        self.inverse_null(value, false)
    }
    pub fn tau_v(&self, value: f64) -> Result<f64, Error> {
        self.inverse_null(value, true)
    }

    /// Ship time of a light-cone intersection with an arbitrary event. `along`
    /// is its coordinate along the trip axis and `perp` its perpendicular
    /// distance. `send=false` receives its light; `send=true` sends the ping.
    pub fn intersection(&self, t: f64, along: f64, perp: f64, send: bool) -> Result<f64, Error> {
        if !t.is_finite() || !along.is_finite() || !perp.is_finite() || perp < 0.0 {
            return Err(Error::InvalidInput);
        }
        let line = |tt: f64| -> Result<f64, Error> {
            let a = self.tau_u(tt - along)?;
            let b = self.tau_v(tt + along)?;
            Ok(if send { a.min(b) } else { a.max(b) })
        };
        let nearest = line(t)?;
        if perp <= 1e-13 * (1.0 + along.abs()) {
            return Ok(nearest);
        }
        let (lo, hi) = if send {
            (line(t - perp)?, nearest)
        } else {
            (nearest, line(t + perp)?)
        };
        newton(
            |tau| {
                // Do not use the display endpoint snap inside a root query.
                let s = self.exact_state(tau)?;
                let r = libm::hypot(along - s.x, perp);
                let rapidity = s.dir * s.phi;
                let sign = if send { 1.0 } else { -1.0 };
                let value = s.t + sign * r - t;
                let slope = libm::cosh(rapidity) + sign * (s.x - along) * libm::sinh(rapidity) / r;
                Ok((value, slope))
            },
            lo,
            hi,
            nearest,
            1e-15,
        )
    }

    // Root solving must not inherit a UI slider's near-arrival snapping.
    fn exact_state(&self, tau: f64) -> Result<State, Error> {
        self.state_impl(tau, false)
    }
}

/// Conservative go/flip/brake reach at a heading relative to the initial
/// velocity. The burn direction stays fixed aboard; the final brake is along
/// the resulting velocity. Returns `None` when the budget cannot stop the ship.
pub fn reach(phi: f64, a: f64, budget: f64, heading: f64) -> Result<Option<[f64; 2]>, Error> {
    if !phi.is_finite()
        || phi < 0.0
        || !a.is_finite()
        || a <= 0.0
        || !budget.is_finite()
        || budget < 0.0
        || !heading.is_finite()
    {
        return Err(Error::InvalidInput);
    }
    let (g0, gb0) = (libm::cosh(phi), libm::sinh(phi));
    let (c, sn) = (libm::cos(heading), libm::sin(heading));
    let fd = |t1: f64| {
        let (sh, ch) = (libm::sinh(a * t1), libm::cosh(a * t1));
        let ua = g0 * sh * c + gb0 * ch;
        let us = sh * sn;
        let um = libm::hypot(ua, us);
        let dum = if um > 0.0 {
            a * (ua * (g0 * ch * c + gb0 * sh) + us * ch * sn) / um
        } else {
            0.0
        };
        (
            t1 + libm::asinh(um) / a - budget,
            1.0 + dum / libm::hypot(1.0, um) / a,
        )
    };
    if fd(0.0).0 > 1e-12 {
        return Ok(None);
    }
    let (mut lo, mut hi) = (0.0, budget);
    let mut t1 = ((budget - phi / a) / 2.0).max(0.0).min(budget);
    for _ in 0..80 {
        let (f, d) = fd(t1);
        finite(&[f, d])?;
        if f < 0.0 {
            lo = t1;
        } else {
            hi = t1;
        }
        if f.abs() < 1e-12 * (1.0 + budget) || hi - lo < 1e-13 * (1.0 + budget) {
            break;
        }
        let candidate = if d > 1e-9 { t1 - f / d } else { f64::NAN };
        t1 = if candidate > lo && candidate < hi {
            candidate
        } else {
            (lo + hi) / 2.0
        };
    }
    let (sh, ch) = (libm::sinh(a * t1), libm::cosh(a * t1));
    let k = (g0 - 1.0) * sh * c + gb0 * ch;
    let tm = sh / a;
    let half = libm::sinh(a * t1 / 2.0);
    let hm = 2.0 * half * half / a;
    let (ua, us) = (sh * c + k, sh * sn);
    let um = libm::hypot(ua, us);
    let db = (um / (libm::hypot(1.0, um) + 1.0)) * um / a;
    let out = [
        g0 * hm * c + gb0 * tm + if um > 0.0 { ua * db / um } else { 0.0 },
        hm * sn + if um > 0.0 { us * db / um } else { 0.0 },
    ];
    finite(&out)?;
    Ok(Some(out))
}

/// Emission time from a test traveller that waits at the origin before t=0,
/// then coasts at `beta` along `heading`. Solves the light cone quadratically,
/// without a fixed lookback window or a numerical simulation step.
pub fn traveller_emission(
    beta: f64,
    heading: Vec3,
    observation_t: f64,
    observer: Vec3,
) -> Result<f64, Error> {
    if !beta.is_finite()
        || !(0.0..1.0).contains(&beta)
        || !observation_t.is_finite()
        || !observer.is_finite()
    {
        return Err(Error::InvalidInput);
    }
    let direction = heading.normalized()?;
    let rest_emission = observation_t - observer.norm();
    if !rest_emission.is_finite() {
        return Err(Error::Overflow);
    }
    if rest_emission <= 0.0 {
        return Ok(rest_emission);
    }
    let velocity = direction * beta;
    let w = observer - velocity * observation_t;
    let ww = w.norm_squared();
    let wv = w.dot(velocity);
    let one_minus_vv = (1.0 - beta) * (1.0 + beta);
    let radical = libm::sqrt(wv * wv + one_minus_vv * ww);
    let delay = if ww == 0.0 {
        0.0
    } else if wv >= 0.0 {
        (wv + radical) / one_minus_vv
    } else {
        ww / (radical - wv)
    };
    let emission = observation_t - delay;
    if emission.is_finite() {
        Ok(emission)
    } else {
        Err(Error::Overflow)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Orbit {
    pub r: f64,
    pub a: f64,
    pub gamma: f64,
    pub beta: f64,
    pub phi: f64,
    pub lap_t: f64,
    pub lap_tau: f64,
    pub omega: f64,
    pub sagnac: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct OrbitState {
    pub t: f64,
    pub position: Vec3,
    pub velocity: Vec3,
    pub heading: Vec3,
    pub angle: f64,
}

impl Orbit {
    pub fn new(r: f64, g: f64) -> Result<Self, Error> {
        if !r.is_finite() || r <= 0.0 || !g.is_finite() || g <= 0.0 {
            return Err(Error::InvalidInput);
        }
        let a = g * G_LY;
        let gb = libm::sqrt(a) * libm::sqrt(r);
        let gamma = libm::hypot(1.0, gb);
        let beta = gb / gamma;
        let lap_t = 2.0 * core::f64::consts::PI * r / beta;
        let out = Self {
            r,
            a,
            gamma,
            beta,
            phi: libm::asinh(gb),
            lap_t,
            lap_tau: lap_t / gamma,
            omega: beta / r,
            sagnac: 2.0 * core::f64::consts::PI * r * beta * gamma,
        };
        finite(&[
            a,
            gamma,
            beta,
            out.phi,
            lap_t,
            out.lap_tau,
            out.omega,
            out.sagnac,
        ])?;
        Ok(out)
    }

    pub fn state(&self, a0: f64, tau: f64) -> Result<OrbitState, Error> {
        if !a0.is_finite() || !tau.is_finite() {
            return Err(Error::InvalidInput);
        }
        let t = self.gamma * tau;
        let angle = a0 + self.omega * t;
        finite(&[t, angle])?;
        let (s, c) = (libm::sin(angle), libm::cos(angle));
        let heading = Vec3::new(-s, c, 0.0);
        Ok(OrbitState {
            t,
            position: Vec3::new(self.r * c, self.r * s, 0.0),
            velocity: heading * self.beta,
            heading,
            angle,
        })
    }

    pub fn intersection(&self, a0: f64, t: f64, event: Vec3, send: bool) -> Result<f64, Error> {
        if !a0.is_finite() || !t.is_finite() || !event.is_finite() {
            return Err(Error::InvalidInput);
        }
        let d = event.norm();
        if d == 0.0 {
            return Ok((t + if send { -self.r } else { self.r }) / self.gamma);
        }
        let near = (d - self.r).max(0.0);
        let far = d + self.r;
        let (lo, hi) = if send {
            ((t - far) / self.gamma, (t - near) / self.gamma)
        } else {
            ((t + near) / self.gamma, (t + far) / self.gamma)
        };
        finite(&[lo, hi])?;
        let guess = (t + if send { -d } else { d }) / self.gamma;
        newton(
            |tau| {
                let state = self.state(a0, tau)?;
                let delta = state.position - event;
                let r = delta.norm();
                let sign = if send { 1.0 } else { -1.0 };
                let slope = self.gamma
                    * (1.0
                        + if r == 0.0 {
                            0.0
                        } else {
                            sign * state.velocity.dot(delta) / r
                        });
                Ok((state.t + sign * r - t, slope))
            },
            lo,
            hi,
            guess,
            1e-12,
        )
    }
}

fn newton<F>(mut f: F, mut lo: f64, mut hi: f64, guess: f64, tolerance: f64) -> Result<f64, Error>
where
    F: FnMut(f64) -> Result<(f64, f64), Error>,
{
    if !lo.is_finite() || !hi.is_finite() || lo > hi {
        return Err(Error::Overflow);
    }
    let mut x = guess.max(lo).min(hi);
    for _ in 0..100 {
        let (value, slope) = f(x)?;
        if !value.is_finite() || !slope.is_finite() {
            return Err(Error::Overflow);
        }
        if value == 0.0 {
            return Ok(x);
        }
        if value > 0.0 {
            hi = x;
        } else {
            lo = x;
        }
        let candidate = x - value / slope;
        let next = if candidate > lo && candidate < hi {
            candidate
        } else {
            lo + (hi - lo) * 0.5
        };
        if (next - x).abs() <= tolerance * x.abs().max(1.0)
            || hi - lo <= 2.0 * tolerance * x.abs().max(1.0)
        {
            return Ok(next);
        }
        x = next;
    }
    Err(Error::PrecisionLoss)
}
