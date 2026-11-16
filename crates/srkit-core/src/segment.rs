//! Closed-form worldline segments.
//!
//! A segment gives position, velocity and elapsed proper time as explicit functions of coordinate time `t`, and
//! the inverse (coordinate time after a given proper time). Seven kinds:
//!
//! | Kind | Motion | Proper time |
//! |---|---|---|
//! | [`Rest`] | `x = p` | `dt` |
//! | [`Inertial`] | `x = p0 + v (t - t0)` | `dt / gamma` |
//! | [`Hyperbolic`] | constant proper acceleration along the line of motion | `ln`, exact |
//! | [`Steered`] | constant proper acceleration at any angle, the direction fixed aboard | `ln`, exact |
//! | [`Burn`] | constant thrust while the mass falls (a real rocket) | bracketed root of a closed form |
//! | [`Circular`] | circle at constant speed | `dt / gamma` |
//! | [`Null`] | light-like run (a signal, or a teleport that moves at `c`) | none |
//!
//! Each segment covers `[t0, t1]` in coordinate time; `t1` may be infinite. Evaluation does not clamp: a segment
//! is its own analytic continuation outside its interval, and [`crate::worldline::Worldline`] decides which
//! segment a time belongs to.
//!
//! Natural units: `c = 1`, velocities are fractions of `c`, accelerations are rates of change of rapidity.

use crate::{
    kinematics::Event,
    roots::{bisect, light_meets_line, newton},
    vec3::Vec3,
    Error,
};

/// Position and velocity at a coordinate time.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct State {
    pub x: Vec3,
    /// Velocity as a fraction of `c`. On a [`Null`] segment it is the unit direction of travel.
    pub v: Vec3,
}

/// The state after a given proper time.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TauState {
    pub t: f64,
    pub x: Vec3,
    pub v: Vec3,
    /// Rapidity: signed along the segment's axis for [`Hyperbolic`] and [`Burn`], the speed's rapidity otherwise.
    pub phi: f64,
    /// Proper acceleration (signed along the axis for [`Hyperbolic`]; the magnitude otherwise).
    pub acc: f64,
    /// True for the state after a [`Null`] jump, which takes no proper time.
    pub jump: bool,
}

/// `sinh(p + d) - sinh(p)` without cancellation.
pub fn d_sinh(p: f64, d: f64) -> f64 {
    2.0 * libm::cosh(p + 0.5 * d) * libm::sinh(0.5 * d)
}

/// `cosh(p + d) - cosh(p)` without cancellation.
pub fn d_cosh(p: f64, d: f64) -> f64 {
    2.0 * libm::sinh(p + 0.5 * d) * libm::sinh(0.5 * d)
}

/// At rest at `p`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rest {
    pub t0: f64,
    pub t1: f64,
    pub p: Vec3,
}

/// Constant velocity `v` from `p0` at `t0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Inertial {
    pub t0: f64,
    pub t1: f64,
    pub p0: Vec3,
    pub v: Vec3,
    /// `1 / gamma = sqrt(1 - v²)`.
    pub ig: f64,
}

/// Constant proper acceleration `alpha` (signed) along the unit vector `n`, starting at `p0` with rapidity `phi0`
/// along `n`. After proper time `s` the rapidity is `phi0 + alpha s`; the velocity stays along `n`.
///
/// Rapidity is stored directly, so a ride at gamma 10⁶ loses no precision to a velocity rounded near `c`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hyperbolic {
    pub t0: f64,
    pub t1: f64,
    pub p0: Vec3,
    pub n: Vec3,
    pub alpha: f64,
    pub phi0: f64,
    /// `sinh(phi0)`.
    pub sh0: f64,
    /// `cosh(phi0)`.
    pub ch0: f64,
}

/// Constant proper acceleration `alpha > 0` in a direction fixed aboard, for any angle between thrust and motion.
///
/// The four-velocity is `U(q) = U0 cosh q + N sinh q` with `q = alpha s`, where `U0 = (gamma0, u0)` is the
/// starting four-velocity and `N` the unit thrust four-vector: the thrust direction in the ship's starting rest
/// frame, boosted to this frame. The path is a hyperbola in the plane of `U0` and `N`, so steering in two or
/// three dimensions needs no integration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Steered {
    pub t0: f64,
    pub t1: f64,
    pub p0: Vec3,
    /// Starting proper velocity `gamma v`.
    pub u0: Vec3,
    /// Starting Lorentz factor.
    pub gamma0: f64,
    /// Unit thrust direction in the starting rest frame.
    pub n: Vec3,
    pub alpha: f64,
    /// Time component of the thrust four-vector: `u0·n`.
    pub nt: f64,
    /// Space part of the thrust four-vector.
    pub ns: Vec3,
    /// `gamma0 + nt` and `gamma0 - nt`, each formed without cancellation.
    pub plus: f64,
    pub minus: f64,
}

/// Constant thrust along (or against) the unit vector `n` while the ship's mass falls linearly in its own time:
/// `m(s) = m0 (1 - s / t_burn)`. With exhaust speed `w` (a fraction of `c`; 1 for a photon rocket) the rapidity
/// is `phi0 - sigma w ln(1 - s / t_burn)`, the relativistic rocket equation, and the null coordinates
/// `t ± x` are integrals of powers of `1 - s / t_burn`: a closed form in proper time. The inverse (proper time
/// from coordinate time) is the unique root of that increasing function.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Burn {
    pub t0: f64,
    pub t1: f64,
    pub p0: Vec3,
    pub n: Vec3,
    pub phi0: f64,
    /// Exhaust speed as a fraction of `c`, `0 < w <= 1`.
    pub w: f64,
    /// Proper time in which the whole starting mass would be spent. The segment must end before it.
    pub t_burn: f64,
    /// `+1` thrust along `n`, `-1` against it.
    pub sigma: f64,
}

/// A circle of radius `r` about `centre` in the plane of the orthonormal pair `(e1, e2)`, at angle
/// `a0 + omega (t - t0)` from `e1` towards `e2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circular {
    pub t0: f64,
    pub t1: f64,
    pub centre: Vec3,
    pub r: f64,
    pub omega: f64,
    pub e1: Vec3,
    pub e2: Vec3,
    pub a0: f64,
    /// `1 / gamma = sqrt(1 - r² omega²)`.
    pub ig: f64,
}

/// A light-like run from `p0` along the unit vector `n`, covering `t1 - t0` of distance. No proper time passes.
/// `p1` is the end point; it is stored so that a jump can be pinned exactly on its target.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Null {
    pub t0: f64,
    pub t1: f64,
    pub p0: Vec3,
    pub n: Vec3,
    pub p1: Vec3,
}

/// One closed-form piece of a worldline.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Segment {
    Rest(Rest),
    Inertial(Inertial),
    Hyperbolic(Hyperbolic),
    Steered(Steered),
    Burn(Burn),
    Circular(Circular),
    Null(Null),
}

impl Rest {
    pub const fn new(t0: f64, t1: f64, p: Vec3) -> Self {
        Self { t0, t1, p }
    }
}

impl Inertial {
    pub fn new(t0: f64, t1: f64, p0: Vec3, v: Vec3) -> Self {
        let ig = libm::sqrt((1.0 - v.norm_squared()).max(0.0));
        Self { t0, t1, p0, v, ig }
    }
}

impl Hyperbolic {
    pub fn new(t0: f64, t1: f64, p0: Vec3, n: Vec3, alpha: f64, phi0: f64) -> Self {
        Self {
            t0,
            t1,
            p0,
            n,
            alpha,
            phi0,
            sh0: libm::sinh(phi0),
            ch0: libm::cosh(phi0),
        }
    }
    /// `sinh` and `cosh` of the rapidity at coordinate time `t`, and the change in each since `t0`.
    fn hyp(&self, t: f64) -> (f64, f64, f64, f64) {
        let dsh = self.alpha * (t - self.t0);
        let sh = self.sh0 + dsh;
        // Past 1e150 the square would overflow, and cosh equals |sinh| to every digit.
        let ch = if sh.abs() > 1e150 {
            sh.abs()
        } else {
            libm::sqrt(1.0 + sh * sh)
        };
        let dch = dsh * (sh + self.sh0) / (ch + self.ch0);
        (sh, ch, dsh, dch)
    }
    /// The event the hyperbola is centred on: the corner of the ship's Rindler wedge. Light from events beyond
    /// the plane through it never reaches a ship that keeps accelerating. Returns `(t, x)` of the pivot.
    pub fn pivot(&self) -> Event {
        Event {
            t: self.t0 - self.sh0 / self.alpha,
            x: self.p0 - self.n * (self.ch0 / self.alpha),
        }
    }
}

impl Steered {
    /// A burn from `p0` at `t0`, starting at velocity `v0`, thrusting at proper acceleration `alpha > 0` along
    /// `thrust`, a direction in the ship's starting rest frame.
    pub fn new(
        t0: f64,
        t1: f64,
        p0: Vec3,
        v0: Vec3,
        thrust: Vec3,
        alpha: f64,
    ) -> Result<Self, Error> {
        let gamma0 = crate::kinematics::gamma(v0)?;
        Self::from_proper_velocity(t0, t1, p0, v0 * gamma0, thrust, alpha)
    }
    /// As [`Steered::new`] from the starting proper velocity `u0 = gamma v`, which keeps its precision near `c`.
    pub fn from_proper_velocity(
        t0: f64,
        t1: f64,
        p0: Vec3,
        u0: Vec3,
        thrust: Vec3,
        alpha: f64,
    ) -> Result<Self, Error> {
        if !(alpha.is_finite() && alpha > 0.0) || !u0.is_finite() || !p0.is_finite() {
            return Err(Error::InvalidInput);
        }
        let n = thrust.normalized()?;
        let gamma0 = libm::hypot(1.0, u0.norm());
        let nt = u0.dot(n);
        let ns = n + u0 * (nt / (gamma0 + 1.0));
        // gamma0² - nt² = 1 + |u0 perpendicular to n|², so the smaller of gamma0 ± nt comes from the larger.
        let perp2 = (u0 - n * nt).norm_squared();
        let (plus, minus) = if nt >= 0.0 {
            let plus = gamma0 + nt;
            (plus, (1.0 + perp2) / plus)
        } else {
            let minus = gamma0 - nt;
            ((1.0 + perp2) / minus, minus)
        };
        Ok(Self {
            t0,
            t1,
            p0,
            u0,
            gamma0,
            n,
            alpha,
            nt,
            ns,
            plus,
            minus,
        })
    }
    /// The rapidity parameter `q = alpha s` reached at coordinate time `t`.
    fn q_at(&self, t: f64) -> f64 {
        let big_t = self.alpha * (t - self.t0);
        if big_t == 0.0 {
            return 0.0;
        }
        // plus E² - 2 (T + nt) E - minus = 0 for E = e^q; the root that is 1 at T = 0.
        let tb = big_t + self.nt;
        let r = libm::hypot(tb, libm::sqrt(self.plus * self.minus));
        let r_plus_tb = if tb >= 0.0 {
            r + tb
        } else {
            self.plus * self.minus / (r - tb)
        };
        let e_minus_1 = big_t * (r_plus_tb + self.plus) / ((r + self.gamma0) * self.plus);
        if e_minus_1 < 1.0 {
            libm::log1p(e_minus_1)
        } else {
            libm::log(r_plus_tb / self.plus)
        }
    }
    fn at_q(&self, q: f64) -> (f64, Vec3, Vec3) {
        let (ep, em) = (libm::expm1(q), libm::expm1(-q));
        // gamma0 sinh q + nt (cosh q - 1), as two terms of one sign
        let dt = 0.5 * (self.plus * ep - self.minus * em) / self.alpha;
        let sh = libm::sinh(q);
        let half = libm::sinh(0.5 * q);
        let cm1 = 2.0 * half * half;
        let x = self.p0 + (self.u0 * sh + self.ns * cm1) / self.alpha;
        let gamma = 0.5 * (self.plus * libm::exp(q) + self.minus * libm::exp(-q));
        let u = self.u0 * (1.0 + cm1) + self.ns * sh;
        (dt, x, u / gamma)
    }
}

impl Burn {
    /// A burn from `p0` at `t0` with starting rapidity `phi0` along the unit vector `n`.
    ///
    /// `alpha0` is the proper acceleration at the start (negative to thrust against `n`), `w` the exhaust speed
    /// as a fraction of `c`. The mass would run out after proper time `w / |alpha0|`; the segment must end
    /// before then.
    pub fn new(
        t0: f64,
        t1: f64,
        p0: Vec3,
        n: Vec3,
        phi0: f64,
        alpha0: f64,
        w: f64,
    ) -> Result<Self, Error> {
        if !(alpha0.is_finite() && alpha0 != 0.0) || !(w > 0.0 && w <= 1.0) || !phi0.is_finite() {
            return Err(Error::InvalidInput);
        }
        Ok(Self {
            t0,
            t1,
            p0,
            n: n.normalized()?,
            phi0,
            w,
            t_burn: w / alpha0.abs(),
            sigma: if alpha0 > 0.0 { 1.0 } else { -1.0 },
        })
    }
    /// `[1 - (1 - y)^p] / p` given `l = ln(1 - y)`, continuous through `p = 0` (where it is `-l`).
    fn power_integral(p: f64, l: f64) -> f64 {
        if p == 0.0 {
            -l
        } else {
            -libm::expm1(p * l) / p
        }
    }
    /// Coordinate time and distance along `n` covered in proper time `s`, and the rapidity then.
    fn at_s(&self, s: f64) -> (f64, f64, f64) {
        let l = libm::log1p(-s / self.t_burn);
        let k = self.sigma * self.w;
        let v = libm::exp(self.phi0) * self.t_burn * Self::power_integral(1.0 - k, l);
        let u = libm::exp(-self.phi0) * self.t_burn * Self::power_integral(1.0 + k, l);
        (0.5 * (u + v), 0.5 * (v - u), self.phi0 - k * l)
    }
    /// Proper time elapsed at coordinate time `t`.
    fn s_at(&self, t: f64) -> f64 {
        let dt = t - self.t0;
        if dt <= 0.0 || dt.is_nan() {
            return if dt == 0.0 { 0.0 } else { f64::NAN };
        }
        // dt/ds = cosh(phi) >= 1, so s <= dt; and s stays below the burn-out time.
        let top = self.t_burn * (1.0 - f64::EPSILON);
        let hi = dt.min(top);
        if self.at_s(hi).0 <= dt {
            return hi;
        }
        newton(
            |s| {
                let (t_s, _, phi) = self.at_s(s);
                (t_s - dt, libm::cosh(phi))
            },
            0.0,
            hi,
            dt / libm::cosh(self.phi0),
        )
    }
    /// The fraction of the starting mass left after proper time `s`.
    pub fn mass_fraction(&self, s: f64) -> f64 {
        1.0 - s / self.t_burn
    }
    /// Proper acceleration after proper time `s` (it rises as the ship lightens).
    pub fn acceleration(&self, s: f64) -> f64 {
        self.w / (self.t_burn - s)
    }
}

impl Circular {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        t0: f64,
        t1: f64,
        centre: Vec3,
        r: f64,
        omega: f64,
        e1: Vec3,
        e2: Vec3,
        a0: f64,
    ) -> Self {
        let ig = libm::sqrt((1.0 - r * r * omega * omega).max(0.0));
        Self {
            t0,
            t1,
            centre,
            r,
            omega,
            e1,
            e2,
            a0,
            ig,
        }
    }
}

impl Null {
    pub fn new(t0: f64, t1: f64, p0: Vec3, n: Vec3) -> Self {
        Self {
            t0,
            t1,
            p0,
            n,
            p1: p0 + n * (t1 - t0),
        }
    }
    /// A jump from `p0` at `t0` to the point `p1`, arriving after the light-travel time between them.
    pub fn between(t0: f64, p0: Vec3, p1: Vec3) -> Self {
        let d = p0.dist(p1);
        let n = if d > 0.0 { (p1 - p0) / d } else { Vec3::X };
        Self {
            t0,
            t1: t0 + d,
            p0,
            n,
            p1,
        }
    }
}

impl Segment {
    /// At rest at `p` over `[t0, t1]`.
    pub fn rest(t0: f64, t1: f64, p: Vec3) -> Self {
        Self::Rest(Rest::new(t0, t1, p))
    }
    /// Constant velocity `v` from `p0`.
    pub fn inertial(t0: f64, t1: f64, p0: Vec3, v: Vec3) -> Self {
        Self::Inertial(Inertial::new(t0, t1, p0, v))
    }
    /// Constant proper acceleration `alpha` along the unit vector `n` from rapidity `phi0` along `n`.
    pub fn hyperbolic(t0: f64, t1: f64, p0: Vec3, n: Vec3, alpha: f64, phi0: f64) -> Self {
        Self::Hyperbolic(Hyperbolic::new(t0, t1, p0, n, alpha, phi0))
    }
    /// A circle; see [`Circular`].
    #[allow(clippy::too_many_arguments)]
    pub fn circular(
        t0: f64,
        t1: f64,
        centre: Vec3,
        r: f64,
        omega: f64,
        e1: Vec3,
        e2: Vec3,
        a0: f64,
    ) -> Self {
        Self::Circular(Circular::new(t0, t1, centre, r, omega, e1, e2, a0))
    }
    /// A light-like run from `p0` along the unit vector `n`.
    pub fn null(t0: f64, t1: f64, p0: Vec3, n: Vec3) -> Self {
        Self::Null(Null::new(t0, t1, p0, n))
    }

    /// Start of the segment in coordinate time.
    pub fn t0(&self) -> f64 {
        match self {
            Self::Rest(g) => g.t0,
            Self::Inertial(g) => g.t0,
            Self::Hyperbolic(g) => g.t0,
            Self::Steered(g) => g.t0,
            Self::Burn(g) => g.t0,
            Self::Circular(g) => g.t0,
            Self::Null(g) => g.t0,
        }
    }
    /// End of the segment in coordinate time (may be infinite).
    pub fn t1(&self) -> f64 {
        match self {
            Self::Rest(g) => g.t1,
            Self::Inertial(g) => g.t1,
            Self::Hyperbolic(g) => g.t1,
            Self::Steered(g) => g.t1,
            Self::Burn(g) => g.t1,
            Self::Circular(g) => g.t1,
            Self::Null(g) => g.t1,
        }
    }
    /// Set the end of the segment.
    pub fn set_t1(&mut self, t1: f64) {
        match self {
            Self::Rest(g) => g.t1 = t1,
            Self::Inertial(g) => g.t1 = t1,
            Self::Hyperbolic(g) => g.t1 = t1,
            Self::Steered(g) => g.t1 = t1,
            Self::Burn(g) => g.t1 = t1,
            Self::Circular(g) => g.t1 = t1,
            Self::Null(g) => g.t1 = t1,
        }
    }
    pub fn is_null(&self) -> bool {
        matches!(self, Self::Null(_))
    }

    /// Position and velocity at coordinate time `t`.
    pub fn at(&self, t: f64) -> State {
        match self {
            Self::Rest(g) => State {
                x: g.p,
                v: Vec3::ZERO,
            },
            Self::Inertial(g) => State {
                x: g.p0 + g.v * (t - g.t0),
                v: g.v,
            },
            Self::Hyperbolic(g) => {
                let (sh, ch, dsh, _) = g.hyp(t);
                let x = dsh * (sh + g.sh0) / (ch + g.ch0) / g.alpha;
                State {
                    x: g.p0 + g.n * x,
                    v: g.n * (sh / ch),
                }
            }
            Self::Steered(g) => {
                let (_, x, v) = g.at_q(g.q_at(t));
                State { x, v }
            }
            Self::Burn(g) => {
                let (_, dx, phi) = g.at_s(g.s_at(t));
                State {
                    x: g.p0 + g.n * dx,
                    v: g.n * libm::tanh(phi),
                }
            }
            Self::Circular(g) => {
                let a = g.a0 + g.omega * (t - g.t0);
                let (s, c) = (libm::sin(a), libm::cos(a));
                let rw = g.r * g.omega;
                State {
                    x: g.centre + (g.e1 * c + g.e2 * s) * g.r,
                    v: (g.e2 * c - g.e1 * s) * rw,
                }
            }
            Self::Null(g) => State {
                x: g.p0 + g.n * (t - g.t0),
                v: g.n,
            },
        }
    }

    /// Proper time elapsed from `t0` to coordinate time `t`.
    pub fn dtau(&self, t: f64) -> f64 {
        match self {
            Self::Rest(g) => t - g.t0,
            Self::Inertial(g) => (t - g.t0) * g.ig,
            Self::Circular(g) => (t - g.t0) * g.ig,
            Self::Hyperbolic(g) => {
                let (_, _, dsh, dch) = g.hyp(t);
                libm::log1p((dsh + dch) / (g.sh0 + g.ch0)) / g.alpha
            }
            Self::Steered(g) => g.q_at(t) / g.alpha,
            Self::Burn(g) => g.s_at(t),
            Self::Null(_) => 0.0,
        }
    }

    /// Coordinate time after proper time `s` from `t0`. On a [`Null`] segment: its end.
    pub fn t_of_tau(&self, s: f64) -> f64 {
        match self {
            Self::Rest(g) => g.t0 + s,
            Self::Inertial(g) => g.t0 + s / g.ig,
            Self::Circular(g) => g.t0 + s / g.ig,
            Self::Hyperbolic(g) => g.t0 + d_sinh(g.phi0, g.alpha * s) / g.alpha,
            Self::Steered(g) => g.t0 + g.at_q(g.alpha * s).0,
            Self::Burn(g) => g.t0 + g.at_s(s).0,
            Self::Null(g) => g.t1,
        }
    }

    /// The state after proper time `s` from `t0`. On a [`Null`] segment: the state after the jump, at rest.
    pub fn at_tau(&self, s: f64) -> TauState {
        match self {
            Self::Rest(g) => TauState {
                t: g.t0 + s,
                x: g.p,
                ..TauState::default()
            },
            Self::Inertial(g) => {
                let dt = s / g.ig;
                TauState {
                    t: g.t0 + dt,
                    x: g.p0 + g.v * dt,
                    v: g.v,
                    phi: libm::atanh(g.v.norm().min(1.0 - 1e-16)),
                    acc: 0.0,
                    jump: false,
                }
            }
            Self::Hyperbolic(g) => {
                let d = g.alpha * s;
                let phi = g.phi0 + d;
                TauState {
                    t: g.t0 + d_sinh(g.phi0, d) / g.alpha,
                    x: g.p0 + g.n * (d_cosh(g.phi0, d) / g.alpha),
                    v: g.n * libm::tanh(phi),
                    phi,
                    acc: g.alpha,
                    jump: false,
                }
            }
            Self::Steered(g) => {
                let (dt, x, v) = g.at_q(g.alpha * s);
                TauState {
                    t: g.t0 + dt,
                    x,
                    v,
                    phi: libm::atanh(v.norm().min(1.0 - 1e-16)),
                    acc: g.alpha,
                    jump: false,
                }
            }
            Self::Burn(g) => {
                let (dt, dx, phi) = g.at_s(s);
                TauState {
                    t: g.t0 + dt,
                    x: g.p0 + g.n * dx,
                    v: g.n * libm::tanh(phi),
                    phi,
                    acc: g.sigma * g.acceleration(s),
                    jump: false,
                }
            }
            Self::Circular(g) => {
                let t = g.t0 + s / g.ig;
                let st = self.at(t);
                let b = g.r * g.omega.abs();
                TauState {
                    t,
                    x: st.x,
                    v: st.v,
                    phi: libm::atanh(b),
                    acc: b * b / (g.ig * g.ig) / g.r,
                    jump: false,
                }
            }
            Self::Null(g) => TauState {
                t: g.t1,
                x: g.p1,
                jump: true,
                ..TauState::default()
            },
        }
    }

    /// Proper time the whole segment takes (infinite if it never ends, zero for a [`Null`] segment).
    pub fn proper_duration(&self) -> f64 {
        self.dtau(self.t1())
    }

    /// The largest speed reached on `[t0, t1]`, as a fraction of `c`. An accelerating segment that never ends
    /// approaches 1.
    pub fn vmax(&self) -> f64 {
        match self {
            Self::Rest(_) => 0.0,
            Self::Inertial(g) => g.v.norm(),
            Self::Circular(g) => (g.r * g.omega).abs(),
            Self::Null(_) => 1.0,
            Self::Hyperbolic(g) => {
                // The speed has no interior maximum (it can only dip through zero), so it peaks at an end.
                let start = libm::tanh(g.phi0).abs();
                if g.t1.is_finite() {
                    start.max(self.at(g.t1).v.norm())
                } else {
                    1.0
                }
            }
            Self::Steered(_) | Self::Burn(_) => {
                let start = self.at(self.t0()).v.norm();
                if self.t1().is_finite() {
                    start.max(self.at(self.t1()).v.norm())
                } else {
                    1.0
                }
            }
        }
    }

    /// For a segment of constant proper acceleration continued for ever, how far `event` lies inside the region
    /// that can signal the ship (`past = false`: the future horizon) or that the ship's signals reach
    /// (`past = true`: the past horizon). Positive is inside; zero or negative is on or beyond the horizon.
    /// `None` for the other kinds, which have no horizon.
    ///
    /// The worldline is `x_c + (U0 sinh q + N cosh q) / alpha`. It approaches the null plane through `x_c` with
    /// normal `k = U0 + N` in the future, and the one with normal `U0 - N` in the past.
    pub fn horizon_margin(&self, event: Event, past: bool) -> Option<f64> {
        match self {
            Self::Hyperbolic(g) => {
                let pivot = g.pivot();
                let sign = if g.alpha >= 0.0 { 1.0 } else { -1.0 };
                let along = (event.x - pivot.x).dot(g.n) * sign;
                let dt = event.t - pivot.t;
                Some(if past { along + dt } else { along - dt })
            }
            Self::Steered(g) => {
                let dt = event.t - (g.t0 - g.nt / g.alpha);
                let dx = event.x - (g.p0 - g.ns / g.alpha);
                Some(if past {
                    // (U0 - N)·(E - x_c) > 0
                    g.minus * dt - (g.u0 - g.ns).dot(dx)
                } else {
                    // (U0 + N)·(E - x_c) < 0
                    (g.u0 + g.ns).dot(dx) - g.plus * dt
                })
            }
            _ => None,
        }
    }

    /// When light from `event` reaches this segment, searching `t` in `[lo, hi]`.
    ///
    /// `f(t) = t - event.t - |x(t) - event.x|` increases along any worldline slower than light, so the root is
    /// unique: closed form for [`Rest`], [`Inertial`] and [`Null`], a bracketed Newton iteration otherwise.
    /// `None` if the light does not arrive within the range. If the event's light has already passed at `lo`,
    /// the answer is `lo`.
    pub fn recv(&self, event: Event, lo: f64, hi: f64) -> Option<f64> {
        let lo = lo.max(event.t);
        // `!(hi >= lo)` also rejects a NaN bound.
        if !(hi >= lo) {
            return None;
        }
        let f = |t: f64| {
            let s = self.at(t);
            let d = s.x - event.x;
            let r = d.norm();
            (
                t - event.t - r,
                1.0 - if r > 0.0 { s.v.dot(d) / r } else { 0.0 },
            )
        };
        if f(lo).0 >= 0.0 {
            return Some(lo);
        }
        if hi.is_finite() && f(hi).0 < 0.0 {
            return None;
        }
        match self {
            Self::Rest(g) => Some(event.t + g.p.dist(event.x)),
            Self::Inertial(g) => {
                let w0 = self.at(event.t).x - event.x;
                Some(event.t + light_meets_line(w0.norm_squared(), w0.dot(g.v), g.v.norm_squared()))
            }
            Self::Null(g) => null_root(self, g, event, lo, hi, false),
            _ => {
                let mut b = hi;
                if !b.is_finite() {
                    // A ship that accelerates for ever outruns light from behind its horizon. Decide that in
                    // closed form: the search below would otherwise chase a root that rounding invents.
                    if self.horizon_margin(event, false).is_some_and(|m| m <= 0.0) {
                        return None;
                    }
                    let r0 = self.at(lo).x.dist(event.x);
                    let mut step = r0 * 2.0 + 1.0;
                    b = lo + step;
                    while f(b).0 < 0.0 {
                        step *= 2.0;
                        b = lo + step;
                        if !b.is_finite() {
                            return None;
                        }
                    }
                }
                Some(newton(f, lo, b, lo + self.at(lo).x.dist(event.x)))
            }
        }
    }

    /// The time on this segment, in `[lo, hi]`, of the event whose light reaches `event`: what an observer at
    /// `event` sees of this body. `None` if that moment is not within the range.
    ///
    /// `f(t) = t + |event.x - x(t)| - event.t` increases along any worldline slower than light.
    pub fn send(&self, event: Event, lo: f64, hi: f64) -> Option<f64> {
        let hi = hi.min(event.t);
        if !(hi >= lo) {
            return None;
        }
        let f = |t: f64| {
            let s = self.at(t);
            let d = event.x - s.x;
            let r = d.norm();
            (
                t + r - event.t,
                1.0 - if r > 0.0 { s.v.dot(d) / r } else { 0.0 },
            )
        };
        let at_hi = f(hi).0;
        if at_hi <= 0.0 {
            return if at_hi == 0.0 { Some(hi) } else { None };
        }
        if lo.is_finite() && f(lo).0 > 0.0 {
            return None;
        }
        match self {
            Self::Rest(g) => Some(event.t - event.x.dist(g.p)),
            Self::Inertial(g) => {
                let w0 = event.x - self.at(event.t).x;
                // Looking back along the line: the body recedes from the event at -v in reversed time.
                Some(event.t - light_meets_line(w0.norm_squared(), w0.dot(g.v), g.v.norm_squared()))
            }
            Self::Null(g) => null_root(self, g, event, lo, hi, true),
            _ => {
                let mut a = lo;
                if !a.is_finite() {
                    if self.horizon_margin(event, true).is_some_and(|m| m <= 0.0) {
                        return None;
                    }
                    let r1 = event.x.dist(self.at(hi).x);
                    let mut step = r1 * 2.0 + 1.0;
                    a = hi - step;
                    while f(a).0 > 0.0 {
                        step *= 2.0;
                        a = hi - step;
                        if !a.is_finite() {
                            return None;
                        }
                    }
                }
                Some(newton(f, a, hi, hi - event.x.dist(self.at(hi).x)))
            }
        }
    }
}

/// Light-cone roots on a light-like segment.
///
/// `f` is non-decreasing there and can be flat (light running alongside the jump), so `recv` takes the first `t`
/// with `f >= 0` and `send` the last with `f <= 0`. With `w0 = x(tE) - xE` and a unit `n`: for `recv`,
/// `|w0 + n u| = u` gives `u = -|w0|² / (2 w0·n)` when `w0·n < 0`; for `send`, `u = |w0|² / (2 w0·n)` when
/// `w0·n > 0`. Bisection covers the degenerate cases.
fn null_root(seg: &Segment, g: &Null, event: Event, lo: f64, hi: f64, send: bool) -> Option<f64> {
    let w0 = seg.at(event.t).x - event.x;
    let ww = w0.norm_squared();
    let wv = w0.dot(g.n);
    let mut t = if ww == 0.0 {
        event.t
    } else if send {
        if wv > 0.0 {
            event.t - ww / (2.0 * wv)
        } else {
            f64::NAN
        }
    } else if wv < 0.0 {
        event.t + ww / (-2.0 * wv)
    } else {
        f64::NAN
    };
    if !(t >= lo && t <= hi) || !t.is_finite() {
        if !lo.is_finite() || !hi.is_finite() {
            return None;
        }
        t = if send {
            bisect(|t| t + event.x.dist(seg.at(t).x) - event.t, lo, hi, false)
        } else {
            bisect(|t| t - event.t - seg.at(t).x.dist(event.x), lo, hi, true)
        };
    }
    Some(t)
}
