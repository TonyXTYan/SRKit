//! Events, four-vectors, boosts, velocity composition, Doppler factors and aberration.
//!
//! Natural units (`c = 1`). The metric signature is `(+ - - -)`: a timelike interval is positive.

use core::ops::{Add, Mul, Sub};

use crate::{vec3::Vec3, Error};

/// A point in spacetime: coordinate time `t` and position `x` in one inertial frame.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Event {
    pub t: f64,
    pub x: Vec3,
}

/// How two events are related, decided in floating point. For an exact decision on stored integer times use
/// [`crate::causal`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Separation {
    /// One event can influence the other with something slower than light.
    Timelike,
    /// Light from one just reaches the other.
    Null,
    /// Neither can influence the other.
    Spacelike,
}

impl Event {
    pub const fn new(t: f64, x: Vec3) -> Self {
        Self { t, x }
    }
    pub fn is_finite(self) -> bool {
        self.t.is_finite() && self.x.is_finite()
    }
    /// The invariant interval squared to `other`: `dt² - |dx|²`.
    pub fn interval2(self, other: Self) -> f64 {
        let dt = other.t - self.t;
        let dx = other.x.dist(self.x);
        (dt - dx) * (dt + dx)
    }
    /// Classify the separation to `other`.
    pub fn separation(self, other: Self) -> Separation {
        let s2 = self.interval2(other);
        if s2 > 0.0 {
            Separation::Timelike
        } else if s2 < 0.0 {
            Separation::Spacelike
        } else {
            Separation::Null
        }
    }
    /// Whether light (or anything slower) leaving this event can reach `other`.
    pub fn can_reach(self, other: Self) -> bool {
        other.t - self.t >= other.x.dist(self.x)
    }
    /// The displacement four-vector from this event to `other`.
    pub fn to(self, other: Self) -> FourVector {
        FourVector {
            t: other.t - self.t,
            x: other.x - self.x,
        }
    }
}

/// A four-vector `(t, x)` with both parts in the same units.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FourVector {
    pub t: f64,
    pub x: Vec3,
}

impl FourVector {
    pub const fn new(t: f64, x: Vec3) -> Self {
        Self { t, x }
    }
    /// The four-velocity `(gamma, gamma v)` of a body moving at velocity `v`.
    pub fn velocity(v: Vec3) -> Result<Self, Error> {
        let g = gamma(v)?;
        Ok(Self { t: g, x: v * g })
    }
    /// The Minkowski product, signature `(+ - - -)`.
    pub fn dot(self, other: Self) -> f64 {
        self.t * other.t - self.x.dot(other.x)
    }
    /// The same vector in a frame moving at velocity `beta` relative to this one.
    ///
    /// Uses `gamma² / (gamma + 1)` for `(gamma - 1) / beta²`, which does not cancel at low speed. The inverse is
    /// `boost(-beta)`.
    pub fn boost(self, beta: Vec3) -> Result<Self, Error> {
        if !self.t.is_finite() || !self.x.is_finite() {
            return Err(Error::InvalidInput);
        }
        let g = gamma(beta)?;
        let projection = beta.dot(self.x);
        let coefficient = g / (1.0 + 1.0 / g);
        let out = Self {
            t: g * (self.t - projection),
            x: self.x + beta * (coefficient * projection - g * self.t),
        };
        if out.t.is_finite() && out.x.is_finite() {
            Ok(out)
        } else {
            Err(Error::Overflow)
        }
    }
    /// The same vector in a frame whose velocity points along `rapidity` with that rapidity (see
    /// [`boost_rapidity`]).
    pub fn boost_rapidity(self, rapidity: Vec3) -> Self {
        let (t, x) = boost_rapidity(rapidity, self.t, self.x);
        Self { t, x }
    }
}

impl Add for FourVector {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.t + rhs.t, self.x + rhs.x)
    }
}
impl Sub for FourVector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.t - rhs.t, self.x - rhs.x)
    }
}
impl Mul<f64> for FourVector {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.t * rhs, self.x * rhs)
    }
}

/// Lorentz factor of a velocity (a fraction of `c`). A speed at or above `c` is an error.
pub fn gamma(beta: Vec3) -> Result<f64, Error> {
    if !beta.is_finite() {
        return Err(Error::InvalidInput);
    }
    gamma_of_speed(beta.norm())
}

/// Lorentz factor of a speed `0 <= b < 1`, as `1 / sqrt((1 - b)(1 + b))` so that it stays accurate near `c`.
pub fn gamma_of_speed(b: f64) -> Result<f64, Error> {
    if !(0.0..1.0).contains(&b) {
        return Err(Error::InvalidInput);
    }
    Ok(1.0 / libm::sqrt((1.0 - b) * (1.0 + b)))
}

/// Rapidity of a speed: `atanh(b)`.
pub fn rapidity_of_speed(b: f64) -> f64 {
    libm::atanh(b)
}

/// Speed of a rapidity: `tanh(phi)`.
pub fn speed_of_rapidity(phi: f64) -> f64 {
    libm::tanh(phi)
}

/// `1 - beta` for a rapidity, without the cancellation of `1 - tanh(phi)`: `2 / (e^(2 phi) + 1)`.
pub fn speed_shortfall(phi: f64) -> f64 {
    2.0 / (libm::exp(2.0 * phi) + 1.0)
}

/// Kinetic energy per unit rest energy, `gamma - 1`, from a rapidity, without cancellation: `2 sinh²(phi / 2)`.
pub fn kinetic_per_rest_energy(phi: f64) -> f64 {
    let s = libm::sinh(0.5 * phi);
    2.0 * s * s
}

/// Boost a displacement `(dt, dx)` into the frame whose velocity points along `rapidity` and whose rapidity is
/// that vector's length.
///
/// Works on the null components `dt ± dx∥`, each scaled by `e^(∓phi)`, so light running along the boost keeps
/// full precision at any gamma. Total: non-finite input gives non-finite output.
pub fn boost_rapidity(rapidity: Vec3, dt: f64, dx: Vec3) -> (f64, Vec3) {
    let phi = rapidity.norm();
    if phi == 0.0 || phi.is_nan() {
        return (dt, dx);
    }
    let axis = rapidity / phi;
    let parallel = axis.dot(dx);
    // Halve before adding so that the sum cannot overflow where the result does not.
    let minus = (0.5 * dt - 0.5 * parallel) * libm::exp(phi);
    let plus = (0.5 * dt + 0.5 * parallel) * libm::exp(-phi);
    (minus + plus, (dx - axis * parallel) + axis * (plus - minus))
}

/// Relativistic composition of velocities: the velocity, in the lab, of a body that moves at `u` in a frame
/// which itself moves at `v` through the lab.
pub fn compose_velocities(v: Vec3, u: Vec3) -> Result<Vec3, Error> {
    let g = gamma(v)?;
    gamma(u)?;
    let vu = v.dot(u);
    let out = (v * (1.0 + g * vu / (1.0 + g)) + u / g) / (1.0 + vu);
    if out.is_finite() {
        Ok(out)
    } else {
        Err(Error::Overflow)
    }
}

/// The velocity of a body moving at `u` in the lab, as measured in a frame moving at `v` through the lab.
pub fn relative_velocity(v: Vec3, u: Vec3) -> Result<Vec3, Error> {
    let four = FourVector::velocity(u)?.boost(v)?;
    Ok(four.x / four.t)
}

/// Received frequency divided by emitted frequency.
///
/// Both velocities are in one inertial frame. `photon_direction` points along the light's travel, from the
/// emitter to the receiver, in that frame, and need not be a unit vector. The factor is
/// `gamma_r (1 - beta_r·n) / [gamma_e (1 - beta_e·n)]`.
pub fn doppler_factor(
    emitter_velocity: Vec3,
    receiver_velocity: Vec3,
    photon_direction: Vec3,
) -> Result<f64, Error> {
    let n = photon_direction.normalized()?;
    let ge = gamma(emitter_velocity)?;
    let gr = gamma(receiver_velocity)?;
    let out = gr * (1.0 - receiver_velocity.dot(n)) / (ge * (1.0 - emitter_velocity.dot(n)));
    if out.is_finite() && out > 0.0 {
        Ok(out)
    } else {
        Err(Error::OutOfRange)
    }
}

/// Head-on Doppler factor for a rapidity: `e^phi` approaching, `e^(-phi)` receding (pass a negative rapidity).
pub fn doppler_of_rapidity(phi: f64) -> f64 {
    libm::exp(phi)
}

/// The direction a photon travels, as measured in a frame moving at velocity `beta`.
///
/// The input is the photon's direction of travel in the original frame. The direction on the sky towards the
/// source is the negative of both input and output.
pub fn aberrate(direction: Vec3, beta: Vec3) -> Result<Vec3, Error> {
    let photon = FourVector::new(1.0, direction.normalized()?).boost(beta)?;
    (photon.x / photon.t).normalized()
}

/// The angle from the direction of motion at which a source appears to an observer of rapidity `phi`, given the
/// angle `theta` at which it lies in the rest frame: `tan(theta'/2) = e^(-phi) tan(theta/2)`.
pub fn aberrate_angle(theta: f64, phi: f64) -> f64 {
    2.0 * libm::atan(libm::exp(-phi) * libm::tan(0.5 * theta))
}
