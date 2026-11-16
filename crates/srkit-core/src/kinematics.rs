//! Three-dimensional Minkowski kinematics in SI units, with signature `(+---)`.
//!
//! Boost velocities are dimensionless beta vectors; ordinary velocities are m/s.
//! These local floating-point constructions do not decide integer-time causality.

use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::Error;

/// Speed of light in vacuum, exactly defined in m/s.
pub const C: f64 = 299_792_458.0;
/// Julian year, exactly 365.25 days, in seconds.
pub const JULIAN_YEAR: f64 = 31_557_600.0;
/// Standard gravity in m/s² (not a local gravitational-field model).
pub const G0: f64 = 9.80665;

/// Cartesian vector. Components have the units of the quantity being represented.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);

    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean norm, using scaled hypot to avoid intermediate square overflow.
    pub fn norm(self) -> f64 {
        libm::hypot(libm::hypot(self.x, self.y), self.z)
    }

    pub fn normalized(self) -> Result<Self, Error> {
        if !self.is_finite() {
            return Err(Error::InvalidInput);
        }
        // Normalize a scaled vector, not the original norm. A subnormal norm
        // can round by a large fraction (e.g. hypot(MIN_SUBNORMAL, ...));
        // a finite vector's norm may also exceed f64::MAX.
        let scale = self.x.abs().max(self.y.abs()).max(self.z.abs());
        if scale == 0.0 {
            return Err(Error::InvalidInput);
        }
        let scaled = self / scale;
        Ok(scaled / scaled.norm())
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

/// Local event: coordinate time in seconds and position in metres.
///
/// Use exact timestamps and an explicit local origin for deep-time storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Event {
    pub t: f64,
    pub position: Vec3,
}

impl Event {
    pub fn new(t: f64, position: Vec3) -> Result<Self, Error> {
        if !t.is_finite() || !position.is_finite() {
            return Err(Error::InvalidInput);
        }
        Ok(Self { t, position })
    }
}

/// A four-vector with the same units in temporal and spatial components.
///
/// For a displacement use `(c * dt, dx)` in metres; for a dimensionless
/// four-velocity use `(gamma, gamma * v / c)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FourVector {
    pub temporal: f64,
    pub spatial: Vec3,
}

impl FourVector {
    pub fn new(temporal: f64, spatial: Vec3) -> Result<Self, Error> {
        if !temporal.is_finite() || !spatial.is_finite() {
            return Err(Error::InvalidInput);
        }
        Ok(Self { temporal, spatial })
    }

    /// Minkowski product with signature `(+---)`.
    pub fn dot(self, other: Self) -> f64 {
        self.temporal * other.temporal - self.spatial.dot(other.spatial)
    }

    /// Coordinates in a frame moving at `beta * c` relative to the input frame.
    ///
    /// Uses `gamma² / (gamma + 1)` for `(gamma - 1) / beta²`, avoiding
    /// cancellation at small beta. The inverse boost is `boost(-beta)`.
    pub fn boost(self, beta: Vec3) -> Result<Self, Error> {
        if !self.temporal.is_finite() || !self.spatial.is_finite() {
            return Err(Error::InvalidInput);
        }
        let gamma = gamma_from_beta(beta)?;
        let projection = beta.dot(self.spatial);
        let coefficient = gamma / (1.0 + 1.0 / gamma);
        let result = Self {
            temporal: gamma * (self.temporal - projection),
            spatial: self.spatial + beta * (coefficient * projection - gamma * self.temporal),
        };
        if !result.temporal.is_finite() || !result.spatial.is_finite() {
            return Err(Error::OutOfRange);
        }
        Ok(result)
    }
}

impl Add for FourVector {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            temporal: self.temporal + rhs.temporal,
            spatial: self.spatial + rhs.spatial,
        }
    }
}

impl Sub for FourVector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            temporal: self.temporal - rhs.temporal,
            spatial: self.spatial - rhs.spatial,
        }
    }
}

impl Mul<f64> for FourVector {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self {
            temporal: self.temporal * rhs,
            spatial: self.spatial * rhs,
        }
    }
}

impl Div<f64> for FourVector {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self {
            temporal: self.temporal / rhs,
            spatial: self.spatial / rhs,
        }
    }
}

fn gamma_from_beta(beta: Vec3) -> Result<f64, Error> {
    if !beta.is_finite() {
        return Err(Error::InvalidInput);
    }
    let speed = beta.norm();
    if speed >= 1.0 {
        return Err(Error::InvalidInput);
    }
    Ok(1.0 / libm::sqrt((1.0 - speed) * (1.0 + speed)))
}

/// Lorentz factor for a subluminal velocity in m/s.
pub fn lorentz_factor(velocity: Vec3) -> Result<f64, Error> {
    gamma_from_beta(velocity / C)
}

/// Received frequency divided by emitted rest-frame frequency.
///
/// Both velocities are in the same lab frame (m/s). `photon_direction` points
/// along propagation **from emitter to receiver** in that lab frame. It may
/// have any nonzero finite length; it is normalized internally. Formula:
/// `gamma_r (1 - beta_r · n) / [gamma_e (1 - beta_e · n)]`.
/// This is a local frequency transformation; it does not establish visibility.
pub fn doppler_factor(
    emitter_velocity: Vec3,
    receiver_velocity: Vec3,
    photon_direction: Vec3,
) -> Result<f64, Error> {
    let direction = photon_direction.normalized()?;
    let emitter_beta = emitter_velocity / C;
    let receiver_beta = receiver_velocity / C;
    let emitter_gamma = gamma_from_beta(emitter_beta)?;
    let receiver_gamma = gamma_from_beta(receiver_beta)?;
    let result = receiver_gamma * (1.0 - receiver_beta.dot(direction))
        / (emitter_gamma * (1.0 - emitter_beta.dot(direction)));
    if !result.is_finite() || result <= 0.0 {
        return Err(Error::OutOfRange);
    }
    Ok(result)
}

/// Photon propagation direction seen in a frame moving at dimensionless `beta`.
///
/// The input is the propagation direction, rather than the direction on the
/// observer's sky toward the source (which is its negative).
pub fn aberration(direction: Vec3, beta: Vec3) -> Result<Vec3, Error> {
    let photon = FourVector::new(1.0, direction.normalized()?)?.boost(beta)?;
    (photon.spatial / photon.temporal).normalized()
}
