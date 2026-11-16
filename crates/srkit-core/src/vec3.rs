//! Cartesian three-vectors. Components carry the units of the quantity they stand for.

use core::ops::{Add, AddAssign, Div, Index, Mul, Neg, Sub, SubAssign};

use crate::Error;

/// A Cartesian vector `(x, y, z)`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);
    pub const X: Self = Self::new(1.0, 0.0, 0.0);
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub const fn from_array(a: [f64; 3]) -> Self {
        Self::new(a[0], a[1], a[2])
    }
    pub const fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }
    /// Euclidean length, by scaled `hypot` so that the squares cannot overflow or underflow.
    pub fn norm(self) -> f64 {
        libm::hypot(libm::hypot(self.x, self.y), self.z)
    }
    /// Distance to another point.
    pub fn dist(self, other: Self) -> f64 {
        (self - other).norm()
    }
    /// The unit vector along this one. A zero or non-finite vector is an error.
    pub fn normalized(self) -> Result<Self, Error> {
        if !self.is_finite() {
            return Err(Error::InvalidInput);
        }
        // Scale first: the norm of a subnormal vector rounds badly and that of a huge one can overflow.
        let scale = self.x.abs().max(self.y.abs()).max(self.z.abs());
        if scale == 0.0 {
            return Err(Error::InvalidInput);
        }
        let scaled = self / scale;
        Ok(scaled / scaled.norm())
    }
    /// The unit vector along this one, or `fallback` for a zero or non-finite vector.
    pub fn normalized_or(self, fallback: Self) -> Self {
        self.normalized().unwrap_or(fallback)
    }
    /// A unit vector at right angles to `self` (assumed a unit vector), the one nearest the `z` axis; for a
    /// vector along `z`, the one nearest `x`.
    pub fn perpendicular(self) -> Self {
        let a = if self.z.abs() < 0.9 { Self::Z } else { Self::X };
        self.cross(a).cross(self).normalized_or(Self::X)
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}
impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}
impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}
impl SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
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
impl Index<usize> for Vec3 {
    type Output = f64;
    fn index(&self, i: usize) -> &f64 {
        match i {
            0 => &self.x,
            1 => &self.y,
            _ => &self.z,
        }
    }
}
impl From<[f64; 3]> for Vec3 {
    fn from(a: [f64; 3]) -> Self {
        Self::from_array(a)
    }
}
