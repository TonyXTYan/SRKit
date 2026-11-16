//! Circular orbits flown under thrust: speed, period, the clocks' rates and what the pilot feels.

use core::f64::consts::PI;

use crate::{kinematics::Event, segment::Segment, vec3::Vec3};

/// A circle of radius `r` flown at constant speed `beta`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularOrbit {
    pub r: f64,
    /// Speed as a fraction of `c`.
    pub beta: f64,
    pub gamma: f64,
    /// Rapidity of the speed.
    pub phi: f64,
    /// Angular velocity in coordinate time, `beta / r`.
    pub omega: f64,
    /// Coordinate time for one lap.
    pub period: f64,
    /// Proper time for one lap.
    pub lap_tau: f64,
    /// Proper acceleration felt aboard: `gamma² beta² / r`.
    pub acc: f64,
}

impl CircularOrbit {
    /// From the radius and the speed.
    pub fn from_radius_speed(r: f64, beta: f64) -> Self {
        let gamma = 1.0 / libm::sqrt(1.0 - beta * beta);
        let period = 2.0 * PI * r / beta;
        Self {
            r,
            beta,
            gamma,
            phi: libm::atanh(beta),
            omega: beta / r,
            period,
            lap_tau: period / gamma,
            acc: gamma * gamma * beta * beta / r,
        }
    }
    /// From the radius and the proper acceleration held aboard: `gamma beta = sqrt(acc r)`.
    ///
    /// The square roots are taken separately so that `acc * r` cannot overflow or underflow first.
    pub fn from_radius_acceleration(r: f64, acc: f64) -> Self {
        let gb = libm::sqrt(acc) * libm::sqrt(r);
        let gamma = libm::hypot(1.0, gb);
        let beta = gb / gamma;
        let period = 2.0 * PI * r / beta;
        Self {
            r,
            beta,
            gamma,
            phi: libm::asinh(gb),
            omega: beta / r,
            period,
            lap_tau: period / gamma,
            acc,
        }
    }
    /// How far one clock sent each way round the loop would disagree at the meeting point, per lap, in the
    /// orbiting frame's proper time over `gamma`: the Sagnac gap `2 pi r beta gamma` (coordinate circumference
    /// times `beta gamma`).
    pub fn sagnac(&self) -> f64 {
        2.0 * PI * self.r * self.beta * self.gamma
    }
    /// Thomas precession per lap, in radians: a gyroscope carried round turns by `2 pi (gamma - 1)` against the
    /// stars.
    pub fn thomas_precession_per_lap(&self) -> f64 {
        2.0 * PI * self.beta * self.beta * self.gamma * self.gamma / (self.gamma + 1.0)
    }
    /// The smallest and largest Doppler factors a distant source in the orbit's plane shows over a lap:
    /// `gamma (1 - beta)` and `gamma (1 + beta)`.
    pub fn doppler_range(&self) -> (f64, f64) {
        (libm::exp(-self.phi), libm::exp(self.phi))
    }
    /// The state after proper time `tau` on a circle about the origin in the `xy` plane, starting at angle `a0`
    /// at `t = 0`: coordinate time, position, velocity, the unit heading and the angle.
    pub fn state(&self, a0: f64, tau: f64) -> OrbitState {
        let t = self.gamma * tau;
        let angle = a0 + self.omega * t;
        let (s, c) = (libm::sin(angle), libm::cos(angle));
        let heading = Vec3::new(-s, c, 0.0);
        OrbitState {
            t,
            x: Vec3::new(self.r * c, self.r * s, 0.0),
            v: heading * self.beta,
            heading,
            angle,
        }
    }
    /// The orbit as a segment about `centre` in the plane `(e1, e2)`, starting at angle `a0` at time `t0`.
    pub fn segment(&self, t0: f64, t1: f64, centre: Vec3, e1: Vec3, e2: Vec3, a0: f64) -> Segment {
        Segment::circular(t0, t1, centre, self.r, self.omega, e1, e2, a0)
    }
    /// Proper time at which the orbiting ship (about the origin in the `xy` plane, at angle `a0` at `t = 0`)
    /// receives light from `event` (`send = false`), or emits the light that reaches `event` (`send = true`).
    pub fn light_cone_tau(&self, a0: f64, event: Event, send: bool) -> f64 {
        let seg = self.segment(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, Vec3::Y, a0);
        let t = if send {
            seg.send(event, f64::NEG_INFINITY, f64::INFINITY)
        } else {
            seg.recv(event, f64::NEG_INFINITY, f64::INFINITY)
        };
        t.map_or(f64::NAN, |t| t / self.gamma)
    }
}

/// A point on a [`CircularOrbit`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrbitState {
    pub t: f64,
    pub x: Vec3,
    pub v: Vec3,
    pub heading: Vec3,
    pub angle: f64,
}
