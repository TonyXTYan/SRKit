//! Closed-form symmetric accelerate/decelerate trips, from rest to rest.
use crate::{kinematics::C, Error};

#[derive(Clone, Copy, Debug)]
pub struct FlipBurn {
    pub distance: f64,
    pub acceleration: f64,
    pub coordinate_time: f64,
    pub proper_time: f64,
    pub peak_rapidity: f64,
    pub peak_gamma: f64,
    /// Coordinate duration minus the light-travel time, evaluated directly to
    /// avoid cancellation when `coordinate_time` and `distance/C` round equal.
    pub light_time_excess: f64,
}

impl FlipBurn {
    pub fn new(distance: f64, acceleration: f64) -> Result<Self, Error> {
        if !distance.is_finite()
            || distance < 0.0
            || !acceleration.is_finite()
            || acceleration <= 0.0
        {
            return Err(Error::InvalidInput);
        }
        if distance == 0.0 {
            return Ok(Self {
                distance,
                acceleration,
                coordinate_time: 0.0,
                proper_time: 0.0,
                peak_rapidity: 0.0,
                peak_gamma: 1.0,
                light_time_excess: 0.0,
            });
        }
        let root_distance = libm::sqrt(distance);
        let root_acceleration = libm::sqrt(acceleration);
        let light_time = distance / C;
        // t² = (d/c)² + 4d/a. Taking square roots before division avoids
        // overflowing c/a or underflowing a*d for nonrelativistic trips.
        let newtonian_time = (2.0 * root_distance) / root_acceleration;
        let coordinate_time = libm::hypot(light_time, newtonian_time);
        // acosh(1+x) = 2 asinh(sqrt(x/2)), with x = a*d/(2c²).
        let half_sinh = (root_acceleration / (2.0 * C)) * root_distance;
        let asinh = libm::asinh(half_sinh);
        let peak_rapidity = 2.0 * asinh;
        let peak_gamma = 1.0 + (2.0 * half_sinh) * half_sinh;
        // tau = t_Newtonian * asinh(y)/y. The y=0 limit is exact; it
        // also handles positive products below the subnormal f64 range.
        let proper_time = newtonian_time
            * if half_sinh == 0.0 {
                1.0
            } else {
                (asinh / half_sinh).min(1.0)
            };
        // Rationalize t-d/c without squaring either term or subtracting
        // nearly equal times. This remains valid when 2/x would overflow.
        let light_time_excess = newtonian_time
            * ((newtonian_time / coordinate_time) / (1.0 + light_time / coordinate_time));
        if !coordinate_time.is_finite()
            || !proper_time.is_finite()
            || !light_time_excess.is_finite()
            || !peak_gamma.is_finite()
        {
            return Err(Error::Overflow);
        }
        Ok(Self {
            distance,
            acceleration,
            coordinate_time,
            proper_time,
            peak_rapidity,
            peak_gamma,
            light_time_excess,
        })
    }

    /// Ideal photon rocket rest-mass ratio for the two burns. This assumes
    /// perfect conversion to collimated photons and excludes payload margins.
    pub fn photon_rocket_mass_ratio(&self) -> Result<f64, Error> {
        let ratio = libm::exp(2.0 * self.peak_rapidity);
        if ratio.is_finite() {
            Ok(ratio)
        } else {
            Err(Error::Overflow)
        }
    }
}

/// Constant-speed circular gamma from laboratory radius and proper acceleration.
pub fn circular_gamma(radius: f64, proper_acceleration: f64) -> Result<f64, Error> {
    if !radius.is_finite()
        || radius <= 0.0
        || !proper_acceleration.is_finite()
        || proper_acceleration < 0.0
    {
        return Err(Error::InvalidInput);
    }
    // sqrt(a*r)/c can be finite even when a*r/c² overflows.
    let scaled = (libm::sqrt(proper_acceleration) / C) * libm::sqrt(radius);
    let gamma = libm::hypot(1.0, scaled);
    if gamma.is_finite() {
        Ok(gamma)
    } else {
        Err(Error::Overflow)
    }
}
