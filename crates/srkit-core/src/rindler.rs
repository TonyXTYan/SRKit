//! The horizon of an accelerating observer.
//!
//! A ship that holds proper acceleration `alpha` for ever outruns light from part of spacetime: the region
//! beyond a plane that trails it at distance `1 / alpha` (`c² / alpha` in SI) and closes on it at the speed of
//! light. Signals sent from beyond that plane never arrive unless the ship eases off. The same closed form tells
//! a braking ship the latest moment a message from ahead can still catch it before it stops.

use crate::{kinematics::Event, segment::Hyperbolic, vec3::Vec3};

/// The Rindler horizon of a hyperbolic segment, extended for ever.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RindlerHorizon {
    /// The event the hyperbola is centred on.
    pub pivot: Event,
    /// Unit vector along the proper acceleration (the direction the ship is being pushed).
    pub n: Vec3,
    /// Distance from the ship to the horizon in the ship's own rest frame: `1 / |alpha|`.
    pub distance: f64,
}

impl RindlerHorizon {
    /// The horizon of `segment` if its acceleration never stopped.
    pub fn of(segment: &Hyperbolic) -> Self {
        let sign = if segment.alpha >= 0.0 { 1.0 } else { -1.0 };
        Self {
            pivot: segment.pivot(),
            n: segment.n * sign,
            distance: 1.0 / segment.alpha.abs(),
        }
    }
    /// How far `event` lies inside the region whose light still reaches the ship: `(x - x_p)·n - (t - t_p)`.
    /// Positive: the light arrives. Zero: on the horizon. Negative: it never arrives while the ship keeps
    /// accelerating.
    pub fn margin(&self, event: Event) -> f64 {
        (event.x - self.pivot.x).dot(self.n) - (event.t - self.pivot.t)
    }
    /// Whether light from `event` ever reaches the ship.
    pub fn can_signal(&self, event: Event) -> bool {
        self.margin(event) > 0.0
    }
    /// The last coordinate time at which a signal sent from the point `x` still reaches the ship. After it, the
    /// sender has dropped behind the horizon.
    pub fn last_signal_time(&self, x: Vec3) -> f64 {
        self.pivot.t + (x - self.pivot.x).dot(self.n)
    }
    /// How far `event` lies inside the region the ship's own signals can reach: `(x - x_p)·n + (t - t_p)`.
    /// Positive means some signal from the ship arrives there.
    pub fn margin_to_be_seen(&self, event: Event) -> f64 {
        (event.x - self.pivot.x).dot(self.n) + (event.t - self.pivot.t)
    }
}

/// The distance to the Rindler horizon for proper acceleration `alpha`: `1 / alpha` (`c² / alpha` in SI).
pub fn horizon_distance(alpha: f64) -> f64 {
    1.0 / alpha.abs()
}
