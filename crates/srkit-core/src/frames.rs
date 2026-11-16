//! Charts an observer can draw of spacetime: null coordinates, the momentary rest frame, radar coordinates and
//! light-cone ("as seen") coordinates.

use crate::{
    kinematics::{boost_rapidity, Event},
    vec3::Vec3,
    worldline::Worldline,
};

/// Null coordinates along one axis: `u = t - x`, `v = t + x`.
pub fn null_coordinates(t: f64, x: f64) -> (f64, f64) {
    (t - x, t + x)
}

/// Back from null coordinates: `t = (u + v) / 2`, `x = (v - u) / 2`.
pub fn from_null_coordinates(u: f64, v: f64) -> (f64, f64) {
    (0.5 * (u + v), 0.5 * (v - u))
}

/// A momentary rest frame: an origin event and the rapidity vector of the frame's velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RestFrame {
    pub origin: Event,
    /// Points along the frame's velocity; its length is the rapidity.
    pub rapidity: Vec3,
}

impl RestFrame {
    /// The coordinates of `event` in this frame, relative to its origin.
    pub fn coordinates(&self, event: Event) -> Event {
        let (t, x) = boost_rapidity(
            self.rapidity,
            event.t - self.origin.t,
            event.x - self.origin.x,
        );
        Event { t, x }
    }
    /// The event with the given coordinates in this frame, back in the original frame.
    pub fn event(&self, local: Event) -> Event {
        let (t, x) = boost_rapidity(-self.rapidity, local.t, local.x);
        Event {
            t: self.origin.t + t,
            x: self.origin.x + x,
        }
    }
}

/// Radar coordinates of an event for an observer: the observer's proper times when a ping would have to leave
/// to reach the event, and when its echo comes back.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Radar {
    /// Proper time at which the observer's light reaches the event if sent then.
    pub tau_send: f64,
    /// Proper time at which the event's light reaches the observer.
    pub tau_recv: f64,
}

impl Radar {
    /// Radar time: the midpoint, `(tau_recv + tau_send) / 2`. On the observer's own worldline it is proper time,
    /// and it is smooth through accelerations.
    pub fn time(&self) -> f64 {
        0.5 * (self.tau_recv + self.tau_send)
    }
    /// Radar distance: half the round trip, `(tau_recv - tau_send) / 2`.
    pub fn distance(&self) -> f64 {
        0.5 * (self.tau_recv - self.tau_send)
    }
}

/// Radar coordinates of `event` for the observer on `worldline`. A side is infinite where the light cone never
/// meets the worldline (an event behind a Rindler horizon has no echo).
pub fn radar(worldline: &Worldline<'_>, event: Event) -> Radar {
    let t_send = worldline.send(event);
    let t_recv = worldline.recv(event);
    Radar {
        tau_send: if t_send.is_finite() {
            worldline.tau(t_send)
        } else {
            f64::NEG_INFINITY
        },
        tau_recv: if t_recv.is_finite() {
            worldline.tau(t_recv)
        } else {
            f64::INFINITY
        },
    }
}

/// Light-cone coordinates: the observer's proper time when the event is seen, and the event's distance in the
/// observer's rest frame at that moment (its light age there).
pub fn light_cone_coordinates(worldline: &Worldline<'_>, event: Event) -> Option<(f64, f64)> {
    let t_recv = worldline.recv(event);
    if !t_recv.is_finite() {
        return None;
    }
    let state = worldline.at(t_recv);
    let phi = libm::atanh(state.v.norm().min(1.0 - 1e-16));
    let rapidity = if phi > 0.0 {
        state.v / state.v.norm() * phi
    } else {
        Vec3::ZERO
    };
    let (dt, _) = boost_rapidity(rapidity, event.t - t_recv, event.x - state.x);
    Some((worldline.tau(t_recv), -dt))
}
