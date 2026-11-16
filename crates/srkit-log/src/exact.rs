//! Arrivals on the integer clock.
//!
//! The worldline gives an arrival time in floating point. The engine stores time as whole milliseconds and must
//! never deliver a signal a tick early. [`arrival_tick`] takes the floating-point answer as an estimate and
//! settles it with the exact integer check `c dt >= |dx|` from [`srkit_core::causal`].

use srkit_core::{
    causal::{self, GridVec},
    kinematics::Event,
    time::{Duration, Timestamp},
    units::Units,
    worldline::Worldline,
    Error,
};

/// How a floating-point scene maps onto the integer clock and the metre grid.
#[derive(Clone, Copy, Debug)]
pub struct Frame<const N: usize> {
    /// The stored tick at which the scene's coordinate time is `t_epoch`.
    pub epoch: Timestamp<N>,
    /// The scene's coordinate time at the epoch tick.
    pub t_epoch: f64,
    /// The scene's units (its unit of time in seconds).
    pub units: Units,
}

impl<const N: usize> Frame<N> {
    /// The scene's coordinate time at a tick. The difference from the epoch is taken in integers first.
    pub fn time_at(&self, tick: Timestamp<N>) -> Result<f64, Error> {
        Ok(self.t_epoch
            + self
                .units
                .time_from_si(tick.local_seconds_since(self.epoch)?))
    }
    /// The tick at or after a scene time (rounded up to a whole millisecond).
    pub fn tick_at(&self, t: f64) -> Result<Timestamp<N>, Error> {
        let ms = libm_ceil(self.units.time_to_si(t - self.t_epoch) * 1e3);
        if !ms.is_finite() || ms.abs() > 9.0e15 {
            return Err(Error::PrecisionLoss);
        }
        self.epoch.checked_add(Duration::from_millis(ms as i64))
    }
    /// A scene position on the metre grid.
    pub fn grid(&self, x: srkit_core::Vec3) -> Result<GridVec, Error> {
        let m = self.units.length_m();
        GridVec::from_metres(x.x * m, x.y * m, x.z * m)
    }
}

fn libm_ceil(x: f64) -> f64 {
    x.ceil()
}

/// The first tick at which light from `event` has reached `receiver`, exactly.
///
/// The floating-point root is the estimate; the result is the first tick at which the receiver's position,
/// rounded to the metre grid, lies inside the event's future light cone by the integer check. `guard` bounds how
/// many ticks the estimate may be off (a few are plenty: the estimate is good to rounding). Returns
/// [`Error::NoIntersection`] if the light never arrives.
pub fn arrival_tick<const N: usize>(
    frame: &Frame<N>,
    event: Event,
    receiver: &Worldline<'_>,
    guard: u32,
) -> Result<Timestamp<N>, Error> {
    let t = receiver.recv(event);
    if !t.is_finite() {
        return Err(Error::NoIntersection);
    }
    // The emission is the last tick at or before the event, so that rounding the emission cannot make the signal
    // early either.
    let emission = {
        let up = frame.tick_at(event.t)?;
        if frame.time_at(up)? > event.t {
            up.checked_sub(Duration::from_unsigned_millis(1))?
        } else {
            up
        }
    };
    let source = frame.grid(event.x)?;
    causal::settle_arrival(emission, frame.tick_at(t)?, guard, |tick| {
        source.to(frame.grid(receiver.at(frame.time_at(tick)?).x)?)
    })
}
