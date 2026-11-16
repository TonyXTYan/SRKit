//! # srkit-core
//!
//! Special relativity in 3+1 Minkowski spacetime, built from closed forms.
//!
//! * **Worldlines are piecewise analytic.** A [`segment::Segment`] is rest, inertial motion, constant proper
//!   acceleration (along the motion, or at any angle fixed aboard), a constant-thrust burn that spends mass, a
//!   circle at constant speed, or a light-like jump. Position, velocity and proper time are evaluated directly at
//!   any coordinate time; nothing is integrated step by step.
//! * **What an observer sees is a light-cone query.** [`worldline::Worldline::send`] finds the event on a
//!   worldline whose light reaches an observer, and [`worldline::Worldline::recv`] the event where light from
//!   somewhere reaches the worldline. Both are unique roots of a monotone function, found in closed form where
//!   one exists and by a bracketed Newton iteration otherwise.
//! * **Stored time is an exact integer.** [`time::Timestamp`] counts milliseconds in 384 or 512 bits, and
//!   [`causal`] decides "has this signal arrived by tick n?" in integer arithmetic.
//!
//! ## Units
//!
//! Geometry uses natural units with `c = 1`: pick any unit of time and measure lengths in the distance light
//! covers in it (seconds and light-seconds, Julian years and light-years). Velocities are fractions of `c`,
//! rapidities are dimensionless and proper accelerations are in inverse time. [`units`] holds the exact SI
//! constants and [`units::Units`] converts to and from SI.
//!
//! ## Numerical contract
//!
//! Every transcendental goes through the pure-Rust `libm`, there is no `unsafe`, no threads and no hash-map
//! iteration, so a given build gives the same bits on every platform that implements IEEE 754 binary64.
//! Geometry is floating point: an intersection found here is a construction, not an exact predicate. Exact
//! decisions belong to [`time`] and [`causal`].
#![no_std]
#![forbid(unsafe_code)]

#[cfg(feature = "alloc")]
extern crate alloc;

// pub mod causal;
// pub mod frames;
// pub mod kepler;
// pub mod kinematics;
// pub mod optics;
// pub mod orbit;
// pub mod plan;
// pub mod reach;
// pub mod rindler;
// pub mod rocket;
// pub mod roots;
// pub mod segment;
pub mod time;
// pub mod trip;
pub mod units;
pub mod vec3;
// pub mod worldline;

// pub use kinematics::Event;
// pub use segment::Segment;
pub use vec3::Vec3;
// pub use worldline::Worldline;

/// Why a constructor or an exact operation refused its input.
///
/// Evaluating a segment or a worldline never fails: those functions are total and follow IEEE 754 (a NaN in
/// gives a NaN out). Errors come from constructors that validate, from integer arithmetic that would wrap and
/// from queries that have no answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    /// A non-finite number, a speed at or above `c`, a zero direction, or another input outside the physics.
    InvalidInput,
    /// A query outside the domain where the function is defined.
    OutOfRange,
    /// Two worldline segments do not join in position or velocity.
    Discontinuous,
    /// The light cone does not meet the worldline in the range asked for.
    NoIntersection,
    /// The result does not fit the representation.
    Overflow,
    /// The answer cannot be given to the precision the function promises.
    PrecisionLoss,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => "invalid or non-finite physical input",
            Self::OutOfRange => "query outside the supported domain",
            Self::Discontinuous => "worldline segments do not join continuously",
            Self::NoIntersection => "no intersection in the supplied interval",
            Self::Overflow => "result exceeds the representation",
            Self::PrecisionLoss => "the result cannot be given to the promised precision",
        })
    }
}
