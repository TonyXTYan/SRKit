//! Experimental SR primitives. Coordinates use SI units and local epochs.
//!
//! Integer clocks and floating-point local geometry have separate contracts:
//! a bounded floating-point intersection is not an exact causal predicate.
#![no_std]
#![forbid(unsafe_code)]

pub mod kinematics;
pub mod lightcone;
pub mod navigator;
mod navigator_optics;
pub mod segment;
pub mod time;
pub mod trip;

/// Errors are explicit rather than propagating NaNs into causal decisions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    InvalidInput,
    OutOfRange,
    Discontinuous,
    NoIntersection,
    Overflow,
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
            Self::PrecisionLoss => "local interval exceeds floating-point precision contract",
        })
    }
}
