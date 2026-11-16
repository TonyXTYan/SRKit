//! Newtonian gravity helpers: the non-relativistic formulas the physics oracle uses.
//!
//! Point masses, two bodies, no relativity. They are here because the same game systems that need the
//! relativistic scalars also need these, and because they are the limits the relativistic results are
//! checked against. The fractional error from leaving relativity out is of order `GM/(r c²)`: `1e-8` for a
//! planet's orbit round the Sun.
//!
//! Masses are in kg, lengths in metres, times in seconds, speeds in m/s and energies in joules. Functions
//! return NaN outside the stated domain. "Within a few ULP" is the accuracy with which the stated formula is
//! evaluated.
//!
//! References: R. R. Bate, D. D. Mueller and J. E. White, *Fundamentals of Astrodynamics* (1971), chapters 1,
//! 3 and 8; W. Hohmann, *Die Erreichbarkeit der Himmelskörper* (1925); S. Chandrasekhar, *An Introduction to
//! the Study of Stellar Structure* (1939), chapter IV.

use core::f64::consts::PI;
use libm::{fabs, sqrt};
use srkit_core::units::G_NEWTON;

/// The gravitational binding energy of a uniform sphere, `U = (3/5) G M²/R`, in joules.
///
/// * `mass`: kg, `≥ 0`.
/// * `radius`: m, `> 0`.
///
/// The energy needed to take a sphere of uniform density apart and move every piece far away. `2.24e32 J`
/// for an Earth of uniform density. Real bodies are denser towards the centre and are bound more tightly:
/// the Earth by about 11% more, and Jupiter, as an `n = 1` polytrope with `(3/4) G M²/R`, by 25% more. Within
/// a few ULP.
///
/// Reference: Chandrasekhar (1939), chapter IV (the potential energy of a polytrope, `3GM²/((5 − n)R)`, at
/// index `n = 0`); the oracle rows SR-47 to SR-49.
pub fn binding_energy_uniform_sphere(mass: f64, radius: f64) -> f64 {
    if mass >= 0.0 && radius > 0.0 {
        0.6 * G_NEWTON * mass * mass / radius
    } else {
        f64::NAN
    }
}

/// The speed on a circular orbit, `v = √(GM/r)`, in m/s.
///
/// * `mass`: the central mass, kg, `≥ 0`; the orbiting body's mass is neglected.
/// * `r`: radius of the orbit, m, `> 0`.
///
/// 29.78 km/s for the Earth round the Sun. Within a few ULP.
///
/// Reference: Bate, Mueller and White (1971), chapter 1.
pub fn circular_speed(mass: f64, r: f64) -> f64 {
    if mass >= 0.0 && r > 0.0 {
        sqrt(G_NEWTON * mass / r)
    } else {
        f64::NAN
    }
}

/// The escape speed, `v = √(2GM/r)`, in m/s.
///
/// * `mass`: kg, `≥ 0`.
/// * `r`: distance from the centre of the mass, m, `> 0`.
///
/// The speed at which a body at `r` has just enough energy to reach infinity; `√2` times
/// [`circular_speed`]. 11.19 km/s from the Earth's surface. It equals `c` at the Schwarzschild radius.
/// Within a few ULP.
///
/// Reference: Bate, Mueller and White (1971), chapter 1.
pub fn escape_speed(mass: f64, r: f64) -> f64 {
    if mass >= 0.0 && r > 0.0 {
        sqrt(2.0 * G_NEWTON * mass / r)
    } else {
        f64::NAN
    }
}

/// The period of an orbit, `T = 2π √(a³/(GM))`, in seconds (Kepler's third law).
///
/// * `mass`: the central mass, kg, `> 0`; the orbiting body's mass is neglected.
/// * `r`: radius of a circular orbit, or the semi-major axis `a` of an elliptical one, m, `≥ 0`.
///
/// 365.2569 days for 1 au round the Sun. The same formula is exact for a circular orbit in the Schwarzschild
/// metric when `r` is the areal radius and the period is measured by a static clock far away. Within a few
/// ULP.
///
/// Reference: Bate, Mueller and White (1971), chapter 1.
pub fn circular_period(mass: f64, r: f64) -> f64 {
    if mass > 0.0 && r >= 0.0 {
        2.0 * PI * r * sqrt(r / (G_NEWTON * mass))
    } else {
        f64::NAN
    }
}

/// The synodic period of two bodies on circular orbits round the same centre,
/// `1/|1/T₁ − 1/T₂|`, in the unit of the two periods.
///
/// * `period_a`: orbital period of one body, any unit of time, `> 0`.
/// * `period_b`: orbital period of the other body, same unit, `> 0`.
///
/// The time between successive identical alignments of the two bodies and the centre (opposition to
/// opposition), and so the interval between launch windows for a transfer between them. 779.9 days for the
/// Earth and Mars. `+∞` for equal periods. Evaluated as `T₁ T₂/|T₂ − T₁|`; within a few ULP unless the
/// two periods are nearly equal.
///
/// Reference: Bate, Mueller and White (1971), chapter 8.
pub fn synodic_period(period_a: f64, period_b: f64) -> f64 {
    if period_a > 0.0 && period_b > 0.0 {
        period_a * (period_b / fabs(period_b - period_a))
    } else {
        f64::NAN
    }
}

/// The flight time of a Hohmann transfer, `π √(a³/(GM))` with `a = (r₁ + r₂)/2`, in seconds.
///
/// * `mass`: the central mass, kg, `> 0`.
/// * `r_from`: radius of the starting circular orbit, m, `> 0`.
/// * `r_to`: radius of the target circular orbit, m, `> 0`.
///
/// Half the period of the ellipse that touches both circular, coplanar orbits. 258.9 days from the Earth's
/// orbit to that of Mars. The same in both directions. Within a few ULP.
///
/// Reference: Hohmann (1925); Bate, Mueller and White (1971), chapter 3.
pub fn hohmann_transfer_time(mass: f64, r_from: f64, r_to: f64) -> f64 {
    if mass > 0.0 && r_from > 0.0 && r_to > 0.0 {
        let a = 0.5 * (r_from + r_to);
        PI * a * sqrt(a / (G_NEWTON * mass))
    } else {
        f64::NAN
    }
}

/// The two burns of a Hohmann transfer between circular, coplanar orbits.
///
/// Each value is the change of speed along the direction of motion, in m/s: positive speeds the craft up
/// (a prograde burn), negative slows it down. Both are positive for an outward transfer and both negative
/// for an inward one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HohmannDeltaV {
    /// The burn that leaves the starting orbit, m/s:
    /// `√(GM/r₁) (√(2 r₂/(r₁ + r₂)) − 1)`.
    pub departure: f64,
    /// The burn that enters the target orbit, m/s:
    /// `√(GM/r₂) (1 − √(2 r₁/(r₁ + r₂)))`.
    pub arrival: f64,
}

impl HohmannDeltaV {
    /// The total change of speed the craft must supply, `|departure| + |arrival|`, in m/s.
    pub fn total(self) -> f64 {
        fabs(self.departure) + fabs(self.arrival)
    }
}

/// The changes of speed for a Hohmann transfer between two circular, coplanar orbits, in m/s.
///
/// * `mass`: the central mass, kg, `≥ 0`.
/// * `r_from`: radius of the starting circular orbit, `r₁`, m, `> 0`.
/// * `r_to`: radius of the target circular orbit, `r₂`, m, `> 0`.
///
/// From the Earth's orbit to that of Mars round the Sun: 2.94 km/s to leave and 2.65 km/s to arrive,
/// 5.59 km/s in total. These are heliocentric changes of speed between circular orbits; climbing out of and
/// into the planets' own gravity is not included. The Hohmann transfer has the lowest total of any two-burn
/// transfer for `r₂/r₁` below 11.94.
///
/// Each burn is evaluated without subtracting nearly equal speeds, so it keeps its relative accuracy for
/// neighbouring orbits: within a few ULP, limited only by how well `r_to − r_from` is known.
///
/// Reference: Hohmann (1925); Bate, Mueller and White (1971), chapter 3.
pub fn hohmann_delta_v(mass: f64, r_from: f64, r_to: f64) -> HohmannDeltaV {
    if mass >= 0.0 && r_from > 0.0 && r_to > 0.0 {
        let gm = G_NEWTON * mass;
        let sum = r_from + r_to;
        let gap = r_to - r_from;
        // √(2 r₂/(r₁+r₂)) − 1 = ((r₂ − r₁)/(r₁+r₂)) / (√(2 r₂/(r₁+r₂)) + 1), and likewise for arrival.
        let k_from = sqrt(2.0 * r_to / sum);
        let k_to = sqrt(2.0 * r_from / sum);
        HohmannDeltaV {
            departure: sqrt(gm / r_from) * (gap / sum) / (k_from + 1.0),
            arrival: sqrt(gm / r_to) * (gap / sum) / (1.0 + k_to),
        }
    } else {
        HohmannDeltaV {
            departure: f64::NAN,
            arrival: f64::NAN,
        }
    }
}
