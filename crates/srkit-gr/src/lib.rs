//! # srkit-gr
//!
//! Closed-form general-relativity scalars: numbers that have an exact formula, evaluated directly.
//!
//! * [`schwarzschild`]: a non-rotating hole or any spherical body seen from outside. Characteristic radii,
//!   the clock rate of a static observer, the redshift between two static radii, and the clock rate and
//!   angular velocity on a circular orbit.
//! * [`kerr`]: a rotating hole. Horizons, the ergosphere, frame dragging, the clock rates of static and
//!   zero-angular-momentum observers, equatorial circular orbits, the innermost stable circular orbit
//!   (Bardeen, Press and Teukolsky 1972) with its energy, angular momentum, clock rate and accretion
//!   efficiency, and the rotational energy that can be extracted.
//! * [`weak_field`]: first post-Newtonian results for the Solar System and for lensing. Shapiro delay, light
//!   deflection, the Einstein angle, and clock-rate offsets, including an orbiting clock against a surface
//!   clock.
//! * [`hawking`]: temperature, entropy, power, mass-loss rate and lifetime of an evaporating hole. Power and
//!   lifetime depend on which particles are emitted, so every such function takes a
//!   [`hawking::EmissionModel`].
//! * [`newtonian`]: the Newtonian helpers the physics oracle uses. Binding energy of a uniform sphere,
//!   circular and escape speeds, orbital and synodic periods, Hohmann transfers.
//!
//! There is no geodesic integration, no ray tracing and no metric algebra here. Those belong to EinsteinPy and
//! the Black Hole Perturbation Toolkit; see the crate README.
//!
//! ## Units
//!
//! Every argument and result is in SI units: masses in kilograms, lengths in metres, times in seconds, angles
//! in radians, angular velocities in radians per second, energies in joules, powers in watts and temperatures
//! in kelvin. Each function states the units of its arguments and of its result. Ratios such as `dt/dτ`, an
//! efficiency or a redshift are dimensionless. The spin of a Kerr hole is the dimensionless
//! `a* = J c / (G M²)`, between −1 and 1.
//!
//! Radii are coordinate radii: the Schwarzschild areal radius (a sphere at `r` has area `4πr²`) or the
//! Boyer-Lindquist radius of the Kerr metric. They are not proper distances.
//!
//! In the Schwarzschild, Kerr, weak-field and orbit functions a mass enters only through `G M`. For a body
//! whose gravitational parameter is known better than its mass (the Sun, the planets), pass
//! [`mass_from_gravitational_parameter`]`(GM)`; the functions multiply by the same `G` again, so the tabulated
//! `GM` is recovered to one rounding error. [`gravitational_radius`] and [`gravitational_time`] give the
//! length `GM/c²` and the time `GM/c³` that turn a result into the geometrised units (`G = c = 1`, lengths
//! in `M`) of the literature.
//!
//! ## Numerical contract
//!
//! * `no_std`, no allocation, no `unsafe`. The only operations are IEEE 754 binary64 arithmetic and the
//!   functions of the pure-Rust `libm`, so the results do not depend on the platform's maths library and
//!   should be the same bits on every architecture. That was checked once by hand (16 Nov 2026): the 1,400
//!   reference inputs gave identical bits natively on x86-64 and as `wasm32-unknown-unknown` under V8. It is
//!   not yet an automated test, and no other target has been run.
//! * Functions are total. Outside the stated domain, and for a NaN argument, the result is NaN. At the edge of
//!   a domain the IEEE limit is returned: a clock rate of 0 on a horizon, `+∞` for the redshift of light from
//!   it.
//! * Each function states its formula, its domain, a reference and its accuracy. "Within a few ULP" means that
//!   the relative error against 120-digit references (`tools/crosscheck.py`, `tests/crosscheck.rs`) is below
//!   `2e-15`; the largest measured where the result is well conditioned is `7.4e-16`. Formulas are arranged
//!   so that small differences (a redshift near zero, a clock-rate offset, the distance to a horizon) are not
//!   formed by subtracting nearly equal numbers. Where a result is ill conditioned in its arguments, the
//!   function says so.
//!
//! ## Example
//!
//! ```
//! use srkit_core::units::{JULIAN_YEAR, M_SUN};
//! use srkit_gr::hawking::{self, EmissionModel};
//! use srkit_gr::kerr::{self, Orbit};
//! use srkit_gr::gravitational_radius;
//!
//! // A station on the prograde ISCO of a hole of four million solar masses spinning at a* = 0.998.
//! let mass = 4.0e6 * M_SUN; // kg
//! let radius = kerr::isco_radius(mass, 0.998, Orbit::Prograde); // m
//! assert!((radius / gravitational_radius(mass) - 1.237).abs() < 1e-3); // 1.237 GM/c²
//!
//! // 10.8 years pass far away for each year aboard, and matter that spirals in to it releases 32.1% of mc².
//! assert!((kerr::isco_dt_dtau(0.998, Orbit::Prograde) - 10.79).abs() < 0.01);
//! assert!((kerr::accretion_efficiency(0.998, Orbit::Prograde) - 0.321).abs() < 1e-3);
//!
//! // The hole's evaporation time, with the emission model named: 7.4e86 years.
//! let years = hawking::lifetime(mass, EmissionModel::PagePhotonsGravitons) / JULIAN_YEAR;
//! assert!((years / 7.43e86 - 1.0).abs() < 0.01);
//! ```
#![no_std]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

pub mod hawking;
pub mod kerr;
pub mod newtonian;
pub mod schwarzschild;
pub mod weak_field;

use srkit_core::units::{C, G_NEWTON};

/// `c²`, m²/s².
pub(crate) const C2: f64 = C * C;
/// `c³`, m³/s³.
pub(crate) const C3: f64 = C * C * C;

/// The gravitational parameter `μ = G M`, in m³/s².
///
/// * `mass`: kg.
///
/// Inverse of [`mass_from_gravitational_parameter`]. Exact to one rounding error.
pub fn gravitational_parameter(mass: f64) -> f64 {
    G_NEWTON * mass
}

/// The mass, in kg, that this crate's `G` assigns to a gravitational parameter `GM`.
///
/// * `gm`: the gravitational parameter, m³/s².
///
/// `GM` is measured far better than `G` or `M` separately (the Sun's to about 1 part in 10¹⁰, `G` to 2.2 parts
/// in 10⁵). Passing `mass_from_gravitational_parameter(gm)` as a mass makes every function here use the
/// measured `GM`, because each multiplies by the same `G` again.
pub fn mass_from_gravitational_parameter(gm: f64) -> f64 {
    gm / G_NEWTON
}

/// The gravitational radius `GM/c²`, in metres: the length called `M` in geometrised units.
///
/// * `mass`: kg.
///
/// Half the Schwarzschild radius. 1,476.6 m for the Sun, 4.435 mm for the Earth. Divide a radius by this to
/// get it "in units of `M`". Exact to a few ULP.
pub fn gravitational_radius(mass: f64) -> f64 {
    G_NEWTON * mass / C2
}

/// The gravitational time `GM/c³`, in seconds: the light-crossing time of one gravitational radius.
///
/// * `mass`: kg.
///
/// 4.925 µs for the Sun. Multiply an angular velocity by this to get it in units of `1/M`. Exact to a few
/// ULP.
pub fn gravitational_time(mass: f64) -> f64 {
    G_NEWTON * mass / C3
}

/// `1 − inner/r`, arranged so that the difference is exact when `r` is close to `inner`.
///
/// For `inner/r > 1/2` the subtraction `r − inner` is exact (Sterbenz), so the result keeps full relative
/// accuracy all the way down to a horizon. For smaller ratios `1 − inner/r` is already accurate, and it also
/// handles `r = ∞`.
pub(crate) fn one_minus_ratio(inner: f64, r: f64) -> f64 {
    let x = inner / r;
    if x > 0.5 {
        (r - inner) / r
    } else {
        1.0 - x
    }
}
