//! Weak-field results: the first post-Newtonian corrections around a slowly moving mass.
//!
//! Everything here is first order in the potential `GM/(r c²)`, which is `2.1e-6` at the Sun's surface,
//! `7.0e-10` at the Earth's and about `1e-6` across the Galaxy. The neglected terms are smaller than the
//! result by a factor of that order. Close to a black hole or a neutron star use [`crate::schwarzschild`] or
//! [`crate::kerr`] instead.
//!
//! Light is treated as moving along the straight line it would follow without the mass, in isotropic (or
//! harmonic) coordinates, with general relativity's value `γ = 1` of the post-Newtonian parameter. The
//! mass is taken as a point: rotation, oblateness and motion of the body are left out.
//!
//! Clock rates are returned as **offsets**, `dτ/dt − 1`, because the rates themselves differ from 1 by parts
//! in a billion and would lose most of their digits if returned as `1 − small`. A negative offset is a slow
//! clock. Multiply an offset by an interval to get the accumulated difference: `4.4647e-10 × 86,400 s =
//! 38.6 µs` per day for a GPS clock.
//!
//! Masses are in kg, lengths in metres, times in seconds and angles in radians. Functions return NaN outside
//! the stated domain. "Within a few ULP" below is the accuracy with which the stated first-order formula is
//! evaluated, not the accuracy of the first-order approximation.
//!
//! References: C. M. Will, "The confrontation between general relativity and experiment", *Living Rev.
//! Relativ.* 17 (2014) 4, sections "The deflection of light" and "The time delay of light" (cited below as
//! Will 2014); N. Ashby, "Relativity in the Global Positioning System", *Living Rev. Relativ.* 6 (2003) 1,
//! section "Relativistic effects on satellite clocks" (Ashby 2003); G. Petit and B. Luzum (eds.), *IERS
//! Conventions (2010)*, IERS Technical Note 36, chapter 11 (IERS 2010).

use crate::{C2, C3};
use libm::{asinh, atanh, sqrt};
use srkit_core::units::G_NEWTON;

/// The Shapiro delay along one leg of a light path, `(2GM/c³) arcosh(r/b)`, in seconds.
///
/// * `mass`: kg, `≥ 0`.
/// * `r`: distance of the end of the leg from the mass, m, `r ≥ b`.
/// * `impact_parameter`: `b`, the closest distance of the straight path to the mass, m, `> 0`.
///
/// The extra coordinate time light takes between the point of closest approach and a point at distance `r`,
/// compared with the Euclidean path length over `c`. `arcosh(r/b) = ln((r + √(r² − b²))/b)`, which tends to
/// `ln(2r/b)` for `r ≫ b`. A path whose closest approach lies between its ends is two such legs, see
/// [`shapiro_delay`]; if both ends are on the same side of the closest approach, the delay is the difference
/// of two legs, or use [`shapiro_delay_between`].
///
/// Valid for `b ≫ GM/c²`. Evaluated as `arsinh(√((r − b)(r + b))/b)`, in which `r − b` is exact, so the
/// result is within a few ULP for a short leg (`r` close to `b`) as well as a long one.
///
/// Reference: I. I. Shapiro, *Phys. Rev. Lett.* 13 (1964) 789; Will 2014.
pub fn shapiro_delay_leg(mass: f64, r: f64, impact_parameter: f64) -> f64 {
    if mass >= 0.0 && impact_parameter > 0.0 && r >= impact_parameter {
        let b = impact_parameter;
        // The length of the leg in units of b; the two factors are kept apart so that neither overflows.
        let leg = sqrt((r - b) / b) * sqrt((r + b) / b);
        2.0 * G_NEWTON * mass / C3 * asinh(leg)
    } else {
        f64::NAN
    }
}

/// The one-way Shapiro delay of light passing a mass, in seconds:
/// `(2GM/c³) [arcosh(r_emitter/b) + arcosh(r_receiver/b)]`.
///
/// * `mass`: kg, `≥ 0`.
/// * `r_emitter`: distance of the emitter from the mass, m, `≥ b`.
/// * `r_receiver`: distance of the receiver from the mass, m, `≥ b`.
/// * `impact_parameter`: `b`, the closest distance of the straight path to the mass, m, `> 0`.
///
/// The extra coordinate time, over the Euclidean distance divided by `c`, for a path whose closest approach
/// to the mass lies between the emitter and the receiver (a superior conjunction): the sum of two
/// [`shapiro_delay_leg`]s. For `b` much smaller than both distances it is the familiar
/// `(2GM/c³) ln(4 r_emitter r_receiver/b²)`. A radar echo takes twice this.
///
/// For a signal between the Earth and Mars that grazes the Sun the one-way delay is 124 µs and the round
/// trip 247 µs. Valid for `b ≫ GM/c²`. Within a few ULP.
///
/// Reference: Shapiro (1964); Will 2014.
pub fn shapiro_delay(mass: f64, r_emitter: f64, r_receiver: f64, impact_parameter: f64) -> f64 {
    shapiro_delay_leg(mass, r_emitter, impact_parameter)
        + shapiro_delay_leg(mass, r_receiver, impact_parameter)
}

/// The one-way Shapiro delay between two points given their separation, in seconds:
/// `(2GM/c³) ln((r_emitter + r_receiver + R)/(r_emitter + r_receiver − R))`.
///
/// * `mass`: kg, `≥ 0`.
/// * `r_emitter`: distance of the emitter from the mass, m, `> 0`.
/// * `r_receiver`: distance of the receiver from the mass, m, `> 0`.
/// * `separation`: `R`, the straight-line distance between emitter and receiver, m, with
///   `|r_emitter − r_receiver| ≤ R < r_emitter + r_receiver`.
///
/// The same first-order delay as [`shapiro_delay`], written for any geometry: the closest approach need not
/// lie between the two points, and no impact parameter is needed. When it does lie between them the two
/// forms are equal. Evaluated as `(4GM/c³) artanh(R/(r_emitter + r_receiver))`, which stays accurate for
/// nearby points.
///
/// For a long path that passes close to the mass, `r_emitter + r_receiver − R` is a small difference of
/// large numbers and the result is only as good as that difference; [`shapiro_delay`] with the impact
/// parameter is the better-conditioned form there. Otherwise within a few ULP.
///
/// Reference: IERS 2010, the light-time equation for ranging; Will 2014.
pub fn shapiro_delay_between(mass: f64, r_emitter: f64, r_receiver: f64, separation: f64) -> f64 {
    let sum = r_emitter + r_receiver;
    let gap = r_emitter - r_receiver;
    if mass >= 0.0
        && r_emitter > 0.0
        && r_receiver > 0.0
        && separation >= gap
        && separation >= -gap
        && separation < sum
    {
        4.0 * G_NEWTON * mass / C3 * atanh(separation / sum)
    } else {
        f64::NAN
    }
}

/// The deflection angle of light passing a mass, `α = 4GM/(b c²)`, in radians.
///
/// * `mass`: kg, `≥ 0`.
/// * `impact_parameter`: `b`, the closest distance of the undeflected path to the mass, m, `> 0`.
///
/// The total angle through which a ray from far away is bent towards the mass: twice the Newtonian value.
/// 1.75 arcseconds at the Sun's limb. Valid for `b ≫ GM/c²`: the next term is smaller by a factor of order
/// `GM/(b c²)`, and near a black hole's photon sphere the true deflection diverges. Within a few ULP.
///
/// Reference: A. Einstein, *Ann. Phys.* 49 (1916) 769; Will 2014.
pub fn deflection_angle(mass: f64, impact_parameter: f64) -> f64 {
    if mass >= 0.0 && impact_parameter > 0.0 {
        4.0 * G_NEWTON * mass / (impact_parameter * C2)
    } else {
        f64::NAN
    }
}

/// The Einstein angle of a point lens, `θ_E = √((4GM/c²) D_ls/(D_l D_s))`, in radians.
///
/// * `mass`: kg, `≥ 0`.
/// * `d_lens`: distance from the observer to the lens, `D_l`, m, `> 0`.
/// * `d_source`: distance from the observer to the source, `D_s`, m, `> 0`.
/// * `d_lens_source`: distance from the lens to the source, `D_ls`, m, `≥ 0`.
///
/// The angular radius of the ring into which a source exactly behind a point mass is imaged, and the angular
/// scale of lensing by that mass otherwise. In flat space `d_lens_source = d_source − d_lens`. Across
/// cosmological distances pass angular-diameter distances, for which that difference does not hold. The
/// radius of the ring in the lens plane is `θ_E d_lens`.
///
/// 1.0 milliarcsecond for a solar mass halfway to a source 8 kpc away. Valid when the ring radius
/// `θ_E d_lens` is much larger than `GM/c²` and all angles are small. Within a few ULP.
///
/// Reference: A. Einstein, *Science* 84 (1936) 506; P. Schneider, J. Ehlers and E. E. Falco, *Gravitational
/// Lenses* (1992), chapter 2.
pub fn einstein_angle(mass: f64, d_lens: f64, d_source: f64, d_lens_source: f64) -> f64 {
    if mass >= 0.0 && d_lens > 0.0 && d_source > 0.0 && d_lens_source >= 0.0 {
        sqrt(4.0 * G_NEWTON * mass / C2 * (d_lens_source / d_source) / d_lens)
    } else {
        f64::NAN
    }
}

/// The clock-rate offset of a static clock in a weak field, `dτ/dt − 1 = −GM/(r c²)`, dimensionless.
///
/// * `mass`: kg, `≥ 0`.
/// * `r`: distance of the clock from the centre of the mass, m, `> 0`.
///
/// `τ` is the proper time of a clock at rest at distance `r` and `t` the time of a clock at rest far away.
/// The offset is negative: the lower clock runs slow. On the surface of a body of radius `R` it is
/// `−GM/(R c²)`: `−2.12e-6` for the Sun, `−6.96e-10` for the Earth. The redshift of light climbing from `r`
/// to infinity is the same number with the sign reversed.
///
/// First post-Newtonian order: the exact Schwarzschild value, `√(1 − 2GM/(r c²)) − 1`, differs by about
/// `(GM/(r c²))²/2`. Within a few ULP.
///
/// Reference: Ashby 2003; Will 2014.
pub fn static_rate_offset(mass: f64, r: f64) -> f64 {
    if mass >= 0.0 && r > 0.0 {
        -(G_NEWTON * mass) / (r * C2)
    } else {
        f64::NAN
    }
}

/// The clock-rate offset of a clock on a circular orbit in a weak field,
/// `dτ/dt − 1 = −3GM/(2 r c²)`, dimensionless.
///
/// * `mass`: kg, `≥ 0`.
/// * `r`: radius of the orbit, m, `> 0`.
///
/// Relative to a clock at rest far away. Two thirds of the offset is gravitational, `−GM/(r c²)`, and one
/// third is the time dilation of the orbital speed, `−v²/(2c²)` with `v² = GM/r`.
///
/// First post-Newtonian order: the exact Schwarzschild value, `√(1 − 3GM/(r c²)) − 1`, differs by about
/// `(9/8)(GM/(r c²))²`. Within a few ULP.
///
/// Reference: Ashby 2003.
pub fn circular_orbit_rate_offset(mass: f64, r: f64) -> f64 {
    if mass >= 0.0 && r > 0.0 {
        -1.5 * (G_NEWTON * mass) / (r * C2)
    } else {
        f64::NAN
    }
}

/// The rate at which a clock on a circular orbit gains on a clock on the body's surface, dimensionless:
/// `GM/(R c²) + v²/(2c²) − 3GM/(2 r c²)`.
///
/// * `mass`: kg, `≥ 0`.
/// * `orbit_radius`: `r`, radius of the circular orbit, m, `> 0`.
/// * `surface_radius`: `R`, distance of the surface clock from the centre of the body, m, `> 0`.
/// * `surface_speed`: `v`, speed of the surface clock in the non-rotating frame of the body (its rotation
///   speed there), m/s, `≥ 0`; pass 0 for a non-rotating body.
///
/// This is `(dτ_orbit − dτ_surface)/dt`, the sum used for satellite navigation clocks: the orbiting clock is
/// higher in the potential, which makes it fast, and moves faster, which makes it slow. A positive result
/// means the orbiting clock runs fast. The result is zero for an orbit at `r = 3R/2` above a non-rotating
/// body, negative below that and positive above.
///
/// For a GPS satellite (`r` = 26,561.75 km) against a clock on the Earth's equator (`R` = 6,378.137 km,
/// `v` = 465.1 m/s) the result is `4.461e-10`, which is +38.5 µs per day: +45.7 µs from the potential and
/// −7.2 µs from the orbital speed, with +0.1 µs from the rotation of the surface clock. The published figure
/// of `4.4647e-10` (+38.6 µs per day) uses the potential of the geoid, `W₀/c² = 6.9693e-10`, which also
/// includes the Earth's oblateness. That is outside a point-mass model; to reproduce it, add `W₀/c²` to
/// [`circular_orbit_rate_offset`] instead of calling this function.
///
/// First post-Newtonian order, point mass, circular orbit. Within a few ULP unless the terms nearly cancel
/// (an orbit close to `3R/2`).
///
/// Reference: Ashby 2003.
pub fn orbit_vs_surface_rate_offset(
    mass: f64,
    orbit_radius: f64,
    surface_radius: f64,
    surface_speed: f64,
) -> f64 {
    if mass >= 0.0 && orbit_radius > 0.0 && surface_radius > 0.0 && surface_speed >= 0.0 {
        let surface =
            static_rate_offset(mass, surface_radius) - 0.5 * surface_speed * surface_speed / C2;
        circular_orbit_rate_offset(mass, orbit_radius) - surface
    } else {
        f64::NAN
    }
}
