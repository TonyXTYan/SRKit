//! The Schwarzschild metric: the field outside a non-rotating, uncharged, spherical mass.
//!
//! In Schwarzschild coordinates `(t, r, θ, φ)` the line element is
//! `ds² = −(1 − r_s/r) c² dt² + dr²/(1 − r_s/r) + r² dΩ²` with `r_s = 2GM/c²`. `t` is the time of a static
//! clock far away and `r` is the areal radius. By Birkhoff's theorem the same metric holds outside any
//! spherical body, so these functions apply above the surface of a star or a planet as well as outside a
//! hole.
//!
//! Masses are in kg and radii in metres. Functions return NaN outside the stated domain.
//!
//! **Accuracy near a horizon.** A clock rate close to the horizon, or a time dilation close to the photon
//! sphere, depends on the small difference between `r` and that radius. The functions form the difference
//! exactly, but the horizon radius itself, `2GM/c²`, is only known to its last bit. The relative error is
//! therefore a few ULP times `r_s/(r − r_s)` (or `3GM/(c² r − 3GM)`), which is a few ULP anywhere outside
//! `2 r_s`. That is a property of asking the question in terms of `M` and `r`, not of the evaluation.
//!
//! References: K. Schwarzschild, *Sitzungsber. Preuss. Akad. Wiss.* (1916) 189; C. W. Misner, K. S. Thorne and
//! J. A. Wheeler, *Gravitation* (1973), chapters 25 and 31 (cited below as MTW); J. B. Hartle, *Gravity: An
//! Introduction to Einstein's General Relativity* (2003), chapter 9; J. M. Bardeen, W. H. Press and
//! S. A. Teukolsky, *Astrophys. J.* 178 (1972) 347 (BPT), whose Kerr formulas reduce to these at zero spin.

use crate::{gravitational_radius, one_minus_ratio};
use libm::sqrt;
use srkit_core::units::G_NEWTON;

/// The Schwarzschild radius `r_s = 2GM/c²`, in metres: the radius of the event horizon.
///
/// * `mass`: kg, `≥ 0`.
///
/// 2,953 m for the Sun, 8.87 mm for the Earth. Within a few ULP.
///
/// Reference: MTW, chapter 31; Hartle, section 9.1.
pub fn schwarzschild_radius(mass: f64) -> f64 {
    scaled_radius(mass, 2.0)
}

/// The radius of the photon sphere, `3GM/c²`, in metres.
///
/// * `mass`: kg, `≥ 0`.
///
/// Light can circle the hole here, unstably. No circular orbit of a massive body exists at or inside this
/// radius. Within a few ULP.
///
/// Reference: MTW, section 25.6; Hartle, section 9.4.
pub fn photon_sphere_radius(mass: f64) -> f64 {
    scaled_radius(mass, 3.0)
}

/// The radius of the marginally bound circular orbit, `4GM/c²`, in metres.
///
/// * `mass`: kg, `≥ 0`.
///
/// The circular orbit whose energy equals the rest energy. A body that falls from rest far away with
/// specific angular momentum `4GM/c` approaches it without ever turning back out. Circular orbits between
/// `4GM/c²` and `6GM/c²` are bound but unstable. Within a few ULP.
///
/// Reference: BPT, eq. (2.19) at zero spin.
pub fn marginally_bound_radius(mass: f64) -> f64 {
    scaled_radius(mass, 4.0)
}

/// The radius of the innermost stable circular orbit (ISCO), `6GM/c²`, in metres.
///
/// * `mass`: kg, `≥ 0`.
///
/// Circular orbits are stable only at or outside this radius. On it `dt/dτ = √2` and the specific energy is
/// `√(8/9)`, so matter spiralling in through a thin disc gives up `1 − √(8/9) = 5.72%` of its rest energy.
/// Within a few ULP.
///
/// Reference: MTW, section 25.5; Hartle, section 9.3; BPT, eq. (2.21) at zero spin.
pub fn isco_radius(mass: f64) -> f64 {
    scaled_radius(mass, 6.0)
}

/// The clock rate of a static observer, `dτ/dt = √(1 − r_s/r)`, dimensionless.
///
/// * `mass`: kg, `≥ 0`.
/// * `r`: areal radius of the observer, m, `r ≥ r_s` and `r > 0`; `+∞` is allowed.
///
/// `τ` is the proper time of an observer who hovers at fixed `r` and `t` is the time of a static clock far
/// away. The result is 1 at infinity and 0 on the horizon; inside the horizon nothing can hover and the
/// result is NaN.
///
/// Far from the mass the result is close to 1 and its difference from 1 is poorly resolved; use
/// [`crate::weak_field::static_rate_offset`] for that difference. Within a few ULP outside `2 r_s`; see the
/// module note for the approach to the horizon.
///
/// Reference: Hartle, section 9.2; MTW, chapter 25.
pub fn static_dtau_dt(mass: f64, r: f64) -> f64 {
    let rs = schwarzschild_radius(mass);
    if mass >= 0.0 && r > 0.0 && r >= rs {
        sqrt(one_minus_ratio(rs, r))
    } else {
        f64::NAN
    }
}

/// The gravitational redshift `z` between two static observers, dimensionless.
///
/// `1 + z = λ_observed/λ_emitted = √((1 − r_s/r_observer)/(1 − r_s/r_emitter))`.
///
/// * `mass`: kg, `≥ 0`.
/// * `r_emitter`: areal radius of the static emitter, m, `≥ r_s` and `> 0`; `+∞` is allowed.
/// * `r_observer`: areal radius of the static observer, m, `≥ r_s` and `> 0`; `+∞` is allowed.
///
/// Positive when the light climbs (`r_observer > r_emitter`), negative (a blueshift, `z > −1`) when it
/// falls, and `+∞` for an emitter on the horizon. It depends only on the two radii, not on the directions
/// of the two observers from the mass. The received frequency is `ν_emitted/(1 + z)`.
///
/// The result is computed as `z` itself, not as `(1 + z) − 1`, so a weak-field redshift such as the Sun's
/// `2.12e-6`, or the `1.1e-16` of a one-metre climb on the Earth, keeps full precision. Within a few ULP for
/// radii outside `2 r_s`; see the module note for the approach to the horizon.
///
/// Reference: Hartle, section 9.2; MTW, chapter 25.
pub fn gravitational_redshift(mass: f64, r_emitter: f64, r_observer: f64) -> f64 {
    let rs = schwarzschild_radius(mass);
    if !(mass >= 0.0 && r_emitter > 0.0 && r_emitter >= rs && r_observer > 0.0 && r_observer >= rs)
    {
        return f64::NAN;
    }
    // Lapse squared at each end, and their difference without cancellation:
    // (1 − r_s/r_o) − (1 − r_s/r_e) = r_s (r_o − r_e)/(r_e r_o).
    let lapse2_e = one_minus_ratio(rs, r_emitter);
    let lapse2_o = one_minus_ratio(rs, r_observer);
    let difference = if r_emitter.is_infinite() || r_observer.is_infinite() {
        rs / r_emitter - rs / r_observer
    } else {
        (rs / r_emitter) * ((r_observer - r_emitter) / r_observer)
    };
    // z = √(A_o/A_e) − 1 = (A_o − A_e)/(A_e + √(A_e A_o)), which is +∞ for an emitter on the horizon.
    difference / (lapse2_e + sqrt(lapse2_e * lapse2_o))
}

/// The time dilation on a circular orbit, `dt/dτ = 1/√(1 − 3GM/(r c²))`, dimensionless.
///
/// * `mass`: kg, `≥ 0`.
/// * `r`: areal radius of the orbit, m, `r ≥ 3GM/c²` and `r > 0`; `+∞` is allowed.
///
/// `τ` is the proper time of the freely falling orbiter and `t` the time of a static clock far away. The
/// factor combines the gravitational term with the time dilation of the orbital speed. It is `√2` on the
/// ISCO and `+∞` on the photon sphere; inside the photon sphere there is no circular orbit and the result
/// is NaN. Orbits between `3GM/c²` and `6GM/c²` exist but are unstable.
///
/// For the small offset of a distant orbit use [`crate::weak_field::circular_orbit_rate_offset`]. Within a
/// few ULP outside `6GM/c²`; see the module note for the approach to the photon sphere.
///
/// Reference: Hartle, section 9.3; BPT, eqs. (2.12) and (2.16) at zero spin.
pub fn circular_orbit_dt_dtau(mass: f64, r: f64) -> f64 {
    let photon = photon_sphere_radius(mass);
    if mass >= 0.0 && r > 0.0 && r >= photon {
        1.0 / sqrt(one_minus_ratio(photon, r))
    } else {
        f64::NAN
    }
}

/// The angular velocity on a circular orbit, `Ω = dφ/dt = √(GM/r³)`, in rad/s.
///
/// * `mass`: kg, `≥ 0`.
/// * `r`: areal radius of the orbit, m, `r ≥ 3GM/c²` and `r > 0`.
///
/// Measured in the time `t` of a static clock far away, for which Kepler's third law holds exactly in
/// Schwarzschild coordinates. Divide by `2π` for the orbital frequency in Hz; the period seen from far away
/// is `2π/Ω`, the same as [`crate::newtonian::circular_period`]. The orbiter's own clock measures a period
/// shorter by the factor [`circular_orbit_dt_dtau`]. At `r = 3GM/c²` this is the angular velocity of light on
/// the photon sphere.
///
/// Within a few ULP.
///
/// Reference: Hartle, section 9.3; BPT, eq. (2.16) at zero spin.
pub fn circular_orbit_angular_velocity(mass: f64, r: f64) -> f64 {
    if mass >= 0.0 && r > 0.0 && r >= photon_sphere_radius(mass) {
        sqrt(G_NEWTON * mass / r) / r
    } else {
        f64::NAN
    }
}

/// `k GM/c²` for a non-negative mass, NaN otherwise.
fn scaled_radius(mass: f64, k: f64) -> f64 {
    if mass >= 0.0 {
        k * gravitational_radius(mass)
    } else {
        f64::NAN
    }
}
