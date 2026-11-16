//! The Kerr metric: the field of a rotating, uncharged black hole.
//!
//! The hole is described by its mass `M` and its dimensionless spin `a* = J c / (G M²)`, between −1 and 1,
//! where `J` is its angular momentum. `a* = 0` is the Schwarzschild hole and `|a*| = 1` is the extremal hole.
//! Radii are Boyer-Lindquist radii `r` and `θ` is the polar angle measured from the spin axis, so the
//! equatorial plane is `θ = π/2`. With `m = GM/c²` and `a = a* m` the metric functions are
//! `Σ = r² + a² cos²θ`, `Δ = r² − 2 m r + a²` and `A = (r² + a²)² − a² Δ sin²θ`.
//!
//! **Sign of the spin.** The sign of `a*` is the sense in which the hole turns about its axis. It changes only
//! the sign of the angular velocities [`horizon_angular_velocity`] and [`zamo_angular_velocity`]; every other
//! function depends on `|a*|`.
//!
//! **Sense of an orbit.** Orbit functions take an [`Orbit`]: [`Orbit::Prograde`] turns the same way as the
//! hole and [`Orbit::Retrograde`] the opposite way, whatever the sign of `a*`. Their angular velocities and
//! angular momenta are returned as magnitudes, measured in the direction of the orbital motion. At `a* = 0`
//! the two senses coincide and every function reduces to its [`crate::schwarzschild`] counterpart.
//!
//! Masses are in kg and radii in metres. Functions return NaN outside the stated domain, in particular for
//! `|a*| > 1`.
//!
//! **Accuracy near a horizon and near the extremal limit.** A clock rate just outside the horizon or the
//! static limit depends on the small difference between `r` and that radius; the difference is formed
//! exactly, but the radius itself is known only to its last bit, so the relative error is a few ULP times
//! `r/(r − r₊)`. In the same way a spin close to 1 is resolved only to `1.1e-16`, while the prograde ISCO
//! depends on `(1 − a*)^(1/3)`. Each function says where this matters.
//!
//! References: R. P. Kerr, *Phys. Rev. Lett.* 11 (1963) 237; R. H. Boyer and R. W. Lindquist, *J. Math. Phys.*
//! 8 (1967) 265; J. M. Bardeen, W. H. Press and S. A. Teukolsky, *Astrophys. J.* 178 (1972) 347 (cited below
//! as BPT): section II for the metric and eqs. (2.12) to (2.21) for circular orbits, section III for the
//! zero-angular-momentum observers ("locally nonrotating frames"); C. W. Misner, K. S. Thorne and
//! J. A. Wheeler, *Gravitation* (1973), chapter 33 (MTW); D. Christodoulou, *Phys. Rev. Lett.* 25 (1970) 1596
//! for the irreducible mass.

use crate::{gravitational_radius, one_minus_ratio, C2};
use libm::{acos, cbrt, cos, fabs, sin, sqrt};
use srkit_core::units::{C, G_NEWTON};

/// `2/(3√3)`, the factor in the angular momentum on the ISCO.
const TWO_OVER_3_ROOT_3: f64 = 0.384_900_179_459_750_5;

/// The sense of an equatorial orbit relative to the rotation of the hole.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orbit {
    /// The orbit turns the same way as the hole (co-rotating).
    Prograde,
    /// The orbit turns against the hole (counter-rotating).
    Retrograde,
}

impl Orbit {
    /// The spin as the orbit sees it: `+|a*|` for a prograde orbit, `−|a*|` for a retrograde one.
    fn effective_spin(self, spin: f64) -> f64 {
        match self {
            Self::Prograde => fabs(spin),
            Self::Retrograde => -fabs(spin),
        }
    }
}

/// The radius of the outer (event) horizon, `r₊ = (GM/c²)(1 + √(1 − a*²))`, in metres.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
///
/// The larger root of `Δ = 0`. `2GM/c²` at `a* = 0`, falling to `GM/c²` for the extremal hole. Within a few
/// ULP.
///
/// Reference: BPT, section II; MTW, section 33.2.
pub fn outer_horizon_radius(mass: f64, spin: f64) -> f64 {
    if mass >= 0.0 && spin_ok(spin) {
        gravitational_radius(mass) * (1.0 + root_one_minus_square(spin))
    } else {
        f64::NAN
    }
}

/// The radius of the inner (Cauchy) horizon, `r₋ = (GM/c²)(1 − √(1 − a*²))`, in metres.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
///
/// The smaller root of `Δ = 0`. 0 at `a* = 0`, rising to `GM/c²` for the extremal hole, where the two
/// horizons meet. Evaluated as `(GM/c²) a*²/(1 + √(1 − a*²))`, so it is within a few ULP for small spins
/// too, where it tends to `a*² GM/(2c²)`.
///
/// Reference: BPT, section II; MTW, section 33.2.
pub fn inner_horizon_radius(mass: f64, spin: f64) -> f64 {
    if mass >= 0.0 && spin_ok(spin) {
        gravitational_radius(mass) * one_minus_root_one_minus_square(spin)
    } else {
        f64::NAN
    }
}

/// The radius of the static limit, the outer edge of the ergosphere,
/// `r_E(θ) = (GM/c²)(1 + √(1 − a*² cos²θ))`, in metres.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `polar_angle`: `θ`, rad, measured from the spin axis; `π/2` is the equatorial plane.
///
/// Between the outer horizon and this surface nothing can stay at rest relative to distant observers. In the
/// equatorial plane the static limit is at `2GM/c²` for every spin; on the axis it touches the horizon.
/// Within a few ULP.
///
/// Reference: MTW, section 33.2 and Box 33.2.
pub fn ergosphere_radius(mass: f64, spin: f64, polar_angle: f64) -> f64 {
    if mass >= 0.0 && spin_ok(spin) && polar_angle.is_finite() {
        gravitational_radius(mass) * (1.0 + root_one_minus_square(spin * cos(polar_angle)))
    } else {
        f64::NAN
    }
}

/// The angular velocity of the horizon, `Ω_H = a* c³ / (2 G M (1 + √(1 − a*²)))`, in rad/s.
///
/// * `mass`: kg, `> 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
///
/// The rate at which the horizon turns as seen from far away, and the limit of the frame-dragging rate
/// [`zamo_angular_velocity`] at the horizon. It has the sign of the spin, is 0 at `a* = 0` and reaches
/// `c³/(2GM)` for the extremal hole (`1.0e5` rad/s for one solar mass). Within a few ULP.
///
/// Reference: BPT, section III; MTW, section 33.4.
pub fn horizon_angular_velocity(mass: f64, spin: f64) -> f64 {
    if mass > 0.0 && spin_ok(spin) {
        spin * C * C2 / (2.0 * G_NEWTON * mass * (1.0 + root_one_minus_square(spin)))
    } else {
        f64::NAN
    }
}

/// The irreducible mass, `M_irr = M √((1 + √(1 − a*²))/2)`, in kg.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
///
/// The mass of the Schwarzschild hole with the same horizon area, `A = 16π (G M_irr/c²)²`. No classical
/// process can lower it; what lies above it, `M − M_irr`, is rotational energy that can in principle be
/// taken out. Equal to `M` at `a* = 0` and `M/√2` for the extremal hole. Within a few ULP.
///
/// Reference: D. Christodoulou, *Phys. Rev. Lett.* 25 (1970) 1596; D. Christodoulou and R. Ruffini,
/// *Phys. Rev. D* 4 (1971) 3552; MTW, section 33.8.
pub fn irreducible_mass(mass: f64, spin: f64) -> f64 {
    if mass >= 0.0 && spin_ok(spin) {
        mass * sqrt(0.5 * (1.0 + root_one_minus_square(spin)))
    } else {
        f64::NAN
    }
}

/// The fraction of the hole's mass-energy that is extractable rotational energy,
/// `1 − M_irr/M = 1 − √((1 + √(1 − a*²))/2)`, dimensionless.
///
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
///
/// 0 at `a* = 0`, rising to `1 − 1/√2 = 29.29%` for the extremal hole. This is the most that the Penrose
/// process, superradiance or the Blandford-Znajek mechanism can take out in total; no real process is this
/// efficient. Evaluated without cancellation: within a few ULP for small spins too, where it tends to
/// `a*²/8`.
///
/// Reference: Christodoulou (1970); MTW, sections 33.7 and 33.8.
pub fn extractable_energy_fraction(spin: f64) -> f64 {
    if spin_ok(spin) {
        // With h = (1 + √(1 − a*²))/2 = (M_irr/M)²: 1 − √h = (1 − h)/(1 + √h), and
        // 1 − h = (1 − √(1 − a*²))/2.
        let h = 0.5 * (1.0 + root_one_minus_square(spin));
        0.5 * one_minus_root_one_minus_square(spin) / (1.0 + sqrt(h))
    } else {
        f64::NAN
    }
}

/// The clock rate of a static observer, `dτ/dt = √(1 − 2 m r/Σ)`, dimensionless.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `r`: Boyer-Lindquist radius of the observer, m, at or outside the static limit [`ergosphere_radius`]
///   and `> 0`; `+∞` is allowed.
/// * `polar_angle`: `θ`, rad, measured from the spin axis.
///
/// A static observer keeps fixed `r`, `θ` and `φ`, at rest relative to the distant stars. `τ` is that
/// observer's proper time and `t` the time of a static clock far away. The rate is 0 on the static limit;
/// inside the ergosphere no static observer exists and the result is NaN. At `a* = 0` this is
/// [`crate::schwarzschild::static_dtau_dt`].
///
/// Within a few ULP at twice the static-limit radius and beyond; closer in, a few ULP times
/// `r/(r − r_E)` (see the module note).
///
/// Reference: MTW, sections 33.2 and 33.4.
pub fn static_dtau_dt(mass: f64, spin: f64, r: f64, polar_angle: f64) -> f64 {
    if !(mass >= 0.0 && spin_ok(spin) && r > 0.0 && polar_angle.is_finite()) {
        return f64::NAN;
    }
    let m = gravitational_radius(mass);
    let ac = spin * cos(polar_angle);
    // Σ − 2 m r = (r − e₊)(r − e₋) with e± = m (1 ± √(1 − a*² cos²θ)).
    let outer = m * (1.0 + root_one_minus_square(ac));
    let inner = m * one_minus_root_one_minus_square(ac);
    if r < outer {
        return f64::NAN;
    }
    let u = m / r;
    let sigma = 1.0 + ac * ac * u * u;
    sqrt(one_minus_ratio(outer, r) * one_minus_ratio(inner, r) / sigma)
}

/// The frame-dragging angular velocity, `ω = −g_tφ/g_φφ = 2 m a r c / A`, in rad/s.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `r`: Boyer-Lindquist radius, m, `r ≥ r₊` ([`outer_horizon_radius`]) and `> 0`; `+∞` is allowed.
/// * `polar_angle`: `θ`, rad, measured from the spin axis.
///
/// The angular velocity `dφ/dt`, seen from far away, of a zero-angular-momentum observer (ZAMO) at that
/// place: an observer with no angular momentum is still carried round by the hole. It has the sign of the
/// spin, equals [`horizon_angular_velocity`] on the horizon and falls off as `2 G J/(c² r³)` far away (the
/// Lense-Thirring rate). Within a few ULP.
///
/// Reference: BPT, section III; J. M. Bardeen, *Astrophys. J.* 162 (1970) 71.
pub fn zamo_angular_velocity(mass: f64, spin: f64, r: f64, polar_angle: f64) -> f64 {
    match zamo_parts(mass, spin, r, polar_angle) {
        Some(parts) => C * 2.0 * spin * parts.u * parts.u / (r * parts.a_hat),
        None => f64::NAN,
    }
}

/// The clock rate of a zero-angular-momentum observer, `dτ/dt = √(Σ Δ / A)`, dimensionless.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `r`: Boyer-Lindquist radius, m, `r ≥ r₊` ([`outer_horizon_radius`]) and `> 0`; `+∞` is allowed.
/// * `polar_angle`: `θ`, rad, measured from the spin axis.
///
/// The observer keeps fixed `r` and `θ` and circles at the frame-dragging rate [`zamo_angular_velocity`].
/// Unlike a static observer, a ZAMO exists everywhere outside the horizon, the ergosphere included. `τ` is
/// its proper time and `t` the time of a static clock far away; the ratio is the lapse function. It is 0 on
/// the horizon and 1 at infinity, it is never below the static rate at the same place, and at `a* = 0` it
/// is [`crate::schwarzschild::static_dtau_dt`].
///
/// The redshift between two observers in the Kerr field depends on the angular momentum of the light, so it
/// is not a function of two radii alone. For light with zero angular momentum about the axis (along the
/// axis, for instance) `1 + z` is the ratio of the two ZAMO clock rates.
///
/// Within a few ULP at `2 r₊` and beyond; closer in, a few ULP times `r/(r − r₊)` (see the module note).
///
/// Reference: BPT, section III.
pub fn zamo_dtau_dt(mass: f64, spin: f64, r: f64, polar_angle: f64) -> f64 {
    match zamo_parts(mass, spin, r, polar_angle) {
        Some(parts) => sqrt(parts.sigma_hat * parts.delta_hat / parts.a_hat),
        None => f64::NAN,
    }
}

/// The radius of the equatorial circular photon orbit,
/// `r_ph = 2 (GM/c²) (1 + cos((2/3) arccos(∓|a*|)))`, in metres (upper sign prograde).
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
///
/// Circular orbits of massive bodies exist only outside this radius. `3GM/c²` at `a* = 0`; for the extremal
/// hole `GM/c²` prograde and `4GM/c²` retrograde. Those three values are returned exactly; otherwise within
/// a few ULP.
///
/// Reference: BPT, eq. (2.18).
pub fn photon_orbit_radius(mass: f64, spin: f64, orbit: Orbit) -> f64 {
    if mass >= 0.0 && spin_ok(spin) {
        gravitational_radius(mass) * photon_x(orbit.effective_spin(spin))
    } else {
        f64::NAN
    }
}

/// The radius of the marginally bound equatorial circular orbit,
/// `r_mb = (GM/c²)(2 ∓ |a*| + 2 √(1 ∓ |a*|))`, in metres (upper sign prograde).
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
///
/// The circular orbit whose energy equals the rest energy. `4GM/c²` at `a* = 0`; for the extremal hole
/// `GM/c²` prograde and `(3 + 2√2) GM/c²` retrograde. Within a few ULP.
///
/// Reference: BPT, eq. (2.19).
pub fn marginally_bound_radius(mass: f64, spin: f64, orbit: Orbit) -> f64 {
    if mass >= 0.0 && spin_ok(spin) {
        let a = orbit.effective_spin(spin);
        gravitational_radius(mass) * (2.0 - a + 2.0 * sqrt(1.0 - a))
    } else {
        f64::NAN
    }
}

/// The radius of the innermost stable circular orbit (ISCO), in metres.
///
/// With `Z₁ = 1 + (1 − a*²)^(1/3) ((1 + |a*|)^(1/3) + (1 − |a*|)^(1/3))` and `Z₂ = √(3 a*² + Z₁²)`,
/// `r_isco = (GM/c²) (3 + Z₂ ∓ √((3 − Z₁)(3 + Z₁ + 2 Z₂)))` (upper sign prograde).
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
///
/// `6GM/c²` at `a* = 0`. For the extremal hole it is `GM/c²` prograde and `9GM/c²` retrograde, and those two
/// limits are returned exactly. The radius falls monotonically as the spin grows for a prograde orbit and
/// rises for a retrograde one.
///
/// The formula is rearranged so that neither `3 − Z₁` (which vanishes as `a*²` at small spin) nor the
/// prograde difference (which vanishes at the extremal limit) is formed by subtraction. Within a few ULP for
/// every spin.
///
/// Reference: BPT, eq. (2.21).
pub fn isco_radius(mass: f64, spin: f64, orbit: Orbit) -> f64 {
    if mass >= 0.0 && spin_ok(spin) {
        gravitational_radius(mass) * isco_x(orbit.effective_spin(spin))
    } else {
        f64::NAN
    }
}

/// The specific energy on the ISCO, `E = √(1 − 2GM/(3 r_isco c²))`, dimensionless.
///
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
///
/// The conserved energy per unit rest energy, `E/(μc²)`, of a body on the ISCO, as measured from far away.
/// `√(8/9)` at `a* = 0`; for the extremal hole `1/√3` prograde and `√(25/27)` retrograde. Within a few ULP.
///
/// Reference: BPT, eqs. (2.12) and (2.21); the closed form in terms of `r_isco` follows from them.
pub fn isco_specific_energy(spin: f64, orbit: Orbit) -> f64 {
    if spin_ok(spin) {
        let x = isco_x(orbit.effective_spin(spin));
        sqrt((3.0 * x - 2.0) / (3.0 * x))
    } else {
        f64::NAN
    }
}

/// The specific angular momentum on the ISCO,
/// `L = (2/(3√3)) (GM/c) (1 + 2 √(3 r_isco c²/(GM) − 2))`, in m²/s.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
///
/// The conserved angular momentum per unit rest mass about the spin axis, as a magnitude in the direction of
/// the orbital motion. In units of `GM/c` it is `2√3` at `a* = 0`; for the extremal hole `2/√3` prograde and
/// `22/(3√3)` retrograde. Within a few ULP.
///
/// Reference: BPT, eqs. (2.13) and (2.21); the closed form in terms of `r_isco` follows from them.
pub fn isco_specific_angular_momentum(mass: f64, spin: f64, orbit: Orbit) -> f64 {
    if mass >= 0.0 && spin_ok(spin) {
        let x = isco_x(orbit.effective_spin(spin));
        G_NEWTON * mass / C * TWO_OVER_3_ROOT_3 * (1.0 + 2.0 * sqrt(3.0 * x - 2.0))
    } else {
        f64::NAN
    }
}

/// The time dilation on the ISCO, `dt/dτ`, dimensionless.
///
/// With `x = r_isco c²/(GM)` and `s = √(3x − 2)`, `dt/dτ = E (s² − s + 6)/(s (s − 1))`, where `E` is
/// [`isco_specific_energy`]. This is BPT's `u^t = (x^(3/2) + a)/(x^(3/4) √(x^(3/2) − 3 x^(1/2) + 2a))` on the
/// ISCO, with `a` eliminated so that the vanishing denominator is not formed by subtraction.
///
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
///
/// `τ` is the proper time of a body on the ISCO and `t` the time of a static clock far away. For a prograde
/// orbit it is `√2` at `a* = 0`, 2.70 at `a* = 0.9`, 10.8 at `a* = 0.998`, 145 at `a* = 1 − 1e-6`, and `+∞`
/// for the extremal hole, growing as `(1 − a*)^(−1/3)`. For a retrograde orbit it falls from `√2` to
/// `13/(6√3) = 1.251`.
///
/// Accuracy: relative error below `4e-15` for `|a*| ≤ 1 − 1e-6`. Closer to the extremal prograde limit the
/// error rises as `2e-16 GM/(c² r_isco − GM)`, to about `2e-11` at the last representable spins below 1. That
/// is far below the effect of the spin's own resolution there: `1 − 1e-14` is not a binary64 number, and the
/// nearest one, `1 − 9.992e-15`, gives a result 0.03% higher.
///
/// Reference: BPT, eqs. (2.12), (2.16) and (2.21).
pub fn isco_dt_dtau(spin: f64, orbit: Orbit) -> f64 {
    if spin_ok(spin) {
        let x = isco_x(orbit.effective_spin(spin));
        let s = sqrt(3.0 * x - 2.0);
        // s − 1 = 3 (x − 1)/(s + 1); x − 1 is exact for x ≤ 2.
        let s_minus_1 = 3.0 * (x - 1.0) / (s + 1.0);
        sqrt((3.0 * x - 2.0) / (3.0 * x)) * (3.0 * x + 4.0 - s) / (s * s_minus_1)
    } else {
        f64::NAN
    }
}

/// The accretion efficiency of a thin disc, `η = 1 − E_isco`, dimensionless.
///
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense in which the disc orbits.
///
/// The fraction of its rest energy that matter radiates while spiralling slowly inwards through circular
/// orbits to the ISCO, after which it falls in without radiating more (the Novikov-Thorne disc).
/// `1 − √(8/9) = 5.72%` at `a* = 0`, 15.6% at `a* = 0.9`, 32.1% at `a* = 0.998` and `1 − 1/√3 = 42.26%` for
/// the extremal prograde disc; `1 − √(25/27) = 3.77%` for the extremal retrograde disc. The radiation that
/// the hole itself captures is not subtracted; with it, accretion cannot spin a hole beyond `a* ≈ 0.998`
/// (Thorne 1974). Within a few ULP.
///
/// Reference: I. D. Novikov and K. S. Thorne, in *Black Holes* (Les Houches 1972), eds. C. DeWitt and
/// B. S. DeWitt (1973); BPT, eqs. (2.12) and (2.21); K. S. Thorne, *Astrophys. J.* 191 (1974) 507.
pub fn accretion_efficiency(spin: f64, orbit: Orbit) -> f64 {
    if spin_ok(spin) {
        let x = isco_x(orbit.effective_spin(spin));
        // 1 − √(1 − y) = y/(1 + √(1 − y)) with y = 2/(3x).
        let energy = sqrt((3.0 * x - 2.0) / (3.0 * x));
        2.0 / (3.0 * x) / (1.0 + energy)
    } else {
        f64::NAN
    }
}

/// The time dilation on an equatorial circular orbit, dimensionless.
///
/// With `w = √(GM/(r c²))` and `a = ±|a*|` (plus for prograde),
/// `dt/dτ = (1 + a w³)/√(1 − 3 w² + 2 a w³)`.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
/// * `r`: Boyer-Lindquist radius of the orbit, m, outside the photon orbit [`photon_orbit_radius`]; `+∞`
///   is allowed.
///
/// `τ` is the proper time of the freely falling orbiter and `t` the time of a static clock far away. The
/// orbit is a geodesic for every radius outside the photon orbit, but it is stable only at or outside
/// [`isco_radius`]. The result tends to `+∞` at the photon orbit and is NaN inside it. At `a* = 0` this is
/// [`crate::schwarzschild::circular_orbit_dt_dtau`].
///
/// Accuracy: the denominator is a difference of nearly equal terms near the photon orbit, and on the ISCO of
/// a nearly extremal hole; the relative error is a few ULP divided by `1 − 3 w² + 2 a w³`. On the ISCO use
/// [`isco_dt_dtau`], which does not lose accuracy.
///
/// Reference: BPT, eqs. (2.12) to (2.16): the time component `u^t` of the four-velocity of the circular
/// geodesic.
pub fn circular_orbit_dt_dtau(mass: f64, spin: f64, orbit: Orbit, r: f64) -> f64 {
    match circular_parts(mass, spin, orbit, r) {
        Some(c) => (1.0 + c.a * c.w3) / c.root,
        None => f64::NAN,
    }
}

/// The specific energy on an equatorial circular orbit, dimensionless.
///
/// With `w = √(GM/(r c²))` and `a = ±|a*|` (plus for prograde),
/// `E = (1 − 2 w² + a w³)/√(1 − 3 w² + 2 a w³)`.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
/// * `r`: Boyer-Lindquist radius of the orbit, m, outside the photon orbit [`photon_orbit_radius`]; `+∞`
///   is allowed.
///
/// The conserved energy per unit rest energy as measured from far away; `1 − E` is the binding energy
/// released in reaching the orbit. It is below 1 outside [`marginally_bound_radius`], has its minimum on
/// the ISCO and tends to `+∞` at the photon orbit. Accuracy as for [`circular_orbit_dt_dtau`].
///
/// Reference: BPT, eq. (2.12).
pub fn circular_orbit_specific_energy(mass: f64, spin: f64, orbit: Orbit, r: f64) -> f64 {
    match circular_parts(mass, spin, orbit, r) {
        Some(c) => (1.0 - 2.0 * c.w2 + c.a * c.w3) / c.root,
        None => f64::NAN,
    }
}

/// The specific angular momentum on an equatorial circular orbit, in m²/s.
///
/// With `w = √(GM/(r c²))` and `a = ±|a*|` (plus for prograde),
/// `L = √(G M r) (1 − 2 a w³ + a² w⁴)/√(1 − 3 w² + 2 a w³)`.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
/// * `r`: Boyer-Lindquist radius of the orbit, m, outside the photon orbit [`photon_orbit_radius`].
///
/// The conserved angular momentum per unit rest mass about the spin axis, as a magnitude in the direction of
/// the orbital motion. It tends to the Newtonian `√(G M r)` far away and has its minimum on the ISCO.
/// Accuracy as for [`circular_orbit_dt_dtau`].
///
/// Reference: BPT, eq. (2.13).
pub fn circular_orbit_specific_angular_momentum(mass: f64, spin: f64, orbit: Orbit, r: f64) -> f64 {
    match circular_parts(mass, spin, orbit, r) {
        Some(c) => {
            sqrt(G_NEWTON * mass * r) * (1.0 - 2.0 * c.a * c.w3 + c.a * c.a * c.w2 * c.w2) / c.root
        }
        None => f64::NAN,
    }
}

/// The angular velocity on an equatorial circular orbit, in rad/s.
///
/// With `w = √(GM/(r c²))` and `a = ±|a*|` (plus for prograde), `Ω = √(GM/r³)/(1 + a w³)`.
///
/// * `mass`: kg, `≥ 0`.
/// * `spin`: `a*`, dimensionless, in `[−1, 1]`.
/// * `orbit`: the sense of the orbit.
/// * `r`: Boyer-Lindquist radius of the orbit, m, at or outside the photon orbit [`photon_orbit_radius`].
///
/// `|dφ/dt|` as seen from far away, a magnitude in the direction of the orbital motion. A prograde orbit is
/// slower than the Kepler rate `√(GM/r³)` and a retrograde one faster. The period seen from far away is
/// `2π/Ω`; the orbiter's own clock measures it shorter by [`circular_orbit_dt_dtau`]. Within a few ULP.
///
/// Reference: BPT, eq. (2.16).
pub fn circular_orbit_angular_velocity(mass: f64, spin: f64, orbit: Orbit, r: f64) -> f64 {
    if !(mass >= 0.0 && spin_ok(spin) && r > 0.0) {
        return f64::NAN;
    }
    let a = orbit.effective_spin(spin);
    let m = gravitational_radius(mass);
    if r < m * photon_x(a) {
        return f64::NAN;
    }
    let w2 = m / r;
    sqrt(G_NEWTON * mass / r) / r / (1.0 + a * w2 * sqrt(w2))
}

/// `|a*| ≤ 1`, false for NaN.
fn spin_ok(spin: f64) -> bool {
    (-1.0..=1.0).contains(&spin)
}

/// `√(1 − x²)` for `|x| ≤ 1`, with `1 − x²` formed as `(1 − x)(1 + x)` so it stays accurate near `|x| = 1`.
fn root_one_minus_square(x: f64) -> f64 {
    sqrt((1.0 - x) * (1.0 + x))
}

/// `1 − √(1 − x²)` for `|x| ≤ 1`, as `x²/(1 + √(1 − x²))` so it stays accurate for small `x` and near
/// `|x| = 1`.
fn one_minus_root_one_minus_square(x: f64) -> f64 {
    x * x / (1.0 + root_one_minus_square(x))
}

/// Photon orbit radius in units of `GM/c²` for an effective spin in `[−1, 1]` (positive: prograde).
///
/// BPT eq. (2.18). The three spins with an exact answer return it exactly, so that the Schwarzschild and
/// extremal limits hold to the last bit.
fn photon_x(a: f64) -> f64 {
    if a == 0.0 {
        3.0
    } else if a == 1.0 {
        1.0
    } else if a == -1.0 {
        4.0
    } else {
        2.0 * (1.0 + cos(2.0 / 3.0 * acos(-a)))
    }
}

/// ISCO radius in units of `GM/c²` for an effective spin in `[−1, 1]` (positive: prograde).
///
/// BPT eq. (2.21), `x = 3 + Z₂ ∓ √((3 − Z₁)(3 + Z₁ + 2 Z₂))`, rearranged in two places.
///
/// With `p = (1 + a)^(1/3)`, `q = (1 − a)^(1/3)` and `s = p + q`, so that `Z₁ = 1 + p q s` and `p³ + q³ = 2`:
/// `3 − Z₁ = (8 − s³)/3` and `2 − s = (1 − p) + (1 − q)`, which factorises to
/// `3 − Z₁ = 2 a² (1 + s)(4 + 2 s + s²) / (3 (p² + p q + q²)(1 + p + p²)(1 + q + q²))`.
/// Every factor is positive, so nothing cancels at small spin.
///
/// For the prograde root, `(3 + Z₂)² − (3 − Z₁)(3 + Z₁ + 2 Z₂) = (Z₁ + Z₂)²`, so
/// `x = (Z₁ + Z₂)²/(3 + Z₂ + √(…))`, which does not cancel near the extremal limit.
fn isco_x(a_eff: f64) -> f64 {
    let a = fabs(a_eff);
    // The extremal limits are exact; returning them as such keeps `x ≥ 1`, which `isco_dt_dtau` relies on.
    if a == 1.0 {
        return if a_eff > 0.0 { 1.0 } else { 9.0 };
    }
    let p = cbrt(1.0 + a);
    let q = cbrt(1.0 - a);
    let s = p + q;
    let pq = p * q;
    let z1 = 1.0 + pq * s;
    let z2 = sqrt(3.0 * a * a + z1 * z1);
    let three_minus_z1_over_a2 = 2.0 * (1.0 + s) * (4.0 + 2.0 * s + s * s)
        / (3.0 * (p * p + pq + q * q) * (1.0 + p + p * p) * (1.0 + q + q * q));
    let root = a * sqrt(three_minus_z1_over_a2 * (3.0 + z1 + 2.0 * z2));
    if a_eff >= 0.0 {
        (z1 + z2) * (z1 + z2) / (3.0 + z2 + root)
    } else {
        3.0 + z2 + root
    }
}

/// The dimensionless pieces of the metric at a point outside the outer horizon.
struct ZamoParts {
    /// `GM/(r c²)`.
    u: f64,
    /// `Σ/r²`.
    sigma_hat: f64,
    /// `Δ/r²`.
    delta_hat: f64,
    /// `A/r⁴`.
    a_hat: f64,
}

/// Metric pieces for the ZAMO functions, or `None` outside their domain.
fn zamo_parts(mass: f64, spin: f64, r: f64, polar_angle: f64) -> Option<ZamoParts> {
    if !(mass >= 0.0 && spin_ok(spin) && r > 0.0 && polar_angle.is_finite()) {
        return None;
    }
    let m = gravitational_radius(mass);
    let outer = m * (1.0 + root_one_minus_square(spin));
    let inner = m * one_minus_root_one_minus_square(spin);
    if r < outer {
        return None;
    }
    let u = m / r;
    let au = spin * u;
    let (sine, cosine) = (sin(polar_angle), cos(polar_angle));
    // Δ = (r − r₊)(r − r₋), formed from exact differences near the horizon.
    let delta_hat = one_minus_ratio(outer, r) * one_minus_ratio(inner, r);
    let sigma_hat = 1.0 + au * au * cosine * cosine;
    let one_plus = 1.0 + au * au;
    let a_hat = one_plus * one_plus - au * au * delta_hat * sine * sine;
    Some(ZamoParts {
        u,
        sigma_hat,
        delta_hat,
        a_hat,
    })
}

/// The pieces shared by the circular-orbit functions.
struct CircularParts {
    /// Effective spin: `+|a*|` prograde, `−|a*|` retrograde.
    a: f64,
    /// `w² = GM/(r c²)`.
    w2: f64,
    /// `w³`.
    w3: f64,
    /// `√(1 − 3 w² + 2 a w³)`.
    root: f64,
}

/// Pieces for the circular-orbit functions, or `None` outside their domain.
///
/// The radius is tested against the photon orbit explicitly. The sign of `1 − 3 w² + 2 a w³` is not enough:
/// the expression is positive again deep inside the hole, and for a nearly extremal hole it is smaller than
/// its own rounding error everywhere between the horizon and the photon orbit.
fn circular_parts(mass: f64, spin: f64, orbit: Orbit, r: f64) -> Option<CircularParts> {
    if !(mass >= 0.0 && spin_ok(spin) && r > 0.0) {
        return None;
    }
    let a = orbit.effective_spin(spin);
    let m = gravitational_radius(mass);
    if r < m * photon_x(a) {
        return None;
    }
    let w2 = m / r;
    let w3 = w2 * sqrt(w2);
    let root = sqrt(1.0 - 3.0 * w2 + 2.0 * a * w3);
    Some(CircularParts { a, w2, w3, root })
}
