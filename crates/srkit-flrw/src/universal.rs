//! The dimensionless background shared by every flat matter + Λ cosmology.
//!
//! With `a_eq = (Ω_m/Ω_Λ)^(1/3)` (the scale factor at which matter and Λ have equal density),
//! `H_Λ = H0 √Ω_Λ` and `θ = (3/2) H_Λ t`, the scale factor is `a = a_eq sinh^(2/3) θ`, and in the variable
//!
//! ```text
//! v = ln(a / a_eq) = (2/3) ln sinh θ
//! ```
//!
//! every background integral is a universal function of `v` times a power of `a_eq` and `H_Λ`:
//!
//! ```text
//! η            = eta_hat(v) / (a_eq H_Λ)        conformal time
//! η_∞ − η      = rem_hat(v) / (a_eq H_Λ)        conformal time left
//! ∫₀ᵗ a dt     = (a_eq / H_Λ) e^v drift_hat(v)
//! ```
//!
//! Writing `x = e^(3v) = (Ω_Λ/Ω_m) a³`, the closed forms are (see [`crate::series`] for the Pfaff step):
//!
//! ```text
//! eta_hat(v)   = 2 e^(v/2) (1 + x)^(-1/2) 2F1(1, 1/2; 7/6; x/(1 + x))            used for v ≤ 0
//! rem_hat(v)   = e^(-v) (1 + 1/x)^(-1/2) 2F1(1, 1/2; 4/3; 1/(1 + x))             used for v ≥ 0
//! eta_hat(v) + rem_hat(v) = Γ(1/6) Γ(1/3) / (3 √π)
//! drift_hat(v) = (2/5) (x/(1 + x))^(1/2) 2F1(1, 1/2; 11/6; x/(1 + x))            used for v ≤ 0
//!              = (x/(1 + x))^(1/2) 2F1(1, 1/2; 2/3; 1/(1 + x)) − K e^(-v)        used for v ≥ 0
//! K            = −(2/5) Γ(11/6) Γ(−1/3) / √π
//! ```
//!
//! In terms of the Hubble rate these read `η = (2/(aH)) 2F1(1, 1/2; 7/6; Ω_Λ(a))` and
//! `η_∞ − η = (1/(aH)) 2F1(1, 1/2; 4/3; Ω_m(a))`: the matter-era and de Sitter results times a correction that
//! tends to 1 in the respective era.
//!
//! Working in `v` keeps everything finite far beyond the point where `a` overflows binary64, and it makes
//! `ln η` and `ln(η_∞ − η)` concave in `v`, which is what lets the inverse functions use Newton's method with a
//! guarantee of convergence.

use crate::series::{f_11_6, f_2_3, f_4_3, f_7_6};
use core::f64::consts::LN_2;
use libm::{asinh, exp, fabs, log, log1p, sinh, sqrt};

/// `Γ(1/6) Γ(1/3) / (3 √π)`: the total conformal time in units of `1/(a_eq H_Λ)`.
pub(crate) const ETA_HAT_TOTAL: f64 = 2.804_364_210_650_908_4;

/// `rem_hat(0) = 2F1(1, 1/2; 4/3; 1/2) / √2`: the conformal time left at matter-Λ equality.
const REM_HAT_AT_EQUALITY: f64 = 0.909_604_242_638_895_7;

/// `eta_hat(0) = √2 2F1(1, 1/2; 7/6; 1/2)`: the conformal time elapsed at matter-Λ equality.
const ETA_HAT_AT_EQUALITY: f64 = 1.894_759_968_012_012_9;

/// `−(2/5) Γ(11/6) Γ(−1/3) / √π`: the constant in the late form of `∫ a dt`.
const DRIFT_CONSTANT: f64 = 0.862_369_853_076_596_9;

/// Above this `θ` (and this `v`) the decaying exponential is below one part in 2^57 and is dropped.
const LARGE: f64 = 20.0;

/// Newton steps allowed in the inverse functions. Convergence is quadratic from a one-sided start and takes
/// three to seven steps; the cap bounds the work for subnormal arguments, where rounding can keep the
/// iteration from settling.
const MAX_NEWTON_STEPS: usize = 40;

/// `v = (2/3) ln sinh θ` for `θ ≥ 0`.
pub(crate) fn v_of_theta(theta: f64) -> f64 {
    const TWO_THIRDS: f64 = 2.0 / 3.0;
    if theta < LARGE {
        TWO_THIRDS * log(sinh(theta))
    } else {
        TWO_THIRDS * (theta - LN_2)
    }
}

/// `θ = asinh(e^(3v/2))`, the inverse of [`v_of_theta`].
pub(crate) fn theta_of_v(v: f64) -> f64 {
    if v < LARGE {
        asinh(exp(1.5 * v))
    } else {
        1.5 * v + LN_2
    }
}

/// `(eta_hat(v), F)` from the early form, where `F = 2F1(1, 1/2; 7/6; x/(1 + x))`. Needs `v ≤ 0`.
fn early(v: f64) -> (f64, f64) {
    let x = exp(3.0 * v);
    let f = f_7_6(x / (1.0 + x));
    (2.0 * exp(0.5 * v) / sqrt(1.0 + x) * f, f)
}

/// `(rem_hat(v), F)` from the late form, where `F = 2F1(1, 1/2; 4/3; 1/(1 + x))`. Needs `v ≥ 0`.
fn late(v: f64) -> (f64, f64) {
    let w = exp(-3.0 * v);
    let f = f_4_3(w / (1.0 + w));
    (exp(-v) / sqrt(1.0 + w) * f, f)
}

/// Conformal time in units of `1/(a_eq H_Λ)`.
pub(crate) fn eta_hat(v: f64) -> f64 {
    if v <= 0.0 {
        early(v).0
    } else {
        ETA_HAT_TOTAL - late(v).0
    }
}

/// Conformal time left before `η_∞`, in units of `1/(a_eq H_Λ)`.
pub(crate) fn rem_hat(v: f64) -> f64 {
    if v <= 0.0 {
        ETA_HAT_TOTAL - early(v).0
    } else {
        late(v).0
    }
}

/// `eta_hat(v2) − eta_hat(v1)`, taking the difference inside one closed form where both epochs allow it.
pub(crate) fn eta_hat_between(v1: f64, v2: f64) -> f64 {
    if v1 <= 0.0 && v2 <= 0.0 {
        early(v2).0 - early(v1).0
    } else if v1 > 0.0 && v2 > 0.0 {
        late(v1).0 - late(v2).0
    } else {
        eta_hat(v2) - eta_hat(v1)
    }
}

/// `e^v rem_hat(v)`: the proper distance to the event horizon in units of `1/H_Λ`. Finite for every `v`
/// and equal to `(1 + 1/x)^(-1/2) 2F1(1, 1/2; 4/3; 1/(1 + x))` in the late form.
pub(crate) fn proper_rem_hat(v: f64) -> f64 {
    if v <= 0.0 {
        exp(v) * (ETA_HAT_TOTAL - early(v).0)
    } else {
        let w = exp(-3.0 * v);
        f_4_3(w / (1.0 + w)) / sqrt(1.0 + w)
    }
}

/// `H_Λ ∫₀ᵗ a dt / a`: the peculiar momentum per unit proper acceleration, in units of `1/H_Λ`, that a
/// rocket thrusting since the Big Bang would hold. Tends to 1 as `v → ∞`.
pub(crate) fn drift_hat(v: f64) -> f64 {
    if v <= 0.0 {
        let x = exp(3.0 * v);
        let y = x / (1.0 + x);
        0.4 * sqrt(y) * f_11_6(y)
    } else {
        let w = exp(-3.0 * v);
        f_2_3(w / (1.0 + w)) / sqrt(1.0 + w) - DRIFT_CONSTANT * exp(-v)
    }
}

/// Solve `eta_hat(v) = e` for `v`, given `0 < e < ETA_HAT_TOTAL`.
pub(crate) fn v_of_eta_hat(e: f64) -> f64 {
    if e <= ETA_HAT_AT_EQUALITY {
        solve_early(e)
    } else {
        solve_late(ETA_HAT_TOTAL - e)
    }
}

/// Solve `rem_hat(v) = r` for `v`, given `0 < r < ETA_HAT_TOTAL`.
pub(crate) fn v_of_rem_hat(r: f64) -> f64 {
    if r <= REM_HAT_AT_EQUALITY {
        solve_late(r)
    } else {
        solve_early(ETA_HAT_TOTAL - r)
    }
}

/// Newton's method on `ln eta_hat(v) − ln e`, for `0 < e ≤ eta_hat(0)`.
///
/// `d ln eta_hat / dv = 1/(2F)` falls from 1/2 in the matter era to 0.37 at equality, so the function is
/// increasing and concave. The matter-era guess `2 e^(v/2) = e` is never above the root, and from there every
/// Newton step moves up without passing the root.
fn solve_early(e: f64) -> f64 {
    let mut v = 2.0 * log(0.5 * e);
    for _ in 0..MAX_NEWTON_STEPS {
        let (value, f) = early(v);
        let step = 2.0 * f * log1p((value - e) / e);
        v = (v - step).min(0.0);
        if converged(step, v) {
            break;
        }
    }
    v
}

/// Newton's method on `ln rem_hat(v) − ln r`, for `0 < r ≤ rem_hat(0)`.
///
/// `d ln rem_hat / dv = −1/F` falls from −0.78 at equality towards −1, so the function is decreasing and
/// concave. The de Sitter guess `e^(-v) = r` is never below the root, and from there every Newton step moves
/// down without passing the root.
fn solve_late(r: f64) -> f64 {
    let mut v = -log(r);
    for _ in 0..MAX_NEWTON_STEPS {
        let (value, f) = late(v);
        let step = f * log1p((value - r) / r);
        v = (v + step).max(0.0);
        if converged(step, v) {
            break;
        }
    }
    v
}

/// True once a Newton step no longer moves `v` by more than rounding, or is not a number.
fn converged(step: f64, v: f64) -> bool {
    step.is_nan() || fabs(step) <= 4.0 * f64::EPSILON * (1.0 + fabs(v))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs()
    }

    #[test]
    fn the_two_forms_meet_at_equality() {
        assert!(close(early(0.0).0, ETA_HAT_AT_EQUALITY, 2.0 * f64::EPSILON));
        assert!(close(late(0.0).0, REM_HAT_AT_EQUALITY, 2.0 * f64::EPSILON));
        assert!(close(
            early(0.0).0 + late(0.0).0,
            ETA_HAT_TOTAL,
            2.0 * f64::EPSILON
        ));
        // drift_hat is continuous across v = 0
        let below = drift_hat(0.0);
        let above = drift_hat(1e-300);
        assert!(close(below, above, 8.0 * f64::EPSILON));
    }

    #[test]
    fn theta_and_v_are_inverse() {
        for &theta in &[1e-12, 1e-3, 0.3, 0.88, 1.0, 5.0, 19.9, 20.1, 300.0, 1e4] {
            let back = theta_of_v(v_of_theta(theta));
            assert!(close(back, theta, 8.0 * f64::EPSILON), "{theta} -> {back}");
        }
    }

    /// The inverse functions return within a few units in the last place of the forward functions, and never
    /// need more than a handful of Newton steps.
    #[test]
    fn newton_inverts_both_branches() {
        let mut v = -80.0;
        while v < 600.0 {
            let e = eta_hat(v);
            let r = rem_hat(v);
            if v <= 0.0 {
                let back = v_of_eta_hat(e);
                assert!((back - v).abs() <= 1e-14 * (1.0 + v.abs()), "{v} -> {back}");
            } else {
                let back = v_of_rem_hat(r);
                assert!((back - v).abs() <= 1e-14 * (1.0 + v.abs()), "{v} -> {back}");
            }
            v += 0.37;
        }
    }

    #[test]
    fn newton_needs_few_steps() {
        fn steps_late(r: f64) -> usize {
            let mut v = -log(r);
            for n in 0..MAX_NEWTON_STEPS {
                let (value, f) = late(v);
                let step = f * log1p((value - r) / r);
                v = (v + step).max(0.0);
                if converged(step, v) {
                    return n + 1;
                }
            }
            MAX_NEWTON_STEPS
        }
        fn steps_early(e: f64) -> usize {
            let mut v = 2.0 * log(0.5 * e);
            for n in 0..MAX_NEWTON_STEPS {
                let (value, f) = early(v);
                let step = 2.0 * f * log1p((value - e) / e);
                v = (v - step).min(0.0);
                if converged(step, v) {
                    return n + 1;
                }
            }
            MAX_NEWTON_STEPS
        }
        let mut worst = 0;
        let mut k = 0.0;
        while k < 1.0 {
            worst = worst.max(steps_late(REM_HAT_AT_EQUALITY * (1.0 - k)));
            worst = worst.max(steps_late(REM_HAT_AT_EQUALITY * exp(-700.0 * k)));
            worst = worst.max(steps_early(ETA_HAT_AT_EQUALITY * (1.0 - k)));
            worst = worst.max(steps_early(ETA_HAT_AT_EQUALITY * exp(-700.0 * k)));
            k += 1.0 / 1024.0;
        }
        assert!(worst <= 8, "worst case took {worst} Newton steps");
    }
}
