//! Hubble drag: how expansion wears down motion relative to the comoving frame.
//!
//! A body moving at peculiar velocity `β` (measured by the comoving observer it is passing) has proper velocity
//! `u = γβ = sinh(rapidity)`, its momentum per unit rest energy. On a geodesic that momentum redshifts like a
//! photon's, `u ∝ 1/a` (Weinberg, *Cosmology*, 2008, section 1.1). A rocket holding constant proper
//! acceleration `α` along its motion obeys
//!
//! ```text
//! du/dt = α − H u
//! ```
//!
//! which is linear, so `u(t) = [a₀ u₀ + α ∫ a dt] / a(t)`, and `u` saturates at `α/H_Λ` once Λ dominates
//! (`docs/engine/deep-time-engine.md` section 4). Motion and thrust lie along one line; `u` and `α` are signed
//! along it.

use crate::cosmology::Cosmology;
use crate::universal::drift_hat;
use core::f64::consts::LN_2;
use libm::{asinh, atanh, copysign, exp, expm1, fabs, log1p, sinh, sqrt};

/// `ln cosh x`, without overflow.
fn ln_cosh(x: f64) -> f64 {
    let x = fabs(x);
    x + log1p(exp(-2.0 * x)) - LN_2
}

impl Cosmology {
    /// Proper velocity `u = γβ` at `t_to` of a free particle that had proper velocity `u` at `t_from`.
    ///
    /// ```text
    /// u(t_to) = u · a(t_from) / a(t_to)
    /// ```
    ///
    /// Domain: both times `≥ 0`. Exact for geodesic motion in any FLRW universe; the scale-factor ratio is
    /// formed from logarithms and does not overflow. Relative error about `(1 + |ln ratio|)` units in the
    /// last place.
    pub fn coasting_proper_velocity(&self, u: f64, t_from: f64, t_to: f64) -> f64 {
        u * exp(self.v(t_from) - self.v(t_to))
    }

    /// Rapidity at `t_to` of a free particle that had rapidity `rapidity` at `t_from`:
    /// `asinh(sinh(rapidity) · a(t_from)/a(t_to))`.
    ///
    /// Domain: both times `≥ 0`, any rapidity (no overflow for large ones).
    pub fn coasting_rapidity(&self, rapidity: f64, t_from: f64, t_to: f64) -> f64 {
        let ln_ratio = self.v(t_from) - self.v(t_to);
        let phi = fabs(rapidity);
        let out = if phi < 20.0 {
            asinh(sinh(phi) * exp(ln_ratio))
        } else {
            // sinh φ = e^φ / 2 to rounding, and asinh(e^m / 2) = m to rounding once m > 20
            let m = phi + ln_ratio;
            if m > 20.0 {
                m
            } else {
                asinh(0.5 * exp(m))
            }
        };
        copysign(out, rapidity)
    }

    /// Lorentz factor at `t_to` of a free particle that had Lorentz factor `gamma` at `t_from`:
    /// `√(1 + (γ² − 1) (a(t_from)/a(t_to))²)`.
    ///
    /// Domain: both times `≥ 0`, `gamma ≥ 1` (NaN otherwise).
    pub fn coasting_gamma(&self, gamma: f64, t_from: f64, t_to: f64) -> f64 {
        if gamma.is_nan() || gamma < 1.0 {
            return f64::NAN;
        }
        let ratio = exp(self.v(t_from) - self.v(t_to));
        sqrt(1.0 + (gamma - 1.0) * (gamma + 1.0) * ratio * ratio)
    }

    /// `∫₀ᵗ a(t′) dt′`, in units of time: the third background function, beside `a(t)` and `η(t)`.
    ///
    /// ```text
    /// ∫₀ᵗ a dt = (2/5) (a/H) · 2F1(1, 1/2; 11/6; Ω_Λ(a))                                          x ≤ 1
    ///          = (a/H) · 2F1(1, 1/2; 2/3; Ω_m(a)) − K a_eq/H_Λ,   K = −(2/5) Γ(11/6) Γ(−1/3)/√π    x ≥ 1
    /// ```
    ///
    /// with `x = (Ω_Λ/Ω_m) a³` and `K` = 0.862 369 853 076 596 8. Both are
    /// `(2/5) a^(5/2) / (H0 √Ω_m) · 2F1(1/2, 5/6; 11/6; −x)` after a Pfaff transformation, the second through
    /// the connection formula for large argument (DLMF 15.8.2), one of whose two series ends after its first
    /// term. Domain: `t ≥ 0` (NaN otherwise). Overflows with `a`. Relative error below 1e-14 for
    /// `H_Λ t < 30`.
    pub fn scale_factor_integral(&self, t: f64) -> f64 {
        let v = self.v(t);
        self.scale_factor_at_equality() / self.h_lambda() * exp(v) * drift_hat(v)
    }

    /// Proper velocity `u = γβ` at `t_to` of a rocket that had proper velocity `u` at `t_from` and holds
    /// constant proper acceleration `alpha` (rate of change of rapidity per unit of time, signed along the
    /// motion) in between.
    ///
    /// ```text
    /// u(t_to) = [a(t_from) u + α ∫ a dt] / a(t_to),   the integral from t_from to t_to
    /// ```
    ///
    /// the exact solution of `du/dt = α − H u`. Domain: both times `≥ 0`. Nothing overflows for late times.
    /// The integral is a difference of two values of [`Cosmology::scale_factor_integral`], so a burn much
    /// shorter than the Hubble time loses digits: the error in the thrust term is about
    /// `1e-16 / (H Δt)` of it. For such burns the drag is negligible and flat-space kinematics apply.
    pub fn thrust_proper_velocity(&self, u: f64, alpha: f64, t_from: f64, t_to: f64) -> f64 {
        let v_from = self.v(t_from);
        let v_to = self.v(t_to);
        let ratio = exp(v_from - v_to);
        u * ratio + alpha / self.h_lambda() * (drift_hat(v_to) - drift_hat(v_from) * ratio)
    }

    /// The proper velocity `u = γβ` at which Hubble drag balances a constant proper acceleration `alpha` once Λ
    /// dominates: `α/H_Λ`. Dimensionless. For 1 g it is 1.81e10 (oracle row SR-20).
    pub fn drag_limited_proper_velocity(&self, alpha: f64) -> f64 {
        alpha / self.h_lambda()
    }

    /// The Lorentz factor at which Hubble drag balances a constant proper acceleration `alpha` once Λ
    /// dominates: `√(1 + (α/H_Λ)²)`. For 1 g it is 1.81e10 (oracle row SR-20).
    pub fn drag_limited_gamma(&self, alpha: f64) -> f64 {
        let u = self.drag_limited_proper_velocity(alpha);
        sqrt(1.0 + u * u)
    }

    /// Proper velocity after thrusting for cosmic time `dt` in the de Sitter limit (`H = H_Λ` throughout),
    /// starting from proper velocity `u`.
    ///
    /// ```text
    /// u(dt) = U + (u − U) e^(−H_Λ dt),   U = α/H_Λ
    /// ```
    ///
    /// Domain: `dt ≥ 0`. This is [`Cosmology::thrust_proper_velocity`] with matter neglected; it is accurate
    /// for burns of any length, including short ones.
    pub fn de_sitter_thrust_proper_velocity(&self, u: f64, alpha: f64, dt: f64) -> f64 {
        let h = self.h_lambda();
        // two terms of one sign when u and α agree in sign, so nothing cancels
        u * exp(-h * dt) - alpha / h * expm1(-h * dt)
    }

    /// Proper time that passes aboard while thrusting for cosmic time `dt` in the de Sitter limit, starting
    /// from proper velocity `u`: `τ = ∫ dt/γ` along `u(t) = U + (u − U) e^(−H_Λ t)`.
    ///
    /// ```text
    /// τ = 2 / (H_Λ γ_U) · atanh( γ_U s / (γ₁ + γ₀ (1 − s)) ),          s = 1 − e^(−H_Λ dt)
    ///   = 1 / (H_Λ γ_U) · [ H_Λ dt + 2 ln cosh((φ_U + φ₁)/2) − 2 ln cosh((φ_U + φ₀)/2) ]
    /// ```
    ///
    /// where `U = α/H_Λ`, `γ_U = √(1 + U²)`, `φ_U = asinh U`, and the subscripts 0 and 1 mark the start and
    /// the end of the burn (`γ = √(1 + u²)`, `φ = asinh u`). The two forms are equal; the first is used while
    /// its argument is at most 1/2 and the second after that, so neither cancels. For large `dt` this is
    /// `dt/γ_U + φ-terms`: from rest, `τ → dt/γ_U + (1/(H_Λ γ_U)) ln(2 γ_U²/(γ_U + 1))`, which for 1 g is
    /// `dt/1.81e10 + 23.55 yr` (oracle row SR-21).
    ///
    /// Domain: `dt ≥ 0` (NaN otherwise). Relative error below 1e-14, checked against mpmath quadrature.
    pub fn de_sitter_thrust_proper_time(&self, u: f64, alpha: f64, dt: f64) -> f64 {
        if dt.is_nan() || dt < 0.0 {
            return f64::NAN;
        }
        let h = self.h_lambda();
        let big_u = alpha / h;
        let gamma_u = sqrt(1.0 + big_u * big_u);
        let decay = exp(-h * dt);
        let s = -expm1(-h * dt);
        let u_end = u * decay + big_u * s;
        let gamma_start = sqrt(1.0 + u * u);
        let gamma_end = sqrt(1.0 + u_end * u_end);
        let rho = gamma_u * s / (gamma_end + gamma_start * decay);
        if rho <= 0.5 {
            2.0 * atanh(rho) / (h * gamma_u)
        } else {
            let phi_u = asinh(big_u);
            let x_end = 0.5 * (phi_u + asinh(u_end));
            let x_start = 0.5 * (phi_u + asinh(u));
            (h * dt + 2.0 * (ln_cosh(x_end) - ln_cosh(x_start))) / (h * gamma_u)
        }
    }
}
