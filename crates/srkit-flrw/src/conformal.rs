//! Conformal time and everything read off the conformal chart: comoving distances, the particle and event
//! horizons, reach, last contact and light cones.
//!
//! Flat FLRW space is conformally flat, `ds² = a²(η)(−dη² + dχ²)`, so light moves on straight 45° lines in
//! conformal time `η = ∫ dt/a` and comoving distance `χ`. The integral converges as `t → ∞`: conformal time has
//! a finite limit `η_∞`, and that is the whole reason an event horizon exists (Rindler 1956, MNRAS 116, 662;
//! Davis & Lineweaver 2004, arXiv:astro-ph/0310808).

use crate::cosmology::{scale_factor_at_redshift, Cosmology};
use crate::universal::{
    eta_hat, eta_hat_between, proper_rem_hat, rem_hat, v_of_eta_hat, v_of_rem_hat, ETA_HAT_TOTAL,
};
use libm::log;
use srkit_core::units::{C, G_NEWTON};

/// Where and when a light signal reaches a comoving observer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Reception {
    /// Conformal time of reception, `η_emitted + χ`, in units of time.
    pub conformal_time: f64,
    /// Cosmic time of reception, in units of time.
    pub time: f64,
}

impl Cosmology {
    /// `1/(a_eq H_Λ) = 1/(H0 Ω_m^(1/3) Ω_Λ^(1/6))`: the unit in which conformal time is a universal function.
    fn conformal_unit(&self) -> f64 {
        1.0 / (self.scale_factor_at_equality() * self.h_lambda())
    }

    /// Conformal time at cosmic time `t`: `η(t) = ∫₀ᵗ dt′/a(t′)`, in units of time (equally, the comoving
    /// distance light has covered since the Big Bang, in units of length).
    ///
    /// ```text
    /// η = 2/(a H) · 2F1(1, 1/2; 7/6; Ω_Λ(a))                Ω_Λ(a) = x/(1 + x),  x = (Ω_Λ/Ω_m) a³
    ///   = 2 √a / (H0 √Ω_m) · 2F1(1/6, 1/2; 7/6; −x)
    /// ```
    ///
    /// The first form is the second after a Pfaff transformation and is summed directly for `x ≤ 1`; for
    /// `x > 1` the result is `η_∞` minus [`Cosmology::conformal_time_remaining`]. Domain: `t ≥ 0` (NaN
    /// otherwise). Relative error below 1e-15 (5 units in the last place), measured against mpmath
    /// quadrature. Late in the Λ era `η` is within rounding of `η_∞`, so ask for the remaining conformal time
    /// there. Reference: Baes, Camps & Van De Putte 2017, MNRAS 468, 927.
    pub fn conformal_time(&self, t: f64) -> f64 {
        self.conformal_unit() * eta_hat(self.v(t))
    }

    /// Conformal time at scale factor `a`: [`Cosmology::conformal_time`] as a function of `a ≥ 0`.
    pub fn conformal_time_at_scale_factor(&self, a: f64) -> f64 {
        self.conformal_unit() * eta_hat(self.v_at_scale_factor(a))
    }

    /// The limit of conformal time as `t → ∞`, in units of time.
    ///
    /// ```text
    /// η_∞ = Γ(1/6) Γ(1/3) / (3 √π) / (H0 Ω_m^(1/3) Ω_Λ^(1/6)) = 2.804 364 210 650 908 5 / (a_eq H_Λ)
    /// ```
    ///
    /// As a length it is the comoving radius of everything that can ever be seen from here: 63.7 billion
    /// light-years for Planck 2018.
    pub fn conformal_time_limit(&self) -> f64 {
        self.conformal_unit() * ETA_HAT_TOTAL
    }

    /// Conformal time left at cosmic time `t`: `η_∞ − η(t) = ∫ₜ^∞ dt′/a(t′)`, in units of time.
    ///
    /// ```text
    /// η_∞ − η = 1/(a H) · 2F1(1, 1/2; 4/3; Ω_m(a))            Ω_m(a) = 1/(1 + x),  x = (Ω_Λ/Ω_m) a³
    ///         = 1/(a H_Λ) · 2F1(1/3, 1/2; 4/3; −1/x)
    /// ```
    ///
    /// summed directly for `x ≥ 1`, and `η_∞ − η(t)` for `x < 1`. It tends to `1/(a H_Λ)`, which is
    /// `e^(−H_Λ t)` up to a constant, so it stays accurate when `η` itself no longer resolves the difference.
    /// This is the quantity to store. Domain: `t ≥ 0` (NaN otherwise). It underflows to zero at `H_Λ t` =
    /// 745. Relative error about `(1 + H_Λ t)` units in the last place: below 1e-15 until the Λ era and
    /// 1.5e-13 at `H_Λ t` = 600, because a result that falls as `e^(−H_Λ t)` magnifies the rounding of `t`.
    pub fn conformal_time_remaining(&self, t: f64) -> f64 {
        self.conformal_unit() * rem_hat(self.v(t))
    }

    /// Cosmic time at conformal time `eta`: the inverse of [`Cosmology::conformal_time`].
    ///
    /// `None` unless `0 ≤ eta < η_∞`. Solved by Newton's method on `ln η` against `ln a`, which is concave, from
    /// a start on the safe side of the root: at most 40 steps, in practice 3 to 7, with no dependence on
    /// anything but the arguments. Relative error about `(1 + η/(η_∞ − η))` units in the last place: when
    /// `eta` is within rounding of `η_∞` the answer is not determined by it, and
    /// [`Cosmology::time_at_conformal_time_remaining`] is the function to use.
    pub fn time_at_conformal_time(&self, eta: f64) -> Option<f64> {
        let e = eta / self.conformal_unit();
        if e == 0.0 {
            Some(0.0)
        } else if e > 0.0 && e < ETA_HAT_TOTAL {
            Some(self.time_of_v(v_of_eta_hat(e)))
        } else {
            None
        }
    }

    /// Cosmic time at which the conformal time left is `remaining`: the inverse of
    /// [`Cosmology::conformal_time_remaining`].
    ///
    /// `None` unless `0 < remaining ≤ η_∞`. Same method and bounds as [`Cosmology::time_at_conformal_time`].
    /// Relative error a few units in the last place from matter-Λ equality on; before it, `remaining` sits
    /// close to `η_∞` and the error grows as `3 η_∞/η` units in the last place (1e-13 at 1,400 years).
    pub fn time_at_conformal_time_remaining(&self, remaining: f64) -> Option<f64> {
        let unit = self.conformal_unit();
        if remaining > 0.0 && remaining <= unit * ETA_HAT_TOTAL {
            let r = remaining / unit;
            if r < ETA_HAT_TOTAL {
                Some(self.time_of_v(v_of_rem_hat(r)))
            } else {
                Some(0.0)
            }
        } else {
            None
        }
    }

    /// Comoving distance light covers between cosmic times `t_from` and `t_to`: `η(t_to) − η(t_from)`, in units
    /// of length. Negative if `t_to < t_from`.
    ///
    /// The difference is taken inside one closed form when both times lie on the same side of matter-Λ
    /// equality. It still cancels when the two times are close: the relative error is about
    /// `1e-16 · η / |Δη|`.
    pub fn comoving_distance(&self, t_from: f64, t_to: f64) -> f64 {
        self.conformal_unit() * eta_hat_between(self.v(t_from), self.v(t_to))
    }

    /// Comoving distance today to a source seen at redshift `z`, in units of length (the line-of-sight comoving
    /// distance of Hogg 1999, arXiv:astro-ph/9905116). Domain: `z > −1`.
    pub fn comoving_distance_at_redshift(&self, z: f64) -> f64 {
        let v_then = self.v_at_scale_factor(scale_factor_at_redshift(z));
        let v_now = self.v_at_scale_factor(1.0);
        self.conformal_unit() * eta_hat_between(v_then, v_now)
    }

    /// Comoving particle horizon at cosmic time `t`: the comoving distance to the farthest source whose light
    /// has had time to arrive. Equal to [`Cosmology::conformal_time`], in units of length.
    pub fn particle_horizon(&self, t: f64) -> f64 {
        self.conformal_time(t)
    }

    /// Proper distance to the particle horizon at cosmic time `t`: `a(t) η(t)`, in units of length.
    pub fn particle_horizon_proper(&self, t: f64) -> f64 {
        self.scale_factor(t) * self.conformal_time(t)
    }

    /// Comoving cosmic event horizon at cosmic time `t`: light emitted farther away than this at `t` never
    /// arrives here, and light sent from here at `t` never gets farther than this. Equal to
    /// [`Cosmology::conformal_time_remaining`], in units of length. 16.7 billion light-years today for
    /// Planck 2018.
    pub fn event_horizon(&self, t: f64) -> f64 {
        self.conformal_time_remaining(t)
    }

    /// Proper distance to the event horizon at cosmic time `t`: `a(t) (η_∞ − η(t))`, in units of length.
    ///
    /// ```text
    /// a (η_∞ − η) = 1/H · 2F1(1, 1/2; 4/3; Ω_m(a))
    /// ```
    ///
    /// It grows from 0 towards the de Sitter horizon `1/H_Λ` and stays finite after `a` overflows.
    pub fn event_horizon_proper(&self, t: f64) -> f64 {
        proper_rem_hat(self.v(t)) / self.h_lambda()
    }

    /// Comoving distance that light (or anything slower) leaving here at cosmic time `t` can ever cover:
    /// the reach. It is the same number as [`Cosmology::event_horizon`], read as a statement about the future
    /// light cone.
    pub fn reach(&self, t: f64) -> f64 {
        self.conformal_time_remaining(t)
    }

    /// Share of the observable universe at cosmic time `t` that can still be reached from here at `t`:
    /// `((η_∞ − η)/η)³`, the ratio of the comoving volumes inside the event horizon and the particle horizon.
    /// Matter is uniform in comoving volume, so it is also the share by mass. Dimensionless.
    ///
    /// It is exactly 1 at `a = a_eq/2` (0.386 for Planck 2018), where the two horizons cross because
    /// `η(a_eq/2) = η_∞/2` in every flat matter + Λ universe. It is 4.5% today for Planck 2018. Radiation,
    /// which this model leaves out, shrinks the particle horizon by 1.9% and so raises today's share by 6%
    /// of itself: `tools/physics-checks/lcdm_reach.py` gets 4.6% with radiation where this function gives
    /// 4.3% for the same `H0` and `Ω_m`.
    pub fn reachable_share_of_observable(&self, t: f64) -> f64 {
        let v = self.v(t);
        let ratio = rem_hat(v) / eta_hat(v);
        ratio * ratio * ratio
    }

    /// Comoving volume still reachable at cosmic time `t` as a fraction of that reachable at `t_reference`:
    /// `(reach(t)/reach(t_reference))³`. Dimensionless.
    pub fn reach_volume_ratio(&self, t: f64, t_reference: f64) -> f64 {
        let ratio = rem_hat(self.v(t)) / rem_hat(self.v(t_reference));
        ratio * ratio * ratio
    }

    /// Mass of matter inside the event horizon at cosmic time `t`, in kilograms: all the matter a traveller
    /// leaving at `t` could ever reach.
    ///
    /// ```text
    /// M = ρ_m a³ · (4π/3) χ³ = c³ F³ / (2 G H_Λ),      F = a_eq H_Λ (η_∞ − η(t))
    /// ```
    ///
    /// `F` depends on `H_Λ t` alone, so the reach mass does not depend on `Ω_m` at a given `H_Λ` and `H_Λ t`.
    /// `F(0) = Γ(1/6) Γ(1/3)/(3 √π)` = 2.8044. The mass counts dark and baryonic matter (`Ω_m`). Reference:
    /// `tools/physics-checks/reach_mass_closed_form.py` and `docs/research/constants-gameplay.md` section 6.
    /// The value carries the 2.2e-5 relative uncertainty of `G`.
    pub fn reach_mass(&self, t: f64) -> f64 {
        let f = rem_hat(self.v(t));
        let h_lambda_si = self.h_lambda() / self.time_unit_seconds();
        C * C * C * f * f * f / (2.0 * G_NEWTON * h_lambda_si)
    }

    /// Latest cosmic time at which a signal sent from here still reaches a comoving target at comoving distance
    /// `chi`: the root of `η_∞ − η(t) = chi`. By symmetry it is also the latest time the target can send a
    /// signal that reaches here.
    ///
    /// `None` unless `0 ≤ chi ≤ η_∞`; a target farther than `η_∞` was never in causal contact. `chi = 0` gives
    /// infinity. A signal sent exactly at this time arrives only as `t → ∞`. The root is found as in
    /// [`Cosmology::time_at_conformal_time`]: `ln(η_∞ − η)` is monotone and concave in `ln a`, Newton's method
    /// starts from the de Sitter estimate on the safe side, and takes at most 40 steps (3 to 7 in practice).
    /// References: Loeb 2002, Phys. Rev. D 65, 047301; Busha et al. 2003, ApJ 596, 713.
    pub fn last_contact_time(&self, chi: f64) -> Option<f64> {
        if chi == 0.0 {
            Some(f64::INFINITY)
        } else {
            self.time_at_conformal_time_remaining(chi)
        }
    }

    /// Emission time of the last light that will ever arrive from a comoving source at comoving distance
    /// `chi`: the source's clock reading at which its image freezes. Equal to
    /// [`Cosmology::last_contact_time`]; the source is seen to approach this moment ever more slowly and ever
    /// more redshifted, and is never seen to pass it.
    ///
    /// This is the horizon definition. The photon-count definition of Last Light in
    /// `docs/engine/deep-time-engine.md` section 4 depends on the telescope; see
    /// [`Cosmology::last_light_gain`] for its aperture scaling.
    pub fn last_light_time(&self, chi: f64) -> Option<f64> {
        self.last_contact_time(chi)
    }

    /// Latest cosmic time at which a signal sent from here to a comoving target at comoving distance `chi`
    /// can still be answered: the root of `η_∞ − η(t) = 2 chi`, for a reply sent the moment the signal
    /// arrives.
    ///
    /// `None` unless `0 ≤ 2 chi ≤ η_∞`.
    pub fn last_reply_time(&self, chi: f64) -> Option<f64> {
        self.last_contact_time(2.0 * chi)
    }

    /// Reception of a light signal emitted at conformal time `eta_emitted` towards a comoving observer at
    /// comoving distance `chi`.
    ///
    /// In the conformal chart the light cone is `η_received = η_emitted + χ`, exactly as in Minkowski space.
    /// Returns the reception conformal time and the cosmic time it corresponds to, or `None` if the signal
    /// never arrives (`η_emitted + χ ≥ η_∞`) or the arguments are outside `eta_emitted ≥ 0`, `chi ≥ 0`.
    ///
    /// Late in the Λ era conformal times crowd against `η_∞` and this form loses precision;
    /// [`Cosmology::light_arrival_time`] works from the conformal time left and does not.
    pub fn light_reception(&self, eta_emitted: f64, chi: f64) -> Option<Reception> {
        if eta_emitted >= 0.0 && chi >= 0.0 {
            let conformal_time = eta_emitted + chi;
            let time = self.time_at_conformal_time(conformal_time)?;
            Some(Reception {
                conformal_time,
                time,
            })
        } else {
            None
        }
    }

    /// Cosmic time at which light emitted at cosmic time `t_emitted` reaches a comoving observer at comoving
    /// distance `chi`.
    ///
    /// `None` if the light never arrives (`chi ≥ η_∞ − η(t_emitted)`, the target is beyond the event horizon)
    /// or the arguments are outside `t_emitted ≥ 0`, `chi ≥ 0`. The subtraction is done on the conformal time
    /// left, so the answer keeps its precision at any epoch; the delay `t_received − t_emitted` has relative
    /// error about `1e-16 (η_∞ − η)/χ`.
    pub fn light_arrival_time(&self, t_emitted: f64, chi: f64) -> Option<f64> {
        if chi == 0.0 && t_emitted >= 0.0 {
            return Some(t_emitted);
        }
        let left = self.conformal_time_remaining(t_emitted) - chi;
        if chi > 0.0 && left > 0.0 {
            self.time_at_conformal_time_remaining(left)
        } else {
            None
        }
    }

    /// Cosmic time at which light received at cosmic time `t_received` left a comoving source at comoving
    /// distance `chi`: the past light cone.
    ///
    /// `None` if the source is beyond the particle horizon (`chi > η(t_received)`) or the arguments are outside
    /// `t_received ≥ 0`, `chi ≥ 0`.
    pub fn light_emission_time(&self, t_received: f64, chi: f64) -> Option<f64> {
        if chi == 0.0 && t_received >= 0.0 {
            return Some(t_received);
        }
        let v = self.v(t_received);
        let unit = self.conformal_unit();
        if chi > 0.0 && v <= 0.0 {
            self.time_at_conformal_time(unit * eta_hat(v) - chi)
        } else if chi > 0.0 && v > 0.0 {
            self.time_at_conformal_time_remaining(unit * rem_hat(v) + chi)
        } else {
            None
        }
    }

    /// `v = ln(a / a_eq)` for a scale factor, or NaN outside `a ≥ 0`.
    fn v_at_scale_factor(&self, a: f64) -> f64 {
        if a >= 0.0 {
            log(a / self.scale_factor_at_equality())
        } else {
            f64::NAN
        }
    }
}
