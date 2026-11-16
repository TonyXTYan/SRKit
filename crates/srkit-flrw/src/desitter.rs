//! The de Sitter limit: what the universe tends to once Λ dominates, and the temperatures that set the far
//! future's heat sink.

use crate::cosmology::Cosmology;
use core::f64::consts::PI;
use libm::{exp, log};
use srkit_core::units::{BOLTZMANN, HBAR};

/// Temperature of the cosmic microwave background today, in kelvin: 2.7255 K (Fixsen 2009, ApJ 707, 916, who
/// gives 2.72548 ± 0.00057 K).
pub const T_CMB0: f64 = 2.7255;

/// Temperature of the cosmic microwave background at scale factor `a`, in kelvin: `T = t_cmb0 / a`, where
/// `t_cmb0` is its temperature today in kelvin (see [`T_CMB0`]). Domain: `a > 0`.
///
/// A blackbody stays a blackbody under expansion and its temperature redshifts like its photons.
pub fn cmb_temperature_at_scale_factor(a: f64, t_cmb0: f64) -> f64 {
    t_cmb0 / a
}

impl Cosmology {
    /// The de Sitter rate `H_Λ` in inverse seconds. 1.808e-18 /s for Planck 2018.
    pub fn h_lambda_si(&self) -> f64 {
        self.h_lambda() / self.time_unit_seconds()
    }

    /// The de Sitter horizon `c/H_Λ`, in units of length: the proper distance the event horizon tends to.
    /// 17.53 billion light-years for Planck 2018.
    pub fn de_sitter_horizon(&self) -> f64 {
        1.0 / self.h_lambda()
    }

    /// The Gibbons-Hawking temperature of the de Sitter horizon, in kelvin.
    ///
    /// ```text
    /// T_dS = ħ H_Λ / (2π k_B)
    /// ```
    ///
    /// 2.20e-30 K (oracle row DT-15). It is the coldest heat sink the far future offers. Reference: Gibbons &
    /// Hawking 1977, Phys. Rev. D 15, 2738.
    pub fn de_sitter_temperature(&self) -> f64 {
        HBAR * self.h_lambda_si() / (2.0 * PI * BOLTZMANN)
    }

    /// Temperature of the cosmic microwave background at cosmic time `t`, in kelvin: `T = t_cmb0 / a(t)`, where
    /// `t_cmb0` is its temperature today in kelvin (see [`T_CMB0`]).
    ///
    /// Domain: `t ≥ 0`. Formed from `ln a`, so it falls smoothly to zero and does not fail when `a` overflows.
    /// The radiation is a tracer here; it does not act on the expansion.
    pub fn cmb_temperature(&self, t: f64, t_cmb0: f64) -> f64 {
        t_cmb0 * exp(-self.ln_scale_factor(t))
    }

    /// Scale factor at which the microwave background has cooled to the de Sitter temperature:
    /// `a = t_cmb0 / T_dS`. Dimensionless; 1.24e30, or 69.3 e-folds from today (oracle rows DT-16 and DT-17).
    pub fn cmb_crossover_scale_factor(&self, t_cmb0: f64) -> f64 {
        t_cmb0 / self.de_sitter_temperature()
    }

    /// Cosmic time at which the microwave background has cooled to the de Sitter temperature, in units of
    /// time: 1.23e12 yr (oracle row DT-18). After it, waiting no longer makes the sky colder.
    pub fn cmb_crossover_time(&self, t_cmb0: f64) -> f64 {
        self.time_at_scale_factor(self.cmb_crossover_scale_factor(t_cmb0))
    }

    /// Temperature of the coldest available heat sink at cosmic time `t`, in kelvin:
    /// `max(T_CMB(t), T_dS)` (`docs/engine/deep-time-engine.md` section 5).
    pub fn sink_temperature(&self, t: f64, t_cmb0: f64) -> f64 {
        let cmb = self.cmb_temperature(t, t_cmb0);
        if cmb.is_nan() {
            cmb
        } else {
            cmb.max(self.de_sitter_temperature())
        }
    }

    /// Outward tidal acceleration that Λ produces at distance `r` from the centre of a bound system:
    /// `H_Λ² r`, as a rate of change of rapidity per unit of time (`r` in units of length). It equals
    /// `Λ c² r / 3`.
    pub fn lambda_tidal_acceleration(&self, r: f64) -> f64 {
        let h = self.h_lambda();
        h * h * r
    }

    /// Outward tidal acceleration that Λ produces at distance `r_metres`, in m/s²: `H_Λ² r`. At 50,000
    /// light-years it is 1.55e-15 m/s² (oracle row DT-25), about 10⁵ times less than the galaxy's own gravity
    /// there, which is why bound systems do not expand.
    pub fn lambda_tidal_acceleration_si(&self, r_metres: f64) -> f64 {
        let h = self.h_lambda_si();
        h * h * r_metres
    }

    /// Recession speed of a comoving point at proper distance `d` at cosmic time `t`, as a fraction of `c`:
    /// the Hubble law `H(t) d`. It exceeds 1 beyond the Hubble radius; that is not a signal speed.
    pub fn recession_speed(&self, t: f64, d: f64) -> f64 {
        self.hubble(t) * d
    }

    /// Recession speed across proper distance `d` in the de Sitter limit, as a fraction of `c`: `H_Λ d`.
    /// Across 100,000 light-years it is 5.7e-6 (oracle row DT-26), the size of the error made by treating a
    /// galaxy as a patch of flat spacetime.
    pub fn de_sitter_recession_speed(&self, d: f64) -> f64 {
        self.h_lambda() * d
    }

    /// Extra time, in units of time, for which a receding source keeps delivering photons when the collecting
    /// area grows by `area_factor`: `ln(area_factor) / (3 H_Λ)`.
    ///
    /// Once Λ dominates, the photon rate from a comoving source into a fixed aperture falls as `e^(−3 H_Λ t)`
    /// (one factor from time dilation, two from the aperture's shrinking share of the emitted sphere), so the
    /// moment the expected number of photons still to come drops below one moves later by this much. Ten times
    /// the area buys 13.5 Gyr (oracle row DT-27). This is the photon-count definition of Last Light in
    /// `docs/engine/deep-time-engine.md` section 4; the background-limited variant is not settled (R37.5) and
    /// is not implemented.
    pub fn last_light_gain(&self, area_factor: f64) -> f64 {
        log(area_factor) / (3.0 * self.h_lambda())
    }

    /// Cosmic time at which the scale factor reaches the largest finite binary64 value, in units of time:
    /// `H_Λ t` = 710.5, `t` = 1.2e13 yr (oracle rows DT-28 and DT-29). Beyond it
    /// [`Cosmology::scale_factor`] returns infinity and [`Cosmology::ln_scale_factor`] must be used.
    pub fn scale_factor_overflow_time(&self) -> f64 {
        self.time_at_ln_scale_factor(log(f64::MAX))
    }
}
