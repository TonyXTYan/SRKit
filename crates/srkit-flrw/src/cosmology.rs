//! The [`Cosmology`] type, its unit system, and the background functions of cosmic time: scale factor, Hubble
//! rate, age and redshift.

use crate::universal::{theta_of_v, v_of_theta};
use core::f64::consts::LN_2;
use libm::{asinh, cbrt, exp, expm1, log, sinh, sqrt, tanh};
use srkit_core::units::{C, JULIAN_YEAR, PARSEC};

/// One kilometre per second per megaparsec, in inverse seconds: `1e3 / (1e6 pc)` = 3.2408e-20 /s.
pub const KM_S_MPC: f64 = 1.0e3 / (1.0e6 * PARSEC);

/// Scale factor of a source seen today at redshift `z`: `a = 1/(1 + z)`. Dimensionless; needs `z > −1`.
pub fn scale_factor_at_redshift(z: f64) -> f64 {
    1.0 / (1.0 + z)
}

/// Redshift today of light emitted at scale factor `a`: `z = 1/a − 1`. Dimensionless; needs `a > 0`.
pub fn redshift_at_scale_factor(a: f64) -> f64 {
    1.0 / a - 1.0
}

/// A spatially flat universe that holds pressureless matter and a cosmological constant (flat ΛCDM without
/// radiation).
///
/// It is fixed by the Hubble rate today `H0` and the matter density parameter `Ω_m`; `Ω_Λ = 1 − Ω_m`. The scale
/// factor is 1 today and cosmic time `t` counts from the Big Bang of this model.
///
/// # Units
///
/// A cosmology carries its unit of time. Every time is in that unit, every length is the distance light covers
/// in that time (so `c = 1`), every rate is per unit of time, and accelerations are rates of change of rapidity
/// per unit of time. [`Cosmology::per_year`] gives Julian years and light-years, [`Cosmology::per_second`] gives
/// seconds and light-seconds. Functions whose names end in `_si`, temperatures and masses return SI values
/// (kelvin, kilograms, m/s²).
///
/// Comoving distances are quoted as proper distances today (`a = 1`).
///
/// # Radiation is left out
///
/// The closed forms hold for matter and Λ only. Radiation matters before `a ≈ 0.01`: with the measured
/// radiation density the age today is shorter by 0.04%, and the particle horizon today by 1.9%. Everything that
/// looks forward from today (the event horizon, reach, last contact, arrival times, the de Sitter limit) changes
/// by less than one part in 10⁵. Use this type for the future and for epochs after `a ≈ 0.01`.
///
/// ```
/// use srkit_flrw::Cosmology;
/// let u = Cosmology::PLANCK_2018;
/// // 13.8 billion years old, and light sent today reaches 16.7 billion light-years (comoving)
/// assert!((u.age() / 1e9 - 13.80).abs() < 0.01);
/// assert!((u.event_horizon(u.age()) / 1e9 - 16.68).abs() < 0.01);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cosmology {
    /// Hubble rate today, per unit of time.
    h0: f64,
    /// Matter density parameter today.
    omega_m: f64,
    /// Seconds in one unit of time.
    time_unit_s: f64,
}

impl Cosmology {
    /// Planck 2018 base ΛCDM, rounded as in the abstract of Planck Collaboration VI (2020, A&A 641, A6,
    /// arXiv:1807.06209): `H0` = 67.4 km/s/Mpc, `Ω_m` = 0.315. In **Julian years and light-years**:
    /// `H0` = 6.893 079 992 407 984e-11 per year, `1/H0` = 14.507 Gyr, `1/H_Λ` = 17.528 Gyr.
    pub const PLANCK_2018: Self = Self {
        h0: 67.4 * KM_S_MPC * JULIAN_YEAR,
        omega_m: 0.315,
        time_unit_s: JULIAN_YEAR,
    };

    /// The same Planck 2018 cosmology as [`Cosmology::PLANCK_2018`], in **seconds and light-seconds**:
    /// `H0` = 2.184 285 241 085 502e-18 per second.
    pub const PLANCK_2018_SI: Self = Self {
        h0: 67.4 * KM_S_MPC,
        omega_m: 0.315,
        time_unit_s: 1.0,
    };

    /// A cosmology with `h0` per unit of time, where one unit of time is `time_unit_s` seconds.
    ///
    /// Returns `None` unless `h0` and `time_unit_s` are positive and finite and `0 < omega_m < 1`. The
    /// end points are excluded because the closed forms divide by both densities: `omega_m = 1` is the
    /// Einstein-de Sitter universe (no event horizon) and `omega_m = 0` is pure de Sitter space (no Big Bang).
    pub const fn with_time_unit(h0: f64, omega_m: f64, time_unit_s: f64) -> Option<Self> {
        let positive = h0 > 0.0 && h0 < f64::INFINITY;
        let unit = time_unit_s > 0.0 && time_unit_s < f64::INFINITY;
        if positive && unit && omega_m > 0.0 && omega_m < 1.0 {
            Some(Self {
                h0,
                omega_m,
                time_unit_s,
            })
        } else {
            None
        }
    }

    /// A cosmology in Julian years and light-years, from `H0` in inverse Julian years.
    ///
    /// `None` on the same conditions as [`Cosmology::with_time_unit`].
    pub const fn per_year(h0_per_year: f64, omega_m: f64) -> Option<Self> {
        Self::with_time_unit(h0_per_year, omega_m, JULIAN_YEAR)
    }

    /// A cosmology in seconds and light-seconds, from `H0` in inverse seconds.
    ///
    /// `None` on the same conditions as [`Cosmology::with_time_unit`].
    pub const fn per_second(h0_per_second: f64, omega_m: f64) -> Option<Self> {
        Self::with_time_unit(h0_per_second, omega_m, 1.0)
    }

    /// A cosmology in Julian years and light-years, from `H0` in km/s/Mpc as astronomers quote it.
    ///
    /// `None` on the same conditions as [`Cosmology::with_time_unit`].
    pub const fn from_km_s_mpc(h0_km_s_mpc: f64, omega_m: f64) -> Option<Self> {
        Self::per_year(h0_km_s_mpc * KM_S_MPC * JULIAN_YEAR, omega_m)
    }

    /// A cosmology in Julian years and light-years, from the de Sitter rate `H_Λ = H0 √Ω_Λ` in inverse Julian
    /// years. The repository's physics oracle fixes `1/H_Λ` = 17.53 Gyr this way.
    ///
    /// `None` on the same conditions as [`Cosmology::with_time_unit`].
    pub fn from_h_lambda_per_year(h_lambda_per_year: f64, omega_m: f64) -> Option<Self> {
        Self::per_year(h_lambda_per_year / sqrt(1.0 - omega_m), omega_m)
    }

    /// The same universe expressed with a different unit of time, `time_unit_s` seconds long.
    ///
    /// `None` unless `time_unit_s` is positive and finite.
    pub const fn in_time_unit(self, time_unit_s: f64) -> Option<Self> {
        Self::with_time_unit(
            self.h0 / self.time_unit_s * time_unit_s,
            self.omega_m,
            time_unit_s,
        )
    }

    /// Hubble rate today `H0`, per unit of time.
    pub const fn h0(&self) -> f64 {
        self.h0
    }

    /// Matter density parameter today `Ω_m` (dark and baryonic matter together). Dimensionless.
    pub const fn omega_m(&self) -> f64 {
        self.omega_m
    }

    /// Cosmological-constant density parameter today `Ω_Λ = 1 − Ω_m`. Dimensionless.
    pub const fn omega_lambda(&self) -> f64 {
        1.0 - self.omega_m
    }

    /// Seconds in one unit of time.
    pub const fn time_unit_seconds(&self) -> f64 {
        self.time_unit_s
    }

    /// Metres in one unit of length: the distance light covers in one unit of time.
    pub const fn length_unit_metres(&self) -> f64 {
        C * self.time_unit_s
    }

    /// A proper acceleration in m/s² as a rate of change of rapidity per unit of time. One standard gravity
    /// is 1.0323 per Julian year.
    pub const fn acceleration_from_si(&self, metres_per_second2: f64) -> f64 {
        metres_per_second2 / C * self.time_unit_s
    }

    /// The de Sitter rate `H_Λ = H0 √Ω_Λ`, per unit of time: the value the Hubble rate falls to as `t → ∞`.
    /// It is also `c √(Λ/3)`.
    pub fn h_lambda(&self) -> f64 {
        self.h0 * sqrt(self.omega_lambda())
    }

    /// Scale factor at which matter and Λ have equal density, `a_eq = (Ω_m/Ω_Λ)^(1/3)`. Dimensionless; 0.772
    /// for Planck 2018. It is the prefactor of [`Cosmology::scale_factor`].
    pub fn scale_factor_at_equality(&self) -> f64 {
        cbrt(self.omega_m / self.omega_lambda())
    }

    /// `θ = (3/2) H_Λ t`, or NaN outside `t ≥ 0`.
    pub(crate) fn theta(&self, t: f64) -> f64 {
        if t >= 0.0 {
            1.5 * self.h_lambda() * t
        } else {
            f64::NAN
        }
    }

    /// `v = ln(a(t) / a_eq)`, or NaN outside `t ≥ 0`.
    pub(crate) fn v(&self, t: f64) -> f64 {
        v_of_theta(self.theta(t))
    }

    /// Cosmic time at `v = ln(a / a_eq)`.
    pub(crate) fn time_of_v(&self, v: f64) -> f64 {
        theta_of_v(v) / (1.5 * self.h_lambda())
    }

    /// Scale factor at cosmic time `t`. Dimensionless, 1 today.
    ///
    /// ```text
    /// a(t) = (Ω_m/Ω_Λ)^(1/3) sinh^(2/3)((3/2) √Ω_Λ H0 t)
    /// ```
    ///
    /// Domain: `t ≥ 0` (NaN otherwise). The result overflows to infinity near `H_Λ t` = 710 (1.2e13 yr for
    /// Planck 2018); use [`Cosmology::ln_scale_factor`] beyond that. Relative error about `2 (1 + H_Λ t)` units
    /// in the last place: below 1e-15 until the Λ era, 9e-15 at `H_Λ t` = 24. An exponential magnifies the
    /// rounding of its argument, and no binary64 formula can avoid that. Reference: the matter + Λ solution
    /// of the Friedmann equation, e.g. Ryden, *Introduction to Cosmology*, section "Matter + Lambda".
    pub fn scale_factor(&self, t: f64) -> f64 {
        const TWO_THIRDS: f64 = 2.0 / 3.0;
        let theta = self.theta(t);
        if theta < 20.0 {
            let root = cbrt(sinh(theta));
            self.scale_factor_at_equality() * (root * root)
        } else {
            // sinh θ = e^θ / 2 to rounding. One exponential of the whole logarithm, so that the result
            // overflows only when a itself does.
            exp(log(self.omega_m / self.omega_lambda()) / 3.0 + TWO_THIRDS * (theta - LN_2))
        }
    }

    /// Natural logarithm of the scale factor at cosmic time `t`: `ln a = ln a_eq + (2/3) ln sinh θ`.
    ///
    /// Domain: `t ≥ 0` (NaN otherwise; −∞ at `t = 0`). Finite for every finite `t`, which is why stored state
    /// should hold `ln a` and not `a`. Absolute error below `4e-16 (1 + |ln a|)`.
    pub fn ln_scale_factor(&self, t: f64) -> f64 {
        log(self.omega_m / self.omega_lambda()) / 3.0 + self.v(t)
    }

    /// Cosmic time at which the scale factor is `a`: the inverse of [`Cosmology::scale_factor`].
    ///
    /// ```text
    /// t(a) = 2 / (3 H_Λ) · asinh((a / a_eq)^(3/2))
    /// ```
    ///
    /// Domain: `a ≥ 0` (NaN otherwise). Relative error below 4 units in the last place.
    pub fn time_at_scale_factor(&self, a: f64) -> f64 {
        if a.is_nan() || a < 0.0 {
            return f64::NAN;
        }
        let r = a / self.scale_factor_at_equality();
        let theta = if r < 1.0e100 {
            asinh(r * sqrt(r))
        } else {
            1.5 * log(r) + LN_2
        };
        theta / (1.5 * self.h_lambda())
    }

    /// Cosmic time at which the scale factor is `e^ln_a`: the inverse of [`Cosmology::ln_scale_factor`], for
    /// epochs where `a` itself would overflow.
    pub fn time_at_ln_scale_factor(&self, ln_a: f64) -> f64 {
        self.time_of_v(ln_a - log(self.omega_m / self.omega_lambda()) / 3.0)
    }

    /// Hubble rate at cosmic time `t`, per unit of time.
    ///
    /// ```text
    /// H(t) = H_Λ coth((3/2) H_Λ t)
    /// ```
    ///
    /// Domain: `t ≥ 0` (NaN otherwise; infinite at `t = 0`). Tends to `2/(3t)` early and to `H_Λ` late.
    /// Relative error below 3 units in the last place.
    pub fn hubble(&self, t: f64) -> f64 {
        self.h_lambda() / tanh(self.theta(t))
    }

    /// Hubble rate at scale factor `a`, per unit of time.
    ///
    /// ```text
    /// H(a) = H0 √(Ω_m a⁻³ + Ω_Λ)
    /// ```
    ///
    /// Domain: `a > 0`. Relative error below 3 units in the last place.
    pub fn hubble_at_scale_factor(&self, a: f64) -> f64 {
        self.h0 * sqrt(self.omega_m / (a * a * a) + self.omega_lambda())
    }

    /// Age of the universe today (`a = 1`), in units of time.
    ///
    /// ```text
    /// t0 = 2 / (3 H_Λ) · asinh(√(Ω_Λ/Ω_m))
    /// ```
    ///
    /// 13.80 Gyr for Planck 2018 (13.79 Gyr once radiation is included, which this model leaves out).
    pub fn age(&self) -> f64 {
        asinh(sqrt(self.omega_lambda() / self.omega_m)) / (1.5 * self.h_lambda())
    }

    /// Cosmic time at which light seen today with redshift `z` was emitted: `t(a = 1/(1 + z))`.
    ///
    /// Domain: `z > −1`; negative `z` gives a time in the future.
    pub fn time_at_redshift(&self, z: f64) -> f64 {
        self.time_at_scale_factor(scale_factor_at_redshift(z))
    }

    /// Lookback time to redshift `z`: `t0 − t(z)`, in units of time.
    pub fn lookback_time(&self, z: f64) -> f64 {
        self.age() - self.time_at_redshift(z)
    }

    /// Redshift `z = a(t_received)/a(t_emitted) − 1` of light exchanged between two comoving observers.
    ///
    /// Domain: both times `≥ 0`. Negative (a blueshift) if `t_received < t_emitted`. The ratio is formed from
    /// logarithms, so it stays finite after `a` overflows; it overflows itself when the scale factor grows by
    /// more than e⁷⁰⁹ between the two times. The playback rate of a comoving source is `1/(1 + z)`.
    pub fn redshift_between(&self, t_emitted: f64, t_received: f64) -> f64 {
        expm1(self.v(t_received) - self.v(t_emitted))
    }
}
