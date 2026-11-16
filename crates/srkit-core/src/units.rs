//! Exact SI constants and the conversion between SI and the crate's natural units.
//!
//! Constants carry their full official precision (CODATA 2018, the 2019 SI redefinition, IAU 2012 and 2015
//! resolutions). Round them only when you show them.

/// Speed of light in vacuum, m/s. Exact by the definition of the metre.
pub const C: f64 = 299_792_458.0;
/// Julian year, s: exactly 365.25 days of 86,400 s.
pub const JULIAN_YEAR: f64 = 31_557_600.0;
/// Day, s.
pub const DAY: f64 = 86_400.0;
/// Light-year, m: `C * JULIAN_YEAR` (IAU), exactly 9,460,730,472,580,800 m.
pub const LIGHT_YEAR: f64 = 9_460_730_472_580_800.0;
/// Astronomical unit, m. Exact (IAU 2012 Resolution B2).
pub const AU: f64 = 149_597_870_700.0;
/// Parsec, m: `648000 / pi` astronomical units (IAU 2015 Resolution B2).
pub const PARSEC: f64 = 3.085_677_581_491_367e16;
/// Standard gravity, m/s². Exact (a conventional value, not a local field).
pub const G0: f64 = 9.806_65;
/// Planck constant, J s. Exact (2019 SI).
pub const PLANCK: f64 = 6.626_070_15e-34;
/// Reduced Planck constant, J s.
pub const HBAR: f64 = 1.054_571_817_646_156_5e-34;
/// Boltzmann constant, J/K. Exact (2019 SI).
pub const BOLTZMANN: f64 = 1.380_649e-23;
/// Newtonian constant of gravitation, m³ kg⁻¹ s⁻² (CODATA 2018; relative uncertainty 2.2e-5).
pub const G_NEWTON: f64 = 6.674_30e-11;
/// Stefan-Boltzmann constant, W m⁻² K⁻⁴ (exact in the 2019 SI; rounded to binary64).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419_184_429_5e-8;
/// Nominal solar luminosity, W (IAU 2015 Resolution B3).
pub const L_SUN: f64 = 3.828e26;
/// Nominal solar radius, m (IAU 2015 Resolution B3).
pub const R_SUN: f64 = 6.957e8;
/// Nominal solar mass parameter `GM`, m³/s² (IAU 2015 Resolution B3).
pub const GM_SUN: f64 = 1.327_124_4e20;
/// Solar mass, kg: `GM_SUN / G_NEWTON`.
pub const M_SUN: f64 = GM_SUN / G_NEWTON;

/// Standard gravity in natural units of years: `G0 * JULIAN_YEAR / C` = 1.0323 per year (ly/yr²).
pub const G0_PER_YEAR: f64 = G0 * JULIAN_YEAR / C;
/// Astronomical unit in light-years.
pub const AU_LY: f64 = AU / LIGHT_YEAR;
/// Parsec in light-years.
pub const PARSEC_LY: f64 = PARSEC / LIGHT_YEAR;

/// A natural unit system: one unit of time, and the distance light covers in it as the unit of length.
///
/// ```
/// use srkit_core::units::Units;
/// let yr = Units::YEARS;
/// // 1 g as a rate of change of rapidity per year
/// assert!((yr.accel_from_si(9.80665) - 1.0323).abs() < 1e-4);
/// // 4.2465 light-years in metres
/// assert!((yr.length_to_si(4.2465) / 4.0175e16 - 1.0).abs() < 1e-4);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Units {
    /// Seconds in one unit of time.
    pub time_s: f64,
}

impl Units {
    /// Seconds and light-seconds.
    pub const SECONDS: Self = Self { time_s: 1.0 };
    /// Julian years and light-years.
    pub const YEARS: Self = Self {
        time_s: JULIAN_YEAR,
    };
    /// Metres of length and metres of light-travel time (geometrised).
    pub const METRES: Self = Self { time_s: 1.0 / C };

    /// A system whose unit of time is `time_s` seconds.
    pub const fn new(time_s: f64) -> Self {
        Self { time_s }
    }
    /// Metres in one unit of length.
    pub fn length_m(self) -> f64 {
        C * self.time_s
    }
    /// Seconds to units of time.
    pub fn time_from_si(self, seconds: f64) -> f64 {
        seconds / self.time_s
    }
    /// Units of time to seconds.
    pub fn time_to_si(self, time: f64) -> f64 {
        time * self.time_s
    }
    /// Metres to units of length.
    pub fn length_from_si(self, metres: f64) -> f64 {
        metres / C / self.time_s
    }
    /// Units of length to metres.
    pub fn length_to_si(self, length: f64) -> f64 {
        length * self.time_s * C
    }
    /// m/s to a fraction of `c`.
    pub fn speed_from_si(self, metres_per_second: f64) -> f64 {
        metres_per_second / C
    }
    /// A fraction of `c` to m/s.
    pub fn speed_to_si(self, beta: f64) -> f64 {
        beta * C
    }
    /// m/s² to a rate of change of rapidity per unit of time.
    pub fn accel_from_si(self, metres_per_second2: f64) -> f64 {
        metres_per_second2 / C * self.time_s
    }
    /// A rate of change of rapidity per unit of time to m/s².
    pub fn accel_to_si(self, alpha: f64) -> f64 {
        alpha / self.time_s * C
    }
}
