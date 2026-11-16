//! `srkit.units`: the exact SI constants and the conversion between SI and natural units.

use pyo3::prelude::*;
use srkit_core::units as core;

use crate::convert::{fl, positive};

/// A natural unit system: one unit of time, and the distance light covers in it as the unit of length.
///
/// `Units(time_s)` makes a system whose unit of time is `time_s` seconds. Three are ready made:
/// `Units.SECONDS` (seconds and light-seconds), `Units.YEARS` (Julian years and light-years) and
/// `Units.METRES` (metres, and metres of light-travel time).
///
/// In natural units speeds are fractions of c, rapidities have no unit and a proper acceleration is a rate of
/// change of rapidity per unit of time.
#[pyclass(frozen, eq, module = "srkit.units", name = "Units")]
#[derive(Clone, Copy, PartialEq)]
pub struct Units {
    inner: core::Units,
}

#[pymethods]
impl Units {
    #[new]
    fn new(time_s: f64) -> PyResult<Self> {
        Ok(Self {
            inner: core::Units::new(positive("time_s", time_s)?),
        })
    }

    /// Seconds and light-seconds.
    #[classattr]
    #[pyo3(name = "SECONDS")]
    fn seconds() -> Self {
        Self {
            inner: core::Units::SECONDS,
        }
    }
    /// Julian years and light-years.
    #[classattr]
    #[pyo3(name = "YEARS")]
    fn years() -> Self {
        Self {
            inner: core::Units::YEARS,
        }
    }
    /// Metres of length and metres of light-travel time (geometrised).
    #[classattr]
    #[pyo3(name = "METRES")]
    fn metres() -> Self {
        Self {
            inner: core::Units::METRES,
        }
    }

    /// Seconds in one unit of time.
    #[getter]
    fn time_s(&self) -> f64 {
        self.inner.time_s
    }
    /// Metres in one unit of length.
    #[getter]
    fn length_m(&self) -> f64 {
        self.inner.length_m()
    }
    /// Seconds to units of time.
    fn time_from_si(&self, seconds: f64) -> f64 {
        self.inner.time_from_si(seconds)
    }
    /// Units of time to seconds.
    fn time_to_si(&self, time: f64) -> f64 {
        self.inner.time_to_si(time)
    }
    /// Metres to units of length.
    fn length_from_si(&self, metres: f64) -> f64 {
        self.inner.length_from_si(metres)
    }
    /// Units of length to metres.
    fn length_to_si(&self, length: f64) -> f64 {
        self.inner.length_to_si(length)
    }
    /// Metres per second to a fraction of c.
    fn speed_from_si(&self, metres_per_second: f64) -> f64 {
        self.inner.speed_from_si(metres_per_second)
    }
    /// A fraction of c to metres per second.
    fn speed_to_si(&self, beta: f64) -> f64 {
        self.inner.speed_to_si(beta)
    }
    /// Metres per second squared to a rate of change of rapidity per unit of time.
    fn accel_from_si(&self, metres_per_second2: f64) -> f64 {
        self.inner.accel_from_si(metres_per_second2)
    }
    /// A rate of change of rapidity per unit of time to metres per second squared.
    fn accel_to_si(&self, alpha: f64) -> f64 {
        self.inner.accel_to_si(alpha)
    }

    fn __repr__(&self) -> String {
        let t = self.inner.time_s;
        if self.inner == core::Units::SECONDS {
            "Units.SECONDS".to_owned()
        } else if self.inner == core::Units::YEARS {
            "Units.YEARS".to_owned()
        } else if self.inner == core::Units::METRES {
            "Units.METRES".to_owned()
        } else {
            format!("Units(time_s={})", fl(t))
        }
    }
}

/// The constants, in the order the module lists them: name, value, one line of documentation.
pub(crate) const CONSTANTS: &[(&str, f64, &str)] = &[
    ("C", core::C, "Speed of light in vacuum, m/s. Exact."),
    (
        "JULIAN_YEAR",
        core::JULIAN_YEAR,
        "Julian year, s: exactly 365.25 days of 86,400 s.",
    ),
    ("DAY", core::DAY, "Day, s."),
    (
        "LIGHT_YEAR",
        core::LIGHT_YEAR,
        "Light-year, m: C * JULIAN_YEAR, exactly 9,460,730,472,580,800 m.",
    ),
    ("AU", core::AU, "Astronomical unit, m. Exact (IAU 2012)."),
    (
        "PARSEC",
        core::PARSEC,
        "Parsec, m: 648000 / pi astronomical units (IAU 2015).",
    ),
    ("G0", core::G0, "Standard gravity, m/s^2. Exact."),
    ("PLANCK", core::PLANCK, "Planck constant, J s. Exact."),
    ("HBAR", core::HBAR, "Reduced Planck constant, J s."),
    (
        "BOLTZMANN",
        core::BOLTZMANN,
        "Boltzmann constant, J/K. Exact.",
    ),
    (
        "G_NEWTON",
        core::G_NEWTON,
        "Newtonian constant of gravitation, m^3 kg^-1 s^-2 (CODATA 2018).",
    ),
    (
        "STEFAN_BOLTZMANN",
        core::STEFAN_BOLTZMANN,
        "Stefan-Boltzmann constant, W m^-2 K^-4.",
    ),
    (
        "L_SUN",
        core::L_SUN,
        "Nominal solar luminosity, W (IAU 2015).",
    ),
    ("R_SUN", core::R_SUN, "Nominal solar radius, m (IAU 2015)."),
    (
        "GM_SUN",
        core::GM_SUN,
        "Nominal solar mass parameter GM, m^3/s^2 (IAU 2015).",
    ),
    ("M_SUN", core::M_SUN, "Solar mass, kg: GM_SUN / G_NEWTON."),
    (
        "G0_PER_YEAR",
        core::G0_PER_YEAR,
        "Standard gravity in natural units of years: 1.0323 per year (ly/yr^2).",
    ),
    ("AU_LY", core::AU_LY, "Astronomical unit in light-years."),
    ("PARSEC_LY", core::PARSEC_LY, "Parsec in light-years."),
];

const DOC: &str = "\
Exact SI constants and the conversion between SI and natural units.

Geometry in srkit uses natural units with c = 1: pick any unit of time and measure lengths in the distance
light covers in it (seconds and light-seconds, Julian years and light-years). Constants carry their full
official precision; round them only when you show them.

Constants (floats, SI unless the name says otherwise):

";

pub(crate) fn module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let m = PyModule::new(py, "srkit.units")?;
    let mut doc = DOC.to_owned();
    let mut all = vec!["Units"];
    for (name, value, text) in CONSTANTS {
        m.add(*name, *value)?;
        doc.push_str(&format!("{name:<17}{text}\n"));
        all.push(name);
    }
    m.setattr("__doc__", doc)?;
    m.add_class::<Units>()?;
    m.setattr("__all__", all)?;
    Ok(m)
}
