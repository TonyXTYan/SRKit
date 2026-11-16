//! Hawking radiation from a non-rotating, uncharged black hole.
//!
//! The temperature and the entropy follow from the mass alone. The **power**, and with it the mass-loss rate
//! and the lifetime, depend on which particle species the hole can emit and on how much of each wave gets
//! out through the potential barrier around the hole (the greybody factors). There is no single closed form
//! for that, so every function that needs it takes an [`EmissionModel`], and the caller sees which model a
//! number comes from.
//!
//! All models here write the power as `P = α ħ c⁶/(G² M²)` with a dimensionless coefficient `α`
//! ([`EmissionModel::power_coefficient`]). The mass then falls as `dM/dt = −α ħ c⁴/(G² M²)` and, **if `α`
//! stays constant**, the hole is gone after `t = G² M³/(3 α ħ c⁴)`. In reality `α` grows as the hole shrinks
//! and heats up, because heavier particles come within reach, so a lifetime computed with the initial `α`
//! is an upper bound. It is a close one for a large hole, because most of the lifetime is spent near the
//! initial mass.
//!
//! The physics oracle (`oracle/oracle.csv`) quotes powers and lifetimes from
//! [`EmissionModel::PhotonBlackbody`], the textbook formula, and `oracle_check.py` marks those rows as
//! ambiguous because the rows do not name the model. For a hole of 10⁹ kg that model gives `3.6e14 W` and
//! 2,700 years, while the oracle's estimate for all Standard Model species
//! (`EmissionModel::MacGibbon { f: 12.0 }`) gives `5.8e16 W` and 16 years: for a small hole the model
//! matters by more than two orders of magnitude. A stellar or supermassive hole is far too cold to emit
//! anything but photons and gravitons, and there [`EmissionModel::PagePhotonsGravitons`] is the realistic
//! choice; it differs from the textbook formula by a factor of 1.8.
//!
//! Masses are in kg, temperatures in kelvin, powers in watts and times in seconds. Functions return NaN
//! outside the stated domain. The formulas are semiclassical and hold for masses far above the Planck mass
//! (`2.2e-8` kg). "Within a few ULP" below is the accuracy with which the stated formula is evaluated for
//! the given `α`.
//!
//! References: S. W. Hawking, *Nature* 248 (1974) 30 and *Commun. Math. Phys.* 43 (1975) 199; D. N. Page,
//! *Phys. Rev. D* 13 (1976) 198; J. H. MacGibbon and B. R. Webber, *Phys. Rev. D* 41 (1990) 3052;
//! J. H. MacGibbon, *Phys. Rev. D* 44 (1991) 376; J. D. Bekenstein, *Phys. Rev. D* 7 (1973) 2333.

use crate::{C2, C3};
use core::f64::consts::PI;
use libm::cbrt;
use srkit_core::units::{BOLTZMANN, C, G_NEWTON, HBAR};

/// Power coefficient `α = 1/(15360π) = 2.0723e-5` of [`EmissionModel::PhotonBlackbody`].
pub const PHOTON_BLACKBODY_COEFFICIENT: f64 = 1.0 / (15360.0 * PI);

/// Total power coefficient `α = 2.011e-4` that Page (1976) computed for a hole emitting only massless
/// particles: photons, gravitons and the two neutrino species (four kinds, counting antiparticles) then
/// thought to be massless.
pub const PAGE_1976_TOTAL_COEFFICIENT: f64 = 2.011e-4;

/// Fraction of Page's (1976) total power carried by photons, 16.7%.
pub const PAGE_1976_PHOTON_FRACTION: f64 = 0.167;

/// Fraction of Page's (1976) total power carried by gravitons, 1.9%.
pub const PAGE_1976_GRAVITON_FRACTION: f64 = 0.019;

/// Power coefficient `α = 3.740e-5` of [`EmissionModel::PagePhotonsGravitons`]: Page's total times the
/// photon and graviton fractions, the same product as in `oracle_check.py`.
pub const PAGE_PHOTONS_GRAVITONS_COEFFICIENT: f64 =
    PAGE_1976_TOTAL_COEFFICIENT * (PAGE_1976_PHOTON_FRACTION + PAGE_1976_GRAVITON_FRACTION);

/// MacGibbon's mass-loss normalisation, in kg³/s: `dM/dt = −5.34e16 f(M)/M²` kg/s with `M` in kg
/// (`5.34e25 f(M)/M²` g/s with `M` in grams in the original).
pub const MACGIBBON_MASS_LOSS: f64 = 5.34e16;

/// Which particles a hole is taken to emit, and how: the model behind a power, a mass-loss rate or a
/// lifetime.
///
/// Every model is a dimensionless coefficient `α` in `P = α ħ c⁶/(G² M²)`; see
/// [`EmissionModel::power_coefficient`].
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum EmissionModel {
    /// Photons only, radiated as by a perfect black body with the area of the horizon, `A = 16π G² M²/c⁴`, at
    /// the Hawking temperature: the Stefan-Boltzmann law `P = σ A T⁴`. No greybody factors and no other
    /// species. `α = 1/(15360π) = 2.072e-5`, which gives the familiar lifetime `5120 π G² M³/(ħ c⁴)`.
    ///
    /// This is the textbook estimate and the model behind the oracle rows SR-64, SR-65, DT-21 and DT-22. It
    /// is not a physical emission rate. The horizon area is not the emitting area that geometric optics
    /// gives (the capture cross-section for light, `27π G² M²/c⁴`, is 27/4 times larger), most of the
    /// spectrum has wavelengths comparable with the hole and is partly reflected by the potential barrier,
    /// and photons are never the only species.
    PhotonBlackbody,
    /// Photons and gravitons with Page's (1976) numerically computed greybody factors:
    /// `α = 2.011e-4 × (0.167 + 0.019) = 3.740e-5`, 1.80 times [`EmissionModel::PhotonBlackbody`].
    ///
    /// The realistic model for a hole too cold to emit any massive particle, with `k T` far below the
    /// lightest neutrino mass: any stellar or supermassive hole. If the lightest neutrino is massless it is
    /// emitted as well; Page's 81.4% for neutrinos is for two species with their antiparticles, so one
    /// species adds 40.7% of `2.011e-4`. Pass the total through [`EmissionModel::Coefficient`].
    PagePhotonsGravitons,
    /// MacGibbon's (1991) parametrisation `dM/dt = −5.34e16 f(M)/M²` kg/s (`M` in kg), where `f(M)` counts
    /// the emitted species: `α = 5.34e16 f G²/(ħ c⁴) = 2.793e-4 f`.
    ///
    /// `f` is normalised to 1 for a hole that emits only massless particles and rises as the hole shrinks
    /// and more species come within reach. This crate does not compute `f(M)`; the caller supplies it
    /// (dimensionless, `> 0`). The oracle's note for SR-64 uses the estimate `f ≈ 12` at 10⁹ kg. A lifetime
    /// from this model holds `f` at the given value.
    MacGibbon {
        /// The species factor `f(M)`, dimensionless, `> 0`.
        f: f64,
    },
    /// A caller-supplied coefficient `α` in `P = α ħ c⁶/(G² M²)`, dimensionless, `> 0`. Use it for any other
    /// emission calculation.
    Coefficient(f64),
}

impl EmissionModel {
    /// The dimensionless coefficient `α` in `P = α ħ c⁶/(G² M²)`.
    ///
    /// `2.072e-5` for [`EmissionModel::PhotonBlackbody`], `3.740e-5` for
    /// [`EmissionModel::PagePhotonsGravitons`], `2.793e-4 f` for [`EmissionModel::MacGibbon`], and the given
    /// value for [`EmissionModel::Coefficient`]. The functions of this module return NaN unless it is
    /// positive.
    pub fn power_coefficient(self) -> f64 {
        match self {
            Self::PhotonBlackbody => PHOTON_BLACKBODY_COEFFICIENT,
            Self::PagePhotonsGravitons => PAGE_PHOTONS_GRAVITONS_COEFFICIENT,
            Self::MacGibbon { f } => {
                MACGIBBON_MASS_LOSS * f * G_NEWTON * G_NEWTON / (HBAR * C2 * C2)
            }
            Self::Coefficient(alpha) => alpha,
        }
    }
}

/// The Hawking temperature, `T = ħ c³/(8π G M k_B)`, in kelvin.
///
/// * `mass`: kg, `> 0`.
///
/// The temperature of the thermal radiation as seen from far away. `6.17e-8 K` for one solar mass and
/// `1.23e14 K` for 10⁹ kg; it is inversely proportional to the mass. It does not depend on the emission
/// model. Within a few ULP.
///
/// Reference: Hawking (1974, 1975).
pub fn temperature(mass: f64) -> f64 {
    if mass > 0.0 {
        HBAR * C3 / (8.0 * PI * G_NEWTON * mass * BOLTZMANN)
    } else {
        f64::NAN
    }
}

/// The Bekenstein-Hawking entropy, `S = k_B c³ A/(4 G ħ) = 4π G M² k_B/(ħ c)`, in J/K.
///
/// * `mass`: kg, `≥ 0`.
///
/// One quarter of the horizon area `A = 16π G² M²/c⁴` in Planck units, times `k_B`. `1.45e54 J/K` for one
/// solar mass. Divide by `k_B ln 2` for bits. It does not depend on the emission model. Within a few ULP.
///
/// Reference: Bekenstein (1973); Hawking (1975).
pub fn entropy(mass: f64) -> f64 {
    if mass >= 0.0 {
        4.0 * PI * G_NEWTON * mass * mass * BOLTZMANN / (HBAR * C)
    } else {
        f64::NAN
    }
}

/// The radiated power, `P = α ħ c⁶/(G² M²)`, in watts.
///
/// * `mass`: kg, `> 0`.
/// * `model`: the emission model that sets `α`; see [`EmissionModel`].
///
/// For 10⁹ kg: `3.56e14 W` with [`EmissionModel::PhotonBlackbody`] (`P = ħ c⁶/(15360π G² M²)`, the oracle's
/// SR-64) and `5.76e16 W` with `MacGibbon { f: 12.0 }`. The power is inversely proportional to the square of
/// the mass. Within a few ULP.
///
/// Reference: Hawking (1975); Page (1976); MacGibbon (1991).
pub fn power(mass: f64, model: EmissionModel) -> f64 {
    let alpha = model.power_coefficient();
    if mass > 0.0 && alpha > 0.0 {
        let gm = G_NEWTON * mass;
        alpha * HBAR * C3 * C3 / (gm * gm)
    } else {
        f64::NAN
    }
}

/// The rate of mass loss, `−dM/dt = P/c² = α ħ c⁴/(G² M²)`, in kg/s (a positive number).
///
/// * `mass`: kg, `> 0`.
/// * `model`: the emission model that sets `α`; see [`EmissionModel`].
///
/// Within a few ULP.
///
/// Reference: Page (1976); MacGibbon (1991).
pub fn mass_loss_rate(mass: f64, model: EmissionModel) -> f64 {
    let alpha = model.power_coefficient();
    if mass > 0.0 && alpha > 0.0 {
        let gm = G_NEWTON * mass;
        alpha * HBAR * C2 * C2 / (gm * gm)
    } else {
        f64::NAN
    }
}

/// The time for a hole to evaporate completely at constant `α`, `t = G² M³/(3 α ħ c⁴)`, in seconds.
///
/// * `mass`: the initial mass, kg, `≥ 0`.
/// * `model`: the emission model that sets `α`; see [`EmissionModel`].
///
/// With [`EmissionModel::PhotonBlackbody`] this is `5120 π G² M³/(ħ c⁴)`: `2.10e67` years for one solar mass
/// and 2,670 years for 10⁹ kg (the oracle's DT-22 and SR-65). With [`EmissionModel::PagePhotonsGravitons`]
/// the solar-mass figure is `1.16e67` years. The lifetime is proportional to the cube of the mass.
///
/// `α` is held constant, which overestimates the lifetime once the hole is hot enough to emit further
/// species; see the module documentation. Accretion is ignored, and so is the cosmic background radiation,
/// which today is hotter than any hole above `4.5e22` kg and feeds it faster than it evaporates. Within a
/// few ULP.
///
/// Reference: Hawking (1975); Page (1976).
pub fn lifetime(mass: f64, model: EmissionModel) -> f64 {
    let alpha = model.power_coefficient();
    if mass >= 0.0 && alpha > 0.0 {
        let gm = G_NEWTON * mass;
        gm * gm * mass / (3.0 * alpha * HBAR * C2 * C2)
    } else {
        f64::NAN
    }
}

/// The initial mass whose lifetime at constant `α` is a given time, `M = (3 α ħ c⁴ t/G²)^(1/3)`, in kg.
///
/// * `lifetime`: s, `≥ 0`.
/// * `model`: the emission model that sets `α`; see [`EmissionModel`].
///
/// The inverse of [`lifetime`]. With [`EmissionModel::PhotonBlackbody`], a hole that lasts one year starts
/// at `7.2e7` kg and one that lasts 13.8 billion years at `1.7e11` kg. The same caveat applies: with
/// constant `α` the mass needed for a given lifetime is underestimated for small, hot holes. Within a few
/// ULP.
///
/// Reference: as for [`lifetime`].
pub fn mass_for_lifetime(lifetime: f64, model: EmissionModel) -> f64 {
    let alpha = model.power_coefficient();
    if lifetime >= 0.0 && alpha > 0.0 {
        cbrt(3.0 * alpha * HBAR * C2 * C2 * lifetime / (G_NEWTON * G_NEWTON))
    } else {
        f64::NAN
    }
}
