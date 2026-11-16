//! The relativistic rocket equation, light sails and what a trip costs in mass and energy.

/// Mass ratio (start over end) to change rapidity by `dphi` with exhaust speed `w` (a fraction of `c`; 1 for a
/// photon rocket): `exp(dphi / w)`. Thrust in any direction counts: a circle held at proper acceleration `a` for
/// proper time `s` costs `dphi = a s`, the same as a straight burn.
pub fn mass_ratio(dphi: f64, w: f64) -> f64 {
    libm::exp(dphi / w)
}

/// The change of rapidity a mass ratio buys: `w ln(MR)`.
pub fn rapidity_of_mass_ratio(mass_ratio: f64, w: f64) -> f64 {
    w * libm::log(mass_ratio)
}

/// Propellant mass per unit of dry mass: `MR - 1`, accurate for small burns.
pub fn propellant_per_dry_mass(dphi: f64, w: f64) -> f64 {
    libm::expm1(dphi / w)
}

/// Energy to make the propellant, in units of the dry mass's rest energy: `(MR - 1) / efficiency`, where
/// `efficiency` is the fraction of input energy that ends up as propellant mass-energy. Multiply by
/// `m_dry c²` for joules.
pub fn energy_per_dry_rest_energy(dphi: f64, w: f64, efficiency: f64) -> f64 {
    libm::expm1(dphi / w) / efficiency
}

/// Rapidity of a perfectly reflecting sail of rest energy `m c²` pushed from rest by beam energy `e` (measured
/// in the frame it started in), whatever the beam's power history: `ln(1 + 2 e / (m c²))`. Pass the ratio
/// `e / (m c²)`.
pub fn sail_rapidity(energy_per_rest_energy: f64) -> f64 {
    libm::log1p(2.0 * energy_per_rest_energy)
}

/// Kinetic energy per unit rest energy at rapidity `phi`: `gamma - 1`.
pub fn kinetic_per_rest_energy(phi: f64) -> f64 {
    crate::kinematics::kinetic_per_rest_energy(phi)
}
