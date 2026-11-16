"""Type stubs for srkit.rocket: the relativistic rocket equation, light sails and what a trip costs."""

def mass_ratio(dphi: float, w: float = 1.0) -> float:
    """Mass ratio (start over end) to change rapidity by dphi with exhaust speed w: exp(dphi / w)."""

def rapidity_of_mass_ratio(mass_ratio: float, w: float = 1.0) -> float:
    """The change of rapidity a mass ratio buys: w ln(mass_ratio)."""

def propellant_per_dry_mass(dphi: float, w: float = 1.0) -> float:
    """Propellant mass per unit of dry mass: the mass ratio minus 1."""

def energy_per_dry_rest_energy(dphi: float, w: float = 1.0, efficiency: float = 1.0) -> float:
    """Energy to make the propellant, in units of the dry mass's rest energy."""

def sail_rapidity(energy_per_rest_energy: float) -> float:
    """Rapidity of a perfectly reflecting sail pushed from rest by a beam: ln(1 + 2 e)."""

def kinetic_per_rest_energy(phi: float) -> float:
    """Kinetic energy per unit rest energy at rapidity phi: gamma - 1."""

__all__ = [
    "energy_per_dry_rest_energy",
    "kinetic_per_rest_energy",
    "mass_ratio",
    "propellant_per_dry_mass",
    "rapidity_of_mass_ratio",
    "sail_rapidity",
]
