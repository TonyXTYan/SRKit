"""Type stubs for srkit.units: exact SI constants and the conversion between SI and natural units."""

from typing import ClassVar, final

C: float
"""Speed of light in vacuum, m/s. Exact."""
JULIAN_YEAR: float
"""Julian year, s: exactly 365.25 days of 86,400 s."""
DAY: float
"""Day, s."""
LIGHT_YEAR: float
"""Light-year, m: C * JULIAN_YEAR."""
AU: float
"""Astronomical unit, m. Exact (IAU 2012)."""
PARSEC: float
"""Parsec, m: 648000 / pi astronomical units (IAU 2015)."""
G0: float
"""Standard gravity, m/s^2. Exact."""
PLANCK: float
"""Planck constant, J s. Exact."""
HBAR: float
"""Reduced Planck constant, J s."""
BOLTZMANN: float
"""Boltzmann constant, J/K. Exact."""
G_NEWTON: float
"""Newtonian constant of gravitation, m^3 kg^-1 s^-2 (CODATA 2018)."""
STEFAN_BOLTZMANN: float
"""Stefan-Boltzmann constant, W m^-2 K^-4."""
L_SUN: float
"""Nominal solar luminosity, W (IAU 2015)."""
R_SUN: float
"""Nominal solar radius, m (IAU 2015)."""
GM_SUN: float
"""Nominal solar mass parameter GM, m^3/s^2 (IAU 2015)."""
M_SUN: float
"""Solar mass, kg: GM_SUN / G_NEWTON."""
G0_PER_YEAR: float
"""Standard gravity in natural units of years: 1.0323 per year (ly/yr^2)."""
AU_LY: float
"""Astronomical unit in light-years."""
PARSEC_LY: float
"""Parsec in light-years."""

@final
class Units:
    """A natural unit system: one unit of time, and the distance light covers in it as the unit of length."""

    SECONDS: ClassVar[Units]
    YEARS: ClassVar[Units]
    METRES: ClassVar[Units]
    def __new__(cls, time_s: float) -> Units: ...
    @property
    def time_s(self) -> float: ...
    @property
    def length_m(self) -> float: ...
    def time_from_si(self, seconds: float) -> float: ...
    def time_to_si(self, time: float) -> float: ...
    def length_from_si(self, metres: float) -> float: ...
    def length_to_si(self, length: float) -> float: ...
    def speed_from_si(self, metres_per_second: float) -> float: ...
    def speed_to_si(self, beta: float) -> float: ...
    def accel_from_si(self, metres_per_second2: float) -> float: ...
    def accel_to_si(self, alpha: float) -> float: ...
    def __eq__(self, other: object, /) -> bool: ...
    def __ne__(self, other: object, /) -> bool: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]

__all__ = [
    "Units",
    "C",
    "JULIAN_YEAR",
    "DAY",
    "LIGHT_YEAR",
    "AU",
    "PARSEC",
    "G0",
    "PLANCK",
    "HBAR",
    "BOLTZMANN",
    "G_NEWTON",
    "STEFAN_BOLTZMANN",
    "L_SUN",
    "R_SUN",
    "GM_SUN",
    "M_SUN",
    "G0_PER_YEAR",
    "AU_LY",
    "PARSEC_LY",
]
