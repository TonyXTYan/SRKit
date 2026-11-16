"""Type stubs for srkit.optics: what a moving observer sees, and telescope and beam formulas."""

from typing import Protocol, SupportsFloat, final

from .kinematics import aberrate as aberrate
from .kinematics import aberrate_angle as aberrate_angle
from .kinematics import boost_rapidity as boost_rapidity
from .kinematics import doppler_factor as doppler_factor
from .kinematics import doppler_of_rapidity as doppler_of_rapidity

_Vec3 = tuple[float, float, float]

class _Vec3Like(Protocol):
    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> SupportsFloat: ...

@final
class Seen:
    """How a source looks to a moving observer."""

    @property
    def doppler(self) -> float: ...
    @property
    def theta(self) -> float: ...
    @property
    def theta_seen(self) -> float: ...
    @property
    def forward(self) -> float: ...
    @property
    def right(self) -> float: ...
    @property
    def up(self) -> float: ...

def observe(phi: float, a: float, br: float, bu: float, aberration: bool = True) -> Seen:
    """Aberration and Doppler shift for an observer moving at rapidity phi along its forward axis."""

def playback(v: _Vec3Like, n: _Vec3Like, v_source: _Vec3Like) -> float:
    """The rate at which an observer moving at velocity v sees a source's clock run."""

def basis(v: _Vec3Like) -> tuple[_Vec3, _Vec3, _Vec3]:
    """An observer's axes from the unit vector along its motion, as (forward, up, right)."""

def rayleigh_resolution(aperture: float, wavelength: float, distance: float) -> float:
    """Smallest detail an aperture resolves at a distance (Rayleigh criterion)."""

def friis_received_power(
    power: float, d_t: float, d_r: float, wavelength: float, distance: float
) -> float:
    """Power received between two circular apertures in the far field."""

def photon_bit_rate(power_w: float, wavelength_m: float, bits_per_photon: float = 1.0) -> float:
    """Bits per second carried by received power (watts) at a wavelength (metres)."""

def black_body_luminosity(radius_m: float, temperature_k: float) -> float:
    """Luminosity in watts of a black-body sphere."""

__all__ = [
    "Seen",
    "aberrate",
    "aberrate_angle",
    "basis",
    "black_body_luminosity",
    "boost_rapidity",
    "doppler_factor",
    "doppler_of_rapidity",
    "friis_received_power",
    "observe",
    "photon_bit_rate",
    "playback",
    "rayleigh_resolution",
]
