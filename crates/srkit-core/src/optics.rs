//! What a moving observer sees: where a source appears, how its light is shifted and how fast its clock seems to
//! run; and the diffraction and link-budget formulas for telescopes and beams.

use core::f64::consts::PI;

use crate::{units, vec3::Vec3, Error};

/// How a source looks to an observer moving at rapidity `phi` along the forward axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Seen {
    /// Doppler factor of the source's light and clock rate: `gamma (1 + beta cos theta)`, where `theta` is the
    /// rest-frame angle of the source from the direction of motion. Above 1 is a blueshift.
    pub doppler: f64,
    /// Angle of the source from the forward axis in the rest frame of the stars.
    pub theta: f64,
    /// Angle from the forward axis as the moving observer sees it: `tan(theta'/2) = e^(-phi) tan(theta/2)`.
    pub theta_seen: f64,
    /// The seen unit direction towards the source in the observer's forward, right and up axes.
    pub forward: f64,
    pub right: f64,
    pub up: f64,
}

/// Aberration and Doppler shift in a ship-aligned basis.
///
/// `(a, br, bu)` is the direction **towards the source** in the rest frame of the stars, in the observer's
/// forward, right and up axes; it need not be a unit vector. With `aberration = false` the seen direction is left
/// as the rest-frame one (the Doppler factor is unchanged), for a display that turns aberration off.
///
/// The Doppler factor is formed as `e^phi cos²(theta/2) + e^(-phi) sin²(theta/2)`, which is exact in form and
/// has no cancellation at any rapidity.
pub fn observe(phi: f64, a: f64, br: f64, bu: f64, aberration: bool) -> Result<Seen, Error> {
    if !phi.is_finite() || !a.is_finite() || !br.is_finite() || !bu.is_finite() {
        return Err(Error::InvalidInput);
    }
    let scale = a.abs().max(br.abs()).max(bu.abs());
    if scale == 0.0 {
        return Err(Error::InvalidInput);
    }
    let (a, br, bu) = (a / scale, br / scale, bu / scale);
    let transverse = libm::hypot(br, bu);
    let length = libm::hypot(a, transverse);
    let (cosine, sine) = (a / length, transverse / length);
    let theta = libm::atan2(transverse, a);
    // cos(theta/2) and sin(theta/2); the small one comes from sin(theta) rather than from 1 -+ cos(theta)
    let (ch, sh) = if a >= 0.0 {
        let ch = libm::sqrt(0.5 * (1.0 + cosine));
        (ch, sine / (2.0 * ch))
    } else {
        let sh = libm::sqrt(0.5 * (1.0 - cosine));
        (sine / (2.0 * sh), sh)
    };
    let weighted_exp = |exponent: f64, weight: f64| {
        if weight == 0.0 {
            0.0
        } else {
            libm::exp(exponent + 2.0 * libm::log(weight))
        }
    };
    let doppler = weighted_exp(phi, ch) + weighted_exp(-phi, sh);
    if !doppler.is_finite() {
        return Err(Error::Overflow);
    }
    if doppler == 0.0 {
        return Err(Error::PrecisionLoss);
    }
    let theta_seen = if !aberration || transverse == 0.0 {
        theta
    } else if phi >= 0.0 {
        2.0 * libm::atan2(libm::exp(-phi) * sh, ch)
    } else {
        2.0 * libm::atan2(sh, libm::exp(phi) * ch)
    };
    let forward = if transverse == 0.0 {
        cosine
    } else {
        libm::cos(theta_seen)
    };
    let projected = libm::sin(theta_seen);
    let (right, up) = if transverse > 0.0 {
        (br / transverse * projected, bu / transverse * projected)
    } else {
        (0.0, 0.0)
    };
    Ok(Seen {
        doppler,
        theta,
        theta_seen,
        forward,
        right,
        up,
    })
}

/// The rate at which an observer moving at velocity `v` sees a source's clock run, the source having had
/// velocity `v_source` when the light left, and `n` the unit vector from observer to source:
/// `sqrt((1 - v_source²) / (1 - v²)) (1 + n·v) / (1 + n·v_source)`.
pub fn playback(v: Vec3, n: Vec3, v_source: Vec3) -> f64 {
    libm::sqrt((1.0 - v_source.norm_squared()) / (1.0 - v.norm_squared())) * (1.0 + n.dot(v))
        / (1.0 + n.dot(v_source))
}

/// An observer's forward, up and right axes from its velocity: forward along the motion (`+x` at rest), up as
/// near `+z` as it can be while perpendicular to forward, right completing the set.
pub fn basis(v: Vec3) -> (Vec3, Vec3, Vec3) {
    let f = if v.norm() > 0.0 { v } else { Vec3::X };
    let mut u = Vec3::Z - f * f.z;
    if u.norm() < 1e-9 {
        u = Vec3::X;
    }
    let n = u.norm();
    let u = if n > 0.0 { u / n } else { u };
    (f, u, f.cross(u))
}

/// Smallest detail an aperture resolves at a distance (Rayleigh criterion): `1.22 lambda / D * distance`.
/// All lengths in the same unit.
pub fn rayleigh_resolution(aperture: f64, wavelength: f64, distance: f64) -> f64 {
    1.22 * wavelength / aperture * distance
}

/// Power received between two circular apertures of diameters `d_t` and `d_r` in the far field (Friis, with the
/// aperture gains): `P_t (pi d_t d_r / (4 lambda distance))²`. Lengths in the same unit; power in the unit given.
pub fn friis_received_power(power: f64, d_t: f64, d_r: f64, wavelength: f64, distance: f64) -> f64 {
    let x = PI * d_t * d_r / (4.0 * wavelength * distance);
    power * x * x
}

/// Bits per second carried by received power `power_w` (watts) at wavelength `wavelength_m` (metres) at
/// `bits_per_photon` bits for each photon: `P / (h c / lambda)` photons per second.
pub fn photon_bit_rate(power_w: f64, wavelength_m: f64, bits_per_photon: f64) -> f64 {
    power_w / (units::PLANCK * units::C / wavelength_m) * bits_per_photon
}

/// Luminosity in watts of a black-body sphere of radius `radius_m` (metres) at temperature `temperature_k`
/// (kelvin): `4 pi R² sigma T⁴`.
pub fn black_body_luminosity(radius_m: f64, temperature_k: f64) -> f64 {
    let t2 = temperature_k * temperature_k;
    4.0 * PI * radius_m * radius_m * units::STEFAN_BOLTZMANN * t2 * t2
}
