//! Rapidity-based optical and frame transforms for the Navigator.
//!
//! These local floating-point constructions use `c = 1`: time and distance
//! coordinates must have matching units (the page uses Julian years and ly).
use crate::{kinematics::Vec3, Error};

/// Doppler and aberration in a ship-aligned orthonormal sky basis.
///
/// The input direction points **toward the source**: `a` is forward, `br`
/// right and `bu` up. Returns `[D, theta, theta_ship, forward, right, up]`.
/// Turning aberration off changes the displayed direction only, not Doppler.
pub fn observe(phi: f64, a: f64, br: f64, bu: f64, aberr: bool) -> Result<[f64; 6], Error> {
    if !phi.is_finite() || !a.is_finite() || !br.is_finite() || !bu.is_finite() {
        return Err(Error::InvalidInput);
    }
    // Scaling also keeps subnormal directions and finite very large directions
    // well conditioned. The page normally supplies a unit direction.
    let scale = a.abs().max(br.abs()).max(bu.abs());
    if scale == 0.0 {
        return Err(Error::InvalidInput);
    }
    let a = a / scale;
    let br = br / scale;
    let bu = bu / scale;
    let transverse = libm::hypot(br, bu);
    let length = libm::hypot(a, transverse);
    let cosine = a / length;
    let sine = transverse / length;
    let theta = libm::atan2(transverse, a);
    // Recover the small half-angle component from sin(theta), rather than
    // subtracting a cosine near ±1 or evaluating cos(pi/2) approximately.
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
    let theta_ship = if !aberr || transverse == 0.0 {
        theta
    } else if phi >= 0.0 {
        2.0 * libm::atan2(libm::exp(-phi) * sh, ch)
    } else {
        2.0 * libm::atan2(sh, libm::exp(phi) * ch)
    };
    let forward = if transverse == 0.0 {
        cosine
    } else {
        libm::cos(theta_ship)
    };
    let projected = libm::sin(theta_ship);
    let (right, up) = if transverse > 0.0 {
        (br / transverse * projected, bu / transverse * projected)
    } else {
        (0.0, 0.0)
    };
    Ok([doppler, theta, theta_ship, forward, right, up])
}

/// Boost a local displacement into the frame specified by a rapidity vector.
///
/// The vector points along the frame's velocity, and its norm is rapidity.
/// Null-coordinate scaling avoids subtracting gamma-sized terms for light
/// travelling exactly along the boost. This is an inexact construction, not a
/// certified causal predicate or an arbitrary-precision transformation.
pub fn boost(rapidity: Vec3, dt: f64, dx: Vec3) -> Result<(f64, Vec3), Error> {
    if !rapidity.is_finite() || !dt.is_finite() || !dx.is_finite() {
        return Err(Error::InvalidInput);
    }
    let phi = rapidity.norm();
    if !phi.is_finite() {
        return Err(Error::Overflow);
    }
    if phi == 0.0 {
        return Ok((dt, dx));
    }
    let axis = rapidity / phi;
    let parallel = axis.dot(dx);
    // Half before adding/subtracting to avoid needless intermediate overflow.
    let minus = (0.5 * dt - 0.5 * parallel) * libm::exp(phi);
    let plus = (0.5 * dt + 0.5 * parallel) * libm::exp(-phi);
    let time = minus + plus;
    let projected = plus - minus;
    let position = (dx - axis * parallel) + axis * projected;
    if !time.is_finite() || !position.is_finite() {
        return Err(Error::Overflow);
    }
    Ok((time, position))
}
