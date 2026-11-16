//! The function modules: `srkit.kinematics`, `srkit.optics`, `srkit.rocket` and `srkit.kepler`.

use pyo3::prelude::*;
use srkit_core::{
    frames, kepler as core_kepler,
    kinematics::{self as core_kin, Event, FourVector, Separation},
    optics as core_optics, rocket as core_rocket,
};

use crate::{
    convert::{fl, invalid, repr_fields, subluminal, tup, IntoPy, Tup3, V3},
    trip::exhaust,
};

// ------------------------------------------------------------------------------------------- kinematics

/// Lorentz factor of a velocity `beta` (a vector, as a fraction of c). Raises `ValueError` at or above c.
#[pyfunction]
fn gamma(beta: V3) -> PyResult<f64> {
    core_kin::gamma(beta.0).py()
}

/// Lorentz factor of a speed `0 <= b < 1`, as `1 / sqrt((1 - b)(1 + b))` so that it stays accurate near c.
#[pyfunction]
fn gamma_of_speed(b: f64) -> PyResult<f64> {
    core_kin::gamma_of_speed(b).py()
}

/// Rapidity of a speed: `atanh(b)`.
#[pyfunction]
fn rapidity_of_speed(b: f64) -> f64 {
    core_kin::rapidity_of_speed(b)
}

/// Speed of a rapidity, as a fraction of c: `tanh(phi)`.
#[pyfunction]
fn speed_of_rapidity(phi: f64) -> f64 {
    core_kin::speed_of_rapidity(phi)
}

/// `1 - beta` for a rapidity, without the cancellation of `1 - tanh(phi)`: `2 / (exp(2 phi) + 1)`.
#[pyfunction]
fn speed_shortfall(phi: f64) -> f64 {
    core_kin::speed_shortfall(phi)
}

/// Kinetic energy per unit rest energy, `gamma - 1`, from a rapidity, without cancellation:
/// `2 sinh(phi / 2)^2`. Multiply by `m c^2` for joules.
#[pyfunction]
fn kinetic_per_rest_energy(phi: f64) -> f64 {
    core_kin::kinetic_per_rest_energy(phi)
}

/// Boost a displacement `(dt, dx)` into the frame moving at velocity `beta` (a fraction of c) relative to this
/// one. Returns `(dt', dx')`. The inverse is the boost by `-beta`.
#[pyfunction]
fn boost(beta: V3, dt: f64, dx: V3) -> PyResult<(f64, Tup3)> {
    let out = FourVector::new(dt, dx.0).boost(beta.0).py()?;
    Ok((out.t, tup(out.x)))
}

/// Boost a displacement `(dt, dx)` into the frame whose velocity points along `rapidity` and whose rapidity is
/// that vector's length. Returns `(dt', dx')`.
///
/// Works on the null components, so light running along the boost keeps full precision at any Lorentz factor.
#[pyfunction]
fn boost_rapidity(rapidity: V3, dt: f64, dx: V3) -> (f64, Tup3) {
    let (t, x) = core_kin::boost_rapidity(rapidity.0, dt, dx.0);
    (t, tup(x))
}

/// Relativistic composition of velocities: the velocity, in the lab, of a body that moves at `u` in a frame
/// which itself moves at `v` through the lab. All are fractions of c.
#[pyfunction]
fn compose_velocities(v: V3, u: V3) -> PyResult<Tup3> {
    core_kin::compose_velocities(v.0, u.0).py().map(tup)
}

/// The velocity of a body moving at `u` in the lab, as measured in a frame moving at `v` through the lab.
#[pyfunction]
fn relative_velocity(v: V3, u: V3) -> PyResult<Tup3> {
    core_kin::relative_velocity(v.0, u.0).py().map(tup)
}

/// Received frequency divided by emitted frequency.
///
/// Both velocities are fractions of c in one inertial frame. `photon_direction` points along the light's
/// travel, from the emitter to the receiver, and need not be a unit vector.
#[pyfunction]
fn doppler_factor(
    emitter_velocity: V3,
    receiver_velocity: V3,
    photon_direction: V3,
) -> PyResult<f64> {
    core_kin::doppler_factor(emitter_velocity.0, receiver_velocity.0, photon_direction.0).py()
}

/// Head-on Doppler factor for a rapidity: `exp(phi)` approaching, `exp(-phi)` receding (pass a negative
/// rapidity).
#[pyfunction]
fn doppler_of_rapidity(phi: f64) -> f64 {
    core_kin::doppler_of_rapidity(phi)
}

/// The unit direction a photon travels in, as measured in a frame moving at velocity `beta`.
///
/// `direction` is the photon's direction of travel in the original frame. The direction on the sky towards the
/// source is the negative of both.
#[pyfunction]
fn aberrate(direction: V3, beta: V3) -> PyResult<Tup3> {
    core_kin::aberrate(direction.0, beta.0).py().map(tup)
}

/// The angle from the direction of motion at which a source appears to an observer of rapidity `phi`, given
/// the angle `theta` (radians) at which it lies in the rest frame: `tan(theta'/2) = exp(-phi) tan(theta/2)`.
#[pyfunction]
fn aberrate_angle(theta: f64, phi: f64) -> f64 {
    core_kin::aberrate_angle(theta, phi)
}

/// The invariant interval squared between two events, `dt^2 - |dx|^2` (positive when timelike).
#[pyfunction]
fn interval2(t_a: f64, x_a: V3, t_b: f64, x_b: V3) -> f64 {
    Event::new(t_a, x_a.0).interval2(Event::new(t_b, x_b.0))
}

/// How two events are related, decided in floating point: 'timelike', 'null' or 'spacelike'. For an exact
/// decision on stored integer times use `srkit.causal.classify`.
#[pyfunction]
fn separation(t_a: f64, x_a: V3, t_b: f64, x_b: V3) -> &'static str {
    match Event::new(t_a, x_a.0).separation(Event::new(t_b, x_b.0)) {
        Separation::Timelike => "timelike",
        Separation::Null => "null",
        Separation::Spacelike => "spacelike",
    }
}

/// Whether light (or anything slower) leaving event a can reach event b.
#[pyfunction]
fn can_reach(t_a: f64, x_a: V3, t_b: f64, x_b: V3) -> bool {
    Event::new(t_a, x_a.0).can_reach(Event::new(t_b, x_b.0))
}

/// Null coordinates along one axis: `(u, v) = (t - x, t + x)`.
#[pyfunction]
fn null_coordinates(t: f64, x: f64) -> (f64, f64) {
    frames::null_coordinates(t, x)
}

/// Back from null coordinates: `(t, x) = ((u + v) / 2, (v - u) / 2)`.
#[pyfunction]
fn from_null_coordinates(u: f64, v: f64) -> (f64, f64) {
    frames::from_null_coordinates(u, v)
}

const KINEMATICS_DOC: &str = "\
Boosts, velocity composition, rapidity, Doppler factors and aberration.

Natural units, c = 1: velocities are fractions of c, as `(x, y, z)`. The metric signature is (+ - - -): a
timelike interval is positive.
";

pub(crate) fn kinematics_module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let m = PyModule::new(py, "srkit.kinematics")?;
    m.setattr("__doc__", KINEMATICS_DOC)?;
    m.add_function(wrap_pyfunction!(gamma, &m)?)?;
    m.add_function(wrap_pyfunction!(gamma_of_speed, &m)?)?;
    m.add_function(wrap_pyfunction!(rapidity_of_speed, &m)?)?;
    m.add_function(wrap_pyfunction!(speed_of_rapidity, &m)?)?;
    m.add_function(wrap_pyfunction!(speed_shortfall, &m)?)?;
    m.add_function(wrap_pyfunction!(kinetic_per_rest_energy, &m)?)?;
    m.add_function(wrap_pyfunction!(boost, &m)?)?;
    m.add_function(wrap_pyfunction!(boost_rapidity, &m)?)?;
    m.add_function(wrap_pyfunction!(compose_velocities, &m)?)?;
    m.add_function(wrap_pyfunction!(relative_velocity, &m)?)?;
    m.add_function(wrap_pyfunction!(doppler_factor, &m)?)?;
    m.add_function(wrap_pyfunction!(doppler_of_rapidity, &m)?)?;
    m.add_function(wrap_pyfunction!(aberrate, &m)?)?;
    m.add_function(wrap_pyfunction!(aberrate_angle, &m)?)?;
    m.add_function(wrap_pyfunction!(interval2, &m)?)?;
    m.add_function(wrap_pyfunction!(separation, &m)?)?;
    m.add_function(wrap_pyfunction!(can_reach, &m)?)?;
    m.add_function(wrap_pyfunction!(null_coordinates, &m)?)?;
    m.add_function(wrap_pyfunction!(from_null_coordinates, &m)?)?;
    m.setattr(
        "__all__",
        [
            "aberrate",
            "aberrate_angle",
            "boost",
            "boost_rapidity",
            "can_reach",
            "compose_velocities",
            "doppler_factor",
            "doppler_of_rapidity",
            "from_null_coordinates",
            "gamma",
            "gamma_of_speed",
            "interval2",
            "kinetic_per_rest_energy",
            "null_coordinates",
            "rapidity_of_speed",
            "relative_velocity",
            "separation",
            "speed_of_rapidity",
            "speed_shortfall",
        ],
    )?;
    Ok(m)
}

// ----------------------------------------------------------------------------------------------- optics

/// How a source looks to a moving observer: what `observe` returns.
#[pyclass(frozen, get_all, module = "srkit.optics", name = "Seen")]
pub struct Seen {
    /// Doppler factor of the source's light and clock rate. Above 1 is a blueshift.
    doppler: f64,
    /// Angle of the source from the forward axis in the rest frame of the stars, in radians.
    theta: f64,
    /// Angle from the forward axis as the moving observer sees it, in radians.
    theta_seen: f64,
    /// Forward component of the seen unit direction towards the source.
    forward: f64,
    /// Rightward component of the seen unit direction towards the source.
    right: f64,
    /// Upward component of the seen unit direction towards the source.
    up: f64,
}

#[pymethods]
impl Seen {
    fn __repr__(&self, py: Python<'_>) -> String {
        repr_fields!(py, "Seen", self, doppler, theta, theta_seen, forward, right, up)
    }
}

/// Aberration and Doppler shift for an observer moving at rapidity `phi` along its forward axis.
///
/// `(a, br, bu)` is the direction towards the source in the rest frame of the stars, in the observer's
/// forward, right and up axes; it need not be a unit vector. With `aberration=False` the seen direction is
/// left as the rest-frame one (the Doppler factor is unchanged). Returns a `Seen`.
#[pyfunction]
#[pyo3(signature = (phi, a, br, bu, aberration=true))]
fn observe(phi: f64, a: f64, br: f64, bu: f64, aberration: bool) -> PyResult<Seen> {
    let s = core_optics::observe(phi, a, br, bu, aberration).py()?;
    Ok(Seen {
        doppler: s.doppler,
        theta: s.theta,
        theta_seen: s.theta_seen,
        forward: s.forward,
        right: s.right,
        up: s.up,
    })
}

/// The rate at which an observer moving at velocity `v` sees a source's clock run. `v_source` is the source's
/// velocity when the light left and `n` the unit vector from the observer to the source. Velocities are
/// fractions of c.
#[pyfunction]
fn playback(v: V3, n: V3, v_source: V3) -> PyResult<f64> {
    Ok(core_optics::playback(
        subluminal("v", v.0)?,
        n.0,
        subluminal("v_source", v_source.0)?,
    ))
}

/// An observer's axes from the unit vector `v` along its motion, as `(forward, up, right)`: forward along the
/// motion (the x axis for a zero vector), up as near the z axis as it can be while at right angles to forward,
/// right completing the set.
#[pyfunction]
fn basis(v: V3) -> (Tup3, Tup3, Tup3) {
    let (f, u, r) = core_optics::basis(v.0);
    (tup(f), tup(u), tup(r))
}

/// Smallest detail an aperture resolves at a distance (Rayleigh criterion): `1.22 wavelength / aperture *
/// distance`. All lengths in the same unit.
#[pyfunction]
fn rayleigh_resolution(aperture: f64, wavelength: f64, distance: f64) -> f64 {
    core_optics::rayleigh_resolution(aperture, wavelength, distance)
}

/// Power received between two circular apertures of diameters `d_t` and `d_r` in the far field (Friis, with
/// the aperture gains): `power (pi d_t d_r / (4 wavelength distance))^2`. Lengths in the same unit; the result
/// is in the unit of `power`.
#[pyfunction]
fn friis_received_power(power: f64, d_t: f64, d_r: f64, wavelength: f64, distance: f64) -> f64 {
    core_optics::friis_received_power(power, d_t, d_r, wavelength, distance)
}

/// Bits per second carried by received power `power_w` (watts) at wavelength `wavelength_m` (metres), at
/// `bits_per_photon` bits for each photon.
#[pyfunction]
#[pyo3(signature = (power_w, wavelength_m, bits_per_photon=1.0))]
fn photon_bit_rate(power_w: f64, wavelength_m: f64, bits_per_photon: f64) -> f64 {
    core_optics::photon_bit_rate(power_w, wavelength_m, bits_per_photon)
}

/// Luminosity in watts of a black-body sphere of radius `radius_m` (metres) at temperature `temperature_k`
/// (kelvin): `4 pi R^2 sigma T^4`.
#[pyfunction]
fn black_body_luminosity(radius_m: f64, temperature_k: f64) -> f64 {
    core_optics::black_body_luminosity(radius_m, temperature_k)
}

const OPTICS_DOC: &str = "\
What a moving observer sees: where a source appears, how its light is shifted and how fast its clock seems to
run; and the diffraction and link-budget formulas for telescopes and beams.

`boost_rapidity`, `doppler_factor`, `doppler_of_rapidity`, `aberrate` and `aberrate_angle` are the functions
of `srkit.kinematics`, listed here as well. Natural units, c = 1, except where an argument's name carries an
SI unit (`_m`, `_w`, `_k`).
";

pub(crate) fn optics_module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let m = PyModule::new(py, "srkit.optics")?;
    m.setattr("__doc__", OPTICS_DOC)?;
    m.add_class::<Seen>()?;
    m.add_function(wrap_pyfunction!(observe, &m)?)?;
    m.add_function(wrap_pyfunction!(playback, &m)?)?;
    m.add_function(wrap_pyfunction!(basis, &m)?)?;
    m.add_function(wrap_pyfunction!(rayleigh_resolution, &m)?)?;
    m.add_function(wrap_pyfunction!(friis_received_power, &m)?)?;
    m.add_function(wrap_pyfunction!(photon_bit_rate, &m)?)?;
    m.add_function(wrap_pyfunction!(black_body_luminosity, &m)?)?;
    m.setattr(
        "__all__",
        [
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
        ],
    )?;
    Ok(m)
}

/// The names `srkit.optics` shares with `srkit.kinematics`.
pub(crate) const SHARED_WITH_OPTICS: &[&str] = &[
    "aberrate",
    "aberrate_angle",
    "boost_rapidity",
    "doppler_factor",
    "doppler_of_rapidity",
];

// ----------------------------------------------------------------------------------------------- rocket

/// Check a mass ratio.
fn ratio(mass_ratio: f64) -> PyResult<f64> {
    if mass_ratio > 0.0 {
        Ok(mass_ratio)
    } else {
        Err(invalid(format!(
            "mass_ratio must be positive (start mass over end mass), got {}",
            fl(mass_ratio)
        )))
    }
}

/// Mass ratio (start over end) to change rapidity by `dphi` with exhaust speed `w`, a fraction of c (1 for a
/// photon rocket): `exp(dphi / w)`. Thrust in any direction counts.
#[pyfunction]
#[pyo3(signature = (dphi, w=1.0))]
fn mass_ratio(dphi: f64, w: f64) -> PyResult<f64> {
    Ok(core_rocket::mass_ratio(dphi, exhaust(w)?))
}

/// The change of rapidity a mass ratio buys: `w ln(mass_ratio)`.
#[pyfunction]
#[pyo3(signature = (mass_ratio, w=1.0))]
fn rapidity_of_mass_ratio(mass_ratio: f64, w: f64) -> PyResult<f64> {
    Ok(core_rocket::rapidity_of_mass_ratio(
        ratio(mass_ratio)?,
        exhaust(w)?,
    ))
}

/// Propellant mass per unit of dry mass: the mass ratio minus 1, accurate for small burns.
#[pyfunction]
#[pyo3(signature = (dphi, w=1.0))]
fn propellant_per_dry_mass(dphi: f64, w: f64) -> PyResult<f64> {
    Ok(core_rocket::propellant_per_dry_mass(dphi, exhaust(w)?))
}

/// Energy to make the propellant, in units of the dry mass's rest energy, where `efficiency` is the fraction
/// of input energy that ends up as propellant mass-energy. Multiply by `m_dry c^2` for joules.
#[pyfunction]
#[pyo3(signature = (dphi, w=1.0, efficiency=1.0))]
fn energy_per_dry_rest_energy(dphi: f64, w: f64, efficiency: f64) -> PyResult<f64> {
    if !(efficiency > 0.0 && efficiency <= 1.0) {
        return Err(invalid(format!(
            "efficiency must be above 0 and at most 1, got {}",
            fl(efficiency)
        )));
    }
    Ok(core_rocket::energy_per_dry_rest_energy(
        dphi,
        exhaust(w)?,
        efficiency,
    ))
}

/// Rapidity of a perfectly reflecting sail pushed from rest by a beam, whatever the beam's power history:
/// `ln(1 + 2 e)`, where `e` is the beam energy (measured in the frame the sail started in) over the sail's
/// rest energy.
#[pyfunction]
fn sail_rapidity(energy_per_rest_energy: f64) -> PyResult<f64> {
    if !(energy_per_rest_energy >= 0.0) {
        return Err(invalid(format!(
            "energy_per_rest_energy must not be negative, got {}",
            fl(energy_per_rest_energy)
        )));
    }
    Ok(core_rocket::sail_rapidity(energy_per_rest_energy))
}

/// Kinetic energy per unit rest energy at rapidity `phi`: `gamma - 1`, without cancellation.
#[pyfunction]
#[pyo3(name = "kinetic_per_rest_energy")]
fn rocket_kinetic_per_rest_energy(phi: f64) -> f64 {
    core_rocket::kinetic_per_rest_energy(phi)
}

const ROCKET_DOC: &str = "\
The relativistic rocket equation, light sails and what a trip costs in mass and energy.

Rapidities have no unit; exhaust speeds are fractions of c. Energies come out per unit of rest energy:
multiply by `m c^2` (with `srkit.units.C`) for joules.
";

pub(crate) fn rocket_module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let m = PyModule::new(py, "srkit.rocket")?;
    m.setattr("__doc__", ROCKET_DOC)?;
    m.add_function(wrap_pyfunction!(mass_ratio, &m)?)?;
    m.add_function(wrap_pyfunction!(rapidity_of_mass_ratio, &m)?)?;
    m.add_function(wrap_pyfunction!(propellant_per_dry_mass, &m)?)?;
    m.add_function(wrap_pyfunction!(energy_per_dry_rest_energy, &m)?)?;
    m.add_function(wrap_pyfunction!(sail_rapidity, &m)?)?;
    m.add_function(wrap_pyfunction!(rocket_kinetic_per_rest_energy, &m)?)?;
    m.setattr(
        "__all__",
        [
            "energy_per_dry_rest_energy",
            "kinetic_per_rest_energy",
            "mass_ratio",
            "propellant_per_dry_mass",
            "rapidity_of_mass_ratio",
            "sail_rapidity",
        ],
    )?;
    Ok(m)
}

// ----------------------------------------------------------------------------------------------- kepler

/// Check an eccentricity.
fn eccentricity(e: f64) -> PyResult<f64> {
    if (0.0..1.0).contains(&e) {
        Ok(e)
    } else {
        Err(invalid(format!(
            "e must be at least 0 and below 1 (a bound orbit), got {}",
            fl(e)
        )))
    }
}

/// Eccentric anomaly `E` solving Kepler's equation `E - e sin E = M` for mean anomaly `M` (radians) and
/// eccentricity `0 <= e < 1`.
#[pyfunction]
fn eccentric_anomaly(mean_anomaly: f64, e: f64) -> PyResult<f64> {
    Ok(core_kepler::eccentric_anomaly(
        mean_anomaly,
        eccentricity(e)?,
    ))
}

/// Position in the orbital plane, in units of the semi-major axis, with the focus at the origin and periapsis
/// along +x: `(cos E - e, sqrt(1 - e^2) sin E)`.
#[pyfunction]
fn position_in_plane(eccentric_anomaly: f64, e: f64) -> PyResult<(f64, f64)> {
    Ok(core_kepler::position_in_plane(
        eccentric_anomaly,
        eccentricity(e)?,
    ))
}

/// True anomaly (radians) from eccentric anomaly.
#[pyfunction]
fn true_anomaly(eccentric_anomaly: f64, e: f64) -> PyResult<f64> {
    Ok(core_kepler::true_anomaly(
        eccentric_anomaly,
        eccentricity(e)?,
    ))
}

const KEPLER_DOC: &str = "\
Kepler's equation, for bound two-body orbits drawn or used as landmarks. Angles are in radians.
";

pub(crate) fn kepler_module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let m = PyModule::new(py, "srkit.kepler")?;
    m.setattr("__doc__", KEPLER_DOC)?;
    m.add_function(wrap_pyfunction!(eccentric_anomaly, &m)?)?;
    m.add_function(wrap_pyfunction!(position_in_plane, &m)?)?;
    m.add_function(wrap_pyfunction!(true_anomaly, &m)?)?;
    m.setattr(
        "__all__",
        ["eccentric_anomaly", "position_in_plane", "true_anomaly"],
    )?;
    Ok(m)
}
