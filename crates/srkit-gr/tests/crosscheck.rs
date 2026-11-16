//! Every public function against independent 120-digit references.
//!
//! `tests/data/crosscheck.csv` is written by `tools/crosscheck.py` (mpmath). Each row names a function, gives
//! its arguments as exact binary64 values, the reference value and a relative tolerance. The tolerance is
//! `2e-15` ("a few ULP") unless the quantity is ill conditioned at that point, in which case the script
//! multiplies it by the condition number and says why.

mod common;

use common::{csv_fields, relative_error};
use srkit_core::units;
use srkit_gr::hawking::{self, EmissionModel};
use srkit_gr::kerr::{self, Orbit};
use srkit_gr::{newtonian, schwarzschild, weak_field};
use std::collections::{BTreeMap, BTreeSet};

const REFERENCE: &str = include_str!("data/crosscheck.csv");

/// Every function the reference file must cover: the whole public API.
const FUNCTIONS: &[&str] = &[
    "gravitational_parameter",
    "gravitational_radius",
    "gravitational_time",
    "mass_from_gravitational_parameter",
    "schwarzschild::schwarzschild_radius",
    "schwarzschild::photon_sphere_radius",
    "schwarzschild::marginally_bound_radius",
    "schwarzschild::isco_radius",
    "schwarzschild::static_dtau_dt",
    "schwarzschild::gravitational_redshift",
    "schwarzschild::circular_orbit_dt_dtau",
    "schwarzschild::circular_orbit_angular_velocity",
    "kerr::outer_horizon_radius",
    "kerr::inner_horizon_radius",
    "kerr::ergosphere_radius",
    "kerr::horizon_angular_velocity",
    "kerr::irreducible_mass",
    "kerr::extractable_energy_fraction",
    "kerr::static_dtau_dt",
    "kerr::zamo_angular_velocity",
    "kerr::zamo_dtau_dt",
    "kerr::photon_orbit_radius",
    "kerr::marginally_bound_radius",
    "kerr::isco_radius",
    "kerr::isco_specific_energy",
    "kerr::isco_specific_angular_momentum",
    "kerr::isco_dt_dtau",
    "kerr::accretion_efficiency",
    "kerr::circular_orbit_dt_dtau",
    "kerr::circular_orbit_specific_energy",
    "kerr::circular_orbit_specific_angular_momentum",
    "kerr::circular_orbit_angular_velocity",
    "weak_field::shapiro_delay_leg",
    "weak_field::shapiro_delay",
    "weak_field::shapiro_delay_between",
    "weak_field::deflection_angle",
    "weak_field::einstein_angle",
    "weak_field::static_rate_offset",
    "weak_field::circular_orbit_rate_offset",
    "weak_field::orbit_vs_surface_rate_offset",
    "hawking::power_coefficient",
    "hawking::temperature",
    "hawking::entropy",
    "hawking::power",
    "hawking::mass_loss_rate",
    "hawking::lifetime",
    "hawking::mass_for_lifetime",
    "newtonian::binding_energy_uniform_sphere",
    "newtonian::circular_speed",
    "newtonian::escape_speed",
    "newtonian::circular_period",
    "newtonian::synodic_period",
    "newtonian::hohmann_transfer_time",
    "newtonian::hohmann_departure",
    "newtonian::hohmann_arrival",
];

fn orbit(code: f64) -> Orbit {
    if code == 1.0 {
        Orbit::Prograde
    } else if code == -1.0 {
        Orbit::Retrograde
    } else {
        panic!("orbit code {code}")
    }
}

fn model(code: f64, parameter: f64) -> EmissionModel {
    match code as i32 {
        0 => EmissionModel::PhotonBlackbody,
        1 => EmissionModel::PagePhotonsGravitons,
        2 => EmissionModel::MacGibbon { f: parameter },
        3 => EmissionModel::Coefficient(parameter),
        _ => panic!("model code {code}"),
    }
}

fn evaluate(function: &str, x: &[f64; 4]) -> f64 {
    match function {
        "gravitational_parameter" => srkit_gr::gravitational_parameter(x[0]),
        "gravitational_radius" => srkit_gr::gravitational_radius(x[0]),
        "gravitational_time" => srkit_gr::gravitational_time(x[0]),
        "mass_from_gravitational_parameter" => srkit_gr::mass_from_gravitational_parameter(x[0]),

        "schwarzschild::schwarzschild_radius" => schwarzschild::schwarzschild_radius(x[0]),
        "schwarzschild::photon_sphere_radius" => schwarzschild::photon_sphere_radius(x[0]),
        "schwarzschild::marginally_bound_radius" => schwarzschild::marginally_bound_radius(x[0]),
        "schwarzschild::isco_radius" => schwarzschild::isco_radius(x[0]),
        "schwarzschild::static_dtau_dt" => schwarzschild::static_dtau_dt(x[0], x[1]),
        "schwarzschild::gravitational_redshift" => {
            schwarzschild::gravitational_redshift(x[0], x[1], x[2])
        }
        "schwarzschild::circular_orbit_dt_dtau" => {
            schwarzschild::circular_orbit_dt_dtau(x[0], x[1])
        }
        "schwarzschild::circular_orbit_angular_velocity" => {
            schwarzschild::circular_orbit_angular_velocity(x[0], x[1])
        }

        "kerr::outer_horizon_radius" => kerr::outer_horizon_radius(x[0], x[1]),
        "kerr::inner_horizon_radius" => kerr::inner_horizon_radius(x[0], x[1]),
        "kerr::ergosphere_radius" => kerr::ergosphere_radius(x[0], x[1], x[2]),
        "kerr::horizon_angular_velocity" => kerr::horizon_angular_velocity(x[0], x[1]),
        "kerr::irreducible_mass" => kerr::irreducible_mass(x[0], x[1]),
        "kerr::extractable_energy_fraction" => kerr::extractable_energy_fraction(x[0]),
        "kerr::static_dtau_dt" => kerr::static_dtau_dt(x[0], x[1], x[2], x[3]),
        "kerr::zamo_angular_velocity" => kerr::zamo_angular_velocity(x[0], x[1], x[2], x[3]),
        "kerr::zamo_dtau_dt" => kerr::zamo_dtau_dt(x[0], x[1], x[2], x[3]),
        "kerr::photon_orbit_radius" => kerr::photon_orbit_radius(x[0], x[1], orbit(x[2])),
        "kerr::marginally_bound_radius" => kerr::marginally_bound_radius(x[0], x[1], orbit(x[2])),
        "kerr::isco_radius" => kerr::isco_radius(x[0], x[1], orbit(x[2])),
        "kerr::isco_specific_energy" => kerr::isco_specific_energy(x[0], orbit(x[1])),
        "kerr::isco_specific_angular_momentum" => {
            kerr::isco_specific_angular_momentum(x[0], x[1], orbit(x[2]))
        }
        "kerr::isco_dt_dtau" => kerr::isco_dt_dtau(x[0], orbit(x[1])),
        "kerr::accretion_efficiency" => kerr::accretion_efficiency(x[0], orbit(x[1])),
        "kerr::circular_orbit_dt_dtau" => {
            kerr::circular_orbit_dt_dtau(x[0], x[1], orbit(x[2]), x[3])
        }
        "kerr::circular_orbit_specific_energy" => {
            kerr::circular_orbit_specific_energy(x[0], x[1], orbit(x[2]), x[3])
        }
        "kerr::circular_orbit_specific_angular_momentum" => {
            kerr::circular_orbit_specific_angular_momentum(x[0], x[1], orbit(x[2]), x[3])
        }
        "kerr::circular_orbit_angular_velocity" => {
            kerr::circular_orbit_angular_velocity(x[0], x[1], orbit(x[2]), x[3])
        }

        "weak_field::shapiro_delay_leg" => weak_field::shapiro_delay_leg(x[0], x[1], x[2]),
        "weak_field::shapiro_delay" => weak_field::shapiro_delay(x[0], x[1], x[2], x[3]),
        "weak_field::shapiro_delay_between" => {
            weak_field::shapiro_delay_between(x[0], x[1], x[2], x[3])
        }
        "weak_field::deflection_angle" => weak_field::deflection_angle(x[0], x[1]),
        "weak_field::einstein_angle" => weak_field::einstein_angle(x[0], x[1], x[2], x[3]),
        "weak_field::static_rate_offset" => weak_field::static_rate_offset(x[0], x[1]),
        "weak_field::circular_orbit_rate_offset" => {
            weak_field::circular_orbit_rate_offset(x[0], x[1])
        }
        "weak_field::orbit_vs_surface_rate_offset" => {
            weak_field::orbit_vs_surface_rate_offset(x[0], x[1], x[2], x[3])
        }

        "hawking::power_coefficient" => model(x[0], x[1]).power_coefficient(),
        "hawking::temperature" => hawking::temperature(x[0]),
        "hawking::entropy" => hawking::entropy(x[0]),
        "hawking::power" => hawking::power(x[0], model(x[1], x[2])),
        "hawking::mass_loss_rate" => hawking::mass_loss_rate(x[0], model(x[1], x[2])),
        "hawking::lifetime" => hawking::lifetime(x[0], model(x[1], x[2])),
        "hawking::mass_for_lifetime" => hawking::mass_for_lifetime(x[0], model(x[1], x[2])),

        "newtonian::binding_energy_uniform_sphere" => {
            newtonian::binding_energy_uniform_sphere(x[0], x[1])
        }
        "newtonian::circular_speed" => newtonian::circular_speed(x[0], x[1]),
        "newtonian::escape_speed" => newtonian::escape_speed(x[0], x[1]),
        "newtonian::circular_period" => newtonian::circular_period(x[0], x[1]),
        "newtonian::synodic_period" => newtonian::synodic_period(x[0], x[1]),
        "newtonian::hohmann_transfer_time" => newtonian::hohmann_transfer_time(x[0], x[1], x[2]),
        "newtonian::hohmann_departure" => newtonian::hohmann_delta_v(x[0], x[1], x[2]).departure,
        "newtonian::hohmann_arrival" => newtonian::hohmann_delta_v(x[0], x[1], x[2]).arrival,

        _ => panic!("the reference file names an unknown function: {function}"),
    }
}

fn constant(name: &str) -> f64 {
    match name {
        "C" => units::C,
        "G_NEWTON" => units::G_NEWTON,
        "HBAR" => units::HBAR,
        "BOLTZMANN" => units::BOLTZMANN,
        "GM_SUN" => units::GM_SUN,
        "M_SUN" => units::M_SUN,
        "AU" => units::AU,
        "R_SUN" => units::R_SUN,
        "JULIAN_YEAR" => units::JULIAN_YEAR,
        _ => panic!("the reference file names an unknown constant: {name}"),
    }
}

#[test]
fn every_function_matches_its_high_precision_reference() {
    let mut worst: BTreeMap<String, (f64, usize)> = BTreeMap::new();
    let mut failures = Vec::new();
    let mut constants = 0;
    for record in REFERENCE
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .skip(1)
    {
        let fields = csv_fields(record);
        assert_eq!(fields.len(), 7, "bad reference record: {record}");
        if let Some(name) = fields[0].strip_prefix("constant:") {
            // The references were computed with these constants; they must be the crate's, bit for bit.
            let value: f64 = fields[1].parse().unwrap();
            assert_eq!(
                value.to_bits(),
                constant(name).to_bits(),
                "srkit_core::units::{name} differs from the value the references were computed with: \
                 update and rerun tools/crosscheck.py"
            );
            constants += 1;
            continue;
        }
        let mut x = [f64::NAN; 4];
        for (slot, field) in x.iter_mut().zip(&fields[1..5]) {
            if !field.is_empty() {
                *slot = field
                    .parse()
                    .unwrap_or_else(|_| panic!("bad argument in: {record}"));
            }
        }
        let want: f64 = fields[5]
            .parse()
            .unwrap_or_else(|_| panic!("bad reference in: {record}"));
        let tol: f64 = fields[6]
            .parse()
            .unwrap_or_else(|_| panic!("bad tolerance in: {record}"));
        let got = evaluate(&fields[0], &x);
        // A reference of exactly zero must be reproduced exactly; relative_error is then absolute.
        let error = relative_error(got, want);
        let allowed = if want == 0.0 { 0.0 } else { tol };
        if error.is_nan() || error > allowed {
            failures.push(format!("{record}\n    got {got:e}, error {error:e}"));
        }
        let entry = worst.entry(fields[0].clone()).or_insert((0.0, 0));
        entry.0 = entry.0.max(if want == 0.0 { 0.0 } else { error / tol });
        entry.1 += 1;
    }
    assert_eq!(constants, 9, "the reference file should pin nine constants");

    println!("{:<50} {:>5}  worst error / tolerance", "function", "rows");
    for (function, (ratio, rows)) in &worst {
        println!("{function:<50} {rows:>5}  {ratio:.3}");
    }
    assert!(
        failures.is_empty(),
        "{} reference rows failed:\n{}",
        failures.len(),
        failures.join("\n")
    );

    let covered: BTreeSet<&str> = worst.keys().map(String::as_str).collect();
    let expected: BTreeSet<&str> = FUNCTIONS.iter().copied().collect();
    assert_eq!(
        covered, expected,
        "the reference file does not cover exactly the public API"
    );
}
