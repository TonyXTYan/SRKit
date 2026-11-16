//! Golden values: every public function against 50-digit mpmath quadrature, and the background against
//! astropy. The reference file is written by `tools/crosscheck.py`; nothing in it comes from this crate.

mod common;

use common::{csv_fields, relative_error};
use srkit_flrw::Cosmology;
use std::collections::BTreeMap;

const GOLDEN: &str = include_str!("golden.csv");

/// Evaluate one named quantity. `x` holds the row's arguments in order.
fn evaluate(u: &Cosmology, quantity: &str, x: &[f64]) -> f64 {
    let some = |value: Option<f64>| value.unwrap_or(f64::NAN);
    match quantity {
        "age" => u.age(),
        "h_lambda" => u.h_lambda(),
        "scale_factor_at_equality" => u.scale_factor_at_equality(),
        "scale_factor" => u.scale_factor(x[0]),
        "ln_scale_factor" => u.ln_scale_factor(x[0]),
        "time_at_scale_factor" => u.time_at_scale_factor(x[0]),
        "time_at_ln_scale_factor" => u.time_at_ln_scale_factor(x[0]),
        "hubble" => u.hubble(x[0]),
        "hubble_at_scale_factor" => u.hubble_at_scale_factor(x[0]),
        "time_at_redshift" => u.time_at_redshift(x[0]),
        "lookback_time" => u.lookback_time(x[0]),
        "redshift_between" => u.redshift_between(x[0], x[1]),
        "conformal_time" => u.conformal_time(x[0]),
        "conformal_time_at_scale_factor" => u.conformal_time_at_scale_factor(x[0]),
        "conformal_time_limit" => u.conformal_time_limit(),
        "conformal_time_remaining" => u.conformal_time_remaining(x[0]),
        "time_at_conformal_time" => some(u.time_at_conformal_time(x[0])),
        "time_at_conformal_time_remaining" => some(u.time_at_conformal_time_remaining(x[0])),
        "comoving_distance" => u.comoving_distance(x[0], x[1]),
        "comoving_distance_at_redshift" => u.comoving_distance_at_redshift(x[0]),
        "particle_horizon_proper" => u.particle_horizon_proper(x[0]),
        "event_horizon_proper" => u.event_horizon_proper(x[0]),
        "reachable_share_of_observable" => u.reachable_share_of_observable(x[0]),
        "reach_volume_ratio" => u.reach_volume_ratio(x[0], x[1]),
        "reach_mass" => u.reach_mass(x[0]),
        "last_contact_time" => some(u.last_contact_time(x[0])),
        "last_reply_time" => some(u.last_reply_time(x[0])),
        "light_arrival_time" => some(u.light_arrival_time(x[0], x[1])),
        "light_emission_time" => some(u.light_emission_time(x[0], x[1])),
        "light_reception" => some(u.light_reception(x[0], x[1]).map(|r| r.time)),
        "coasting_proper_velocity" => u.coasting_proper_velocity(x[0], x[1], x[2]),
        "coasting_rapidity" => u.coasting_rapidity(x[0], x[1], x[2]),
        "coasting_gamma" => u.coasting_gamma(x[0], x[1], x[2]),
        "scale_factor_integral" => u.scale_factor_integral(x[0]),
        "thrust_proper_velocity" => u.thrust_proper_velocity(x[0], x[1], x[2], x[3]),
        "drag_limited_gamma" => u.drag_limited_gamma(x[0]),
        "de_sitter_thrust_proper_velocity" => u.de_sitter_thrust_proper_velocity(x[0], x[1], x[2]),
        "de_sitter_thrust_proper_time" => u.de_sitter_thrust_proper_time(x[0], x[1], x[2]),
        "de_sitter_temperature" => u.de_sitter_temperature(),
        "cmb_temperature" => u.cmb_temperature(x[0], x[1]),
        "cmb_crossover_scale_factor" => u.cmb_crossover_scale_factor(x[0]),
        "cmb_crossover_time" => u.cmb_crossover_time(x[0]),
        "lambda_tidal_acceleration" => u.lambda_tidal_acceleration(x[0]),
        "lambda_tidal_acceleration_si" => u.lambda_tidal_acceleration_si(x[0]),
        "last_light_gain" => u.last_light_gain(x[0]),
        "scale_factor_overflow_time" => u.scale_factor_overflow_time(),
        other => panic!("golden.csv names a quantity this test does not know: {other}"),
    }
}

/// Worst errors seen for one quantity from one source.
#[derive(Default)]
struct Worst {
    rows: usize,
    /// Worst relative error among the rows held to 1e-12: ordinary conditioning.
    error: f64,
    /// Rows whose tolerance follows a documented loss of digits (close epochs, the far future, a target at
    /// the edge of the horizon), and the worst error among them as a fraction of that tolerance.
    conditioned: usize,
    conditioned_ratio: f64,
}

#[test]
fn every_function_matches_the_reference() {
    let mut worst: BTreeMap<(String, String), Worst> = BTreeMap::new();
    let mut failures = Vec::new();
    let mut rows = 0;
    for line in GOLDEN.lines() {
        if line.is_empty() || line.starts_with('#') || line.starts_with("source,") {
            continue;
        }
        let f = csv_fields(line);
        assert_eq!(f.len(), 10, "malformed golden row: {line}");
        let number = |s: &str| -> f64 {
            s.parse()
                .unwrap_or_else(|_| panic!("bad number in: {line}"))
        };
        let u = Cosmology::per_year(number(&f[2]), number(&f[3])).expect("valid cosmology");
        let x: Vec<f64> = f[4..8]
            .iter()
            .filter(|s| !s.is_empty())
            .map(|s| number(s))
            .collect();
        let expected = number(&f[8]);
        let tolerance = number(&f[9]);
        let got = evaluate(&u, &f[1], &x);
        // A logarithm is compared on an absolute scale near zero: ln a = 0 today.
        let error = if f[1] == "ln_scale_factor" {
            (got - expected).abs() / (1.0 + expected.abs())
        } else {
            relative_error(got, expected)
        };
        let entry = worst.entry((f[0].clone(), f[1].clone())).or_default();
        entry.rows += 1;
        if tolerance == 1e-12 {
            entry.error = entry.error.max(error);
        } else {
            entry.conditioned += 1;
            entry.conditioned_ratio = entry.conditioned_ratio.max(error / tolerance);
        }
        if error.is_nan() || error > tolerance {
            failures.push(format!(
                "{line}\n    got {got:e}, relative error {error:.3e}"
            ));
        }
        rows += 1;
    }
    println!("{rows} golden rows");
    println!(
        "{:<8} {:<34} {:>5} {:>14} {:>12} {:>14}",
        "source", "quantity", "rows", "worst at 1e-12", "conditioned", "worst / tol"
    );
    let mut overall: f64 = 0.0;
    for ((source, quantity), w) in &worst {
        println!(
            "{:<8} {:<34} {:>5} {:>14.2e} {:>12} {:>14.2}",
            source, quantity, w.rows, w.error, w.conditioned, w.conditioned_ratio
        );
        overall = overall.max(w.error);
    }
    println!("worst relative error among rows held to 1e-12: {overall:.2e}");
    assert!(
        failures.is_empty(),
        "{} of {rows} golden rows failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
    assert!(rows > 2000, "golden.csv looks truncated: {rows} rows");
}

/// The presets hold exactly the numbers their documentation states, and the reference file was generated
/// for the same binary64 values.
#[test]
fn presets_are_the_documented_numbers() {
    let years = Cosmology::PLANCK_2018;
    assert_eq!(years.h0(), 6.893079992407984e-11);
    assert_eq!(years.omega_m(), 0.315);
    assert_eq!(years.omega_lambda(), 1.0 - 0.315);
    assert_eq!(years.time_unit_seconds(), 31_557_600.0);
    assert_eq!(years.length_unit_metres(), 9_460_730_472_580_800.0);
    assert_eq!(Some(years), Cosmology::from_km_s_mpc(67.4, 0.315));

    let seconds = Cosmology::PLANCK_2018_SI;
    assert_eq!(seconds.h0(), 2.184285241085502e-18);
    assert_eq!(seconds.time_unit_seconds(), 1.0);
    assert_eq!(seconds.length_unit_metres(), 299_792_458.0);

    assert!(GOLDEN.contains("mpmath,age,6.893079992407984e-11,0.315,"));
}

/// The two presets describe one universe: every quantity agrees once converted between the unit systems.
#[test]
fn the_unit_systems_agree() {
    let years = Cosmology::PLANCK_2018;
    let seconds = Cosmology::PLANCK_2018_SI;
    let year = years.time_unit_seconds();
    let close = |a: f64, b: f64| relative_error(a, b) < 1e-13;
    let now_y = years.age();
    let now_s = seconds.age();
    assert!(close(now_y * year, now_s));
    for m in [0.01, 0.5, 1.0, 3.0, 20.0] {
        let (ty, ts) = (now_y * m, now_s * m);
        assert!(close(years.scale_factor(ty), seconds.scale_factor(ts)));
        assert!(close(years.hubble(ty), seconds.hubble(ts) * year));
        assert!(close(
            years.conformal_time(ty) * year,
            seconds.conformal_time(ts)
        ));
        assert!(close(
            years.event_horizon(ty) * year,
            seconds.event_horizon(ts)
        ));
        assert!(close(years.reach_mass(ty), seconds.reach_mass(ts)));
        assert!(close(
            years.scale_factor_integral(ty) * year,
            seconds.scale_factor_integral(ts)
        ));
    }
    assert!(close(
        years.de_sitter_temperature(),
        seconds.de_sitter_temperature()
    ));
    assert!(close(years.h_lambda_si(), seconds.h_lambda()));
    assert!(close(years.h_lambda_si(), seconds.h_lambda_si()));
    assert!(close(
        years.lambda_tidal_acceleration_si(1e20),
        seconds.lambda_tidal_acceleration_si(1e20)
    ));
    assert!(close(
        years.acceleration_from_si(9.80665) / year,
        seconds.acceleration_from_si(9.80665)
    ));
    // converting a cosmology between units gives the other preset to rounding
    let converted = years.in_time_unit(1.0).unwrap();
    assert!(close(converted.h0(), seconds.h0()));
    assert_eq!(converted.time_unit_seconds(), 1.0);
}

/// Names that are documented as the same number are the same number.
#[test]
fn aliases_agree() {
    let u = Cosmology::PLANCK_2018;
    for m in [1e-3, 0.4, 1.0, 5.0, 60.0] {
        let t = u.age() * m;
        assert_eq!(u.particle_horizon(t), u.conformal_time(t));
        assert_eq!(u.event_horizon(t), u.conformal_time_remaining(t));
        assert_eq!(u.reach(t), u.conformal_time_remaining(t));
        let chi = 0.3 * u.event_horizon(t);
        assert_eq!(u.last_light_time(chi), u.last_contact_time(chi));
        assert_eq!(u.last_reply_time(chi), u.last_contact_time(2.0 * chi));
        assert_eq!(u.recession_speed(t, 7.0), u.hubble(t) * 7.0);
    }
    assert_eq!(u.de_sitter_horizon(), 1.0 / u.h_lambda());
    assert_eq!(u.de_sitter_recession_speed(3.0), u.h_lambda() * 3.0);
    assert_eq!(u.drag_limited_proper_velocity(2.0), 2.0 / u.h_lambda());
    assert_eq!(srkit_flrw::scale_factor_at_redshift(1.0), 0.5);
    assert_eq!(srkit_flrw::redshift_at_scale_factor(0.25), 3.0);
    assert_eq!(
        srkit_flrw::KM_S_MPC,
        1.0e3 / (1.0e6 * srkit_core::units::PARSEC)
    );
    assert_eq!(srkit_flrw::T_CMB0, 2.7255);
    assert_eq!(
        srkit_flrw::cmb_temperature_at_scale_factor(0.5, 2.7255),
        5.451
    );
    let today = u.cmb_temperature(u.age(), 2.7255);
    assert!((today - 2.7255).abs() < 1e-14);
}
