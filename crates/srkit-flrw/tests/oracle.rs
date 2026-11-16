//! The cosmology rows of the shared physics oracle, `oracle/oracle.csv`.
//!
//! Each row is computed twice: with the oracle's own convention (`1/H_Λ` = 17.53 Gyr, `Ω_m` = 0.315, as in
//! `oracle_check.py`) and with the crate's Planck 2018 preset (`1/H_Λ` = 17.528 Gyr). Both must agree with the
//! stated value under the oracle's rule: within the row's relative tolerance (MATCH), or equal to the stated
//! value once rounded to the significant figures it is quoted with (ROUNDED), for the rows `oracle_check.py`
//! itself classes as ROUNDED.

mod common;

use common::csv_fields;
use srkit_core::units::{BOLTZMANN, C, G0, LIGHT_YEAR};
use srkit_flrw::{Cosmology, T_CMB0};
use std::collections::BTreeMap;

const ORACLE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../oracle/oracle.csv"
));

/// Solar mass in kilograms, as `oracle_check.py` has it since R71.10: the IAU nominal GM☉ over CODATA G.
const M_SUN: f64 = srkit_core::units::M_SUN;

/// None since R71.10 (16 Nov 2026), when SR-68, DT-18 and DT-29 were given a third figure. Before that:
/// rows whose stated value is quoted with too few digits to meet its own 1% tolerance. `oracle_check.py`
/// reports the same three as ROUNDED: 1.2287e30 against 1.2e30, 1.2273e12 yr against 1.2e12 yr, and
/// 1.2455e13 yr against 1.2e13 yr.
const ROUNDED: [&str; 0] = [];

/// Number of significant figures a stated value is quoted with (the rule in `oracle_check.py`).
fn significant_figures(stated: &str) -> usize {
    let mantissa = stated
        .to_ascii_lowercase()
        .trim_start_matches(['+', '-'])
        .split('e')
        .next()
        .unwrap()
        .to_string();
    let mut digits = mantissa
        .replace('.', "")
        .trim_start_matches('0')
        .to_string();
    if !mantissa.contains('.') {
        digits = digits.trim_end_matches('0').to_string();
    }
    digits.len().max(1)
}

/// `x` rounded to `n` significant figures.
fn round_to(x: f64, n: usize) -> f64 {
    format!("{:.*e}", n - 1, x).parse().unwrap()
}

/// The cosmology rows, computed with one cosmology in Julian years and light-years.
fn compute(u: &Cosmology) -> BTreeMap<&'static str, f64> {
    let g = u.acceleration_from_si(G0);
    let t_ds = u.de_sitter_temperature();
    let crossover = u.cmb_crossover_scale_factor(T_CMB0);
    let overflow = u.scale_factor_overflow_time();
    BTreeMap::from([
        // Hubble drag on a 1 g rocket: gamma saturates at sqrt(1 + (g / (c H_Λ))²)
        ("SR-20", u.drag_limited_gamma(g)),
        // ... and its proper time over 1e100 yr of cosmic time, from rest, in de Sitter space
        ("SR-21", u.de_sitter_thrust_proper_time(0.0, g, 1e100)),
        // Landauer gain from a 2.7 K sky to the de Sitter temperature
        ("SR-68", u.cmb_crossover_scale_factor(2.7)),
        ("DT-15", t_ds),
        ("DT-16", crossover),
        ("DT-17", crossover.ln()),
        ("DT-18", u.cmb_crossover_time(T_CMB0)),
        // an island's lifetime budget: 1e12 solar masses of rest energy, erased at the de Sitter temperature
        (
            "DT-23",
            1e12 * M_SUN * C * C / (BOLTZMANN * t_ds * core::f64::consts::LN_2),
        ),
        ("DT-25", u.lambda_tidal_acceleration_si(5.0e4 * LIGHT_YEAR)),
        ("DT-26", u.de_sitter_recession_speed(1.0e5)),
        ("DT-27", u.last_light_gain(10.0) / 1e9),
        ("DT-28", u.h_lambda() * overflow),
        ("DT-29", overflow),
    ])
}

#[test]
fn cosmology_rows_of_the_oracle() {
    let mut table = BTreeMap::new();
    for record in ORACLE.lines().skip(1).filter(|line| !line.is_empty()) {
        let fields = csv_fields(record);
        assert_eq!(fields.len(), 6, "malformed oracle record: {record}");
        assert!(
            table.insert(fields[0].clone(), fields).is_none(),
            "repeated oracle id in: {record}"
        );
    }

    let oracle_convention = Cosmology::from_h_lambda_per_year(1.0 / 17.53e9, 0.315).unwrap();
    let cases = [
        ("1/H_Λ = 17.53 Gyr", oracle_convention),
        ("Planck 2018 preset", Cosmology::PLANCK_2018),
    ];
    let mut checked = 0;
    for (label, u) in cases {
        println!("{label}");
        for (id, value) in compute(&u) {
            let row = table
                .get(id)
                .unwrap_or_else(|| panic!("oracle.csv has no row {id}"));
            let stated: f64 = row[3].parse().unwrap();
            let tolerance: f64 = row[5].parse().unwrap();
            let relative = (value - stated) / stated;
            let status = if relative.abs() <= tolerance {
                "MATCH"
            } else if round_to(value, significant_figures(&row[3])) == stated {
                "ROUNDED"
            } else {
                "MISMATCH"
            };
            println!(
                "  {status:<8} {id}  {:<42} stated {:>9}  computed {value:<12.5e} {:+.2}%",
                row[2],
                row[3],
                100.0 * relative
            );
            assert!(value.is_finite(), "{id}: not finite");
            assert!(
                status == "MATCH" || (status == "ROUNDED" && ROUNDED.contains(&id)),
                "{id} ({label}): computed {value:e}, stated {stated:e}, tolerance {tolerance}: {status}"
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 2 * 13);
}

/// The oracle's own computed values (`oracle_check.py`, 1/H_Λ = 17.53 Gyr), to the digits it prints.
#[test]
fn agrees_with_oracle_check_py_to_its_printed_digits() {
    let u = Cosmology::from_h_lambda_per_year(1.0 / 17.53e9, 0.315).unwrap();
    let computed = compute(&u);
    let printed = [
        ("SR-20", 1.8096e10),
        ("SR-21", 5.526e89),
        ("SR-68", 1.2287e30),
        ("DT-15", 2.1975e-30),
        ("DT-16", 1.2403e30),
        ("DT-17", 69.293),
        ("DT-18", 1.2273e12),
        ("DT-23", 8.4979e111), // 8.4982e111 before R71.10 moved the solar mass to GM☉ / G
        ("DT-25", 1.5457e-15),
        ("DT-26", 5.7045e-6),
        ("DT-27", 13.455),
        ("DT-28", 710.5),
        ("DT-29", 1.2455e13),
    ];
    for (id, value) in printed {
        let digits = if id == "DT-28" || id == "SR-21" { 4 } else { 5 };
        assert_eq!(
            round_to(computed[id], digits),
            value,
            "{id}: {:e}",
            computed[id]
        );
    }
    // SR-21 is T / gamma_sat plus a transient the oracle integrates numerically: 23.55 yr
    let g = u.acceleration_from_si(G0);
    let gamma = u.drag_limited_gamma(g);
    let transient = u.de_sitter_thrust_proper_time(0.0, g, 1.0e13) - 1.0e13 / gamma;
    assert!((transient - 23.55).abs() < 0.005, "transient {transient}");
}

/// DT-28 and DT-29 say where `a` overflows binary64; the scale factor does overflow there and nowhere before.
#[test]
fn the_scale_factor_overflows_where_the_oracle_says() {
    for u in [Cosmology::PLANCK_2018, Cosmology::PLANCK_2018_SI] {
        let t = u.scale_factor_overflow_time();
        assert!(u.scale_factor(t * (1.0 - 1e-9)).is_finite());
        assert!(u.scale_factor(t * (1.0 + 1e-9)).is_infinite());
        assert!((u.ln_scale_factor(t) - f64::MAX.ln()).abs() < 1e-9);
        // the logarithm carries on
        assert!(u.ln_scale_factor(1e6 * t).is_finite());
    }
}
