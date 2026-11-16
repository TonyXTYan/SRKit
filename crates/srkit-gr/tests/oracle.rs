//! The general-relativity and gravity rows of the shared physics oracle.
//!
//! `tools/physics-checks/oracle.csv` is the single source of the stated values; this test reads it and
//! computes each covered row with this crate's public functions, in the way `oracle_check.py` defines the row.
//!
//! A row passes if the computed value is within the row's own tolerance (column 6, ±1%). Nine of the covered
//! rows are quoted to two significant figures and the true value is more than 1% from the quoted one;
//! `oracle_check.py` classifies those as ROUNDED ("equal to the stated value at its quoted significant
//! figures, but outside ±1%"). For exactly those rows, listed in [`ROUNDED`], this test accepts a value that
//! rounds to the stated one at the stated number of figures, as the Python check does. No tolerance is
//! widened.

mod common;

use common::csv_fields;
use srkit_core::units::{AU, C, JULIAN_YEAR, M_SUN};
use srkit_gr::hawking::{self, EmissionModel};
use srkit_gr::kerr::{self, Orbit};
use srkit_gr::newtonian;
use std::collections::BTreeMap;

const ORACLE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../../tools/physics-checks/oracle.csv"
));

/// None since R71.10 (16 Nov 2026), when these nine rows were given a third figure. Before that:
/// rows that `oracle_check.py` reports as ROUNDED, or as AMBIGUOUS with a computed value outside ±1%: the
/// stated value has too few figures for its tolerance. Computed values, 16 Nov 2026:
/// SR-40 42.52 (stated 43), SR-41 2.732 (2.7), SR-42 30.62 (31), SR-43 2.135 (2.1), SR-48 2.242e32 (2.2e32),
/// SR-49 2.064e36 (2.1e36), SR-63 2.031e25 (2.0e25), SR-64 3.562e14 (3.6e14), SR-65 2665 (2700).
const ROUNDED: &[&str] = &[];

/// The rows this crate is responsible for.
const COVERED: &[&str] = &[
    "SR-40", "SR-41", "SR-42", "SR-43", "SR-47", "SR-48", "SR-49", "SR-53", "SR-54", "SR-55",
    "SR-56", "SR-57", "SR-58", "SR-59", "SR-60", "SR-61", "SR-62", "SR-63", "SR-64", "SR-65",
    "DT-19", "DT-20", "DT-21", "DT-22",
];

struct Row {
    quantity: String,
    stated: String,
    unit: String,
    tolerance: f64,
}

fn oracle() -> BTreeMap<String, Row> {
    let mut rows = BTreeMap::new();
    let mut lines = ORACLE.lines().filter(|line| !line.is_empty());
    assert_eq!(
        lines.next(),
        Some("id,group,quantity,stated,unit,tolerance"),
        "unexpected oracle header"
    );
    for record in lines {
        let fields = csv_fields(record);
        assert_eq!(fields.len(), 6, "bad oracle record: {record}");
        let row = Row {
            quantity: fields[2].clone(),
            stated: fields[3].clone(),
            unit: fields[4].clone(),
            tolerance: fields[5]
                .parse()
                .unwrap_or_else(|_| panic!("bad tolerance: {record}")),
        };
        assert!(
            rows.insert(fields[0].clone(), row).is_none(),
            "repeated oracle id {}",
            fields[0]
        );
    }
    rows
}

/// Significant figures of a stated value, as `oracle_check.py` counts them: trailing zeros of an integer
/// are not significant ("2700" has two figures, "2.10e67" three).
fn significant_figures(stated: &str) -> usize {
    let lower = stated.to_lowercase();
    let mantissa = lower
        .trim_start_matches(['+', '-'])
        .split('e')
        .next()
        .unwrap();
    let mut digits = mantissa.replace('.', "");
    digits = digits.trim_start_matches('0').to_string();
    if !mantissa.contains('.') {
        digits = digits.trim_end_matches('0').to_string();
    }
    digits.len().max(1)
}

/// `x` rounded to `figures` significant figures.
fn round_to_figures(x: f64, figures: usize) -> f64 {
    format!("{:.*e}", figures - 1, x).parse().unwrap()
}

#[test]
fn csv_fields_handle_quotes_and_commas() {
    assert_eq!(csv_fields("a,b,,c"), ["a", "b", "", "c"]);
    assert_eq!(
        csv_fields(r#"SR-53,spine,"Kerr ISCO dt/dtau, a=0",1.41,,0.01"#),
        [
            "SR-53",
            "spine",
            "Kerr ISCO dt/dtau, a=0",
            "1.41",
            "",
            "0.01"
        ]
    );
    assert_eq!(
        csv_fields(r#"x,"say ""hi"", twice",y"#),
        ["x", r#"say "hi", twice"#, "y"]
    );
    assert_eq!(csv_fields(""), [""]);
}

#[test]
fn significant_figures_follow_the_python_check() {
    for (stated, figures) in [
        ("43", 2),
        ("2.7", 2),
        ("2700", 2),
        ("2.10e67", 3),
        ("6.8e4", 2),
        ("145", 3),
    ] {
        assert_eq!(significant_figures(stated), figures, "{stated}");
    }
    assert_eq!(round_to_figures(2665.4, 2), 2700.0);
    assert_eq!(round_to_figures(3.5616e14, 2), 3.6e14);
    assert_eq!(round_to_figures(2.7316, 2), 2.7);
}

#[test]
fn gravity_rows_of_the_physics_oracle() {
    let reference = oracle();
    let mut computed: BTreeMap<&str, f64> = BTreeMap::new();
    let mut add = |id: &'static str, value: f64| {
        assert!(value.is_finite(), "{id}: non-finite result {value}");
        assert!(computed.insert(id, value).is_none(), "{id} computed twice");
    };

    // SR-40 to SR-43: Hohmann transfers from the Earth's orbit and the Mars synodic period, shown at the
    // game's rate of one year per minute. "real" seconds are 60 times the years; "real" minutes are the years.
    let years = |seconds: f64| seconds / JULIAN_YEAR;
    add(
        "SR-40",
        60.0 * years(newtonian::hohmann_transfer_time(M_SUN, AU, 1.523679 * AU)),
    );
    add(
        "SR-41",
        years(newtonian::hohmann_transfer_time(M_SUN, AU, 5.2038 * AU)),
    );
    add(
        "SR-42",
        years(newtonian::hohmann_transfer_time(M_SUN, AU, 30.07 * AU)),
    );
    // sidereal periods of the Earth and Mars in days; the result in Julian years
    add(
        "SR-43",
        newtonian::synodic_period(365.256363, 686.980) / 365.25,
    );

    // SR-47 to SR-49: uniform-sphere binding energy of Mercury, the Earth and Jupiter (mass in kg, mean
    // radius in m, as in oracle_check.py).
    add(
        "SR-47",
        newtonian::binding_energy_uniform_sphere(3.3011e23, 2.4397e6),
    );
    add(
        "SR-48",
        newtonian::binding_energy_uniform_sphere(5.9722e24, 6.3710e6),
    );
    add(
        "SR-49",
        newtonian::binding_energy_uniform_sphere(1.89813e27, 6.9911e7),
    );

    // SR-53 to SR-62: dt/dtau on the prograde ISCO and the thin-disc efficiency, in percent.
    for (spin, dt_id, efficiency_id) in [
        (0.0, "SR-53", "SR-54"),
        (0.9, "SR-55", "SR-56"),
        (0.998, "SR-57", "SR-58"),
        (1.0 - 1e-6, "SR-59", "SR-60"),
        (1.0 - 1e-14, "SR-61", "SR-62"),
    ] {
        add(dt_id, kerr::isco_dt_dtau(spin, Orbit::Prograde));
        add(
            efficiency_id,
            100.0 * kerr::accretion_efficiency(spin, Orbit::Prograde),
        );
    }

    // SR-63: the kinetic energy (gamma - 1) m c^2 of 1e6 kg at gamma = 227. This is special relativity, not
    // gravity: the "spin-up" is of a loop, not of a hole. No function of this crate computes it; it is here
    // so that the block SR-53 to SR-65 has no gap.
    add("SR-63", (227.0 - 1.0) * 1e6 * C * C);

    // SR-64, SR-65, DT-19 to DT-22: Hawking radiation. The power and lifetime rows are the photon-only black
    // body of the horizon area (`oracle_check.py` marks them AMBIGUOUS because the row does not say so).
    let photons = EmissionModel::PhotonBlackbody;
    add("SR-64", hawking::power(1e9, photons));
    add("SR-65", years(hawking::lifetime(1e9, photons)));
    add("DT-19", hawking::temperature(1e11 * M_SUN));
    add("DT-20", hawking::temperature(M_SUN));
    add("DT-21", years(hawking::lifetime(1e11 * M_SUN, photons)));
    add("DT-22", years(hawking::lifetime(M_SUN, photons)));

    assert_eq!(computed.keys().copied().collect::<Vec<_>>(), {
        let mut ids = COVERED.to_vec();
        ids.sort_unstable();
        ids
    });

    let mut failures = Vec::new();
    println!(
        "{:<6} {:<34} {:>10} {:>14} {:>8}  status",
        "id", "quantity", "stated", "computed", "diff"
    );
    for (id, &value) in &computed {
        let row = reference
            .get(*id)
            .unwrap_or_else(|| panic!("{id} is not in oracle.csv"));
        let stated: f64 = row.stated.replace(',', "").parse().unwrap();
        let difference = (value - stated) / stated;
        let within = difference.abs() <= row.tolerance;
        let rounds =
            (round_to_figures(value, significant_figures(&row.stated)) / stated - 1.0).abs() < 1e-9;
        let status = if within {
            "MATCH"
        } else if rounds && ROUNDED.contains(id) {
            "ROUNDED"
        } else {
            failures.push(format!(
                "{id} ({}): computed {value:e} {}, stated {}, difference {:+.2}% against a tolerance of {}%",
                row.quantity,
                row.unit,
                row.stated,
                100.0 * difference,
                100.0 * row.tolerance
            ));
            "MISMATCH"
        };
        let quantity: String = row.quantity.chars().take(34).collect();
        println!(
            "{id:<6} {quantity:<34} {:>10} {value:>14.5e} {:>+7.2}%  {status}",
            row.stated,
            100.0 * difference
        );
    }
    assert!(
        failures.is_empty(),
        "oracle rows not reproduced:\n{}",
        failures.join("\n")
    );
}

/// `tools/physics-checks/hawking_scaling.py`: the photon-only lifetime of one solar mass, its cubic scaling
/// with mass, and the hole whose lifetime overflows a 384-bit count of milliseconds.
#[test]
fn hawking_scaling_script() {
    let photons = EmissionModel::PhotonBlackbody;
    let solar_years = hawking::lifetime(M_SUN, photons) / JULIAN_YEAR;
    assert!(
        (solar_years / 2.10e67 - 1.0).abs() < 0.005,
        "t_ev(1 Msun) = {solar_years:e} yr"
    );

    // t is proportional to M^3 (and to G^2, which this crate does not vary).
    for factor in [10.0, 1e-3, 7.5] {
        let ratio = hawking::lifetime(factor * M_SUN, photons) / hawking::lifetime(M_SUN, photons);
        assert!(
            (ratio / (factor * factor * factor) - 1.0).abs() < 1e-14,
            "scaling at {factor}"
        );
    }

    // U384 milliseconds: 2^384 ms = 1.25e105 yr; the script prints "overflow above 3.91e+12 Msun".
    let u384_seconds = 2f64.powi(384) / 1e3;
    assert!((u384_seconds / JULIAN_YEAR / 1.25e105 - 1.0).abs() < 0.005);
    let overflow_mass = hawking::mass_for_lifetime(u384_seconds, photons) / M_SUN;
    assert!(
        (overflow_mass / 3.91e12 - 1.0).abs() < 0.005,
        "overflow above {overflow_mass:e} Msun"
    );
}

/// The other emission models quoted in the oracle's notes for the ambiguous rows.
#[test]
fn hawking_models_in_the_oracle_notes() {
    // DT-22 note: "Photons+gravitons with greybody factors (Page 1976): 1.16e+67 yr".
    let page = hawking::lifetime(M_SUN, EmissionModel::PagePhotonsGravitons) / JULIAN_YEAR;
    assert!(
        (page / 1.16e67 - 1.0).abs() < 0.005,
        "Page lifetime {page:e} yr"
    );

    // SR-64 and SR-65 notes: MacGibbon with f ~ 12 gives "~6e+16 W" and "~16 yr" for 1e9 kg.
    let standard_model = EmissionModel::MacGibbon { f: 12.0 };
    let power = hawking::power(1e9, standard_model);
    let lifetime = hawking::lifetime(1e9, standard_model) / JULIAN_YEAR;
    assert!(
        (power / 5.76e16 - 1.0).abs() < 0.005,
        "MacGibbon power {power:e} W"
    );
    assert!(
        (lifetime / 16.5 - 1.0).abs() < 0.005,
        "MacGibbon lifetime {lifetime} yr"
    );

    // The model changes the 1e9 kg numbers by a factor of about 160.
    let ratio = power / hawking::power(1e9, EmissionModel::PhotonBlackbody);
    assert!((150.0..170.0).contains(&ratio), "ratio {ratio}");
}
