//! The implemented SR subset of the shared physics oracle. GR, FLRW, gameplay
//! time scaling and wide-clock rows are deliberately not claimed by this test.
use srkit_core::kinematics::{doppler_factor, lorentz_factor, Event, Vec3, C, G0, JULIAN_YEAR};
use srkit_core::segment::Segment;
use srkit_core::trip::{circular_gamma, FlipBurn};
use std::collections::{BTreeMap, BTreeSet};

const ORACLE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../../tools/physics-checks/oracle.csv"
));
const LY: f64 = C * JULIAN_YEAR;

// The source has quoted quantities containing commas. Parse CSV fields rather
// than splitting columns on commas; doubled quotes inside quoted fields are
// accepted too. The repository's records do not contain embedded newlines.
fn fields(record: &str) -> Vec<String> {
    let mut values = Vec::new();
    let mut value = String::new();
    let mut quoted = false;
    let mut chars = record.chars().peekable();
    while let Some(ch) = chars.next() {
        match ch {
            '"' if quoted && chars.peek() == Some(&'"') => {
                value.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            ',' if !quoted => {
                values.push(std::mem::take(&mut value));
            }
            _ => value.push(ch),
        }
    }
    assert!(!quoted, "unterminated CSV field");
    values.push(value);
    values
}

#[test]
fn implemented_sr_oracle_subset() {
    let mut reference = BTreeMap::new();
    for record in ORACLE.lines().skip(1).filter(|line| !line.is_empty()) {
        let columns = fields(record);
        assert_eq!(columns.len(), 6, "invalid oracle record: {record}");
        assert!(reference.insert(columns[0].clone(), columns).is_none());
    }
    let mut computed = BTreeMap::new();
    let mut add = |id: &str, value: f64| {
        assert!(value.is_finite(), "nonfinite result for {id}");
        assert!(computed.insert(id.to_string(), value).is_none());
    };

    add("SR-01", C / 10.0 / JULIAN_YEAR);
    add("SR-02", C * C / 10.0 / LY);
    let trip = FlipBurn::new(1e5 * LY, G0).unwrap();
    add("SR-03", trip.proper_time / JULIAN_YEAR);
    add("SR-04", trip.coordinate_time / JULIAN_YEAR);
    add("SR-05", trip.peak_gamma);
    add("SR-06", trip.photon_rocket_mass_ratio().unwrap());
    add("SR-07", 2.0 * trip.peak_rapidity);

    let origin = Event::new(0.0, Vec3::ZERO).unwrap();
    let burn = Segment::accelerated(
        origin,
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 1.0),
        G0,
        2.0 * JULIAN_YEAR,
    )
    .unwrap();
    let time_to_point9c = C / G0 * 0.9 / libm::sqrt(1.0 - 0.9 * 0.9);
    add(
        "SR-08",
        burn.proper_time(time_to_point9c).unwrap() / JULIAN_YEAR,
    );
    for (id, radius) in [
        ("SR-09", 1.0),
        ("SR-10", 100.0),
        ("SR-11", 1e4),
        ("SR-12", 5e4),
    ] {
        add(id, circular_gamma(radius * LY, G0).unwrap());
    }

    let radius = 100.0 * LY;
    let gamma = circular_gamma(radius, G0).unwrap();
    let beta = libm::sqrt(1.0 - 1.0 / (gamma * gamma));
    let period = 2.0 * core::f64::consts::PI * radius / (beta * C);
    let orbit = Segment::circular(
        Event::new(0.0, Vec3::new(radius, 0.0, 0.0)).unwrap(),
        Vec3::ZERO,
        Vec3::new(radius, 0.0, 0.0),
        Vec3::new(0.0, radius, 0.0),
        beta * C / radius,
        period,
    )
    .unwrap();
    add("SR-13", orbit.duration() / JULIAN_YEAR);
    add("SR-14", orbit.proper_duration() / JULIAN_YEAR);
    let receiver_velocity = Vec3::new(beta * C, 0.0, 0.0);
    add(
        "SR-15",
        doppler_factor(Vec3::ZERO, receiver_velocity, Vec3::new(1.0, 0.0, 0.0)).unwrap(),
    );
    add(
        "SR-16",
        doppler_factor(Vec3::ZERO, receiver_velocity, Vec3::new(-1.0, 0.0, 0.0)).unwrap(),
    );

    // Uniform coordinate-time samples are also uniform in proper time at
    // constant gamma. Each source is stationary and outside the circular path.
    // Average received/emitted frequency over a closed lap is exactly gamma.
    let mut worst_mean_error: f64 = 0.0;
    for (r, height) in [(1e3, 0.0), (300.0, 50.0), (1e6, 0.0), (150.0, -20.0)] {
        let source = Vec3::new(r * LY, 0.0, height * LY);
        let mut total = 0.0;
        for i in 0..4096 {
            let state = orbit.state(period * (i as f64 / 4096.0)).unwrap();
            total +=
                doppler_factor(Vec3::ZERO, state.velocity, state.event.position - source).unwrap();
        }
        worst_mean_error = worst_mean_error.max((total / 4096.0 / gamma - 1.0).abs());
    }
    assert!(
        worst_mean_error < 1e-12,
        "closed-lap Doppler error {worst_mean_error}"
    );
    add("SR-17", 1.0 + worst_mean_error);
    let nonrelativistic_gap = 2.0 * core::f64::consts::PI * radius * beta / C;
    add("SR-18", nonrelativistic_gap * gamma / JULIAN_YEAR);
    add("SR-19", nonrelativistic_gap * gamma * gamma / JULIAN_YEAR);

    add(
        "DT-11",
        doppler_factor(
            Vec3::ZERO,
            Vec3::new(C * libm::tanh(5.0), 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .unwrap(),
    );

    for (beta, energy_id, gamma_id, blue_id) in [
        (0.9999, "SR-27", "SR-28", "SR-29"),
        (0.99999, "SR-30", "SR-31", "SR-32"),
    ] {
        let velocity = Vec3::new(0.0, beta * C, 0.0);
        let gamma = lorentz_factor(velocity).unwrap();
        add(energy_id, (gamma - 1.0) * C * C / 4.184e18);
        add(gamma_id, gamma);
        add(
            blue_id,
            doppler_factor(Vec3::ZERO, velocity, Vec3::new(0.0, -1.0, 0.0)).unwrap(),
        );
    }
    for (id, distance) in [("SR-44", 4.2465), ("SR-45", 1e5)] {
        let ratio = FlipBurn::new(distance * LY, G0)
            .unwrap()
            .photon_rocket_mass_ratio()
            .unwrap();
        add(id, (ratio - 1.0) * 1e6 / 2.0);
    }
    for (id, beta) in [("SR-50", 0.01), ("SR-51", 0.1)] {
        add(
            id,
            (lorentz_factor(Vec3::new(0.0, 0.0, beta * C)).unwrap() - 1.0) * C * C / 4.184e12,
        );
    }
    add("SR-63", (227.0 - 1.0) * 1e6 * C * C);
    add("SR-72", 220_000.0 / C);
    add(
        "SR-73",
        lorentz_factor(Vec3::new(0.0, 220_000.0, 0.0)).unwrap() - 1.0,
    );

    // These two rows have only two stated significant figures and are already
    // classified ROUNDED by oracle_check.py: 1.0741 -> 1.1 kt, and
    // 2.0312e25 -> 2.0e25 J. Accept their rounding cells, not a wider blanket
    // relative tolerance. All other covered rows must meet the CSV tolerance.
    let rounding_cells = BTreeMap::from([("SR-50", (1.1, 0.05)), ("SR-63", (2.0e25, 0.05e25))]);
    for (id, value) in &computed {
        let columns = reference
            .get(id)
            .unwrap_or_else(|| panic!("missing row {id}"));
        let stated: f64 = columns[3].parse().unwrap();
        let tolerance: f64 = columns[5].parse().unwrap();
        let error = (value - stated).abs();
        if let Some(&(expected_stated, half_cell)) = rounding_cells.get(id.as_str()) {
            assert_eq!(
                stated, expected_stated,
                "review rounding exemption for {id}"
            );
            assert!(
                error <= half_cell,
                "{id}: {value} does not round to {stated}"
            );
        } else {
            assert!(
                error <= stated.abs() * tolerance,
                "{id}: {value} versus {stated}, tolerance {tolerance}"
            );
        }
    }
    let covered: BTreeSet<_> = computed.keys().map(String::as_str).collect();
    let expected: BTreeSet<_> = [
        "SR-01", "SR-02", "SR-03", "SR-04", "SR-05", "SR-06", "SR-07", "SR-08", "SR-09", "SR-10",
        "SR-11", "SR-12", "SR-13", "SR-14", "SR-15", "SR-16", "SR-17", "SR-18", "SR-19", "SR-27",
        "SR-28", "SR-29", "SR-30", "SR-31", "SR-32", "SR-44", "SR-45", "SR-50", "SR-51", "SR-63",
        "SR-72", "SR-73", "DT-11",
    ]
    .into_iter()
    .collect();
    assert_eq!(covered, expected);
    assert_eq!(reference.len(), 111, "review coverage after oracle changes");
    assert_eq!(covered.len(), 33);
    assert_eq!(
        reference.len() - covered.len(),
        78,
        "unimplemented oracle rows"
    );
    println!("Shared oracle: 33 covered rows; 78 explicitly outside this SR milestone");
}
