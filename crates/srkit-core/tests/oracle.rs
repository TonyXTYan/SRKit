//! The physics oracle (`tools/physics-checks/oracle.csv`): every row that is special relativity, optics or
//! exact timekeeping, recomputed through this crate's public API. The cosmology rows are checked by `srkit-flrw`
//! and the gravity rows by `srkit-gr`; rows that are game design (real-time pacing, engineering estimates) are not
//! this library's to claim and are left to `oracle_check.py`.
mod common;

use srkit_core::{
    causal::{self, GridVec},
    kinematics::{self, Event},
    optics,
    orbit::CircularOrbit,
    rocket,
    segment::Segment,
    time::{Time384, Time512},
    trip::{FlipBurn, Trip},
    units::{self, Units, AU, C, G0, G0_PER_YEAR, JULIAN_YEAR, LIGHT_YEAR},
    Vec3,
};

const J_PER_KT: f64 = 4.184e12;
const D_PROXIMA: f64 = 4.2465; // ly

#[test]
fn special_relativity_rows() {
    let yr = Units::YEARS;
    let g = G0_PER_YEAR;
    let mut rows: Vec<(&str, f64)> = Vec::new();

    // the scales of a 10 m/s² drive
    rows.push(("SR-01", 1.0 / yr.accel_from_si(10.0)));
    rows.push(("SR-02", yr.length_from_si(C * C / 10.0)));

    // 1 g flip-and-burn across the galaxy
    let trip = FlipBurn::new(1e5, g).unwrap();
    rows.push(("SR-03", trip.tau));
    rows.push(("SR-04", trip.t));
    rows.push(("SR-05", trip.gamma_peak));
    rows.push(("SR-06", trip.mass_ratio(1.0)));
    rows.push(("SR-07", 2.0 * trip.phi_peak));
    // the same ride as a plan of legs must agree with the closed form
    let plan = Trip::new(1e5, g, f64::INFINITY, 0, false).unwrap();
    assert!((plan.tau / trip.tau - 1.0).abs() < 1e-14 && (plan.lag / trip.lag - 1.0).abs() < 1e-12);

    // proper time to reach 0.9 c at 1 g: evaluate the burn where its speed is 0.9 c
    let burn = Segment::hyperbolic(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, g, 0.0);
    let t_09 = 0.9 / (1.0 - 0.81_f64).sqrt() / g; // gamma beta / a
    assert!((burn.at(t_09).v.x - 0.9).abs() < 1e-15);
    rows.push(("SR-08", burn.dtau(t_09)));

    // loops flown at 1 g
    for (id, r) in [
        ("SR-09", 1.0),
        ("SR-10", 100.0),
        ("SR-11", 1e4),
        ("SR-12", 5e4),
    ] {
        rows.push((id, CircularOrbit::from_radius_acceleration(r, g).gamma));
    }
    let ring = CircularOrbit::from_radius_acceleration(100.0, g);
    rows.push(("SR-13", ring.period));
    rows.push(("SR-14", ring.lap_tau));
    let (d_min, d_max) = ring.doppler_range();
    rows.push(("SR-15", d_min));
    rows.push(("SR-16", d_max));
    // lap average of the Doppler factor of a distant source in the plane, over gamma
    let n = 3600;
    let mean: f64 = (0..n)
        .map(|k| {
            let th = 2.0 * std::f64::consts::PI * (k as f64 + 0.5) / n as f64;
            optics::observe(ring.phi, th.cos(), th.sin(), 0.0, true)
                .unwrap()
                .doppler
        })
        .sum::<f64>()
        / n as f64;
    rows.push(("SR-17", mean / ring.gamma));
    rows.push(("SR-18", ring.sagnac()));
    rows.push(("SR-19", ring.sagnac() * ring.gamma));

    // 1 g held for 1e100 years of coordinate time, in flat space
    rows.push(("SR-22", burn.dtau(1e100)));
    rows.push(("DT-12", burn.dtau(1e100) * g));

    // warning time: light from a launch 20 light-minutes away against the projectile itself
    let minutes = Units::new(60.0);
    for (id, beta, unit_s) in [
        ("SR-23", 0.01, 3600.0),
        ("SR-24", 0.5, 60.0),
        ("SR-25", 0.9999, 1.0),
        ("SR-26", 0.99999, 1e-3),
    ] {
        let d = 20.0;
        let shot = Segment::inertial(
            0.0,
            f64::INFINITY,
            Vec3::new(d, 0.0, 0.0),
            Vec3::new(-beta, 0.0, 0.0),
        );
        let target = Segment::rest(f64::NEG_INFINITY, f64::INFINITY, Vec3::ZERO);
        let seen = target
            .recv(
                Event::new(0.0, Vec3::new(d, 0.0, 0.0)),
                f64::NEG_INFINITY,
                f64::INFINITY,
            )
            .unwrap();
        // when the shot reaches the target: its light, sent from the target's position, "arrives" at once
        let hit = shot
            .send(Event::new(d / beta, Vec3::ZERO), 0.0, f64::INFINITY)
            .unwrap();
        assert!((hit - d / beta).abs() <= 1e-9 * d / beta);
        // the shot covers the last (1 - beta) of the light's lead: d / beta - d = d (1 - beta) / beta
        let warning = d * kinematics::speed_shortfall(kinematics::rapidity_of_speed(beta)) / beta;
        assert!((warning - (hit - seen)).abs() <= 1e-9 * hit);
        rows.push((id, minutes.time_to_si(warning) / unit_s));
    }

    // kinetic energy, gamma and head-on blueshift near c
    for (beta, ke, gam, blue) in [
        (0.9999, "SR-27", "SR-28", "SR-29"),
        (0.99999, "SR-30", "SR-31", "SR-32"),
    ] {
        let phi = kinematics::rapidity_of_speed(beta);
        rows.push((
            ke,
            rocket::kinetic_per_rest_energy(phi) * C * C / (J_PER_KT * 1e6),
        ));
        rows.push((gam, kinematics::gamma_of_speed(beta).unwrap()));
        rows.push((blue, kinematics::doppler_of_rapidity(phi)));
    }
    for (id, beta) in [("SR-50", 0.01), ("SR-51", 0.1)] {
        rows.push((
            id,
            rocket::kinetic_per_rest_energy(kinematics::rapidity_of_speed(beta)) * C * C / J_PER_KT,
        ));
    }
    // a slow ship overtaken by a later, faster one: the gap between their arrivals at Proxima
    let (slow, fast) = (
        Segment::inertial(0.0, f64::INFINITY, Vec3::ZERO, Vec3::new(0.01, 0.0, 0.0)),
        Segment::inertial(0.0, f64::INFINITY, Vec3::ZERO, Vec3::new(0.1, 0.0, 0.0)),
    );
    let at = |seg: &Segment| D_PROXIMA / seg.at(0.0).v.x;
    assert!((slow.at(at(&slow)).x.x - D_PROXIMA).abs() < 1e-12);
    rows.push(("SR-52", at(&slow) - at(&fast)));
    rows.push((
        "SR-63",
        (CircularOrbit::from_radius_acceleration(5e4, g)
            .gamma
            .round()
            - 1.0)
            * 1e6
            * C
            * C,
    ));

    // beams and apertures
    rows.push((
        "SR-69",
        optics::rayleigh_resolution(1e6, 1e-6, 100.0 * LIGHT_YEAR) / 1.22 / 1e3,
    ));
    rows.push(("SR-70", 1e6 * 1e6 / 1e-6 / LIGHT_YEAR));
    rows.push(("SR-87", optics::rayleigh_resolution(AU, 1e-6, 1.0) / 1.22));

    // the galaxy's rotation
    let beta_rot = Units::SECONDS.speed_from_si(220e3);
    rows.push(("SR-72", beta_rot));
    rows.push((
        "SR-73",
        rocket::kinetic_per_rest_energy(kinematics::rapidity_of_speed(beta_rot)),
    ));

    // a 1,000 t ship to Proxima at 1 g: a beamed sail to the turnover, and a photon rocket for the whole ride
    let prox = FlipBurn::new(D_PROXIMA, g).unwrap();
    let rest_energy = 1e6 * C * C;
    // a sail held at 1 g receives its beam redshifted twice over (D², D = e^-phi): E = (m c²/2)(e^phi - 1)
    let sail_energy = 0.5 * rest_energy * libm_expm1(prox.phi_peak);
    assert!((rocket::sail_rapidity(sail_energy / rest_energy) - prox.phi_peak).abs() < 1e-12);
    rows.push(("SR-80", sail_energy));
    rows.push(("SR-81", 0.5 * 1e6 * G0 * C * (2.0 * prox.phi_peak).exp()));
    rows.push((
        "SR-82",
        rocket::propellant_per_dry_mass(2.0 * prox.phi_peak, 1.0) * rest_energy,
    ));

    // playback rate of the launch site's clock at rapidity 5 on the way out
    let far = Trip::new(1e6, g, f64::INFINITY, 0, false).unwrap();
    rows.push(("DT-11", far.state(5.0 / g).doppler_home));

    // exact clocks: ranges of the stored integer, in Julian years
    let ms_per_year = JULIAN_YEAR * 1e3;
    let as_f64 = |text: String| text.parse::<f64>().unwrap();
    rows.push(("DT-02", as_f64(Time384::MAX.to_string()) / ms_per_year));
    assert!(as_f64(Time512::MAX.to_string()) / ms_per_year > 1e140);
    // an orbit's phase at 1e20 years by integer modular reduction is exact: the same phase as at the remainder
    let period_ms = 31_557_600_000_u64;
    let late = Time384::from_decimal(&format!("3155760000000{}", "0".repeat(17))).unwrap();
    let phase = late
        .checked_add(srkit_core::time::Duration384::from_millis(123_456_789))
        .unwrap()
        .modulo_millis(period_ms)
        .unwrap();
    rows.push(("DT-13", (phase as f64 - 123_456_789.0) / period_ms as f64));

    // causal rounding: adversarial arrivals just above a tick, settled by the exact integer check
    rows.push(("DT-14", causal_failures() as f64));

    let total = common::check(&rows);
    println!(
        "srkit-core reproduces {} of the oracle's {total} rows",
        rows.len()
    );
    assert_eq!(rows.len(), 47);
    let _ = units::PARSEC;
}

fn libm_expm1(x: f64) -> f64 {
    x.exp_m1()
}

/// The oracle's DT-14 experiment: a receiver moving at up to 0.9 c, placed so that the signal arrives a hair
/// after a whole tick. `ceil` of the floating-point closed form is sometimes a tick early; settling it with
/// [`causal::settle_arrival`] must never be.
fn causal_failures() -> usize {
    let c_m_per_ms = 299_792.458_f64;
    let mut state = 0x2545_F491_4F6C_DD1D_u64; // xorshift, fixed seed
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let (mut failures, mut naive_early) = (0, 0);
    let emission = Time384::from_decimal("1000000000000000000000000000000000000000").unwrap();
    for _ in 0..100_000 {
        let k = (1e6 + next() * (1e13 - 1e6)).floor();
        let v = (next() * 1.8 - 0.9) * c_m_per_ms; // receiver velocity along the line, m per ms
        let target = k + next() * 4.0 * (k * f64::EPSILON);
        let x0 = target * (c_m_per_ms - v); // receiver's offset at emission, m
        let estimate_ms = (x0 / (c_m_per_ms - v)).ceil();
        let displacement = |tick: Time384| -> Result<GridVec, srkit_core::Error> {
            let dt = tick.signed_duration_since(emission).local_seconds()? * 1e3;
            GridVec::from_metres(x0 + v * dt, 0.0, 0.0)
        };
        let estimate = emission
            .checked_add(srkit_core::time::Duration384::from_unsigned_millis(
                estimate_ms as u64,
            ))
            .unwrap();
        if !causal::arrived(emission, estimate, displacement(estimate).unwrap()) {
            naive_early += 1;
        }
        let settled = causal::settle_arrival(emission, estimate, 8, displacement).unwrap();
        let ok = causal::arrived(emission, settled, displacement(settled).unwrap());
        let earlier = settled
            .checked_sub(srkit_core::time::Duration384::from_unsigned_millis(1))
            .unwrap();
        let minimal = !causal::arrived(emission, earlier, displacement(earlier).unwrap());
        if !(ok && minimal) {
            failures += 1;
        }
    }
    println!("causal rounding: plain ceil() of the float was a tick early {naive_early} times in 100,000; settled exactly: {failures} failures");
    failures
}
