//! Known limits, textbook values and identities between the functions.

mod common;

use common::assert_close;
use core::f64::consts::{FRAC_PI_2, PI};
use srkit_core::units::{
    AU, BOLTZMANN, C, DAY, GM_SUN, G_NEWTON, HBAR, JULIAN_YEAR, M_SUN, PARSEC, R_SUN,
    STEFAN_BOLTZMANN,
};
use srkit_gr::hawking::{self, EmissionModel};
use srkit_gr::kerr::{self, Orbit};
use srkit_gr::{
    gravitational_parameter, gravitational_radius, gravitational_time,
    mass_from_gravitational_parameter, newtonian, schwarzschild, weak_field,
};

/// "Within a few ULP".
const ULPS: f64 = 2e-15;
/// Seconds of arc per radian.
const ARCSEC: f64 = 180.0 * 3600.0 / PI;
/// The Earth's gravitational parameter, m³/s² (IERS Conventions 2010, TCG-compatible value).
const GM_EARTH: f64 = 3.986_004_418e14;
/// A stellar-mass hole for the strong-field tests, kg.
const HOLE: f64 = 10.0 * M_SUN;
const SENSES: [Orbit; 2] = [Orbit::Prograde, Orbit::Retrograde];

// ---------------------------------------------------------------- scales

#[test]
fn solar_scales() {
    assert_close(
        gravitational_parameter(M_SUN),
        GM_SUN,
        2.3e-16,
        "GM round trip",
    );
    assert_close(
        mass_from_gravitational_parameter(GM_SUN),
        M_SUN,
        0.0,
        "M_SUN",
    );
    assert_close(
        gravitational_radius(M_SUN),
        1476.625,
        1e-6,
        "GM/c^2 of the Sun",
    );
    assert_close(
        gravitational_time(M_SUN),
        4.925_49e-6,
        1e-6,
        "GM/c^3 of the Sun",
    );
    assert_close(
        schwarzschild::schwarzschild_radius(M_SUN),
        2953.25,
        1e-6,
        "r_s of the Sun",
    );
    assert_close(
        schwarzschild::schwarzschild_radius(mass_from_gravitational_parameter(GM_EARTH)),
        8.870_06e-3,
        1e-6,
        "r_s of the Earth",
    );
}

// ---------------------------------------------------------------- Schwarzschild

#[test]
fn schwarzschild_characteristic_radii() {
    let m = gravitational_radius(HOLE);
    assert_eq!(schwarzschild::schwarzschild_radius(HOLE), 2.0 * m);
    assert_eq!(schwarzschild::photon_sphere_radius(HOLE), 3.0 * m);
    assert_eq!(schwarzschild::marginally_bound_radius(HOLE), 4.0 * m);
    assert_eq!(schwarzschild::isco_radius(HOLE), 6.0 * m);
    assert_eq!(schwarzschild::schwarzschild_radius(0.0), 0.0);
}

#[test]
fn schwarzschild_static_clock() {
    let rs = schwarzschild::schwarzschild_radius(HOLE);
    assert_eq!(schwarzschild::static_dtau_dt(HOLE, rs), 0.0);
    assert_eq!(schwarzschild::static_dtau_dt(HOLE, f64::INFINITY), 1.0);
    assert_eq!(schwarzschild::static_dtau_dt(0.0, 1.0), 1.0);
    assert!(schwarzschild::static_dtau_dt(HOLE, 0.99 * rs).is_nan());
    // sqrt(1 - 1/2), sqrt(1 - 1/4)
    assert_close(
        schwarzschild::static_dtau_dt(HOLE, 2.0 * rs),
        0.5f64.sqrt(),
        ULPS,
        "2 r_s",
    );
    assert_close(
        schwarzschild::static_dtau_dt(HOLE, 4.0 * rs),
        0.75f64.sqrt(),
        ULPS,
        "4 r_s",
    );
}

#[test]
fn schwarzschild_redshift() {
    let rs = schwarzschild::schwarzschild_radius(HOLE);
    // 1 + z from r to infinity is the reciprocal of the static clock rate.
    for x in [1.001, 1.5, 3.0, 50.0, 1e6] {
        let r = x * rs;
        let z = schwarzschild::gravitational_redshift(HOLE, r, f64::INFINITY);
        assert_close(
            1.0 + z,
            1.0 / schwarzschild::static_dtau_dt(HOLE, r),
            4.0 * ULPS,
            "1 + z to infinity",
        );
    }
    // Light sent back the other way is blueshifted by the reciprocal factor, and redshifts compose.
    let (a, b, c) = (1.3 * rs, 4.0 * rs, 90.0 * rs);
    let z_ab = schwarzschild::gravitational_redshift(HOLE, a, b);
    let z_ba = schwarzschild::gravitational_redshift(HOLE, b, a);
    let z_bc = schwarzschild::gravitational_redshift(HOLE, b, c);
    let z_ac = schwarzschild::gravitational_redshift(HOLE, a, c);
    assert!(z_ab > 0.0 && z_ba < 0.0 && z_ba > -1.0);
    assert_close(
        (1.0 + z_ab) * (1.0 + z_ba),
        1.0,
        4.0 * ULPS,
        "there and back",
    );
    assert_close(
        (1.0 + z_ab) * (1.0 + z_bc),
        1.0 + z_ac,
        4.0 * ULPS,
        "composition",
    );
    assert_eq!(schwarzschild::gravitational_redshift(HOLE, b, b), 0.0);
    assert_eq!(
        schwarzschild::gravitational_redshift(HOLE, rs, b),
        f64::INFINITY
    );
    assert!(schwarzschild::gravitational_redshift(HOLE, 0.5 * rs, b).is_nan());

    // The Sun's surface seen from far away: z = 2.12e-6, the weak-field GM/(R c^2) to first order.
    let z_sun = schwarzschild::gravitational_redshift(M_SUN, R_SUN, f64::INFINITY);
    assert_close(z_sun, 2.1225e-6, 1e-4, "solar redshift");
    assert_close(
        z_sun,
        -weak_field::static_rate_offset(M_SUN, R_SUN),
        4e-6,
        "first order",
    );
    // A small redshift is resolved in full: a 1 m climb at the Earth's surface is g h / c^2 = 1.09e-16.
    let earth = mass_from_gravitational_parameter(GM_EARTH);
    let z_metre = schwarzschild::gravitational_redshift(earth, 6.371e6, 6.371e6 + 1.0);
    let g_h = GM_EARTH / (6.371e6 * 6.371e6) / (C * C);
    assert_close(z_metre, g_h, 1e-6, "1 m at the Earth's surface");
}

#[test]
fn schwarzschild_circular_orbits() {
    let m = gravitational_radius(HOLE);
    assert_close(
        schwarzschild::circular_orbit_dt_dtau(HOLE, 6.0 * m),
        2f64.sqrt(),
        ULPS,
        "ISCO",
    );
    assert_eq!(
        schwarzschild::circular_orbit_dt_dtau(HOLE, 3.0 * m),
        f64::INFINITY
    );
    assert!(schwarzschild::circular_orbit_dt_dtau(HOLE, 2.9 * m).is_nan());
    assert_eq!(
        schwarzschild::circular_orbit_dt_dtau(HOLE, f64::INFINITY),
        1.0
    );
    for x in [3.2, 6.0, 10.0, 1e4] {
        let r = x * m;
        // dt/dtau is the static factor times the Lorentz factor of the locally measured orbital speed,
        // v^2 = GM/(r - r_s).
        let local_speed2 = G_NEWTON * HOLE / (r - 2.0 * m);
        let gamma = 1.0 / (1.0 - local_speed2 / (C * C)).sqrt();
        let expected = gamma / schwarzschild::static_dtau_dt(HOLE, r);
        assert_close(
            schwarzschild::circular_orbit_dt_dtau(HOLE, r),
            expected,
            40.0 * ULPS,
            "static x Lorentz",
        );
        // Kepler's third law holds exactly in the far-away time.
        let omega = schwarzschild::circular_orbit_angular_velocity(HOLE, r);
        assert_close(
            omega * omega * r * r * r,
            G_NEWTON * HOLE,
            2.0 * ULPS,
            "Kepler III",
        );
        assert_close(
            2.0 * PI / omega,
            newtonian::circular_period(HOLE, r),
            2.0 * ULPS,
            "period",
        );
    }
    assert!(schwarzschild::circular_orbit_angular_velocity(HOLE, 2.9 * m).is_nan());
}

// ---------------------------------------------------------------- Kerr at zero spin

#[test]
fn kerr_reduces_to_schwarzschild_at_zero_spin() {
    let m = gravitational_radius(HOLE);
    let rs = schwarzschild::schwarzschild_radius(HOLE);
    assert_eq!(kerr::outer_horizon_radius(HOLE, 0.0), rs);
    assert_eq!(kerr::inner_horizon_radius(HOLE, 0.0), 0.0);
    assert_eq!(kerr::horizon_angular_velocity(HOLE, 0.0), 0.0);
    assert_eq!(kerr::irreducible_mass(HOLE, 0.0), HOLE);
    assert_eq!(kerr::extractable_energy_fraction(0.0), 0.0);
    for theta in [0.0, 0.7, FRAC_PI_2, 3.0] {
        assert_eq!(kerr::ergosphere_radius(HOLE, 0.0, theta), rs);
        for x in [2.0, 2.000_001, 2.7, 6.0, 1e5] {
            let r = x * m;
            let exact = schwarzschild::static_dtau_dt(HOLE, r);
            assert_close(
                kerr::static_dtau_dt(HOLE, 0.0, r, theta),
                exact,
                ULPS,
                "static clock",
            );
            assert_close(
                kerr::zamo_dtau_dt(HOLE, 0.0, r, theta),
                exact,
                ULPS,
                "ZAMO clock",
            );
            assert_eq!(kerr::zamo_angular_velocity(HOLE, 0.0, r, theta), 0.0);
        }
    }
    for orbit in SENSES {
        assert_eq!(
            kerr::photon_orbit_radius(HOLE, 0.0, orbit),
            schwarzschild::photon_sphere_radius(HOLE)
        );
        assert_eq!(
            kerr::marginally_bound_radius(HOLE, 0.0, orbit),
            schwarzschild::marginally_bound_radius(HOLE)
        );
        assert_eq!(
            kerr::isco_radius(HOLE, 0.0, orbit),
            schwarzschild::isco_radius(HOLE)
        );
        assert_close(
            kerr::isco_specific_energy(0.0, orbit),
            (8f64 / 9.0).sqrt(),
            ULPS,
            "E_isco",
        );
        assert_close(
            kerr::accretion_efficiency(0.0, orbit),
            1.0 - (8f64 / 9.0).sqrt(),
            4.0 * ULPS,
            "efficiency",
        );
        assert_close(
            kerr::isco_dt_dtau(0.0, orbit),
            2f64.sqrt(),
            ULPS,
            "dt/dtau on the ISCO",
        );
        assert_close(
            kerr::isco_specific_angular_momentum(HOLE, 0.0, orbit),
            2.0 * 3f64.sqrt() * G_NEWTON * HOLE / C,
            ULPS,
            "L_isco",
        );
        for x in [3.0, 3.5, 6.0, 20.0, 1e6] {
            let r = x * m;
            assert_close(
                kerr::circular_orbit_angular_velocity(HOLE, 0.0, orbit, r),
                schwarzschild::circular_orbit_angular_velocity(HOLE, r),
                ULPS,
                "angular velocity",
            );
            if x > 3.0 {
                assert_close(
                    kerr::circular_orbit_dt_dtau(HOLE, 0.0, orbit, r),
                    schwarzschild::circular_orbit_dt_dtau(HOLE, r),
                    8.0 * ULPS,
                    "dt/dtau",
                );
                // E = (1 - 2m/r)/sqrt(1 - 3m/r), L = sqrt(G M r)/sqrt(1 - 3m/r)
                let root = (1.0 - 3.0 / x).sqrt();
                assert_close(
                    kerr::circular_orbit_specific_energy(HOLE, 0.0, orbit, r),
                    (1.0 - 2.0 / x) / root,
                    8.0 * ULPS,
                    "E",
                );
                assert_close(
                    kerr::circular_orbit_specific_angular_momentum(HOLE, 0.0, orbit, r),
                    (G_NEWTON * HOLE * r).sqrt() / root,
                    8.0 * ULPS,
                    "L",
                );
            }
        }
    }
}

// ---------------------------------------------------------------- Kerr at the extremal limit

#[test]
fn kerr_extremal_limits() {
    let m = gravitational_radius(HOLE);
    let (pro, retro) = (Orbit::Prograde, Orbit::Retrograde);
    for spin in [1.0, -1.0] {
        // ISCO at 1 M prograde and 9 M retrograde, exactly.
        assert_eq!(kerr::isco_radius(HOLE, spin, pro), m);
        assert_eq!(kerr::isco_radius(HOLE, spin, retro), 9.0 * m);
        assert_eq!(kerr::photon_orbit_radius(HOLE, spin, pro), m);
        assert_eq!(kerr::photon_orbit_radius(HOLE, spin, retro), 4.0 * m);
        assert_eq!(kerr::marginally_bound_radius(HOLE, spin, pro), m);
        assert_close(
            kerr::marginally_bound_radius(HOLE, spin, retro),
            (3.0 + 2.0 * 2f64.sqrt()) * m,
            ULPS,
            "r_mb retrograde",
        );
        // Efficiency 1 - 1/sqrt(3) = 42.26% prograde, 1 - sqrt(25/27) = 3.77% retrograde.
        assert_close(
            kerr::accretion_efficiency(spin, pro),
            1.0 - 1.0 / 3f64.sqrt(),
            ULPS,
            "efficiency",
        );
        assert_close(
            kerr::accretion_efficiency(spin, retro),
            1.0 - (25f64 / 27.0).sqrt(),
            8.0 * ULPS,
            "retro",
        );
        assert_close(
            kerr::isco_specific_energy(spin, pro),
            1.0 / 3f64.sqrt(),
            ULPS,
            "E prograde",
        );
        assert_close(
            kerr::isco_specific_energy(spin, retro),
            5.0 / (3.0 * 3f64.sqrt()),
            ULPS,
            "E retrograde",
        );
        let unit = G_NEWTON * HOLE / C;
        assert_close(
            kerr::isco_specific_angular_momentum(HOLE, spin, pro),
            2.0 / 3f64.sqrt() * unit,
            ULPS,
            "L",
        );
        assert_close(
            kerr::isco_specific_angular_momentum(HOLE, spin, retro),
            22.0 / (3.0 * 3f64.sqrt()) * unit,
            ULPS,
            "L retrograde",
        );
        // The prograde ISCO reaches the horizon and its clock stops; the retrograde one stays mild.
        assert_eq!(kerr::isco_dt_dtau(spin, pro), f64::INFINITY);
        assert_close(
            kerr::isco_dt_dtau(spin, retro),
            13.0 / (6.0 * 3f64.sqrt()),
            ULPS,
            "dt/dtau retrograde",
        );
        // The two horizons meet at M, and 1 - 1/sqrt(2) of the hole's mass is spin energy.
        assert_eq!(kerr::outer_horizon_radius(HOLE, spin), m);
        assert_eq!(kerr::inner_horizon_radius(HOLE, spin), m);
        assert_close(
            kerr::irreducible_mass(HOLE, spin),
            HOLE / 2f64.sqrt(),
            ULPS,
            "M_irr",
        );
        assert_close(
            kerr::extractable_energy_fraction(spin),
            1.0 - 1.0 / 2f64.sqrt(),
            ULPS,
            "29.29%",
        );
        assert_close(
            kerr::horizon_angular_velocity(HOLE, spin),
            spin * 0.5 / gravitational_time(HOLE),
            ULPS,
            "Omega_H",
        );
    }
    assert_close(
        kerr::extractable_energy_fraction(1.0),
        0.292_893_218_813_452_5,
        ULPS,
        "29.29%",
    );
    assert_close(
        kerr::accretion_efficiency(1.0, pro),
        0.422_649_730_810_374_2,
        ULPS,
        "42.26%",
    );
}

#[test]
fn kerr_tabulated_values() {
    // Thorne's limit a* = 0.998: r_isco = 1.2370 M, efficiency 32.1%; a* = 0.9: 2.3209 M, 15.58%.
    let m = gravitational_radius(HOLE);
    assert_close(
        kerr::isco_radius(HOLE, 0.998, Orbit::Prograde) / m,
        1.236_97,
        2e-6,
        "r_isco(0.998)",
    );
    assert_close(
        kerr::accretion_efficiency(0.998, Orbit::Prograde),
        0.320_99,
        2e-5,
        "efficiency(0.998)",
    );
    assert_close(
        kerr::isco_radius(HOLE, 0.9, Orbit::Prograde) / m,
        2.320_88,
        2e-6,
        "r_isco(0.9)",
    );
    assert_close(
        kerr::isco_radius(HOLE, 0.9, Orbit::Retrograde) / m,
        8.717_35,
        2e-6,
        "r_isco(0.9) retrograde",
    );
    assert_close(
        kerr::accretion_efficiency(0.9, Orbit::Prograde),
        0.155_75,
        2e-5,
        "efficiency(0.9)",
    );
    // Equatorial static limit at 2M for every spin; on the axis it meets the horizon.
    for spin in [0.0, 0.3, -0.9, 1.0] {
        assert_eq!(kerr::ergosphere_radius(HOLE, spin, FRAC_PI_2), 2.0 * m);
        assert_eq!(
            kerr::ergosphere_radius(HOLE, spin, 0.0),
            kerr::outer_horizon_radius(HOLE, spin)
        );
    }
    // Small spin: extractable fraction a^2/8, inner horizon a^2 M/2, without loss of digits.
    assert_close(
        kerr::extractable_energy_fraction(1e-8),
        1e-16 / 8.0,
        1e-15,
        "a^2/8",
    );
    assert_close(
        kerr::inner_horizon_radius(HOLE, 1e-8),
        0.5e-16 * m,
        1e-15,
        "a^2 M/2",
    );
}

#[test]
fn kerr_sign_of_the_spin() {
    let m = gravitational_radius(HOLE);
    let r = 3.7 * m;
    for spin in [0.2, 0.7, 0.998] {
        // Reversing the hole reverses the angular velocities and nothing else.
        assert_eq!(
            kerr::horizon_angular_velocity(HOLE, -spin),
            -kerr::horizon_angular_velocity(HOLE, spin)
        );
        assert_eq!(
            kerr::zamo_angular_velocity(HOLE, -spin, r, 1.1),
            -kerr::zamo_angular_velocity(HOLE, spin, r, 1.1)
        );
        assert!(kerr::horizon_angular_velocity(HOLE, spin) > 0.0);
        assert_eq!(
            kerr::outer_horizon_radius(HOLE, -spin),
            kerr::outer_horizon_radius(HOLE, spin)
        );
        assert_eq!(
            kerr::zamo_dtau_dt(HOLE, -spin, r, 1.1),
            kerr::zamo_dtau_dt(HOLE, spin, r, 1.1)
        );
        assert_eq!(
            kerr::static_dtau_dt(HOLE, -spin, r, 1.1),
            kerr::static_dtau_dt(HOLE, spin, r, 1.1)
        );
        assert_eq!(
            kerr::extractable_energy_fraction(-spin),
            kerr::extractable_energy_fraction(spin)
        );
        for orbit in SENSES {
            // "Prograde" always means "with the hole", whichever way the hole turns.
            assert_eq!(
                kerr::isco_radius(HOLE, -spin, orbit),
                kerr::isco_radius(HOLE, spin, orbit)
            );
            assert_eq!(
                kerr::isco_dt_dtau(-spin, orbit),
                kerr::isco_dt_dtau(spin, orbit)
            );
            assert_eq!(
                kerr::accretion_efficiency(-spin, orbit),
                kerr::accretion_efficiency(spin, orbit)
            );
        }
        assert!(kerr::isco_radius(HOLE, spin, Orbit::Prograde) < 6.0 * m);
        assert!(kerr::isco_radius(HOLE, spin, Orbit::Retrograde) > 6.0 * m);
    }
}

#[test]
fn kerr_frame_dragging() {
    let m = gravitational_radius(HOLE);
    for spin in [0.3, 0.9, 0.998, 1.0] {
        let horizon = kerr::outer_horizon_radius(HOLE, spin);
        let omega_h = kerr::horizon_angular_velocity(HOLE, spin);
        for theta in [0.0, 0.4, FRAC_PI_2] {
            // On the horizon every ZAMO turns at the horizon's rate and its clock stops.
            assert_close(
                kerr::zamo_angular_velocity(HOLE, spin, horizon, theta),
                omega_h,
                4.0 * ULPS,
                "on the horizon",
            );
            assert_eq!(kerr::zamo_dtau_dt(HOLE, spin, horizon, theta), 0.0);
            assert!(kerr::zamo_dtau_dt(HOLE, spin, 0.999 * horizon, theta).is_nan());
            assert!(kerr::zamo_angular_velocity(HOLE, spin, 0.999 * horizon, theta).is_nan());
            // Far away: the Lense-Thirring rate 2 G J/(c^2 r^3) with J = a* G M^2/c.
            let r = 1e7 * m;
            let j = spin * G_NEWTON * HOLE * HOLE / C;
            let lense_thirring = 2.0 * G_NEWTON * j / (C * C * r * r * r);
            assert_close(
                kerr::zamo_angular_velocity(HOLE, spin, r, theta),
                lense_thirring,
                1e-6,
                "Lense-Thirring",
            );
            assert_eq!(
                kerr::zamo_angular_velocity(HOLE, spin, f64::INFINITY, theta),
                0.0
            );
            assert_eq!(kerr::zamo_dtau_dt(HOLE, spin, f64::INFINITY, theta), 1.0);
            // The static clock stops on the static limit and does not exist inside it.
            let limit = kerr::ergosphere_radius(HOLE, spin, theta);
            assert_eq!(kerr::static_dtau_dt(HOLE, spin, limit, theta), 0.0);
            assert!(kerr::static_dtau_dt(HOLE, spin, 0.999 * limit, theta).is_nan());
        }
        // Inside the ergosphere a ZAMO still exists and its clock runs.
        let inside = 0.5 * (horizon + 2.0 * m);
        assert!(kerr::static_dtau_dt(HOLE, spin, inside, FRAC_PI_2).is_nan());
        assert!(kerr::zamo_dtau_dt(HOLE, spin, inside, FRAC_PI_2) > 0.0);
    }
}

#[test]
fn kerr_isco_closed_forms_match_the_general_orbit_formulas() {
    for spin in [0.0, 0.1, 0.5, 0.9, 0.998, 0.999_999] {
        for orbit in SENSES {
            let r = kerr::isco_radius(HOLE, spin, orbit);
            // The general formulas lose digits where 1 - 3w^2 + 2aw^3 is small (near-extremal prograde).
            let slack = if spin > 0.99 && orbit == Orbit::Prograde {
                1e-9
            } else {
                40.0 * ULPS
            };
            assert_close(
                kerr::circular_orbit_specific_energy(HOLE, spin, orbit, r),
                kerr::isco_specific_energy(spin, orbit),
                slack,
                "E",
            );
            assert_close(
                kerr::circular_orbit_specific_angular_momentum(HOLE, spin, orbit, r),
                kerr::isco_specific_angular_momentum(HOLE, spin, orbit),
                slack,
                "L",
            );
            assert_close(
                kerr::circular_orbit_dt_dtau(HOLE, spin, orbit, r),
                kerr::isco_dt_dtau(spin, orbit),
                slack,
                "dt/dtau",
            );
            // The ISCO is where the energy and the angular momentum of circular orbits are least.
            let e_isco = kerr::isco_specific_energy(spin, orbit);
            let l_isco = kerr::isco_specific_angular_momentum(HOLE, spin, orbit);
            for factor in [0.98, 1.02, 1.5, 4.0] {
                let other = factor * r;
                if other <= kerr::photon_orbit_radius(HOLE, spin, orbit) {
                    continue;
                }
                assert!(kerr::circular_orbit_specific_energy(HOLE, spin, orbit, other) > e_isco);
                assert!(
                    kerr::circular_orbit_specific_angular_momentum(HOLE, spin, orbit, other)
                        > l_isco
                );
            }
            // The marginally bound orbit has E = 1, and the radii are ordered.
            let r_mb = kerr::marginally_bound_radius(HOLE, spin, orbit);
            assert_close(
                kerr::circular_orbit_specific_energy(HOLE, spin, orbit, r_mb),
                1.0,
                1e-9,
                "E(r_mb)",
            );
            let r_ph = kerr::photon_orbit_radius(HOLE, spin, orbit);
            assert!(kerr::outer_horizon_radius(HOLE, spin) <= r_ph && r_ph < r_mb && r_mb < r);
            // No circular orbit inside the photon orbit; just outside it the clock nearly stops.
            assert!(kerr::circular_orbit_dt_dtau(HOLE, spin, orbit, 0.999 * r_ph).is_nan());
            assert!(kerr::circular_orbit_dt_dtau(HOLE, spin, orbit, 1.000_001 * r_ph) > 100.0);
        }
    }
}

// ---------------------------------------------------------------- weak field

#[test]
fn shapiro_delay_at_superior_conjunction() {
    // A signal between the Earth and Mars that grazes the Sun: about 250 microseconds for the round trip.
    let (r_earth, r_mars) = (AU, 1.523_679 * AU);
    let one_way = weak_field::shapiro_delay(M_SUN, r_earth, r_mars, R_SUN);
    let round_trip = 2.0 * one_way;
    assert!(
        (240e-6..255e-6).contains(&round_trip),
        "round trip {round_trip:e} s"
    );
    assert_close(round_trip, 247.2e-6, 5e-4, "round trip");
    // The familiar form (4GM/c^3) ln(4 r_e r_m / b^2) for the round trip, good to b^2/r^2.
    let familiar = 4.0 * GM_SUN / (C * C * C) * (4.0 * r_earth * r_mars / (R_SUN * R_SUN)).ln();
    assert_close(round_trip, familiar, 2e-6, "log form");
    // Two legs equal the two-point form when the closest approach lies between the ends.
    let separation =
        (r_earth * r_earth - R_SUN * R_SUN).sqrt() + (r_mars * r_mars - R_SUN * R_SUN).sqrt();
    let between = weak_field::shapiro_delay_between(M_SUN, r_earth, r_mars, separation);
    assert_close(between, one_way, 1e-9, "two-point form");
    // A leg of no length has no delay, and the delay falls as the path moves away from the mass.
    assert_eq!(weak_field::shapiro_delay_leg(M_SUN, R_SUN, R_SUN), 0.0);
    assert!(weak_field::shapiro_delay(M_SUN, r_earth, r_mars, 2.0 * R_SUN) < one_way);
    assert!(weak_field::shapiro_delay_leg(M_SUN, 0.5 * R_SUN, R_SUN).is_nan());
    // The separation must make a triangle with the two distances.
    assert!(weak_field::shapiro_delay_between(M_SUN, r_earth, r_mars, 3.0 * AU).is_nan());
    assert!(weak_field::shapiro_delay_between(M_SUN, r_earth, r_mars, 0.4 * AU).is_nan());
    assert_eq!(
        weak_field::shapiro_delay_between(M_SUN, r_earth, r_earth, 0.0),
        0.0
    );
}

#[test]
fn light_deflection_at_the_solar_limb() {
    let arcsec = weak_field::deflection_angle(M_SUN, R_SUN) * ARCSEC;
    assert_close(arcsec, 1.75, 2e-3, "deflection at the limb");
    assert_close(
        arcsec,
        1.751_2,
        5e-5,
        "deflection at the limb, five figures",
    );
    // Inversely proportional to the impact parameter; twice the Newtonian value 2GM/(b c^2).
    assert_close(
        weak_field::deflection_angle(M_SUN, 3.0 * R_SUN),
        arcsec / ARCSEC / 3.0,
        ULPS,
        "1/b",
    );
    assert_close(
        arcsec / ARCSEC,
        2.0 * (2.0 * GM_SUN / (R_SUN * C * C)),
        ULPS,
        "twice Newtonian",
    );
    assert!(weak_field::deflection_angle(M_SUN, 0.0).is_nan());
}

#[test]
fn einstein_ring() {
    // One solar mass halfway to a source 8 kpc away: about 1 milliarcsecond.
    let kpc = 1e3 * PARSEC;
    let (d_lens, d_source) = (4.0 * kpc, 8.0 * kpc);
    let theta = weak_field::einstein_angle(M_SUN, d_lens, d_source, d_source - d_lens);
    assert_close(theta * ARCSEC * 1e3, 1.009, 1e-3, "Einstein angle in mas");
    // The ring is where the deflection brings a ray from a source behind the lens back to the observer:
    // theta_E d_s = alpha(theta_E d_l) d_ls.
    let bend = weak_field::deflection_angle(M_SUN, theta * d_lens);
    assert_close(
        bend * (d_source - d_lens) / d_source,
        theta,
        4.0 * ULPS,
        "lens equation",
    );
    // It grows as the square root of the mass, and vanishes for a source at the lens.
    assert_close(
        weak_field::einstein_angle(4.0 * M_SUN, d_lens, d_source, d_source - d_lens),
        2.0 * theta,
        ULPS,
        "sqrt(M)",
    );
    assert_eq!(weak_field::einstein_angle(M_SUN, d_lens, d_lens, 0.0), 0.0);
    // The solar gravitational lens: light grazing the Sun from a distant source converges at
    // R_sun/theta = R_sun^2 c^2/(4GM) = 548 au.
    let focus = R_SUN / weak_field::deflection_angle(M_SUN, R_SUN);
    assert_close(focus / AU, 547.8, 1e-3, "solar focal distance");
}

#[test]
fn weak_field_clock_offsets_are_the_first_order_of_schwarzschild() {
    let earth = mass_from_gravitational_parameter(GM_EARTH);
    for r in [6.371e6, 2.656_175e7, 4.2164e7, 3.844e8] {
        // The exact rates are 1 - 7e-10 and so on, and binary64 resolves a number near 1 only to 1.1e-16.
        // That, not the second-order term (below 1e-18 here), limits the comparison.
        let resolution = 4.0 * f64::EPSILON;
        let exact_static = schwarzschild::static_dtau_dt(earth, r) - 1.0;
        assert!(
            (weak_field::static_rate_offset(earth, r) - exact_static).abs() < resolution,
            "static at {r}"
        );
        let exact_orbit = 1.0 / schwarzschild::circular_orbit_dt_dtau(earth, r) - 1.0;
        assert!(
            (weak_field::circular_orbit_rate_offset(earth, r) - exact_orbit).abs() < resolution,
            "orbit at {r}"
        );
        // gravitational two thirds, velocity one third
        let v = newtonian::circular_speed(earth, r);
        let sum = weak_field::static_rate_offset(earth, r) - 0.5 * v * v / (C * C);
        assert_close(
            weak_field::circular_orbit_rate_offset(earth, r),
            sum,
            2.0 * ULPS,
            "potential + speed",
        );
    }
    assert_close(
        weak_field::static_rate_offset(M_SUN, R_SUN),
        -2.1225e-6,
        1e-4,
        "solar surface",
    );
    assert_close(
        weak_field::static_rate_offset(earth, 6.371e6),
        -6.961e-10,
        1e-4,
        "Earth's surface",
    );
}

#[test]
fn gps_clock_offset() {
    let earth = mass_from_gravitational_parameter(GM_EARTH);
    let orbit = 2.656_175e7; // GPS semi-major axis, m
    let equator = 6.378_137e6; // WGS 84 equatorial radius, m
    let spin_speed = 465.1; // equatorial rotation speed, m/s
    let per_day = DAY * 1e6; // fractional offset to microseconds per day

    // Point-mass Earth, clock on the rotating equator: +38.5 us/day.
    let net = weak_field::orbit_vs_surface_rate_offset(earth, orbit, equator, spin_speed);
    assert_close(net, 4.461e-10, 2e-4, "fractional offset");
    assert_close(net * per_day, 38.54, 5e-4, "microseconds per day");
    assert!(
        (net * per_day - 38.6).abs() < 0.1,
        "net {} us/day",
        net * per_day
    );

    // The two classic pieces for a clock at rest on the surface: +45.7 us/day from the potential and
    // -7.2 us/day from the orbital speed.
    let potential = weak_field::static_rate_offset(earth, orbit)
        - weak_field::static_rate_offset(earth, equator);
    let speed = weak_field::circular_orbit_rate_offset(earth, orbit)
        - weak_field::static_rate_offset(earth, orbit);
    assert_close(potential * per_day, 45.65, 5e-4, "gravitational part");
    assert_close(speed * per_day, -7.21, 1e-3, "velocity part");
    assert_close(
        weak_field::orbit_vs_surface_rate_offset(earth, orbit, equator, 0.0),
        potential + speed,
        4.0 * ULPS,
        "sum of the parts",
    );

    // With the potential of the geoid, W0 = 62,636,856.0 m^2/s^2 (which adds the Earth's oblateness and
    // rotation), the sum is the published 4.4647e-10: +38.6 us/day, the offset built into GPS clocks.
    let geoid = 62_636_856.0 / (C * C);
    let published = weak_field::circular_orbit_rate_offset(earth, orbit) + geoid;
    assert_close(published, 4.4647e-10, 2e-5, "Ashby's factory offset");
    assert_close(published * per_day, 38.6, 1e-3, "+38.6 us/day");

    // A clock in orbit at 3R/2 keeps time with one on a non-rotating surface.
    assert!(
        weak_field::orbit_vs_surface_rate_offset(earth, 1.5 * equator, equator, 0.0).abs() < 1e-24
    );
    assert!(weak_field::orbit_vs_surface_rate_offset(earth, 1.4 * equator, equator, 0.0) < 0.0);
    assert!(weak_field::orbit_vs_surface_rate_offset(earth, 1.6 * equator, equator, 0.0) > 0.0);
}

// ---------------------------------------------------------------- Hawking

#[test]
fn hawking_temperature_and_entropy() {
    // T M is a constant: 1.2270e23 K kg.
    let t_m = HBAR * C * C * C / (8.0 * PI * G_NEWTON * BOLTZMANN);
    for mass in [1e5, 1e9, M_SUN, 1e11 * M_SUN] {
        assert_close(hawking::temperature(mass) * mass, t_m, ULPS, "T M");
    }
    assert_close(t_m, 1.227e23, 1e-3, "T M in K kg");
    assert_close(
        hawking::temperature(M_SUN),
        6.170e-8,
        1e-3,
        "solar-mass temperature",
    );
    // S/k_B = 4 pi (M/m_Planck)^2, 1.05e77 for the Sun; S = 1.45e54 J/K.
    let planck_mass = (HBAR * C / G_NEWTON).sqrt();
    let s = hawking::entropy(M_SUN);
    assert_close(
        s / BOLTZMANN,
        4.0 * PI * (M_SUN / planck_mass).powi(2),
        4.0 * ULPS,
        "S/k_B",
    );
    assert_close(s, 1.448e54, 1e-3, "solar-mass entropy");
    // First law: dE = T dS, checked with a small step in mass.
    let (mass, step) = (1e12, 1e6);
    let d_entropy = hawking::entropy(mass + step) - hawking::entropy(mass - step);
    let temperature = hawking::temperature(mass);
    assert_close(
        temperature * d_entropy,
        2.0 * step * C * C,
        1e-9,
        "dE = T dS",
    );
    assert!(hawking::temperature(0.0).is_nan() && hawking::temperature(-1.0).is_nan());
}

#[test]
fn hawking_photon_blackbody_is_the_stefan_boltzmann_law() {
    for mass in [1e9, 3.3e14, M_SUN] {
        let rs = schwarzschild::schwarzschild_radius(mass);
        let area = 4.0 * PI * rs * rs;
        let stefan_boltzmann = STEFAN_BOLTZMANN * area * hawking::temperature(mass).powi(4);
        assert_close(
            hawking::power(mass, EmissionModel::PhotonBlackbody),
            stefan_boltzmann,
            8.0 * ULPS,
            "sigma A T^4",
        );
        // and the familiar closed forms
        let closed = HBAR * C.powi(6) / (15360.0 * PI * (G_NEWTON * mass).powi(2));
        assert_close(
            hawking::power(mass, EmissionModel::PhotonBlackbody),
            closed,
            4.0 * ULPS,
            "P",
        );
        let life = 5120.0 * PI * G_NEWTON.powi(2) * mass.powi(3) / (HBAR * C.powi(4));
        assert_close(
            hawking::lifetime(mass, EmissionModel::PhotonBlackbody),
            life,
            4.0 * ULPS,
            "t",
        );
    }
    assert_close(
        hawking::lifetime(M_SUN, EmissionModel::PhotonBlackbody) / JULIAN_YEAR,
        2.0955e67,
        1e-4,
        "yr",
    );
}

#[test]
fn hawking_models_are_consistent() {
    let models = [
        EmissionModel::PhotonBlackbody,
        EmissionModel::PagePhotonsGravitons,
        EmissionModel::MacGibbon { f: 1.0 },
        EmissionModel::MacGibbon { f: 12.0 },
        EmissionModel::Coefficient(2.011e-4),
    ];
    for model in models {
        let alpha = model.power_coefficient();
        assert!(alpha > 0.0);
        for mass in [1e6, 1e9, 1e20, M_SUN] {
            let power = hawking::power(mass, model);
            let rate = hawking::mass_loss_rate(mass, model);
            let life = hawking::lifetime(mass, model);
            assert_close(power, rate * C * C, 2.0 * ULPS, "P = -dM/dt c^2");
            // dM/dt is proportional to 1/M^2, so t = M/(3 |dM/dt|) and the mean power is three times
            // the initial power.
            assert_close(life, mass / (3.0 * rate), 2.0 * ULPS, "t = M/(3 dM/dt)");
            assert_close(
                hawking::mass_for_lifetime(life, model),
                mass,
                2.0 * ULPS,
                "inverse",
            );
            // every model scales the photon black body by the ratio of the coefficients
            let photons = hawking::power(mass, EmissionModel::PhotonBlackbody);
            assert_close(
                power / photons,
                alpha / hawking::PHOTON_BLACKBODY_COEFFICIENT,
                4.0 * ULPS,
                "ratio",
            );
        }
    }
    assert_close(
        hawking::PHOTON_BLACKBODY_COEFFICIENT,
        2.0723e-5,
        1e-4,
        "1/(15360 pi)",
    );
    assert_close(
        hawking::PAGE_PHOTONS_GRAVITONS_COEFFICIENT,
        3.740e-5,
        1e-3,
        "Page photons + gravitons",
    );
    assert_close(
        EmissionModel::MacGibbon { f: 1.0 }.power_coefficient(),
        2.793e-4,
        1e-3,
        "MacGibbon f = 1",
    );
    assert_eq!(EmissionModel::Coefficient(0.25).power_coefficient(), 0.25);
    // A model without a positive coefficient has no power.
    assert!(hawking::power(1e9, EmissionModel::Coefficient(0.0)).is_nan());
    assert!(hawking::lifetime(1e9, EmissionModel::MacGibbon { f: -1.0 }).is_nan());
    assert!(hawking::power(0.0, EmissionModel::PhotonBlackbody).is_nan());
    assert_eq!(hawking::lifetime(0.0, EmissionModel::PhotonBlackbody), 0.0);
    // Examples quoted in the documentation: a hole that lasts a year, and one that lasts 13.8 Gyr.
    let photons = EmissionModel::PhotonBlackbody;
    assert_close(
        hawking::mass_for_lifetime(JULIAN_YEAR, photons),
        7.21e7,
        2e-3,
        "one year",
    );
    assert_close(
        hawking::mass_for_lifetime(13.8e9 * JULIAN_YEAR, photons),
        1.73e11,
        3e-3,
        "13.8 Gyr",
    );
    // As warm as the microwave background today (2.7255 K): 4.5e22 kg.
    assert_close(hawking::temperature(4.50e22), 2.7255, 3e-3, "T = T_CMB");
}

// ---------------------------------------------------------------- Newtonian

#[test]
fn newtonian_orbits() {
    // The Earth round the Sun.
    assert_close(
        newtonian::circular_speed(M_SUN, AU),
        29_784.7,
        1e-5,
        "orbital speed",
    );
    assert_close(
        newtonian::circular_period(M_SUN, AU) / DAY,
        365.2569,
        1e-6,
        "Gaussian year",
    );
    let earth = mass_from_gravitational_parameter(GM_EARTH);
    assert_close(
        newtonian::escape_speed(earth, 6.371e6),
        11_186.0,
        1e-4,
        "escape from the Earth",
    );
    for (mass, r) in [(M_SUN, AU), (earth, 6.778e6), (7.342e22, 1.7374e6)] {
        let v = newtonian::circular_speed(mass, r);
        assert_close(
            newtonian::escape_speed(mass, r),
            2f64.sqrt() * v,
            2.0 * ULPS,
            "sqrt(2) v_circ",
        );
        assert_close(
            newtonian::circular_period(mass, r) * v,
            2.0 * PI * r,
            2.0 * ULPS,
            "T v = 2 pi r",
        );
        // binding energy of a uniform sphere = 0.3 M v_esc^2
        let escape = newtonian::escape_speed(mass, r);
        assert_close(
            newtonian::binding_energy_uniform_sphere(mass, r),
            0.3 * mass * escape * escape,
            4.0 * ULPS,
            "0.3 M v_esc^2",
        );
    }
    // The Newtonian escape speed reaches c at the Schwarzschild radius.
    assert_close(
        newtonian::escape_speed(M_SUN, schwarzschild::schwarzschild_radius(M_SUN)),
        C,
        2.0 * ULPS,
        "c",
    );
    assert!(newtonian::circular_speed(M_SUN, 0.0).is_nan());
    assert!(newtonian::circular_period(0.0, AU).is_nan());
}

#[test]
fn newtonian_synodic_period() {
    let (earth, mars) = (365.256_363, 686.980);
    assert_close(
        newtonian::synodic_period(earth, mars),
        779.936,
        1e-6,
        "Earth and Mars, days",
    );
    assert_close(
        newtonian::synodic_period(earth, mars),
        newtonian::synodic_period(mars, earth),
        ULPS,
        "symmetric",
    );
    assert_eq!(newtonian::synodic_period(earth, earth), f64::INFINITY);
    // a very slow outer body: the synodic period tends to the inner period
    assert_close(
        newtonian::synodic_period(earth, 1e12),
        earth,
        1e-9,
        "distant body",
    );
    assert!(newtonian::synodic_period(0.0, mars).is_nan());
}

#[test]
fn hohmann_transfer_from_the_earth_to_mars() {
    let (r_earth, r_mars) = (AU, 1.523_679 * AU);
    let time = newtonian::hohmann_transfer_time(M_SUN, r_earth, r_mars);
    assert_close(time / DAY, 258.86, 1e-4, "transfer time in days");
    assert_eq!(
        time,
        newtonian::hohmann_transfer_time(M_SUN, r_mars, r_earth)
    );
    // between half the inner period and half the outer one
    assert!(time > 0.5 * newtonian::circular_period(M_SUN, r_earth));
    assert!(time < 0.5 * newtonian::circular_period(M_SUN, r_mars));

    let out = newtonian::hohmann_delta_v(M_SUN, r_earth, r_mars);
    assert_close(out.departure, 2_945.0, 1e-3, "departure burn");
    assert_close(out.arrival, 2_649.0, 1e-3, "arrival burn");
    assert_close(out.total(), 5_594.0, 1e-3, "total");
    // The return trip is the same manoeuvre reversed: both burns retrograde, swapped.
    let back = newtonian::hohmann_delta_v(M_SUN, r_mars, r_earth);
    assert!(back.departure < 0.0 && back.arrival < 0.0);
    assert_close(back.departure, -out.arrival, 4.0 * ULPS, "reversed arrival");
    assert_close(
        back.arrival,
        -out.departure,
        4.0 * ULPS,
        "reversed departure",
    );
    assert_close(back.total(), out.total(), 4.0 * ULPS, "same total");

    // The burns follow from the vis-viva equation on the transfer ellipse.
    let gm = G_NEWTON * M_SUN;
    let a = 0.5 * (r_earth + r_mars);
    let perihelion = (gm * (2.0 / r_earth - 1.0 / a)).sqrt();
    let aphelion = (gm * (2.0 / r_mars - 1.0 / a)).sqrt();
    assert_close(
        out.departure,
        perihelion - newtonian::circular_speed(M_SUN, r_earth),
        1e-12,
        "vis-viva",
    );
    assert_close(
        out.arrival,
        newtonian::circular_speed(M_SUN, r_mars) - aphelion,
        1e-12,
        "vis-viva",
    );

    // No transfer, no burn; and a climb of 1 km from a 400 km Earth orbit keeps its digits.
    let none = newtonian::hohmann_delta_v(M_SUN, r_earth, r_earth);
    assert_eq!(
        (none.departure, none.arrival, none.total()),
        (0.0, 0.0, 0.0)
    );
    let earth = mass_from_gravitational_parameter(GM_EARTH);
    let small = newtonian::hohmann_delta_v(earth, 6.778e6, 6.779e6);
    // each burn is v dr/(4r) to first order
    let first_order = newtonian::circular_speed(earth, 6.778e6) * 1e3 / (4.0 * 6.778e6);
    assert_close(small.departure, first_order, 2e-4, "small climb");
    assert_close(small.arrival, first_order, 4e-4, "small climb");
    assert!(newtonian::hohmann_delta_v(M_SUN, -1.0, AU)
        .departure
        .is_nan());
    assert!(newtonian::hohmann_transfer_time(0.0, AU, AU).is_nan());
}

#[test]
fn uniform_sphere_binding_energy() {
    // The Earth as a uniform sphere: 2.24e32 J.
    assert_close(
        newtonian::binding_energy_uniform_sphere(5.9722e24, 6.3710e6),
        2.242e32,
        1e-3,
        "Earth",
    );
    assert_eq!(newtonian::binding_energy_uniform_sphere(0.0, 1.0), 0.0);
    assert!(newtonian::binding_energy_uniform_sphere(1.0, 0.0).is_nan());
}

// ---------------------------------------------------------------- domains

#[test]
fn out_of_domain_arguments_give_nan() {
    let m = gravitational_radius(HOLE);
    for bad_spin in [1.000_000_1, -1.5, f64::NAN, f64::INFINITY] {
        assert!(kerr::outer_horizon_radius(HOLE, bad_spin).is_nan());
        assert!(kerr::inner_horizon_radius(HOLE, bad_spin).is_nan());
        assert!(kerr::ergosphere_radius(HOLE, bad_spin, 1.0).is_nan());
        assert!(kerr::horizon_angular_velocity(HOLE, bad_spin).is_nan());
        assert!(kerr::irreducible_mass(HOLE, bad_spin).is_nan());
        assert!(kerr::extractable_energy_fraction(bad_spin).is_nan());
        assert!(kerr::static_dtau_dt(HOLE, bad_spin, 10.0 * m, 1.0).is_nan());
        assert!(kerr::zamo_angular_velocity(HOLE, bad_spin, 10.0 * m, 1.0).is_nan());
        assert!(kerr::zamo_dtau_dt(HOLE, bad_spin, 10.0 * m, 1.0).is_nan());
        for orbit in SENSES {
            assert!(kerr::photon_orbit_radius(HOLE, bad_spin, orbit).is_nan());
            assert!(kerr::marginally_bound_radius(HOLE, bad_spin, orbit).is_nan());
            assert!(kerr::isco_radius(HOLE, bad_spin, orbit).is_nan());
            assert!(kerr::isco_specific_energy(bad_spin, orbit).is_nan());
            assert!(kerr::isco_specific_angular_momentum(HOLE, bad_spin, orbit).is_nan());
            assert!(kerr::isco_dt_dtau(bad_spin, orbit).is_nan());
            assert!(kerr::accretion_efficiency(bad_spin, orbit).is_nan());
            assert!(kerr::circular_orbit_dt_dtau(HOLE, bad_spin, orbit, 10.0 * m).is_nan());
            assert!(kerr::circular_orbit_specific_energy(HOLE, bad_spin, orbit, 10.0 * m).is_nan());
            assert!(kerr::circular_orbit_specific_angular_momentum(
                HOLE,
                bad_spin,
                orbit,
                10.0 * m
            )
            .is_nan());
            assert!(
                kerr::circular_orbit_angular_velocity(HOLE, bad_spin, orbit, 10.0 * m).is_nan()
            );
        }
    }
    for bad_mass in [-1.0, f64::NAN] {
        assert!(schwarzschild::schwarzschild_radius(bad_mass).is_nan());
        assert!(schwarzschild::static_dtau_dt(bad_mass, 1e9).is_nan());
        assert!(schwarzschild::gravitational_redshift(bad_mass, 1e9, 2e9).is_nan());
        assert!(schwarzschild::circular_orbit_dt_dtau(bad_mass, 1e9).is_nan());
        assert!(schwarzschild::circular_orbit_angular_velocity(bad_mass, 1e9).is_nan());
        assert!(kerr::isco_radius(bad_mass, 0.5, Orbit::Prograde).is_nan());
        assert!(kerr::zamo_dtau_dt(bad_mass, 0.5, 1e9, 1.0).is_nan());
        assert!(weak_field::shapiro_delay(bad_mass, 2.0, 3.0, 1.0).is_nan());
        assert!(weak_field::deflection_angle(bad_mass, 1.0).is_nan());
        assert!(weak_field::einstein_angle(bad_mass, 1.0, 2.0, 1.0).is_nan());
        assert!(weak_field::static_rate_offset(bad_mass, 1.0).is_nan());
        assert!(weak_field::orbit_vs_surface_rate_offset(bad_mass, 2.0, 1.0, 0.0).is_nan());
        assert!(hawking::temperature(bad_mass).is_nan());
        assert!(hawking::entropy(bad_mass).is_nan());
        assert!(hawking::lifetime(bad_mass, EmissionModel::PhotonBlackbody).is_nan());
        assert!(newtonian::binding_energy_uniform_sphere(bad_mass, 1.0).is_nan());
        assert!(newtonian::circular_speed(bad_mass, 1.0).is_nan());
        assert!(newtonian::hohmann_transfer_time(bad_mass, 1.0, 2.0).is_nan());
    }
    // Circular orbits do not exist inside the photon orbit, including far inside the horizon, where the
    // expression under the square root is positive again.
    for orbit in SENSES {
        for x in [0.05, 0.5, 1.2] {
            assert!(
                kerr::circular_orbit_dt_dtau(HOLE, 0.5, orbit, x * m).is_nan(),
                "x = {x}"
            );
            assert!(kerr::circular_orbit_specific_energy(HOLE, 0.5, orbit, x * m).is_nan());
            assert!(kerr::circular_orbit_angular_velocity(HOLE, 0.5, orbit, x * m).is_nan());
        }
    }
    // Inside the inner horizon Delta is positive again; the ZAMO functions still refuse.
    assert!(kerr::zamo_dtau_dt(HOLE, 0.9, 0.2 * m, 1.0).is_nan());
}
