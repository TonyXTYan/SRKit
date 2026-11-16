//! Properties that hold for every argument: monotonicity, bounds, identities and totality.

use core::f64::consts::PI;
use proptest::prelude::*;
use srkit_core::units::{C, G_NEWTON, M_SUN};
use srkit_gr::hawking::{self, EmissionModel};
use srkit_gr::kerr::{self, Orbit};
use srkit_gr::{gravitational_radius, newtonian, schwarzschild, weak_field};

/// Relative slack for comparisons that hold exactly in real arithmetic: a few ULP.
const SLACK: f64 = 2e-15;
const HOLE: f64 = 10.0 * M_SUN;

/// `a ≤ b` up to rounding.
fn le(a: f64, b: f64) -> bool {
    a <= b + SLACK * b.abs().max(a.abs())
}

/// `a = b` to a relative tolerance.
fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs())
}

/// `10^e` for `e` uniform in a range: a log-uniform positive number.
fn decades(low: f64, high: f64) -> impl Strategy<Value = f64> {
    (low..high).prop_map(|e| 10f64.powf(e))
}

/// A spin magnitude in `[0, 1]`: uniform, close to zero, close to one, and the three special values.
fn spin() -> impl Strategy<Value = f64> {
    prop_oneof![
        4 => 0.0..=1.0f64,
        2 => decades(-16.0, 0.0).prop_map(|deficit| 1.0 - deficit),
        2 => decades(-18.0, 0.0),
        1 => prop_oneof![Just(0.0), Just(1.0), Just(0.998)],
    ]
}

fn orbit() -> impl Strategy<Value = Orbit> {
    prop_oneof![Just(Orbit::Prograde), Just(Orbit::Retrograde)]
}

fn model() -> impl Strategy<Value = EmissionModel> {
    prop_oneof![
        Just(EmissionModel::PhotonBlackbody),
        Just(EmissionModel::PagePhotonsGravitons),
        (0.5..20.0f64).prop_map(|f| EmissionModel::MacGibbon { f }),
        decades(-7.0, -2.0).prop_map(EmissionModel::Coefficient),
    ]
}

proptest! {
    // ------------------------------------------------------------ Kerr: the ISCO

    /// The ISCO moves inwards as a prograde spin grows and outwards as a retrograde spin grows.
    #[test]
    fn isco_radius_is_monotonic_in_spin(a in spin(), b in spin()) {
        let (low, high) = if a <= b { (a, b) } else { (b, a) };
        let pro = |s| kerr::isco_radius(HOLE, s, Orbit::Prograde);
        let retro = |s| kerr::isco_radius(HOLE, s, Orbit::Retrograde);
        prop_assert!(le(pro(high), pro(low)), "prograde: {} at {high} > {} at {low}", pro(high), pro(low));
        prop_assert!(le(retro(low), retro(high)), "retrograde: {} at {low} > {} at {high}", retro(low), retro(high));
        // Strictly, once the spins differ by more than rounding can hide.
        if high - low > 1e-9 {
            prop_assert!(pro(high) < pro(low));
            prop_assert!(retro(high) > retro(low));
        }
    }

    /// 1 M ≤ prograde ISCO ≤ 6 M ≤ retrograde ISCO ≤ 9 M, and the characteristic radii are ordered.
    #[test]
    fn kerr_radii_are_bounded_and_ordered(a in spin(), sense in orbit(), negative in any::<bool>()) {
        let m = gravitational_radius(HOLE);
        let signed = if negative { -a } else { a };
        let isco = kerr::isco_radius(HOLE, signed, sense);
        match sense {
            Orbit::Prograde => prop_assert!(m <= isco && isco <= 6.0 * m, "prograde ISCO {}", isco / m),
            Orbit::Retrograde => prop_assert!(6.0 * m <= isco && isco <= 9.0 * m, "retrograde ISCO {}", isco / m),
        }
        let inner = kerr::inner_horizon_radius(HOLE, signed);
        let outer = kerr::outer_horizon_radius(HOLE, signed);
        let photon = kerr::photon_orbit_radius(HOLE, signed, sense);
        let bound = kerr::marginally_bound_radius(HOLE, signed, sense);
        prop_assert!(0.0 <= inner && inner <= m && m <= outer && outer <= 2.0 * m);
        prop_assert!(le(outer, photon) && le(photon, bound) && le(bound, isco),
            "r+ {} r_ph {} r_mb {} r_isco {}", outer / m, photon / m, bound / m, isco / m);
        // The horizons are the roots of r^2 - 2 m r + a^2 m^2.
        prop_assert!(close(inner + outer, 2.0 * m, SLACK));
        prop_assert!(close(inner * outer, a * a * m * m, 2.0 * SLACK));
    }

    /// The thin-disc efficiency stays between its Schwarzschild and extremal values, grows with a prograde
    /// spin, and is one minus the ISCO energy.
    #[test]
    fn accretion_efficiency_is_bounded_and_monotonic(a in spin(), b in spin()) {
        let (low, high) = if a <= b { (a, b) } else { (b, a) };
        let schwarzschild = 1.0 - (8f64 / 9.0).sqrt();
        for s in [low, high] {
            let pro = kerr::accretion_efficiency(s, Orbit::Prograde);
            let retro = kerr::accretion_efficiency(s, Orbit::Retrograde);
            prop_assert!(le(schwarzschild, pro) && le(pro, 1.0 - 1.0 / 3f64.sqrt()), "prograde {pro} at {s}");
            prop_assert!(le(1.0 - (25f64 / 27.0).sqrt(), retro) && le(retro, schwarzschild), "retrograde {retro} at {s}");
            for sense in [Orbit::Prograde, Orbit::Retrograde] {
                let energy = kerr::isco_specific_energy(s, sense);
                prop_assert!(energy > 0.57 && energy < 0.97);
                prop_assert!((kerr::accretion_efficiency(s, sense) - (1.0 - energy)).abs() < 3e-16);
            }
        }
        prop_assert!(le(kerr::accretion_efficiency(low, Orbit::Prograde), kerr::accretion_efficiency(high, Orbit::Prograde)));
        prop_assert!(le(kerr::accretion_efficiency(high, Orbit::Retrograde), kerr::accretion_efficiency(low, Orbit::Retrograde)));
    }

    /// Time on the ISCO always runs slow, more so the faster a prograde hole spins.
    #[test]
    fn isco_time_dilation_is_bounded_and_monotonic(a in spin(), b in spin()) {
        let (low, high) = if a <= b { (a, b) } else { (b, a) };
        let pro = |s| kerr::isco_dt_dtau(s, Orbit::Prograde);
        let retro = |s| kerr::isco_dt_dtau(s, Orbit::Retrograde);
        prop_assert!(le(2f64.sqrt(), pro(low)) && le(pro(low), pro(high)), "{} then {}", pro(low), pro(high));
        prop_assert!(le(retro(high), retro(low)) && le(retro(low), 2f64.sqrt()));
        prop_assert!(le(13.0 / (6.0 * 3f64.sqrt()), retro(high)));
        prop_assert!(pro(high) > 0.0 && (pro(high).is_finite() || high == 1.0));
    }

    /// On the ISCO the general circular-orbit formulas give the ISCO closed forms, and the ISCO has the
    /// least energy of all circular orbits.
    #[test]
    fn isco_is_the_least_bound_energy_of_circular_orbits(
        a in 0.0..0.99f64, sense in orbit(), factor in decades(-3.0, 3.0),
    ) {
        let isco = kerr::isco_radius(HOLE, a, sense);
        let photon = kerr::photon_orbit_radius(HOLE, a, sense);
        let e_isco = kerr::isco_specific_energy(a, sense);
        prop_assert!(close(kerr::circular_orbit_specific_energy(HOLE, a, sense, isco), e_isco, 1e-12));
        prop_assert!(close(kerr::circular_orbit_dt_dtau(HOLE, a, sense, isco), kerr::isco_dt_dtau(a, sense), 1e-12));
        prop_assert!(close(
            kerr::circular_orbit_specific_angular_momentum(HOLE, a, sense, isco),
            kerr::isco_specific_angular_momentum(HOLE, a, sense),
            1e-12
        ));
        // any other radius outside the photon orbit
        let r = photon + (isco - photon) * factor;
        let energy = kerr::circular_orbit_specific_energy(HOLE, a, sense, r);
        prop_assert!(energy >= e_isco * (1.0 - 1e-13), "E({}) = {energy} < E_isco = {e_isco}", r / isco);
        prop_assert!(kerr::circular_orbit_dt_dtau(HOLE, a, sense, r) >= 1.0);
        // A prograde orbit is slower than the Kepler rate sqrt(GM/r^3), a retrograde one faster.
        let kepler = (G_NEWTON * HOLE / r).sqrt() / r;
        let omega = kerr::circular_orbit_angular_velocity(HOLE, a, sense, r);
        match sense {
            Orbit::Prograde => prop_assert!(le(omega, kepler)),
            Orbit::Retrograde => prop_assert!(le(kepler, omega)),
        }
    }

    /// No circular orbit exists inside the photon orbit, for any spin up to the extremal limit, where the
    /// photon orbit and the horizon nearly coincide.
    #[test]
    fn no_circular_orbit_inside_the_photon_orbit(a in spin(), sense in orbit(), fraction in 0.0..1.0f64) {
        let photon = kerr::photon_orbit_radius(HOLE, a, sense);
        let r = photon * fraction;
        prop_assume!(r > 0.0 && r < photon);
        prop_assert!(kerr::circular_orbit_dt_dtau(HOLE, a, sense, r).is_nan());
        prop_assert!(kerr::circular_orbit_specific_energy(HOLE, a, sense, r).is_nan());
        prop_assert!(kerr::circular_orbit_specific_angular_momentum(HOLE, a, sense, r).is_nan());
        prop_assert!(kerr::circular_orbit_angular_velocity(HOLE, a, sense, r).is_nan());
        // Just outside it, the orbit exists and its clock is slow.
        let outside = kerr::circular_orbit_dt_dtau(HOLE, a, sense, photon * 1.001);
        prop_assert!(outside > 1.0, "dt/dtau = {outside} just outside the photon orbit at spin {a}");
    }

    // ------------------------------------------------------------ Kerr: the hole and its observers

    /// Christodoulou's mass formula: M² = M_irr² + (J c/(2 G M_irr))², with J = a* G M²/c.
    #[test]
    fn irreducible_mass_satisfies_the_mass_formula(a in spin(), mass in decades(20.0, 42.0)) {
        let irreducible = kerr::irreducible_mass(mass, a);
        prop_assert!(le(mass / 2f64.sqrt(), irreducible) && le(irreducible, mass));
        let j = a * G_NEWTON * mass * mass / C;
        let rotational = j * C / (2.0 * G_NEWTON * irreducible);
        prop_assert!(close(mass * mass, irreducible * irreducible + rotational * rotational, 4.0 * SLACK));
        let fraction = kerr::extractable_energy_fraction(a);
        prop_assert!((0.0..=1.0 - 1.0 / 2f64.sqrt() + 1e-16).contains(&fraction));
        prop_assert!(close(mass * (1.0 - fraction), irreducible, 2.0 * SLACK));
        // the horizon area is 16 pi (G M_irr/c^2)^2 = 4 pi (r+^2 + a^2 m^2)
        let m = gravitational_radius(mass);
        let outer = kerr::outer_horizon_radius(mass, a);
        let area = 4.0 * PI * (outer * outer + a * a * m * m);
        prop_assert!(close(area, 16.0 * PI * gravitational_radius(irreducible).powi(2), 4.0 * SLACK));
    }

    /// A zero-angular-momentum observer's clock runs between 0 and 1, never slower than a static clock at
    /// the same place, and it is dragged round no faster than the horizon turns.
    #[test]
    fn zamo_clock_and_frame_dragging_are_bounded(
        a in spin(), negative in any::<bool>(), k in decades(0.0, 6.0), theta in 0.0..=PI,
    ) {
        let signed = if negative { -a } else { a };
        let horizon = kerr::outer_horizon_radius(HOLE, signed);
        let r = horizon * k;
        let lapse = kerr::zamo_dtau_dt(HOLE, signed, r, theta);
        prop_assert!((0.0..=1.0).contains(&lapse), "lapse {lapse}");
        let omega = kerr::zamo_angular_velocity(HOLE, signed, r, theta);
        let omega_h = kerr::horizon_angular_velocity(HOLE, signed);
        prop_assert!(omega == 0.0 || omega.signum() == signed.signum(), "omega {omega} for spin {signed}");
        prop_assert!(le(omega.abs(), omega_h.abs()), "omega {omega} beyond the horizon rate {omega_h}");
        let limit = kerr::ergosphere_radius(HOLE, signed, theta);
        prop_assert!(le(horizon, limit) && le(limit, 2.0 * gravitational_radius(HOLE)));
        let fixed = kerr::static_dtau_dt(HOLE, signed, r, theta);
        if r >= limit {
            prop_assert!(le(fixed, lapse), "static {fixed} above ZAMO {lapse}");
        } else {
            prop_assert!(fixed.is_nan());
        }
    }

    /// Frame dragging weakens outwards in the equatorial plane, and the ZAMO clock speeds up.
    #[test]
    fn frame_dragging_falls_with_radius(a in spin(), k1 in decades(0.0, 5.0), k2 in decades(0.0, 5.0)) {
        let horizon = kerr::outer_horizon_radius(HOLE, a);
        let (near, far) = if k1 <= k2 { (horizon * k1, horizon * k2) } else { (horizon * k2, horizon * k1) };
        let half_pi = PI / 2.0;
        prop_assert!(le(
            kerr::zamo_angular_velocity(HOLE, a, far, half_pi),
            kerr::zamo_angular_velocity(HOLE, a, near, half_pi)
        ));
        // Close to the horizon the lapse is only as well determined as r - r+, so allow for that there.
        let lapse_near = kerr::zamo_dtau_dt(HOLE, a, near, half_pi);
        let lapse_far = kerr::zamo_dtau_dt(HOLE, a, far, half_pi);
        prop_assert!(lapse_near <= lapse_far * (1.0 + 1e-9), "{lapse_near} at {near} > {lapse_far} at {far}");
    }

    /// At zero spin every Kerr function is its Schwarzschild counterpart.
    #[test]
    fn kerr_at_zero_spin_is_schwarzschild(x in decades(0.31, 8.0), theta in 0.0..=PI, sense in orbit()) {
        let r = x * gravitational_radius(HOLE);
        let exact = schwarzschild::static_dtau_dt(HOLE, r);
        prop_assert!(close(kerr::static_dtau_dt(HOLE, 0.0, r, theta), exact, SLACK));
        prop_assert!(close(kerr::zamo_dtau_dt(HOLE, 0.0, r, theta), exact, SLACK));
        prop_assert_eq!(kerr::zamo_angular_velocity(HOLE, 0.0, r, theta), 0.0);
        if x > 3.0 {
            // both lose digits together as the photon sphere is approached
            let slack = 4.0 * SLACK * x / (x - 3.0);
            prop_assert!(close(
                kerr::circular_orbit_dt_dtau(HOLE, 0.0, sense, r),
                schwarzschild::circular_orbit_dt_dtau(HOLE, r),
                slack
            ));
            prop_assert!(close(
                kerr::circular_orbit_angular_velocity(HOLE, 0.0, sense, r),
                schwarzschild::circular_orbit_angular_velocity(HOLE, r),
                SLACK
            ));
        }
    }

    // ------------------------------------------------------------ Schwarzschild

    /// Redshifts compose: (1 + z_ab)(1 + z_bc) = 1 + z_ac; and a redshift is undone on the way back.
    #[test]
    fn redshift_composes(xa in decades(0.0001, 7.0), xb in decades(0.0001, 7.0), xc in decades(0.0001, 7.0)) {
        let rs = schwarzschild::schwarzschild_radius(HOLE);
        let (a, b, c) = (xa * rs, xb * rs, xc * rs);
        let z = |from, to| schwarzschild::gravitational_redshift(HOLE, from, to);
        prop_assert!(close((1.0 + z(a, b)) * (1.0 + z(b, c)), 1.0 + z(a, c), 1e-13));
        prop_assert!(close((1.0 + z(a, b)) * (1.0 + z(b, a)), 1.0, 1e-13));
        prop_assert!(z(a, b) > -1.0);
        prop_assert_eq!(z(a, b) > 0.0, b > a);
        // against the static clock rates
        let ratio = schwarzschild::static_dtau_dt(HOLE, b) / schwarzschild::static_dtau_dt(HOLE, a);
        prop_assert!(close(1.0 + z(a, b), ratio, 1e-13));
    }

    /// Clocks run faster further out: the static rate rises towards 1 and the orbital dilation falls towards 1.
    #[test]
    fn schwarzschild_clock_rates_are_monotonic(x1 in decades(0.0, 9.0), x2 in decades(0.0, 9.0)) {
        let rs = schwarzschild::schwarzschild_radius(HOLE);
        let (near, far) = if x1 <= x2 { (x1 * rs, x2 * rs) } else { (x2 * rs, x1 * rs) };
        let (slow, fast) = (schwarzschild::static_dtau_dt(HOLE, near), schwarzschild::static_dtau_dt(HOLE, far));
        prop_assert!((0.0..=1.0).contains(&slow) && slow <= fast * (1.0 + 1e-9) && fast <= 1.0);
        let photon = schwarzschild::photon_sphere_radius(HOLE);
        let (inner, outer) = (near.max(photon) * 1.000_001, far.max(photon) * 1.000_002);
        let (deep, shallow) = (
            schwarzschild::circular_orbit_dt_dtau(HOLE, inner),
            schwarzschild::circular_orbit_dt_dtau(HOLE, outer),
        );
        prop_assert!(deep >= shallow && shallow >= 1.0, "{deep} at {inner}, {shallow} at {outer}");
    }

    // ------------------------------------------------------------ weak field

    /// The two-leg and two-point forms of the Shapiro delay agree; the delay grows with the length of the
    /// path, falls as the path moves away from the mass, and depends only on the shape of the triangle.
    #[test]
    fn shapiro_delay_forms_agree(b in decades(8.0, 11.0), ke in decades(0.001, 3.0), kr in decades(0.001, 3.0)) {
        let (r_e, r_r) = (b * ke, b * kr);
        let legs = weak_field::shapiro_delay(M_SUN, r_e, r_r, b);
        let separation = (r_e * r_e - b * b).sqrt() + (r_r * r_r - b * b).sqrt();
        let two_point = weak_field::shapiro_delay_between(M_SUN, r_e, r_r, separation);
        // The two-point form subtracts nearly equal lengths for long, grazing paths: its argument is within
        // 1/(2 k^2) of 1, so rounding is magnified by k^2.
        let k = ke.max(kr);
        prop_assert!(close(legs, two_point, 1e-13 + 1e-15 * k * k), "{legs} against {two_point}");
        prop_assert!(legs > 0.0);
        prop_assert!(weak_field::shapiro_delay(M_SUN, 2.0 * r_e, r_r, b) > legs);
        prop_assert!(weak_field::shapiro_delay(M_SUN, r_e, r_r, 0.5 * b) > legs);
        prop_assert!(close(weak_field::shapiro_delay(M_SUN, 2.0 * r_e, 2.0 * r_r, 2.0 * b), legs, SLACK));
        // the long-path logarithm is an upper bound
        let log_form = 2.0 * G_NEWTON * M_SUN / (C * C * C) * (4.0 * r_e * r_r / (b * b)).ln();
        prop_assert!(le(legs, log_form));
    }

    /// The Einstein ring is where the deflection closes the lens equation.
    #[test]
    fn einstein_angle_solves_the_lens_equation(
        mass in decades(28.0, 42.0), d_lens in decades(16.0, 24.0), ratio in 1.01..50.0f64,
    ) {
        let d_source = d_lens * ratio;
        let d_lens_source = d_source - d_lens;
        let theta = weak_field::einstein_angle(mass, d_lens, d_source, d_lens_source);
        let bend = weak_field::deflection_angle(mass, theta * d_lens);
        prop_assert!(close(bend * d_lens_source / d_source, theta, 4.0 * SLACK));
    }

    /// An orbiting clock gains on a surface clock exactly when it is above 3R/2 (non-rotating body).
    #[test]
    fn orbiting_clock_changes_sign_at_three_halves(radius in decades(5.0, 8.0), k in 1.0..20.0f64) {
        let earth = 5.9722e24;
        let offset = weak_field::orbit_vs_surface_rate_offset(earth, k * radius, radius, 0.0);
        if k > 1.5 + 1e-9 {
            prop_assert!(offset > 0.0);
        } else if k < 1.5 - 1e-9 {
            prop_assert!(offset < 0.0);
        }
        // spinning the surface clock makes it slower, so the orbiting clock gains more
        prop_assert!(weak_field::orbit_vs_surface_rate_offset(earth, k * radius, radius, 400.0) > offset);
    }

    // ------------------------------------------------------------ Hawking

    /// T ∝ 1/M, P ∝ 1/M², t ∝ M³, and the mass for a lifetime inverts the lifetime.
    #[test]
    fn hawking_scalings(mass in decades(-3.0, 45.0), factor in decades(-2.0, 2.0), emission in model()) {
        let other = factor * mass;
        prop_assert!(close(hawking::temperature(mass) * mass, hawking::temperature(other) * other, 2.0 * SLACK));
        prop_assert!(close(
            hawking::power(mass, emission) * mass * mass,
            hawking::power(other, emission) * other * other,
            4.0 * SLACK
        ));
        prop_assert!(close(
            hawking::lifetime(other, emission) * mass.powi(3),
            hawking::lifetime(mass, emission) * other.powi(3),
            4.0 * SLACK
        ));
        let lifetime = hawking::lifetime(mass, emission);
        prop_assert!(close(hawking::mass_for_lifetime(lifetime, emission), mass, 2.0 * SLACK));
        // the energy radiated over the lifetime is the rest energy: mean power = 3 x initial power
        prop_assert!(close(3.0 * hawking::power(mass, emission) * lifetime, mass * C * C, 4.0 * SLACK));
        // entropy grows as the square of the mass
        prop_assert!(close(hawking::entropy(other), factor * factor * hawking::entropy(mass), 4.0 * SLACK));
    }

    // ------------------------------------------------------------ Newtonian

    /// A Hohmann transfer takes between half the inner and half the outer period, costs nothing between
    /// equal orbits, and reverses exactly.
    #[test]
    fn hohmann_transfer_is_consistent(mass in decades(20.0, 32.0), r1 in decades(6.0, 13.0), ratio in 1.0..100.0f64) {
        let r2 = r1 * ratio;
        let time = newtonian::hohmann_transfer_time(mass, r1, r2);
        prop_assert!(le(0.5 * newtonian::circular_period(mass, r1), time));
        prop_assert!(le(time, 0.5 * newtonian::circular_period(mass, r2)));
        let out = newtonian::hohmann_delta_v(mass, r1, r2);
        let back = newtonian::hohmann_delta_v(mass, r2, r1);
        prop_assert!(out.departure >= 0.0 && out.arrival >= 0.0 && back.departure <= 0.0 && back.arrival <= 0.0);
        prop_assert!(close(out.departure, -back.arrival, SLACK) && close(out.arrival, -back.departure, SLACK));
        // neither burn exceeds the gap between circular and escape speed at its end
        let v1 = newtonian::circular_speed(mass, r1);
        prop_assert!(le(out.departure, (2f64.sqrt() - 1.0) * v1));
        prop_assert!(le(out.arrival, newtonian::circular_speed(mass, r2)));
        prop_assert!(close(newtonian::escape_speed(mass, r1), 2f64.sqrt() * v1, SLACK));
    }

    /// The synodic period is symmetric in the two bodies and never shorter than the shorter period.
    #[test]
    fn synodic_period_exceeds_the_shorter_period(a in decades(0.0, 9.0), b in decades(0.0, 9.0)) {
        let synodic = newtonian::synodic_period(a, b);
        prop_assert!(close(synodic, newtonian::synodic_period(b, a), SLACK) || synodic.is_infinite());
        prop_assert!(le(a.min(b), synodic));
        if a != b {
            // the defining difference of rates, which cancels for nearly equal periods
            let tolerance = 1e-15 + 4e-16 * (a + b) / (a - b).abs();
            prop_assert!(close(1.0 / synodic, (1.0 / a - 1.0 / b).abs(), tolerance));
        }
    }

    // ------------------------------------------------------------ totality

    /// Every function accepts any binary64 arguments without panicking, and a spin outside [-1, 1] gives NaN.
    #[test]
    fn functions_are_total(
        mass in any::<f64>(), a in any::<f64>(), r in any::<f64>(), theta in any::<f64>(), other in any::<f64>(),
        sense in orbit(),
    ) {
        let valid_spin = (-1.0..=1.0).contains(&a);
        let kerr_results = [
            kerr::outer_horizon_radius(mass, a),
            kerr::inner_horizon_radius(mass, a),
            kerr::ergosphere_radius(mass, a, theta),
            kerr::horizon_angular_velocity(mass, a),
            kerr::irreducible_mass(mass, a),
            kerr::extractable_energy_fraction(a),
            kerr::static_dtau_dt(mass, a, r, theta),
            kerr::zamo_angular_velocity(mass, a, r, theta),
            kerr::zamo_dtau_dt(mass, a, r, theta),
            kerr::photon_orbit_radius(mass, a, sense),
            kerr::marginally_bound_radius(mass, a, sense),
            kerr::isco_radius(mass, a, sense),
            kerr::isco_specific_energy(a, sense),
            kerr::isco_specific_angular_momentum(mass, a, sense),
            kerr::isco_dt_dtau(a, sense),
            kerr::accretion_efficiency(a, sense),
            kerr::circular_orbit_dt_dtau(mass, a, sense, r),
            kerr::circular_orbit_specific_energy(mass, a, sense, r),
            kerr::circular_orbit_specific_angular_momentum(mass, a, sense, r),
            kerr::circular_orbit_angular_velocity(mass, a, sense, r),
        ];
        if !valid_spin {
            prop_assert!(kerr_results.iter().all(|v| v.is_nan()));
        }
        if mass.is_nan() || mass < 0.0 {
            // every function that takes a mass refuses a negative or NaN mass
            for (index, value) in kerr_results.iter().enumerate() {
                let takes_mass = ![5, 12, 14, 15].contains(&index);
                prop_assert!(!takes_mass || value.is_nan(), "Kerr function {index} accepted mass {mass}");
            }
            prop_assert!(schwarzschild::schwarzschild_radius(mass).is_nan());
            prop_assert!(schwarzschild::static_dtau_dt(mass, r).is_nan());
            prop_assert!(schwarzschild::gravitational_redshift(mass, r, other).is_nan());
            prop_assert!(schwarzschild::circular_orbit_dt_dtau(mass, r).is_nan());
            prop_assert!(weak_field::shapiro_delay(mass, r, other, theta).is_nan());
            prop_assert!(weak_field::shapiro_delay_between(mass, r, other, theta).is_nan());
            prop_assert!(weak_field::deflection_angle(mass, r).is_nan());
            prop_assert!(weak_field::einstein_angle(mass, r, other, theta).is_nan());
            prop_assert!(weak_field::orbit_vs_surface_rate_offset(mass, r, other, theta).is_nan());
            prop_assert!(hawking::temperature(mass).is_nan());
            prop_assert!(hawking::power(mass, EmissionModel::PhotonBlackbody).is_nan());
            prop_assert!(hawking::lifetime(mass, EmissionModel::Coefficient(other)).is_nan());
            prop_assert!(newtonian::binding_energy_uniform_sphere(mass, r).is_nan());
            prop_assert!(newtonian::hohmann_delta_v(mass, r, other).total().is_nan());
        }
        // The rest must simply return: any value, no panic.
        let _ = (
            schwarzschild::static_dtau_dt(mass, r),
            schwarzschild::gravitational_redshift(mass, r, other),
            schwarzschild::circular_orbit_dt_dtau(mass, r),
            schwarzschild::circular_orbit_angular_velocity(mass, r),
            weak_field::shapiro_delay_leg(mass, r, other),
            weak_field::shapiro_delay_between(mass, r, other, theta),
            weak_field::static_rate_offset(mass, r),
            weak_field::circular_orbit_rate_offset(mass, r),
            hawking::entropy(mass),
            hawking::mass_loss_rate(mass, EmissionModel::MacGibbon { f: other }),
            hawking::mass_for_lifetime(r, EmissionModel::Coefficient(other)),
            newtonian::circular_speed(mass, r),
            newtonian::escape_speed(mass, r),
            newtonian::circular_period(mass, r),
            newtonian::synodic_period(r, other),
            newtonian::hohmann_transfer_time(mass, r, other),
        );
        // A static or ZAMO clock rate, where it exists, is a number between 0 and 1.
        for rate in [schwarzschild::static_dtau_dt(mass, r), kerr_results[6], kerr_results[8]] {
            prop_assert!(rate.is_nan() || (0.0..=1.0).contains(&rate), "clock rate {rate}");
        }
    }
}
