//! Property tests: identities and limits that must hold for any flat matter + Λ cosmology and any epoch.
//!
//! Each property runs 512 cases from proptest's deterministic generator, so every run tests the same inputs
//! and a failure reproduces exactly.

use proptest::prelude::*;
use proptest::test_runner::{RngAlgorithm, TestCaseError, TestRng, TestRunner};
use srkit_flrw::Cosmology;

const EPS: f64 = f64::EPSILON;

fn check<S: Strategy>(strategy: S, test: impl Fn(S::Value) -> Result<(), TestCaseError>)
where
    S::Value: core::fmt::Debug,
{
    let config = ProptestConfig {
        cases: 512,
        failure_persistence: None,
        ..ProptestConfig::default()
    };
    let rng = TestRng::deterministic_rng(RngAlgorithm::ChaCha);
    let mut runner = TestRunner::new_with_rng(config, rng);
    if let Err(error) = runner.run(&strategy, test) {
        panic!("{error}");
    }
}

/// Any cosmology with `0.02 < Ω_m < 0.98` and `H0` between 1e-12 and 1e-8 per year.
fn cosmology() -> impl Strategy<Value = Cosmology> {
    (0.02..0.98f64, -12.0..-8.0f64).prop_map(|(omega_m, log_h0)| {
        Cosmology::per_year(10f64.powf(log_h0), omega_m).expect("valid cosmology")
    })
}

/// A cosmology and a cosmic time with `H_Λ t` log-uniform between `10^lo` and `10^hi`.
fn epoch(lo: f64, hi: f64) -> impl Strategy<Value = (Cosmology, f64)> {
    (cosmology(), lo..hi).prop_map(|(u, log_ht)| (u, 10f64.powf(log_ht) / u.h_lambda()))
}

/// Three cosmic times in increasing order, each with `H_Λ t` log-uniform between `10^lo` and `10^hi`.
fn three_epochs(lo: f64, hi: f64) -> impl Strategy<Value = (Cosmology, [f64; 3])> {
    (cosmology(), lo..hi, lo..hi, lo..hi).prop_map(|(u, x, y, z)| {
        let mut times = [x, y, z].map(|log_ht| 10f64.powf(log_ht) / u.h_lambda());
        times.sort_by(f64::total_cmp);
        (u, times)
    })
}

/// `|a − b| ≤ tol · max(|a|, |b|)`.
fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs())
}

/// A step for central differences: small against both the age and the Hubble time, so the truncation error
/// is near 1e-8 and the rounding error stays far smaller.
fn small_step(u: &Cosmology, t: f64) -> f64 {
    1e-4 * t.min(1.0 / u.h_lambda())
}

#[test]
fn time_and_scale_factor_round_trip() {
    check(epoch(-8.0, 2.5), |(u, t)| {
        let a = u.scale_factor(t);
        prop_assert!(a > 0.0 && a.is_finite());
        prop_assert!(close(u.time_at_scale_factor(a), t, 1e-13), "t(a(t))");
        // and through the logarithm, which also works after a overflows
        let ln_a = u.ln_scale_factor(t);
        prop_assert!(
            close(u.time_at_ln_scale_factor(ln_a), t, 1e-13),
            "t(ln a(t))"
        );
        prop_assert!((ln_a - a.ln()).abs() <= 1e-13 * (1.0 + ln_a.abs()));
        // a(t(a)) = a: a relative error in t is amplified by d ln a / d ln t ≤ max(1, H_Λ t)
        let back = u.scale_factor(u.time_at_scale_factor(a));
        prop_assert!(close(back, a, 1e-13 * (1.0 + u.h_lambda() * t)), "a(t(a))");
        Ok(())
    });
}

#[test]
fn the_logarithm_round_trips_far_beyond_overflow() {
    check(epoch(2.5, 12.0), |(u, t)| {
        let ln_a = u.ln_scale_factor(t);
        prop_assert!(ln_a.is_finite());
        prop_assert!(close(u.time_at_ln_scale_factor(ln_a), t, 1e-13));
        Ok(())
    });
}

#[test]
fn the_friedmann_equation_holds() {
    check(epoch(-8.0, 2.0), |(u, t)| {
        // H(t) and H(a(t)) are two closed forms for one number
        let a = u.scale_factor(t);
        let tol = 1e-13 * (1.0 + u.h_lambda() * t);
        prop_assert!(close(u.hubble(t), u.hubble_at_scale_factor(a), tol));
        prop_assert!(u.hubble(t) >= u.h_lambda());
        // H = d ln a / dt, by a central difference
        let h = small_step(&u, t);
        let slope = (u.ln_scale_factor(t + h) - u.ln_scale_factor(t - h)) / (2.0 * h);
        prop_assert!(close(slope, u.hubble(t), 1e-6));
        Ok(())
    });
}

#[test]
fn today_is_scale_factor_one() {
    check(cosmology(), |u| {
        let now = u.age();
        prop_assert!(close(u.scale_factor(now), 1.0, 8.0 * EPS));
        prop_assert!(close(u.hubble(now), u.h0(), 8.0 * EPS));
        prop_assert!(u.ln_scale_factor(now).abs() <= 8.0 * EPS);
        prop_assert!(close(u.time_at_scale_factor(1.0), now, 8.0 * EPS));
        prop_assert!(close(u.time_at_redshift(0.0), now, 8.0 * EPS));
        prop_assert!(u.lookback_time(0.0).abs() <= 8.0 * EPS * now);
        prop_assert!(
            u.comoving_distance_at_redshift(0.0).abs() <= 8.0 * EPS * u.conformal_time(now)
        );
        // matter and Λ are equally dense at a_eq
        let t_eq = u.time_at_scale_factor(u.scale_factor_at_equality());
        prop_assert!(close(
            u.hubble(t_eq),
            u.h_lambda() * core::f64::consts::SQRT_2,
            8.0 * EPS
        ));
        Ok(())
    });
}

#[test]
fn conformal_time_is_monotone_and_bounded() {
    check((epoch(-8.0, 2.5), 1e-6..10.0f64), |((u, t1), step)| {
        let t2 = t1 * (1.0 + step);
        let limit = u.conformal_time_limit();
        let (e1, e2) = (u.conformal_time(t1), u.conformal_time(t2));
        let (r1, r2) = (
            u.conformal_time_remaining(t1),
            u.conformal_time_remaining(t2),
        );
        // monotone, to within the rounding of one evaluation
        prop_assert!(e2 - e1 >= -8.0 * EPS * limit, "η fell: {e1} -> {e2}");
        prop_assert!(r1 - r2 >= -8.0 * EPS * r1, "η_∞ − η rose: {r1} -> {r2}");
        prop_assert!(u.scale_factor(t2) >= u.scale_factor(t1));
        prop_assert!(u.hubble(t2) <= u.hubble(t1));
        // bounded: 0 < η < η_∞ and the two parts sum to the whole
        prop_assert!(e1 > 0.0 && e1 <= limit);
        prop_assert!(r1 >= 0.0 && r1 < limit);
        prop_assert!(close(e1 + r1, limit, 8.0 * EPS));
        // the difference of the two is the comoving distance between the epochs
        let d = u.comoving_distance(t1, t2);
        prop_assert!(d >= 0.0);
        prop_assert!((d - (e2 - e1)).abs() <= 16.0 * EPS * limit);
        prop_assert!((d + u.comoving_distance(t2, t1)).abs() <= 16.0 * EPS * limit);
        Ok(())
    });
}

#[test]
fn conformal_time_inverts() {
    check(epoch(-8.0, 2.5), |(u, t)| {
        let eta = u.conformal_time(t);
        let rem = u.conformal_time_remaining(t);
        let a = u.scale_factor(t);
        // each inverse is as good as its conditioning: d ln t / d ln η = η a / t, and the same for η_∞ − η
        if eta < 0.999 * u.conformal_time_limit() {
            let back = u.time_at_conformal_time(eta).unwrap();
            prop_assert!(close(
                back,
                t,
                1e-13 * (1.0 + eta * a / t) + 1e-13 * eta / rem
            ));
        }
        let back = u.time_at_conformal_time_remaining(rem).unwrap();
        prop_assert!(close(back, t, 1e-13 * (1.0 + rem * a / t)));
        // d η / dt = 1 / a, by a central difference
        let h = small_step(&u, t);
        let slope = u.comoving_distance(t - h, t + h) / (2.0 * h);
        prop_assert!(close(slope, 1.0 / a, 1e-6));
        Ok(())
    });
}

#[test]
fn conformal_time_agrees_between_time_and_scale_factor() {
    check(epoch(-8.0, 1.5), |(u, t)| {
        let a = u.scale_factor(t);
        let tol = 1e-13 * (1.0 + u.h_lambda() * t);
        prop_assert!(close(
            u.conformal_time_at_scale_factor(a),
            u.conformal_time(t),
            tol
        ));
        prop_assert!(close(
            u.particle_horizon_proper(t),
            a * u.conformal_time(t),
            tol
        ));
        prop_assert!(close(
            u.event_horizon_proper(t),
            a * u.conformal_time_remaining(t),
            tol
        ));
        Ok(())
    });
}

#[test]
fn the_early_universe_is_matter_dominated() {
    // for θ = (3/2) H_Λ t ≪ 1: a = a_eq θ^(2/3) (1 + θ²/9 ...), H = 2/(3t), η = 3t/a
    check(epoch(-8.0, -3.0), |(u, t)| {
        let theta = 1.5 * u.h_lambda() * t;
        let tol = theta * theta + 1e-13;
        let a = u.scale_factor(t);
        prop_assert!(close(
            a,
            u.scale_factor_at_equality() * theta.powf(2.0 / 3.0),
            tol
        ));
        prop_assert!(close(u.hubble(t), 2.0 / (3.0 * t), tol));
        prop_assert!(close(u.conformal_time(t), 3.0 * t / a, tol));
        prop_assert!(close(u.scale_factor_integral(t), 0.6 * a * t, tol));
        Ok(())
    });
}

#[test]
fn the_late_universe_is_de_sitter() {
    // for H_Λ t ≥ 15 every matter correction is below e^(−45)
    check((epoch(1.2, 2.4), 0.01..40.0f64), |((u, t), h_dt)| {
        let h = u.h_lambda();
        prop_assert!(close(u.hubble(t), h, 2.0 * EPS));
        // exponential expansion
        let dt = h_dt / h;
        let growth = u.ln_scale_factor(t + dt) - u.ln_scale_factor(t);
        prop_assert!((growth - h_dt).abs() <= 1e-12 * (1.0 + h * (t + dt)));
        // the event horizon sits at the de Sitter radius, and the conformal time left is 1/(a H_Λ)
        prop_assert!(close(
            u.event_horizon_proper(t),
            u.de_sitter_horizon(),
            4.0 * EPS
        ));
        let ln_product = u.conformal_time_remaining(t).ln() + u.ln_scale_factor(t) + h.ln();
        prop_assert!(
            ln_product.abs() <= 1e-12 * (1.0 + h * t),
            "ln(a H_Λ (η_∞ − η)) = {ln_product}"
        );
        // η itself is within e^(−H_Λ t) of its limit (and within rounding of it once that is below 1e-16)
        let limit = u.conformal_time_limit();
        let short = limit - u.conformal_time(t);
        prop_assert!((0.0..=((-h * t).exp() + 4.0 * EPS) * limit).contains(&short));
        // light obeys the de Sitter light cone: a target at a fraction f of the horizon is reached after
        // −ln(1 − f)/H_Λ
        let f = (h_dt / 40.0).clamp(0.1, 0.99);
        let chi = f * u.conformal_time_remaining(t);
        let arrival = u.light_arrival_time(t, chi).unwrap();
        prop_assert!(close(arrival - t, -(-f).ln_1p() / h, 1e-9));
        // Hubble drag: thrust in full ΛCDM equals the de Sitter closed form
        let alpha = 3.0 * h;
        let full = u.thrust_proper_velocity(0.5, alpha, t, t + dt);
        let limit = u.de_sitter_thrust_proper_velocity(0.5, alpha, dt);
        prop_assert!(close(full, limit, 1e-9), "{full} against {limit}");
        Ok(())
    });
}

#[test]
fn light_cones_are_causal_and_invertible() {
    check((epoch(-3.0, 1.5), 1e-6..0.99f64), |((u, t), f)| {
        let left = u.conformal_time_remaining(t);
        let chi = f * left;
        let arrival = u
            .light_arrival_time(t, chi)
            .expect("inside the event horizon");
        prop_assert!(arrival > t);
        // expansion only ever lengthens the trip: a_emitted χ ≤ Δt ≤ a_received χ
        let (a_e, a_r) = (u.scale_factor(t), u.scale_factor(arrival));
        let dt = arrival - t;
        let slack = 1e-9 * dt + 64.0 * EPS * arrival;
        prop_assert!(dt >= a_e * chi - slack, "{dt} < {}", a_e * chi);
        prop_assert!(dt <= a_r * chi + slack, "{dt} > {}", a_r * chi);
        // the past light cone undoes the future one
        let emission = u
            .light_emission_time(arrival, chi)
            .expect("inside the particle horizon");
        prop_assert!(close(emission, t, 1e-10 * (1.0 + left * a_e / t)));
        // the redshift along the way is the stretch of the scale factor
        prop_assert!(close(
            1.0 + u.redshift_between(t, arrival),
            a_r / a_e,
            1e-12
        ));
        // beyond the event horizon nothing arrives
        prop_assert!(u.light_arrival_time(t, left * (1.0 + 1e-9)).is_none());
        prop_assert!(u
            .light_arrival_time(t, 2.0 * u.conformal_time_limit())
            .is_none());
        // and beyond the particle horizon nothing has been received
        prop_assert!(u
            .light_emission_time(t, u.conformal_time(t) * (1.0 + 1e-9))
            .is_none());
        Ok(())
    });
}

#[test]
fn the_conformal_chart_is_minkowski() {
    check((epoch(-3.0, 0.5), 1e-3..0.9f64), |((u, t), f)| {
        let eta = u.conformal_time(t);
        let chi = f * u.conformal_time_remaining(t);
        let reception = u
            .light_reception(eta, chi)
            .expect("inside the event horizon");
        // the cone is a straight line: η_received = η_emitted + χ, exactly as computed
        prop_assert_eq!(reception.conformal_time, eta + chi);
        // and the cosmic time is the one the remaining-time form gives
        let arrival = u.light_arrival_time(t, chi).unwrap();
        let rem = u.conformal_time_limit() - reception.conformal_time;
        prop_assert!(close(
            reception.time,
            arrival,
            1e-11 * u.conformal_time_limit() / rem
        ));
        // past η_∞ the signal never arrives
        prop_assert!(u.light_reception(eta, u.conformal_time_limit()).is_none());
        prop_assert!(u.light_reception(eta, -chi).is_none());
        Ok(())
    });
}

#[test]
fn last_contact_is_the_edge_of_the_event_horizon() {
    check((cosmology(), -6.0..-0.001f64), |(u, log_fraction)| {
        // any target that was ever inside the horizon
        let chi = 10f64.powf(log_fraction) * u.conformal_time_limit();
        let last = u.last_contact_time(chi).expect("was once in contact");
        prop_assert!(last > 0.0 && last.is_finite());
        // the event horizon at that moment is exactly the target's distance
        let a = u.scale_factor(last);
        let tol = 1e-13 * (1.0 + last / (a * chi));
        prop_assert!(close(u.event_horizon(last), chi, tol));
        // a signal sent a little earlier arrives; one sent a little later never does
        let margin = 1e-6;
        prop_assert!(u.light_arrival_time(last * (1.0 - margin), chi).is_some());
        prop_assert!(u.light_arrival_time(last * (1.0 + margin), chi).is_none());
        // an answer needs twice the conformal time, so the last moment to ask comes earlier
        if let Some(reply) = u.last_reply_time(chi) {
            prop_assert!(reply < last);
            let there = u.light_arrival_time(reply * (1.0 - margin), chi).unwrap();
            prop_assert!(u.light_arrival_time(there, chi).is_some());
        } else {
            prop_assert!(2.0 * chi > u.conformal_time_limit());
        }
        // a nearer target stays in contact longer
        prop_assert!(u.last_contact_time(0.5 * chi).unwrap() > last);
        Ok(())
    });
}

#[test]
fn reach_shrinks_and_the_horizons_cross_at_half_equality() {
    check((epoch(-4.0, 1.5), 1e-3..3.0f64), |((u, t1), step)| {
        let t2 = t1 * (1.0 + step);
        prop_assert!(u.reach(t2) < u.reach(t1));
        prop_assert!(u.reach_mass(t2) < u.reach_mass(t1));
        prop_assert!(u.reachable_share_of_observable(t2) < u.reachable_share_of_observable(t1));
        prop_assert!(u.reach_volume_ratio(t2, t1) < 1.0);
        prop_assert!(close(u.reach_volume_ratio(t1, t1), 1.0, 4.0 * EPS));
        let ratio = u.reach(t2) / u.reach(t1);
        prop_assert!(close(
            u.reach_volume_ratio(t2, t1),
            ratio * ratio * ratio,
            1e-13
        ));
        prop_assert!(close(
            u.reach_mass(t2) / u.reach_mass(t1),
            ratio * ratio * ratio,
            1e-13
        ));
        // the reach mass depends on H_Λ and H_Λ t alone: at equal H_Λ another Ω_m gives the same mass
        let other = Cosmology::from_h_lambda_per_year(u.h_lambda(), 0.5).unwrap();
        prop_assert!(close(other.reach_mass(t1), u.reach_mass(t1), 1e-12));
        // the particle horizon equals the event horizon when a = a_eq / 2, for every cosmology
        let crossing = u.time_at_scale_factor(0.5 * u.scale_factor_at_equality());
        prop_assert!(close(u.reachable_share_of_observable(crossing), 1.0, 1e-13));
        Ok(())
    });
}

#[test]
fn peculiar_momentum_decays_as_one_over_a() {
    check(
        (three_epochs(-6.0, 2.3), -3.0..9.0f64),
        |((u, [t1, t2, t3]), log_u)| {
            let p = 10f64.powf(log_u);
            let (a1, a2) = (u.scale_factor(t1), u.scale_factor(t2));
            let tol = 1e-12 * (1.0 + u.h_lambda() * t3);
            // p a is conserved
            let p2 = u.coasting_proper_velocity(p, t1, t2);
            prop_assert!(close(p2 * a2, p * a1, tol));
            prop_assert!(p2 <= p);
            // coasting in two legs is coasting in one
            let p3 = u.coasting_proper_velocity(p2, t2, t3);
            prop_assert!(close(p3, u.coasting_proper_velocity(p, t1, t3), tol));
            // rapidity and Lorentz factor say the same thing
            let phi2 = u.coasting_rapidity(p.asinh(), t1, t2);
            prop_assert!(close(phi2.sinh(), p2, tol));
            let gamma2 = u.coasting_gamma(p.hypot(1.0), t1, t2);
            prop_assert!(gamma2 >= 1.0);
            prop_assert!(close(gamma2, p2.hypot(1.0), tol));
            // direction is kept
            prop_assert_eq!(u.coasting_rapidity(-p.asinh(), t1, t2), -phi2);
            prop_assert_eq!(u.coasting_proper_velocity(-p, t1, t2), -p2);
            // coasting is thrust with the engine off
            prop_assert!(close(u.thrust_proper_velocity(p, 0.0, t1, t2), p2, tol));
            Ok(())
        },
    );
}

#[test]
fn thrust_obeys_its_equation_of_motion() {
    let strategy = (
        epoch(-2.0, 1.5),
        0.05..5.0f64,
        0.05..5.0f64,
        -2.0..2.0f64,
        0.0..5.0f64,
    );
    check(strategy, |((u, t0), step1, step2, log_alpha, u0)| {
        let t1 = t0 * (1.0 + step1);
        let t2 = t1 * (1.0 + step2);
        let alpha = 10f64.powf(log_alpha) * u.h_lambda();
        // burning in two legs is burning in one
        let u1 = u.thrust_proper_velocity(u0, alpha, t0, t1);
        let u2 = u.thrust_proper_velocity(u1, alpha, t1, t2);
        prop_assert!(close(
            u2,
            u.thrust_proper_velocity(u0, alpha, t0, t2),
            1e-11
        ));
        // d(∫ a dt)/dt = a, by a central difference
        let h = small_step(&u, t1);
        let slope = (u.scale_factor_integral(t1 + h) - u.scale_factor_integral(t1 - h)) / (2.0 * h);
        prop_assert!(close(slope, u.scale_factor(t1), 1e-6));
        // du/dt = α − H u, by a central difference around t1
        let ahead = u.thrust_proper_velocity(u1, alpha, t1, t1 + h);
        let behind = u.thrust_proper_velocity(u1, alpha, t1, t1 - h);
        let rate = (ahead - behind) / (2.0 * h);
        let expected = alpha - u.hubble(t1) * u1;
        prop_assert!(
            (rate - expected).abs() <= 1e-5 * (alpha + u.hubble(t1) * u1.abs()),
            "du/dt = {rate}, α − H u = {expected}"
        );
        // drag in ΛCDM is at least the de Sitter drag, so a rocket from rest is never faster than in de
        // Sitter space, and never faster than the saturation value α/H_Λ
        let from_rest = u.thrust_proper_velocity(0.0, alpha, t0, t2);
        let de_sitter = u.de_sitter_thrust_proper_velocity(0.0, alpha, t2 - t0);
        prop_assert!(from_rest > 0.0);
        prop_assert!(from_rest <= de_sitter * (1.0 + 1e-9));
        prop_assert!(de_sitter <= u.drag_limited_proper_velocity(alpha));
        Ok(())
    });
}

#[test]
fn de_sitter_proper_time_is_additive_and_bounded() {
    let strategy = (
        cosmology(),
        -8.0..2.5f64,
        -8.0..2.5f64,
        -4.0..11.0f64,
        0.0..1e3f64,
        prop::bool::ANY,
    );
    check(strategy, |(u, log_h1, log_h2, log_alpha, u0, brake)| {
        let h = u.h_lambda();
        let (dt1, dt2) = (10f64.powf(log_h1) / h, 10f64.powf(log_h2) / h);
        let alpha = 10f64.powf(log_alpha) * h * if brake { -1.0 } else { 1.0 };
        let tau1 = u.de_sitter_thrust_proper_time(u0, alpha, dt1);
        let u1 = u.de_sitter_thrust_proper_velocity(u0, alpha, dt1);
        let tau2 = u.de_sitter_thrust_proper_time(u1, alpha, dt2);
        let total = u.de_sitter_thrust_proper_time(u0, alpha, dt1 + dt2);
        // two legs make one burn
        prop_assert!(
            close(tau1 + tau2, total, 1e-11),
            "{tau1} + {tau2} against {total}"
        );
        let u2 = u.de_sitter_thrust_proper_velocity(u1, alpha, dt2);
        let direct = u.de_sitter_thrust_proper_velocity(u0, alpha, dt1 + dt2);
        prop_assert!((u2 - direct).abs() <= 1e-11 * (u2.abs() + u0 + alpha.abs() / h));
        // the clock aboard runs slow by the Lorentz factor, which lies between its values at the ends
        // unless the velocity changes sign on the way, when it dips to 1
        let (g0, g1) = (u0.hypot(1.0), u1.hypot(1.0));
        let slowest = if u1 < 0.0 { 1.0 } else { g0.min(g1) };
        prop_assert!(tau1 <= dt1 / slowest * (1.0 + 1e-12));
        prop_assert!(tau1 >= dt1 / g0.max(g1) * (1.0 - 1e-12));
        // the velocity heads for α/H_Λ and never overshoots it
        let limit = u.drag_limited_proper_velocity(alpha);
        prop_assert!((u1 - limit).abs() <= (u0 - limit).abs() * (1.0 + 4.0 * EPS));
        Ok(())
    });
}

#[test]
fn de_sitter_thrust_reduces_to_flat_space_for_short_burns() {
    // H_Λ dt ≪ 1: τ = (asinh(u0 + α dt) − asinh(u0)) / α, the flat-space rocket
    let strategy = (cosmology(), -12.0..-7.0f64, 0.0..14.0f64, 0.0..30.0f64);
    check(strategy, |(u, log_h_dt, log_alpha, u0)| {
        let h = u.h_lambda();
        let h_dt = 10f64.powf(log_h_dt);
        let dt = h_dt / h;
        let alpha = 10f64.powf(log_alpha) * h;
        let tau = u.de_sitter_thrust_proper_time(u0, alpha, dt);
        // the rapidity gained, asinh(u1) − asinh(u0) = ln((u1 + γ1)/(u0 + γ0)), without the cancellation
        let u1 = u0 + alpha * dt;
        let (g0, g1) = (u0.hypot(1.0), u1.hypot(1.0));
        let flat = (alpha * dt * (1.0 + (u1 + u0) / (g1 + g0)) / (u0 + g0)).ln_1p() / alpha;
        // drag changes the velocity by a fraction H dt of itself over the burn
        let tol = 4.0 * h_dt * (1.0 + u0) + 1e-9;
        prop_assert!(close(tau, flat, tol), "{tau} against flat {flat}");
        Ok(())
    });
}

#[test]
fn the_sky_cools_to_the_de_sitter_temperature() {
    check((cosmology(), 0.1..10.0f64), |(u, t_cmb0)| {
        let t_ds = u.de_sitter_temperature();
        prop_assert!(t_ds > 0.0);
        let crossover = u.cmb_crossover_time(t_cmb0);
        prop_assert!(close(u.cmb_temperature(crossover, t_cmb0), t_ds, 1e-12));
        prop_assert!(close(u.cmb_temperature(u.age(), t_cmb0), t_cmb0, 8.0 * EPS));
        // before the crossover the sink is the microwave background; after it, the horizon
        prop_assert!(close(
            u.sink_temperature(0.5 * crossover, t_cmb0),
            u.cmb_temperature(0.5 * crossover, t_cmb0),
            EPS
        ));
        prop_assert_eq!(u.sink_temperature(2.0 * crossover, t_cmb0), t_ds);
        prop_assert_eq!(u.sink_temperature(1e9 * crossover, t_cmb0), t_ds);
        // the tidal acceleration is the same number in both unit systems
        let r = 3.0e20;
        let natural = u.lambda_tidal_acceleration(r / u.length_unit_metres());
        let si = natural * u.length_unit_metres() / (u.time_unit_seconds() * u.time_unit_seconds());
        prop_assert!(close(u.lambda_tidal_acceleration_si(r), si, 1e-13));
        Ok(())
    });
}

#[test]
fn out_of_domain_inputs_are_refused() {
    assert!(Cosmology::per_year(0.0, 0.3).is_none());
    assert!(Cosmology::per_year(-1.0, 0.3).is_none());
    assert!(Cosmology::per_year(f64::NAN, 0.3).is_none());
    assert!(Cosmology::per_year(f64::INFINITY, 0.3).is_none());
    assert!(Cosmology::per_year(7e-11, 0.0).is_none());
    assert!(Cosmology::per_year(7e-11, 1.0).is_none());
    assert!(Cosmology::per_year(7e-11, f64::NAN).is_none());
    assert!(Cosmology::per_second(2e-18, 0.3).is_some());
    assert!(Cosmology::with_time_unit(1.0, 0.3, 0.0).is_none());
    assert!(Cosmology::from_km_s_mpc(70.0, 1.5).is_none());
    assert!(Cosmology::from_h_lambda_per_year(5.7e-11, 1.0).is_none());
    assert!(Cosmology::PLANCK_2018.in_time_unit(-1.0).is_none());

    let u = Cosmology::PLANCK_2018;
    let now = u.age();
    for bad in [-1.0, f64::NAN] {
        assert!(u.scale_factor(bad).is_nan());
        assert!(u.ln_scale_factor(bad).is_nan());
        assert!(u.hubble(bad).is_nan());
        assert!(u.time_at_scale_factor(bad).is_nan());
        assert!(u.conformal_time(bad).is_nan());
        assert!(u.conformal_time_remaining(bad).is_nan());
        assert!(u.conformal_time_at_scale_factor(bad).is_nan());
        assert!(u.event_horizon_proper(bad).is_nan());
        assert!(u.reach_mass(bad).is_nan());
        assert!(u.scale_factor_integral(bad).is_nan());
        assert!(u.cmb_temperature(bad, 2.7).is_nan());
        assert!(u.sink_temperature(bad, 2.7).is_nan());
        assert!(u.coasting_proper_velocity(1.0, bad, now).is_nan());
        assert!(u.thrust_proper_velocity(1.0, 1.0, now, bad).is_nan());
        assert!(u.de_sitter_thrust_proper_time(0.0, 1.0, bad).is_nan());
        assert!(u.time_at_conformal_time(bad).is_none());
        assert!(u.time_at_conformal_time_remaining(bad).is_none());
        assert!(u.last_contact_time(bad).is_none());
        assert!(u.light_arrival_time(bad, 1.0).is_none());
        assert!(u.light_arrival_time(now, bad).is_none());
        assert!(u.light_emission_time(bad, 1.0).is_none());
        assert!(u.light_emission_time(now, bad).is_none());
        assert!(u.light_reception(bad, 1.0).is_none());
        assert!(u.light_reception(1.0, bad).is_none());
    }
    assert!(u.coasting_gamma(0.5, now, 2.0 * now).is_nan());

    // the ends of the domain
    assert_eq!(u.scale_factor(0.0), 0.0);
    assert_eq!(u.ln_scale_factor(0.0), f64::NEG_INFINITY);
    assert_eq!(u.hubble(0.0), f64::INFINITY);
    assert_eq!(u.conformal_time(0.0), 0.0);
    assert_eq!(u.conformal_time_remaining(0.0), u.conformal_time_limit());
    assert_eq!(u.event_horizon_proper(0.0), 0.0);
    assert_eq!(u.scale_factor_integral(0.0), 0.0);
    assert_eq!(u.time_at_scale_factor(0.0), 0.0);
    assert_eq!(u.time_at_scale_factor(f64::INFINITY), f64::INFINITY);
    assert_eq!(u.time_at_conformal_time(0.0), Some(0.0));
    assert_eq!(u.time_at_conformal_time(u.conformal_time_limit()), None);
    assert_eq!(u.time_at_conformal_time_remaining(0.0), None);
    assert_eq!(
        u.time_at_conformal_time_remaining(u.conformal_time_limit()),
        Some(0.0)
    );
    assert_eq!(
        u.time_at_conformal_time_remaining(2.0 * u.conformal_time_limit()),
        None
    );
    assert_eq!(u.last_contact_time(0.0), Some(f64::INFINITY));
    assert_eq!(u.light_arrival_time(now, 0.0), Some(now));
    assert_eq!(u.light_emission_time(now, 0.0), Some(now));
    assert_eq!(u.hubble(f64::INFINITY), u.h_lambda());
    assert_eq!(u.conformal_time(f64::INFINITY), u.conformal_time_limit());
    assert_eq!(u.conformal_time_remaining(f64::INFINITY), 0.0);
    assert_eq!(u.conformal_time_remaining(1e30), 0.0);
    assert_eq!(u.event_horizon_proper(f64::INFINITY), u.de_sitter_horizon());
    assert_eq!(u.de_sitter_thrust_proper_time(0.0, 0.0, 5.0), 5.0);
    assert_eq!(u.de_sitter_thrust_proper_time(3.0, 1.0, 0.0), 0.0);
    assert_eq!(u.cmb_temperature(1e30, 2.7255), 0.0);
}
