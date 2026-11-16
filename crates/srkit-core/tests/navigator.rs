use srkit_core::{
    kinematics::{Vec3, C, G0, JULIAN_YEAR},
    navigator::{reach, traveller_emission, Orbit, Trip, G_LY},
    trip::FlipBurn,
    Error,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}, tolerance {tolerance}"
    );
}

#[test]
fn one_way_trip_matches_the_independent_si_core_and_proxima_reference() {
    for distance in [4.2465, 1e5, 2.5e6] {
        let trip = Trip::new(distance, 1.0, f64::INFINITY, 0, false).unwrap();
        let si = FlipBurn::new(distance * C * JULIAN_YEAR, G0).unwrap();
        near(trip.tau, si.proper_time / JULIAN_YEAR, 2e-14);
        near(trip.t_end, si.coordinate_time / JULIAN_YEAR, 2e-14);
        near(trip.g_peak, si.peak_gamma, 2e-14);
        near(trip.lag, si.light_time_excess / JULIAN_YEAR, 2e-14);
    }
    let proxima = Trip::new(4.2465, 1.0, f64::INFINITY, 0, false).unwrap();
    near(proxima.tau, 3.542, 2e-4);
    near(proxima.t_end, 5.873, 1e-4);
}

#[test]
fn stopped_capped_and_round_trips_preserve_time_and_null_inverses() {
    for distance in [4.2465, 1e5, 2.5e6] {
        for cap in [2.0, 100.0, f64::INFINITY] {
            for stops in [0, 2] {
                for roundtrip in [false, true] {
                    let trip = Trip::new(distance, 0.1, cap, stops, roundtrip).unwrap();
                    assert!(trip.g_peak <= cap);
                    for fraction in [0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.99, 1.0] {
                        let tau = trip.tau * fraction;
                        let state = trip.state(tau).unwrap();
                        let tu = trip.tau_u(state.u).unwrap();
                        let tv = trip.tau_v(state.v).unwrap();
                        // Null coordinates accumulated across a return trip
                        // can be millions of years while their derivative is
                        // exp(-phi). Inverting their rounded f64 values is
                        // conditioned by that derivative, even analytically.
                        let u_condition =
                            8.0 * f64::EPSILON * state.u.abs() * libm::exp(state.dir * state.phi);
                        let v_condition =
                            8.0 * f64::EPSILON * state.v.abs() * libm::exp(-state.dir * state.phi);
                        assert!((tu-tau).abs() <= 2e-8*tau.abs().max(1.0)+u_condition, "u: d={distance},cap={cap},stops={stops},rt={roundtrip},f={fraction}: {tu}!={tau}");
                        assert!((tv-tau).abs() <= 2e-8*tau.abs().max(1.0)+v_condition, "v: d={distance},cap={cap},stops={stops},rt={roundtrip},f={fraction}: {tv}!={tau}");
                        assert!(state.gamma >= 1.0);
                        assert!(state.t >= state.tau - 1e-12);
                    }
                    let end = trip.state(trip.tau).unwrap();
                    assert_eq!(end.phi, 0.0);
                    assert_eq!(end.x, if roundtrip { 0.0 } else { distance });
                }
            }
        }
    }
}

#[test]
fn return_leg_mirrors_outbound_and_waiting_extends_both_clocks() {
    let trip = Trip::new(100_000.0, 1.0, 1e5, 1, true).unwrap();
    let half = trip.tau / 2.0;
    for fraction in [0.01, 0.1, 0.37, 0.6, 0.9] {
        let outward = trip.state(half * fraction).unwrap();
        let inward = trip.state(half * (1.0 + fraction)).unwrap();
        assert!((inward.x - (trip.d - outward.x)).abs() < 4e-15 * trip.d);
        near(inward.t, trip.t_out + outward.t, 2e-13);
        near(outward.phi, inward.phi, 2e-13);
        assert_eq!(inward.dir, -1.0);
    }
    let before = trip.state(-1e6).unwrap();
    assert_eq!((before.t, before.x, before.gamma), (-1e6, 0.0, 1.0));
    let late = trip.state(trip.tau + 1e6).unwrap();
    assert_eq!(late.x, 0.0);
    near(late.t, trip.t_end + 1e6, 1e-15);
    near(trip.tau_u(late.u).unwrap(), late.tau, 1e-15);
    near(trip.tau_v(late.v).unwrap(), late.tau, 1e-15);
}

#[test]
fn off_axis_received_light_and_ping_satisfy_the_null_equation() {
    for (distance, cap) in [(4.2465, f64::INFINITY), (1e5, 1e4)] {
        let trip = Trip::new(distance, 1.0, cap, 0, true).unwrap();
        for (t, along, perp) in [
            (-20.0, 1.0, 2.0),
            (1.0, 0.5, 1.0),
            (50.0, 2.0, 3.0),
            (1e6, 1.0, 1.0),
        ] {
            let receive = trip.intersection(t, along, perp, false).unwrap();
            let ping = trip.intersection(t, along, perp, true).unwrap();
            assert!(ping <= receive);
            for (tau, send) in [(receive, false), (ping, true)] {
                let state = trip.state(tau).unwrap();
                let radius = libm::hypot(state.x - along, perp);
                near(state.t + if send { radius } else { -radius }, t, 2e-9);
            }
        }
    }
}

#[test]
fn orbit_clock_and_center_contact_have_closed_form_solutions() {
    for radius in [1.0, 100.0, 50_000.0] {
        let orbit = Orbit::new(radius, 1.0).unwrap();
        near(orbit.gamma, libm::sqrt(1.0 + G_LY * radius), 2e-15);
        near(orbit.lap_tau, orbit.lap_t / orbit.gamma, 2e-15);
        for t in [-1e6, 0.0, 10.0, 1e6] {
            near(
                orbit.intersection(0.3, t, Vec3::ZERO, false).unwrap(),
                (t + radius) / orbit.gamma,
                2e-13,
            );
            near(
                orbit.intersection(0.3, t, Vec3::ZERO, true).unwrap(),
                (t - radius) / orbit.gamma,
                2e-13,
            );
        }
        let state = orbit.state(0.3, orbit.lap_tau).unwrap();
        near(state.position.x, radius * libm::cos(0.3), 2e-14);
        near(state.position.y, radius * libm::sin(0.3), 2e-14);
        near(state.velocity.norm(), orbit.beta, 2e-15);
    }
}

#[test]
fn off_plane_orbit_contacts_are_null_and_ordered() {
    let orbit = Orbit::new(100.0, 1.0).unwrap();
    let event = Vec3::new(30.0, 70.0, 50.0);
    let (ping, receive) = (
        orbit.intersection(0.5, 123.0, event, true).unwrap(),
        orbit.intersection(0.5, 123.0, event, false).unwrap(),
    );
    assert!(ping < receive);
    for (tau, send) in [(ping, true), (receive, false)] {
        let state = orbit.state(0.5, tau).unwrap();
        let radius = (state.position - event).norm();
        near(state.t + if send { radius } else { -radius }, 123.0, 3e-12);
    }
}

#[test]
fn reach_from_rest_is_radial_and_brake_now_uses_the_stopping_distance() {
    let a = G_LY;
    let budget = 3.0;
    let radius = 4.0 / a * libm::sinh(a * budget / 4.0).powi(2);
    for heading in [0.0, 0.3, 1.5, 3.0] {
        let result = reach(0.0, a, budget, heading).unwrap().unwrap();
        near(result[0], radius * libm::cos(heading), 2e-12);
        near(result[1], radius * libm::sin(heading), 2e-12);
    }
    let phi = 1.5;
    assert_eq!(reach(phi, a, phi / a - 0.1, 0.0).unwrap(), None);
    let stopped = reach(phi, a, phi / a, 0.3).unwrap().unwrap();
    near(stopped[0], (libm::cosh(phi) - 1.0) / a, 2e-12);
    near(stopped[1], 0.0, 2e-12);
}

#[test]
fn zero_distance_and_invalid_inputs_have_explicit_outcomes() {
    let zero = Trip::new(0.0, 1.0, 1.0, 0, false).unwrap();
    assert_eq!(zero.tau, 0.0);
    assert_eq!(zero.state(20.0).unwrap().t, 20.0);
    assert_eq!(zero.tau_u(12.0), Ok(12.0));
    assert!(matches!(
        Trip::new(1.0, 1.0, 1.0, 0, false),
        Err(Error::InvalidInput)
    ));
    assert!(Trip::new(-1.0, 1.0, f64::INFINITY, 0, false).is_err());
    assert!(Trip::new(1.0, 0.0, f64::INFINITY, 0, false).is_err());
    assert!(Trip::new(1.0, 1.0, f64::INFINITY, u32::MAX, false).is_err());
    assert!(Orbit::new(0.0, 1.0).is_err());
    assert!(reach(0.0, 0.0, 1.0, 0.0).is_err());
}

#[test]
fn test_traveller_contact_crosses_launch_without_a_lookback_limit() {
    let heading = Vec3::new(0.0, 0.0, 1.0);
    assert_eq!(
        traveller_emission(0.9, heading, 1.0, Vec3::new(0.0, 0.0, 5.0)),
        Ok(-4.0)
    );
    near(
        traveller_emission(0.9, heading, 10.0, Vec3::ZERO).unwrap(),
        10.0 / 1.9,
        2e-15,
    );
    near(
        traveller_emission(0.9, heading, 1e9, Vec3::ZERO).unwrap(),
        1e9 / 1.9,
        2e-15,
    );
    assert!(traveller_emission(1.0, heading, 10.0, Vec3::ZERO).is_err());
}

#[test]
fn off_axis_test_traveller_emission_is_on_the_observers_past_light_cone() {
    let heading = Vec3::new(1.0, 2.0, 3.0).normalized().unwrap();
    let position = Vec3::new(3.0, -5.0, 4.0);
    for beta in [0.0, 0.5, 0.9, 0.99] {
        let emission = traveller_emission(beta, heading, 123.0, position).unwrap();
        assert!((0.0..=123.0).contains(&emission));
        near(
            emission + (position - heading * (beta * emission)).norm(),
            123.0,
            3e-15,
        );
    }
}
