use srkit_core::{
    kinematics::{lorentz_factor, Event, Vec3, C, G0},
    segment::Segment,
    Error,
};

fn origin() -> Event {
    Event::new(0.0, Vec3::ZERO).unwrap()
}

#[test]
fn underflowed_rapidity_does_not_silently_erase_the_burn() {
    assert!(matches!(
        Segment::accelerated(origin(), Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0), 1e-320, 1.0),
        Err(Error::PrecisionLoss)
    ));
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(f64::MIN_POSITIVE),
        "{actual} != {expected} (relative tolerance {tolerance})"
    );
}

#[test]
fn arbitrarily_small_burn_retains_its_acceleration_displacement() {
    // cosh(q)-1 and expm1(q)-expm1(-q) constructions both lose this term
    // if evaluated by subtracting approximately equal unit-sized numbers.
    let acceleration = 1e-20;
    let segment = Segment::accelerated(
        origin(),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 1.0),
        acceleration,
        1.0,
    )
    .unwrap();
    let state = segment.state(segment.coord_time(1.0).unwrap()).unwrap();
    relative(state.event.position.z, 0.5 * acceleration, 1e-14);
    relative(state.velocity.z, acceleration, 1e-14);
    relative(state.proper_time, 1.0, 1e-14);
}

#[test]
fn near_inertial_limit_does_not_require_finite_c_over_acceleration() {
    let acceleration = 1e-300;
    assert!((C / acceleration).is_infinite());
    let segment = Segment::accelerated(
        origin(),
        Vec3::ZERO,
        Vec3::new(0.0, 1.0, 0.0),
        acceleration,
        1.0,
    )
    .unwrap();
    let state = segment.state(segment.coord_time(1.0).unwrap()).unwrap();
    relative(state.event.position.y, 0.5 * acceleration, 1e-13);
    relative(state.velocity.y, acceleration, 1e-13);
    assert_eq!(state.proper_time, 1.0);
}

#[test]
fn transverse_burn_matches_boosted_rest_frame_solution() {
    let gamma = 1.25;
    let q = 1.5;
    let tau = q * C / G0;
    let segment = Segment::accelerated(
        origin(),
        Vec3::new(0.6 * C, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        G0,
        tau,
    )
    .unwrap();
    for fraction in [1e-15, 0.01, 0.5, 1.0] {
        let proper = fraction * tau;
        let rapidity = fraction * q;
        let time = segment.coord_time(proper).unwrap();
        let state = segment.state(time).unwrap();
        relative(time, gamma * (C / G0) * libm::sinh(rapidity), 2e-14);
        relative(state.event.position.x, 0.6 * C * time, 2e-14);
        relative(
            state.event.position.z,
            2.0 * C * C / G0 * libm::sinh(rapidity / 2.0).powi(2),
            2e-14,
        );
        relative(state.velocity.x, 0.6 * C, 2e-14);
        relative(state.velocity.z, C / gamma * libm::tanh(rapidity), 2e-14);
        relative(state.proper_time, proper, 2e-14);
    }
}

#[test]
fn high_gamma_counterthrust_survives_the_turnaround() {
    // Previous gamma*cosh + b0*sinh formulas lost ~1e-5 in the inverse
    // already at q=1, and could produce a materially superluminal velocity.
    let initial_velocity = Vec3::new((1.0 - 1e-12) * C, 0.0, 0.0);
    let beta = initial_velocity.x / C;
    let eta = libm::atanh(beta);
    let gamma = lorentz_factor(initial_velocity).unwrap();
    let max_q = eta + 3.0;
    let segment = Segment::accelerated(
        origin(),
        initial_velocity,
        Vec3::new(-1.0, 0.0, 0.0),
        G0,
        max_q * C / G0,
    )
    .unwrap();
    for q in [0.0, 1e-10, 0.1, 1.0, eta / 2.0, eta, eta + 1.0, max_q] {
        let tau = q * C / G0;
        let time = segment.coord_time(tau).unwrap();
        let state = segment.state(time).unwrap();
        assert!(
            (state.proper_time - tau).abs() <= 0.1,
            "q={q}: {} != {tau}",
            state.proper_time
        );
        assert!(
            (state.velocity.x - C * libm::tanh(eta - q)).abs() <= 1.0,
            "q={q}: velocity {}",
            state.velocity.x
        );
        assert!(state.velocity.norm() < C);
        if q >= 0.1 {
            let expected_t = C / G0 * (gamma * beta - libm::sinh(eta - q));
            let expected_x = C * C / G0 * (gamma - libm::cosh(eta - q));
            relative(time, expected_t, 2e-12);
            relative(state.event.position.x, expected_x, 2e-12);
        }
    }
}

#[test]
fn unrepresentable_massive_velocity_is_rejected() {
    let result = Segment::accelerated(
        origin(),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 1.0),
        G0,
        100.0 * C / G0,
    );
    assert!(matches!(result, Err(Error::PrecisionLoss)));
}

#[test]
fn initial_velocity_is_preserved_exactly() {
    let velocity = Vec3::new(123_456.7, -234_567.8, 345_678.9);
    let segment =
        Segment::accelerated(origin(), velocity, Vec3::new(1.0, 2.0, 3.0), G0, 1.0).unwrap();
    assert_eq!(segment.state(0.0).unwrap().velocity, velocity);
}

#[test]
fn overflowing_circular_radius_is_rejected() {
    let axis = Vec3::new(f64::MAX, f64::MAX, 0.0);
    let result = Segment::circular(
        Event::new(0.0, axis).unwrap(),
        Vec3::ZERO,
        axis,
        Vec3::new(0.0, 0.0, f64::MAX),
        0.0,
        1.0,
    );
    assert!(matches!(result, Err(Error::InvalidInput)));
}
