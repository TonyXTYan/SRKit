use srkit_core::{
    kinematics::{Event, Vec3, C},
    lightcone::{future_intersection, retarded_state},
    segment::Segment,
    Error,
};

#[test]
fn accelerated_retarded_clock_matches_null_coordinate_closed_form() {
    let a = C / 10.0;
    let receiver = Event::new(10.0, Vec3::ZERO).unwrap();
    let source = Segment::accelerated(
        Event::new(0.0, Vec3::ZERO).unwrap(),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 1.0),
        a,
        20.0,
    )
    .unwrap();
    let state = retarded_state(&source, receiver).unwrap();
    let expected_tau = C / a * libm::log1p(a * receiver.t / C);
    assert!((state.proper_time - expected_tau).abs() < 1e-12);
    assert!(state.event.t <= receiver.t);
}

#[test]
fn accelerating_receiver_has_a_horizon_and_no_false_late_contact() {
    let a = C / 10.0;
    let receiver = Segment::accelerated(
        Event::new(0.0, Vec3::ZERO).unwrap(),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 1.0),
        a,
        40.0,
    )
    .unwrap();
    let emission = Event::new(0.0, Vec3::new(0.0, 0.0, -5.0 * C)).unwrap();
    let state = future_intersection(&receiver, emission).unwrap();
    // ct - z = c²/a (1-exp(-a*tau/c)); 5 light-seconds < c/a.
    let expected_tau = -C / a * libm::log1p(-0.5);
    assert!((state.proper_time - expected_tau).abs() < 1e-12);
    let outside = Event::new(0.0, Vec3::new(0.0, 0.0, -11.0 * C)).unwrap();
    assert!(matches!(
        future_intersection(&receiver, outside),
        Err(Error::NoIntersection)
    ));
}

#[test]
fn circular_lightcone_query_preserves_the_numerical_causal_bracket() {
    let radius = 0.5 * C;
    let orbit = Segment::circular(
        Event::new(0.0, Vec3::new(radius, 0.0, 0.0)).unwrap(),
        Vec3::ZERO,
        Vec3::new(radius, 0.0, 0.0),
        Vec3::new(0.0, 0.0, radius),
        0.8,
        20.0,
    )
    .unwrap();
    let observer = Event::new(10.0, Vec3::new(0.0, C, C)).unwrap();
    let state = retarded_state(&orbit, observer).unwrap();
    let residual =
        observer.t - state.event.t - (observer.position - state.event.position).norm() / C;
    assert!(residual >= 0.0);
    assert!(residual < 1e-13);
}
