use proptest::prelude::*;
use srkit_core::{
    kinematics::{Event, Vec3, C, G0},
    segment::{Segment, Worldline},
    trip::FlipBurn,
    Error,
};

fn origin() -> Event {
    Event::new(0.0, Vec3::ZERO).unwrap()
}
fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

#[test]
fn analytic_acceleration_from_rest_and_small_time() {
    let segment =
        Segment::accelerated(origin(), Vec3::ZERO, Vec3::new(0., 0., 1.), G0, 1e7).unwrap();
    for tau in [0., 1e-9, 0.01, 100., 1e6, 1e7] {
        let state = segment.state(segment.coord_time(tau).unwrap()).unwrap();
        near(state.proper_time, tau, 1e-12);
        near(
            state.event.position.z,
            2. * (C * C / G0) * libm::sinh(G0 * tau / (2. * C)).powi(2),
            1e-12,
        );
        near(state.velocity.z, C * libm::tanh(G0 * tau / C), 1e-12);
    }
}

#[test]
fn noncollinear_acceleration_matches_lorentz_transformed_rest_burn() {
    use srkit_core::kinematics::FourVector;
    let beta = Vec3::new(0.4, -0.2, 0.3);
    let direction = Vec3::new(1., 2., -1.).normalized().unwrap();
    let segment = Segment::accelerated(origin(), beta * C, direction, G0, 1e7).unwrap();
    let tau = 5e6;
    let q = G0 * tau / C;
    let rest_displacement = FourVector::new(
        C * C / G0 * libm::sinh(q),
        direction * (2. * C * C / G0 * libm::sinh(q / 2.).powi(2)),
    )
    .unwrap();
    let lab = rest_displacement.boost(-beta).unwrap();
    let state = segment.state(lab.temporal / C).unwrap();
    near((state.event.position - lab.spatial).norm(), 0., 0.5);
    near(state.proper_time, tau, 1e-12);
}

#[test]
fn circle_is_3d_with_constant_clock_rate() {
    let r = C;
    let segment = Segment::circular(
        Event::new(2., Vec3::new(0., 0., r)).unwrap(),
        Vec3::ZERO,
        Vec3::new(0., 0., r),
        Vec3::new(0., r, 0.),
        0.6,
        10.,
    )
    .unwrap();
    let state = segment
        .state(2. + core::f64::consts::FRAC_PI_2 / 0.6)
        .unwrap();
    near(state.event.position.y, r, 1e-12);
    near(state.event.position.z, 0., 1e-6);
    near(state.velocity.norm(), 0.6 * C, 1e-12);
    near(state.proper_time, (state.event.t - 2.) * 0.8, 1e-12);
    assert!(Segment::circular(
        origin(),
        Vec3::ZERO,
        Vec3::new(r, 0., 0.),
        Vec3::new(r, 0., 0.),
        0.6,
        10.
    )
    .is_err());
}

#[test]
fn worldline_joins_and_null_clock_domain() {
    let v = Vec3::new(0., 0., 0.3 * C);
    let a = Segment::inertial(origin(), v, 10.).unwrap();
    let b = Segment::inertial(a.end(), v, 20.).unwrap();
    let segments = [a, b];
    let worldline = Worldline::new(&segments, 1e-6, 1e-6).unwrap();
    let state = worldline.state(25.).unwrap();
    near(worldline.coord_time(state.proper_time).unwrap(), 25., 1e-12);
    let invalid = [a, Segment::inertial(a.end(), -v, 20.).unwrap()];
    assert!(matches!(
        Worldline::new(&invalid, 1e-6, 1e-6),
        Err(Error::Discontinuous)
    ));
    let null = Segment::null(origin(), Vec3::new(0., 0., 1.), 10.).unwrap();
    assert_eq!(null.proper_time(5.), Ok(0.));
    assert_eq!(null.coord_time(0.), Err(Error::OutOfRange));
}

#[test]
fn trip_preserves_small_light_time_excess_at_deep_time() {
    let trip = FlipBurn::new(C * 1e100, G0).unwrap();
    near(trip.light_time_excess, 2. * C / G0, 1e-12);
    assert_eq!(trip.coordinate_time, trip.distance / C);
    assert!(trip.photon_rocket_mass_ratio().is_ok());
    assert!(FlipBurn::new(-1., G0).is_err());
    assert!(FlipBurn::new(1., 0.).is_err());
    assert_eq!(FlipBurn::new(0., G0).unwrap().proper_time, 0.);
}

#[test]
fn nonfinite_and_unrepresentable_intervals_are_rejected() {
    assert!(Segment::inertial(origin(), Vec3::new(C, 0., 0.), 1.).is_err());
    assert!(Segment::inertial(origin(), Vec3::ZERO, 0.).is_err());
    assert!(Segment::inertial(
        Event {
            t: f64::NAN,
            position: Vec3::ZERO
        },
        Vec3::ZERO,
        1.
    )
    .is_err());
    assert!(Segment::inertial(
        Event {
            t: 1e100,
            position: Vec3::ZERO
        },
        Vec3::ZERO,
        1.
    )
    .is_err());
    assert!(Segment::accelerated(origin(), Vec3::ZERO, Vec3::ZERO, G0, 1.).is_err());
}

proptest! {
    #[test]
    fn acceleration_clock_inverse_is_monotone(beta in -0.8f64..0.8, tau in 1f64..1e7) {
        let segment = Segment::accelerated(origin(), Vec3::new(beta*C,0.,0.), Vec3::new(0.,1.,1.), G0, 2e7).unwrap();
        let t = segment.coord_time(tau).unwrap();
        let state = segment.state(t).unwrap();
        prop_assert!((state.proper_time/tau-1.).abs() < 1e-10);
        prop_assert!(state.proper_time <= t*(1.+1e-12));
        prop_assert!(state.velocity.norm() < C);
        prop_assert!(segment.coord_time(tau+1.).unwrap() > t);
    }
}
