use srkit_core::kinematics::{Event, Vec3, C};
use srkit_core::segment::Segment;
use srkit_core::trip::{circular_gamma, FlipBurn};
use srkit_core::Error;

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn normalization_scales_subnormal_and_near_maximum_vectors() {
    for scale in [f64::from_bits(1), f64::MAX] {
        let direction = Vec3::new(scale, scale, 0.0).normalized().unwrap();
        relative(direction.norm(), 1.0, 3e-16);
        relative(direction.x, core::f64::consts::FRAC_1_SQRT_2, 3e-16);
        let ray = Segment::null(
            Event::new(0.0, Vec3::ZERO).unwrap(),
            Vec3::new(scale, scale, 0.0),
            1.0,
        )
        .unwrap();
        relative(ray.state(1.0).unwrap().velocity.norm(), C, 3e-16);
    }
}

#[test]
fn tiny_distance_has_finite_nonrelativistic_clocks_and_excess() {
    let trip = FlipBurn::new(1e-300, 1.0).unwrap();
    relative(trip.coordinate_time, 2e-150, 1e-14);
    relative(trip.proper_time, trip.coordinate_time, 1e-14);
    relative(trip.light_time_excess, trip.coordinate_time, 1e-14);
    assert!(trip.proper_time <= trip.coordinate_time);
}

#[test]
fn tiny_acceleration_does_not_overflow_an_intermediate_timescale() {
    let trip = FlipBurn::new(1.0, 1e-300).unwrap();
    relative(trip.coordinate_time, 2e150, 1e-14);
    relative(trip.proper_time, trip.coordinate_time, 1e-14);
    relative(trip.light_time_excess, trip.coordinate_time, 1e-14);
    assert!(trip.proper_time <= trip.coordinate_time);
}

#[test]
fn minimum_positive_inputs_preserve_representable_clock_results() {
    let tiny = f64::from_bits(1);
    let trip = FlipBurn::new(tiny, tiny).unwrap();
    relative(trip.coordinate_time, 2.0, 1e-14);
    relative(trip.proper_time, 2.0, 1e-14);
    relative(trip.light_time_excess, 2.0, 1e-14);
    for acceleration in [tiny, 1.0, f64::MAX] {
        let zero_trip = FlipBurn::new(0.0, acceleration).unwrap();
        assert_eq!(zero_trip.coordinate_time, 0.0);
        assert_eq!(zero_trip.proper_time, 0.0);
        assert_eq!(zero_trip.peak_gamma, 1.0);
    }
}

#[test]
fn large_circular_gamma_avoids_squaring_the_result() {
    relative(circular_gamma(1e308, 1e308).unwrap(), 1e308 / C, 1e-14);
}

#[test]
fn genuinely_unrepresentable_trip_outputs_still_return_overflow() {
    assert!(matches!(
        FlipBurn::new(f64::MAX, f64::from_bits(1)),
        Err(Error::Overflow)
    ));
    assert!(matches!(
        FlipBurn::new(f64::MAX, f64::MAX),
        Err(Error::Overflow)
    ));
}
