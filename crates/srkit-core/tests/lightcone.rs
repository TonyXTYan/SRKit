use proptest::prelude::*;
use srkit_core::kinematics::{Event, Vec3, C};
use srkit_core::lightcone::{future_intersection, retarded_state};
use srkit_core::segment::Segment;
use srkit_core::Error;

fn event(t: f64, x: f64, y: f64, z: f64) -> Event {
    Event::new(t, Vec3::new(x, y, z)).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn stationary_source_retarded_and_future_light_times() {
    let source = Segment::inertial(event(0.0, 0.0, 0.0, 0.0), Vec3::ZERO, 20.0).unwrap();
    let observer = event(12.0, 3.0 * C, 4.0 * C, 0.0);
    let state = retarded_state(&source, observer).unwrap();
    close(state.event.t, 7.0);
    close(state.proper_time, 7.0);

    let receiver = Segment::inertial(event(0.0, 3.0 * C, 4.0 * C, 0.0), Vec3::ZERO, 20.0).unwrap();
    let arrival = future_intersection(&receiver, event(2.0, 0.0, 0.0, 0.0)).unwrap();
    close(arrival.event.t, 7.0);
}

#[test]
fn moving_receiver_and_source_use_z_dimension() {
    let moving = Segment::inertial(
        event(0.0, 0.0, 0.0, 2.0 * C),
        Vec3::new(0.0, 0.0, 0.5 * C),
        20.0,
    )
    .unwrap();
    let arrival = future_intersection(&moving, event(0.0, 0.0, 0.0, 0.0)).unwrap();
    close(arrival.event.t, 4.0);
    close(arrival.event.position.z, 4.0 * C);

    let seen = retarded_state(&moving, event(11.0, 0.0, 0.0, 0.0)).unwrap();
    close(seen.event.t, 6.0);
    close(seen.event.position.z, 5.0 * C);
}

#[test]
fn intersections_at_segment_endpoints_are_included() {
    let receiver = Segment::inertial(event(2.0, 2.0 * C, 0.0, 0.0), Vec3::ZERO, 3.0).unwrap();
    close(
        future_intersection(&receiver, event(0.0, 0.0, 0.0, 0.0))
            .unwrap()
            .event
            .t,
        2.0,
    );
    close(
        retarded_state(&receiver, event(7.0, 0.0, 0.0, 0.0))
            .unwrap()
            .event
            .t,
        5.0,
    );
}

#[test]
fn roots_outside_segment_are_not_clamped_into_false_intersections() {
    let receiver = Segment::inertial(event(0.0, 10.0 * C, 0.0, 0.0), Vec3::ZERO, 5.0).unwrap();
    assert!(matches!(
        future_intersection(&receiver, event(0.0, 0.0, 0.0, 0.0)),
        Err(Error::NoIntersection)
    ));
    assert!(matches!(
        retarded_state(&receiver, event(9.0, 0.0, 0.0, 0.0)),
        Err(Error::NoIntersection)
    ));
    assert!(matches!(
        retarded_state(&receiver, event(-1.0, 0.0, 0.0, 0.0)),
        Err(Error::NoIntersection)
    ));
    assert!(matches!(
        future_intersection(&receiver, event(6.0, 0.0, 0.0, 0.0)),
        Err(Error::NoIntersection)
    ));
}

#[test]
fn null_intervals_choose_latest_retarded_and_earliest_future_root() {
    let ray = Segment::null(event(0.0, 0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), 10.0).unwrap();
    close(
        retarded_state(&ray, event(15.0, 15.0 * C, 0.0, 0.0))
            .unwrap()
            .event
            .t,
        10.0,
    );
    close(
        future_intersection(&ray, event(0.0, 0.0, 0.0, 0.0))
            .unwrap()
            .event
            .t,
        0.0,
    );
}

#[test]
fn parallel_null_motion_off_the_cone_has_no_root() {
    for offset in [-2.0, 2.0] {
        let ray = Segment::null(
            event(0.0, offset * C, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            10.0,
        )
        .unwrap();
        assert!(matches!(
            retarded_state(&ray, event(15.0, 15.0 * C, 0.0, 0.0)),
            Err(Error::NoIntersection)
        ));
        if offset > 0.0 {
            assert!(matches!(
                future_intersection(&ray, event(0.0, 0.0, 0.0, 0.0)),
                Err(Error::NoIntersection)
            ));
        } else {
            // The left-going flash meets a ray that starts behind the origin.
            close(
                future_intersection(&ray, event(0.0, 0.0, 0.0, 0.0))
                    .unwrap()
                    .event
                    .t,
                1.0,
            );
        }
    }
}

proptest! {
    #[test]
    fn inertial_arrivals_match_closed_form_and_are_on_causal_side(
        beta in -0.9_f64..0.9,
        distance_light_seconds in 0.1_f64..10.0,
    ) {
        let receiver = Segment::inertial(
            event(0.0, 0.0, 0.0, distance_light_seconds * C),
            Vec3::new(0.0, 0.0, beta * C),
            100.0,
        ).unwrap();
        let arrival = future_intersection(&receiver, event(0.0, 0.0, 0.0, 0.0)).unwrap();
        let expected = distance_light_seconds / (1.0 - beta);
        prop_assert!((arrival.event.t - expected).abs() <= 1e-12 * expected.max(1.0));
        let residual = arrival.event.t - arrival.event.position.norm() / C;
        prop_assert!(residual >= 0.0);
    }
}

#[test]
fn invalid_observation_event_is_rejected_even_when_public_fields_bypass_constructor() {
    let source = Segment::inertial(event(0.0, 0.0, 0.0, 0.0), Vec3::ZERO, 10.0).unwrap();
    let invalid = Event {
        t: f64::NAN,
        position: Vec3::ZERO,
    };
    assert!(matches!(
        retarded_state(&source, invalid),
        Err(Error::InvalidInput)
    ));
    let invalid = Event {
        t: 1.0,
        position: Vec3::new(0.0, f64::INFINITY, 0.0),
    };
    assert!(matches!(
        future_intersection(&source, invalid),
        Err(Error::InvalidInput)
    ));
}
