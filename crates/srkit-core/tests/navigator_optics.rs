use srkit_core::{
    kinematics::Vec3,
    navigator::{boost, observe},
    Error,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual} != {expected}"
    );
}

#[test]
fn optical_doppler_convention_and_exact_axes() {
    let phi = 5.0;
    let ahead = observe(phi, 1.0, 0.0, 0.0, true).unwrap();
    let behind = observe(phi, -1.0, 0.0, 0.0, true).unwrap();
    let transverse = observe(phi, 0.0, 1.0, 0.0, true).unwrap();
    near(ahead[0], libm::exp(phi), 2e-15);
    near(behind[0], libm::exp(-phi), 2e-15);
    near(transverse[0], libm::cosh(phi), 2e-15);
    assert_eq!(ahead[3], 1.0);
    assert_eq!(behind[3], -1.0);
    near(transverse[3], libm::tanh(phi), 2e-15);
    let high_behind = observe(100.0, -1.0, 0.0, 0.0, true).unwrap();
    near(high_behind[0], libm::exp(-100.0), 2e-15);
}

#[test]
fn aberration_toggle_retains_frequency_shift_and_direction_azimuth() {
    let (a, right, up) = (0.6, 0.48, 0.64);
    let shifted = observe(2.0, a, right, up, true).unwrap();
    let static_sky = observe(2.0, a, right, up, false).unwrap();
    assert_eq!(shifted[0], static_sky[0]);
    near(static_sky[3], a, 2e-15);
    near(static_sky[4], right, 2e-15);
    near(static_sky[5], up, 2e-15);
    near(shifted[4] / shifted[5], right / up, 2e-15);
    near(
        shifted[3] * shifted[3] + shifted[4] * shifted[4] + shifted[5] * shifted[5],
        1.0,
        2e-15,
    );
    assert!(shifted[2] < shifted[1]);
}

#[test]
fn full_three_dimensional_boost_matches_invariant_and_inverse() {
    let rapidity = Vec3::new(0.3, -0.4, 0.7);
    let position = Vec3::new(2.0, -3.0, 1.5);
    let time = 5.0;
    let (shifted_time, shifted_position) = boost(rapidity, time, position).unwrap();
    near(
        shifted_time * shifted_time - shifted_position.norm_squared(),
        time * time - position.norm_squared(),
        2e-14,
    );
    let (recovered_time, recovered_position) =
        boost(-rapidity, shifted_time, shifted_position).unwrap();
    near(recovered_time, time, 2e-14);
    near((recovered_position - position).norm() + 1.0, 1.0, 2e-14);
}

#[test]
fn null_coordinate_boost_keeps_tiny_forward_light_displacement() {
    let (time, position) = boost(Vec3::new(60.0, 0.0, 0.0), 1.0, Vec3::new(1.0, 0.0, 0.0)).unwrap();
    near(time, libm::exp(-60.0), 2e-15);
    assert_eq!(time, position.x);
    assert_eq!(position.y, 0.0);
    let (opposite_time, opposite_position) =
        boost(Vec3::new(60.0, 0.0, 0.0), 1.0, Vec3::new(-1.0, 0.0, 0.0)).unwrap();
    near(opposite_time, libm::exp(60.0), 2e-15);
    assert_eq!(opposite_time, -opposite_position.x);
}

#[test]
fn invalid_and_unrepresentable_optics_are_explicit_errors() {
    assert_eq!(observe(0.0, 0.0, 0.0, 0.0, true), Err(Error::InvalidInput));
    assert_eq!(
        observe(f64::NAN, 1.0, 0.0, 0.0, true),
        Err(Error::InvalidInput)
    );
    assert_eq!(observe(1000.0, 1.0, 0.0, 0.0, true), Err(Error::Overflow));
    assert_eq!(
        boost(Vec3::new(1000.0, 0.0, 0.0), 1.0, Vec3::ZERO),
        Err(Error::Overflow)
    );
}
