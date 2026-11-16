use proptest::prelude::*;
use srkit_core::kinematics::{
    aberration, doppler_factor, lorentz_factor, Event, FourVector, Vec3, C,
};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b} (tolerance {tolerance})"
    );
}

#[test]
fn gamma_and_invalid_domains() {
    close(
        lorentz_factor(Vec3::new(0.6 * C, 0.0, 0.0)).unwrap(),
        1.25,
        1e-15,
    );
    close(
        lorentz_factor(Vec3::new(0.0, 0.0, 0.8 * C)).unwrap(),
        5.0 / 3.0,
        1e-15,
    );
    for v in [
        Vec3::new(C, 0.0, 0.0),
        Vec3::new(0.0, -C, 0.0),
        Vec3::new(f64::NAN, 0.0, 0.0),
        Vec3::new(f64::INFINITY, 0.0, 0.0),
    ] {
        assert!(lorentz_factor(v).is_err());
    }
    assert!(Event::new(f64::NAN, Vec3::ZERO).is_err());
    assert!(Event::new(0.0, Vec3::new(0.0, 0.0, f64::INFINITY)).is_err());
    assert!(Vec3::ZERO.normalized().is_err());
    let large = Vec3::new(1e250, 0.0, 1e250).normalized().unwrap();
    close(large.norm(), 1.0, 2e-16);
}

#[test]
fn boost_convention_and_tiny_beta() {
    let rest = FourVector::new(1.0, Vec3::ZERO).unwrap();
    let boosted = rest.boost(Vec3::new(0.0, 0.0, 0.6)).unwrap();
    close(boosted.temporal, 1.25, 1e-15);
    close(boosted.spatial.z, -0.75, 1e-15);
    close(boosted.dot(boosted), 1.0, 1e-15);
    let vector = FourVector::new(2.0, Vec3::new(3.0, 4.0, 5.0)).unwrap();
    let tiny = vector.boost(Vec3::new(1e-20, 0.0, 0.0)).unwrap();
    assert_eq!(tiny, vector);
    assert!(vector.boost(Vec3::new(1.0, 0.0, 0.0)).is_err());
}

#[test]
fn doppler_conventions_and_transverse_effect() {
    let n = Vec3::new(1.0, 0.0, 0.0);
    let velocity = Vec3::new(0.6 * C, 0.0, 0.0);
    // Receiver pursues the photon: redshift. Emitter pursues receiver: blueshift.
    close(doppler_factor(Vec3::ZERO, velocity, n).unwrap(), 0.5, 1e-15);
    close(doppler_factor(velocity, Vec3::ZERO, n).unwrap(), 2.0, 1e-15);
    close(doppler_factor(velocity, velocity, n).unwrap(), 1.0, 1e-15);
    close(
        doppler_factor(Vec3::new(0.0, 0.6 * C, 0.0), Vec3::ZERO, n).unwrap(),
        0.8,
        1e-15,
    );
    close(
        doppler_factor(Vec3::ZERO, Vec3::new(0.0, 0.0, 0.6 * C), n).unwrap(),
        1.25,
        1e-15,
    );
    assert!(doppler_factor(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO).is_err());
}

#[test]
fn photon_aberration_handles_z() {
    let direction = aberration(Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.6, 0.0, 0.0)).unwrap();
    close(direction.x, -0.6, 1e-15);
    close(direction.y, 0.0, 1e-15);
    close(direction.z, 0.8, 1e-15);
    let original = aberration(direction, Vec3::new(-0.6, 0.0, 0.0)).unwrap();
    close(original.x, 0.0, 1e-15);
    close(original.z, 1.0, 1e-15);
}

proptest! {
    #[test]
    fn arbitrary_boost_preserves_interval_and_is_invertible(
        temporal in -10.0f64..10.0,
        x in -10.0f64..10.0, y in -10.0f64..10.0, z in -10.0f64..10.0,
        bx in -0.5f64..0.5, by in -0.5f64..0.5, bz in -0.5f64..0.5,
    ) {
        let vector = FourVector::new(temporal, Vec3::new(x, y, z)).unwrap();
        let beta = Vec3::new(bx, by, bz);
        let transformed = vector.boost(beta).unwrap();
        let restored = transformed.boost(-beta).unwrap();
        prop_assert!((transformed.dot(transformed) - vector.dot(vector)).abs() < 2e-11);
        prop_assert!((restored.temporal - temporal).abs() < 1e-12);
        prop_assert!((restored.spatial - vector.spatial).norm() < 1e-12);
    }
}
