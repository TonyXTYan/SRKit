//! Exact causal predicates: integer arithmetic gives the same answer at any epoch and any distance.
use proptest::prelude::*;
use srkit_core::{
    causal::{self, Causal, GridVec},
    time::{Duration384, Duration512, Time384, Time512},
};

const C: u128 = 299_792_458;

#[test]
fn light_ticks_are_the_first_that_pass() {
    // exactly one light-second: 1000 ticks, not 999
    let one = GridVec::new(C as i128, 0, 0);
    assert_eq!(
        causal::light_ticks::<6>(one).unwrap(),
        Duration384::from_unsigned_millis(1000)
    );
    // one metre more needs one tick more
    assert_eq!(
        causal::light_ticks::<6>(GridVec::new(C as i128 + 1, 0, 0)).unwrap(),
        Duration384::from_unsigned_millis(1001)
    );
    assert_eq!(
        causal::light_ticks::<6>(GridVec::ZERO).unwrap(),
        Duration384::ZERO
    );
    // a 3-4-5 triangle of light-seconds
    let d = GridVec::new(3 * C as i128, -4 * C as i128, 0);
    assert_eq!(
        causal::light_ticks::<8>(d).unwrap(),
        Duration512::from_unsigned_millis(5000)
    );
}

#[test]
fn the_same_answer_at_any_epoch() {
    let d = GridVec::new(123_456_789_012, -98_765_432_109, 5_555_555_555);
    let ticks = causal::light_ticks::<6>(d).unwrap();
    for epoch in ["0", "1000000000000", &format!("1{}", "0".repeat(100))] {
        let e = Time384::from_decimal(epoch).unwrap();
        let arrival = causal::first_arrival(e, d).unwrap();
        assert_eq!(arrival.signed_duration_since(e), ticks);
        assert!(causal::arrived(e, arrival, d));
        assert!(!causal::arrived(
            e,
            arrival
                .checked_sub(Duration384::from_unsigned_millis(1))
                .unwrap(),
            d
        ));
    }
    // the far future in 512 bits, across the observable universe (8.8e26 m)
    let far = GridVec::new(880_000_000_000_000_000_000_000_000, 0, 0);
    let e = Time512::from_decimal(&format!("7{}", "0".repeat(140))).unwrap();
    let arrival = causal::first_arrival(e, far).unwrap();
    assert!(causal::arrived(e, arrival, far));
    assert!(!causal::arrived(
        e,
        arrival
            .checked_sub(Duration512::from_unsigned_millis(1))
            .unwrap(),
        far
    ));
    // before the emission nothing has arrived, even at no distance
    assert!(!causal::arrived(arrival, e, GridVec::ZERO));
    assert!(causal::arrived(e, e, GridVec::ZERO));
}

#[test]
fn classification_is_exact_on_the_cone() {
    // c * 1000 ms = 1 light-second exactly: null, and one metre either way tips it
    let dt = Duration384::from_unsigned_millis(1000);
    assert_eq!(
        causal::classify(dt, GridVec::new(C as i128, 0, 0)),
        Causal::Null
    );
    assert_eq!(
        causal::classify(dt, GridVec::new(C as i128 - 1, 0, 0)),
        Causal::Timelike
    );
    assert_eq!(
        causal::classify(dt, GridVec::new(C as i128, 1, 0)),
        Causal::Spacelike
    );
    // the sign of the interval does not matter
    assert_eq!(
        causal::classify(dt.negated(), GridVec::new(0, 0, -(C as i128))),
        Causal::Null
    );
    // a wait far longer than any distance is timelike without overflowing
    let ages = Time512::MAX.signed_duration_since(Time512::ZERO);
    assert_eq!(
        causal::classify(ages, GridVec::new(i128::MAX, i128::MIN + 1, i128::MAX)),
        Causal::Timelike
    );
}

#[test]
fn grid_rounding_and_overflow() {
    assert_eq!(
        GridVec::from_metres(1.4, -2.6, 0.5).unwrap(),
        GridVec::new(1, -3, 1)
    );
    assert!(GridVec::from_metres(f64::NAN, 0.0, 0.0).is_err());
    assert!(GridVec::from_metres(1e39, 0.0, 0.0).is_err());
    assert_eq!(
        GridVec::new(1, 2, 3).to(GridVec::new(4, 6, 8)).unwrap(),
        GridVec::new(3, 4, 5)
    );
    assert!(GridVec::new(i128::MIN, 0, 0)
        .to(GridVec::new(i128::MAX, 0, 0))
        .is_err());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    /// `light_ticks` is the smallest tick count that passes the exact check, for any displacement.
    #[test]
    fn light_ticks_is_minimal(x in any::<i64>(), y in any::<i64>(), z in any::<i64>(), scale in 0u32..60) {
        let d = GridVec::new((x as i128) << scale, (y as i128) << scale, (z as i128) << scale);
        let e = Time384::from_millis(12_345);
        let ticks = causal::light_ticks::<6>(d).unwrap();
        let arrival = e.checked_add(ticks).unwrap();
        prop_assert!(causal::arrived(e, arrival, d));
        if ticks != Duration384::ZERO {
            prop_assert!(!causal::arrived(e, arrival.checked_sub(Duration384::from_unsigned_millis(1)).unwrap(), d));
        }
        // along one axis it is the plain ceiling
        let axis = GridVec::new(x as i128, 0, 0);
        let expected = (x.unsigned_abs() as u128 * 1000).div_ceil(C);
        prop_assert_eq!(causal::light_ticks::<6>(axis).unwrap(), Duration384::from_unsigned_millis(expected as u64));
    }

    /// Settling from any estimate within the guard band lands on the same first tick.
    #[test]
    fn settle_does_not_depend_on_the_estimate(x in 1i64..i64::MAX / 4, off in -6i64..6) {
        let d = GridVec::new(x as i128, x as i128 / 3, 7);
        let e = Time384::from_millis(1 << 40);
        let exact = causal::first_arrival(e, d).unwrap();
        let estimate = if off >= 0 { exact.checked_add(Duration384::from_millis(off)).unwrap() } else { exact.checked_sub(Duration384::from_millis(-off)).unwrap_or(e) };
        prop_assert_eq!(causal::settle_arrival(e, estimate, 16, |_| Ok(d)).unwrap(), exact);
    }
}

mod finer_ticks {
    //! R71.4: a 1,024-bit clock with microsecond and nanosecond ticks.
    use proptest::prelude::*;
    use srkit_core::{
        causal::{self, Causal, GridVec, Resolution},
        time::{Duration1024, Duration512, Time1024, Time384, Time512},
    };

    const C: i128 = 299_792_458;

    #[test]
    fn a_nanosecond_needs_a_millimetre() {
        let (ms, us, ns) = (
            Resolution::MILLISECOND_METRE,
            Resolution::MICROSECOND_METRE,
            Resolution::NANOSECOND_MILLIMETRE,
        );
        assert!((ms.light_steps_per_tick() - 299_792.458).abs() < 1e-9);
        assert!((us.light_steps_per_tick() - 299.792458).abs() < 1e-12);
        assert!((ns.light_steps_per_tick() - 299.792458).abs() < 1e-12);
        // a nanosecond tick on a metre grid: light covers less than one grid step a tick, so rounding a
        // position to the grid moves an arrival by more than a tick
        assert!(
            Resolution::new(1_000_000_000, 1)
                .unwrap()
                .light_steps_per_tick()
                < 1.0
        );
        assert!(Resolution::new(0, 1).is_err() && Resolution::new(1, 0).is_err());
        assert!(
            Resolution::new(10_u64.pow(13), 1).is_err() && Resolution::new(1, 10_000_000).is_err()
        );
    }

    #[test]
    fn one_light_second_at_every_resolution() {
        // metres: 1e6 microseconds, and one metre more is one tick more (light crosses 300 m a microsecond: +1)
        let us = Resolution::MICROSECOND_METRE;
        assert_eq!(
            causal::light_ticks_at::<16>(us, GridVec::new(C, 0, 0)).unwrap(),
            Duration1024::from_unsigned_ticks(1_000_000)
        );
        assert_eq!(
            causal::light_ticks_at::<16>(us, GridVec::new(C + 1, 0, 0)).unwrap(),
            Duration1024::from_unsigned_ticks(1_000_001)
        );
        // millimetres: 1e9 nanoseconds
        let ns = Resolution::NANOSECOND_MILLIMETRE;
        let one = GridVec::new(C * 1000, 0, 0);
        assert_eq!(
            causal::light_ticks_at::<16>(ns, one).unwrap(),
            Duration1024::from_unsigned_ticks(1_000_000_000)
        );
        assert_eq!(
            causal::light_ticks_at::<16>(ns, GridVec::new(C * 1000 + 1, 0, 0)).unwrap(),
            Duration1024::from_unsigned_ticks(1_000_000_001)
        );
        assert_eq!(
            causal::classify_at(ns, Duration1024::from_unsigned_ticks(1_000_000_000), one),
            Causal::Null
        );
        assert_eq!(
            causal::classify_at(ns, Duration1024::from_unsigned_ticks(999_999_999), one),
            Causal::Spacelike
        );
        // a 3-4-5 triangle of light-seconds in millimetres
        let d = GridVec::new(3 * C * 1000, -4 * C * 1000, 0);
        assert_eq!(
            causal::light_ticks_at::<16>(ns, d).unwrap(),
            Duration1024::from_unsigned_ticks(5_000_000_000)
        );
        // the millisecond functions are the millisecond-metre resolution
        let ms = Resolution::MILLISECOND_METRE;
        let x = GridVec::new(123_456_789_012, -98_765_432_109, 5_555_555_555);
        assert_eq!(
            causal::light_ticks_at::<8>(ms, x).unwrap(),
            causal::light_ticks::<8>(x).unwrap()
        );
    }

    #[test]
    fn what_each_width_reaches() {
        // 1e100 years in nanoseconds is 3.16e116: 388 bits. It overflows 384 bits and fits 512 and 1,024.
        let ns_in_1e100_years = format!("31557600{}", "0".repeat(109));
        assert!(Time384::from_decimal(&ns_in_1e100_years).is_err());
        assert!(Time512::from_decimal(&ns_in_1e100_years).is_ok());
        // 1,024 bits at a nanosecond: an epoch of 1e290 years, exact to the tick, across the observable universe
        let epoch = Time1024::from_decimal(&format!("31557600{}", "0".repeat(299))).unwrap();
        let ns = Resolution::NANOSECOND_MILLIMETRE;
        let far = GridVec::new(880_000_000_000_000_000_000_000_000_000, 3, -7); // 8.8e26 m in mm
        let arrival = causal::first_arrival_at(ns, epoch, far).unwrap();
        assert!(causal::arrived_at(ns, epoch, arrival, far));
        assert!(!causal::arrived_at(
            ns,
            epoch,
            arrival
                .checked_sub(Duration1024::from_unsigned_ticks(1))
                .unwrap(),
            far
        ));
        assert_eq!(epoch.to_string().len(), 307);
        // and the same check in 512 bits at 1e100 years
        let e512 = Time512::from_decimal(&ns_in_1e100_years).unwrap();
        let a512 = causal::first_arrival_at(ns, e512, far).unwrap();
        assert_eq!(
            a512.signed_duration_since(e512).to_string(),
            arrival.signed_duration_since(epoch).to_string()
        );
        assert!(!causal::arrived_at(
            ns,
            e512,
            a512.checked_sub(Duration512::from_unsigned_ticks(1))
                .unwrap(),
            far
        ));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(2000))]

        /// At every resolution `light_ticks_at` is the smallest tick count that passes the exact check, and a
        /// moving receiver's arrival settles on the first tick that passes.
        #[test]
        fn minimal_at_every_resolution(x in any::<i64>(), y in any::<i64>(), z in any::<i64>(), scale in 0u32..50, which in 0usize..3, off in -5i64..5) {
            let res = [Resolution::MILLISECOND_METRE, Resolution::MICROSECOND_METRE, Resolution::NANOSECOND_MILLIMETRE][which];
            let d = GridVec::new((x as i128) << scale, (y as i128) << scale, (z as i128) << scale);
            let e = Time1024::from_decimal(&format!("7{}", "0".repeat(250))).unwrap();
            let ticks = causal::light_ticks_at::<16>(res, d).unwrap();
            let arrival = e.checked_add(ticks).unwrap();
            prop_assert!(causal::arrived_at(res, e, arrival, d));
            if ticks != Duration1024::ZERO {
                prop_assert!(!causal::arrived_at(res, e, arrival.checked_sub(Duration1024::from_unsigned_ticks(1)).unwrap(), d));
            }
            let estimate = if off >= 0 { arrival.checked_add(Duration1024::from_millis(off)).unwrap() } else { arrival.checked_sub(Duration1024::from_millis(-off)).unwrap_or(e) };
            prop_assert_eq!(causal::settle_arrival_at(res, e, estimate.max(e), 16, |_| Ok(d)).unwrap(), arrival);
        }
    }
}
