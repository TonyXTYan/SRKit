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

mod any_clock {
    //! R71.4: four widths, each at any tick and any grid spacing, behind one `Clock`.
    use proptest::prelude::*;
    use srkit_core::{
        causal::{Causal, Clock, Clock1024, Clock384, Clock512, Clock768, GridVec, Resolution},
        time::{Duration, Duration768, Time768, Timestamp},
        Error,
    };

    const C: i128 = 299_792_458;

    /// Ticks to the second and grid steps to the metre: the three named resolutions, a picosecond on a
    /// micrometre grid (the finest allowed), a second on a metre grid (the coarsest), and two that are not
    /// powers of ten (a 60 Hz frame on a centimetre grid; 1,024 ticks on a 3-to-the-metre grid).
    const RESOLUTIONS: [(u64, u64); 7] = [
        (1_000, 1),
        (1_000_000, 1),
        (1_000_000_000, 1_000),
        (1_000_000_000_000, 1_000_000),
        (1, 1),
        (60, 100),
        (1_024, 3),
    ];

    /// The first-arrival tick is the first that passes the exact check, and is the same number of ticks at
    /// every width.
    fn first_tick_passes<const N: usize>(res: Resolution, epoch: &str, dx: GridVec) -> String {
        let clock = Clock::<N>::new(res);
        let sent = Timestamp::<N>::from_decimal(epoch).unwrap();
        let heard = clock.first_arrival(sent, dx).unwrap();
        assert!(clock.arrived(sent, heard, dx));
        let wait = heard.duration_since(sent).unwrap();
        assert_eq!(wait, clock.light_ticks(dx).unwrap());
        if wait != Duration::<N>::ZERO {
            let before = heard.checked_sub(Duration::from_unsigned_ticks(1)).unwrap();
            assert!(!clock.arrived(sent, before, dx));
            assert_eq!(
                clock.classify(before.duration_since(sent).unwrap(), dx),
                Causal::Spacelike
            );
        }
        assert_ne!(clock.classify(wait, dx), Causal::Spacelike);
        wait.to_string()
    }

    #[test]
    fn one_light_second_on_every_clock() {
        for (ticks, grid) in RESOLUTIONS {
            let res = Resolution::new(ticks, grid).unwrap();
            let one = GridVec::new(C * grid as i128, 0, 0);
            let epoch = "3000000000000000000000000000000000000000000000000000000000000";
            let waits = [
                first_tick_passes::<6>(res, epoch, one),
                first_tick_passes::<8>(res, epoch, one),
                first_tick_passes::<12>(res, epoch, one),
                first_tick_passes::<16>(res, epoch, one),
            ];
            // exactly one second of ticks, at every width
            assert!(
                waits.iter().all(|w| *w == ticks.to_string()),
                "{ticks} {grid}: {waits:?}"
            );
            // and it is null, not merely inside
            assert_eq!(
                Clock768::new(res).classify(Duration768::from_unsigned_ticks(ticks), one),
                Causal::Null
            );
            // one grid step further is strictly later whenever the grid is finer than a tick of light
            let further = GridVec::new(C * grid as i128 + 1, 0, 0);
            let later = Clock768::new(res).light_ticks(further).unwrap();
            assert_eq!(later, Duration768::from_unsigned_ticks(ticks + 1));
        }
    }

    #[test]
    fn what_each_clock_reaches() {
        let ms = Resolution::MILLISECOND_METRE;
        let ns = Resolution::NANOSECOND_MILLIMETRE;
        let ps = Resolution::new(1_000_000_000_000, 1_000_000).unwrap();
        let near = |got: f64, want: f64| (got - want).abs() < 0.05;
        // 2^384 ms = 1.2e105 yr; at a nanosecond 1.2e99; 2^512 ns = 4.2e137; 2^768 ns = 4.9e214; 2^1024 ns = 5.7e291
        assert!(near(Clock384::new(ms).reach_log10_years(), 105.1));
        assert!(near(Clock384::new(ns).reach_log10_years(), 99.1));
        assert!(near(Clock512::new(ns).reach_log10_years(), 137.6));
        assert!(near(Clock768::new(ns).reach_log10_years(), 214.7));
        assert!(near(Clock768::new(ps).reach_log10_years(), 211.7));
        assert!(near(Clock1024::new(ns).reach_log10_years(), 291.8));
        assert_eq!(
            (
                Clock384::BITS,
                Clock512::BITS,
                Clock768::BITS,
                Clock1024::BITS
            ),
            (384, 512, 768, 1024)
        );
        assert_eq!(Clock768::millisecond_metre().resolution(), ms);
        // 768 bits hold 1e200 years of picoseconds exactly; 512 bits do not
        let epoch = format!("31557600{}", "0".repeat(212));
        assert!(Timestamp::<8>::from_decimal(&epoch).is_err());
        let t = Time768::from_decimal(&epoch).unwrap();
        let clock = Clock768::new(ps);
        let across = GridVec::new(8_800_000_000_000_000_000_000_000_000_000, 3, -7); // 8.8e24 m in micrometres
        let heard = clock.first_arrival(t, across).unwrap();
        assert!(clock.arrived(t, heard, across));
        assert!(!clock.arrived(
            t,
            heard
                .checked_sub(Duration768::from_unsigned_ticks(1))
                .unwrap(),
            across
        ));
        assert_eq!(Time768::MAX.to_string().len(), 232);
    }

    #[test]
    fn seconds_and_phase_follow_the_tick() {
        let ns = Clock768::new(Resolution::NANOSECOND_MILLIMETRE);
        let t0 = Time768::from_decimal(&format!("7{}", "0".repeat(150))).unwrap();
        let t1 = t0.checked_add(ns.seconds(90).unwrap()).unwrap();
        assert_eq!(ns.local_seconds_between(t1, t0).unwrap(), 90.0);
        assert_eq!(ns.local_seconds_between(t0, t1).unwrap(), -90.0);
        assert_eq!(ns.seconds(90).unwrap().to_string(), "90000000000");
        // a Julian year of nanoseconds is 3.2e16 ticks: past 2^53, so its phase is an exact integer first
        let year = 31_557_600 * 1_000_000_000_u64;
        let phase = t1.modulo_ticks(year).unwrap();
        assert_eq!(
            t1.checked_sub(Duration768::from_unsigned_ticks(phase))
                .unwrap()
                .modulo_ticks(year)
                .unwrap(),
            0
        );
        assert_eq!(
            ns.phase_seconds(Time768::from_ticks(1_500_000_000), 1_000_000_000)
                .unwrap(),
            0.5
        );
        // the default methods read a tick as a millisecond
        let ms = Clock768::millisecond_metre();
        let d = Duration768::from_ticks(-2_500);
        assert_eq!(ms.local_seconds(d).unwrap(), d.local_seconds().unwrap());
        assert_eq!(ms.local_seconds(d).unwrap(), -2.5);
        // beyond 2^53 ticks a float cannot hold the count
        let long = Duration768::from_unsigned_ticks(1 << 60);
        assert_eq!(ns.local_seconds(long), Err(Error::PrecisionLoss));
    }

    #[test]
    fn grid_from_metres_and_from_text() {
        let mm = Clock512::new(Resolution::NANOSECOND_MILLIMETRE);
        assert_eq!(
            mm.grid(1.0, -0.0004, 2.5).unwrap(),
            GridVec::new(1000, 0, 2500)
        );
        assert_eq!(
            mm.grid(0.0005, -0.0005, 0.0).unwrap(),
            GridVec::new(1, -1, 0)
        ); // ties away from zero
        assert!(mm.grid(f64::NAN, 0.0, 0.0).is_err() && mm.grid(1e40, 0.0, 0.0).is_err());
        assert_eq!(
            GridVec::from_decimal("40175000000000000000,-3,77").unwrap(),
            GridVec::new(40_175_000_000_000_000_000, -3, 77)
        );
        assert_eq!(GridVec::from_decimal("0,0,0").unwrap(), GridVec::ZERO);
        for bad in [
            "", "1,2", "1,2,3,4", "1, 2,3", "+1,2,3", "1,,3", "1,-,3", "1.5,2,3", "a,b,c",
        ] {
            assert_eq!(
                GridVec::from_decimal(bad),
                Err(Error::InvalidInput),
                "{bad:?}"
            );
        }
        let too_big = format!("{},0,0", "9".repeat(40));
        assert_eq!(GridVec::from_decimal(&too_big), Err(Error::Overflow));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(160))]
        /// At any resolution and any displacement, every width agrees and the first tick is the first that
        /// passes.
        #[test]
        fn widths_agree_at_any_resolution(pick in 0usize..7, x in -1.0e15f64..1.0e15, y in -1.0e15f64..1.0e15, z in -1.0e12f64..1.0e12, lead in 0u64..1_000_000) {
            let (ticks, grid) = RESOLUTIONS[pick];
            let res = Resolution::new(ticks, grid).unwrap();
            let dx = GridVec::from_metres(x, y, z).unwrap();
            let epoch = format!("{}{}", lead + 1, "0".repeat(90));
            let a = first_tick_passes::<6>(res, &epoch, dx);
            prop_assert_eq!(&a, &first_tick_passes::<8>(res, &epoch, dx));
            prop_assert_eq!(&a, &first_tick_passes::<12>(res, &epoch, dx));
            prop_assert_eq!(&a, &first_tick_passes::<16>(res, &epoch, dx));
        }

        /// A moving receiver's arrival, settled from an estimate that is off by a few ticks, is the exact one.
        #[test]
        fn settle_at_any_resolution(pick in 0usize..7, x in 1.0e3f64..1.0e13, off in -6i64..7) {
            let (ticks, grid) = RESOLUTIONS[pick];
            let clock = Clock768::new(Resolution::new(ticks, grid).unwrap());
            let dx = GridVec::from_metres(x, x / 3.0, -x / 7.0).unwrap();
            let sent = Time768::from_decimal(&format!("5{}", "0".repeat(200))).unwrap();
            let exact = clock.first_arrival(sent, dx).unwrap();
            let estimate = if off >= 0 { exact.checked_add(Duration768::from_ticks(off)).unwrap() } else { exact.checked_sub(Duration768::from_ticks(-off)).unwrap_or(sent).max(sent) };
            prop_assert_eq!(clock.settle_arrival(sent, estimate, 16, |_| Ok(dx)).unwrap(), exact);
        }
    }
}
