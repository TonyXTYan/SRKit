use bnum::BUint;
use proptest::prelude::*;
use srkit_core::{
    time::{Duration384, Duration512, Time384, Time512, MAX_LOCAL_MILLIS},
    Error,
};

#[test]
fn both_widths_preserve_huge_epochs_and_millisecond_differences() {
    let epoch = "10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    let t384 = Time384::from_decimal(epoch).unwrap();
    let t512 = Time512::from_decimal(epoch).unwrap();
    assert_eq!(t384.to_string(), epoch);
    assert_eq!(t512.to_string(), epoch);
    let later384 = t384.checked_add(Duration384::from_millis(1)).unwrap();
    let later512 = t512.checked_add(Duration512::from_millis(1)).unwrap();
    assert_eq!(later384.local_seconds_since(t384), Ok(0.001));
    assert_eq!(later512.local_seconds_since(t512), Ok(0.001));
    assert_eq!(t384.local_seconds_since(later384), Ok(-0.001));
    assert_eq!(t512.local_seconds_since(later512), Ok(-0.001));
}

#[test]
fn timestamp_arithmetic_rejects_wrapping() {
    assert_eq!(
        Time384::MAX.checked_add(Duration384::from_millis(1)),
        Err(Error::Overflow)
    );
    assert_eq!(
        Time512::MAX.checked_add(Duration512::from_millis(1)),
        Err(Error::Overflow)
    );
    assert_eq!(
        Time384::ZERO.checked_sub(Duration384::from_millis(1)),
        Err(Error::Overflow)
    );
    assert_eq!(
        Time512::ZERO.checked_add(Duration512::from_millis(-1)),
        Err(Error::Overflow)
    );
    let max_difference = Time512::MAX.signed_duration_since(Time512::ZERO);
    assert_eq!(Time512::ZERO.checked_add(max_difference), Ok(Time512::MAX));
    assert_eq!(
        max_difference.checked_add(Duration512::from_millis(1)),
        Err(Error::Overflow)
    );
}

#[test]
fn widths_remain_distinct_and_full_width_durations_round_trip() {
    let wider_epoch = format!("1{}", "0".repeat(120));
    assert_eq!(Time384::from_decimal(&wider_epoch), Err(Error::Overflow));
    let epoch512 = Time512::from_decimal(&wider_epoch).unwrap();
    let negative = Time512::ZERO.signed_duration_since(epoch512);
    assert!(negative.is_negative());
    assert_eq!(negative.to_string(), format!("-{wider_epoch}"));
    assert_eq!(epoch512.checked_add(negative), Ok(Time512::ZERO));
}

#[test]
fn signed_durations_normalize_zero_and_support_i64_min() {
    assert_eq!(Duration384::from_decimal("-0").unwrap(), Duration384::ZERO);
    assert!(!Duration384::ZERO.negated().is_negative());
    assert_eq!(
        Duration384::from_millis(i64::MIN).to_string(),
        i64::MIN.to_string()
    );
    let neg = Duration384::from_millis(-11);
    let pos = Duration384::from_millis(7);
    assert_eq!(neg.checked_add(pos).unwrap().to_string(), "-4");
    assert_eq!(pos.checked_add(neg).unwrap().to_string(), "-4");
    assert_eq!(neg.checked_sub(pos).unwrap().to_string(), "-18");
    assert_eq!(pos.checked_sub(neg).unwrap().to_string(), "18");
    assert_eq!(neg.checked_add(neg.negated()), Ok(Duration384::ZERO));
}

#[test]
fn invalid_decimal_and_overflow_are_distinct() {
    for text in ["", "-1", "+1", " 1", "1 ", "1.0", "1_000", "１２"] {
        assert_eq!(Time384::from_decimal(text), Err(Error::InvalidInput));
    }
    for text in ["-", "--1", "+1", "- 1"] {
        assert_eq!(Duration512::from_decimal(text), Err(Error::InvalidInput));
    }
    let too_big = "9".repeat(160);
    assert_eq!(Time384::from_decimal(&too_big), Err(Error::Overflow));
    assert_eq!(Time512::from_decimal(&too_big), Err(Error::Overflow));
    assert_eq!("00012".parse::<Time512>().unwrap().to_string(), "12");
}

#[test]
fn local_conversion_has_an_explicit_precision_limit() {
    let edge = Time512::from_millis(MAX_LOCAL_MILLIS);
    assert!(edge.local_seconds_since(Time512::ZERO).is_ok());
    let outside = edge.checked_add(Duration512::from_millis(1)).unwrap();
    assert_eq!(
        outside.local_seconds_since(Time512::ZERO),
        Err(Error::PrecisionLoss)
    );
    assert_eq!(
        Time512::ZERO.local_seconds_since(outside),
        Err(Error::PrecisionLoss)
    );
    assert_eq!(
        Time512::ZERO.duration_since(Time512::from_millis(1)),
        Err(Error::OutOfRange)
    );
}

#[test]
fn periodic_phase_reduces_exactly_before_conversion() {
    // 10^3 = -1 mod 1001, so 10^100 = -10 = 991 mod 1001.
    // This epoch cannot be represented in f64 exactly.
    let epoch = format!("1{}", "0".repeat(100));
    let timestamp = Time512::from_decimal(&epoch).unwrap();
    assert_eq!(timestamp.modulo_millis(1001), Ok(991));
    assert_eq!(timestamp.phase_seconds(1001), Ok(0.991));
    assert_eq!(timestamp.modulo_millis(0), Err(Error::InvalidInput));
    assert_eq!(Time512::MAX.modulo_millis(1), Ok(0));
}

#[test]
fn stationary_signals_round_up_and_never_arrive_a_tick_early() {
    let epoch = Time384::from_decimal(&format!("1{}", "0".repeat(100))).unwrap();
    // 299,792 metres is slightly less than one millisecond of light travel.
    // 299,793 metres is slightly more and must take two ticks.
    for (distance, expected_millis) in [
        (0, 0),
        (1, 1),
        (299_792, 1),
        (299_793, 2),
        (299_792_458, 1000),
    ] {
        let metres = BUint::<6>::from_digit(distance);
        let arrival = epoch.stationary_arrival_metres(metres).unwrap();
        assert_eq!(
            arrival.signed_duration_since(epoch),
            Duration384::from_millis(expected_millis)
        );
        assert_eq!(arrival.stationary_signal_arrived(epoch, metres), Ok(true));
        if expected_millis > 0 {
            let before = arrival.checked_sub(Duration384::from_millis(1)).unwrap();
            assert_eq!(before.stationary_signal_arrived(epoch, metres), Ok(false));
        }
    }
    assert_eq!(
        Time384::ZERO.stationary_signal_arrived(epoch, BUint::<6>::ZERO),
        Ok(false)
    );
}

#[test]
fn full_width_stationary_distance_avoids_multiplication_overflow() {
    let distance = BUint::<8>::MAX;
    let delay = Duration512::light_travel_time_metres(distance).unwrap();
    assert!(!delay.is_negative());
    assert!(delay.magnitude_millis() < distance);
    // Validate the direct causal inequality independently in a wider test
    // integer. The production path stays within the selected clock width.
    let distance_wide = BUint::<9>::from_str_radix(&distance.to_string(), 10).unwrap();
    let delay_wide = BUint::<9>::from_str_radix(&delay.magnitude_millis().to_string(), 10).unwrap();
    let c = BUint::<9>::from_digit(299_792_458);
    let scaled_distance = distance_wide * BUint::from_digit(1000);
    assert!(delay_wide * c >= scaled_distance);
    assert!((delay_wide - BUint::ONE) * c < scaled_distance);
    assert_eq!(
        Time512::MAX.stationary_arrival_metres(BUint::from_digit(1)),
        Err(Error::Overflow)
    );
    assert_eq!(
        Time512::MAX.stationary_signal_arrived(Time512::MAX, BUint::from_digit(1)),
        Ok(false)
    );
}

proptest! {
    #[test]
    fn local_round_trip_never_loses_a_tick(epoch in any::<u64>(), delta in -1000000_i64..1000000) {
        let base = Time512::from_millis(epoch).checked_add(Duration512::from_millis(1000000)).unwrap();
        let later = base.checked_add(Duration512::from_millis(delta)).unwrap();
        prop_assert_eq!(later.signed_duration_since(base), Duration512::from_millis(delta));
        prop_assert_eq!(later.checked_sub(Duration512::from_millis(delta)), Ok(base));
    }

    #[test]
    fn phase_agrees_with_integer_remainder(epoch in any::<u64>(), period in 1_u64..u64::MAX) {
        prop_assert_eq!(Time384::from_millis(epoch).modulo_millis(period), Ok(epoch % period));
        prop_assert_eq!(Time512::from_millis(epoch).modulo_millis(period), Ok(epoch % period));
    }

    #[test]
    fn causal_delay_is_the_first_valid_integer_tick(distance_metres in any::<u64>()) {
        let delay = Duration512::light_travel_time_metres(BUint::<8>::from_digit(distance_metres)).unwrap();
        let ticks = delay.magnitude_millis().digits()[0] as u128;
        let distance_times_1000 = distance_metres as u128 * 1000;
        prop_assert!(ticks * 299_792_458 >= distance_times_1000);
        if ticks > 0 {
            prop_assert!((ticks - 1) * 299_792_458 < distance_times_1000);
        }
    }
}
