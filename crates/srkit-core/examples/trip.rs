use srkit_core::{
    kinematics::{Event, Vec3, C, G0, JULIAN_YEAR},
    lightcone::retarded_state,
    segment::Segment,
    time::Time512,
    trip::FlipBurn,
};

fn main() -> Result<(), srkit_core::Error> {
    let trip = FlipBurn::new(4.2465 * C * JULIAN_YEAR, G0)?;
    println!(
        "Proxima rest-to-rest at 1 g: {:.3} coordinate years, {:.3} proper years",
        trip.coordinate_time / JULIAN_YEAR,
        trip.proper_time / JULIAN_YEAR
    );

    // A source stationary one light-second away is seen one second in the past.
    let source = Segment::inertial(Event::new(0.0, Vec3::new(0.0, 0.0, C))?, Vec3::ZERO, 10.0)?;
    let observed = retarded_state(&source, Event::new(5.0, Vec3::ZERO)?)?;
    println!("Source emission time: {:.3} s", observed.event.t);

    // Exact arrival at a huge epoch; this restricted predicate takes an integer
    // stationary distance in metres, not a floating-point moving worldline.
    let epoch = Time512::from_decimal("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000")?;
    let arrival = epoch.stationary_arrival_metres(bnum::BUint::from_digit(C as u64))?;
    println!(
        "Exact stationary signal delay: {} ms",
        arrival.duration_since(epoch)?
    );
    Ok(())
}
