//! A 1 g ride to Proxima Centauri: the plan, the worldline, what home sees, and an exact arrival tick.
//!
//!     cargo run --example trip
use srkit_core::{
    causal::{self, GridVec},
    kinematics::Event,
    plan::{self, Start},
    segment::Segment,
    time::Time512,
    units::{Units, G0_PER_YEAR, LIGHT_YEAR},
    worldline::WorldlineBuf,
    Vec3,
};

fn main() -> Result<(), srkit_core::Error> {
    // Julian years and light-years; 1 g is 1.0323 per year
    let proxima = Vec3::new(4.2465, 0.0, 0.0);
    let start = Start {
        t: 2100.0,
        tau: 0.0,
        x: Vec3::ZERO,
    };
    let ride = plan::line(
        start,
        proxima,
        G0_PER_YEAR,
        f64::INFINITY,
        0.0,
        0.0,
        G0_PER_YEAR,
        None,
    )
    .expect("a rest-to-rest run always has a plan");
    println!(
        "Proxima at 1 g: {:.3} yr on Earth's calendar, {:.3} yr aboard, peak gamma {:.3}, photon-rocket mass ratio {:.1}",
        ride.end.t - start.t,
        ride.end.tau,
        ride.phi_peak.cosh(),
        srkit_core::rocket::mass_ratio(ride.dphi, 1.0)
    );

    // the worldline: wait, boost, brake, stay
    let mut segments = vec![Segment::rest(2000.0, start.t, Vec3::ZERO)];
    segments.extend(ride.pieces.iter().map(|p| p.seg));
    segments.push(Segment::rest(ride.end.t, f64::INFINITY, proxima));
    let ship = WorldlineBuf::new(segments, start.t, 0.0)?;
    let ship = ship.view();
    ship.check_continuity(1e-12, 1e-12)?;

    // what an observer on Earth sees of the ship, year by year: a past-light-cone query, never the ship's "now"
    for year in [2101.0, 2103.0, 2105.0, 2107.0, 2110.0, 2111.0] {
        let (t_emitted, seen) = ship
            .seen_from(Event::new(year, Vec3::ZERO))
            .expect("the ship is always visible from home");
        println!(
            "  Earth {year}: sees the ship {:.3} ly out at {:.1}% of c, its clock reading {:.3} yr (light left in {:.2})",
            seen.x.x,
            100.0 * seen.v.norm(),
            ship.tau(t_emitted),
            t_emitted
        );
    }

    // "we have arrived", sent on landing: when does Earth hear it?
    let landing = Event::new(ride.end.t, proxima);
    let home = WorldlineBuf::new(
        [Segment::rest(2000.0, f64::INFINITY, Vec3::ZERO)],
        2000.0,
        0.0,
    )?;
    println!(
        "  the landing message reaches Earth in {:.3}",
        home.view().recv(landing)
    );

    // the same delay on the engine's exact clock, at an epoch no double could hold: whole milliseconds, never a tick early
    let epoch = Time512::from_decimal(&format!("3{}", "0".repeat(110)))?; // about 1e100 years
    let metres = Units::YEARS.length_to_si(4.2465);
    let distance = GridVec::from_metres(metres, 0.0, 0.0)?;
    let arrival = causal::first_arrival(epoch, distance)?;
    println!(
        "  exact: {} m ({:.4} ly) is {} ms of light time; arrived a tick earlier? {}",
        distance.x,
        distance.x as f64 / LIGHT_YEAR,
        arrival.signed_duration_since(epoch),
        causal::arrived(
            epoch,
            arrival.checked_sub(srkit_core::time::Duration512::from_unsigned_millis(1))?,
            distance
        )
    );
    Ok(())
}
