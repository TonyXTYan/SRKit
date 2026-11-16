//! Itineraries: runs, brakes and orbit entries join continuously and end where they were asked to.
use proptest::prelude::*;
use srkit_core::{
    plan::{self, EntryOptions, Phase, Start},
    reach,
    segment::Segment,
    trip::{FlipBurn, Trip},
    worldline::WorldlineBuf,
    Vec3,
};

fn worldline(pieces: &plan::Pieces, tail: Option<Segment>) -> WorldlineBuf {
    let first = pieces.get(0).unwrap();
    let mut segs: Vec<Segment> = pieces.iter().map(|p| p.seg).collect();
    segs.extend(tail);
    WorldlineBuf::new(segs, first.seg.t0(), first.tau0).unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(300))]

    /// A straight run reaches its target, joins continuously and spends the rapidity it reports.
    #[test]
    fn line_joins_and_arrives(
        target in (-50.0..50.0f64, -50.0..50.0f64, -5.0..5.0f64),
        a in 0.05..3.0f64, brake_ratio in 0.3..3.0f64,
        cap in prop_oneof![Just(f64::INFINITY), 0.2..3.0f64],
        phi0 in prop_oneof![Just(0.0), 0.0..0.15f64], phi1 in prop_oneof![Just(0.0), 0.0..0.15f64],
    ) {
        let x1 = Vec3::new(target.0, target.1, target.2);
        prop_assume!(x1.norm() > 0.5);
        let start = Start { t: 2000.0, tau: 7.0, x: Vec3::ZERO };
        let Some(line) = plan::line(start, x1, a, cap, phi0, phi1, a * brake_ratio, None) else {
            // refused: only when the run is too short for the speeds asked
            prop_assert!(phi0 > 0.0 || phi1 > 0.0);
            return Ok(());
        };
        prop_assert!(!line.pieces.is_empty());
        let w = worldline(&line.pieces, None);
        w.view().check_continuity(1e-9 * (1.0 + x1.norm()), 1e-9).unwrap();
        let last = line.pieces.get(line.pieces.len() - 1).unwrap().end();
        prop_assert!(last.x.dist(x1) <= 1e-9 * (1.0 + x1.norm()), "ends {:e} from the target", last.x.dist(x1));
        prop_assert_eq!(line.end.x, x1);
        prop_assert!((last.phi - phi1).abs() < 1e-9);
        prop_assert!(line.phi_peak <= cap + 1e-12 && line.phi_peak >= phi0.max(phi1) - 1e-12);
        // proper time adds up, and the phases come in order
        let total: f64 = line.pieces.iter().map(|p| p.dtau).sum();
        prop_assert!((line.end.tau - 7.0 - total).abs() < 1e-9 * (1.0 + total));
        let phases: Vec<Phase> = line.pieces.iter().map(|p| p.phase).collect();
        prop_assert!(phases.windows(2).all(|w| w[0] != w[1]));
        prop_assert!(!phases.contains(&Phase::Coast) || line.coast_length > 0.0);
    }

    /// From rest to rest, a line is the flip-and-burn (or the capped leg) of the trip planner.
    #[test]
    fn line_agrees_with_the_trip_planner(d in 0.01..1e5f64, a in 0.05..3.0f64, cap in prop_oneof![Just(f64::INFINITY), 1.05..50.0f64]) {
        let start = Start { t: 0.0, tau: 0.0, x: Vec3::ZERO };
        let line = plan::line(start, Vec3::new(d, 0.0, 0.0), a, cap.acosh(), 0.0, 0.0, a, None).unwrap();
        let trip = Trip::new(d, a, cap, 0, false).unwrap();
        prop_assert!((line.end.tau / trip.tau - 1.0).abs() < 1e-9, "{} vs {}", line.end.tau, trip.tau);
        prop_assert!((line.end.t / trip.t_end - 1.0).abs() < 1e-9);
        prop_assert!((line.dphi - 2.0 * trip.leg.phc).abs() < 1e-9 * (1.0 + line.dphi));
        if cap.is_infinite() {
            let fb = FlipBurn::new(d, a).unwrap();
            prop_assert!((trip.tau / fb.tau - 1.0).abs() < 1e-13 && (trip.t_end / fb.t - 1.0).abs() < 1e-13);
            prop_assert!((reach::flip_burn_tau(d, a) / fb.tau - 1.0).abs() < 1e-13);
            prop_assert!((reach::flip_burn_distance(fb.tau, a) / d - 1.0).abs() < 1e-9);
        }
    }

    /// Entering an orbit: the velocity is continuous where the run meets the circle.
    #[test]
    fn orbit_entry_is_smooth(from in (-30.0..30.0f64, -30.0..30.0f64, -3.0..3.0f64), r in 0.5..8.0f64, beta in 0.05..0.9f64, a in 0.2..3.0f64) {
        let start = Start { t: 100.0, tau: 0.0, x: Vec3::new(from.0, from.1, from.2) };
        let centre = Vec3::new(1.0, -2.0, 0.5);
        let Some(entry) = plan::orbit_entry(start, centre, r, beta, a, f64::INFINITY, EntryOptions::default()) else { return Ok(()); };
        let circle = plan::circ_piece(entry.end.start(), &entry.orb, f64::INFINITY);
        let w = worldline(&entry.pieces, Some(circle.seg));
        w.view().check_continuity(1e-8 * (1.0 + start.x.norm()), 1e-8).unwrap();
        // on the circle: at the radius, at the orbit's speed, counter-clockwise in (e1, e2)
        let s = circle.seg.at(circle.seg.t0() + 0.3);
        prop_assert!((s.x.dist(centre) - r).abs() < 1e-9 * (1.0 + r));
        prop_assert!((s.v.norm() - beta).abs() < 1e-12);
        let turn = (s.x - centre).cross(s.v).dot(entry.orb.e1.cross(entry.orb.e2));
        prop_assert!(turn > 0.0);
        // and leaving it brakes to rest
        let (piece, end, dphi) = plan::orbit_exit(Start { t: 0.0, tau: 0.0, x: s.x }, s.v, a);
        prop_assert!(piece.is_some() && end.v.norm() < 1e-12 && (dphi - beta.atanh()).abs() < 1e-12);
    }
}

#[test]
fn short_hops_keep_their_digits() {
    // One astronomical unit at 1 g: the peak rapidity is 4e-3, where acosh(1 + x) would lose half its digits.
    let (au, g) = (1.5812507409819728e-5, 1.0322952409169414);
    let start = Start {
        t: 0.0,
        tau: 0.0,
        x: Vec3::ZERO,
    };
    let line = plan::line(
        start,
        Vec3::new(au, 0.0, 0.0),
        g,
        f64::INFINITY,
        0.0,
        0.0,
        g,
        None,
    )
    .unwrap();
    // Newtonian limit: tau = 2 sqrt(d / a), corrected by 1 - a d / 24 to leading order
    let newton = 2.0 * (au / g).sqrt();
    let expected = newton * (1.0 - g * au / 24.0);
    assert!(
        (line.end.tau / expected - 1.0).abs() < 1e-10,
        "{} vs {expected}",
        line.end.tau
    );
    let last = line.pieces.get(1).unwrap().end();
    assert!(
        last.x.dist(Vec3::new(au, 0.0, 0.0)) < 1e-15 * au.max(1e-3),
        "missed by {:e} ly",
        last.x.dist(Vec3::new(au, 0.0, 0.0))
    );
}

#[test]
fn a_run_of_no_length_is_empty_and_a_jump_is_light_like() {
    let start = Start {
        t: 5.0,
        tau: 1.0,
        x: Vec3::new(1.0, 2.0, 3.0),
    };
    let line = plan::line(start, start.x, 1.0, f64::INFINITY, 0.0, 0.0, 1.0, None).unwrap();
    assert!(line.pieces.is_empty() && line.dphi == 0.0 && line.end.x == start.x);
    let jump = plan::null_piece(start, Vec3::new(1.0, 2.0, 7.0));
    assert_eq!(jump.dtau, 0.0);
    assert_eq!(jump.seg.t1(), 9.0);
    assert_eq!(jump.end().x, Vec3::new(1.0, 2.0, 7.0));
    let rest = plan::rest(start, 2.5);
    assert_eq!(rest.end().tau, 3.5);
    let (none, end, spent) = plan::brake(start, Vec3::X, 0.0, 1.0);
    assert!(none.is_none() && spent == 0.0 && end.x == start.x);
}
