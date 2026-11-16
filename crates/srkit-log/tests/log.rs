use proptest::prelude::*;
use srkit_core::{
    causal,
    kinematics::Event,
    plan::{self, Start},
    segment::Segment,
    time::{Duration384, Time384},
    units::Units,
    worldline::WorldlineBuf,
    Vec3,
};
use srkit_log::{exact, EventLog, Inbox};

fn at_rest(x: Vec3) -> WorldlineBuf {
    WorldlineBuf::new([Segment::rest(0.0, f64::INFINITY, x)], 0.0, 0.0).unwrap()
}

/// A ship that leaves the origin at t = 10 for a point 8 away at a = 1, and stays.
fn traveller() -> WorldlineBuf {
    let start = Start {
        t: 10.0,
        tau: 10.0,
        x: Vec3::ZERO,
    };
    let line = plan::line(
        start,
        Vec3::new(8.0, 0.0, 0.0),
        1.0,
        f64::INFINITY,
        0.0,
        0.0,
        1.0,
        None,
    )
    .unwrap();
    let mut segs = vec![Segment::rest(0.0, 10.0, Vec3::ZERO)];
    segs.extend(line.pieces.iter().map(|p| p.seg));
    segs.push(Segment::rest(line.end.t, f64::INFINITY, line.end.x));
    WorldlineBuf::new(segs, 0.0, 0.0).unwrap()
}

#[test]
fn news_arrives_in_light_cone_order_not_log_order() {
    let mut log = EventLog::new();
    let far = log.append(Event::new(1.0, Vec3::new(10.0, 0.0, 0.0)), "far, early");
    let near = log.append(Event::new(5.0, Vec3::new(1.0, 0.0, 0.0)), "near, late");
    let here = log.append(Event::new(7.0, Vec3::ZERO), "here");
    let home = at_rest(Vec3::ZERO);
    let order: Vec<_> = log
        .known_by(&home.view(), 100.0)
        .iter()
        .map(|a| a.id)
        .collect();
    assert_eq!(order, [near, here, far]);
    // nothing is known before its light arrives, and it is known from the instant it does
    assert!(log.known_by(&home.view(), 5.999).is_empty());
    assert_eq!(log.known_by(&home.view(), 6.0).len(), 1);
    assert_eq!(log.arrival(far, &home.view()), Some(11.0));
    assert_eq!(log.known_by(&home.view(), 10.999).len(), 2);
    let a = log.known_by(&home.view(), 11.0)[2];
    assert_eq!(
        (a.t_emitted, a.t_arrival, a.light_age(), a.tau_arrival),
        (1.0, 11.0, 10.0, 11.0)
    );
}

#[test]
fn an_inbox_delivers_each_event_once_and_catches_up_in_one_call() {
    let mut log = EventLog::new();
    let home = at_rest(Vec3::ZERO);
    let mut tick_by_tick = Inbox::new();
    let mut away = Inbox::new();
    let mut seen = Vec::new();
    for step in 0..200 {
        let t = step as f64 * 0.5;
        if step % 3 == 0 {
            // events appear all over, as the world goes on
            let x = Vec3::new(
                ((step * 7) % 23) as f64 - 11.0,
                ((step * 5) % 13) as f64 - 6.0,
                0.0,
            );
            log.append(Event::new(t, x), step);
        }
        seen.extend_from_slice(tick_by_tick.read(&log, &home.view(), t));
        assert_eq!(tick_by_tick.read(&log, &home.view(), t), &[]); // reading again gives nothing new
    }
    // one read after a long absence returns the same news in the same order
    let all = away.read(&log, &home.view(), 99.5).to_vec();
    assert_eq!(all, seen);
    assert_eq!(away.history(), &seen[..]);
    assert!(all.windows(2).all(|w| w[0].t_arrival <= w[1].t_arrival));
    let ids: std::collections::BTreeSet<_> = all.iter().map(|a| a.id).collect();
    assert_eq!(ids.len(), all.len());
    // reading backwards in time gives nothing and does not rewind
    assert_eq!(away.read(&log, &home.view(), 3.0), &[]);
    assert_eq!(away.read_until(), 99.5);
    // what is still on its way is next
    let next = away.peek(&log, &home.view());
    assert!(next.is_none_or(|a| a.t_arrival > 99.5));
}

#[test]
fn a_moving_seat_gets_the_news_later_or_never() {
    let mut log = EventLog::new();
    let ship = traveller();
    let w = ship.view();
    // home broadcasts every unit of time; the ship flees, then stops
    for k in 0..40 {
        log.append(Event::new(k as f64, Vec3::ZERO), k);
    }
    let arrivals = log.known_by(&w, 1e6);
    assert_eq!(arrivals.len(), 40);
    // before launch it hears each broadcast at once; in flight they pile up behind it; at rest again, 8 late
    assert_eq!(arrivals[5].t_arrival, 5.0);
    let last = arrivals.last().unwrap();
    assert!((last.t_arrival - (39.0 + 8.0)).abs() < 1e-9);
    assert!(arrivals.windows(2).all(|p| p[1].t_arrival > p[0].t_arrival));
    // the seat's own clock at each arrival runs slow in flight
    let stretched = arrivals
        .iter()
        .filter(|a| a.t_arrival - a.tau_arrival > 1e-6)
        .count();
    assert!(stretched > 10);

    // a ship that accelerates for ever stops hearing home: nothing sent after 1/a arrives
    let fleeing = WorldlineBuf::new(
        [Segment::hyperbolic(
            0.0,
            f64::INFINITY,
            Vec3::ZERO,
            Vec3::X,
            1.0,
            0.0,
        )],
        0.0,
        0.0,
    )
    .unwrap();
    let heard = log.known_by(&fleeing.view(), f64::MAX);
    assert_eq!(
        heard.len(),
        1,
        "only the broadcast at t = 0 is inside the horizon"
    );
}

#[test]
fn a_changed_flight_plan_keeps_what_was_heard() {
    let mut log = EventLog::new();
    for k in 0..30 {
        log.append(Event::new(k as f64, Vec3::new(20.0, 0.0, 0.0)), k);
    }
    let staying = at_rest(Vec3::ZERO);
    let mut inbox = Inbox::new();
    let heard_before = inbox.read(&log, &staying.view(), 25.0).len();
    assert_eq!(heard_before, 6); // sent 0..=5, twenty away
                                 // the seat now flies towards the source: the rest arrives sooner, and nothing is delivered twice
    let start = Start {
        t: 25.0,
        tau: 25.0,
        x: Vec3::ZERO,
    };
    let line = plan::line(
        start,
        Vec3::new(15.0, 0.0, 0.0),
        1.0,
        f64::INFINITY,
        0.0,
        0.0,
        1.0,
        None,
    )
    .unwrap();
    let mut segs = vec![Segment::rest(0.0, 25.0, Vec3::ZERO)];
    segs.extend(line.pieces.iter().map(|p| p.seg));
    segs.push(Segment::rest(line.end.t, f64::INFINITY, line.end.x));
    let flying = WorldlineBuf::new(segs, 0.0, 0.0).unwrap();
    inbox.replan(&log, &flying.view());
    let rest = inbox.read(&log, &flying.view(), 1e6).to_vec();
    assert_eq!(heard_before + rest.len(), 30);
    assert!(rest[0].t_arrival > 25.0);
    assert!(
        rest.last().unwrap().t_arrival < 29.0 + 20.0,
        "closer to the source, the last broadcast arrives sooner"
    );
    let ids: std::collections::BTreeSet<_> = inbox.history().iter().map(|a| a.id).collect();
    assert_eq!(ids.len(), 30);
}

#[test]
fn a_saved_log_replays_to_the_same_history() {
    let mut log = EventLog::new();
    log.append(
        Event::new(0.1 + 0.2, Vec3::new(1.0 / 3.0, -2.5e-17, 1e300)),
        String::from("a b c"),
    );
    log.append(Event::new(-7.0, Vec3::ZERO), String::new());
    let text = log.to_text(|p| p.clone());
    let back = EventLog::from_text(&text, |s| Some(s.to_string())).unwrap();
    assert_eq!(back, log);
    assert!(EventLog::<String>::from_text("1.0 2.0 nonsense", |s| Some(s.to_string())).is_err());
    assert!(EventLog::new()
        .try_append(Event::new(f64::NAN, Vec3::ZERO), ())
        .is_err());
}

#[test]
fn exact_arrival_is_never_a_tick_early() {
    // a scene in seconds and light-seconds, at an epoch far beyond what a double can hold
    let frame = exact::Frame::<6> {
        epoch: Time384::from_decimal(&format!("1{}", "0".repeat(60))).unwrap(),
        t_epoch: 0.0,
        units: Units::SECONDS,
    };
    let one = Duration384::from_unsigned_millis(1);
    // a receiver at rest 3 light-seconds away: exactly 3000 ticks, and not 2999
    let rx = at_rest(Vec3::new(3.0, 0.0, 0.0));
    let tick = exact::arrival_tick(&frame, Event::new(0.0, Vec3::ZERO), &rx.view(), 8).unwrap();
    assert_eq!(
        tick.signed_duration_since(frame.epoch),
        Duration384::from_unsigned_millis(3000)
    );
    // receivers moving at many speeds: the settled tick passes the exact check and the one before does not
    for k in 0..400 {
        let beta = -0.95 + 1.9 * (k as f64 / 399.0);
        let x0 = 0.5 + (k % 17) as f64 * 1.37;
        let rx = WorldlineBuf::new(
            [Segment::inertial(
                0.0,
                f64::INFINITY,
                Vec3::new(x0, 0.3, 0.0),
                Vec3::new(beta, 0.0, 0.0),
            )],
            0.0,
            0.0,
        )
        .unwrap();
        let event = Event::new(0.25, Vec3::ZERO);
        let tick = exact::arrival_tick(&frame, event, &rx.view(), 8).unwrap();
        let emission = frame.tick_at(event.t).unwrap();
        let dx = |tick| {
            frame
                .grid(rx.view().at(frame.time_at(tick).unwrap()).x)
                .unwrap()
        };
        assert!(causal::arrived(emission, tick, dx(tick)), "beta {beta}");
        let before = tick.checked_sub(one).unwrap();
        assert!(
            !causal::arrived(emission, before, dx(before)),
            "beta {beta}: the tick before had not arrived"
        );
    }
    // behind a horizon there is no arrival to settle
    let fleeing = WorldlineBuf::new(
        [Segment::hyperbolic(
            0.0,
            f64::INFINITY,
            Vec3::ZERO,
            Vec3::X,
            1.0,
            0.0,
        )],
        0.0,
        0.0,
    )
    .unwrap();
    assert_eq!(
        exact::arrival_tick(&frame, Event::new(5.0, Vec3::ZERO), &fleeing.view(), 8),
        Err(srkit_core::Error::NoIntersection)
    );
}

proptest! {
    /// However the same events are appended and however the reads are split, a receiver is given the same history.
    #[test]
    fn history_does_not_depend_on_how_it_is_read(
        events in prop::collection::vec((0.0..50.0f64, -20.0..20.0f64, -20.0..20.0f64), 1..40),
        cuts in prop::collection::vec(0.0..120.0f64, 0..12),
    ) {
        let mut log = EventLog::new();
        for (t, x, y) in &events {
            log.append(Event::new(*t, Vec3::new(*x, *y, 0.0)), ());
        }
        let ship = traveller();
        let w = ship.view();
        let whole = log.known_by(&w, 200.0);
        let mut inbox = Inbox::new();
        let mut pieces = Vec::new();
        let mut cuts = cuts.clone();
        cuts.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for c in cuts {
            pieces.extend_from_slice(inbox.read(&log, &w, c));
        }
        pieces.extend_from_slice(inbox.read(&log, &w, 200.0));
        prop_assert_eq!(&pieces, &whole);
        // causality: every arrival is on or inside the event's future light cone
        for a in &whole {
            let e = log.get(a.id).unwrap().at();
            let r = w.at(a.t_arrival).x.dist(e.x);
            prop_assert!(a.t_arrival - e.t >= r - 1e-9 * (1.0 + r));
        }
    }
}
