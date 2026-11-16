//! # srkit-log
//!
//! An append-only log of events, read through light cones.
//!
//! An event is recorded where and when it happens. Nobody reads the log directly: a receiver (a worldline) asks
//! what has **arrived**, and gets events in the order their light reaches it. That one rule is what makes a
//! simulation causal: nothing can act on news it has not received, and the same log replays to the same history
//! for every observer.
//!
//! * [`EventLog`] holds the events. Appending never changes what any receiver has already been given.
//! * [`Inbox`] is one receiver's view. It works out each event's arrival once, keeps them sorted, and hands over
//!   what is new since it was last read: the catch-up after time away is the same call as a single tick.
//! * [`exact`] settles an arrival on the engine's integer clock, so "has it arrived by tick n?" is never off by
//!   one.
//!
//! Arrival order is total and deterministic: by arrival time, then by emission time, then by the order of
//! appending. Events that never arrive (behind a receiver's horizon) are never delivered.
//!
//! ```
//! use srkit_core::{kinematics::Event, segment::Segment, worldline::WorldlineBuf, Vec3};
//! use srkit_log::{EventLog, Inbox};
//!
//! let mut log = EventLog::new();
//! log.append(Event::new(2030.0, Vec3::new(4.0, 0.0, 0.0)), "colony founded");
//! log.append(Event::new(2031.0, Vec3::new(0.0, 0.0, 0.0)), "election at home");
//!
//! // a seat at rest at the origin
//! let home = WorldlineBuf::new([Segment::rest(2000.0, f64::INFINITY, Vec3::ZERO)], 2000.0, 0.0).unwrap();
//! let mut inbox = Inbox::new();
//! let news: Vec<_> = inbox.read(&log, &home.view(), 2033.0).iter().map(|a| *log.get(a.id).unwrap().payload()).collect();
//! assert_eq!(news, ["election at home"]);           // the colony's light is still a year away
//! let later: Vec<_> = inbox.read(&log, &home.view(), 2034.0).iter().map(|a| a.t_arrival).collect();
//! assert_eq!(later, [2034.0]);                       // founded 2030, 4 light-years off
//! ```
#![forbid(unsafe_code)]

pub mod exact;

use srkit_core::{kinematics::Event, worldline::Worldline};

/// The position of an event in the log: the order it was appended in.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EventId(pub u64);

/// An event with what happened there.
#[derive(Clone, Debug, PartialEq)]
pub struct Logged<T> {
    id: EventId,
    at: Event,
    payload: T,
}

impl<T> Logged<T> {
    pub fn id(&self) -> EventId {
        self.id
    }
    /// Where and when it happened.
    pub fn at(&self) -> Event {
        self.at
    }
    pub fn payload(&self) -> &T {
        &self.payload
    }
}

/// Why an append was refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LogError {
    /// The event's time or position is not a finite number.
    NotFinite,
    /// A line of a saved log could not be read.
    Malformed,
}

/// The events of a world, in the order they were recorded.
#[derive(Clone, Debug, PartialEq)]
pub struct EventLog<T> {
    events: Vec<Logged<T>>,
}

impl<T> Default for EventLog<T> {
    fn default() -> Self {
        Self { events: Vec::new() }
    }
}

impl<T> EventLog<T> {
    pub fn new() -> Self {
        Self::default()
    }
    /// Record an event. Events may be appended in any order of time; an [`Inbox`] sorts them by arrival.
    pub fn append(&mut self, at: Event, payload: T) -> EventId {
        self.try_append(at, payload)
            .expect("event coordinates must be finite")
    }
    /// As [`EventLog::append`], refusing an event that is not at a finite time and place.
    pub fn try_append(&mut self, at: Event, payload: T) -> Result<EventId, LogError> {
        if !at.is_finite() {
            return Err(LogError::NotFinite);
        }
        let id = EventId(self.events.len() as u64);
        self.events.push(Logged { id, at, payload });
        Ok(id)
    }
    pub fn len(&self) -> usize {
        self.events.len()
    }
    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }
    pub fn get(&self, id: EventId) -> Option<&Logged<T>> {
        self.events.get(id.0 as usize)
    }
    /// Every event in the order recorded. This is the archive view: it ignores light cones, so game code that
    /// shows a seat its world must go through an [`Inbox`] instead.
    pub fn archive(&self) -> impl Iterator<Item = &Logged<T>> {
        self.events.iter()
    }
    /// When the light of event `id` reaches `receiver`: a coordinate time, or `+inf` if it never does.
    pub fn arrival(&self, id: EventId, receiver: &Worldline<'_>) -> Option<f64> {
        self.get(id).map(|e| receiver.recv(e.at))
    }
    /// Everything `receiver` has received by coordinate time `t`, in order of arrival.
    pub fn known_by(&self, receiver: &Worldline<'_>, t: f64) -> Vec<Arrival> {
        let mut inbox = Inbox::new();
        inbox.read(self, receiver, t).to_vec()
    }
    /// Save the log as text, one event per line: `t x y z payload`, the numbers written so that they read back
    /// to the same bits. `encode` must not produce a line break.
    pub fn to_text(&self, mut encode: impl FnMut(&T) -> String) -> String {
        let mut out = String::new();
        for e in &self.events {
            let p = encode(&e.payload);
            debug_assert!(!p.contains('\n'));
            out.push_str(&format!(
                "{:?} {:?} {:?} {:?} {}\n",
                e.at.t, e.at.x.x, e.at.x.y, e.at.x.z, p
            ));
        }
        out
    }
    /// Load a log saved by [`EventLog::to_text`].
    pub fn from_text(
        text: &str,
        mut decode: impl FnMut(&str) -> Option<T>,
    ) -> Result<Self, LogError> {
        let mut log = Self::new();
        for line in text.lines().filter(|l| !l.is_empty()) {
            let mut parts = line.splitn(5, ' ');
            let mut number = || {
                parts
                    .next()
                    .and_then(|s| s.parse::<f64>().ok())
                    .ok_or(LogError::Malformed)
            };
            let (t, x, y, z) = (number()?, number()?, number()?, number()?);
            let payload = decode(parts.next().unwrap_or("")).ok_or(LogError::Malformed)?;
            log.try_append(Event::new(t, srkit_core::Vec3::new(x, y, z)), payload)?;
        }
        Ok(log)
    }
}

/// An event's arrival at a receiver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Arrival {
    pub id: EventId,
    /// When the light left: the event's own time.
    pub t_emitted: f64,
    /// When it reaches the receiver, in coordinate time.
    pub t_arrival: f64,
    /// The receiver's proper time then.
    pub tau_arrival: f64,
}

impl Arrival {
    /// How old the news is when it arrives, in coordinate time.
    pub fn light_age(&self) -> f64 {
        self.t_arrival - self.t_emitted
    }
    fn key(&self) -> (f64, f64, u64) {
        (self.t_arrival, self.t_emitted, self.id.0)
    }
}

/// One receiver's view of a log.
///
/// The inbox remembers how much of the log it has examined and how far it has been read, so each call does only
/// the new work. Use one inbox per receiver, always with the same worldline: if the receiver's future changes (a
/// new flight plan), call [`Inbox::replan`] so that arrivals not yet delivered are worked out again.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Inbox {
    /// Arrivals of every event examined so far that will arrive at all, sorted.
    pending: Vec<Arrival>,
    /// How many events of the log have been examined.
    examined: usize,
    /// How many of `pending` have been delivered.
    delivered: usize,
    /// The time read up to.
    read_until: f64,
}

impl Inbox {
    pub fn new() -> Self {
        Self {
            read_until: f64::NEG_INFINITY,
            ..Self::default()
        }
    }
    /// Work out the arrival of every event appended since the last call.
    fn sync<T>(&mut self, log: &EventLog<T>, receiver: &Worldline<'_>) {
        for e in &log.events[self.examined.min(log.events.len())..] {
            let t_arrival = receiver.recv(e.at);
            if t_arrival.is_finite() {
                let a = Arrival {
                    id: e.id,
                    t_emitted: e.at.t,
                    t_arrival,
                    tau_arrival: receiver.tau(t_arrival),
                };
                // Keep delivered arrivals where they are: news recorded late that would have arrived earlier is
                // delivered next, never slipped into the past the receiver has already been told.
                let from = self.delivered;
                let offset = self.pending[from..].partition_point(|p| p.key() <= a.key());
                self.pending.insert(from + offset, a);
            }
        }
        self.examined = log.events.len();
    }
    /// Everything that reaches `receiver` up to coordinate time `t` and has not been handed over yet, in order of
    /// arrival. Reading never goes backwards: a `t` earlier than the last read returns nothing new.
    pub fn read<T>(&mut self, log: &EventLog<T>, receiver: &Worldline<'_>, t: f64) -> &[Arrival] {
        self.sync(log, receiver);
        let start = self.delivered;
        let t = t.max(self.read_until);
        let count = self.pending[start..].partition_point(|a| a.t_arrival <= t);
        self.delivered = start + count;
        self.read_until = t;
        &self.pending[start..start + count]
    }
    /// The next arrival not yet delivered, without delivering it: when to wake a sleeping seat.
    pub fn peek<T>(&mut self, log: &EventLog<T>, receiver: &Worldline<'_>) -> Option<Arrival> {
        self.sync(log, receiver);
        self.pending.get(self.delivered).copied()
    }
    /// Everything delivered so far, in the order it was delivered.
    pub fn history(&self) -> &[Arrival] {
        &self.pending[..self.delivered]
    }
    /// The coordinate time this inbox has been read up to.
    pub fn read_until(&self) -> f64 {
        self.read_until
    }
    /// The receiver's future has changed: forget arrivals not yet delivered and work them out again from the log
    /// against the new worldline. What was delivered stays delivered.
    pub fn replan<T>(&mut self, log: &EventLog<T>, receiver: &Worldline<'_>) {
        let done: std::collections::BTreeSet<u64> = self.pending[..self.delivered]
            .iter()
            .map(|a| a.id.0)
            .collect();
        self.pending.truncate(self.delivered);
        for e in log.events.iter().filter(|e| !done.contains(&e.id.0)) {
            let t_arrival = receiver.recv(e.at);
            if t_arrival.is_finite() {
                self.pending.push(Arrival {
                    id: e.id,
                    t_emitted: e.at.t,
                    t_arrival,
                    tau_arrival: receiver.tau(t_arrival),
                });
            }
        }
        let from = self.delivered;
        self.pending[from..].sort_by(|a, b| a.key().partial_cmp(&b.key()).expect("finite keys"));
        self.examined = log.events.len();
    }
}
