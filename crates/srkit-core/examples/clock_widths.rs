//! U384 against U512 as the engine's stored clock (question R23.8): the cost of the operations the engine does on
//! every event. Run in release mode:
//!
//!     cargo run --release --example clock_widths
//!
//! The same workload runs in WebAssembly through `srkit_bench_clock` (see `tools/bench_clocks.js`).
use srkit_core::{
    causal::{self, GridVec},
    time::{Duration, Timestamp},
};
use std::{hint::black_box, time::Instant};

fn bench<const N: usize>(label: &str, n: u32) {
    let epoch = Timestamp::<N>::from_decimal(&format!("3{}", "0".repeat(100))).unwrap();
    let tick = Duration::<N>::from_unsigned_millis(1);
    let dx = GridVec::new(40_175_000_000_000_000, -3_000_000_000, 77);
    let time = |name: &str, mut f: Box<dyn FnMut() -> u64>| {
        let start = Instant::now();
        let mut sum = 0u64;
        for _ in 0..n {
            sum = sum.wrapping_add(f());
        }
        let ns = start.elapsed().as_nanos() as f64 / f64::from(n);
        println!("  {label} {name:<34} {ns:8.1} ns   (checksum {sum})");
    };
    let mut t = epoch;
    time(
        "add one tick",
        Box::new(move || {
            t = black_box(t).checked_add(tick).unwrap();
            t.modulo_millis(7).unwrap()
        }),
    );
    let later = epoch
        .checked_add(Duration::from_unsigned_millis(123_456_789))
        .unwrap();
    time(
        "difference of two times",
        Box::new(move || {
            black_box(later)
                .signed_duration_since(black_box(epoch))
                .magnitude_millis()
                .digits()[0]
        }),
    );
    time(
        "compare two times",
        Box::new(move || u64::from(black_box(later) > black_box(epoch))),
    );
    time(
        "phase in a period (modulo)",
        Box::new(move || black_box(later).modulo_millis(31_557_600_000).unwrap()),
    );
    time(
        "light-travel ticks (integer sqrt)",
        Box::new(move || {
            causal::light_ticks::<N>(black_box(dx))
                .unwrap()
                .magnitude_millis()
                .digits()[0]
        }),
    );
    let arrival = causal::first_arrival(epoch, dx).unwrap();
    time(
        "has the signal arrived? (exact)",
        Box::new(move || u64::from(causal::arrived(black_box(epoch), black_box(arrival), dx))),
    );
    time(
        "to decimal text and back",
        Box::new(move || {
            let s = black_box(later).to_string();
            Timestamp::<N>::from_decimal(&s)
                .unwrap()
                .modulo_millis(1000)
                .unwrap()
        }),
    );
}

fn main() {
    let n = 200_000;
    println!("stored clock, {n} iterations each (native):");
    bench::<6>("U384", n);
    bench::<8>("U512", n);
    println!(
        "sizes: U384 timestamp {} bytes, U512 timestamp {} bytes",
        size_of::<Timestamp<6>>(),
        size_of::<Timestamp<8>>()
    );
}
