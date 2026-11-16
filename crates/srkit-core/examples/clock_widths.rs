//! The engine's stored clock at three widths and three tick lengths (questions R23.8 and R71.4): the cost of
//! the operations the engine does on every event, and how far each combination reaches. Run in release mode:
//!
//!     cargo run --release --example clock_widths
//!
//! The same workload runs in WebAssembly through `srkit_bench_clock` (see `tools/bench_clocks.js`).
use srkit_core::{
    causal::{self, GridVec, Resolution},
    time::{Duration, Timestamp},
};
use std::{hint::black_box, time::Instant};

/// Nanoseconds for one call of `f`: the fastest of five timed batches.
fn time(mut f: impl FnMut() -> u64) -> f64 {
    let n = 100_000;
    let mut best = f64::INFINITY;
    for _ in 0..5 {
        let start = Instant::now();
        let mut sum = 0u64;
        for _ in 0..n {
            sum = sum.wrapping_add(f());
        }
        black_box(sum);
        best = best.min(start.elapsed().as_nanos() as f64 / f64::from(n));
    }
    best
}

fn bench<const N: usize>(label: &str, res: Resolution) {
    // an epoch of 1e100 years where the width reaches it (every case but 384 bits at a nanosecond), in ticks
    let years = if N == 6 && res.ticks_per_second() > 1_000_000 {
        98
    } else {
        100
    };
    let digits = years + (res.ticks_per_second() as f64).log10().round() as usize;
    let epoch = Timestamp::<N>::from_decimal(&format!("31557600{}", "0".repeat(digits))).unwrap();
    let tick = Duration::<N>::from_unsigned_ticks(1);
    // Proxima's distance, a little off the axes, in grid steps
    let g = res.grid_per_metre() as i128;
    let dx = GridVec::new(40_175_000_000_000_000 * g, -3_000_000_000 * g, 77 * g);
    let year = 31_557_600 * res.ticks_per_second();
    let later = epoch
        .checked_add(Duration::from_unsigned_ticks(123_456_789))
        .unwrap();
    let arrival = causal::first_arrival_at(res, epoch, dx).unwrap();
    let mut t = epoch;
    let row = [
        time(|| {
            t = black_box(t).checked_add(tick).unwrap();
            1
        }),
        time(|| {
            black_box(later)
                .signed_duration_since(black_box(epoch))
                .magnitude_millis()
                .digits()[0]
        }),
        time(|| u64::from(black_box(later) > black_box(epoch))),
        time(|| black_box(later).modulo_millis(year).unwrap()),
        time(|| {
            causal::light_ticks_at::<N>(res, black_box(dx))
                .unwrap()
                .magnitude_millis()
                .digits()[0]
        }),
        time(|| {
            u64::from(causal::arrived_at(
                res,
                black_box(epoch),
                black_box(arrival),
                dx,
            ))
        }),
        time(|| {
            let s = black_box(later).to_string();
            Timestamp::<N>::from_decimal(&s)
                .unwrap()
                .modulo_millis(1000)
                .unwrap()
        }),
    ];
    let cells: Vec<String> = row.iter().map(|ns| format!("{ns:7.0}")).collect();
    // reach: 2^(64 N) ticks, in years
    let reach =
        (64.0 * N as f64) * 2f64.log10() - (31_557_600.0 * res.ticks_per_second() as f64).log10();
    println!(
        "{label:<22}{}   {:>3} bytes   1e{:.0} yr",
        cells.join(""),
        size_of::<Timestamp<N>>(),
        reach.floor()
    );
}

fn main() {
    let (ms, us, ns) = (
        Resolution::MILLISECOND_METRE,
        Resolution::MICROSECOND_METRE,
        Resolution::NANOSECOND_MILLIMETRE,
    );
    println!(
        "stored clock, nanoseconds per operation (native, fastest of five batches of 100,000)"
    );
    println!(
        "{:<22}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}   {:>9}   reach",
        "width, tick, grid", "add", "diff", "cmp", "phase", "ticks", "arrived", "text", "size"
    );
    bench::<6>("U384  ms  metre", ms);
    bench::<8>("U512  ms  metre", ms);
    bench::<16>("U1024 ms  metre", ms);
    bench::<6>("U384  us  metre", us);
    bench::<8>("U512  us  metre", us);
    bench::<16>("U1024 us  metre", us);
    bench::<6>("U384  ns  millimetre", ns);
    bench::<8>("U512  ns  millimetre", ns);
    bench::<16>("U1024 ns  millimetre", ns);
    println!("add: one tick later. diff: difference of two times. cmp: which is later. phase: the time modulo a year.");
    println!("ticks: light-travel time over a displacement (integer square root). arrived: the exact causal check.");
    println!("text: to decimal text and back. U384 at a nanosecond reaches 1.2e99 years, so its epoch here is 1e98.");
}
