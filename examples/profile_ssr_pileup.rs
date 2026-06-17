//! CPU-profile the `ssr-pileup` realignment hot path.
//!
//! Runs the synthetic per-locus realignment workload (`bench_harness`) in a
//! tight loop for a fixed wall-clock budget so a sampling profiler can attach to
//! a long-lived, single-threaded process.
//!
//! Host (macOS) sampling profile:
//! ```text
//! cargo build --release --example profile_ssr_pileup
//! ./target/release/examples/profile_ssr_pileup 30 &
//! sample $! 25 -file tmp/ssr_pileup_sample.txt
//! ```
//!
//! Args: `[seconds] [period] [units] [depth] [window] [clip_units]`
//! (defaults: 30s, dinucleotide, units=15, depth=30, window=6, clip=0).

use std::hint::black_box;
use std::time::{Duration, Instant};

use pop_var_caller::ssr_mark1::pileup::bench_harness::{
    analyze_workload, build_synthetic_workload,
};

fn arg<T: std::str::FromStr>(args: &[String], i: usize, default: T) -> T {
    args.get(i).and_then(|s| s.parse().ok()).unwrap_or(default)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let secs: u64 = arg(&args, 1, 30);
    let period: usize = arg(&args, 2, 2);
    let units: u16 = arg(&args, 3, 15);
    let depth: usize = arg(&args, 4, 30);
    let window: u16 = arg(&args, 5, 6);
    let clip: u16 = arg(&args, 6, 0);

    let workload = build_synthetic_workload(period, units, 50, depth, window, clip);
    eprintln!(
        "profiling {secs}s: period={period} units={units} depth={depth} window={window} clip={clip}"
    );

    let deadline = Instant::now() + Duration::from_secs(secs);
    let mut sink = 0.0f64;
    let mut iters: u64 = 0;
    while Instant::now() < deadline {
        // Batch between clock reads so the timer check is negligible.
        for _ in 0..200 {
            sink += analyze_workload(black_box(&workload));
            iters += 1;
        }
    }
    black_box(sink);
    eprintln!(
        "did {iters} workload iterations ({} read-analyses)",
        iters * depth as u64
    );
}
