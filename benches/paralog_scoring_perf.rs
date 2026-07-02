//! Hidden-paralog per-locus scorer throughput.
//!
//! The committed perf surface for the paralog filter's hot kernel,
//! [`score_locus_for_paralogy`](pop_var_caller::paralog::score_locus_for_paralogy).
//! The on-by-default filter's cost is this pure function run once per locus per
//! pass, over 376 k loci × 58 samples, **twice** (calibrate + write) — the
//! +11.6× wall regression measured in
//! `doc/devel/reports/reviews/perf_paralog-filter_2026-07-02.md`. This bench
//! isolates the per-locus cost so the non-parallel levers (memoising the
//! locus-invariant Wright / carrier log-prior tables, hoisting the carrier
//! configs, reusing the per-locus scratch `Vec`s) can be measured before the
//! larger parallelisation change.
//!
//! ## What it measures
//!
//! `paralog_score_locus/N{n}` — one call to `score_locus_for_paralogy` over an
//! `n`-sample synthetic locus with the production
//! [`ParalogModelParams::default`] grids (H1: 200 SFS points; H2: 40
//! carrier-freq points × 7 configs). `n ∈ {20, 58, 200}` brackets the tomato2
//! cohort (58 usable). The per-locus cost is ~linear in `n` and in the grid
//! sizes, so this is the quantity every kernel-level lever moves.
//!
//! The cohort is built **once outside the timed region** and borrowed; only the
//! `score_locus_for_paralogy` call is timed (a fixture rebuild inside `b.iter`
//! would dominate the sample). Inputs and the returned `ParalogScore` are
//! `black_box`'d so the optimiser cannot hoist the call out or constant-fold the
//! marginals.
//!
//! A separate streaming bench over `calibrate` (decode → join → score → fold)
//! is deferred to the parallelisation phase: it needs `pub(crate)` spill
//! internals, and the decode/join-vs-score split it would expose matters for
//! that change, not for the kernel levers this bench gates.

use std::hint::black_box;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

use pop_var_caller::paralog::{
    LocusObservations, ParalogModelParams, ParalogScorePrecompute, SampleObservation,
    score_locus_for_paralogy,
};

/// The one-copy relative-depth SD (σ₀) every synthetic sample carries. Matches
/// the R1 tomato2 fit (σ₀ ≈ 0.28), so the Normal coverage terms sit in a
/// realistic range rather than a degenerate one.
const SIGMA0: f64 = 0.28;

/// Build an `n`-sample synthetic locus with varied per-sample evidence, so the
/// scorer does representative work and the optimiser cannot constant-fold across
/// samples. The pattern deterministically spreads relative copy number, allele
/// depth, VAF, and inbreeding coefficient across the cohort — a mix of
/// single-copy hets, skewed (paralog-like) VAFs, and a few unusable (`None`)
/// samples, matching what a real variant locus feeds the scorer.
fn build_locus(n: usize) -> Vec<Option<SampleObservation>> {
    (0..n)
        .map(|i| {
            // Every 17th sample is unusable (no reads / dropped) — exercises the
            // `usable < n` compaction path without dominating the cohort.
            if i % 17 == 16 {
                return None;
            }
            // Spread copy number 1.0..~2.0 (single-copy through collapsed
            // paralog), total depth 12..40, and VAF between clean-het and
            // skewed. F cycles selfer↔outbred.
            let rel = 1.0 + (i % 5) as f64 * 0.25;
            let total = 12 + (i % 8) as u32 * 4;
            let alt = ((total as f64) * (0.30 + 0.05 * (i % 4) as f64)).round() as u32;
            let f = match i % 3 {
                0 => 0.0,
                1 => 0.3,
                _ => 0.9,
            };
            Some(SampleObservation {
                relative_copy_number: rel,
                alt_reads: alt.min(total),
                total_reads: total,
                inbreeding_coefficient: f,
            })
        })
        .collect()
}

fn bench_score_locus(c: &mut Criterion) {
    let mut group = c.benchmark_group("paralog_score_locus");
    let params = ParalogModelParams::default();

    for &n in &[20usize, 58, 200] {
        // Built once, outside the timed region, then borrowed. The precompute is
        // likewise built once per pass in production, so building it outside
        // `b.iter` measures the true per-locus cost (not the one-time table build).
        let samples = build_locus(n);
        let sds = vec![SIGMA0; n];
        let inbreeding: Vec<f64> = samples
            .iter()
            .map(|o| o.map_or(0.0, |s| s.inbreeding_coefficient))
            .collect();
        let precompute = ParalogScorePrecompute::new(&params, &inbreeding);

        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| {
                let score = score_locus_for_paralogy(
                    black_box(&LocusObservations { samples: &samples }),
                    black_box(&sds),
                    black_box(&precompute),
                );
                black_box(score)
            });
        });
    }

    group.finish();
}

fn config() -> Criterion {
    Criterion::default()
        .sample_size(50)
        .measurement_time(Duration::from_secs(5))
}

criterion_group! {
    name = benches;
    config = config();
    targets = bench_score_locus
}

criterion_main!(benches);
