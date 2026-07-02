//! Heap-profile the hidden-paralog per-locus scorer to quantify its per-locus
//! allocation load (the L3 finding in
//! `doc/devel/reports/reviews/perf_paralog-filter_2026-07-02.md`).
//!
//! Build and run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example dhat_paralog --features dhat-heap
//! ```
//!
//! Produces `dhat-heap.json` (open at
//! <https://nnethercote.github.io/dh_view/dh_view.html>). The precompute is
//! built **once** (as in a real calibrate / write pass), then `N_LOCI` loci are
//! scored in a loop, so every allocation the profile attributes to
//! `score_locus_for_paralogy` / `h1_log_likelihood` / `h2_log_likelihood` is a
//! genuine **per-locus** allocation — the thing L3 would hoist into reusable
//! scratch. The end-of-run summary prints total blocks so
//! `total_blocks / N_LOCI` reads directly as "allocations per scored locus".

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use pop_var_caller::paralog::{
    LocusObservations, ParalogModelParams, ParalogScorePrecompute, SampleObservation,
    score_locus_for_paralogy,
};

/// Samples per synthetic locus (the tomato2 usable cohort).
const N_SAMPLES: usize = 58;
/// Loci scored — enough that per-locus allocations dominate the one-time
/// fixture/precompute allocations in the profile.
const N_LOCI: usize = 50_000;
/// One-copy relative-depth SD (σ₀), matching the R1 tomato2 fit.
const SIGMA0: f64 = 0.28;

/// Build one synthetic locus with varied per-sample evidence (mirrors the
/// `paralog_scoring_perf` bench fixture): a mix of single-copy hets, skewed
/// paralog-like VAFs, and a few unusable (`None`) samples.
fn build_locus(n: usize) -> Vec<Option<SampleObservation>> {
    (0..n)
        .map(|i| {
            if i % 17 == 16 {
                return None;
            }
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

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let params = ParalogModelParams::default();
    let samples = build_locus(N_SAMPLES);
    let sds = vec![SIGMA0; N_SAMPLES];
    let inbreeding: Vec<f64> = samples
        .iter()
        .map(|o| o.map_or(0.0, |s| s.inbreeding_coefficient))
        .collect();

    // Built once, as in a real pass — its allocations are one-time, not per-locus.
    let precompute = ParalogScorePrecompute::new(&params, &inbreeding);

    let observations = LocusObservations { samples: &samples };
    let mut acc = 0.0f64;
    for _ in 0..N_LOCI {
        let score = score_locus_for_paralogy(&observations, &sds, &precompute);
        // Consume the result so the loop is not optimised away.
        acc += score.paralog_log_likelihood_ratio;
    }
    eprintln!(
        "scored {N_LOCI} loci (LR checksum {acc:.3}); total_blocks / {N_LOCI} = allocs/locus"
    );
}
