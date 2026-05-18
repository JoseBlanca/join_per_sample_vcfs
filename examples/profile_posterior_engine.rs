//! Standalone driver for profiling the posterior engine (Stage 6).
//!
//! Builds the biallelic-contam-on fixture from the var_calling bench
//! once, then drains the engine repeatedly so a sampling profiler
//! (`samply` / `perf`) captures self-time inside `run_em_for_record`
//! and friends — not the merger's setup work that contaminates the
//! Criterion-based bench window when profiled directly.
//!
//! Build and run on the host (rootless podman blocks `perf_event_open`,
//! so samply / perf must be invoked outside the container):
//!
//! ```text
//! ./scripts/dev.sh cargo build --release --example profile_posterior_engine
//! samply record -- ./target-container/release/examples/profile_posterior_engine
//! perf record -F 997 -g --call-graph=dwarf -- ./target-container/release/examples/profile_posterior_engine
//! ```
//!
//! Wall-time of each drain is printed to stderr; aggregate totals at
//! the end give an off-Criterion baseline.

use std::sync::Arc;
use std::time::Instant;

use merge_per_sample_vcfs::per_sample_pileup::pileup::{
    AlleleObservation, AlleleSupportStats, ChainId, PileupRecord, RefSeqFetcher,
};
use merge_per_sample_vcfs::var_calling::contamination_estimation::ContaminationEstimates;
use merge_per_sample_vcfs::var_calling::per_group_merger::{
    MergedRecord, PerGroupMerger, PerGroupMergerConfig, SharedRefFetcher,
};
use merge_per_sample_vcfs::var_calling::per_position_merger::PerPositionPileups;
use merge_per_sample_vcfs::var_calling::posterior_engine::backends::{
    ExactMath, InterpUnivariateMath, InterpUnivariateSimdMath, MathBackend,
};
use merge_per_sample_vcfs::var_calling::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig,
};
use merge_per_sample_vcfs::var_calling::variant_grouping::OverlappingVariantGroup;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

const N_SAMPLES: usize = 64;
const N_GROUPS: u32 = 10_000;
const SPACING: u32 = 10;
const N_RUNS: usize = 30;

struct InMemRef {
    seq: Vec<u8>,
    base_offset: u32,
}

impl RefSeqFetcher for InMemRef {
    fn fetch(
        &self,
        _chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, std::io::Error> {
        let start_idx = (start_1based - self.base_offset) as usize;
        let end_idx = start_idx + length as usize;
        if end_idx > self.seq.len() {
            return Err(std::io::Error::other("out of range"));
        }
        Ok(self.seq[start_idx..end_idx].to_vec())
    }
}

fn support(num_obs: u32, q_sum: f64) -> AlleleSupportStats {
    AlleleSupportStats::new(num_obs, q_sum, num_obs / 2, num_obs / 4, num_obs / 8)
}

fn ref_obs(ref_base: u8, num_obs: u32) -> AlleleObservation {
    AlleleObservation::new(vec![ref_base], support(num_obs, 0.0), Vec::new())
}

fn alt_obs(alt_base: u8, num_obs: u32, q_sum: f64, chain_ids: Vec<ChainId>) -> AlleleObservation {
    AlleleObservation::new(vec![alt_base], support(num_obs, q_sum), chain_ids)
}

fn build_biallelic_snp_groups(
    n_groups: u32,
    n_samples: usize,
    spacing: u32,
) -> (Vec<OverlappingVariantGroup>, Vec<u8>) {
    let mut groups = Vec::with_capacity(n_groups as usize);
    let last_pos = (n_groups - 1) * spacing + 100;
    let ref_len = (last_pos - 100 + 1) as usize;
    let ref_seq: Vec<u8> = (0..ref_len).map(|i| BASES[i & 3]).collect();
    for g in 0..n_groups {
        let pos = g * spacing + 100;
        let ref_base = ref_seq[(pos - 100) as usize];
        let alt_base = BASES[((pos - 100) as usize + 1) & 3];
        let mut per_sample: Vec<Option<PileupRecord>> = (0..n_samples).map(|_| None).collect();
        for (s, slot) in per_sample.iter_mut().enumerate().take(n_samples) {
            let (n_ref, n_alt) = if s.is_multiple_of(2) {
                (4, 26)
            } else {
                (24, 6)
            };
            *slot = Some(PileupRecord::new(
                0,
                pos,
                vec![
                    ref_obs(ref_base, n_ref),
                    alt_obs(
                        alt_base,
                        n_alt,
                        -2.0 * n_alt as f64,
                        vec![pos as u64 + s as u64],
                    ),
                ],
            ));
        }
        groups.push(OverlappingVariantGroup {
            chrom_id: 0,
            start: pos,
            end: pos,
            records: vec![PerPositionPileups {
                chrom_id: 0,
                pos,
                per_sample,
            }],
        });
    }
    (groups, ref_seq)
}

fn pre_merge(groups: Vec<OverlappingVariantGroup>, fetcher: SharedRefFetcher) -> Vec<MergedRecord> {
    PerGroupMerger::with_config(
        groups.into_iter().map(Ok),
        fetcher,
        PerGroupMergerConfig::default(),
    )
    .collect::<Result<Vec<_>, _>>()
    .expect("merger fixture produced an error")
}

fn representative_contamination(n_samples: usize) -> ContaminationEstimates {
    ContaminationEstimates::from_user_supplied(
        vec![Some(0.03_f64); n_samples],
        vec![[0.6_f64, 0.3, 0.1]],
        vec![0_usize; n_samples],
    )
    .expect("representative_contamination inputs are valid")
}

fn drain<M: MathBackend + Copy>(
    merged: &[MergedRecord],
    config: &PosteriorEngineConfig,
    math: M,
) -> std::time::Duration {
    let mut total_records: u64 = 0;
    let runs_start = Instant::now();
    for run in 0..N_RUNS {
        let run_start = Instant::now();
        let records = merged.to_vec();
        let engine =
            PosteriorEngine::with_math_backend(records.into_iter().map(Ok), config.clone(), math);
        let mut n = 0_u64;
        for item in engine {
            let r = item.expect("engine produced an error");
            std::hint::black_box(&r);
            n += 1;
        }
        total_records += n;
        eprintln!(
            "run {:>2}/{N_RUNS}: {} records in {:?}",
            run + 1,
            n,
            run_start.elapsed()
        );
    }
    let elapsed = runs_start.elapsed();
    eprintln!(
        "total: {total_records} records across {N_RUNS} runs in {elapsed:?} ({:.2} us/record)",
        elapsed.as_secs_f64() * 1e6 / total_records as f64
    );
    elapsed
}

fn main() {
    eprintln!("building fixture ({N_GROUPS} groups, {N_SAMPLES} samples)…");
    let setup_start = Instant::now();

    let (groups, ref_seq) = build_biallelic_snp_groups(N_GROUPS, N_SAMPLES, SPACING);
    let fetcher: SharedRefFetcher = Arc::new(InMemRef {
        seq: ref_seq,
        base_offset: 100,
    });
    let merged = pre_merge(groups, fetcher);
    let mut config = PosteriorEngineConfig::with_project_defaults();
    config.contamination = Some(representative_contamination(N_SAMPLES));

    eprintln!(
        "fixture built in {:?} ({} merged records)",
        setup_start.elapsed(),
        merged.len()
    );

    // `POSTERIOR_BACKEND` selects the math backend for this run.
    // Defaults to ExactMath so existing scripts keep their meaning.
    let backend = std::env::var("POSTERIOR_BACKEND").unwrap_or_else(|_| "exact".to_string());
    match backend.as_str() {
        "exact" => {
            eprintln!("backend: ExactMath; draining engine {N_RUNS} times…");
            drain(&merged, &config, ExactMath);
        }
        "interp" => {
            eprintln!("backend: InterpUnivariateMath; draining engine {N_RUNS} times…");
            drain(&merged, &config, InterpUnivariateMath);
        }
        "simd" => {
            eprintln!("backend: InterpUnivariateSimdMath; draining engine {N_RUNS} times…");
            drain(&merged, &config, InterpUnivariateSimdMath);
        }
        other => {
            eprintln!("unknown POSTERIOR_BACKEND={other:?} (use 'exact' / 'interp' / 'simd')");
            std::process::exit(2);
        }
    }
}
