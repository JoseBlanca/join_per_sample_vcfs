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

use std::fs::{File, create_dir_all};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use tempfile::TempDir;

use pop_var_caller::per_sample_pileup::ref_fetcher::StreamingChromRefFetcher;
use pop_var_caller::pileup_record::{AlleleObservation, AlleleSupportStats, ChainId, PileupRecord};
use pop_var_caller::var_calling::contamination_estimation::ContaminationEstimates;
use pop_var_caller::var_calling::per_group_merger::{
    MergedRecord, PerGroupMerger, PerGroupMergerConfig, SharedRefFetcher,
};
use pop_var_caller::var_calling::per_position_merger::PerPositionPileups;
use pop_var_caller::var_calling::posterior_engine::backends::{
    ExactMath, InterpUnivariateMath, InterpUnivariateSimdMath, MathBackend,
};
use pop_var_caller::var_calling::posterior_engine::{PosteriorEngine, PosteriorEngineConfig};
use pop_var_caller::var_calling::variant_grouping::OverlappingVariantGroup;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

const N_SAMPLES: usize = 64;
const N_GROUPS: u32 = 10_000;
const SPACING: u32 = 10;
const N_RUNS: usize = 30;

/// Write a single-contig FASTA + `.fai` to `dir`. Contig layout is
/// `99 N` + `seq` so 1-based position 100 maps to `seq[0]`.
fn write_fasta(dir: &Path, contig_name: &str, seq: &[u8]) -> PathBuf {
    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let mut fa = File::create(&fasta_path).expect("create fasta");
    writeln!(fa, ">{contig_name}").expect("fa header");
    let header_len = (contig_name.len() + 2) as u64;
    let prefix = vec![b'N'; 99];
    fa.write_all(&prefix).expect("fa prefix");
    fa.write_all(seq).expect("fa seq");
    fa.write_all(b"\n").expect("fa nl");
    let total_len = prefix.len() + seq.len();
    let mut fai = File::create(&fai_path).expect("create fai");
    writeln!(
        fai,
        "{contig_name}\t{total_len}\t{header_len}\t{total_len}\t{}",
        total_len + 1
    )
    .expect("fai write");
    fasta_path
}

fn build_shared_fetcher(dir: &Path, seq: &[u8]) -> SharedRefFetcher {
    let fasta_path = write_fasta(dir, "chr0", seq);
    let streaming =
        StreamingChromRefFetcher::for_contig(&fasta_path, "chr0").expect("for_contig(chr0)");
    // See cohort_driver.rs:481 — SharedRefFetcher is Arc<dyn ChromRefFetcher + Send>
    // but the concrete StreamingChromRefFetcher is !Sync (RefCell interior mutability).
    #[allow(clippy::arc_with_non_send_sync)]
    let fetcher: SharedRefFetcher = Arc::new(streaming);
    fetcher
}

fn support(num_obs: u32, q_sum: f64) -> AlleleSupportStats {
    AlleleSupportStats::new(num_obs, q_sum, num_obs / 2, num_obs / 4, num_obs / 8, 0, 0)
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
    create_dir_all("tmp").expect("mkdir tmp");
    let dir = TempDir::new_in("tmp").expect("tempdir");
    let fetcher = build_shared_fetcher(dir.path(), &ref_seq);
    let merged = pre_merge(groups, fetcher);
    let config = PosteriorEngineConfig::with_project_defaults()
        .with_contamination(Some(representative_contamination(N_SAMPLES)))
        .expect("with_contamination is infallible today");

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
