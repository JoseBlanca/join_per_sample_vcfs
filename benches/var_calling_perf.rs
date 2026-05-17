//! Variant-calling pipeline end-to-end throughput.
//!
//! Four iterator stages, four bench groups:
//!
//! - `var_calling_merger/*` — `PerPositionMerger` (Stage 3) k-way merge
//!   over per-sample synthetic iterators. Drives the inner linear-scan
//!   over peeked heads.
//! - `var_calling_grouper/*` — `VariantGrouper` (Stage 4) overlap
//!   bundler on top of the merger. Drives the pure-REF drop loop and
//!   the transitive extension chain.
//! - `var_calling_per_group_merger/*` — `PerGroupMerger` (Stage 5)
//!   allele unification + likelihood reconstruction. Drives compound
//!   detection, scalar projection, the closed-form likelihood, and the
//!   rayon-parallel batch path.
//! - `var_calling_posterior_engine/*` — `PosteriorEngine` (Stage 6)
//!   EM-with-HWE prior, optional contamination mixture pre-pass.
//!   Inputs are pre-merged `MergedRecord` vectors built outside the
//!   timed region so the bench measures the engine in isolation.
//!
//! Workloads are synthetic in-memory `Vec`s so each `cargo bench` run
//! is reproducible. Inputs are built outside the timed regions.

// Opt-in mimalloc global allocator (cargo bench --features alloc-mimalloc).
#[cfg(feature = "alloc-mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::hint::black_box;
use std::sync::Arc;
use std::time::Duration;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use merge_per_sample_vcfs::per_sample_pileup::pileup::{
    AlleleObservation, AlleleSupportStats, ChainId, PileupRecord, RefSeqFetcher,
};
use merge_per_sample_vcfs::per_sample_pileup::psp::PspReadError;
use merge_per_sample_vcfs::var_calling::contamination_estimation::ContaminationEstimates;
use merge_per_sample_vcfs::var_calling::per_group_merger::{
    MergedRecord, PerGroupMerger, PerGroupMergerConfig, SharedRefFetcher,
};
use merge_per_sample_vcfs::var_calling::per_position_merger::{
    PerPositionMerger, PerPositionMergerError, PerPositionPileups,
};
use merge_per_sample_vcfs::var_calling::posterior_engine::backends::{
    ExactMath, InterpUnivariateMath, MathBackend,
};
use merge_per_sample_vcfs::var_calling::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig,
};
use merge_per_sample_vcfs::var_calling::variant_grouping::{
    GrouperConfig, OverlappingVariantGroup, VariantGrouper,
};

// ---------------------------------------------------------------------
// Common helpers
// ---------------------------------------------------------------------

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Mock fetcher backed by an in-memory reference buffer. Returns
/// uppercase ASCII per the `RefSeqFetcher` trait contract.
struct InMemRef {
    seq: Vec<u8>,
    /// 1-based position of `seq[0]`.
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

fn shared_fetcher(seq: Vec<u8>, base_offset: u32) -> SharedRefFetcher {
    Arc::new(InMemRef { seq, base_offset })
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

// ---------------------------------------------------------------------
// per_position_merger benches
// ---------------------------------------------------------------------

type MergerIter = std::vec::IntoIter<Result<PileupRecord, PspReadError>>;

/// Build N per-sample iterators each carrying records at every
/// position in `1..=n_positions`. The merger has to walk every head
/// at every output position — worst case for the linear-scan inner
/// loop.
fn build_dense_per_sample_streams(n_samples: usize, n_positions: u32) -> Vec<MergerIter> {
    (0..n_samples)
        .map(|s| {
            let records: Vec<Result<PileupRecord, PspReadError>> = (1..=n_positions)
                .map(|p| {
                    let ref_base = BASES[(p as usize + s) & 3];
                    Ok(PileupRecord::new(0, p, vec![ref_obs(ref_base, 30)]))
                })
                .collect();
            records.into_iter()
        })
        .collect()
}

/// Build N per-sample iterators where each sample only covers half
/// the positions (alternating). Exercises the merger's tied-readers
/// advance path.
fn build_sparse_per_sample_streams(n_samples: usize, n_positions: u32) -> Vec<MergerIter> {
    (0..n_samples)
        .map(|s| {
            let mut records: Vec<Result<PileupRecord, PspReadError>> =
                Vec::with_capacity((n_positions / 2) as usize);
            for p in 1..=n_positions {
                if ((p as usize) + s).is_multiple_of(2) {
                    let ref_base = BASES[(p as usize) & 3];
                    records.push(Ok(PileupRecord::new(0, p, vec![ref_obs(ref_base, 30)])));
                }
            }
            records.into_iter()
        })
        .collect()
}

fn names(n: usize) -> Vec<String> {
    (0..n).map(|i| format!("S{i}")).collect()
}

fn bench_per_position_merger(c: &mut Criterion) {
    let mut group = c.benchmark_group("var_calling_merger");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    const N_SAMPLES: usize = 64;
    const N_POSITIONS: u32 = 200_000;

    group.throughput(Throughput::Elements(N_POSITIONS as u64));

    group.bench_function(
        format!("dense_{N_SAMPLES}_samples_{N_POSITIONS}_positions"),
        |b| {
            b.iter_batched(
                || build_dense_per_sample_streams(N_SAMPLES, N_POSITIONS),
                |streams| {
                    let merger =
                        PerPositionMerger::new(streams, names(N_SAMPLES), Vec::new()).unwrap();
                    let mut count: u64 = 0;
                    for item in merger {
                        let pp = item.unwrap();
                        black_box(&pp);
                        count += 1;
                    }
                    black_box(count)
                },
                criterion::BatchSize::LargeInput,
            );
        },
    );

    group.bench_function(
        format!("sparse_{N_SAMPLES}_samples_{N_POSITIONS}_positions"),
        |b| {
            b.iter_batched(
                || build_sparse_per_sample_streams(N_SAMPLES, N_POSITIONS),
                |streams| {
                    let merger =
                        PerPositionMerger::new(streams, names(N_SAMPLES), Vec::new()).unwrap();
                    let mut count: u64 = 0;
                    for item in merger {
                        let pp = item.unwrap();
                        black_box(&pp);
                        count += 1;
                    }
                    black_box(count)
                },
                criterion::BatchSize::LargeInput,
            );
        },
    );

    group.finish();
}

// ---------------------------------------------------------------------
// variant_grouper benches
// ---------------------------------------------------------------------

type GrouperUpstream = std::vec::IntoIter<Result<PerPositionPileups, PerPositionMergerError>>;

/// SNP-dense per-position stream: every position carries a record;
/// every `snp_every_n`th position has a non-REF allele. Drives the
/// pure-REF drop path and the seed-extend path equally.
fn build_snp_dense_pp_stream(
    n_samples: usize,
    n_positions: u32,
    snp_every_n: u32,
) -> Vec<Result<PerPositionPileups, PerPositionMergerError>> {
    let mut out = Vec::with_capacity(n_positions as usize);
    for p in 1..=n_positions {
        let mut per_sample: Vec<Option<PileupRecord>> = (0..n_samples).map(|_| None).collect();
        // Two samples contribute records at each position: one REF-only,
        // and (every snp_every_n) one with an ALT.
        let ref_base = BASES[(p as usize) & 3];
        let alt_base = BASES[(p as usize + 1) & 3];
        per_sample[0] = Some(PileupRecord::new(0, p, vec![ref_obs(ref_base, 30)]));
        if p % snp_every_n == 0 {
            per_sample[1] = Some(PileupRecord::new(
                0,
                p,
                vec![
                    ref_obs(ref_base, 4),
                    alt_obs(alt_base, 26, -52.0, vec![p as u64]),
                ],
            ));
        } else {
            per_sample[1] = Some(PileupRecord::new(0, p, vec![ref_obs(ref_base, 30)]));
        }
        out.push(Ok(PerPositionPileups {
            chrom_id: 0,
            pos: p,
            per_sample,
        }));
    }
    out
}

/// Builds a stream where every 50 positions contains a deletion-shaped
/// record (ref_span = 5) followed by SNPs that fall inside its
/// footprint, transitively extending the group. Drives the extension
/// chain.
fn build_overlap_extension_pp_stream(
    n_samples: usize,
    n_groups: u32,
    spacing: u32,
) -> Vec<Result<PerPositionPileups, PerPositionMergerError>> {
    let mut out = Vec::with_capacity((n_groups * 6) as usize);
    for g in 0..n_groups {
        let base = g * spacing + 100;
        // Deletion span-5 at base (covers base..base+4).
        let mut per_sample: Vec<Option<PileupRecord>> = (0..n_samples).map(|_| None).collect();
        per_sample[0] = Some(PileupRecord::new(
            0,
            base,
            vec![
                AlleleObservation::new(b"AAAAA".to_vec(), support(2, 0.0), Vec::new()),
                AlleleObservation::new(b"A".to_vec(), support(8, -16.0), vec![base as u64]),
            ],
        ));
        out.push(Ok(PerPositionPileups {
            chrom_id: 0,
            pos: base,
            per_sample,
        }));
        // SNPs at base+1, base+2, base+3, base+4 in other samples.
        for offset in 1..=4_u32 {
            let mut per_sample: Vec<Option<PileupRecord>> = (0..n_samples).map(|_| None).collect();
            let pos = base + offset;
            per_sample[1 % n_samples] = Some(PileupRecord::new(
                0,
                pos,
                vec![ref_obs(b'A', 6), alt_obs(b'T', 24, -48.0, vec![pos as u64])],
            ));
            out.push(Ok(PerPositionPileups {
                chrom_id: 0,
                pos,
                per_sample,
            }));
        }
    }
    out
}

fn bench_variant_grouper(c: &mut Criterion) {
    let mut group = c.benchmark_group("var_calling_grouper");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    const N_SAMPLES: usize = 16;

    // SNP-dense: 200 k positions, ~1 in 200 carries a variant. Stresses
    // the pure-REF drop fast path.
    const N_POSITIONS_SNP: u32 = 200_000;
    group.throughput(Throughput::Elements(N_POSITIONS_SNP as u64));
    group.bench_function(
        format!("snp_dense_{N_SAMPLES}_samples_{N_POSITIONS_SNP}_pos"),
        |b| {
            b.iter_batched(
                || {
                    let pp_stream = build_snp_dense_pp_stream(N_SAMPLES, N_POSITIONS_SNP, 200);
                    pp_stream.into_iter()
                },
                |upstream: GrouperUpstream| {
                    let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
                    let mut count: u64 = 0;
                    for item in grouper {
                        let g = item.unwrap();
                        black_box(&g);
                        count += 1;
                    }
                    black_box(count)
                },
                criterion::BatchSize::LargeInput,
            );
        },
    );

    // Overlap-extension: 5 000 groups × 5 records each.
    const N_GROUPS_EXT: u32 = 5_000;
    group.throughput(Throughput::Elements(N_GROUPS_EXT as u64));
    group.bench_function(
        format!("overlap_extension_{N_SAMPLES}_samples_{N_GROUPS_EXT}_groups"),
        |b| {
            b.iter_batched(
                || {
                    let pp_stream = build_overlap_extension_pp_stream(N_SAMPLES, N_GROUPS_EXT, 100);
                    pp_stream.into_iter()
                },
                |upstream: GrouperUpstream| {
                    let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
                    let mut count: u64 = 0;
                    for item in grouper {
                        let g = item.unwrap();
                        black_box(&g);
                        count += 1;
                    }
                    black_box(count)
                },
                criterion::BatchSize::LargeInput,
            );
        },
    );

    group.finish();
}

// ---------------------------------------------------------------------
// per_group_merger (Stage 5) benches
// ---------------------------------------------------------------------

/// Build N biallelic SNP groups, one record per group, 64 samples.
/// Every sample carries an observation at the SNP locus: half are
/// REF-leaning, half ALT-leaning. No compounds.
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
            // Even samples ALT-leaning, odd samples REF-leaning.
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

/// Build N groups, each spanning two records (positions `pos` and
/// `pos + 2`). A subset of samples in each group carries chain ids
/// linking the two non-REF alleles → cross-record compound. Exercises
/// the compound-detection + scalar-projection + chain-broken-fallback
/// paths.
fn build_compound_groups(
    n_groups: u32,
    n_samples: usize,
    spacing: u32,
    chain_anchored_fraction: f32,
) -> (Vec<OverlappingVariantGroup>, Vec<u8>) {
    let mut groups = Vec::with_capacity(n_groups as usize);
    let last_pos = (n_groups - 1) * spacing + 100 + 2;
    let ref_len = (last_pos - 100 + 1) as usize;
    let ref_seq: Vec<u8> = (0..ref_len).map(|i| BASES[i & 3]).collect();

    for g in 0..n_groups {
        let p1 = g * spacing + 100;
        let p2 = p1 + 2;
        let ref_p1 = ref_seq[(p1 - 100) as usize];
        let ref_p2 = ref_seq[(p2 - 100) as usize];
        let alt_p1 = BASES[((p1 - 100) as usize + 1) & 3];
        let alt_p2 = BASES[((p2 - 100) as usize + 1) & 3];

        let mut per_sample_p1: Vec<Option<PileupRecord>> = (0..n_samples).map(|_| None).collect();
        let mut per_sample_p2: Vec<Option<PileupRecord>> = (0..n_samples).map(|_| None).collect();

        let anchor_cut = (n_samples as f32 * chain_anchored_fraction) as usize;
        for s in 0..n_samples {
            let chain_id_shared = (g as u64) * 1_000_000 + s as u64;
            let chain_id_p1_only = chain_id_shared + 7_000_000_000;
            let chain_id_p2_only = chain_id_shared + 9_000_000_000;
            if s < anchor_cut {
                // Chain-anchored: ALT at p1 and ALT at p2 share a chain id.
                per_sample_p1[s] = Some(PileupRecord::new(
                    0,
                    p1,
                    vec![
                        ref_obs(ref_p1, 4),
                        alt_obs(alt_p1, 8, -16.0, vec![chain_id_shared]),
                    ],
                ));
                per_sample_p2[s] = Some(PileupRecord::new(
                    0,
                    p2,
                    vec![
                        ref_obs(ref_p2, 4),
                        alt_obs(alt_p2, 8, -16.0, vec![chain_id_shared]),
                    ],
                ));
            } else {
                // Chain-broken: different chain ids on each ALT.
                per_sample_p1[s] = Some(PileupRecord::new(
                    0,
                    p1,
                    vec![
                        ref_obs(ref_p1, 4),
                        alt_obs(alt_p1, 8, -16.0, vec![chain_id_p1_only]),
                    ],
                ));
                per_sample_p2[s] = Some(PileupRecord::new(
                    0,
                    p2,
                    vec![
                        ref_obs(ref_p2, 4),
                        alt_obs(alt_p2, 8, -16.0, vec![chain_id_p2_only]),
                    ],
                ));
            }
        }

        groups.push(OverlappingVariantGroup {
            chrom_id: 0,
            start: p1,
            end: p2,
            records: vec![
                PerPositionPileups {
                    chrom_id: 0,
                    pos: p1,
                    per_sample: per_sample_p1,
                },
                PerPositionPileups {
                    chrom_id: 0,
                    pos: p2,
                    per_sample: per_sample_p2,
                },
            ],
        });
    }
    (groups, ref_seq)
}

fn bench_per_group_merger(c: &mut Criterion) {
    let mut group = c.benchmark_group("var_calling_per_group_merger");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(15));

    const N_SAMPLES: usize = 64;
    const N_GROUPS: u32 = 10_000;
    const SPACING: u32 = 10;

    // --- biallelic SNP, no compounds ---
    let (groups_snp, ref_seq_snp) = build_biallelic_snp_groups(N_GROUPS, N_SAMPLES, SPACING);
    let fetcher_snp = shared_fetcher(ref_seq_snp.clone(), 100);
    group.throughput(Throughput::Elements(N_GROUPS as u64));
    group.bench_function(
        format!("biallelic_snp_{N_SAMPLES}_samples_{N_GROUPS}_groups"),
        |b| {
            b.iter_batched(
                || groups_snp.clone(),
                |groups| {
                    let merger = PerGroupMerger::with_config(
                        groups.into_iter().map(Ok),
                        Arc::clone(&fetcher_snp),
                        PerGroupMergerConfig::default(),
                    );
                    let mut count: u64 = 0;
                    for item in merger {
                        let r = item.unwrap();
                        black_box(&r);
                        count += 1;
                    }
                    black_box(count)
                },
                criterion::BatchSize::LargeInput,
            );
        },
    );

    // --- compound-heavy: every sample chain-anchored ---
    let (groups_compound_full, ref_seq_compound) =
        build_compound_groups(N_GROUPS, N_SAMPLES, SPACING, 1.0);
    let fetcher_compound = shared_fetcher(ref_seq_compound.clone(), 100);
    group.bench_function(
        format!("compound_all_anchored_{N_SAMPLES}_samples_{N_GROUPS}_groups"),
        |b| {
            b.iter_batched(
                || groups_compound_full.clone(),
                |groups| {
                    let merger = PerGroupMerger::with_config(
                        groups.into_iter().map(Ok),
                        Arc::clone(&fetcher_compound),
                        PerGroupMergerConfig::default(),
                    );
                    let mut count: u64 = 0;
                    for item in merger {
                        let r = item.unwrap();
                        black_box(&r);
                        count += 1;
                    }
                    black_box(count)
                },
                criterion::BatchSize::LargeInput,
            );
        },
    );

    // --- compound mixed: half the samples chain-broken ⇒ fallback path ---
    let (groups_compound_mixed, _) = build_compound_groups(N_GROUPS, N_SAMPLES, SPACING, 0.5);
    group.bench_function(
        format!("compound_half_chain_broken_{N_SAMPLES}_samples_{N_GROUPS}_groups"),
        |b| {
            b.iter_batched(
                || groups_compound_mixed.clone(),
                |groups| {
                    let merger = PerGroupMerger::with_config(
                        groups.into_iter().map(Ok),
                        Arc::clone(&fetcher_compound),
                        PerGroupMergerConfig::default(),
                    );
                    let mut count: u64 = 0;
                    for item in merger {
                        let r = item.unwrap();
                        black_box(&r);
                        count += 1;
                    }
                    black_box(count)
                },
                criterion::BatchSize::LargeInput,
            );
        },
    );

    // --- single-sample biallelic SNP (no rayon contention) ---
    let (groups_solo, ref_seq_solo) = build_biallelic_snp_groups(N_GROUPS, 1, SPACING);
    let fetcher_solo = shared_fetcher(ref_seq_solo.clone(), 100);
    group.bench_function(format!("biallelic_snp_1_sample_{N_GROUPS}_groups"), |b| {
        b.iter_batched(
            || groups_solo.clone(),
            |groups| {
                let merger = PerGroupMerger::with_config(
                    groups.into_iter().map(Ok),
                    Arc::clone(&fetcher_solo),
                    PerGroupMergerConfig::default(),
                );
                let mut count: u64 = 0;
                for item in merger {
                    let r = item.unwrap();
                    black_box(&r);
                    count += 1;
                }
                black_box(count)
            },
            criterion::BatchSize::LargeInput,
        );
    });

    group.finish();
}

// ---------------------------------------------------------------------
// posterior_engine benches (Stage 6)
// ---------------------------------------------------------------------

/// Drive `PerGroupMerger` to completion and return the collected
/// `MergedRecord` stream. Done **once per fixture**, outside any timed
/// region — Stage 6's per-record cost is what we want to measure, not
/// the merger's.
fn pre_merge(
    groups: Vec<OverlappingVariantGroup>,
    fetcher: SharedRefFetcher,
    merger_config: PerGroupMergerConfig,
) -> Vec<MergedRecord> {
    PerGroupMerger::with_config(groups.into_iter().map(Ok), fetcher, merger_config)
        .collect::<Result<Vec<_>, _>>()
        .expect("merger fixture produced an error")
}

/// Representative single-batch contamination estimates: every sample
/// belongs to batch 0, `c_s = 3 %` for everyone, and the contaminant
/// source distribution `q_b` is `[0.6 REF, 0.3 SnpAlt, 0.1 IndelAlt]`.
/// Modest but visible — exercises the mixture pre-pass without
/// driving every cell to the `mix <= 0` defensive branch.
fn representative_contamination(n_samples: usize) -> ContaminationEstimates {
    let c_s_per_sample = vec![Some(0.03_f64); n_samples];
    let q_b_per_batch = vec![[0.6_f64, 0.3, 0.1]];
    let sample_to_batch = vec![0_usize; n_samples];
    ContaminationEstimates::from_user_supplied(c_s_per_sample, q_b_per_batch, sample_to_batch)
        .expect("representative_contamination inputs are valid")
}

/// Time `PosteriorEngine` on a pre-merged `Vec<MergedRecord>` under
/// an explicit math backend. Used twice per fixture (`exact`,
/// `interp_univariate`) so the bench summary lets us read off the
/// relative speedup side-by-side without flipping defaults in the
/// engine itself.
fn bench_posterior_drain<M: MathBackend + Copy + 'static>(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    name: &str,
    merged: &[MergedRecord],
    config: &PosteriorEngineConfig,
    math: M,
) {
    group.bench_function(name, |b| {
        b.iter_batched(
            || merged.to_vec(),
            |records| {
                let engine = PosteriorEngine::with_math_backend(
                    records.into_iter().map(Ok),
                    config.clone(),
                    math,
                );
                let mut count: u64 = 0;
                for item in engine {
                    let r = item.expect("posterior fixture produced an error");
                    black_box(&r);
                    count += 1;
                }
                black_box(count)
            },
            criterion::BatchSize::LargeInput,
        );
    });
}

fn bench_posterior_engine(c: &mut Criterion) {
    let mut group = c.benchmark_group("var_calling_posterior_engine");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(15));

    const N_SAMPLES: usize = 64;
    const N_GROUPS: u32 = 10_000;
    const SPACING: u32 = 10;

    group.throughput(Throughput::Elements(N_GROUPS as u64));

    let (groups_snp, ref_seq_snp) = build_biallelic_snp_groups(N_GROUPS, N_SAMPLES, SPACING);
    let fetcher_snp = shared_fetcher(ref_seq_snp.clone(), 100);
    let merged_diploid = pre_merge(
        groups_snp.clone(),
        Arc::clone(&fetcher_snp),
        PerGroupMergerConfig::default(),
    );

    let config_no_contam = PosteriorEngineConfig::with_project_defaults();
    let mut config_contam = PosteriorEngineConfig::with_project_defaults();
    config_contam.contamination = Some(representative_contamination(N_SAMPLES));

    // --- biallelic SNP, no contamination (dominant production case) ---
    bench_posterior_drain(
        &mut group,
        &format!("biallelic_snp_no_contam_{N_SAMPLES}_samples_{N_GROUPS}_groups/exact"),
        &merged_diploid,
        &config_no_contam,
        ExactMath,
    );
    bench_posterior_drain(
        &mut group,
        &format!(
            "biallelic_snp_no_contam_{N_SAMPLES}_samples_{N_GROUPS}_groups/interp_univariate"
        ),
        &merged_diploid,
        &config_no_contam,
        InterpUnivariateMath,
    );

    // --- biallelic SNP, contamination on (exercises mixture pre-pass) ---
    bench_posterior_drain(
        &mut group,
        &format!("biallelic_snp_contam_on_{N_SAMPLES}_samples_{N_GROUPS}_groups/exact"),
        &merged_diploid,
        &config_contam,
        ExactMath,
    );
    bench_posterior_drain(
        &mut group,
        &format!(
            "biallelic_snp_contam_on_{N_SAMPLES}_samples_{N_GROUPS}_groups/interp_univariate"
        ),
        &merged_diploid,
        &config_contam,
        InterpUnivariateMath,
    );

    // --- polyploid (ploidy = 3) biallelic SNP ---
    let mut triploid_merger_config = PerGroupMergerConfig::default();
    triploid_merger_config.ploidy = 3;
    let merged_triploid = pre_merge(groups_snp, Arc::clone(&fetcher_snp), triploid_merger_config);
    bench_posterior_drain(
        &mut group,
        &format!("biallelic_snp_ploidy_3_{N_SAMPLES}_samples_{N_GROUPS}_groups/exact"),
        &merged_triploid,
        &config_no_contam,
        ExactMath,
    );
    bench_posterior_drain(
        &mut group,
        &format!(
            "biallelic_snp_ploidy_3_{N_SAMPLES}_samples_{N_GROUPS}_groups/interp_univariate"
        ),
        &merged_triploid,
        &config_no_contam,
        InterpUnivariateMath,
    );

    group.finish();
}


criterion_group! {
    name = benches;
    config = Criterion::default();
    targets =
        bench_per_position_merger,
        bench_variant_grouper,
        bench_per_group_merger,
        bench_posterior_engine
}

criterion_main!(benches);
