//! Shared test fixtures for the paralog-filter wiring (calibrate + write pass).
//!
//! A small synthetic cohort whose samples all share one single-copy coverage
//! model, plus spill-record builders for "normal" (single-copy) and
//! "paralog-like" (2× coverage, skewed VAF) biallelic SNP loci — enough for the
//! calibrate and write-pass tests to exercise the full score → flag → drop path
//! without real data.

#![cfg(test)]

use crate::paralog::CoverageFitConfig;
use crate::pileup_record::AlleleSupportStats;
use crate::sample_summary::{
    CoverageByGcHistogram, HetCounts, SAMPLE_SUMMARY_VERSION, SampleSummary,
};
use crate::var_calling::per_group_merger::MergedAllele;
use crate::var_calling::posterior_engine::{EmDiagnostics, PosteriorRecord, RecordLocus};
use crate::var_calling::types::LocusWindowCoverage;

use super::prepass::ParalogPrePass;
use super::spill::{ParalogSpillRecord, ParalogSpillWriter};

/// A single-copy coverage histogram (scale ≈ 20) every synthetic sample shares,
/// so `relative_copy_number(gc, depth) ≈ depth / 20`.
pub(crate) fn single_copy_summary() -> SampleSummary {
    let gc_bins = 4u32;
    let depth_bins = 200u32;
    let width = 0.5f64;
    let row_stride = depth_bins as usize + 1;
    let mut counts = vec![0u32; gc_bins as usize * row_stride];
    let mut n_tiles = 0u64;
    for gc_bin in 0..gc_bins as usize {
        for (depth, count) in [(19.0, 100u32), (20.0, 300), (21.0, 100)] {
            let db = ((depth / width) as usize).min(depth_bins as usize);
            counts[gc_bin * row_stride + db] += count;
            n_tiles += u64::from(count);
        }
    }
    SampleSummary {
        version: SAMPLE_SUMMARY_VERSION,
        coverage_by_gc: CoverageByGcHistogram {
            window_bp: 500,
            gc_bins,
            depth_bin_width: width,
            depth_bins,
            n_tiles,
            n_skipped_tiles: 0,
            callable_positions: 1_000_000,
            counts,
        },
        heterozygosity: HetCounts {
            n_het_sites: 5000,
            n_hom_alt_sites: 0,
            n_ambiguous_sites: 0,
            n_variant_sites: 5000,
            min_depth: 4,
            error_rate: 0.02,
            lr_margin: std::f64::consts::LN_10,
        },
    }
}

/// A pre-pass over `n_samples` identical single-copy samples.
pub(crate) fn prepass(n_samples: usize) -> ParalogPrePass {
    let summaries: Vec<Option<SampleSummary>> = (0..n_samples)
        .map(|_| Some(single_copy_summary()))
        .collect();
    ParalogPrePass::fit(&summaries, &CoverageFitConfig::default())
}

fn support(num_obs: u32) -> AlleleSupportStats {
    AlleleSupportStats::new(num_obs, 0.0, num_obs / 2, 0, 0, 0, 0)
}

/// A plain (non-compound) merged allele from its sequence bytes.
pub(crate) fn allele(seq: &[u8]) -> MergedAllele {
    MergedAllele {
        seq: seq.to_vec(),
        is_compound: false,
        constituents: Vec::new(),
    }
}

/// A per-sample window-coverage fixture: constant GC `gc` for every sample, and
/// each sample's mean coverage from `depth` (`None` → `NaN`, i.e. the sample has
/// no covered record at the locus and the score skips it).
pub(crate) fn window_coverage(
    n_samples: usize,
    gc: f32,
    depth: &dyn Fn(usize) -> Option<f32>,
) -> LocusWindowCoverage {
    LocusWindowCoverage {
        gc: vec![gc; n_samples],
        coverage: (0..n_samples).map(|s| depth(s).unwrap_or(f32::NAN)).collect(),
    }
}

/// A biallelic-SNP **record-spill** record for `n_samples`, each sample given an
/// `(ref_obs, alt_obs)` AD by `per_sample`, carrying its per-sample centred-
/// window coverage `window` (as the caller now gathers from the psp columns).
pub(crate) fn snp_record(
    pos: u32,
    n_samples: usize,
    per_sample: &dyn Fn(usize) -> (u32, u32),
    window: LocusWindowCoverage,
) -> ParalogSpillRecord {
    let mut scalars = Vec::with_capacity(n_samples * 2);
    for s in 0..n_samples {
        let (ref_obs, alt_obs) = per_sample(s);
        scalars.push(support(ref_obs));
        scalars.push(support(alt_obs));
    }
    ParalogSpillRecord {
        paralog_lr: f64::NAN,
        window_coverage: window,
        record: PosteriorRecord {
            locus: RecordLocus {
                chrom_id: 0,
                start: pos,
                end: pos,
            },
            alleles: vec![allele(b"A"), allele(b"T")],
            ploidy: 2,
            n_samples,
            n_genotypes: 3,
            allele_frequencies: vec![0.5, 0.5],
            compound_frequencies: vec![None, None],
            posteriors: vec![0.0; n_samples * 3],
            best_genotype: vec![1; n_samples],
            gq_phred: vec![50.0; n_samples],
            qual_phred: 100.0,
            scalars,
            other_scalars: Vec::new(),
            chain_anchor_flags: vec![false; n_samples * 2],
            diagnostics: EmDiagnostics {
                iterations: 1,
                final_max_delta_p: 1e-6,
                converged: true,
            },
        },
    }
}

/// A "normal" locus record: every sample a clean het (AD 10/10) — a real
/// single-copy variant — carrying ~1× coverage (depth 20 vs the single-copy
/// scale ≈ 20), GC 0.5. Its LR is `< 0`.
pub(crate) fn normal_locus(pos: u32, n: usize) -> ParalogSpillRecord {
    snp_record(pos, n, &|_| (10, 10), window_coverage(n, 0.5, &|_| Some(20.0)))
}

/// A "paralog-like" locus record: every sample a skewed VAF (~1/3, a
/// collapsed-paralog PSV; AD 20/10) carrying ~2× coverage (depth 40), GC 0.5.
/// Its LR is `> 0`.
pub(crate) fn paralog_locus(pos: u32, n: usize) -> ParalogSpillRecord {
    snp_record(pos, n, &|_| (20, 10), window_coverage(n, 0.5, &|_| Some(40.0)))
}

/// Serialise record-spill records to bytes.
pub(crate) fn write_spill(records: &[ParalogSpillRecord]) -> Vec<u8> {
    let mut w = ParalogSpillWriter::new(std::io::Cursor::new(Vec::new()));
    for r in records {
        w.append(r).unwrap();
    }
    w.finish().unwrap().into_inner()
}
