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

/// A biallelic-SNP spill record for `n_samples`, each sample given a window mean
/// depth and an `(ref_obs, alt_obs)` AD by the `per_sample` closure.
pub(crate) fn snp_spill(
    pos: u32,
    n_samples: usize,
    gc: f32,
    per_sample: &dyn Fn(usize) -> (f32, u32, u32),
) -> ParalogSpillRecord {
    let mut scalars = Vec::with_capacity(n_samples * 2);
    let mut window_mean_depth = Vec::with_capacity(n_samples);
    for s in 0..n_samples {
        let (depth, ref_obs, alt_obs) = per_sample(s);
        scalars.push(support(ref_obs));
        scalars.push(support(alt_obs));
        window_mean_depth.push(Some(depth));
    }
    ParalogSpillRecord {
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
        window_gc: gc,
        window_mean_depth,
    }
}

/// A "normal" locus: every sample at ~1× coverage with a clean het — a real
/// single-copy variant, LR should be `< 0`.
pub(crate) fn normal_locus(pos: u32, n: usize) -> ParalogSpillRecord {
    snp_spill(pos, n, 0.5, &|_| (20.0, 10, 10))
}

/// A "paralog-like" locus: every sample at ~2× coverage with a skewed VAF (~1/3,
/// a collapsed-paralog PSV), LR should be `> 0`.
pub(crate) fn paralog_locus(pos: u32, n: usize) -> ParalogSpillRecord {
    snp_spill(pos, n, 0.5, &|_| (40.0, 20, 10))
}

/// Serialise records to spill bytes.
pub(crate) fn write_spill(records: &[ParalogSpillRecord]) -> Vec<u8> {
    let mut w = ParalogSpillWriter::new(std::io::Cursor::new(Vec::new()));
    for r in records {
        w.append(r).unwrap();
    }
    w.finish().unwrap().into_inner()
}
