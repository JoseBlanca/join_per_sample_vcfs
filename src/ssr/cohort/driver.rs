//! The `ssr-call` driver — open the cohort, run the merge, write the output (arch
//! `doc/devel/architecture/ssr_call_reading.md` §5).
//!
//! **Phase 3: single-threaded.** The producer→queue→worker-pool→writer topology and
//! the genotyping EM are not built here. The worker would exist to overlap *expensive
//! EM work* with decode; with no EM yet (the genotyping doc owns it) and decode
//! parallelism being the separate Phase-5 prefetch pool, a thread pool around a
//! formatting stub would be unverifiable complexity. So the driver streams the merged
//! `CohortLocus`es straight to a **catalog-ordered TSV dump** of the reading layer — a
//! placeholder until the EM + VCF land, and a useful way to inspect a cohort's merged
//! evidence. `--threads` / `--queue-depth` are accepted but reserved.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Seek, Write};
use std::path::PathBuf;

use crate::ssr::cohort::merge::{CohortMerger, SsrMergeError};
use crate::ssr::cohort::param_estimation::{G0FitCfg, G0PseudocountDecay, ParamSet, PerBaseError};
use crate::ssr::cohort::prepass::EstimatedParams;
use crate::ssr::cohort::sample_groups::GroupedParams;
use crate::ssr::cohort::types::CohortLocus;

/// Inputs for an `ssr-call` run.
pub(crate) struct SsrCallConfig {
    /// The shared `.ssr.catalog`.
    pub(crate) catalog: PathBuf,
    /// The per-sample `.ssr.psp` evidence files.
    pub(crate) psp_files: Vec<PathBuf>,
    /// Output path (currently the TSV dump; the VCF lands with the EM).
    pub(crate) output: PathBuf,
    /// Reserved — the EM worker pool is not built yet (Phase 3 is single-threaded).
    pub(crate) threads: usize,
    /// Reserved — the bounded producer→worker queue is not built yet.
    pub(crate) queue_depth: usize,
}

/// Errors from an `ssr-call` run.
#[derive(Debug, thiserror::Error)]
pub(crate) enum SsrCallError {
    /// Opening / merging the cohort failed.
    #[error(transparent)]
    Merge(#[from] SsrMergeError),
    /// Writing the output failed.
    #[error("writing ssr-call output")]
    Write(#[from] std::io::Error),
    /// The pre-pass produced no confident genotype for one or more cohort samples, so
    /// their chemistry cannot be frozen (decision E, arch `ssr_call_driver.md` §3/§7).
    /// A sample that never resolves anywhere is a degenerate input (a blank/near-empty
    /// sample, a mis-supplied file, or a catalog/chemistry mismatch) — fail loud rather
    /// than silently assign it cohort-default chemistry. The driver maps the indices to
    /// sample names for the user-facing message.
    #[error(
        "{} cohort sample(s) produced no confident genotype in the pre-pass \
         (sample indices {samples:?}); cannot estimate their chemistry — check the inputs",
        samples.len()
    )]
    UnresolvedSamples { samples: Vec<u32> },
}

/// Assemble the frozen [`ParamSet`] the genotyping EM consumes from the pre-pass
/// outputs (arch `ssr_call_driver.md` §3). Pure, so it is unit-testable away from any
/// I/O.
///
/// **Decision E:** every cohort sample `0..n_samples` must appear in
/// `grouped.group_of_sample` (i.e. it produced ≥1 confident genotype in the pre-pass).
/// A missing sample is an [`SsrCallError::UnresolvedSamples`] hard error naming the
/// offending indices — never a silent cohort-default group.
///
/// **`G₀` coverage (review Mi1):** the fitted `est.g0_by_period` carries only periods
/// with variable loci. We backfill `g0_cfg.fallback_p` for every other period the
/// pre-pass characterized (those in `est.shape_by_period`), so `g0_cfg.fallback_p` is
/// the single fallback source of truth for every realistically-genotyped period rather
/// than leaving the EM to apply its own independent default. A period the pre-pass
/// never observed at all (no chemistry) is out of scope and falls to the EM default.
pub(crate) fn build_param_set(
    est: &EstimatedParams,
    grouped: &GroupedParams,
    n_samples: usize,
    g0_cfg: &G0FitCfg,
) -> Result<ParamSet, SsrCallError> {
    // Decision E: a dense group-of-sample over the whole cohort, or a hard error.
    let mut group_of_sample = Vec::with_capacity(n_samples);
    let mut unresolved = Vec::new();
    for sample in 0..n_samples as u32 {
        match grouped.group_of_sample.get(&sample) {
            Some(&group) => group_of_sample.push(group),
            None => unresolved.push(sample),
        }
    }
    if !unresolved.is_empty() {
        return Err(SsrCallError::UnresolvedSamples {
            samples: unresolved,
        });
    }

    // G₀: the fitted per-period decays, with every other characterized period backfilled
    // from the configured fallback (single source of truth — review Mi1).
    let mut pseudocount_decay_per_loci_group: HashMap<u8, G0PseudocountDecay> =
        est.g0_by_period.clone();
    for &period in est.shape_by_period.keys() {
        pseudocount_decay_per_loci_group
            .entry(period)
            .or_insert(G0PseudocountDecay {
                p: g0_cfg.fallback_p,
            });
    }

    Ok(ParamSet {
        error_per_sample_group: grouped
            .eps_per_group
            .iter()
            .map(|&e| PerBaseError(e))
            .collect(),
        stutter_shape_parent: est.shape_by_period.clone(),
        stutter_shape_by_cell: grouped.shape_by_group_period.clone(),
        level_seed: grouped.level_per_group.clone(),
        pseudocount_decay_per_loci_group,
        group_of_sample,
        // `F` is estimated and frozen by the burn-in (arch §4), not here.
        f0_seed: 0.0,
    })
}

/// Open the cohort, merge, and write the catalog-ordered dump to `config.output`.
pub(crate) fn run(config: &SsrCallConfig) -> Result<(), SsrCallError> {
    let merger = CohortMerger::open(&config.catalog, &config.psp_files)?;
    let chrom_names = merger.chrom_names();
    let mut out = BufWriter::new(File::create(&config.output)?);
    write_dump(merger, &chrom_names, &mut out)?;
    out.flush()?;
    Ok(())
}

/// Stream every emitted `CohortLocus` to `out` as one TSV row, in catalog order.
/// Generic over the merge sources so it is testable in memory.
fn write_dump<R: Read + Seek, C: Read, W: Write>(
    merger: CohortMerger<R, C>,
    chrom_names: &[String],
    out: &mut W,
) -> Result<(), SsrCallError> {
    writeln!(out, "#chrom\tstart\tend\tmotif\tn_present\tsamples")?;
    for item in merger {
        let (_seq, cohort) = item?;
        writeln!(out, "{}", format_locus(&cohort, chrom_names))?;
    }
    Ok(())
}

/// One TSV row for a locus: coordinates + motif + the present samples, each as
/// `idx:depth=…,distinct=…`. Pure, so the formatting is unit-testable on its own.
fn format_locus(cohort: &CohortLocus, chrom_names: &[String]) -> String {
    let chrom = chrom_names
        .get(cohort.locus.chrom_id as usize)
        .map(String::as_str)
        .unwrap_or("?");
    let motif = String::from_utf8_lossy(cohort.motif.as_bytes());
    let samples = cohort
        .present
        .iter()
        .zip(&cohort.samples)
        .map(|(idx, evidence)| {
            format!(
                "{idx}:depth={},distinct={}",
                evidence.qc.depth,
                evidence.seq_counts.len()
            )
        })
        .collect::<Vec<_>>()
        .join(";");
    format!(
        "{chrom}\t{}\t{}\t{motif}\t{}\t{samples}",
        cohort.locus.start,
        cohort.locus.end,
        cohort.present_count(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::test_support::{
        REF_MD5, catalog_bytes, loc, obs, reader, rec, ssr_header, ssr_psp,
    };
    use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};
    use crate::ssr::types::Motif;
    use std::io::Cursor;

    fn cohort_locus(chrom_id: u32, start: u32, present: &[(u32, u32, usize)]) -> CohortLocus {
        let mut cl = CohortLocus::new(
            LocusId {
                chrom_id,
                start,
                end: start + 6,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGGCACACATTTTTT".as_slice()),
            Box::from(b"CACACA".as_slice()),
        );
        for &(idx, depth, n_alleles) in present {
            cl.push(
                idx,
                SampleEvidence {
                    seq_counts: (0..n_alleles)
                        .map(|i| (vec![b'C', b'A', i as u8].into_boxed_slice(), 1))
                        .collect(),
                    qc: SsrQc {
                        depth,
                        n_filtered: 0,
                        mapped_reads: 0,
                        n_low_quality: 0,
                        n_border_off_end: 0,
                    },
                },
            );
        }
        cl
    }

    #[test]
    fn format_locus_renders_chrom_motif_and_samples() {
        let names = vec!["chr1".to_string()];
        let cl = cohort_locus(0, 16, &[(0, 30, 1), (2, 12, 3)]);
        assert_eq!(
            format_locus(&cl, &names),
            "chr1\t16\t22\tCA\t2\t0:depth=30,distinct=1;2:depth=12,distinct=3"
        );
    }

    #[test]
    fn format_locus_falls_back_when_chrom_id_unknown() {
        let cl = cohort_locus(7, 16, &[(0, 30, 1)]);
        assert!(format_locus(&cl, &[]).starts_with("?\t16\t22\tCA\t1\t"));
    }

    #[test]
    fn write_dump_emits_header_and_one_row_per_locus_in_order() {
        let loci = [loc("chr1", 16), loc("chr1", 60), loc("chr1", 100)];
        let catalog =
            crate::ssr::catalog::io::CatalogReader::new(Cursor::new(catalog_bytes(REF_MD5, &loci)))
                .unwrap();
        // Sample 0 covers 16 + 100; sample 1 covers 60 + 100.
        let a = reader(ssr_psp(
            ssr_header(&["chr1"], REF_MD5),
            &[
                rec(0, 16, obs(&[(b"CACACA", 4)])),
                rec(0, 100, obs(&[(b"CA", 2)])),
            ],
        ));
        let b = reader(ssr_psp(
            ssr_header(&["chr1"], REF_MD5),
            &[
                rec(0, 60, obs(&[(b"CACACACA", 5)])),
                rec(0, 100, obs(&[(b"CA", 3)])),
            ],
        ));
        let merger =
            CohortMerger::from_parts(catalog, vec![("A".into(), a), ("B".into(), b)]).unwrap();
        let names = merger.chrom_names();

        let mut out = Vec::new();
        write_dump(merger, &names, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();

        assert_eq!(
            text,
            "#chrom\tstart\tend\tmotif\tn_present\tsamples\n\
             chr1\t16\t22\tCA\t1\t0:depth=30,distinct=1\n\
             chr1\t60\t66\tCA\t1\t1:depth=30,distinct=1\n\
             chr1\t100\t106\tCA\t2\t0:depth=30,distinct=1;1:depth=30,distinct=1\n"
        );
    }

    #[test]
    fn run_reads_files_and_writes_the_dump() {
        let dir = tempfile::TempDir::new().unwrap();
        let catalog_path = dir.path().join("c.ssr.catalog");
        let a_path = dir.path().join("a.ssr.psp");
        let out_path = dir.path().join("out.tsv");

        std::fs::write(&catalog_path, catalog_bytes(REF_MD5, &[loc("chr1", 16)])).unwrap();
        std::fs::write(
            &a_path,
            ssr_psp(
                ssr_header(&["chr1"], REF_MD5),
                &[rec(0, 16, obs(&[(b"CACACA", 4)]))],
            ),
        )
        .unwrap();

        let config = SsrCallConfig {
            catalog: catalog_path,
            psp_files: vec![a_path],
            output: out_path.clone(),
            threads: 4,
            queue_depth: 4,
        };
        run(&config).unwrap();

        let text = std::fs::read_to_string(&out_path).unwrap();
        assert_eq!(
            text,
            "#chrom\tstart\tend\tmotif\tn_present\tsamples\n\
             chr1\t16\t22\tCA\t1\t0:depth=30,distinct=1\n"
        );
    }

    #[test]
    fn run_errors_on_a_missing_catalog() {
        let config = SsrCallConfig {
            catalog: PathBuf::from("/no/such/catalog.ssr.catalog"),
            psp_files: vec![PathBuf::from("/no/such/a.ssr.psp")],
            output: PathBuf::from("/tmp/unused.tsv"),
            threads: 4,
            queue_depth: 4,
        };
        assert!(matches!(run(&config), Err(SsrCallError::Merge(_))));
    }

    // ── H1: build_param_set ──

    use crate::ssr::cohort::param_estimation::{SampleGroupId, StutterLevel, StutterShape};

    fn shape() -> StutterShape {
        StutterShape {
            up_rate: 1.0,
            down_rate: 2.0,
            decay: 0.1,
        }
    }

    /// An `EstimatedParams` with shape for the given periods and a fitted `G₀` only for
    /// those in `g0_periods` (to exercise the Mi1 backfill).
    fn est_for(shape_periods: &[u8], g0_periods: &[u8]) -> EstimatedParams {
        EstimatedParams {
            eps: 0.004,
            shape_by_period: shape_periods.iter().map(|&p| (p, shape())).collect(),
            level_by_sample: HashMap::new(),
            g0_by_period: g0_periods
                .iter()
                .map(|&p| (p, G0PseudocountDecay { p: 0.3 }))
                .collect(),
        }
    }

    /// A `GroupedParams` placing each of `n` samples in `group_of` (one entry per
    /// sample index, or omit an index to make it unresolved), with `n_groups` groups.
    fn grouped_for(group_of: &[(u32, u16)], n_groups: usize) -> GroupedParams {
        GroupedParams {
            group_of_sample: group_of
                .iter()
                .map(|&(s, g)| (s, SampleGroupId(g)))
                .collect(),
            n_groups,
            eps_per_group: vec![0.004; n_groups],
            level_per_group: vec![
                StutterLevel {
                    baseline: 0.06,
                    slope: 0.0,
                };
                n_groups
            ],
            shape_by_group_period: (0..n_groups)
                .map(|g| ((SampleGroupId(g as u16), 2u8), shape()))
                .collect(),
            eps_freeze_justified: true,
        }
    }

    #[test]
    fn build_param_set_maps_every_field_and_densifies_groups() {
        let est = est_for(&[2], &[2]);
        let grouped = grouped_for(&[(0, 0), (1, 1)], 2);
        let params = build_param_set(&est, &grouped, 2, &G0FitCfg::dev_default()).unwrap();

        assert_eq!(
            params.group_of_sample,
            vec![SampleGroupId(0), SampleGroupId(1)]
        );
        assert_eq!(params.error_per_sample_group, vec![PerBaseError(0.004); 2]);
        assert_eq!(params.stutter_shape_parent[&2], shape());
        assert_eq!(
            params.stutter_shape_by_cell[&(SampleGroupId(0), 2)],
            shape()
        );
        assert_eq!(params.level_seed.len(), 2);
        assert_eq!(params.pseudocount_decay_per_loci_group[&2].p, 0.3); // fitted, not fallback
        assert_eq!(params.f0_seed, 0.0); // F is frozen by the burn-in, not here
    }

    #[test]
    fn build_param_set_hard_errors_on_a_sample_with_no_confident_genotype() {
        // Sample 1 of 2 never resolved (absent from group_of_sample) → decision-E error.
        let est = est_for(&[2], &[2]);
        let grouped = grouped_for(&[(0, 0)], 1);
        match build_param_set(&est, &grouped, 2, &G0FitCfg::dev_default()) {
            Err(SsrCallError::UnresolvedSamples { samples }) => assert_eq!(samples, vec![1]),
            other => panic!("expected UnresolvedSamples, got {other:?}"),
        }
    }

    #[test]
    fn build_param_set_carries_a_g0_only_period() {
        // Period 4 is fitted in g0_by_period but absent from shape_by_period; the
        // clone-first base of the backfill must still carry its fitted value.
        let est = est_for(&[2], &[2, 4]);
        let grouped = grouped_for(&[(0, 0), (1, 0)], 1);
        let params = build_param_set(&est, &grouped, 2, &G0FitCfg::dev_default()).unwrap();
        assert_eq!(params.pseudocount_decay_per_loci_group[&4].p, 0.3);
    }

    #[test]
    fn build_param_set_on_an_empty_cohort_is_ok_and_empty() {
        // Degenerate (the merger rejects an empty cohort upstream): no samples → no
        // unresolved samples → an empty, valid ParamSet.
        let est = est_for(&[2], &[2]);
        let grouped = grouped_for(&[], 0);
        let params = build_param_set(&est, &grouped, 0, &G0FitCfg::dev_default()).unwrap();
        assert!(params.group_of_sample.is_empty());
    }

    #[test]
    fn build_param_set_backfills_g0_fallback_for_a_period_without_a_fit() {
        // Period 3 has chemistry (shape) but no fitted G₀ → it must take fallback_p,
        // so the configured fallback is the single source of truth (Mi1).
        let est = est_for(&[2, 3], &[2]);
        let grouped = grouped_for(&[(0, 0), (1, 0)], 1);
        let cfg = G0FitCfg::dev_default();
        let params = build_param_set(&est, &grouped, 2, &cfg).unwrap();
        assert_eq!(params.pseudocount_decay_per_loci_group[&2].p, 0.3); // fitted
        assert_eq!(
            params.pseudocount_decay_per_loci_group[&3].p,
            cfg.fallback_p // backfilled
        );
    }
}
