//! The `ssr-call` driver — the two-pass streaming pipeline that turns a cohort of
//! `.ssr.psp` files into a VCF (arch `doc/devel/architecture/ssr_call_driver.md`).
//!
//! **Pass 1 (burn-in):** open the cohort, collect a bounded subset of loci, and freeze
//! the cross-locus-pooled parameters — chemistry (ε, stutter shape/level, `G₀`) plus the
//! per-individual inbreeding `F` (arch §4). **Pass 2 (genotyping sweep):** re-open and
//! stream every locus, genotyping each independently on the frozen parameters and
//! emitting a VCF record (emit-iff-variable for PASS, filtered loci kept with their
//! reason — arch §6). The sweep is single-threaded for now; the bounded-queue worker
//! pool is Milestone J. `config.threads` drives the (parallel, byte-identical) burn-in.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Read, Seek, Write};
use std::path::PathBuf;

use crate::ssr::cohort::candidate_set::{Admission, CandidateCfg, assemble_candidates};
use crate::ssr::cohort::em::{EmCfg, run_locus_em_with};
use crate::ssr::cohort::em_init::seed_locus;
use crate::ssr::cohort::inbreeding::{OuterCfg, run_cohort_em};
use crate::ssr::cohort::merge::{CohortMerger, SsrMergeError};
use crate::ssr::cohort::param_estimation::{
    G0FitCfg, G0PseudocountDecay, ParamSet, PerBaseError, StutterLevel,
};
use crate::ssr::cohort::prepass::{EstimatedParams, estimate, run_prepass_stats};
use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
use crate::ssr::cohort::sample_groups::{ClusterCfg, GroupedParams, group_samples};
use crate::ssr::cohort::types::CohortLocus;
use crate::ssr::cohort::vcf_out::{
    FpControlCfg, apply_fp_control, f_is_warning, format_vcf_record, is_variable, site_qual,
    write_vcf_header,
};

/// Inputs for an `ssr-call` run.
pub(crate) struct SsrCallConfig {
    /// The shared `.ssr.catalog`.
    pub(crate) catalog: PathBuf,
    /// The per-sample `.ssr.psp` evidence files.
    pub(crate) psp_files: Vec<PathBuf>,
    /// Where the cohort VCF is written.
    pub(crate) output: PathBuf,
    /// Thread count for the parallel, byte-identical burn-in (the genotyping sweep is
    /// single-threaded until Milestone J).
    pub(crate) threads: usize,
    /// Reserved — the bounded-queue depth of the streaming genotyping sweep (Milestone
    /// J); the current single-threaded sweep does not use it.
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
    /// Two cohort samples reduce to the same VCF sample name — duplicate `#CHROM`
    /// columns would be an invalid VCF. Sample names are the input path basenames
    /// (decision C-1), so two inputs with the same file name collide (H2 review Mi1).
    #[error("duplicate sample name {name:?} — cohort inputs must have distinct file names")]
    DuplicateSampleName { name: String },
    /// Building the burn-in thread pool failed.
    #[error("building the ssr-call thread pool")]
    ThreadPool(#[from] rayon::ThreadPoolBuildError),
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

/// Diploid is the only supported ploidy (decision D); `run_locus_em` asserts it.
const PLOIDY: u8 = 2;

/// Cap on the number of loci the burn-in holds in RAM to estimate and freeze the
/// cross-locus-pooled parameters (arch §4). The genotyping sweep streams *all* loci
/// regardless; only the burn-in is bounded. Loci are taken in catalog-stream order until
/// the cap — a positionally-biased first-`cap` subset; a stratified / reservoir selection
/// is a calibration follow-up (reading Q-R5), not needed for correctness.
const BURN_IN_MAX_LOCI: usize = 20_000;

/// The cross-locus-pooled parameters the burn-in freezes for the genotyping sweep.
struct FrozenParams {
    /// Frozen chemistry (ε, stutter shape parent, per-group level seed, `G₀`, groups).
    params: ParamSet,
    /// Frozen per-(global)-sample inbreeding `F`.
    f_per_sample: Vec<f64>,
    /// Frozen per-group stutter level.
    level_per_group: Vec<StutterLevel>,
}

/// Run `ssr-call`: burn-in over a bounded subset to freeze the pooled parameters, then a
/// streaming genotyping sweep over every locus → VCF at `config.output` (arch §2/§4/§6).
pub(crate) fn run(config: &SsrCallConfig) -> Result<(), SsrCallError> {
    let rung_cfg = RungCfg::dev_default();
    let cand_cfg = CandidateCfg::dev_default();
    let em_cfg = EmCfg::dev_default();
    let outer_cfg = OuterCfg::dev_default();
    let cluster_cfg = ClusterCfg::dev_default();
    let g0_cfg = G0FitCfg::dev_default();
    let fp_cfg = FpControlCfg::dev_default();

    // ── Pass 1: header table + bounded burn-in subset ──
    let merger = CohortMerger::open(&config.catalog, &config.psp_files)?;
    let chromosomes = merger.chromosomes();
    let chrom_names = merger.chrom_names();
    let sample_names = merger.sample_names();
    let n_samples = sample_names.len();
    check_unique_sample_names(&sample_names)?;
    let subset = collect_burn_in_subset(merger)?;

    // Estimate chemistry and freeze `F` + the per-group level on the requested thread
    // pool. The reduces are integer / fixed-point, so the result is byte-identical
    // regardless of thread count (arch §4).
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads.max(1))
        .build()?;
    let frozen = pool.install(|| -> Result<FrozenParams, SsrCallError> {
        let stats = run_prepass_stats(&subset, PLOIDY, &rung_cfg);
        let est = estimate(&stats, &g0_cfg);
        let grouped = group_samples(&stats, &est, &cluster_cfg);
        let params = build_param_set(&est, &grouped, n_samples, &g0_cfg)?;
        let calls = run_cohort_em(
            &subset,
            &params,
            grouped.level_per_group.clone(),
            PLOIDY,
            &em_cfg,
            &rung_cfg,
            &cand_cfg,
            &outer_cfg,
        );
        Ok(FrozenParams {
            params,
            f_per_sample: calls.f_per_sample,
            level_per_group: calls.level_per_group,
        })
    })?;
    drop(subset);

    // Cohort warnings (header + stderr): the apparent-`F_IS` notice (arch §7).
    let mut warnings = Vec::new();
    if let Some(warning) = f_is_warning(&frozen.f_per_sample, &fp_cfg) {
        eprintln!("ssr-call warning: {warning}");
        warnings.push(warning);
    }

    // ── Pass 2: re-open and stream every locus → VCF ──
    let merger = CohortMerger::open(&config.catalog, &config.psp_files)?;
    let mut out = BufWriter::new(File::create(&config.output)?);
    write_vcf_header(&mut out, &chromosomes, &sample_names, &warnings)?;
    for item in merger {
        let (_seq, locus) = item?;
        if let Some(line) = genotype_locus(
            &locus,
            &chrom_names,
            &frozen,
            &rung_cfg,
            &cand_cfg,
            &em_cfg,
            &fp_cfg,
        ) {
            writeln!(out, "{line}")?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Reject a cohort whose sample names are not all distinct (H2 review Mi1): duplicate
/// `#CHROM` columns would be an invalid VCF.
fn check_unique_sample_names(names: &[String]) -> Result<(), SsrCallError> {
    let mut seen = HashSet::with_capacity(names.len());
    for name in names {
        if !seen.insert(name.as_str()) {
            return Err(SsrCallError::DuplicateSampleName { name: name.clone() });
        }
    }
    Ok(())
}

/// Stream the merger and collect up to [`BURN_IN_MAX_LOCI`] loci for the burn-in
/// (bounded memory; the genotyping pass re-reads all loci).
fn collect_burn_in_subset<R: Read + Seek, C: Read>(
    merger: CohortMerger<R, C>,
) -> Result<Vec<CohortLocus>, SsrCallError> {
    let mut subset = Vec::new();
    for item in merger {
        let (_seq, locus) = item?;
        subset.push(locus);
        if subset.len() >= BURN_IN_MAX_LOCI {
            break;
        }
    }
    Ok(subset)
}

/// Genotype one streamed locus on the frozen parameters and format its VCF line, applying
/// the emit policy (arch §6): a filtered locus is emitted with its reason; a PASS locus is
/// emitted iff it is variable; a monomorphic PASS locus is dropped (`None`).
fn genotype_locus(
    locus: &CohortLocus,
    chrom_names: &[String],
    frozen: &FrozenParams,
    rung_cfg: &RungCfg,
    cand_cfg: &CandidateCfg,
    em_cfg: &EmCfg,
    fp_cfg: &FpControlCfg,
) -> Option<String> {
    let rungs = build_rungs(locus, rung_cfg);
    let candidates = assemble_candidates(locus, &rungs, PLOIDY, cand_cfg);
    let seed = seed_locus(
        locus,
        &rungs,
        &candidates,
        &frozen.params,
        PLOIDY,
        rung_cfg.prominence,
    );
    let f_present: Vec<f64> = locus
        .present
        .iter()
        .map(|&g| frozen.f_per_sample[g as usize])
        .collect();
    let mut call = run_locus_em_with(
        locus,
        &rungs,
        &candidates,
        &frozen.params,
        &seed,
        PLOIDY,
        em_cfg,
        &f_present,
        &frozen.level_per_group,
    );
    apply_fp_control(locus, &mut call, fp_cfg);

    // PASS + monomorphic → drop; PASS + variable or any filtered verdict → emit.
    if candidates.admit == Admission::Pass && !is_variable(&call, &candidates) {
        return None;
    }
    let qual = site_qual(&call, &candidates, fp_cfg);
    let chrom = chrom_names
        .get(locus.locus.chrom_id as usize)
        .map(String::as_str)
        .unwrap_or("?");
    Some(format_vcf_record(chrom, locus, &candidates, &call, qual))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::psp::registry_ssr::SsrLocusRecord;
    use crate::ssr::cohort::test_support::{REF_MD5, catalog_bytes, ssr_header, ssr_psp};
    use crate::ssr::types::{Locus, Motif};

    /// `CA` repeated `units` times — one allele's tract sequence.
    fn ca(units: u16) -> Box<[u8]> {
        "CA".repeat(units as usize).into_bytes().into_boxed_slice()
    }

    /// A catalog locus on `chrom` at 0-based `start` with a `units`-unit `CA` reference
    /// tract (`GGGGGG` + `CA`×units + `TTTTTT`, anchored 6 bp upstream).
    fn loc_units(chrom: &str, start: u32, units: u32) -> Locus {
        let mut ref_seq = b"GGGGGG".to_vec();
        ref_seq.extend_from_slice(&"CA".repeat(units as usize).into_bytes());
        ref_seq.extend_from_slice(b"TTTTTT");
        Locus::new(
            chrom.into(),
            start,
            start + 2 * units,
            Motif::new(b"CA").unwrap(),
            1.0,
            ref_seq.into_boxed_slice(),
            start - 6,
        )
        .unwrap()
    }

    /// Per-read evidence for a genotype: a clean peak at each allele length plus light
    /// ±1 stutter (a homozygote gets full depth; a het splits it). An empty allele list
    /// yields a thin record that never resolves confidently.
    fn reads_for(alleles: &[u16]) -> Vec<(Box<[u8]>, u32)> {
        if alleles.is_empty() {
            return vec![(ca(8), 3)]; // below the resolution depth floor
        }
        let peak = if alleles.len() == 1 { 150 } else { 75 };
        let mut counts: Vec<(Box<[u8]>, u32)> = Vec::new();
        for &a in alleles {
            counts.push((ca(a), peak));
            counts.push((ca(a - 1), 4));
            counts.push((ca(a + 1), 4));
        }
        counts.sort_by(|x, y| x.0.cmp(&y.0));
        counts
    }

    /// One `.ssr.psp` record at the `units`-unit locus frame starting at 0-based `start`
    /// (1-based half-open container coords; `end` matches the catalog frame, *not* the
    /// observed allele lengths — the cursor derives the `LocusId` from these). Depth is
    /// the observed read total.
    fn ssr_record(start: u32, units: u32, observed: Vec<(Box<[u8]>, u32)>) -> SsrLocusRecord {
        let depth = observed.iter().map(|(_, c)| *c).sum();
        SsrLocusRecord {
            chrom_id: 0,
            start: start + 1,
            end: start + 2 * units + 1,
            depth,
            n_filtered: 0,
            mapped_reads: depth,
            n_low_quality: 0,
            n_border_off_end: 0,
            observed,
        }
    }

    /// Write a cohort to `dir`: a catalog over `loci` (`(start, ref_units)`) plus one
    /// `.ssr.psp` per sample (`(name, per-locus allele genotype)`), and return the config.
    fn write_cohort(
        dir: &std::path::Path,
        loci: &[(u32, u32)],
        samples: &[(&str, Vec<Vec<u16>>)],
    ) -> SsrCallConfig {
        let catalog_loci: Vec<Locus> = loci.iter().map(|&(s, u)| loc_units("chr1", s, u)).collect();
        let catalog = dir.join("c.ssr.catalog");
        std::fs::write(&catalog, catalog_bytes(REF_MD5, &catalog_loci)).unwrap();
        let mut psp_files = Vec::new();
        for (name, genotypes) in samples {
            let records: Vec<SsrLocusRecord> = loci
                .iter()
                .zip(genotypes)
                .map(|(&(s, u), alleles)| ssr_record(s, u, reads_for(alleles)))
                .collect();
            let path = dir.join(format!("{name}.ssr.psp"));
            std::fs::write(&path, ssr_psp(ssr_header(&["chr1"], REF_MD5), &records)).unwrap();
            psp_files.push(path);
        }
        SsrCallConfig {
            catalog,
            psp_files,
            output: dir.join("out.vcf"),
            threads: 2,
            queue_depth: 4,
        }
    }

    #[test]
    fn run_emits_a_vcf_with_a_pass_variant_and_drops_a_monomorphic_locus() {
        let dir = tempfile::TempDir::new().unwrap();
        // Locus A (start 40): 6 homozygous-ref (8/8) + 6 separated hets (6/10) → variable.
        // Locus B (start 100): every sample 8/8 (= reference) → monomorphic, dropped.
        let mut owned: Vec<(String, Vec<Vec<u16>>)> = Vec::new();
        for i in 0..6 {
            owned.push((format!("hom{i}"), vec![vec![8], vec![8]]));
        }
        for i in 0..6 {
            owned.push((format!("het{i}"), vec![vec![6, 10], vec![8]]));
        }
        let samples: Vec<(&str, Vec<Vec<u16>>)> =
            owned.iter().map(|(n, g)| (n.as_str(), g.clone())).collect();
        let config = write_cohort(dir.path(), &[(40, 8), (100, 8)], &samples);

        run(&config).unwrap();
        let vcf = std::fs::read_to_string(&config.output).unwrap();

        assert!(vcf.starts_with("##fileformat=VCFv4.4"), "{vcf}");
        let chrom_line = vcf.lines().find(|l| l.starts_with("#CHROM")).unwrap();
        assert_eq!(chrom_line.split('\t').count(), 9 + 12, "12 sample columns");

        // Exactly one data record (locus A); the monomorphic locus B is dropped.
        let records: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(
            records.len(),
            1,
            "monomorphic locus B must be dropped\n{vcf}"
        );
        let cols: Vec<&str> = records[0].split('\t').collect();
        assert_eq!(cols[0], "chr1");
        assert_eq!(cols[6], "PASS", "variant locus should PASS: {}", records[0]);
        assert_ne!(cols[4], ".", "should carry an ALT allele");
        // Sample columns are homs (9..15) then hets (15..21); each het calls REPCN 6/10.
        for col in &cols[15..21] {
            let repcn = col.split(':').nth(2).unwrap();
            let mut units: Vec<u16> = repcn.split(',').map(|s| s.parse().unwrap()).collect();
            units.sort_unstable();
            assert_eq!(units, vec![6, 10], "het miscalled: {col}");
        }
    }

    #[test]
    fn run_hard_errors_when_a_sample_never_resolves() {
        let dir = tempfile::TempDir::new().unwrap();
        // Four resolvable homozygotes + one thin sample (index 4) that never resolves.
        let mut owned: Vec<(String, Vec<Vec<u16>>)> = (0..4)
            .map(|i| (format!("ok{i}"), vec![vec![8u16]]))
            .collect();
        owned.push(("thin".to_string(), vec![vec![]]));
        let samples: Vec<(&str, Vec<Vec<u16>>)> =
            owned.iter().map(|(n, g)| (n.as_str(), g.clone())).collect();
        let config = write_cohort(dir.path(), &[(40, 8)], &samples);

        match run(&config) {
            Err(SsrCallError::UnresolvedSamples { samples }) => assert_eq!(samples, vec![4]),
            other => panic!("expected UnresolvedSamples, got {other:?}"),
        }
    }

    #[test]
    fn check_unique_sample_names_rejects_a_collision() {
        let names = vec!["s".to_string(), "t".to_string(), "s".to_string()];
        match check_unique_sample_names(&names) {
            Err(SsrCallError::DuplicateSampleName { name }) => assert_eq!(name, "s"),
            other => panic!("expected DuplicateSampleName, got {other:?}"),
        }
        assert!(check_unique_sample_names(&["a".to_string(), "b".to_string()]).is_ok());
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

    use crate::ssr::cohort::param_estimation::{SampleGroupId, StutterShape};

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
