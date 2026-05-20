//! Integration tests for the cohort CLI subcommands: `var-calling`,
//! `estimate-contamination`, and `var-calling-from-bam`.
//!
//! Each test builds a tiny synthetic FASTA + CRAM(s), drives the
//! relevant `run_*` orchestrator(s), and verifies the resulting
//! artefacts. Shared FASTA / CRAM / SAM-header / read fixture
//! helpers live in `tests/common/mod.rs` (Mi20).

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::{WRONG_MD5, build_cram, build_fasta, fixture_md5, read_record};
use pop_var_caller::per_sample_pileup::baq::{
    DEFAULT_BAQ_CHUNK_SIZE, SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH, SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
    SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
};
use pop_var_caller::per_sample_pileup::cram_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MISMATCH_BQ_FLOOR,
};
use pop_var_caller::per_sample_pileup::pileup::{
    DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_ACTIVE_READS, DEFAULT_MAX_INDEL_COLUMN_DEPTH,
    DEFAULT_MAX_RECORD_SPAN, DEFAULT_MAX_SNP_COLUMN_DEPTH,
};
use pop_var_caller::pop_var_caller::cli::shared_args::{CohortPipelineArgs, Stage1Args};
use pop_var_caller::pop_var_caller::cohort_driver::{
    DEFAULT_MIN_ALT_OBS_PER_SAMPLE, DEFAULT_MIN_QUAL_PHRED,
};
use pop_var_caller::pop_var_caller::estimate_contamination::{
    EstimateContaminationArgs, EstimateContaminationCliError, run_estimate_contamination,
};
use pop_var_caller::pop_var_caller::var_calling::{
    VarCallingArgs, VarCallingCliError, run_var_calling,
};
use pop_var_caller::pop_var_caller::var_calling_from_bam::{
    VarCallingFromBamArgs, VarCallingFromBamCliError, run_var_calling_from_bam,
};
use pop_var_caller::pop_var_caller::{PileupArgs, run_pileup};
use pop_var_caller::var_calling::contamination_estimation::{
    DEFAULT_C_S_INIT, DEFAULT_Q_B_INIT_PER_CLASS,
};
use pop_var_caller::var_calling::dust_filter::{
    DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW,
};
use pop_var_caller::var_calling::per_group_merger::{
    DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY,
};
use pop_var_caller::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT,
};
use pop_var_caller::var_calling::variant_grouping::DEFAULT_MAX_VARIANT_GROUP_SPAN;
use noodles_sam::alignment::record_buf::RecordBuf;
use tempfile::TempDir;

// ---------------------------------------------------------------------
// Cohort-CLI-specific fixture builders (FASTA / CRAM / SAM-header /
// read helpers come from `tests/common/mod.rs`; the helpers below
// build cohort args + drive `run_pileup` to materialise per-sample
// `.psp` files).
// ---------------------------------------------------------------------

fn pileup_args(reference: PathBuf, output: PathBuf, crams: Vec<PathBuf>) -> PileupArgs {
    PileupArgs {
        reference,
        output,
        threads: None,
        crams,
        stage1: Stage1Args {
            min_mapq: DEFAULT_MIN_MAPQ,
            no_baq: true, // skip BAQ — tiny synthetic reads don't need it
            min_read_length: 0,
            keep_qc_fail: false,
            keep_duplicates: false,
            max_read_mismatch_fraction: DEFAULT_MAX_READ_MISMATCH_FRACTION,
            mismatch_bq_floor: DEFAULT_MISMATCH_BQ_FLOOR,
            baq_gap_open_prob: SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
            baq_gap_extend_prob: SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
            baq_band_half_width: SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
            baq_chunk_size: DEFAULT_BAQ_CHUNK_SIZE,
            max_snp_column_depth: DEFAULT_MAX_SNP_COLUMN_DEPTH,
            max_indel_column_depth: DEFAULT_MAX_INDEL_COLUMN_DEPTH,
            max_record_span: DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
            max_active_reads: DEFAULT_MAX_ACTIVE_READS,
        },
    }
}

fn var_calling_args(
    reference: PathBuf,
    output: PathBuf,
    psp_files: Vec<PathBuf>,
    contamination_estimates: Option<PathBuf>,
) -> VarCallingArgs {
    VarCallingArgs {
        reference,
        output,
        threads: None,
        contamination_estimates,
        no_complexity_filter: true, // tiny ref; sdust would mask everything
        psp_files,
        cohort: CohortPipelineArgs {
            ploidy: DEFAULT_PLOIDY,
            complexity_window: DEFAULT_DUST_WINDOW,
            complexity_threshold: DEFAULT_DUST_THRESHOLD,
            var_group_max_span: DEFAULT_MAX_VARIANT_GROUP_SPAN,
            max_alleles_per_var: DEFAULT_MAX_ALLELES_PER_RECORD,
            inbreeding_coefficient: DEFAULT_INBREEDING_COEFFICIENT,
            em_convergence_threshold: DEFAULT_CONVERGENCE_THRESHOLD,
            em_max_iterations: DEFAULT_MAX_ITERATIONS,
            ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
            snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
            indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
            compound_alt_pseudocount: DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
            max_gq_phred: DEFAULT_MAX_GQ_PHRED,
            min_qual_phred: DEFAULT_MIN_QUAL_PHRED,
            min_alt_obs_per_sample: DEFAULT_MIN_ALT_OBS_PER_SAMPLE,
            emit_gp: false,
        },
    }
}

/// Build a `.psp` for one sample by running the real pileup
/// orchestrator on a synthetic CRAM. Returns the .psp path.
fn make_psp_for_sample(
    dir: &Path,
    fasta_path: &Path,
    sample: &str,
    reads: &[RecordBuf],
) -> PathBuf {
    let cram = build_cram(dir, fasta_path, sample, Some(fixture_md5()), reads);
    let psp = dir.join(format!("{sample}.psp"));
    let args = pileup_args(fasta_path.to_path_buf(), psp.clone(), vec![cram]);
    run_pileup(&args).expect("run_pileup OK");
    psp
}

/// Build an `EstimateContaminationArgs` parameterised for tests:
/// relaxed thresholds (block_size = 1, min_depth = 1, stability
/// tolerance = 0.5, stability_blocks = 1) so the side-pass
/// terminates quickly on the tiny synthetic fixtures used in
/// integration tests.
fn estimate_contamination_args(
    reference: PathBuf,
    output: PathBuf,
    psp_files: Vec<PathBuf>,
) -> EstimateContaminationArgs {
    EstimateContaminationArgs {
        reference,
        output,
        threads: None,
        batch_assignment: None,
        psp_files,
        // Relaxed convergence: any one block of stability ends the run.
        // 0.1 is the engine's hard upper bound (STABILITY_TOLERANCE_RANGE_MAX).
        stability_tolerance: 0.1,
        stability_blocks: 1,
        // Tiny blocks; admit any covered site.
        block_size: 1,
        min_depth: 1,
        // Engine requires `> 0.5`; relaxed for the synthetic
        // fixture so depth-1 minor reads alongside depth-3 major
        // reads still qualify the site as hom-major (3/4 = 0.75).
        min_major_fraction: 0.51,
        min_cohort_minor_count: 1,
        min_cohort_minor_fraction: 0.0,
        min_batch_size_for_contamination: 2,
        ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
        snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
        indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
        c_s_init: DEFAULT_C_S_INIT,
        q_b_init_per_class: DEFAULT_Q_B_INIT_PER_CLASS,
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

/// **End-to-end happy path.** Three samples → three `.psp` files →
/// run_var_calling → VCF on disk. Verifies that the whole cohort
/// pipeline glue (PspReader → merger → DustFilter [bypassed] →
/// grouper → per-group merger → posterior engine → VCF writer) wires
/// together without errors and produces a parseable VCF.
#[test]
fn var_calling_happy_path_three_samples() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());

    // Three samples, identical-ish reads covering positions 1..14.
    // Reads are mostly REF (all-A reference); some bases are altered
    // to introduce SNP variants that the grouper will pick up.
    let reads_a = vec![
        read_record("r1", 10, b"AAAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];
    let reads_b = vec![
        read_record("r1", 10, b"ACAAA"), // SNP at pos 11
        read_record("r2", 15, b"AAAAA"),
    ];
    let reads_c = vec![
        read_record("r1", 10, b"AAAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];

    let psp_a = make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads_a);
    let psp_b = make_psp_for_sample(dir.path(), &fasta, "NA00002", &reads_b);
    let psp_c = make_psp_for_sample(dir.path(), &fasta, "NA00003", &reads_c);

    let vcf_out = dir.path().join("cohort.vcf");
    let args = var_calling_args(fasta, vcf_out.clone(), vec![psp_a, psp_b, psp_c], None);
    run_var_calling(&args).expect("run_var_calling OK");

    assert!(vcf_out.exists(), "VCF must exist at the final path");
    let body = fs::read_to_string(&vcf_out).expect("read VCF");
    // Self-checks on the header — we expect a multi-sample VCF.
    assert!(body.contains("##fileformat=VCFv4"), "VCF header missing");
    assert!(body.contains("##source="), "##source missing");
    assert!(
        body.contains("\tNA00001\tNA00002\tNA00003"),
        "sample columns missing"
    );
}

/// **Determinism under the per-chromosome parallel path.** Run the
/// same cohort input through `run_var_calling` twice (different
/// output paths in the same process — the rayon pool is shared, so
/// both runs use the same parallelism). The data body of the
/// resulting VCFs must be byte-identical: per-chrom fragments are
/// assembled in contig-table order regardless of worker finish
/// order, so the final VCF is a deterministic function of the
/// inputs. Plan §test plan #1 (morphed: post-merge the serial
/// path no longer exists, so we assert run-vs-run determinism
/// instead of parallel-vs-serial equivalence).
#[test]
fn var_calling_emits_deterministic_vcf_across_runs() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());

    let reads_a = vec![
        read_record("r1", 10, b"AAAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];
    let reads_b = vec![
        read_record("r1", 10, b"ACAAA"), // SNP at pos 11
        read_record("r2", 15, b"AAAAA"),
    ];
    let psp_a = make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads_a);
    let psp_b = make_psp_for_sample(dir.path(), &fasta, "NA00002", &reads_b);

    let psps = vec![psp_a, psp_b];
    let vcf1 = dir.path().join("cohort_run1.vcf");
    let vcf2 = dir.path().join("cohort_run2.vcf");

    let args1 = var_calling_args(fasta.clone(), vcf1.clone(), psps.clone(), None);
    run_var_calling(&args1).expect("first run OK");

    let args2 = var_calling_args(fasta, vcf2.clone(), psps, None);
    run_var_calling(&args2).expect("second run OK");

    // Bodies must match byte-for-byte once we strip header lines
    // that legitimately differ between runs (##source / ##commandline
    // — the latter includes argv, which would only differ if the
    // tests passed different args, which they don't here, but we
    // strip both to be robust against any future header drift).
    let body1 = fs::read_to_string(&vcf1).expect("read vcf1");
    let body2 = fs::read_to_string(&vcf2).expect("read vcf2");
    let strip_volatile = |s: &str| -> String {
        s.lines()
            .filter(|l| !l.starts_with("##source=") && !l.starts_with("##commandline="))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(
        strip_volatile(&body1),
        strip_volatile(&body2),
        "two runs over identical input must produce byte-identical VCF body \
         (parallel path is deterministic; per-chrom fragments assemble in \
         contig-table order regardless of worker finish order)"
    );
}

/// **var-calling-from-bam happy path.** One CRAM directly → VCF.
/// No `.psp` is written by the test (the helper internally writes
/// the .psp for the cohort path, but `run_var_calling_from_bam`
/// goes BAM → VCF in a single process per the
/// `feedback_no_silent_intermediates` rule).
#[test]
fn var_calling_from_bam_happy_path() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let reads = vec![
        read_record("r1", 10, b"ACAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];
    let cram = build_cram(dir.path(), &fasta, "NA12878", Some(fixture_md5()), &reads);
    let vcf_out = dir.path().join("solo.vcf");

    let args = VarCallingFromBamArgs {
        reference: fasta.clone(),
        output: vcf_out.clone(),
        threads: None,
        crams: vec![cram],
        no_complexity_filter: true,
        stage1: Stage1Args {
            min_mapq: DEFAULT_MIN_MAPQ,
            no_baq: true,
            min_read_length: 0,
            keep_qc_fail: false,
            keep_duplicates: false,
            max_read_mismatch_fraction: DEFAULT_MAX_READ_MISMATCH_FRACTION,
            mismatch_bq_floor: DEFAULT_MISMATCH_BQ_FLOOR,
            baq_gap_open_prob: SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
            baq_gap_extend_prob: SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
            baq_band_half_width: SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
            baq_chunk_size: DEFAULT_BAQ_CHUNK_SIZE,
            max_snp_column_depth: DEFAULT_MAX_SNP_COLUMN_DEPTH,
            max_indel_column_depth: DEFAULT_MAX_INDEL_COLUMN_DEPTH,
            max_record_span: DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
            max_active_reads: DEFAULT_MAX_ACTIVE_READS,
        },
        cohort: CohortPipelineArgs {
            ploidy: DEFAULT_PLOIDY,
            complexity_window: DEFAULT_DUST_WINDOW,
            complexity_threshold: DEFAULT_DUST_THRESHOLD,
            var_group_max_span: DEFAULT_MAX_VARIANT_GROUP_SPAN,
            max_alleles_per_var: DEFAULT_MAX_ALLELES_PER_RECORD,
            inbreeding_coefficient: DEFAULT_INBREEDING_COEFFICIENT,
            em_convergence_threshold: DEFAULT_CONVERGENCE_THRESHOLD,
            em_max_iterations: DEFAULT_MAX_ITERATIONS,
            ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
            snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
            indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
            compound_alt_pseudocount: DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
            max_gq_phred: DEFAULT_MAX_GQ_PHRED,
            min_qual_phred: DEFAULT_MIN_QUAL_PHRED,
            min_alt_obs_per_sample: DEFAULT_MIN_ALT_OBS_PER_SAMPLE,
            emit_gp: false,
        },
    };
    run_var_calling_from_bam(&args).expect("run_var_calling_from_bam OK");

    assert!(vcf_out.exists(), "VCF must exist at the final path");
    let body = fs::read_to_string(&vcf_out).expect("read VCF");
    assert!(body.contains("##fileformat=VCFv4"), "VCF header missing");
    // Single-sample VCF — only NA12878 in the column header.
    assert!(body.contains("\tNA12878\n") || body.contains("\tNA12878\r\n"));
    // No other sample names should appear in the data-header line.
    assert!(!body.contains("NA00001"));
}

/// **Sample-name reconciliation: missing sample is a hard error.**
/// Build a `.psp` for `NA00001` but hand var-calling a contamination
/// artefact that only references `OTHER_SAMPLE`. `run_var_calling`
/// must surface `ContaminationSampleMissing`.
#[test]
fn var_calling_rejects_contamination_artefact_missing_sample() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let reads = vec![read_record("r1", 10, b"AAAAA")];
    let psp = make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads);

    // Hand-craft a contamination artefact missing NA00001. The
    // on-disk TOML is the realistic input path here: `run_var_calling`
    // always loads via `ContaminationArtefact::read`, and writing a
    // literal mirrors what a user editing the file by hand would
    // produce. Avoids depending on the (`#[non_exhaustive]`) artefact
    // struct surface from outside the crate.
    let artefact_toml = r#"[provenance]
tool = "pop_var_caller"
version = "0.1.0"
subcommand = "estimate-contamination"
created = 2026-05-19T12:00:00Z

[provenance.inputs]
reference = "ref.fa"
input_psps = ["OTHER_SAMPLE.psp"]

[parameters]

[[batches]]
id = "lane_1"
contaminant_ref_prob = 0.999
contaminant_snp_alt_prob = 0.0008
contaminant_indel_alt_prob = 0.0002

[[samples]]
name = "OTHER_SAMPLE"
batch = "lane_1"
contamination_fraction = 0.01
"#;
    let artefact_path = dir.path().join("estimates.toml");
    fs::write(&artefact_path, artefact_toml).expect("write artefact");

    let vcf_out = dir.path().join("out.vcf");
    let args = var_calling_args(fasta, vcf_out, vec![psp], Some(artefact_path));
    let err = run_var_calling(&args).expect_err("must error on missing sample");
    let chain = format!("{err}");
    let nested = match &err {
        VarCallingCliError::ContamArtefact(inner) => format!("{inner}"),
        _ => chain.clone(),
    };
    assert!(
        nested.contains("NA00001") || chain.contains("NA00001"),
        "error should mention the missing sample: {chain} / {nested}"
    );
}

/// **M13 multi-invocation pattern.** Two consecutive `run_pileup`
/// calls in the same process must both succeed — the cohort CLI's
/// rayon-pool gate
/// ([`configure_rayon_pool`](pop_var_caller::pop_var_caller))
/// is idempotent, so library consumers and a future
/// multi-subcommand test runner can chain `run_*` helpers freely.
///
/// `args.threads = None` so this exercises the multi-invocation
/// *pattern* without depending on rayon's process-global state
/// being uninitialised by an earlier test (which would be racy
/// under cargo-test's shared-process model). The
/// **Mi23 load-bearing chain test.** Build 2 `.psp`s → run
/// `run_estimate_contamination` to a TOML → run `run_var_calling
/// --contamination-estimates <that>` → assert the resulting VCF
/// parses. Exercises the full estimate-contamination → var-calling
/// integration that the cohort CLI's two subcommands are designed
/// to support (this chain is the **contract between them**; until
/// Wave 5 it lived only in the spec, not in the test suite).
///
/// The synthetic 200-base all-A fixture is far too small for the
/// contamination side-pass to recover meaningful estimates; this
/// test only asserts that the chain runs end-to-end and produces
/// a parseable VCF — the side-pass's accuracy is covered by the
/// per-module `tests/contamination_estimation_integration.rs` with
/// realistic synthetic data.
#[test]
fn estimate_contamination_then_var_calling_chain() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());

    // Two samples, each with three 20bp reads covering positions
    // 10..=29 — all-A reference. Sample B adds five reads each
    // with a single C at a different relative offset (positions
    // 10, 11, 12, 13, 14) — 1/20 = 5% mismatch, under
    // `DEFAULT_MAX_READ_MISMATCH_FRACTION = 0.10`. Each of positions
    // 10..=14 then has the cohort-minor C: sample A hom-A
    // (3/3 = 1.0), sample B hom-major-A (7/8 = 0.875
    // ≥ `min_major_fraction = 0.51`), cohort minor count = 1.
    // Five informative sites — enough for the side-pass's
    // Convergence mode to find a stable block.
    let long_ref = b"AAAAAAAAAAAAAAAAAAAA"; // 20 As
    let reads_a = vec![
        read_record("ra1", 10, long_ref),
        read_record("ra2", 10, long_ref),
        read_record("ra3", 10, long_ref),
    ];
    let reads_b = vec![
        read_record("rb1", 10, long_ref),
        read_record("rb2", 10, long_ref),
        read_record("rb3", 10, long_ref),
        read_record("rb4", 10, b"CAAAAAAAAAAAAAAAAAAA"), // C at pos 10
        read_record("rb5", 10, b"ACAAAAAAAAAAAAAAAAAA"), // C at pos 11
        read_record("rb6", 10, b"AACAAAAAAAAAAAAAAAAA"), // C at pos 12
        read_record("rb7", 10, b"AAACAAAAAAAAAAAAAAAA"), // C at pos 13
        read_record("rb8", 10, b"AAAACAAAAAAAAAAAAAAA"), // C at pos 14
    ];
    let psp_a = make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads_a);
    let psp_b = make_psp_for_sample(dir.path(), &fasta, "NA00002", &reads_b);

    // Estimate-contamination → TOML.
    let toml_out = dir.path().join("contam.toml");
    let est_args = estimate_contamination_args(
        fasta.clone(),
        toml_out.clone(),
        vec![psp_a.clone(), psp_b.clone()],
    );
    run_estimate_contamination(&est_args).expect("run_estimate_contamination OK");
    assert!(toml_out.exists(), "TOML must exist at the final path");
    let toml_body = fs::read_to_string(&toml_out).expect("read TOML");
    assert!(toml_body.contains("[provenance]"));
    assert!(toml_body.contains("\"NA00001\""));
    assert!(toml_body.contains("\"NA00002\""));

    // Var-calling with the artefact wired in.
    let vcf_out = dir.path().join("cohort.vcf");
    let args = var_calling_args(fasta, vcf_out.clone(), vec![psp_a, psp_b], Some(toml_out));
    run_var_calling(&args).expect("run_var_calling OK");

    assert!(vcf_out.exists(), "VCF must exist");
    let body = fs::read_to_string(&vcf_out).expect("read VCF");
    assert!(body.contains("##fileformat=VCFv4"));
    assert!(body.contains("\tNA00001\tNA00002"));
}

/// **M5 follow-up — FASTA → `.psp` MD5 mismatch.** Build a `.psp`
/// whose header `chromosome.md5` is `WRONG_MD5` (forced via the
/// CRAM `@SQ M5`), then run `run_var_calling` against the real
/// FASTA. Verifies that
/// [`verify_fasta_matches_psp_chromosomes`](pop_var_caller::pop_var_caller)
/// catches the "right basename, wrong genome build" failure mode
/// the v1 basename-only check couldn't see, and surfaces a typed
/// [`VarCallingCliError::FastaContigMd5Mismatch`] before the
/// pipeline boots.
#[test]
fn var_calling_rejects_fasta_whose_bytes_dont_match_psp_md5() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let reads = vec![read_record("r1", 10, b"AAAAAAAAAAAAAAAAAAAA")];
    // Build a CRAM with a wrong @SQ M5; run_pileup propagates that
    // value into the .psp header verbatim. The FASTA itself
    // hasn't changed (still 200 A's, real MD5 = `fixture_md5()`),
    // so the post-Wave-5 cross-check catches the discrepancy.
    let cram = build_cram(dir.path(), &fasta, "NA00001", Some(WRONG_MD5), &reads);
    let psp = dir.path().join("NA00001.psp");
    let pa = pileup_args(fasta.clone(), psp.clone(), vec![cram]);
    run_pileup(&pa).expect("run_pileup OK (Stage 1 doesn't cross-check FASTA bytes)");

    let vcf_out = dir.path().join("out.vcf");
    let args = var_calling_args(fasta, vcf_out, vec![psp], None);
    let err = run_var_calling(&args).expect_err("M5 cross-check must fire");
    match err {
        VarCallingCliError::FastaContigMd5Mismatch {
            contig,
            fasta_md5,
            psp_md5,
        } => {
            assert_eq!(contig, "chr1");
            assert_eq!(psp_md5, WRONG_MD5);
            assert_eq!(fasta_md5, fixture_md5());
        }
        other => panic!("expected FastaContigMd5Mismatch, got {other:?}"),
    }
}

/// **Mi23 reference-mismatch coverage** (estimate-contamination
/// side). The `.psp` header's `reference` field is set from the
/// FASTA basename Stage 1 was driven with; when
/// `run_estimate_contamination` is later invoked with a
/// `--reference` whose basename differs, the orchestrator must
/// surface a typed [`EstimateContaminationCliError::ReferenceMismatch`]
/// *before* the side-pass starts (a Stage-3 / Stage-4 surface check
/// — analogous to `run_var_calling`'s).
#[test]
fn estimate_contamination_rejects_reference_mismatch() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path()); // basename "ref.fa"
    let reads = vec![read_record("r1", 10, b"AAAAAAAAAAAAAAAAAAAA")];
    let psp = make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads);

    // Build a second FASTA at a different basename. Identical bytes
    // + identical .fai — only the path differs, which is enough to
    // trip the basename comparison the orchestrator runs.
    let other_fasta = dir.path().join("other.fa");
    let other_fai = dir.path().join("other.fa.fai");
    fs::copy(&fasta, &other_fasta).expect("copy fasta");
    fs::copy(format!("{}.fai", fasta.display()), &other_fai).expect("copy fai");

    let toml_out = dir.path().join("contam.toml");
    let args = estimate_contamination_args(other_fasta, toml_out, vec![psp]);
    let err = run_estimate_contamination(&args).expect_err("must error on reference mismatch");
    match err {
        EstimateContaminationCliError::ReferenceMismatch {
            psp_ref,
            supplied_ref,
            ..
        } => {
            assert_eq!(psp_ref, "ref.fa");
            assert_eq!(supplied_ref, "other.fa");
        }
        other => panic!("expected ReferenceMismatch, got {other:?}"),
    }
}

/// **Mi23 reference-mismatch coverage** (var-calling side). Same
/// invariant as `estimate_contamination_rejects_reference_mismatch`
/// but for `run_var_calling`: the basename comparison happens
/// before the cohort pipeline boots, and the typed
/// [`VarCallingCliError::ReferenceMismatch`] surfaces the .psp +
/// both basenames so the user can fix their CLI invocation.
#[test]
fn var_calling_reports_reference_mismatch_for_psp_with_different_header_basename() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path()); // basename "ref.fa"
    let reads = vec![read_record("r1", 10, b"AAAAAAAAAAAAAAAAAAAA")];
    let psp = make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads);

    let other_fasta = dir.path().join("other.fa");
    let other_fai = dir.path().join("other.fa.fai");
    fs::copy(&fasta, &other_fasta).expect("copy fasta");
    fs::copy(format!("{}.fai", fasta.display()), &other_fai).expect("copy fai");

    let vcf_out = dir.path().join("out.vcf");
    let args = var_calling_args(other_fasta, vcf_out, vec![psp], None);
    let err = run_var_calling(&args).expect_err("must error on reference mismatch");
    match err {
        VarCallingCliError::ReferenceMismatch {
            psp_ref,
            supplied_ref,
            ..
        } => {
            assert_eq!(psp_ref, "ref.fa");
            assert_eq!(supplied_ref, "other.fa");
        }
        other => panic!("expected ReferenceMismatch, got {other:?}"),
    }
}

/// **M1/M2 follow-up — end-to-end walker-error CRAM test.** Synthesize
/// a CRAM with two reads overlapping at the same position, then call
/// `run_var_calling_from_bam` with `max_active_reads = 1` so the
/// second read trips the walker's defensive bound. Verifies the
/// orchestrator's two contracts from the 2026-05-19 review:
///
///   1. Walker errors surface as a typed
///      [`VarCallingFromBamCliError::Walker`] variant (not silently
///      swallowed) — M1/M2 from the review locked the error path in
///      via the
///      [`ErrorSheddingAdapter`](pop_var_caller::pop_var_caller).
///   2. The output VCF does not exist on disk afterwards
///      (publish-then-retract: the orchestrator best-effort removes
///      `<output>` after the rename runs).
#[test]
fn var_calling_from_bam_surfaces_walker_error_on_max_active_reads_trip() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    // Two reads at the same position — both active simultaneously
    // → tripping `max_active_reads = 1` is unavoidable.
    let reads = vec![
        read_record("r1", 10, b"AAAAA"),
        read_record("r2", 10, b"AAAAA"),
    ];
    let cram = build_cram(dir.path(), &fasta, "NA12878", Some(fixture_md5()), &reads);
    let vcf_out = dir.path().join("out.vcf");

    let args = VarCallingFromBamArgs {
        reference: fasta.clone(),
        output: vcf_out.clone(),
        threads: None,
        crams: vec![cram],
        no_complexity_filter: true,
        stage1: Stage1Args {
            min_mapq: DEFAULT_MIN_MAPQ,
            no_baq: true,
            min_read_length: 0,
            keep_qc_fail: false,
            keep_duplicates: false,
            max_read_mismatch_fraction: DEFAULT_MAX_READ_MISMATCH_FRACTION,
            mismatch_bq_floor: DEFAULT_MISMATCH_BQ_FLOOR,
            baq_gap_open_prob: SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
            baq_gap_extend_prob: SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
            baq_band_half_width: SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
            baq_chunk_size: DEFAULT_BAQ_CHUNK_SIZE,
            max_snp_column_depth: DEFAULT_MAX_SNP_COLUMN_DEPTH,
            max_indel_column_depth: DEFAULT_MAX_INDEL_COLUMN_DEPTH,
            max_record_span: DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
            // The defensive bound under test.
            max_active_reads: 1,
        },
        cohort: CohortPipelineArgs {
            ploidy: DEFAULT_PLOIDY,
            complexity_window: DEFAULT_DUST_WINDOW,
            complexity_threshold: DEFAULT_DUST_THRESHOLD,
            var_group_max_span: DEFAULT_MAX_VARIANT_GROUP_SPAN,
            max_alleles_per_var: DEFAULT_MAX_ALLELES_PER_RECORD,
            inbreeding_coefficient: DEFAULT_INBREEDING_COEFFICIENT,
            em_convergence_threshold: DEFAULT_CONVERGENCE_THRESHOLD,
            em_max_iterations: DEFAULT_MAX_ITERATIONS,
            ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
            snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
            indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
            compound_alt_pseudocount: DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
            max_gq_phred: DEFAULT_MAX_GQ_PHRED,
            min_qual_phred: DEFAULT_MIN_QUAL_PHRED,
            min_alt_obs_per_sample: DEFAULT_MIN_ALT_OBS_PER_SAMPLE,
            emit_gp: false,
        },
    };
    let err = run_var_calling_from_bam(&args).expect_err("must trip walker bound");
    assert!(
        matches!(err, VarCallingFromBamCliError::Walker(_)),
        "expected Walker(_), got {err:?}"
    );
    assert!(
        !vcf_out.exists(),
        "publish-then-retract: <output> must not remain on disk after walker error"
    );
}

/// `configure_rayon_pool(None)` path returns `Ok(())` without
/// touching `build_global()`, so both calls are trivially clean.
///
/// `#[serial]` keeps this test from interleaving with other tests
/// in the same binary that might exercise rayon. The unit test
/// `configure_rayon_pool_none_is_always_ok` in
/// `pop_var_caller::common::tests` covers the helper's
/// `None`-path contract directly.
#[serial_test::serial]
#[test]
fn run_pileup_can_be_called_back_to_back() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let reads = vec![read_record("r1", 10, b"AAAAA")];
    let cram_a = build_cram(dir.path(), &fasta, "A", Some(fixture_md5()), &reads);
    let cram_b = build_cram(dir.path(), &fasta, "B", Some(fixture_md5()), &reads);

    let args1 = pileup_args(fasta.clone(), dir.path().join("a.psp"), vec![cram_a]);
    run_pileup(&args1).expect("first run_pileup OK");

    let args2 = pileup_args(fasta, dir.path().join("b.psp"), vec![cram_b]);
    run_pileup(&args2).expect("second run_pileup OK");

    assert!(dir.path().join("a.psp").exists());
    assert!(dir.path().join("b.psp").exists());
}
