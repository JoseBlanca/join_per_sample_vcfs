//! Byte-identity oracle harness for the re-architected cohort pipeline.
//!
//! The re-architecture (`var_calling_new::`) is built in parallel with the
//! shipping `var_calling::` and graded against a **hard byte-identity
//! constraint**: for the same `.psp` cohort, the new pipeline's VCF must match
//! the old one's GT/GQ/AD/AF/AC/FILTER exactly (QUAL excluded — QUAL is
//! non-deterministic since `03e2221`). The retained `var_calling::` is the
//! branch base (`main`), so **new-vs-old is new-vs-`main`** — a far stronger
//! gate than a saved golden file (plan §7).
//!
//! This harness drives **both** pipelines from the same [`VarCallingArgs`] on
//! the same in-process `.psp` fixtures and diffs the VCF bodies (header lines
//! and the QUAL column stripped).
//!
//! ## Status: red until Phase 4
//!
//! The new entry point ([`run_var_calling_new`]) is a Phase-0 stub returning
//! `NotImplemented`, so the diff test is **`#[ignore]`d** — it keeps the tree
//! green (principle 4) while documenting the gate the later phases turn green.
//! Run it on demand with `cargo test --test var_calling_new_oracle -- --ignored`.
//! **Phase 4 removes the `#[ignore]`**, at which point the oracle must pass.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::{build_cram, build_fasta, fixture_md5, read_record};
use noodles_sam::alignment::record_buf::RecordBuf;
use pop_var_caller::bam::alignment_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MISMATCH_BQ_FLOOR,
};
use pop_var_caller::baq::{
    SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH, SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
    SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
};
use pop_var_caller::pileup::per_sample::baq_stream::DEFAULT_BAQ_CHUNK_SIZE;
use pop_var_caller::pileup::walker::{
    DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_ACTIVE_READS, DEFAULT_MAX_INDEL_COLUMN_DEPTH,
    DEFAULT_MAX_RECORD_SPAN, DEFAULT_MAX_SNP_COLUMN_DEPTH,
};
use pop_var_caller::pop_var_caller::cli::shared_args::{CohortPipelineArgs, Stage1Args};
use pop_var_caller::pop_var_caller::var_calling::{VarCallingArgs, run_var_calling};
use pop_var_caller::pop_var_caller::{PileupArgs, run_pileup};
use pop_var_caller::var_calling::dust_filter::{DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW};
use pop_var_caller::var_calling::per_group_merger::{
    DEFAULT_MAX_ALLELES_LH_CALC, DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY,
};
use pop_var_caller::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT,
};
use pop_var_caller::var_calling::variant_grouping::DEFAULT_MAX_VARIANT_GROUP_SPAN;
use pop_var_caller::var_calling::{
    DEFAULT_MIN_ALT_OBS_PER_SAMPLE, DEFAULT_MIN_MAPQ_DIFF_T, DEFAULT_MIN_QUAL_PHRED,
};
use pop_var_caller::var_calling_new::pipeline::run_var_calling as run_var_calling_new;
use tempfile::TempDir;

// ---------------------------------------------------------------------
// Fixture/arg builders (mirrors `cohort_cli_integration.rs`; kept local so the
// oracle is self-driving — these can be factored into `tests/common` later).
// ---------------------------------------------------------------------

fn pileup_args(reference: PathBuf, output: PathBuf, alignment_files: Vec<PathBuf>) -> PileupArgs {
    PileupArgs {
        reference,
        output,
        regions: None,
        build_map_file_index: true,
        threads: None,
        alignment_files,
        block_target_bytes: pop_var_caller::psp::writer::TARGET_BLOCK_BYTES,
        block_window_bp: pop_var_caller::psp::writer::DEFAULT_BLOCK_WINDOW_BP,
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
    }
}

fn var_calling_args(
    reference: PathBuf,
    output: PathBuf,
    psp_files: Vec<PathBuf>,
) -> VarCallingArgs {
    VarCallingArgs {
        reference,
        output,
        regions: None,
        threads: None,
        contamination_estimates: None,
        no_complexity_filter: true, // tiny ref; sdust would mask everything
        target_variants_per_chunk: 0,
        low_memory: false,
        psp_files,
        cohort: CohortPipelineArgs {
            ploidy: DEFAULT_PLOIDY,
            complexity_window: DEFAULT_DUST_WINDOW,
            complexity_threshold: DEFAULT_DUST_THRESHOLD,
            var_group_max_span: DEFAULT_MAX_VARIANT_GROUP_SPAN,
            max_alleles_per_var: DEFAULT_MAX_ALLELES_PER_RECORD,
            max_alleles_lh_calc: DEFAULT_MAX_ALLELES_LH_CALC,
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
            no_mapq_diff_filter: false,
            min_mapq_diff_t: DEFAULT_MIN_MAPQ_DIFF_T,
            emit_gp: false,
        },
    }
}

/// Build a `.psp` for one sample by running the real pileup orchestrator on a
/// synthetic CRAM. `max_read_mismatch_fraction = 0.0` disables the per-read
/// mismatch filter so SNP-bearing synthetic reads survive into the pileup.
fn make_psp_for_sample(dir: &Path, fasta: &Path, sample: &str, reads: &[RecordBuf]) -> PathBuf {
    let cram = build_cram(dir, fasta, sample, Some(fixture_md5()), reads);
    let psp = dir.join(format!("{sample}.psp"));
    let mut args = pileup_args(fasta.to_path_buf(), psp.clone(), vec![cram]);
    args.stage1.max_read_mismatch_fraction = 0.0;
    run_pileup(&args).expect("run_pileup OK");
    psp
}

// ---------------------------------------------------------------------
// VCF normalisation + diff (QUAL excluded; header lines stripped).
// ---------------------------------------------------------------------

/// The QUAL column index in a VCF data line (0-based: CHROM POS ID REF ALT
/// QUAL …). QUAL is excluded from the byte-identity contract (plan §6).
const QUAL_COL: usize = 5;

/// Normalise a VCF for byte-identity comparison: drop `##` meta-header lines,
/// keep the `#CHROM` column header and all data lines, and blank the QUAL
/// column of each data line. Returns the canonical body as a single string.
fn normalise_vcf(path: &Path) -> String {
    let body = fs::read_to_string(path).expect("read VCF");
    let mut out: Vec<String> = Vec::new();
    for line in body.lines() {
        if line.starts_with("##") || line.trim().is_empty() {
            continue;
        }
        if line.starts_with('#') {
            out.push(line.to_string()); // column header (#CHROM ... samples)
            continue;
        }
        let mut cols: Vec<&str> = line.split('\t').collect();
        if cols.len() > QUAL_COL {
            cols[QUAL_COL] = "."; // exclude QUAL from the contract
        }
        out.push(cols.join("\t"));
    }
    out.join("\n")
}

/// Drive both pipelines on the same cohort and assert the normalised VCFs are
/// byte-identical. This is the gate Phase 4 must satisfy.
fn assert_oracle_byte_identical(dir: &Path, fasta: PathBuf, psps: Vec<PathBuf>) {
    let old_vcf = dir.join("old.vcf");
    let new_vcf = dir.join("new.vcf");

    run_var_calling(&var_calling_args(
        fasta.clone(),
        old_vcf.clone(),
        psps.clone(),
    ))
    .expect("OLD var_calling produced a VCF");
    run_var_calling_new(&var_calling_args(fasta, new_vcf.clone(), psps))
        .expect("NEW var_calling_new pipeline produced a VCF");

    let old_norm = normalise_vcf(&old_vcf);
    let new_norm = normalise_vcf(&new_vcf);
    assert_eq!(
        old_norm, new_norm,
        "re-architecture byte-identity gate failed: new VCF differs from old (QUAL excluded)"
    );
}

/// Three-sample cohort with a single SNP — the smallest end-to-end oracle.
///
/// `#[ignore]` until Phase 4 wires the new pipeline (the stub entry returns
/// `NotImplemented`, so this is red by design). Remove the attribute in P4.
#[test]
#[ignore = "re-architecture oracle: green from Phase 4 (new pipeline entry is a P0 stub)"]
fn oracle_three_sample_cohort_byte_identical() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());

    // Sample B carries a SNP (C) at pos 11 on *two* reads, so it clears
    // `min_alt_obs_per_sample = 2` and a real variant record emits — making
    // the eventual P4 byte-identity gate non-trivial (not two empty VCFs).
    let reads_a = vec![
        read_record("r1", 10, b"AAAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];
    let reads_b = vec![
        read_record("r1", 10, b"ACAAA"),
        read_record("r3", 10, b"ACAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];
    let reads_c = vec![
        read_record("r1", 10, b"AAAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];

    let psps = vec![
        make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads_a),
        make_psp_for_sample(dir.path(), &fasta, "NA00002", &reads_b),
        make_psp_for_sample(dir.path(), &fasta, "NA00003", &reads_c),
    ];

    assert_oracle_byte_identical(dir.path(), fasta, psps);
}

/// The oracle scaffolding itself (fixture build + VCF normaliser) is exercised
/// unconditionally so it can't silently rot while the byte-identity test is
/// ignored: build a cohort, run the OLD pipeline, and confirm the normaliser
/// produces a stable, QUAL-blanked body.
#[test]
fn oracle_normaliser_is_stable_on_old_pipeline_output() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let reads = vec![
        read_record("r1", 10, b"ACAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];
    let psp = make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads);

    let vcf = dir.path().join("old.vcf");
    run_var_calling(&var_calling_args(fasta, vcf.clone(), vec![psp])).expect("old run OK");

    let norm = normalise_vcf(&vcf);
    assert!(
        norm.starts_with("#CHROM"),
        "normaliser keeps the column header"
    );
    assert!(
        !norm.lines().any(|l| l.starts_with("##")),
        "normaliser drops ## meta-header lines"
    );
    // Idempotent: re-normalising the same file yields the same body.
    assert_eq!(norm, normalise_vcf(&vcf), "normaliser is deterministic");
}
