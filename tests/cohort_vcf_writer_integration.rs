//! End-to-end tests for the cohort VCF writer.
//!
//! Drives a hand-built sequence of `PosteriorRecord`s through
//! `CohortVcfWriter` to disk and asserts shape on the resulting
//! `.vcf` / `.vcf.gz` file. Field-level round-trip through noodles
//! is already covered by the writer's in-module unit tests; this
//! suite focuses on the file-on-disk contract: text bytes, the
//! atomic-rename behaviour, the bgzf EOF marker, and the
//! out-of-order latch.

use std::io::Read as _;

use tempfile::tempdir;

use pop_var_caller::pileup_record::AlleleSupportStats;
use pop_var_caller::psp::header::ParsedChromosome;
use pop_var_caller::var_calling::per_group_merger::MergedAllele;
use pop_var_caller::var_calling::posterior_engine::{EmDiagnostics, PosteriorRecord, RecordLocus};
use pop_var_caller::vcf::{CohortMetadata, CohortVcfWriter, VcfWriteError, WriterConfig};

fn ref_allele(seq: &[u8]) -> MergedAllele {
    MergedAllele {
        seq: seq.to_vec(),
        is_compound: false,
        constituents: Vec::new(),
    }
}

/// Mi13: alias for `ref_allele` so call sites read REF/ALT intent
/// without anyone having to stop and check whether the fixture is
/// (re)building a REF position.
fn alt_allele(seq: &[u8]) -> MergedAllele {
    ref_allele(seq)
}

fn support(num_obs: u32) -> AlleleSupportStats {
    // `AlleleSupportStats` is `#[non_exhaustive]`; integration tests
    // outside the crate go through the `new` constructor.
    AlleleSupportStats::new(num_obs, 0.0, 0, 0, 0, 0, 0)
}

fn metadata_two_samples() -> CohortMetadata {
    CohortMetadata {
        sample_names: vec!["S0".into(), "S1".into()],
        contigs: vec![ParsedChromosome {
            name: "chr1".into(),
            length: 1_000_000,
            md5: "0123456789abcdef0123456789abcdef".into(),
        }],
        tool_string: "pop_var_caller integration-test".into(),
        command_line: "pop_var_caller cohort --output out.vcf".into(),
        paralog_provenance: String::new(),
    }
}

/// Three-record fixture: biallelic SNP, triallelic SNP, insertion.
fn three_record_fixture() -> Vec<PosteriorRecord> {
    // Biallelic SNP at chr1:100, A→T.
    let biallelic = PosteriorRecord {
        locus: RecordLocus {
            chrom_id: 0,
            start: 100,
            end: 100,
        },
        alleles: vec![ref_allele(b"A"), alt_allele(b"T")],
        ploidy: 2,
        n_samples: 2,
        n_genotypes: 3,
        allele_frequencies: vec![0.75, 0.25],
        compound_frequencies: vec![None, None],
        posteriors: vec![0.98, 0.01, 0.01, 0.05, 0.90, 0.05],
        best_genotype: vec![0, 1],
        gq_phred: vec![60.0, 40.0],
        qual_phred: 150.0,
        scalars: vec![support(20), support(0), support(10), support(10)],
        other_scalars: vec![],
        chain_anchor_flags: vec![false; 4],
        diagnostics: EmDiagnostics {
            iterations: 5,
            final_max_delta_p: 1e-6,
            converged: true,
        },
        paralog_posterior: None,
    };

    // Triallelic SNP at chr1:500, A → T, A → C.
    // genotype enumeration at ploidy=2, n_alleles=3:
    //   [[0,0],[0,1],[1,1],[0,2],[1,2],[2,2]]
    let triallelic = PosteriorRecord {
        locus: RecordLocus {
            chrom_id: 0,
            start: 500,
            end: 500,
        },
        alleles: vec![ref_allele(b"A"), alt_allele(b"T"), alt_allele(b"C")],
        ploidy: 2,
        n_samples: 2,
        n_genotypes: 6,
        allele_frequencies: vec![0.6, 0.25, 0.15],
        compound_frequencies: vec![None, None, None],
        posteriors: vec![
            0.7, 0.1, 0.05, 0.1, 0.025, 0.025, // S0: argmax = 0 (0/0)
            0.05, 0.05, 0.05, 0.7, 0.1, 0.05, // S1: argmax = 3 (0/2)
        ],
        best_genotype: vec![0, 3],
        gq_phred: vec![55.0, 45.0],
        qual_phred: 220.0,
        scalars: vec![
            support(20),
            support(0),
            support(0), // S0
            support(10),
            support(0),
            support(10), // S1
        ],
        other_scalars: vec![],
        chain_anchor_flags: vec![false; 6],
        diagnostics: EmDiagnostics {
            iterations: 4,
            final_max_delta_p: 1e-6,
            converged: true,
        },
        paralog_posterior: None,
    };

    // Insertion at chr1:900, REF=A, ALT=ATT (length-3 inserted bases
    // after the anchor base, expressed as the full ALT string per
    // VCF convention).
    let insertion = PosteriorRecord {
        locus: RecordLocus {
            chrom_id: 0,
            start: 900,
            end: 900,
        },
        alleles: vec![ref_allele(b"A"), alt_allele(b"ATT")],
        ploidy: 2,
        n_samples: 2,
        n_genotypes: 3,
        allele_frequencies: vec![0.5, 0.5],
        compound_frequencies: vec![None, None],
        posteriors: vec![0.05, 0.90, 0.05, 0.01, 0.01, 0.98],
        best_genotype: vec![1, 2],
        gq_phred: vec![60.0, 70.0],
        qual_phred: f64::INFINITY, // exercises the QUAL cap path
        scalars: vec![support(15), support(15), support(0), support(30)],
        other_scalars: vec![],
        chain_anchor_flags: vec![false; 4],
        diagnostics: EmDiagnostics {
            iterations: 3,
            final_max_delta_p: 1e-6,
            converged: true,
        },
        paralog_posterior: None,
    };

    vec![biallelic, triallelic, insertion]
}

#[test]
fn plain_text_path_writes_three_records_with_default_config() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out.vcf");
    let metadata = metadata_two_samples();
    let config = WriterConfig::new(out.clone());

    let mut writer = CohortVcfWriter::new(metadata, config).unwrap();
    for record in three_record_fixture() {
        writer.write_record(&record).unwrap();
    }
    writer.finish().unwrap();

    let text = std::fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let data_lines: Vec<&str> = lines
        .iter()
        .copied()
        .filter(|l| !l.starts_with('#'))
        .collect();
    assert_eq!(data_lines.len(), 3, "expected 3 data lines\n{text}");

    let header_blob = text
        .lines()
        .take_while(|l| l.starts_with('#'))
        .collect::<Vec<_>>()
        .join("\n");
    assert!(header_blob.contains("##fileformat=VCFv4.4"));
    assert!(header_blob.contains("##source=pop_var_caller integration-test"));
    assert!(header_blob.contains("##contig=<ID=chr1"));
    assert!(header_blob.contains("##INFO=<ID=AF"));
    assert!(header_blob.contains("##INFO=<ID=CA"));
    assert!(header_blob.contains("##FORMAT=<ID=AD,Number=R,Type=Integer"));
    assert!(
        !header_blob.contains("##FORMAT=<ID=GP"),
        "GP should not be declared with emit_gp=false"
    );

    // Records are in genomic order.
    for (i, &expected_pos) in [100, 500, 900].iter().enumerate() {
        let fields: Vec<&str> = data_lines[i].split('\t').collect();
        assert_eq!(fields[0], "chr1");
        assert_eq!(fields[1], expected_pos.to_string());
    }

    // The triallelic line has ALT = "T,C".
    assert!(data_lines[1].contains("\tT,C\t"));

    // The insertion line uses the 9999 QUAL cap (its qual_phred was INF).
    let ins_fields: Vec<&str> = data_lines[2].split('\t').collect();
    let qual: f32 = ins_fields[5].parse().unwrap();
    assert_eq!(qual, 9999.0);

    // Default-config FORMAT key string is the four-key default.
    for line in &data_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[8], "GT:GQ:DP:AD");
    }

    // tmp file is gone.
    let tmp = dir.path().join("out.vcf.tmp");
    assert!(!tmp.exists(), "tmp left behind: {tmp:?}");
}

#[test]
fn plain_text_path_with_emit_gp_adds_format_column() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out.vcf");
    let metadata = metadata_two_samples();
    let config = WriterConfig::new(out.clone()).with_emit_gp(true);

    let mut writer = CohortVcfWriter::new(metadata, config).unwrap();
    for record in three_record_fixture() {
        writer.write_record(&record).unwrap();
    }
    writer.finish().unwrap();

    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("##FORMAT=<ID=GP,Number=G,Type=Float"));

    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[8], "GT:GQ:DP:AD:GP");
        // Per-sample cells now carry 5 colon-separated subfields.
        let s0_parts: Vec<&str> = fields[9].split(':').collect();
        let s1_parts: Vec<&str> = fields[10].split(':').collect();
        assert_eq!(s0_parts.len(), 5);
        assert_eq!(s1_parts.len(), 5);
        // GP lists are non-empty.
        assert!(!s0_parts[4].is_empty());
        assert!(!s1_parts[4].is_empty());
    }
}

#[test]
fn bgzf_path_writes_compressed_vcf_with_eof_marker() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out.vcf.gz");
    let metadata = metadata_two_samples();
    let config = WriterConfig::new(out.clone());

    let mut writer = CohortVcfWriter::new(metadata, config).unwrap();
    for record in three_record_fixture() {
        writer.write_record(&record).unwrap();
    }
    writer.finish().unwrap();

    let raw = std::fs::read(&out).unwrap();
    // htslib EOF marker, 28 bytes.
    let eof: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
    assert!(raw.len() >= eof.len());
    assert_eq!(
        &raw[raw.len() - eof.len()..],
        &eof,
        "bgzf file must end with the htslib EOF block"
    );

    // Decompress via noodles_bgzf and assert content shape.
    let file = std::fs::File::open(&out).unwrap();
    let mut bgzf = noodles_bgzf::io::Reader::new(file);
    let mut decompressed = Vec::new();
    bgzf.read_to_end(&mut decompressed).unwrap();
    let text = String::from_utf8(decompressed).unwrap();
    let data_lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(data_lines.len(), 3);
}

/// Helper for the MAPQ-INFO test: builds an `AlleleSupportStats`
/// that pins all five of `num_obs`, `mapq_sum`, `mapq_sum_sq` so
/// the test can predict the cohort-pooled mean / Welch's t.
fn support_with_mapq(num_obs: u32, mapq_sum: u32, mapq_sum_sq: u64) -> AlleleSupportStats {
    AlleleSupportStats::new(num_obs, 0.0, 0, 0, 0, mapq_sum, mapq_sum_sq)
}

/// Phase B: a record where REF reads are clean (all MAPQ=60) and
/// ALT reads are multi-mappers (MAPQ=20 mean, with spread). Verify
/// the four new INFO fields land with the expected values.
///
/// REF: 20 reads × MAPQ=60   → sum=1200, sum_sq=72000  → mean 60.00
/// ALT:  4 reads with MAPQs (0, 20, 30, 30) → sum=80, sum_sq=2200 → mean 20.00
/// Variance (sample, ddof=1):
///   var_ref = 0  (all reads MAPQ=60)
///   var_alt = (2200 - 80²/4) / (4-1) = (2200 - 1600)/3 = 200
/// Welch's t = (20 - 60) / sqrt(200/4 + 0/20)
///           = -40 / sqrt(50)
///           ≈ -5.657
#[test]
fn mapq_info_fields_reflect_cohort_pooled_stats() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out.vcf");
    let metadata = metadata_two_samples();
    let config = WriterConfig::new(out.clone());

    let mut writer = CohortVcfWriter::new(metadata, config).unwrap();
    let record = PosteriorRecord {
        locus: RecordLocus {
            chrom_id: 0,
            start: 100,
            end: 100,
        },
        alleles: vec![ref_allele(b"A"), alt_allele(b"T")],
        ploidy: 2,
        n_samples: 2,
        n_genotypes: 3,
        allele_frequencies: vec![0.8, 0.2],
        compound_frequencies: vec![None, None],
        posteriors: vec![0.98, 0.01, 0.01, 0.05, 0.90, 0.05],
        best_genotype: vec![0, 1],
        gq_phred: vec![60.0, 40.0],
        qual_phred: 150.0,
        // Layout is sample-major: scalars[s][a]. Pool across samples:
        //   REF: S0 = 10 reads × MAPQ 60, S1 = 10 reads × MAPQ 60
        //   ALT: S0 = 0 reads, S1 = 4 reads (mixed)
        // S0 row: REF (n=10, sum=600, ssq=36000), ALT (n=0, 0, 0)
        // S1 row: REF (n=10, sum=600, ssq=36000), ALT (n=4, 80, 2200)
        scalars: vec![
            support_with_mapq(10, 600, 36_000),
            support_with_mapq(0, 0, 0),
            support_with_mapq(10, 600, 36_000),
            support_with_mapq(4, 80, 2_200),
        ],
        other_scalars: vec![],
        chain_anchor_flags: vec![false; 4],
        diagnostics: EmDiagnostics {
            iterations: 5,
            final_max_delta_p: 1e-6,
            converged: true,
        },
        paralog_posterior: None,
    };
    writer.write_record(&record).unwrap();
    writer.finish().unwrap();

    let text = std::fs::read_to_string(&out).unwrap();
    // Header declarations.
    assert!(text.contains("##INFO=<ID=MQRef,Number=1,Type=Float"));
    assert!(text.contains("##INFO=<ID=MQAlt,Number=A,Type=Float"));
    assert!(text.contains("##INFO=<ID=MQDiff,Number=A,Type=Float"));
    assert!(text.contains("##INFO=<ID=MQDiffT,Number=A,Type=Float"));

    let data_lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(data_lines.len(), 1);
    let info = data_lines[0].split('\t').nth(7).unwrap();

    // Parse INFO into a map of key → value string.
    let info_map: std::collections::HashMap<&str, &str> = info
        .split(';')
        .filter_map(|kv| kv.split_once('=').or(Some((kv, ""))))
        .collect();

    // MQRef: cohort mean MAPQ of REF reads = 60.0
    let mq_ref: f32 = info_map["MQRef"].parse().unwrap();
    assert!((mq_ref - 60.0).abs() < 0.01, "MQRef={mq_ref}");

    // MQAlt: cohort mean MAPQ of ALT reads = 80 / 4 = 20.0
    let mq_alt: f32 = info_map["MQAlt"].parse().unwrap();
    assert!((mq_alt - 20.0).abs() < 0.01, "MQAlt={mq_alt}");

    // MQDiff: 20.0 - 60.0 = -40.0
    let mq_diff: f32 = info_map["MQDiff"].parse().unwrap();
    assert!((mq_diff - (-40.0)).abs() < 0.01, "MQDiff={mq_diff}");

    // MQDiffT: -40 / sqrt(200/4 + 0/20) = -40 / sqrt(50) ≈ -5.657
    let mq_t: f32 = info_map["MQDiffT"].parse().unwrap();
    let expected_t = -40.0_f32 / (50.0_f32).sqrt();
    assert!(
        (mq_t - expected_t).abs() < 0.05,
        "MQDiffT={mq_t} expected ≈ {expected_t}"
    );
}

/// Phase C: spot-check the allele-balance drop helper directly (the
/// end-to-end driver path is exercised in `cohort_cli_integration`).
/// Verifies the contract branches of `record_fails_allele_balance`:
///
///  - a deep het at low ALT fraction (≈0.2 VAF) gets the "drop" verdict;
///  - the same skew at shallow depth does NOT (the test is depth-aware
///    and stays agnostic where the data can't support a verdict);
///  - a permissive threshold keeps even the deep-skewed record;
///  - a clean het carrier (or a hom-alt carrier) anywhere rescues the site;
///  - a multiallelic record is skipped (v1 is biallelic-only).
#[test]
fn allele_balance_filter_decision_matches_contract() {
    use pop_var_caller::var_calling::allele_balance::{
        DEFAULT_AB_CONCENTRATION as S, DEFAULT_AB_MIN_LOG_LR as THR,
    };
    use pop_var_caller::var_calling::vcf_writer::record_fails_allele_balance_for_test as fails;

    // AB filter reads only num_obs (allele balance ignores MAPQ).
    let obs = |n: u32| support_with_mapq(n, 0, 0);

    // Single deep het at VAF≈0.2: REF=200, ALT=50.
    let suspect = PosteriorRecord {
        locus: RecordLocus {
            chrom_id: 0,
            start: 100,
            end: 100,
        },
        alleles: vec![ref_allele(b"A"), alt_allele(b"T")],
        ploidy: 2,
        n_samples: 1,
        n_genotypes: 3,
        allele_frequencies: vec![0.8, 0.2],
        compound_frequencies: vec![None, None],
        posteriors: vec![0.05, 0.90, 0.05],
        best_genotype: vec![1], // het
        gq_phred: vec![40.0],
        qual_phred: 150.0,
        scalars: vec![obs(200), obs(50)],
        other_scalars: vec![],
        chain_anchor_flags: vec![false; 2],
        diagnostics: EmDiagnostics {
            iterations: 5,
            final_max_delta_p: 1e-6,
            converged: true,
        },
        paralog_posterior: None,
    };
    assert!(
        fails(&suspect, THR, S),
        "deep 0.2-VAF het should be dropped"
    );
    // A very permissive threshold ( ≈ filter-off ) keeps it.
    assert!(
        !fails(&suspect, -1000.0, S),
        "permissive threshold must keep"
    );

    // Same skew, shallow depth: REF=4, ALT=1. Not enough evidence → keep.
    let shallow = PosteriorRecord {
        scalars: vec![obs(4), obs(1)],
        ..suspect.clone()
    };
    assert!(
        !fails(&shallow, THR, S),
        "shallow 0.2-VAF must NOT be dropped"
    );

    // Two-sample: failing het + a clean balanced het → site rescued.
    let rescued_by_het = PosteriorRecord {
        n_samples: 2,
        best_genotype: vec![1, 1],
        gq_phred: vec![40.0, 60.0],
        scalars: vec![obs(200), obs(50), obs(125), obs(125)],
        chain_anchor_flags: vec![false; 4],
        posteriors: vec![0.05, 0.90, 0.05, 0.02, 0.96, 0.02],
        ..suspect.clone()
    };
    assert!(
        !fails(&rescued_by_het, THR, S),
        "a clean het carrier keeps the site"
    );

    // Two-sample: failing het + a hom-alt carrier → untested support keeps it.
    let rescued_by_hom_alt = PosteriorRecord {
        best_genotype: vec![1, 2], // sample 1 is hom-alt
        scalars: vec![obs(200), obs(50), obs(0), obs(250)],
        ..rescued_by_het.clone()
    };
    assert!(
        !fails(&rescued_by_hom_alt, THR, S),
        "a hom-alt carrier keeps the site"
    );

    // Biallelic INDEL with the same deep skew: SNP-only guard → never dropped
    // (indel allele balance from AD is unreliable).
    let indel = PosteriorRecord {
        alleles: vec![ref_allele(b"AT"), alt_allele(b"A")], // a deletion
        n_samples: 1,
        best_genotype: vec![1],
        gq_phred: vec![40.0],
        scalars: vec![obs(200), obs(50)],
        chain_anchor_flags: vec![false; 2],
        posteriors: vec![0.05, 0.90, 0.05],
        ..suspect.clone()
    };
    assert!(
        !fails(&indel, THR, S),
        "indel records are skipped (SNP-only)"
    );

    // Multiallelic: v1 is biallelic-only → never dropped here.
    let multiallelic = PosteriorRecord {
        alleles: vec![ref_allele(b"A"), alt_allele(b"T"), alt_allele(b"C")],
        n_samples: 1,
        best_genotype: vec![1],
        gq_phred: vec![40.0],
        scalars: vec![obs(200), obs(50), obs(0)],
        chain_anchor_flags: vec![false; 3],
        posteriors: vec![0.0; 6],
        ..suspect.clone()
    };
    assert!(
        !fails(&multiallelic, THR, S),
        "multiallelic records are skipped"
    );
}

#[test]
fn out_of_order_records_surface_a_clear_error() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out.vcf");
    let metadata = metadata_two_samples();
    let config = WriterConfig::new(out);
    let mut writer = CohortVcfWriter::new(metadata, config).unwrap();
    let records = three_record_fixture();
    writer.write_record(&records[2]).unwrap(); // chr1:900
    let err = writer.write_record(&records[0]).unwrap_err(); // chr1:100
    match err {
        VcfWriteError::RecordOutOfOrder {
            chrom_id,
            pos,
            prev_chrom_id,
            prev_pos,
        } => {
            assert_eq!((chrom_id, pos), (0, 100));
            assert_eq!((prev_chrom_id, prev_pos), (0, 900));
        }
        other => panic!("expected RecordOutOfOrder, got {other:?}"),
    }
}
