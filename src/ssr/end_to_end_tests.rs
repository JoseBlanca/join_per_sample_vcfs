//! End-to-end SSR integration test: synthetic reference + per-sample BAMs +
//! catalog → the real `ssr-pileup` driver (one `.ssr.psp` per sample) → the real
//! `ssr-call` driver → a multi-sample VCF.
//!
//! This is the first test that exercises the whole BAM→VCF chain through the two
//! production drivers together, so it guards the cross-stage contracts the
//! per-stage tests cannot see in isolation: the catalog `reference_md5` binding
//! carried from `ssr-pileup` into each `.ssr.psp` and re-checked by `ssr-call`,
//! the chromosome-table reconciliation between the `.ssr.psp` headers and the
//! catalog, and the coordinate frames lining up from the BAM through to the VCF
//! POS column.
//!
//! **Catalog generation (`ssr-catalog` / `trf-mod`) is deliberately skipped** —
//! the catalog is written directly with [`CatalogWriter`] so the test carries no
//! external-tool dependency. The genome FASTA is all-`A`; the SSR pileup realigns
//! each read against the catalog's embedded `ref_seq`, not the genome bases, so
//! the tract content lives in the reads + the catalog, never the FASTA.

#![cfg(test)]

use std::fs::File;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles_sam::alignment::RecordBuf;
use tempfile::TempDir;

use crate::bam::segment_reader::SegmentReadFilter;
use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};
use crate::ssr::catalog::CatalogParams;
use crate::ssr::catalog::io::{CatalogHeader, CatalogWriter};
use crate::ssr::cohort::driver::{self as call_driver, SsrCallConfig};
use crate::ssr::pileup::driver::{self as pileup_driver, SsrPileupConfig};
use crate::ssr::types::{Locus, Motif};

/// Reference contig the whole cohort sits on.
const CONTIG: &str = "chr1";
/// Contig length (all-`A`; only the coordinate frame matters — see the module doc).
const CONTIG_LEN: usize = 256;
/// The reference allele length (in `CA` repeat units) at every catalog locus.
const REF_UNITS: u32 = 8;

/// One spanning read whose repeat tract is `units` copies of `CA`, embedded in
/// `GGGGGG` / `TTTTTT` flanks with short `A` pads, aligned all-Match at `start - 10`
/// (1-based). The SSR pileup delimits the `CA` tract within the read and tallies the
/// observed sequence, so the (naive) all-Match CIGAR over a tract longer/shorter than
/// the reference is exactly the realignment case Stage 1 exists to handle.
fn spanning_read(start: u32, qname: &str, units: u32) -> RecordBuf {
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record::{Flags, MappingQuality};
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

    let mut seq = Vec::new();
    seq.extend_from_slice(b"AAAAA");
    seq.extend_from_slice(b"GGGGGG");
    for _ in 0..units {
        seq.extend_from_slice(b"CA");
    }
    seq.extend_from_slice(b"TTTTTT");
    seq.extend_from_slice(b"AAAAA");
    let len = seq.len();

    RecordBuf::builder()
        .set_name(qname.as_bytes())
        .set_reference_sequence_id(0)
        .set_flags(Flags::default())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_alignment_start(Position::try_from((start - 10) as usize).unwrap())
        .set_cigar([Op::new(Kind::Match, len)].into_iter().collect())
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(vec![40u8; len]))
        .build()
}

/// Reads for one sample's genotype at a locus: a clean peak at each allele length
/// (a homozygote gets full depth, a het splits it) plus light ±1 stutter, mirroring
/// the cohort EM's `reads_for` fixture so the genotype resolves confidently.
fn genotype_reads(start: u32, alleles: &[u32], qprefix: &str, reads: &mut Vec<RecordBuf>) {
    let peak = if alleles.len() == 1 { 40 } else { 20 };
    let mut n = 0u32;
    let push = |units: u32, reads: &mut Vec<RecordBuf>, n: &mut u32| {
        reads.push(spanning_read(start, &format!("{qprefix}_{n}"), units));
        *n += 1;
    };
    for &a in alleles {
        for _ in 0..peak {
            push(a, reads, &mut n);
        }
        for _ in 0..3 {
            push(a - 1, reads, &mut n);
            push(a + 1, reads, &mut n);
        }
    }
}

/// Write a coordinate-sorted single-contig BAM (`@RG SM`, a trusted `@SQ M5`) holding
/// `reads` in the order given (the caller supplies them in coordinate order).
fn write_bam(path: &Path, reads: &[RecordBuf]) {
    use bstr::BString;
    use noodles_bam as bam;
    use noodles_sam as sam;
    use sam::header::record::value::Map;
    use sam::header::record::value::map::header::Version;
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
    use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};

    let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
    hd.other_fields_mut()
        .insert(SORT_ORDER, BString::from("coordinate"));
    let mut sq = Map::<ReferenceSequence>::new(NonZero::new(CONTIG_LEN).unwrap());
    sq.other_fields_mut()
        .insert(MD5_CHECKSUM, BString::from("0".repeat(32)));
    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut()
        .insert(SAMPLE, BString::from("sample0"));
    let header = sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(CONTIG.as_bytes().to_vec(), sq)
        .add_read_group(b"rg0".to_vec(), rg)
        .build();

    let mut writer = bam::io::Writer::new(File::create(path).unwrap());
    use sam::alignment::io::Write as _;
    writer.write_header(&header).unwrap();
    for read in reads {
        writer.write_alignment_record(&header, read).unwrap();
    }
    writer.try_finish().unwrap();
}

/// Write the cohort catalog directly (no `trf-mod`): one `CA(REF_UNITS)` locus per
/// `start`, with the reference tract embedded in `GGGGGG`/`TTTTTT` flanks.
fn write_catalog(path: &Path, starts: &[u32]) {
    let header = CatalogHeader {
        tool_version: "0.1.0".into(),
        reference: "ref.fa".into(),
        reference_md5: "a".repeat(32),
        trf_mod_version: "test".into(),
        params: CatalogParams::default(),
        date: "2026-06-24".into(),
    };
    let mut w = CatalogWriter::new(File::create(path).unwrap(), &header).unwrap();
    for &start in starts {
        let mut ref_bytes = Vec::new();
        ref_bytes.extend_from_slice(b"GGGGGG");
        for _ in 0..REF_UNITS {
            ref_bytes.extend_from_slice(b"CA");
        }
        ref_bytes.extend_from_slice(b"TTTTTT");
        let locus = Locus::new(
            CONTIG.into(),
            start,
            start + 2 * REF_UNITS,
            Motif::new(b"CA").unwrap(),
            1.0,
            ref_bytes.into_boxed_slice(),
            start - 6,
        )
        .unwrap();
        w.write_locus(&locus).unwrap();
    }
    w.finish().unwrap();
}

/// Run the real `ssr-pileup` driver for one sample's BAM against the catalog,
/// returning its `.ssr.psp` path (named `{sample}.ssr.psp` so `ssr-call` derives a
/// distinct sample name from the basename).
fn pileup_sample(
    dir: &Path,
    fasta: &Path,
    catalog: &Path,
    sample: &str,
    reads: &[RecordBuf],
) -> PathBuf {
    let bam = dir.join(format!("{sample}.bam"));
    write_bam(&bam, reads);
    let psp = dir.join(format!("{sample}.ssr.psp"));
    let cfg = SsrPileupConfig {
        alignment_files: vec![bam],
        reference: fasta.to_path_buf(),
        catalog: catalog.to_path_buf(),
        output: psp.clone(),
        filter: SegmentReadFilter::default(),
        cap: 1000,
        build_index_if_missing: true,
        sample: Some(sample.to_string()),
        threads: 1,
    };
    pileup_driver::run(&cfg).expect("ssr-pileup run");
    psp
}

/// The full BAM→VCF chain: synthetic reference + per-sample BAMs + a written catalog
/// → `ssr-pileup` ×N → `ssr-call` → VCF, with one variable locus (a length
/// polymorphism) and one monomorphic locus (dropped).
#[test]
fn bam_to_vcf_emits_a_pass_variant_and_drops_a_monomorphic_locus() {
    let (_fa_dir, fasta) = build_fasta(&[ContigSpec {
        name: CONTIG.into(),
        length: CONTIG_LEN as u64,
    }])
    .unwrap();
    let dir = TempDir::new().unwrap();
    let catalog = dir.path().join("cohort.ssr.catalog");
    // Locus A (start 40): variable. Locus B (start 100): monomorphic.
    let (locus_a, locus_b) = (40u32, 100u32);
    write_catalog(&catalog, &[locus_a, locus_b]);

    // A clean length polymorphism at locus A between the reference allele (CA×8) and a
    // 2-unit-shorter allele (CA×6) — both recovered faithfully by the pileup realignment.
    // Three genotype groups of two samples each: hom-ref 8/8, het 8/6, hom-alt 6/6. Every
    // sample is 8/8 at locus B, so locus B is monomorphic across the cohort and dropped.
    let groups: [(&str, &[u32]); 3] = [("homref", &[8]), ("het", &[6, 8]), ("homalt", &[6])];
    let mut psp_files = Vec::new();
    for (label, alleles) in groups {
        for i in 0..2 {
            let mut reads = Vec::new();
            genotype_reads(locus_a, alleles, &format!("{label}A{i}"), &mut reads);
            genotype_reads(locus_b, &[8], &format!("{label}B{i}"), &mut reads);
            psp_files.push(pileup_sample(
                dir.path(),
                &fasta,
                &catalog,
                &format!("{label}{i}"),
                &reads,
            ));
        }
    }

    let output = dir.path().join("cohort.vcf");
    call_driver::run(&SsrCallConfig {
        catalog: catalog.clone(),
        psp_files,
        output: output.clone(),
        threads: 2,
        queue_depth: 0,
    })
    .expect("ssr-call run");

    let vcf = std::fs::read_to_string(&output).unwrap();
    assert!(vcf.starts_with("##fileformat=VCFv4.4"), "{vcf}");

    // 6 sample columns after the 9 fixed VCF columns.
    let chrom_line = vcf.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    assert_eq!(
        chrom_line.split('\t').count(),
        9 + 6,
        "6 sample columns\n{vcf}"
    );

    // Exactly one data record: the variable locus A; the monomorphic locus B is dropped.
    let records: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(records.len(), 1, "only the variable locus A emits\n{vcf}");
    let cols: Vec<&str> = records[0].split('\t').collect();
    assert_eq!(cols[0], CONTIG);
    assert_eq!(
        cols[1],
        (locus_a + 1).to_string(),
        "VCF POS is 1-based locus start"
    );
    assert_eq!(cols[3], "CACACACACACACACA", "REF is the CA×8 tract");
    assert_eq!(
        cols[6], "PASS",
        "the variable locus should PASS: {}",
        records[0]
    );
    assert!(
        cols[4].split(',').any(|alt| alt == "CACACACACACA"),
        "ALT should carry the CA×6 allele: {}",
        records[0]
    );

    // The six samples' called repeat counts: two hom-ref (8,8), two het (8,6), two
    // hom-alt (6,6). Collect each sample's sorted REPCN and assert the multiset.
    let mut genotypes: Vec<Vec<u16>> = cols[9..]
        .iter()
        .map(|c| {
            let repcn = c.split(':').nth(2).unwrap_or(".");
            let mut units: Vec<u16> = repcn.split(',').filter_map(|s| s.parse().ok()).collect();
            units.sort_unstable();
            units
        })
        .collect();
    genotypes.sort();
    assert_eq!(
        genotypes,
        vec![
            vec![6, 6],
            vec![6, 6],
            vec![6, 8],
            vec![6, 8],
            vec![8, 8],
            vec![8, 8],
        ],
        "expected two each of hom-alt (6,6), het (6,8), hom-ref (8,8)\n{vcf}"
    );
}
