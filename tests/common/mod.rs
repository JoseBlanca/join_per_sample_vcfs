//! Shared fixture helpers for the cohort-CLI integration tests.
//!
//! Both `pileup_cli_integration.rs` and `cohort_cli_integration.rs`
//! used to carry their own near-identical `build_fasta` /
//! `build_sam_header` / `build_cram` / `read_record` helpers. Mi20
//! from the 2026-05-19 cohort CLI review consolidates them here.
//!
//! Each integration-test target lives in a separate binary, but
//! `mod common;` from inside `tests/` is the canonical Cargo-test
//! idiom for sharing helper code across them — the file is compiled
//! once per test target.
//!
//! The fixture is intentionally tiny: a single 200-base all-A
//! reference contig named `chr1`. Stage 1 doesn't cross-check the
//! `.psp`'s `chromosome.md5` against the FASTA bytes (Stage 3 is the
//! one that does, with a real FASTA), so the CRAM-side `@SQ M5` can
//! be a fixed placeholder.

#![allow(dead_code)] // not every integration-test target uses every helper

use std::fs::File;
use std::io::Write;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use bstr::BString;
use md5::{Digest, Md5};
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_csi as csi;
use noodles_fasta as fasta;
use noodles_sam as sam;
use sam::alignment::record::Flags;
use sam::alignment::record::MappingQuality;
use sam::alignment::record::cigar::Op;
use sam::alignment::record::cigar::op::Kind;
use sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};
use sam::header::record::value::Map;
use sam::header::record::value::map::header::Version;
use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};

/// Default contig name used by the helpers.
pub const CONTIG_NAME: &str = "chr1";

/// Default contig length used by the helpers.
pub const CONTIG_LEN: usize = 200;

/// 32-hex MD5 of the all-`A` contig the fixture FASTA produces.
/// Computed lazily on first access so the constant stays
/// in lock-step with [`CONTIG_LEN`] — Wave 5's M5 check
/// (`verify_fasta_matches_psp_chromosomes`) verifies the FASTA
/// bytes against the `.psp` header's md5, so the value here
/// **must** match the real MD5 of `[b'A'; CONTIG_LEN]`.
pub fn fixture_md5() -> &'static str {
    static MD5: OnceLock<String> = OnceLock::new();
    MD5.get_or_init(|| {
        let mut h = Md5::new();
        h.update(vec![b'A'; CONTIG_LEN]);
        let digest = h.finalize();
        let mut out = String::with_capacity(32);
        for b in digest {
            use std::fmt::Write as _;
            let _ = write!(&mut out, "{b:02x}");
        }
        out
    })
}

/// An obviously-wrong MD5 (32 zeroes) that integration tests can
/// inject when they specifically want to trip the FASTA-vs-`.psp`
/// MD5 check (M5 follow-up). For the happy-path fixtures, use
/// [`fixture_md5`] instead.
pub const WRONG_MD5: &str = "00000000000000000000000000000000";

/// Write a minimal FASTA + `.fai` index for the single-contig
/// fixture (`CONTIG_NAME`, `CONTIG_LEN` bases, all `A`). Returns
/// the FASTA path.
pub fn build_fasta(dir: &Path) -> PathBuf {
    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let mut fa = File::create(&fasta_path).expect("create fasta");
    let mut fai = File::create(&fai_path).expect("create fai");
    writeln!(fa, ">{}", CONTIG_NAME).unwrap();
    let header_len = (CONTIG_NAME.len() + 2) as u64; // '>' + name + '\n'
    let seq = vec![b'A'; CONTIG_LEN];
    fa.write_all(&seq).unwrap();
    fa.write_all(b"\n").unwrap();
    writeln!(
        fai,
        "{}\t{}\t{}\t{}\t{}",
        CONTIG_NAME,
        CONTIG_LEN,
        header_len,
        CONTIG_LEN,
        CONTIG_LEN + 1
    )
    .unwrap();
    fasta_path
}

/// Build a single-contig SAM header tagged with `sample` and an
/// `@SQ M5` only when `md5` is `Some`. The `None` path is used by
/// the pileup integration test's "header missing M5" rejection
/// case; cohort tests pass `Some(fixture_md5())` to get the
/// matching MD5 that satisfies Wave 5's M5 follow-up FASTA
/// cross-check, and the M5-mismatch test passes `Some(WRONG_MD5)`
/// to deliberately trip that check.
pub fn build_sam_header(sample: &str, md5: Option<&str>) -> sam::Header {
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
    hd.other_fields_mut()
        .insert(SORT_ORDER, BString::from(b"coordinate".as_ref()));

    let length_nz = NonZero::new(CONTIG_LEN).unwrap();
    let mut ref_seq = Map::<ReferenceSequence>::new(length_nz);
    if let Some(md5_str) = md5 {
        ref_seq
            .other_fields_mut()
            .insert(MD5_CHECKSUM, BString::from(md5_str.as_bytes()));
    }

    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut()
        .insert(SAMPLE, BString::from(sample.as_bytes()));

    sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(CONTIG_NAME.as_bytes().to_vec(), ref_seq)
        .add_read_group(b"rg0".as_ref(), rg)
        .build()
}

/// Build a single-contig CRAM for `sample` with the given records.
/// CRAM filename is `{sample}.cram` under `dir`. `md5` controls the
/// `@SQ M5` tag — see [`build_sam_header`] for the three intended
/// values. Re-reads `fasta_path` internally to drive the CRAM
/// encoder's reference repository.
pub fn build_cram(
    dir: &Path,
    fasta_path: &Path,
    sample: &str,
    md5: Option<&str>,
    records: &[RecordBuf],
) -> PathBuf {
    let cram_path = dir.join(format!("{sample}.cram"));
    let header = build_sam_header(sample, md5);

    let indexed_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .expect("indexed fasta reader");
    let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
    let repository = fasta::Repository::new(adapter);

    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(&cram_path)
        .expect("cram writer");
    writer.write_header(&header).expect("write header");
    for record in records {
        sam::alignment::io::Write::write_alignment_record(&mut writer, &header, record)
            .expect("write record");
    }
    sam::alignment::io::Write::finish(&mut writer, &header).expect("finish cram");
    cram_path
}

/// Build a single-contig BAM for `sample` with the given records.
/// BAM filename is `{sample}.bam` under `dir`. `md5` controls the
/// `@SQ M5` tag — see [`build_sam_header`] for the three intended
/// values. Unlike CRAM, BAM stores sequences inline and does not
/// take a reference repository.
pub fn build_bam(dir: &Path, sample: &str, md5: Option<&str>, records: &[RecordBuf]) -> PathBuf {
    let bam_path = dir.join(format!("{sample}.bam"));
    let header = build_sam_header(sample, md5);

    let file = File::create(&bam_path).expect("create bam");
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).expect("write header");
    for record in records {
        sam::alignment::io::Write::write_alignment_record(&mut writer, &header, record)
            .expect("write record");
    }
    writer.try_finish().expect("finish bam");
    bam_path
}

/// Build a `.csi` next to `bam_path` by scanning the BAM. Returns
/// the index path. Mirrors the build-csi-for-bam helper in
/// `crate::bam::index_preflight` but lives here so the tests don't
/// reach across the integration / library boundary for what is
/// just a one-shot file-writing fixture.
pub fn build_csi(bam_path: &Path) -> PathBuf {
    use csi::binning_index::Indexer;
    use csi::binning_index::index::reference_sequence::bin::Chunk;
    use csi::binning_index::index::reference_sequence::index::BinnedIndex;
    use sam::alignment::Record as _;

    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .expect("open bam for csi build");
    let header = reader.read_header().expect("read header for csi build");

    let mut indexer: Indexer<BinnedIndex> = Indexer::default();
    let mut chunk_start = reader.get_ref().virtual_position();
    let mut record = bam::Record::default();

    while reader.read_record(&mut record).expect("read record") != 0 {
        let chunk_end = reader.get_ref().virtual_position();
        let alignment_context = match (
            record.reference_sequence_id().transpose().expect("ref id"),
            record.alignment_start().transpose().expect("start"),
            record.alignment_end().transpose().expect("end"),
        ) {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };
        let chunk = Chunk::new(chunk_start, chunk_end);
        indexer
            .add_record(alignment_context, chunk)
            .expect("add record");
        chunk_start = chunk_end;
    }

    let index = indexer.build(header.reference_sequences().len());
    let csi_path = csi_path_for(bam_path);
    csi::fs::write(&csi_path, &index).expect("write csi");
    csi_path
}

/// Build a `.bai` next to `bam_path` by scanning the BAM. Returns
/// the index path. Used by the .bai-fallback integration test.
pub fn build_bai(bam_path: &Path) -> PathBuf {
    let index = bam::fs::index(bam_path).expect("build bai");
    let bai_path = bai_path_for(bam_path);
    bam::bai::fs::write(&bai_path, &index).expect("write bai");
    bai_path
}

fn csi_path_for(bam: &Path) -> PathBuf {
    let mut s = bam.as_os_str().to_owned();
    s.push(".csi");
    PathBuf::from(s)
}

fn bai_path_for(bam: &Path) -> PathBuf {
    let mut s = bam.as_os_str().to_owned();
    s.push(".bai");
    PathBuf::from(s)
}

/// Build a single-contig FASTA + `.fai` with a parameterised
/// contig length (all-`A` bases). Returns `(fasta_path, md5_hex)`.
/// Used by tests that need enough records to exercise PSP-writer
/// block boundaries; the default [`build_fasta`] is too small
/// (200 bp).
pub fn build_fasta_with_len(dir: &Path, contig_len: usize) -> (PathBuf, String) {
    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let mut fa = File::create(&fasta_path).expect("create fasta");
    let mut fai = File::create(&fai_path).expect("create fai");
    writeln!(fa, ">{}", CONTIG_NAME).unwrap();
    let header_len = (CONTIG_NAME.len() + 2) as u64;
    let seq = vec![b'A'; contig_len];
    fa.write_all(&seq).unwrap();
    fa.write_all(b"\n").unwrap();
    writeln!(
        fai,
        "{}\t{}\t{}\t{}\t{}",
        CONTIG_NAME,
        contig_len,
        header_len,
        contig_len,
        contig_len + 1
    )
    .unwrap();
    let mut h = Md5::new();
    h.update(&seq);
    let digest = h.finalize();
    let mut md5_hex = String::with_capacity(32);
    for b in digest {
        use std::fmt::Write as _;
        let _ = write!(&mut md5_hex, "{b:02x}");
    }
    (fasta_path, md5_hex)
}

/// Like [`build_sam_header`] but with a parameterised contig length.
pub fn build_sam_header_with_len(
    sample: &str,
    md5: Option<&str>,
    contig_len: usize,
) -> sam::Header {
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
    hd.other_fields_mut()
        .insert(SORT_ORDER, BString::from(b"coordinate".as_ref()));

    let length_nz = NonZero::new(contig_len).unwrap();
    let mut ref_seq = Map::<ReferenceSequence>::new(length_nz);
    if let Some(md5_str) = md5 {
        ref_seq
            .other_fields_mut()
            .insert(MD5_CHECKSUM, BString::from(md5_str.as_bytes()));
    }

    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut()
        .insert(SAMPLE, BString::from(sample.as_bytes()));

    sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(CONTIG_NAME.as_bytes().to_vec(), ref_seq)
        .add_read_group(b"rg0".as_ref(), rg)
        .build()
}

/// Like [`build_cram`] but with a parameterised contig length so it
/// can pair with [`build_fasta_with_len`].
pub fn build_cram_with_len(
    dir: &Path,
    fasta_path: &Path,
    sample: &str,
    md5: Option<&str>,
    records: &[RecordBuf],
    contig_len: usize,
) -> PathBuf {
    let cram_path = dir.join(format!("{sample}.cram"));
    let header = build_sam_header_with_len(sample, md5, contig_len);

    let indexed_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .expect("indexed fasta reader");
    let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
    let repository = fasta::Repository::new(adapter);

    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(&cram_path)
        .expect("cram writer");
    writer.write_header(&header).expect("write header");
    for record in records {
        sam::alignment::io::Write::write_alignment_record(&mut writer, &header, record)
            .expect("write record");
    }
    sam::alignment::io::Write::finish(&mut writer, &header).expect("finish cram");
    cram_path
}

/// Build a synthetic single-end read at `pos` (1-based), all-Match
/// CIGAR, MAPQ 60, BQ 30. `qname` is the read name; `seq` is the
/// read bases.
pub fn read_record(qname: &str, pos: usize, seq: &[u8]) -> RecordBuf {
    let mut rb = RecordBuf::default();
    *rb.name_mut() = Some(BString::from(qname.as_bytes()));
    *rb.flags_mut() = Flags::from(0u16);
    *rb.reference_sequence_id_mut() = Some(0);
    *rb.alignment_start_mut() = noodles_core::Position::new(pos);
    *rb.mapping_quality_mut() = MappingQuality::new(60);
    let cigar_ops = vec![Op::new(Kind::Match, seq.len())];
    *rb.cigar_mut() = Cigar::from(cigar_ops);
    *rb.sequence_mut() = Sequence::from(seq.to_vec());
    *rb.quality_scores_mut() = QualityScores::from(vec![30u8; seq.len()]);
    rb
}
