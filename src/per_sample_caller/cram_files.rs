//! Test fixtures for Group B tests of `cram_input` — synthetic CRAM
//! and FASTA files written to a tempdir, exercising the noodles
//! writer and `CramMergedReader::new`'s file-open / header-parsing
//! path.
//!
//! See `ia/feature_implementation_plans/per_sample_caller_cram_input.md`
//! §"Test-fixture helpers".

use std::fs::File;
use std::io::{self, Write as IoWrite};
use std::num::NonZero;
use std::path::PathBuf;

use bstr::BString;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use sam::alignment::record_buf::RecordBuf;
use sam::header::record::value::Map;
use sam::header::record::value::map::header::Version;
use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};
use tempfile::TempDir;

/// One reference contig used to assemble a synthetic FASTA + `.fai`.
#[derive(Debug, Clone)]
pub(crate) struct ContigSpec {
    pub name: String,
    pub length: u64,
}

/// Knobs that test-driven CRAM writes need to override on the
/// generated SAM header.
#[derive(Debug, Clone, Default)]
pub(crate) struct HeaderOverrides {
    /// Override `@HD SO:`. `None` means use "coordinate".
    pub sort_order: Option<String>,
    /// `(read_group_id, optional SM)` tuples. If empty, defaults to a
    /// single `@RG ID:rg0 SM:sample0`. SM = `None` produces an `@RG`
    /// with no SM tag (used by the missing-SM rejection test).
    pub read_groups: Vec<(String, Option<String>)>,
    /// Per-contig MD5 override map (contig name → 32 hex chars or
    /// raw bytes). The `.fai` itself does not carry MD5; this lets
    /// tests inject `M5` into the CRAM header's `@SQ` entry only.
    pub md5_overrides: Vec<(String, String)>,
    /// Per-contig name override (CRAM `@SQ` SN may disagree with the
    /// FASTA name — used by FASTA mismatch tests).
    pub name_overrides: Vec<(String, String)>,
    /// Per-contig length override (same idea, for length mismatch
    /// tests).
    pub length_overrides: Vec<(String, u64)>,
}

/// Build a FASTA + `.fai` in a freshly allocated tempdir. Returns the
/// FASTA path; the `TempDir` is kept alive by the caller. Each contig
/// gets `length` 'A' bases on a single line.
pub(crate) fn build_fasta(contigs: &[ContigSpec]) -> io::Result<(TempDir, PathBuf)> {
    let dir = tempfile::tempdir()?;
    let fasta_path = dir.path().join("ref.fa");
    let fai_path = dir.path().join("ref.fa.fai");

    let mut fa = File::create(&fasta_path)?;
    let mut fai = File::create(&fai_path)?;
    let mut offset: u64 = 0;
    for contig in contigs {
        // Header line: ">{name}\n"
        let header = format!(">{}\n", contig.name);
        fa.write_all(header.as_bytes())?;
        offset += header.len() as u64;
        // Sequence line: 'A' * length + '\n'
        let seq_len = usize::try_from(contig.length).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "contig length does not fit in usize on this target",
            )
        })?;
        let seq = vec![b'A'; seq_len];
        fa.write_all(&seq)?;
        fa.write_all(b"\n")?;
        // .fai entry: name\tlength\toffset\tline_bases\tline_width
        let line_bases = contig.length;
        let line_width = line_bases + 1;
        writeln!(
            fai,
            "{}\t{}\t{}\t{}\t{}",
            contig.name, contig.length, offset, line_bases, line_width
        )?;
        offset += line_bases + 1;
    }
    Ok((dir, fasta_path))
}

/// Build a coordinate-sorted CRAM (or whatever sort order is set in
/// `overrides.sort_order`) in a freshly allocated tempdir. Returns
/// the CRAM path. The header's `@SQ` list is taken from `contigs`
/// with optional per-contig name/length/MD5 overrides applied.
pub(crate) fn build_cram(
    fasta_path: &std::path::Path,
    contigs: &[ContigSpec],
    overrides: &HeaderOverrides,
    records: &[RecordBuf],
) -> io::Result<(TempDir, PathBuf)> {
    let dir = tempfile::tempdir()?;
    let cram_path = dir.path().join("test.cram");

    let header = build_sam_header(contigs, overrides);
    let repository = repository_for(fasta_path)?;

    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(&cram_path)?;
    writer.write_header(&header)?;
    for record in records {
        sam::alignment::io::Write::write_alignment_record(&mut writer, &header, record)?;
    }
    sam::alignment::io::Write::finish(&mut writer, &header)?;
    Ok((dir, cram_path))
}

fn repository_for(fasta_path: &std::path::Path) -> io::Result<fasta::Repository> {
    let indexed_fasta_reader =
        fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
    let adapter = fasta::repository::adapters::IndexedReader::new(indexed_fasta_reader);
    Ok(fasta::Repository::new(adapter))
}

fn build_sam_header(contigs: &[ContigSpec], overrides: &HeaderOverrides) -> sam::Header {
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    // @HD with sort order
    let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
    let so = overrides
        .sort_order
        .clone()
        .unwrap_or_else(|| "coordinate".into());
    hd.other_fields_mut()
        .insert(SORT_ORDER, BString::from(so.as_bytes()));

    // @SQ list
    let mut builder = sam::Header::builder().set_header(hd);
    for contig in contigs {
        let name = name_for(contig, overrides);
        let length = length_for(contig, overrides);
        let length_usize =
            usize::try_from(length).expect("contig length must fit in usize on this target");
        let length_nz = NonZero::new(length_usize).expect("contig length must be non-zero");
        let mut ref_seq_map = Map::<ReferenceSequence>::new(length_nz);
        if let Some((_, md5)) = overrides
            .md5_overrides
            .iter()
            .find(|(n, _)| n == &contig.name)
        {
            ref_seq_map
                .other_fields_mut()
                .insert(MD5_CHECKSUM, BString::from(md5.as_bytes()));
        }
        builder = builder.add_reference_sequence(name.into_bytes(), ref_seq_map);
    }

    // @RG list — default to one read group with SM:sample0 if the
    // override list is empty.
    let read_groups: Vec<(String, Option<String>)> = if overrides.read_groups.is_empty() {
        vec![("rg0".into(), Some("sample0".into()))]
    } else {
        overrides.read_groups.clone()
    };
    for (rg_id, sm) in &read_groups {
        let mut rg_map = Map::<ReadGroup>::default();
        if let Some(sm_value) = sm {
            rg_map
                .other_fields_mut()
                .insert(SAMPLE, BString::from(sm_value.as_bytes()));
        }
        builder = builder.add_read_group(rg_id.as_bytes(), rg_map);
    }

    builder.build()
}

fn name_for(contig: &ContigSpec, overrides: &HeaderOverrides) -> String {
    overrides
        .name_overrides
        .iter()
        .find(|(n, _)| n == &contig.name)
        .map(|(_, replacement)| replacement.clone())
        .unwrap_or_else(|| contig.name.clone())
}

fn length_for(contig: &ContigSpec, overrides: &HeaderOverrides) -> u64 {
    overrides
        .length_overrides
        .iter()
        .find(|(n, _)| n == &contig.name)
        .map(|(_, replacement)| *replacement)
        .unwrap_or(contig.length)
}

/// Write a 26-byte CRAM file definition with a non-3 major version.
/// Enough to exercise `read_file_definition()` returning a 4.x tag —
/// the rejection happens before any container is decoded, so we do
/// not need a full file body.
///
/// CRAM file definition layout (per the CRAM 3.x spec, identical for
/// 4.x): 4 magic bytes ("CRAM"), 1 major byte, 1 minor byte,
/// 20-byte file ID. Total 26 bytes.
pub(crate) fn build_cram_with_major_version(
    major: u8,
    minor: u8,
) -> io::Result<(TempDir, PathBuf)> {
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("forced_version.cram");
    let mut file = File::create(&path)?;
    file.write_all(b"CRAM")?;
    file.write_all(&[major, minor])?;
    file.write_all(&[0u8; 20])?;
    Ok((dir, path))
}
