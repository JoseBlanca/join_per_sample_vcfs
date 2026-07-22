//! Fixtures shared by this module's tests: a small multi-contig reference, a
//! header builder, and an indexed BAM on disk.
//!
//! They live in their own file because both the gate's tests and the region
//! query's need them, and a `#[cfg(test)] mod tests` block is private to its
//! own module — the alternative is two copies that drift.

use std::fs::File;
use std::num::NonZero;
use std::path::PathBuf;

use noodles_bam as bam;
use noodles_core::Position as RecordPosition;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::{Flags, MappingQuality};
use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
use sam::header::record::value::Map;
use sam::header::record::value::map::{ReadGroup, ReferenceSequence};
use tempfile::TempDir;

use crate::bam::index_preflight::preflight_alignment_indexes;
use crate::ng::reference_info::{ReferenceInfo, ReferenceSource, read_reference_info};
use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};

/// The fixture reference: two contigs of **different lengths**, so a
/// permutation of an `@SQ` list is detectable on `name` and a re-labelling on
/// `length`.
pub(crate) const FIXTURE_CONTIGS: [(&str, usize); 2] = [("chr1", 100), ("chr2", 200)];

/// A header builder taking `(sort order, contigs, read groups)`, so each test
/// states only the field it is about. `None` omits the tag entirely.
pub(crate) fn header(
    sort_order_value: Option<&str>,
    contigs: &[(&str, usize, Option<&str>)],
    read_groups: &[(&str, Option<&str>)],
) -> sam::Header {
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    let mut hd = Map::<sam::header::record::value::map::Header>::default();
    if let Some(value) = sort_order_value {
        hd.other_fields_mut()
            .insert(SORT_ORDER, value.as_bytes().into());
    }

    let mut builder = sam::Header::builder().set_header(hd);

    for (name, length, md5) in contigs {
        let mut sq = Map::<ReferenceSequence>::new(NonZero::new(*length).unwrap());
        if let Some(md5) = md5 {
            sq.other_fields_mut()
                .insert(MD5_CHECKSUM, md5.as_bytes().into());
        }
        builder = builder.add_reference_sequence(*name, sq);
    }

    for (id, sample) in read_groups {
        let mut rg = Map::<ReadGroup>::default();
        if let Some(sample) = sample {
            rg.other_fields_mut()
                .insert(SAMPLE, sample.as_bytes().into());
        }
        builder = builder.add_read_group(*id, rg);
    }

    builder.build()
}

/// A header whose `@SQ` list is `contigs`, `SO:coordinate`, and one read group
/// naming `NA12878` — the shape a file has to have to get past the gate.
pub(crate) fn bam_header(contigs: &[(&str, usize, Option<&str>)]) -> sam::Header {
    header(Some("coordinate"), contigs, &[("rg1", Some("NA12878"))])
}

/// [`FIXTURE_CONTIGS`] in the `@SQ` shape, with no `M5` tags.
pub(crate) fn matching_contigs() -> Vec<(&'static str, usize, Option<&'static str>)> {
    FIXTURE_CONTIGS
        .iter()
        .map(|(name, length)| (*name, *length, None))
        .collect()
}

/// A `ReferenceInfo` over [`FIXTURE_CONTIGS`].
///
/// `with_digests` picks the arm: the `Fasta` arm reads the genome and carries
/// real per-contig MD5s, while the `Fai` arm cannot, so its digests are `None`
/// and the MD5 half of reconciliation is a no-op. Tests that care about the
/// digest comparison need both.
pub(crate) fn fixture_reference(with_digests: bool) -> (TempDir, ReferenceInfo) {
    let specs: Vec<ContigSpec> = FIXTURE_CONTIGS
        .iter()
        .map(|(name, length)| ContigSpec {
            name: (*name).to_string(),
            length: *length as u64,
        })
        .collect();
    let (dir, fasta) = build_fasta(&specs).expect("build fasta");

    let source = if with_digests {
        ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        }
    } else {
        ReferenceSource::Fai(crate::ng::reference_info::sibling_fai_path(&fasta))
    };
    (dir, read_reference_info(source).expect("read reference"))
}

/// A 10 bp perfectly-matching read at `start` on `reference_sequence_id`.
///
/// `qname` matters for the region-query oracle, which identifies reads by name
/// to compare two independently-produced streams.
pub(crate) fn read_named(qname: &str, reference_sequence_id: usize, start: usize) -> RecordBuf {
    read_named_with_length(qname, reference_sequence_id, start, 10)
}

pub(crate) fn read_named_with_length(
    qname: &str,
    reference_sequence_id: usize,
    start: usize,
    length: usize,
) -> RecordBuf {
    RecordBuf::builder()
        .set_name(qname.as_bytes())
        .set_reference_sequence_id(reference_sequence_id)
        // Explicit, because `RecordBuf`'s default flags are `UNMAPPED` — a
        // fixture that left them alone would be silently dropped by filter #1
        // in any test that runs the real filter.
        .set_flags(Flags::empty())
        .set_mapping_quality(MappingQuality::new(60).expect("mapq in range"))
        .set_alignment_start(RecordPosition::try_from(start).unwrap())
        .set_cigar([Op::new(Kind::Match, length)].into_iter().collect())
        .set_sequence(Sequence::from(vec![b'A'; length]))
        .set_quality_scores(QualityScores::from(vec![30u8; length]))
        .build()
}

/// Write a BAM holding `records` and build its index beside it.
///
/// Returns the `TempDir` as well as the path: bind it, or the directory is
/// removed the moment the call returns and the handle points at nothing.
pub(crate) fn indexed_bam(header: &sam::Header, records: &[RecordBuf]) -> (TempDir, PathBuf) {
    let (dir, path) = unindexed_bam(header, records);
    preflight_alignment_indexes(std::slice::from_ref(&path), true).expect("build index");
    (dir, path)
}

/// The same BAM with **no** index beside it — for the gate's index check.
pub(crate) fn unindexed_bam(header: &sam::Header, records: &[RecordBuf]) -> (TempDir, PathBuf) {
    let dir = TempDir::new().expect("tempdir");
    let path = dir.path().join("sample.bam");

    let mut writer = bam::io::Writer::new(File::create(&path).expect("create bam"));
    writer.write_header(header).expect("write header");
    for record in records {
        writer
            .write_alignment_record(header, record)
            .expect("write record");
    }
    writer.try_finish().expect("finish");

    (dir, path)
}

/// The same records written as a **CRAM**, with a `.crai` beside it, over the
/// fixture reference. Returns the CRAM's dir and path plus the FASTA's dir and
/// path — a CRAM cannot be decoded without the reference, so its `ReferenceInfo`
/// has to come from the `Fasta` arm.
///
/// **`records` must all be on one contig.** Not a choice: `noodles_cram::fs::index`
/// decodes *multi-reference* slices with `fasta::Repository::default()` — an
/// empty repository, marked `// TODO` in noodles 0.93 (`src/fs/index.rs:137`) —
/// so indexing a CRAM whose reads are stored as differences from the reference
/// panics with "invalid reference sequence name". Reads spanning two contigs
/// land in one slice, which is exactly that case. Porting noodles' indexer with
/// a real repository is not possible from outside the crate:
/// `ReferenceSequenceContext` and `Slice::header` are private.
///
/// Production's CRAM fixtures are single-contig too, which is why this has not
/// bitten before. The `.crai` contig walk is covered instead by a hand-built
/// index in `region_query`'s tests, which needs no file at all.
pub(crate) fn indexed_cram(records: &[RecordBuf]) -> (TempDir, PathBuf, TempDir, PathBuf) {
    debug_assert!(
        records
            .iter()
            .all(|record| record.reference_sequence_id() == Some(0)),
        "the CRAM fixture must stay on one contig — see this function's docs"
    );
    use crate::pileup::per_sample::cram_files::{HeaderOverrides, build_cram};

    let specs: Vec<ContigSpec> = FIXTURE_CONTIGS
        .iter()
        .map(|(name, length)| ContigSpec {
            name: (*name).to_string(),
            length: *length as u64,
        })
        .collect();
    let (fasta_dir, fasta) = build_fasta(&specs).expect("build fasta");
    let (cram_dir, cram_path) = build_cram(
        &fasta,
        &specs,
        &HeaderOverrides {
            read_groups: vec![("rg1".to_string(), Some("NA12878".to_string()))],
            ..HeaderOverrides::default()
        },
        records,
    )
    .expect("build cram");

    let index = noodles_cram::fs::index(&cram_path).expect(
        "noodles can index a single-reference CRAM; see the doc above for why \
         the fixture must stay single-contig",
    );
    let crai_path = PathBuf::from(format!("{}.crai", cram_path.display()));
    noodles_cram::crai::fs::write(&crai_path, &index).expect("write crai");

    (cram_dir, cram_path, fasta_dir, fasta)
}

/// A single-contig CRAM over a **long** contig, with enough reads to fill
/// several containers.
///
/// noodles writes 10240 records per container, so a fixture that stays under
/// that produces one container and one `.crai` entry — which exercises the
/// container decode but *none* of the `.crai` walk: not the multi-entry loop,
/// not the container-level early stop, not the span skip. Returns the CRAM, its
/// FASTA, and the contig length.
pub(crate) fn multi_container_cram(
    contig_length: usize,
    read_count: usize,
) -> (TempDir, PathBuf, TempDir, PathBuf) {
    use crate::pileup::per_sample::cram_files::{HeaderOverrides, build_cram};

    let specs = vec![ContigSpec {
        name: "chr1".to_string(),
        length: contig_length as u64,
    }];
    let (fasta_dir, fasta) = build_fasta(&specs).expect("build fasta");

    // Spread the reads evenly along the contig and keep them in coordinate
    // order, so successive containers cover successive stretches — which is
    // what makes the container-level early stop observable.
    let step = (contig_length - 40) / read_count.max(1);
    let records: Vec<RecordBuf> = (0..read_count)
        .map(|i| read_named_with_length(&format!("r{i}"), 0, 1 + i * step.max(1), 30))
        .collect();

    let (cram_dir, cram_path) = build_cram(
        &fasta,
        &specs,
        &HeaderOverrides {
            read_groups: vec![("rg1".to_string(), Some("NA12878".to_string()))],
            ..HeaderOverrides::default()
        },
        &records,
    )
    .expect("build cram");

    let index = noodles_cram::fs::index(&cram_path).expect("index a single-reference CRAM");
    let crai_path = PathBuf::from(format!("{}.crai", cram_path.display()));
    noodles_cram::crai::fs::write(&crai_path, &index).expect("write crai");

    (cram_dir, cram_path, fasta_dir, fasta)
}

/// One 10 bp read at the start of the first contig — enough to make a file
/// non-empty for tests that never read it.
pub(crate) fn one_read() -> Vec<RecordBuf> {
    vec![read_named("read-1", 0, 1)]
}
