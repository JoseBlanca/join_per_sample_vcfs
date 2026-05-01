//! Test fixtures for Group A tests of `cram_input` — synthetic
//! `RecordBuf`s assembled from minimal field specs, plus an
//! `OpenCram` factory that wraps a `Vec` in the boxed iterator shape
//! the merge expects. No filesystem, no FASTA, no noodles writer.
//!
//! See `ia/feature_implementation_plans/per_sample_caller_cram_input.md`
//! §"Test-fixture helpers".

use std::io;
use std::path::PathBuf;

use bstr::BString;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::MappingQuality;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record_buf::Cigar;
use noodles_sam::alignment::record_buf::QualityScores;
use noodles_sam::alignment::record_buf::RecordBuf;
use noodles_sam::alignment::record_buf::Sequence;

use super::cram_input::{CigarOp, ContigEntry, ContigList, OpenCram};

/// Spec for a synthetic record built by `record_spec`. Every test
/// composes its records by mutating fields on a base — keeping the
/// builder declarative.
#[derive(Debug, Clone)]
pub(crate) struct RecordSpec {
    pub qname: String,
    pub flag: u16,
    pub ref_id: usize,
    pub pos: u64,
    pub mapq: u8,
    pub cigar_ops: Vec<CigarOp>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub mate_ref_id: Option<usize>,
    pub mate_pos: Option<u64>,
}

pub(crate) fn record_spec(spec: RecordSpec) -> RecordBuf {
    let mut rb = RecordBuf::default();
    *rb.name_mut() = Some(BString::from(spec.qname.as_bytes()));
    *rb.flags_mut() = Flags::from(spec.flag);
    *rb.reference_sequence_id_mut() = Some(spec.ref_id);
    if spec.pos > 0 {
        *rb.alignment_start_mut() = usize::try_from(spec.pos)
            .ok()
            .and_then(noodles_core::Position::new);
    }
    if spec.mapq > 0 {
        *rb.mapping_quality_mut() = MappingQuality::new(spec.mapq);
    }
    let cigar_ops: Vec<Op> = spec
        .cigar_ops
        .iter()
        .map(|op| match *op {
            CigarOp::Match(n) => Op::new(Kind::Match, n as usize),
            CigarOp::Insertion(n) => Op::new(Kind::Insertion, n as usize),
            CigarOp::Deletion(n) => Op::new(Kind::Deletion, n as usize),
            CigarOp::Skip(n) => Op::new(Kind::Skip, n as usize),
            CigarOp::SoftClip(n) => Op::new(Kind::SoftClip, n as usize),
            CigarOp::HardClip(n) => Op::new(Kind::HardClip, n as usize),
            CigarOp::Padding(n) => Op::new(Kind::Pad, n as usize),
            CigarOp::SeqMatch(n) => Op::new(Kind::SequenceMatch, n as usize),
            CigarOp::SeqMismatch(n) => Op::new(Kind::SequenceMismatch, n as usize),
        })
        .collect();
    *rb.cigar_mut() = Cigar::from(cigar_ops);
    *rb.sequence_mut() = Sequence::from(spec.seq);
    *rb.quality_scores_mut() = QualityScores::from(spec.qual);
    *rb.mate_reference_sequence_id_mut() = spec.mate_ref_id;
    *rb.mate_alignment_start_mut() = spec
        .mate_pos
        .and_then(|p| usize::try_from(p).ok())
        .and_then(noodles_core::Position::new);
    rb
}

pub(crate) fn open_cram_from_records(path_for_errors: &str, records: Vec<RecordBuf>) -> OpenCram {
    let iter = records.into_iter().map(Ok::<_, io::Error>);
    OpenCram {
        path_for_errors: PathBuf::from(path_for_errors),
        records: Box::new(iter),
    }
}

pub(crate) fn default_contigs() -> ContigList {
    ContigList {
        entries: vec![ContigEntry {
            name: "chr1".into(),
            length: 100_000,
            md5: None,
        }],
    }
}
