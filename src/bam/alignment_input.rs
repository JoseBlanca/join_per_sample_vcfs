//! Alignment-file input: header validation, peek-and-scan merge,
//! per-read filter cascade. Format-agnostic at the merge/filter
//! seam; per-input record-stream decoding (currently CRAM-only,
//! inlined below) is the format-specific surface.
//!
//! The CRAM-specific design rationale lives in
//! `doc/devel/implementation_plans/per_sample_caller_cram_input.md`;
//! the BAM extension plan that motivated lifting the merger out of
//! the CRAM-named module lives in
//! `doc/devel/implementation_plans/bam_input_support.md`.

use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::bam::bam_input::{BamIndex, open_bam_record_stream, open_indexed_bam_record_stream};
use crate::bam::cram_input::{open_cram_record_stream, open_indexed_cram_record_stream};
use crate::bam::errors::AlignmentInputError;
use crate::bam::index_preflight::{AlignmentFileKind, AlignmentIndex};
use crate::fasta::{ContigEntry, ContigList};
use crate::iter_ext::BufferedPeekable;
use crate::pileup::walker::CigarOp;
use crate::pileup::walker::indel_norm::left_align_indels;

// ---------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------

/// Reads with MAPQ strictly below this are dropped. Matches bcftools'
/// default
///
/// Reads with MAPQ unavailable (SAM 0xFF / noodles `mapping_quality()`
/// returning `None`) are treated as MAPQ 0 and therefore rejected
/// under any non-zero minimum. Matches the bcftools/freebayes
/// convention.
pub const DEFAULT_MIN_MAPQ: u8 = 20;

/// Default for `AlignmentMergedReaderConfig::max_read_mismatch_fraction`.
/// 10% — comfortably outside the ~0.5-1% mismatch rate of normal
/// Illumina at MAPQ 20+, inside the regime where adapter-runthrough,
/// contamination, and chimeric tails live. See finding `F1` in
/// `ia/reviews/pileup_freebayes_comparison_2026-05-08.md`.
pub const DEFAULT_MAX_READ_MISMATCH_FRACTION: f32 = 0.10;

/// Default for `AlignmentMergedReaderConfig::mismatch_bq_floor`. Matches
/// freebayes' `BQL2`. Mismatches whose raw base quality is below this
/// floor do not count toward the mismatch fraction — keeps the filter
/// from firing on genuinely low-quality bases (which the downstream
/// likelihood already de-emphasises). `0` disables the floor entirely.
pub const DEFAULT_MISMATCH_BQ_FLOOR: u8 = 10;

/// Decoded SEQ length below this is dropped. Reads shorter than this
/// rarely contribute reliable alignments at the project's coverage
/// targets.
pub const DEFAULT_MIN_READ_LENGTH: u32 = 30;

/// `BufferedPeekable` look-ahead used per CRAM in the merge. The CRAM
/// decoder underneath already does its own slice-level batching; an
/// extra outer buffer would just delay records without saving work.
const PER_PEEKER_BUFFER_SIZE: usize = 1;

// ---------------------------------------------------------------------
// MappedRead
// ---------------------------------------------------------------------

/// A read decoded out of an alignment file (CRAM or BAM). Every
/// byte (qname, sequence, qualities, CIGAR) is owned by the
/// `MappedRead` itself rather than borrowed from the noodles
/// record it was decoded from. That lets downstream stages hold
/// onto reads without being tied to the source reader's lifetime.
///
/// Every field is public — downstream stages (BAQ, the pileup walker)
/// need to read them directly.
///
/// A borrowing design would save one allocation per read at decode
/// time, but each downstream stage would then have to copy the bytes
/// for itself before using them, which costs more in total than
/// allocating once upfront and handing the read along.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MappedRead {
    pub qname: Vec<u8>,
    pub flag: u16,
    /// Index into the merged reader's `ContigList`.
    pub ref_id: usize,
    /// 1-based leftmost mapped position.
    pub pos: u64,
    pub mapq: u8,
    pub cigar: Vec<CigarOp>,
    /// Uppercase ACGTN.
    pub seq: Vec<u8>,
    /// Phred 0-93, raw BQ before BAQ capping.
    pub qual: Vec<u8>,
    pub mate_ref_id: Option<usize>,
    pub mate_pos: Option<u64>,
    /// Reference position of the first base on this read that lies
    /// inside the mate-pair adaptor (1-based, inclusive on the read's
    /// 3′ side; for reverse-strand reads, inclusive on the 5′ side).
    /// `None` when the boundary cannot be reliably computed (single-end,
    /// mate unmapped, mates on different contigs, TLEN=0, mates on the
    /// same strand, geometry inconsistent, or the molecule is at least
    /// as long as the read so no readthrough is possible). See finding
    /// `G1` in `ia/reviews/pileup_gatk_comparison_2026-05-08.md`.
    pub adaptor_boundary: Option<u32>,
    pub source_file_index: usize,
}

// ---------------------------------------------------------------------
// Filter counts
// ---------------------------------------------------------------------

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct FilterCounts {
    pub unmapped: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub qc_fail: u64,
    pub duplicate: u64,
    pub low_mapq: u64,
    pub too_short: u64,
    /// Reads dropped because their `M`-op mismatch fraction exceeded
    /// `AlignmentMergedReaderConfig::max_read_mismatch_fraction`. See
    /// finding `F1` in
    /// `ia/reviews/pileup_freebayes_comparison_2026-05-08.md`.
    pub high_mismatch_fraction: u64,
    /// Reads dropped because their CIGAR contained an adjacent `I`/`D`
    /// pair, or started/ended with a deletion (after stripping leading
    /// soft/hard clips). See finding `G2` in
    /// `ia/reviews/pileup_gatk_comparison_2026-05-08.md`.
    pub bad_cigar: u64,
    /// Reads dropped because the BAQ stage refused to produce a usable
    /// per-base posterior (HMM overflow, ref window past chrom end,
    /// CIGAR with `N`/no-match, etc.). One bucket per
    /// [`BaqSkipReason`](super::baq::BaqSkipReason) currently does not
    /// exist — every BAQ skip reason rolls up into this counter. Split
    /// later if a deployment cares which reason dominates.
    pub baq_rejected: u64,
}

// ---------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------

/// Tunable thresholds and per-flag toggles for the per-read filter
/// cascade. `min_mapq` and `min_read_length` are `Option`s: `None` is
/// the explicit "no minimum" state. Using `0` as a disable sentinel
/// would make `Some(0)` and `None` redundant representations of the
/// same behaviour; `Option` makes the absent state structurally
/// distinct from any specific threshold.
///
/// **Why no `drop_secondary` / `drop_supplementary` toggles.** Both
/// classes of non-primary alignment are unconditionally dropped — the
/// pileup walker's per-record-closure semantics (`PileupRecord`,
/// `AlleleObservation`, phase-chain slots) are defined only for
/// primary alignments. A secondary alignment is a duplicated
/// projection of a read already represented by its primary, and a
/// supplementary alignment is a chunk of a chimeric read whose other
/// chunks are tracked separately; admitting either would silently
/// double-count alleles or break the one-slot-per-pair invariant.
/// The drops are still surfaced via `FilterCounts.secondary` /
/// `.supplementary` so users can audit how many records were
/// removed, but there is no configuration that lets them through.
/// A future BAM (or other) input path should mirror this policy, or
/// document a strong reason and a redesign of the walker contract
/// before diverging.
#[derive(Debug, Clone, Copy)]
pub struct AlignmentMergedReaderConfig {
    /// `None` = no minimum; `Some(n)` = drop reads with MAPQ < n.
    pub min_mapq: Option<u8>,
    /// `None` = no minimum; `Some(n)` = drop reads with decoded SEQ
    /// length < n.
    pub min_read_length: Option<u32>,
    /// Drop reads with `flag & 0x200` set.
    pub drop_qc_fail: bool,
    /// Drop reads with `flag & 0x400` set.
    pub drop_duplicate: bool,
    /// Default-on defence against contamination, adapter readthrough,
    /// and other "the aligner placed it but it doesn't really belong
    /// here" failure modes. `None` = filter disabled; `Some(x)` = drop
    /// a read whose `M`-op mismatch fraction exceeds `x`.
    ///
    /// Mismatch fraction = `mismatches_with_bq_above_floor / m_op_bases_atgc`,
    /// counting only positions where both the read base and the
    /// reference base are in `{A, C, G, T}` (positions with `N` in
    /// either are skipped from both numerator and denominator).
    /// `mismatch_bq_floor` controls which mismatches count.
    ///
    /// Default `Some(0.10)` (`DEFAULT_MAX_READ_MISMATCH_FRACTION`).
    /// Real BAQ-adjusted Illumina at MAPQ 20+ sits at ~0.5-1% mismatch
    /// rate; 10% is comfortably outside that distribution. See finding
    /// `F1` in
    /// `ia/reviews/pileup_freebayes_comparison_2026-05-08.md`.
    pub max_read_mismatch_fraction: Option<f32>,
    /// BQ floor below which a mismatch does not count toward
    /// `max_read_mismatch_fraction`. Default `10` (`DEFAULT_MISMATCH_BQ_FLOOR`),
    /// matching freebayes' `BQL2`. `0` disables the floor (every
    /// mismatch counts). Only meaningful when
    /// `max_read_mismatch_fraction` is `Some`.
    ///
    /// At the `alignment_input` stage the BQ is the raw value from the
    /// CRAM, *not* BAQ-adjusted (BAQ runs in a later stage). Filtering
    /// on raw BQ is the correct level for this filter — we are
    /// rejecting whole reads, not per-base evidence, and want to do
    /// so before paying the BAQ cost.
    pub mismatch_bq_floor: u8,
}

// SAM/BAM flag bit constants — used both inside this module and by
// tests building synthetic records.
pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE_STRAND: u16 = 0x10;
pub const FLAG_MATE_REVERSE_STRAND: u16 = 0x20;
pub const FLAG_FIRST_OF_PAIR: u16 = 0x40;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_QC_FAIL: u16 = 0x200;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

impl Default for AlignmentMergedReaderConfig {
    fn default() -> Self {
        Self {
            min_mapq: Some(DEFAULT_MIN_MAPQ),
            min_read_length: Some(DEFAULT_MIN_READ_LENGTH),
            drop_qc_fail: true,
            drop_duplicate: true,
            max_read_mismatch_fraction: Some(DEFAULT_MAX_READ_MISMATCH_FRACTION),
            mismatch_bq_floor: DEFAULT_MISMATCH_BQ_FLOOR,
        }
    }
}

// ---------------------------------------------------------------------
// Per-input lazy record stream
// ---------------------------------------------------------------------

/// One alignment file's worth of already-decoded records, plus the
/// path string used for error messages.
///
/// `pub(crate)` on purpose — `RecordBuf` is a noodles type, and
/// exposing it on a public type would re-couple our public API to
/// noodles. Tests live in the same crate, so `pub(crate)` is enough.
pub(crate) struct OpenAlignmentFile {
    pub path: PathBuf,
    pub records: AlignmentRecordsIter,
}

/// Type alias for the per-input record stream the merge consumes.
/// Both [`crate::bam::cram_input`] (today) and the planned
/// `crate::bam::bam_input` sibling produce values of this shape.
pub(super) type AlignmentRecordsIter =
    Box<dyn Iterator<Item = io::Result<sam::alignment::RecordBuf>> + Send>;

// ---------------------------------------------------------------------
// Per-format record-stream decoders live in sibling submodules
// ---------------------------------------------------------------------
// CRAM-side: [`crate::bam::cram_input`] exports
// `open_cram_record_stream` / `open_indexed_cram_record_stream`.
// BAM-side: [`crate::bam::bam_input`] exports
// `open_bam_record_stream` / `open_indexed_bam_record_stream`.
// Both produce an `AlignmentRecordsIter`; this module wraps the
// result into [`OpenAlignmentFile`] inside `new` and `query`.

// ---------------------------------------------------------------------
// Header validation
// ---------------------------------------------------------------------

/// Per-CRAM header summary, built once during pre-flight validation.
struct AlignmentFileHeaderSummary {
    contigs: ContigList,
    sample_name: String,
}

/// Pull `@HD SO`, `@SQ` list, and the (single) `@RG SM` value out of a
/// `sam::Header`. Returns `AlignmentInputError` if any required invariant
/// fails: SO != coordinate, missing SM, multiple distinct SMs in the
/// same file.
fn extract_header(
    path: &Path,
    sam_header: &sam::Header,
) -> Result<AlignmentFileHeaderSummary, AlignmentInputError> {
    // Sort order.
    let sort_order = sort_order_string(sam_header);
    if sort_order.as_deref() != Some("coordinate") {
        return Err(AlignmentInputError::NotCoordinateSorted {
            path: path.to_path_buf(),
            sort_order: sort_order.unwrap_or_else(|| "<missing>".into()),
        });
    }

    // @SQ list.
    let mut entries: Vec<ContigEntry> = Vec::new();
    for (name_bstr, ref_seq_map) in sam_header.reference_sequences() {
        let name = String::from_utf8_lossy(name_bstr.as_ref()).into_owned();
        let length: u64 = usize::from(ref_seq_map.length()) as u64;
        let md5 = md5_from_reference_sequence(path, &name, ref_seq_map)?;
        entries.push(ContigEntry { name, length, md5 });
    }
    let contigs = ContigList { entries };

    // @RG SM. All SMs in this single file must collapse to one value.
    let sample_name = extract_single_sample_name(path, sam_header)?;

    Ok(AlignmentFileHeaderSummary {
        contigs,
        sample_name,
    })
}

fn sort_order_string(sam_header: &sam::Header) -> Option<String> {
    use noodles_sam::header::record::value::map::header::tag::SORT_ORDER;
    let hd = sam_header.header()?;
    let raw = hd.other_fields().get(&SORT_ORDER)?;
    Some(String::from_utf8_lossy(raw.as_ref()).into_owned())
}

fn md5_from_reference_sequence(
    path: &Path,
    contig_name: &str,
    ref_seq_map: &sam::header::record::value::Map<
        sam::header::record::value::map::ReferenceSequence,
    >,
) -> Result<Option<[u8; 16]>, AlignmentInputError> {
    use noodles_sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
    let Some(raw) = ref_seq_map.other_fields().get(&MD5_CHECKSUM) else {
        return Ok(None);
    };
    let bytes = raw.as_ref();
    decode_md5_hex(bytes)
        .map(Some)
        .ok_or_else(|| AlignmentInputError::MalformedMd5 {
            path: path.to_path_buf(),
            contig: contig_name.to_string(),
            detail: if bytes.len() != 32 {
                format!("expected 32 hex chars, got {} bytes", bytes.len())
            } else {
                "non-hex character in M5 value".into()
            },
        })
}

fn decode_md5_hex(hex: &[u8]) -> Option<[u8; 16]> {
    if hex.len() != 32 {
        return None;
    }
    let mut out = [0u8; 16];
    for i in 0..16 {
        let hi = hex_nibble(hex[2 * i])?;
        let lo = hex_nibble(hex[2 * i + 1])?;
        out[i] = (hi << 4) | lo;
    }
    Some(out)
}

fn hex_nibble(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'a'..=b'f' => Some(b - b'a' + 10),
        b'A'..=b'F' => Some(b - b'A' + 10),
        _ => None,
    }
}

fn extract_single_sample_name(
    path: &Path,
    sam_header: &sam::Header,
) -> Result<String, AlignmentInputError> {
    use noodles_sam::header::record::value::map::read_group::tag::SAMPLE;
    let mut current: Option<(String, String)> = None;
    for (rg_id_bstr, rg_map) in sam_header.read_groups() {
        let rg_id = String::from_utf8_lossy(rg_id_bstr.as_ref()).into_owned();
        let sm_raw = rg_map.other_fields().get(&SAMPLE).ok_or_else(|| {
            AlignmentInputError::MissingSampleTag {
                path: path.to_path_buf(),
                read_group_id: rg_id.clone(),
            }
        })?;
        let sm = String::from_utf8_lossy(sm_raw.as_ref()).into_owned();
        match &current {
            None => current = Some((rg_id, sm)),
            Some((existing_rg, existing_sm)) if existing_sm != &sm => {
                return Err(AlignmentInputError::MultipleSampleNamesInFile {
                    path: path.to_path_buf(),
                    rg_a: existing_rg.clone(),
                    sm_a: existing_sm.clone(),
                    rg_b: rg_id,
                    sm_b: sm,
                });
            }
            Some(_) => {}
        }
    }
    current
        .map(|(_, sm)| sm)
        .ok_or_else(|| AlignmentInputError::MissingSampleTag {
            path: path.to_path_buf(),
            read_group_id: "<no @RG entries>".into(),
        })
}

// ---------------------------------------------------------------------
// FASTA agreement
// ---------------------------------------------------------------------

/// Validate that the indexed FASTA's contigs match `contigs`
/// (name + length, in the same order). The alignment file's
/// `@SQ M5` is trusted; we do not recompute it from the FASTA.
fn validate_fasta_agreement(
    fasta_path: &Path,
    canonical_contigs: &ContigList,
    alignment_file_path: &Path,
) -> Result<(), AlignmentInputError> {
    let fai_path = with_fai_extension(fasta_path);
    if !fai_path.exists() {
        return Err(AlignmentInputError::MissingFastaIndex {
            fasta_path: fasta_path.to_path_buf(),
        });
    }
    let index =
        noodles_fasta::fai::fs::read(&fai_path).map_err(|source| AlignmentInputError::Io {
            path: fai_path.clone(),
            source,
        })?;
    let fai_records: &[noodles_fasta::fai::Record] = index.as_ref();

    if fai_records.len() != canonical_contigs.entries.len() {
        return Err(AlignmentInputError::FastaContigMismatch {
            fasta_path: fasta_path.to_path_buf(),
            alignment_file_path: alignment_file_path.to_path_buf(),
            detail: format!(
                "FASTA has {} contigs, alignment file has {}",
                fai_records.len(),
                canonical_contigs.entries.len()
            ),
        });
    }
    for (i, (fai_record, contig)) in fai_records
        .iter()
        .zip(canonical_contigs.entries.iter())
        .enumerate()
    {
        let fai_name = String::from_utf8_lossy(fai_record.name().as_ref()).into_owned();
        if fai_name != contig.name {
            return Err(AlignmentInputError::FastaContigMismatch {
                fasta_path: fasta_path.to_path_buf(),
                alignment_file_path: alignment_file_path.to_path_buf(),
                detail: format!(
                    "name disagreement at index {} ('{}' in FASTA vs '{}' in alignment file)",
                    i, fai_name, contig.name
                ),
            });
        }
        let fai_length: u64 = fai_record.length();
        if fai_length != contig.length {
            return Err(AlignmentInputError::FastaContigMismatch {
                fasta_path: fasta_path.to_path_buf(),
                alignment_file_path: alignment_file_path.to_path_buf(),
                detail: format!(
                    "length disagreement at index {} (contig '{}': {} in FASTA vs {} in alignment file)",
                    i, contig.name, fai_length, contig.length
                ),
            });
        }
    }
    Ok(())
}

fn with_fai_extension(fasta_path: &Path) -> PathBuf {
    let mut buf = fasta_path.as_os_str().to_owned();
    buf.push(".fai");
    PathBuf::from(buf)
}

// ---------------------------------------------------------------------
// AlignmentMergedReader
// ---------------------------------------------------------------------

/// A genome location: index into the merged `ContigList` and 1-based
/// position on that contig.
type Locus = (usize, u64);

/// Identifying fingerprint of a mapped read: the four fields whose
/// equality defines "same read" for duplicate detection — qname,
/// SAM flags, contig-list index, and 1-based position. See
/// `ia/specs/per_sample_pileup.md` §"Duplicate-read detection across
/// CRAMs" for the rationale behind these specific fields.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ReadFingerprint {
    qname: Vec<u8>,
    flag: u16,
    ref_id: usize,
    pos: u64,
}

/// A `ReadFingerprint` paired with the index of the input file it
/// came from. Stored once per accepted read at the current locus so
/// a duplicate error can name *both* source files (the previous
/// acceptance and the colliding one).
#[derive(Debug, Clone)]
struct ReadFingerprintWithSourceFile {
    key: ReadFingerprint,
    source_file_index: usize,
}

type AlignmentPeekable =
    BufferedPeekable<AlignmentRecordsIter, sam::alignment::RecordBuf, io::Error>;

pub struct AlignmentMergedReader {
    record_streams: Vec<AlignmentPeekable>,
    paths: Vec<PathBuf>,
    contigs: ContigList,
    sample_name: String,
    config: AlignmentMergedReaderConfig,
    /// FASTA repository for the F1 mismatch-fraction filter. `None`
    /// disables that filter regardless of `config.max_read_mismatch_fraction`
    /// — only relevant for in-memory test fixtures that do not have a
    /// reference. The production `new()` always passes `Some`.
    repository: Option<fasta::Repository>,
    /// Cache of the most-recently-fetched contig sequence, indexed by
    /// `ref_id`. Avoids one `Repository::get` (and its `RwLock`
    /// acquisition) per accepted read on the steady-state path of a
    /// chromosome. Cleared and refilled when the read's `ref_id`
    /// changes. Held as `Arc<fasta::record::Sequence>` because that's
    /// what `Repository::get` returns and cloning the `Arc` is free.
    cached_contig: Option<(usize, std::sync::Arc<fasta::record::Sequence>)>,
    /// Previous `Locus` accepted from each input, parallel to
    /// `record_streams` and `paths`. `None` until the first record
    /// from that input is accepted. Used to detect within-file
    /// regressions in coordinate order.
    per_file_prev_locus: Vec<Option<Locus>>,
    /// Cross-file duplicate-detection buffer for the current locus:
    /// the `ReadFingerprintWithSourceFile` of every read already
    /// accepted at `current_locus`. A new candidate whose fingerprint
    /// matches any entry here is rejected as a cross-file duplicate.
    /// Cleared whenever the merge advances to a new locus, since
    /// duplicates only matter between records at the same
    /// `(ref_id, pos)`.
    current_locus_read_fingerprints: Vec<ReadFingerprintWithSourceFile>,
    /// The locus the merge is currently processing. `None` before the
    /// first record is examined; otherwise the locus of every entry
    /// in `current_locus_read_fingerprints` (which may be empty if no
    /// record has been accepted at this locus yet).
    current_locus: Option<Locus>,
    filter_counts: FilterCounts,
    /// Set the first time `next()` returns `Some(Err(_))`. After that,
    /// every subsequent call returns `None`. Implements the fuse-on-
    /// error semantics required by the spec's "halts Stage 1"
    /// language so callers can use `for`, `collect`, or `try_fold`
    /// without separate stop-on-first-error bookkeeping.
    fused: bool,
}

impl AlignmentMergedReader {
    pub fn new(
        alignment_files: &[PathBuf],
        fasta: &Path,
        config: AlignmentMergedReaderConfig,
    ) -> Result<Self, AlignmentInputError> {
        if alignment_files.is_empty() {
            return Err(AlignmentInputError::NoInputs);
        }

        // FASTA repository — built once and shared across decoders.
        let fai_path = with_fai_extension(fasta);
        if !fai_path.exists() {
            return Err(AlignmentInputError::MissingFastaIndex {
                fasta_path: fasta.to_path_buf(),
            });
        }
        let indexed_fasta_reader = noodles_fasta::io::indexed_reader::Builder::default()
            .build_from_path(fasta)
            .map_err(|source| AlignmentInputError::Io {
                path: fasta.to_path_buf(),
                source,
            })?;
        let adapter = fasta::repository::adapters::IndexedReader::new(indexed_fasta_reader);
        let repository = fasta::Repository::new(adapter);

        // Reject mixed `.cram` + `.bam` inputs and unknown
        // extensions before touching the filesystem. The
        // `index_preflight::preflight_alignment_indexes` layer
        // enforces the same policy for the indexed-query path;
        // `pileup` does not run through pre-flight so the gate
        // lives here too. Same `display_name` strings on both
        // sides so the error message is consistent.
        // Classify-pass: harvest the per-input kind once. The
        // open-pass below zips this Vec back in, eliminating the
        // need to re-call `AlignmentFileKind::from_path` (and the
        // `.unwrap()` that pattern invited).
        let mut input_kinds: Vec<AlignmentFileKind> = Vec::with_capacity(alignment_files.len());
        let mut first_seen_kind: Option<(usize, AlignmentFileKind)> = None;
        for (idx, input_path) in alignment_files.iter().enumerate() {
            let kind = AlignmentFileKind::from_path(input_path).ok_or_else(|| {
                AlignmentInputError::UnsupportedExtension {
                    path: input_path.clone(),
                }
            })?;
            match first_seen_kind {
                None => first_seen_kind = Some((idx, kind)),
                Some((first_idx, first_kind)) if first_kind != kind => {
                    return Err(AlignmentInputError::MixedAlignmentFileFormats {
                        first_path: alignment_files[first_idx].clone(),
                        first_format: first_kind.display_name(),
                        other_path: input_path.clone(),
                        other_format: kind.display_name(),
                    });
                }
                Some(_) => {}
            }
            input_kinds.push(kind);
        }

        // Open every input, validate per-file invariants, and
        // collect per-file headers.
        let mut open_alignment_files: Vec<OpenAlignmentFile> =
            Vec::with_capacity(alignment_files.len());
        let mut canonical_contigs: Option<ContigList> = None;
        let mut canonical_sample: Option<String> = None;
        let mut reference_input_path: Option<PathBuf> = None;

        for (input_path, kind) in alignment_files.iter().zip(input_kinds.iter().copied()) {
            // Format-specific decoder lives in the corresponding
            // sibling module; each helper opens the file, validates
            // any format-level invariants (CRAM version, BAM magic),
            // reads the SAM header, and returns the header + an
            // owned record iterator. Cross-file checks (contigs
            // identical, single SM tag) are this module's concern
            // and run after.
            let (noodles_sam_header, owned_records) = match kind {
                AlignmentFileKind::Cram => open_cram_record_stream(input_path, repository.clone())?,
                AlignmentFileKind::Bam => open_bam_record_stream(input_path)?,
            };
            let extracted_header = extract_header(input_path, &noodles_sam_header)?;

            // Cross-file checks: contigs and sample name must be
            // identical across every file in the input.
            match &canonical_contigs {
                None => {
                    canonical_contigs = Some(extracted_header.contigs.clone());
                    reference_input_path = Some(input_path.clone());
                }
                Some(existing) => {
                    if let Err(detail) = existing.first_disagreement(&extracted_header.contigs) {
                        return Err(AlignmentInputError::ContigListMismatch {
                            reference_path: reference_input_path.clone().expect(
                                "reference_input_path is set when canonical_contigs is Some",
                            ),
                            other_path: input_path.clone(),
                            detail,
                        });
                    }
                }
            }
            match &canonical_sample {
                None => canonical_sample = Some(extracted_header.sample_name.clone()),
                Some(existing) if existing != &extracted_header.sample_name => {
                    return Err(AlignmentInputError::MultipleSampleNames {
                        path_a: reference_input_path
                            .clone()
                            .expect("reference_input_path is set when canonical_sample is Some"),
                        sm_a: existing.clone(),
                        path_b: input_path.clone(),
                        sm_b: extracted_header.sample_name,
                    });
                }
                _ => {}
            }

            open_alignment_files.push(OpenAlignmentFile {
                path: input_path.clone(),
                records: owned_records,
            });
        }

        let canonical_contigs =
            canonical_contigs.expect("at least one CRAM was opened so canonical_contigs is set");
        let canonical_sample =
            canonical_sample.expect("at least one CRAM was opened so canonical_sample is set");

        // FASTA agreement — only after the canonical contig list is
        // known.
        validate_fasta_agreement(
            fasta,
            &canonical_contigs,
            reference_input_path.as_deref().unwrap_or(fasta),
        )?;

        Self::from_open_alignment_files(
            open_alignment_files,
            canonical_contigs,
            canonical_sample,
            config,
            Some(repository),
        )
    }

    /// Crate-internal core. Takes pre-built per-CRAM record streams
    /// plus the canonical contig list and sample name (the caller is
    /// responsible for cross-file header validation — `new` does this
    /// before calling here). All I/O for opening lives outside this
    /// function, which is why it is the single point of entry tests
    /// drive against.
    ///
    /// `repository` is the FASTA reference for the F1 mismatch-fraction
    /// filter. `None` disables that filter (intended for tests that
    /// build CRAM records in memory and have no reference). The
    /// production `new()` always passes `Some`.
    pub(crate) fn from_open_alignment_files(
        open_alignment_files: Vec<OpenAlignmentFile>,
        contigs: ContigList,
        sample_name: String,
        config: AlignmentMergedReaderConfig,
        repository: Option<fasta::Repository>,
    ) -> Result<Self, AlignmentInputError> {
        let n = open_alignment_files.len();
        let mut record_streams = Vec::with_capacity(n);
        let mut paths = Vec::with_capacity(n);
        for open in open_alignment_files {
            paths.push(open.path);
            record_streams.push(BufferedPeekable::with_buffer_size(
                open.records,
                PER_PEEKER_BUFFER_SIZE,
            ));
        }
        Ok(Self {
            record_streams,
            paths,
            contigs,
            sample_name,
            config,
            repository,
            cached_contig: None,
            per_file_prev_locus: vec![None; n],
            current_locus_read_fingerprints: Vec::new(),
            current_locus: None,
            filter_counts: FilterCounts::default(),
            fused: false,
        })
    }

    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }

    pub fn contigs(&self) -> &ContigList {
        &self.contigs
    }

    pub fn filter_counts(&self) -> &FilterCounts {
        &self.filter_counts
    }

    /// Per-contig random-access variant of [`Self::new`]. Opens each
    /// input CRAM fresh, seeks via the caller-supplied per-input
    /// `Arc<crai::Index>`, and yields only records whose alignment
    /// falls on `contig_name`. Records are k-way-merged across
    /// inputs in coordinate order via the same machinery as
    /// [`Self::new`].
    ///
    /// `contigs` is the canonical contig list shared across all
    /// input CRAMs (the driver harvests this once at startup via
    /// [`Self::new`] or an equivalent validation pass);
    /// `contig_name` is resolved against it to a ref_id.
    /// `sample_name` is propagated to the resulting reader's
    /// [`sample_name`](Self::sample_name) accessor.
    ///
    /// **What this skips vs `new`.** Per-file header parsing
    /// (caller's `Arc<sam::Header>` is reused), FASTA contig
    /// agreement, sample-name reconciliation, and cross-file
    /// contig-list cross-checks. Each per-worker call still opens
    /// fresh `cram::io::Reader<File>` and `fasta::Repository`
    /// instances because both are mutably held by the merge for the
    /// reader's lifetime, but their construction cost is small next
    /// to the per-contig seek + decode that follows.
    ///
    /// # Errors
    ///
    /// - [`AlignmentInputError::NoInputs`] — `crams` is empty.
    /// - [`AlignmentInputError::PerInputHandleCountMismatch`] —
    ///   `headers.len()` or `indexes.len()` disagree with
    ///   `crams.len()`.
    /// - [`AlignmentInputError::MissingFastaIndex`] — `<fasta>.fai` is
    ///   not present next to `fasta`.
    /// - [`AlignmentInputError::ContigNotInList`] — `contig_name` is
    ///   not present in `contigs`.
    /// - [`AlignmentInputError::OpenFailed`] / [`AlignmentInputError::Io`] —
    ///   any underlying CRAM or FASTA I/O failure during open.
    #[allow(clippy::too_many_arguments)]
    pub fn query(
        alignment_files: &[PathBuf],
        fasta: &Path,
        contigs: ContigList,
        sample_name: String,
        headers: &[Arc<sam::Header>],
        indexes: &[AlignmentIndex],
        contig_name: &str,
        config: AlignmentMergedReaderConfig,
    ) -> Result<Self, AlignmentInputError> {
        if alignment_files.is_empty() {
            return Err(AlignmentInputError::NoInputs);
        }
        if headers.len() != alignment_files.len() || indexes.len() != alignment_files.len() {
            return Err(AlignmentInputError::PerInputHandleCountMismatch {
                inputs: alignment_files.len(),
                headers: headers.len(),
                indexes: indexes.len(),
            });
        }

        let target_reference_sequence_id = contigs
            .entries
            .iter()
            .position(|c| c.name == contig_name)
            .ok_or_else(|| AlignmentInputError::ContigNotInList {
                contig: contig_name.to_string(),
                known_contigs: contigs.entries.len(),
            })?;

        // FASTA repository — needed by the slice decoder for the
        // mismatch-fraction filter. Built per-worker because
        // `fasta::Repository`'s internal adapter holds a file handle
        // that is not shareable across threads.
        let fai_path = with_fai_extension(fasta);
        if !fai_path.exists() {
            return Err(AlignmentInputError::MissingFastaIndex {
                fasta_path: fasta.to_path_buf(),
            });
        }
        let indexed_fasta_reader = noodles_fasta::io::indexed_reader::Builder::default()
            .build_from_path(fasta)
            .map_err(|source| AlignmentInputError::Io {
                path: fasta.to_path_buf(),
                source,
            })?;
        let adapter = fasta::repository::adapters::IndexedReader::new(indexed_fasta_reader);
        let repository = fasta::Repository::new(adapter);

        // One indexed iterator per input. Format-specific decoders
        // live in the `cram_input` / `bam_input` sibling modules;
        // this loop's job is just to pick the right opener and
        // match the per-input `AlignmentIndex` variant to the
        // input file's extension.
        //
        // The cohort driver is responsible for loading indexes via
        // [`crate::bam::index_preflight::load_alignment_index`],
        // which projects the on-disk format onto the right enum
        // variant. If the variants disagree with the file
        // extensions here, that's a programmer error in the
        // driver — surfaced as a typed
        // `AlignmentIndexFormatMismatch` rather than a panic.
        let mut open_alignment_files: Vec<OpenAlignmentFile> =
            Vec::with_capacity(alignment_files.len());
        for ((input_path, header), alignment_index) in alignment_files
            .iter()
            .zip(headers.iter())
            .zip(indexes.iter())
        {
            let file_kind = AlignmentFileKind::from_path(input_path).ok_or_else(|| {
                AlignmentInputError::UnsupportedExtension {
                    path: input_path.clone(),
                }
            })?;
            let owned_records = match (file_kind, alignment_index) {
                (AlignmentFileKind::Cram, AlignmentIndex::Crai(idx)) => {
                    open_indexed_cram_record_stream(
                        input_path,
                        Arc::clone(header),
                        Arc::clone(idx),
                        repository.clone(),
                        target_reference_sequence_id,
                    )?
                }
                (AlignmentFileKind::Bam, AlignmentIndex::BamCsi(idx)) => {
                    open_indexed_bam_record_stream(
                        input_path,
                        Arc::clone(header),
                        BamIndex::Csi(Arc::clone(idx)),
                        target_reference_sequence_id,
                    )?
                }
                (AlignmentFileKind::Bam, AlignmentIndex::BamBai(idx)) => {
                    open_indexed_bam_record_stream(
                        input_path,
                        Arc::clone(header),
                        BamIndex::Bai(Arc::clone(idx)),
                        target_reference_sequence_id,
                    )?
                }
                // Genuine driver-bug mismatches — enumerate each
                // explicitly so adding a new AlignmentIndex variant
                // forces a new arm here (the enum is
                // `#[non_exhaustive]`; without enumeration the
                // catch-all would silently absorb the new variant
                // as a "format mismatch", defeating the
                // extensibility intent).
                (
                    AlignmentFileKind::Cram,
                    index @ (AlignmentIndex::BamCsi(_) | AlignmentIndex::BamBai(_)),
                )
                | (AlignmentFileKind::Bam, index @ AlignmentIndex::Crai(_)) => {
                    return Err(AlignmentInputError::AlignmentIndexFormatMismatch {
                        path: input_path.clone(),
                        file_format: file_kind.display_name(),
                        index_format: index.display_name(),
                    });
                }
            };
            open_alignment_files.push(OpenAlignmentFile {
                path: input_path.clone(),
                records: owned_records,
            });
        }

        Self::from_open_alignment_files(
            open_alignment_files,
            contigs,
            sample_name,
            config,
            Some(repository),
        )
    }
}

// ---------------------------------------------------------------------
// Iterator + merge logic
// ---------------------------------------------------------------------

/// Iterates over surviving reads in coordinate order.
///
/// **Fuse-on-error semantics.** Once `next()` returns `Some(Err(_))`,
/// every subsequent call returns `None`. This matches the spec's
/// "halts Stage 1" requirement (`ia/specs/per_sample_pileup.md`
/// §"Errors") and lets callers use the iterator with any consumer
/// (`for`, `collect`, `try_fold`) without separate "stop on first
/// error" bookkeeping.
impl Iterator for AlignmentMergedReader {
    type Item = Result<MappedRead, AlignmentInputError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.fused {
            return None;
        }
        loop {
            // Refill record streams' heads, dropping any head that fails the
            // pre-decode filter cascade. Each dropped head increments
            // its `FilterCounts` bucket.
            if let Err(e) = self.refill_heads() {
                return self.fail(e);
            }

            // Pick the smallest surviving head; `None` means every
            // record stream is exhausted and we are done.
            let chosen_idx = self.argmin_head()?;

            // Validate per-file order against this file's previous
            // accepted `(ref_id, pos)`.
            let head = match self.peek_head_keys(chosen_idx) {
                Ok(Some(keys)) => keys,
                Ok(None) => continue,
                Err(e) => return self.fail(e),
            };
            if let Some((prev_ref, prev_pos)) = self.per_file_prev_locus[chosen_idx]
                && (head.ref_id, head.pos) < (prev_ref, prev_pos)
            {
                return self.fail(AlignmentInputError::OutOfOrderRead {
                    path: self.paths[chosen_idx].clone(),
                    qname: String::from_utf8_lossy(&head.qname).into_owned(),
                    prev_ref_id: prev_ref,
                    prev_pos,
                    this_ref_id: head.ref_id,
                    this_pos: head.pos,
                });
            }

            // Move `current_locus` to the head's locus, clearing the
            // accepted-fingerprints buffer if we've advanced past the
            // previous one.
            self.advance_current_locus_if_needed(head.ref_id, head.pos);

            // Check for cross-file duplicates at the current locus.
            let new_key = ReadFingerprint {
                qname: head.qname.clone(),
                flag: head.flag,
                ref_id: head.ref_id,
                pos: head.pos,
            };
            if let Some(other) = self
                .current_locus_read_fingerprints
                .iter()
                .find(|entry| entry.key == new_key)
            {
                let other_path = self.paths[other.source_file_index].clone();
                return self.fail(AlignmentInputError::DuplicateReadAcrossFiles {
                    qname: String::from_utf8_lossy(&new_key.qname).into_owned(),
                    path_a: other_path,
                    path_b: self.paths[chosen_idx].clone(),
                    ref_id: new_key.ref_id,
                    pos: new_key.pos,
                });
            }

            // Consume the chosen record stream's head.
            let record = match self.record_streams[chosen_idx].next() {
                Some(Ok(rb)) => rb,
                Some(Err(e)) => {
                    return self.fail(AlignmentInputError::MalformedRecord {
                        path: self.paths[chosen_idx].clone(),
                        qname: Some(String::from_utf8_lossy(&head.qname).into_owned()),
                        source: e,
                    });
                }
                None => continue,
            };

            // Min-read-length filter — applied *before* the
            // RecordBuf → MappedRead allocation so short reads do not
            // pay for the qname / cigar / seq / qual clones we are
            // about to throw away. The order/dup checks above have
            // already run, so the dedup buffer at the current locus
            // stays consistent.
            if let Some(min) = self.config.min_read_length
                && (record.sequence().as_ref().len() as u32) < min
            {
                self.filter_counts.too_short += 1;
                continue;
            }

            // Convert RecordBuf -> MappedRead.
            let mut mapped = match record_buf_to_mapped_read(&record, chosen_idx) {
                Ok(m) => m,
                Err(e) => {
                    return self.fail(AlignmentInputError::MalformedRecord {
                        path: self.paths[chosen_idx].clone(),
                        qname: Some(String::from_utf8_lossy(&head.qname).into_owned()),
                        source: e,
                    });
                }
            };

            // G2 — `GoodCigar`-style read rejection. Drops reads
            // whose CIGAR contains an adjacent `I`/`D` pair (no
            // biological event produces this) or whose first/last
            // op (after stripping clips) is a deletion (no flanking
            // evidence on the missing side). Default-on, no opt-out.
            //
            // Runs **before** F3 because F3 can deliberately shift
            // an interior indel to a boundary as part of normal
            // canonicalisation; checking after F3 would punish that.
            // Pre-F3 boundary deletions are the case where the
            // aligner itself produced the suspect alignment.
            if cigar_is_bad(&mapped.cigar) {
                self.filter_counts.bad_cigar += 1;
                continue;
            }

            // F3 + F1 both need the read's reference slice. Fetch
            // once and share. Repository is cheap to clone — its
            // internal state is `Arc<RwLock<...>>` — and cloning
            // releases the immutable borrow on `self.repository`
            // before `fetch_ref_for_read` takes `&mut self` for its
            // cache update.
            let ref_seq_opt = self
                .repository
                .clone()
                .and_then(|repo| self.fetch_ref_for_read(&repo, &mapped));

            // F3 — left-align indels in the CIGAR to their leftmost
            // (canonical) position, so reads carrying the same biological
            // indel converge on one `(anchor, REF, ALT)` and consolidate
            // downstream. Always on, no opt-out. Mutates `mapped.cigar` in
            // place; read seq and qual bytes are unchanged (only their
            // alignment annotation moves). Skipped when no `Repository` is
            // plumbed (test fixtures). Backed by the GATK `leftAlignIndels`
            // port in `pileup::walker::indel_norm` (architecture spec
            // §"Indel normalization (left-alignment)"); see that module for
            // the algorithm.
            if let Some(ref_seq) = ref_seq_opt.as_ref() {
                left_align_indels(&mut mapped.cigar, &mapped.seq, ref_seq);
            }

            // F1 — per-read mismatch-fraction filter. Drops reads
            // whose `M`-op mismatch fraction exceeds the configured
            // threshold; defends against contamination, adapter
            // readthrough, and chimeric tails that pass MAPQ but
            // wouldn't survive base-by-base scrutiny against the
            // reference. Runs *after* F3 so the mismatch fraction
            // is computed against the canonicalised CIGAR. Skipped
            // when the threshold is `None` or when no `Repository`
            // is plumbed.
            if let Some(threshold) = self.config.max_read_mismatch_fraction
                && let Some(ref_seq) = ref_seq_opt.as_ref()
                && read_exceeds_mismatch_fraction(
                    &mapped.cigar,
                    &mapped.seq,
                    &mapped.qual,
                    ref_seq,
                    self.config.mismatch_bq_floor,
                    threshold,
                )
            {
                self.filter_counts.high_mismatch_fraction += 1;
                continue;
            }

            // Accept: record the fingerprint and bump the per-file
            // previous-locus tracker.
            self.current_locus_read_fingerprints
                .push(ReadFingerprintWithSourceFile {
                    key: new_key,
                    source_file_index: chosen_idx,
                });
            self.per_file_prev_locus[chosen_idx] = Some((mapped.ref_id, mapped.pos));
            return Some(Ok(mapped));
        }
    }
}

impl std::iter::FusedIterator for AlignmentMergedReader {}

impl AlignmentMergedReader {
    fn fail(&mut self, e: AlignmentInputError) -> Option<<Self as Iterator>::Item> {
        self.fused = true;
        Some(Err(e))
    }

    /// Fetch the reference slice covering `[mapped.pos, mapped.pos +
    /// cigar_ref_span)` for the F1 mismatch-fraction filter. Returns
    /// `None` if the reference is unavailable or the slice would
    /// extend past the end of the contig — both treated as "no
    /// signal", which causes the F1 caller to skip the read instead
    /// of dropping it.
    ///
    /// Caches one chromosome's `Arc<Sequence>` between calls so
    /// reads contiguous on the same chromosome avoid the
    /// `Repository::get` `RwLock` acquisition per read. The
    /// underlying `noodles_fasta::Repository` already caches whole
    /// contigs unboundedly, so this on-top cache is only an
    /// optimisation against the per-call dispatch cost.
    fn fetch_ref_for_read(
        &mut self,
        repository: &fasta::Repository,
        mapped: &MappedRead,
    ) -> Option<Vec<u8>> {
        let ref_span = cigar_ref_span(&mapped.cigar);
        if ref_span == 0 {
            return None;
        }
        let contig = self.contigs.entries.get(mapped.ref_id)?;

        // Refresh the cache on chrom change.
        let needs_refresh = match &self.cached_contig {
            Some((id, _)) => *id != mapped.ref_id,
            None => true,
        };
        if needs_refresh {
            let seq = repository.get(contig.name.as_bytes())?.ok()?;
            self.cached_contig = Some((mapped.ref_id, seq));
        }
        let seq = &self.cached_contig.as_ref().expect("just set").1;

        // `pos` is 1-based in our `MappedRead`; we want the byte
        // slice covering reference positions `[pos, pos+ref_span)`.
        // `Sequence` implements `AsRef<[u8]>` so we slice the raw
        // byte buffer directly — sidesteps `Position` construction
        // on the hot path.
        let raw: &[u8] = AsRef::<[u8]>::as_ref(seq.as_ref());
        let start = (mapped.pos as usize).checked_sub(1)?;
        let end = start.checked_add(ref_span as usize)?;
        if end > raw.len() {
            return None;
        }
        Some(raw[start..end].to_vec())
    }

    fn advance_current_locus_if_needed(&mut self, ref_id: usize, pos: u64) {
        match self.current_locus {
            Some(existing) if existing != (ref_id, pos) => {
                self.current_locus_read_fingerprints.clear();
                self.current_locus = Some((ref_id, pos));
            }
            None => {
                self.current_locus = Some((ref_id, pos));
            }
            _ => {}
        }
    }

    fn refill_heads(&mut self) -> Result<(), AlignmentInputError> {
        let config = self.config;
        for idx in 0..self.record_streams.len() {
            loop {
                let drop_outcome = match self.record_streams[idx].peek() {
                    Ok(Some(rb)) => classify_pre_decode(&config, rb),
                    Ok(None) => break,
                    Err(e) => {
                        return Err(AlignmentInputError::MalformedRecord {
                            path: self.paths[idx].clone(),
                            qname: None,
                            source: e,
                        });
                    }
                };
                match drop_outcome {
                    PreDecodeOutcome::Keep => break,
                    PreDecodeOutcome::Drop(bucket) => {
                        self.bump_filter_count(bucket);
                        // Consume the dropped head and re-peek.
                        let _ = self.record_streams[idx].next();
                    }
                }
            }
        }
        Ok(())
    }

    /// Inspect a peekable head without consuming it. Returns the
    /// merge-order key plus the few fields we need to perform the
    /// out-of-order and duplicate checks. `Ok(None)` means the
    /// record stream is exhausted (peek is None at this point — should not
    /// happen between `refill_heads` and `argmin_head`, but treat
    /// defensively).
    fn peek_head_keys(&mut self, idx: usize) -> Result<Option<HeadKey>, AlignmentInputError> {
        match self.record_streams[idx].peek() {
            Ok(Some(rb)) => match head_key(rb) {
                Some(key) => Ok(Some(key)),
                None => {
                    // Defence in depth: unmapped reads are now always
                    // filtered upstream by `classify_pre_decode`, so a
                    // record reaching `peek_head_keys` should always
                    // carry a `reference_sequence_id` and an
                    // `alignment_start`. If `head_key` returns `None`
                    // here, the input is genuinely malformed (the CRAM
                    // claims the read is mapped but omits both fields)
                    // and the user is better served by a typed error
                    // than a silent drop.
                    let qname = rb
                        .name()
                        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned());
                    Err(AlignmentInputError::MalformedRecord {
                        path: self.paths[idx].clone(),
                        qname,
                        source: io::Error::new(
                            io::ErrorKind::InvalidData,
                            "record has no reference_sequence_id or alignment_start",
                        ),
                    })
                }
            },
            Ok(None) => Ok(None),
            Err(e) => Err(AlignmentInputError::MalformedRecord {
                path: self.paths[idx].clone(),
                qname: None,
                source: e,
            }),
        }
    }

    fn argmin_head(&mut self) -> Option<usize> {
        let mut best: Option<(usize, Locus)> = None;
        for idx in 0..self.record_streams.len() {
            let Ok(Some(rb)) = self.record_streams[idx].peek() else {
                continue;
            };
            let Some(head) = head_key(rb) else {
                continue;
            };
            let candidate = (head.ref_id, head.pos);
            best = match best {
                None => Some((idx, candidate)),
                Some((_, cur)) if candidate < cur => Some((idx, candidate)),
                Some(b) => Some(b),
            };
        }
        best.map(|(idx, _)| idx)
    }

    fn bump_filter_count(&mut self, bucket: FilterBucket) {
        match bucket {
            FilterBucket::Unmapped => self.filter_counts.unmapped += 1,
            FilterBucket::Secondary => self.filter_counts.secondary += 1,
            FilterBucket::Supplementary => self.filter_counts.supplementary += 1,
            FilterBucket::QcFail => self.filter_counts.qc_fail += 1,
            FilterBucket::Duplicate => self.filter_counts.duplicate += 1,
            FilterBucket::LowMapq => self.filter_counts.low_mapq += 1,
            FilterBucket::BaqRejected => self.filter_counts.baq_rejected += 1,
        }
    }
}

fn classify_pre_decode(
    config: &AlignmentMergedReaderConfig,
    rb: &sam::alignment::RecordBuf,
) -> PreDecodeOutcome {
    let flag = rb.flags().bits();

    // Hit-rate-ordered cascade per
    // `ia/feature_implementation_plans/per_sample_caller_cram_input.md`
    // §"Filter precedence at each pulled record".

    // 1. Duplicate (~10-30%).
    if config.drop_duplicate && (flag & FLAG_DUPLICATE) != 0 {
        return PreDecodeOutcome::Drop(FilterBucket::Duplicate);
    }
    // 2. Low MAPQ (~5-15%).
    if let Some(min) = config.min_mapq {
        let mapq = rb.mapping_quality().map(u8::from).unwrap_or(0);
        if mapq < min {
            return PreDecodeOutcome::Drop(FilterBucket::LowMapq);
        }
    }
    // 3. Supplementary (~1-5%) — unconditionally dropped; see the
    //    "Why no drop_secondary / drop_supplementary toggles" note
    //    on `AlignmentMergedReaderConfig`.
    if (flag & FLAG_SUPPLEMENTARY) != 0 {
        return PreDecodeOutcome::Drop(FilterBucket::Supplementary);
    }
    // 4. Secondary (~1-5%) — unconditionally dropped; same rationale
    //    as supplementary.
    if (flag & FLAG_SECONDARY) != 0 {
        return PreDecodeOutcome::Drop(FilterBucket::Secondary);
    }
    // 5. Unmapped (~0-5%) — always dropped. An unmapped read contributes
    // no allele evidence and would otherwise trip the head-key invariant
    // (no `reference_sequence_id`, no `alignment_start`).
    if (flag & FLAG_UNMAPPED) != 0 {
        return PreDecodeOutcome::Drop(FilterBucket::Unmapped);
    }
    // 6. QC fail (<1%).
    if config.drop_qc_fail && (flag & FLAG_QC_FAIL) != 0 {
        return PreDecodeOutcome::Drop(FilterBucket::QcFail);
    }

    PreDecodeOutcome::Keep
}

#[derive(Debug, Clone)]
struct HeadKey {
    qname: Vec<u8>,
    flag: u16,
    ref_id: usize,
    pos: u64,
}

fn head_key(rb: &sam::alignment::RecordBuf) -> Option<HeadKey> {
    let ref_id = rb.reference_sequence_id()?;
    let pos = rb.alignment_start()?.get() as u64;
    let flag = rb.flags().bits();
    let qname = rb
        .name()
        .map(|n| AsRef::<[u8]>::as_ref(n).to_vec())
        .unwrap_or_default();
    Some(HeadKey {
        qname,
        flag,
        ref_id,
        pos,
    })
}

#[derive(Debug, Clone, Copy)]
enum PreDecodeOutcome {
    Keep,
    Drop(FilterBucket),
}

#[derive(Debug, Clone, Copy)]
pub(super) enum FilterBucket {
    Unmapped,
    Secondary,
    Supplementary,
    QcFail,
    Duplicate,
    LowMapq,
    /// Bucket for reads the BAQ stage refused — incremented by the
    /// pipeline integration in a later commit; no call site in
    /// `alignment_input` itself today.
    #[allow(dead_code)]
    BaqRejected,
}

// ---------------------------------------------------------------------
// Per-record conversion
// ---------------------------------------------------------------------

fn record_buf_to_mapped_read(
    rb: &sam::alignment::RecordBuf,
    source_file_index: usize,
) -> io::Result<MappedRead> {
    let qname = rb
        .name()
        .map(|n| AsRef::<[u8]>::as_ref(n).to_vec())
        .unwrap_or_default();
    let flag = rb.flags().bits();
    let ref_id = rb.reference_sequence_id().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "record has no reference_sequence_id",
        )
    })?;
    let pos = rb
        .alignment_start()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "record has no alignment_start"))?
        .get() as u64;
    let mapq = rb.mapping_quality().map(u8::from).unwrap_or(0);
    let cigar = cigar_to_ops(rb.cigar());
    let seq: Vec<u8> = rb
        .sequence()
        .as_ref()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let qual = rb.quality_scores().as_ref().to_vec();
    let mate_ref_id = rb.mate_reference_sequence_id();
    let mate_pos = rb.mate_alignment_start().map(|p| p.get() as u64);
    let template_length = rb.template_length();
    let adaptor_boundary = compute_adaptor_boundary(
        flag,
        ref_id,
        pos,
        seq.len() as u32,
        mate_ref_id,
        mate_pos,
        template_length,
    );
    Ok(MappedRead {
        qname,
        flag,
        ref_id,
        pos,
        mapq,
        cigar,
        seq,
        qual,
        mate_ref_id,
        mate_pos,
        adaptor_boundary,
        source_file_index,
    })
}

/// Finding `G1` (adaptor-region per-base filter) — pure-logic
/// boundary computation. Returns the 1-based reference position of
/// the first base that lies *inside* the mate-pair adaptor, or
/// `None` when the read is single-end, the mate is unreliable, the
/// fragment geometry is inconsistent, or the molecule is at least
/// as long as the read (so no read base could have been sequenced
/// past the molecule end).
///
/// Mirrors GATK's
/// [`ReadUtils.getAdaptorBoundary`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/read/ReadUtils.java#L527)
/// with one principled deviation: instead of GATK's hardcoded
/// `|tlen| > DEFAULT_ADAPTOR_SIZE = 100` cap (calibrated for 100bp
/// reads), we gate on `|tlen| < seq_len`. The molecule has to be
/// shorter than the *read's own* sequenced length for adaptor
/// readthrough to be possible — that condition is read-aware and
/// stays correct for ancient-DNA short fragments and for modern
/// 150bp+ reads alike.
///
/// On the forward strand the boundary is `read.start + |tlen|`
/// (first base past the molecule's 3′ end). On the reverse strand
/// it is `mate.start - 1` (last base before the molecule's 5′ end);
/// any read base at or before this position is inside the adaptor.
/// The walker's emit sites apply the direction-aware test using
/// `is_reverse_strand`.
pub(crate) fn compute_adaptor_boundary(
    flag: u16,
    ref_id: usize,
    pos: u64,
    seq_len: u32,
    mate_ref_id: Option<usize>,
    mate_pos: Option<u64>,
    template_length: i32,
) -> Option<u32> {
    // Required preconditions, in order of cheapness:
    if flag & FLAG_PAIRED == 0 {
        return None;
    }
    if flag & FLAG_UNMAPPED != 0 || flag & FLAG_MATE_UNMAPPED != 0 {
        return None;
    }
    if template_length == 0 {
        return None;
    }
    let is_reverse = flag & FLAG_REVERSE_STRAND != 0;
    let mate_is_reverse = flag & FLAG_MATE_REVERSE_STRAND != 0;
    // Same-strand pair → not a normal FR/RF pair, geometry is
    // not the simple readthrough case.
    if is_reverse == mate_is_reverse {
        return None;
    }
    let mate_ref_id = mate_ref_id?;
    if mate_ref_id != ref_id {
        return None;
    }
    let mate_pos = mate_pos?;

    let abs_tlen = template_length.unsigned_abs();
    // Per-read gate: only check when the read could physically have
    // run off the end of the molecule. If the molecule is at least
    // as long as the read's sequence, no base sequenced from this
    // read is in adaptor.
    if abs_tlen >= seq_len {
        return None;
    }

    if is_reverse {
        // Reverse-strand read: molecule's 5′ end is at the mate's
        // 1-based start; the first base inside adaptor (on the
        // read's 5′ side, which is the *higher* read offset since
        // the read is reversed) is the position immediately before
        // the mate's start. Skip when ref_pos <= boundary.
        if mate_pos == 0 {
            // Defensive: 1-based positions are non-zero in valid
            // SAM, but a corrupt record could place mate at 0.
            return None;
        }
        Some(mate_pos as u32 - 1)
    } else {
        // Forward-strand read: molecule's 3′ end is at
        // `read.start + |tlen|`. Skip when ref_pos >= boundary.
        Some(pos as u32 + abs_tlen)
    }
}

/// Reference-span sum of `M`-, `=`-, `X`-, `D`-, and `N`-op lengths
/// — i.e. the count of reference positions the read covers. Helper
/// for the F1 mismatch-fraction filter, which needs the slice
/// `[pos, pos + ref_span)` of the reference to compare against.
fn cigar_ref_span(cigar: &[CigarOp]) -> u32 {
    cigar
        .iter()
        .map(|op| match *op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Deletion(n)
            | CigarOp::Skip(n) => n,
            _ => 0,
        })
        .sum()
}

/// Pure-logic implementation of finding `G2` (`GoodCigar`-style
/// read-level rejection). Returns `true` when the CIGAR is
/// "ill-formed" by either of two rules — both sentinel-grade:
///
/// 1. **Adjacent `I`/`D` pair, in either order.** No biological
///    event produces an immediate insertion-then-deletion (or
///    deletion-then-insertion) in a single read; when the aligner
///    emits one, the alignment is genuinely confused and the read's
///    other M-events near that region are also unreliable.
/// 2. **Starts or ends with a deletion** (after stripping leading
///    soft- or hard-clips). A deletion at a CIGAR boundary has no
///    flanking evidence on the missing side; the aligner ran out of
///    read bases at the wrong moment.
///
/// The boundary-deletion rule must be applied to the **original**
/// CIGAR — i.e. before F3 left-alignment runs. F3 deliberately
/// shifts indels through homopolymers/tandem repeats and can in
/// principle land an indel at the read's edge as part of normal
/// canonicalisation; rejecting that case would punish F3 for doing
/// its job. Pre-F3 boundary deletions are the case where the
/// aligner itself produced the suspect alignment.
///
/// Adjacent indels are invariant under F3 (F3 explicitly skips
/// adjacent-indel triplets), so the consecutive-`I`/`D` rule is
/// timing-independent — but for symmetry we run both checks at the
/// same place in the cascade.
fn cigar_is_bad(cigar: &[CigarOp]) -> bool {
    // Rule 1 — adjacent I/D in either order.
    for window in cigar.windows(2) {
        let (a, b) = (window[0], window[1]);
        let a_is_indel = matches!(a, CigarOp::Insertion(_) | CigarOp::Deletion(_));
        let b_is_indel = matches!(b, CigarOp::Insertion(_) | CigarOp::Deletion(_));
        // Adjacent same-kind indels (II, DD) shouldn't occur in
        // canonical CIGARs either, but they are not the GATK
        // GoodCigar pattern — only mixed I/D pairs trigger this rule.
        let a_is_ins = matches!(a, CigarOp::Insertion(_));
        let b_is_ins = matches!(b, CigarOp::Insertion(_));
        if a_is_indel && b_is_indel && a_is_ins != b_is_ins {
            return true;
        }
    }

    // Rule 2 — first or last op (after stripping clips) is a
    // deletion. `iter().find` past the leading clips, and `iter()
    // .rev().find` past the trailing clips, give us the first/last
    // *non-clip* op. A read consisting entirely of clips has no
    // such op and is not a boundary-deletion case.
    let first_non_clip = cigar
        .iter()
        .find(|op| !matches!(op, CigarOp::SoftClip(_) | CigarOp::HardClip(_)));
    if matches!(first_non_clip, Some(CigarOp::Deletion(_))) {
        return true;
    }
    let last_non_clip = cigar
        .iter()
        .rev()
        .find(|op| !matches!(op, CigarOp::SoftClip(_) | CigarOp::HardClip(_)));
    if matches!(last_non_clip, Some(CigarOp::Deletion(_))) {
        return true;
    }

    false
}

/// Pure-logic implementation of finding `F1` (per-read
/// mismatch-fraction filter). Returns `true` when the read should
/// be dropped — i.e. its `M`-op mismatch fraction strictly exceeds
/// `threshold`.
///
/// - **Numerator:** count of `M`-op (or `=`/`X`-op) positions where
///   the read base differs from the reference base **and** the raw
///   base quality clears `bq_floor`. Positions where either base is
///   `N` (or non-ATGC) are skipped from both numerator and
///   denominator — they carry no usable signal.
/// - **Denominator:** count of `M`-op positions counted in the
///   numerator's denominator (i.e. ATGC-on-both-sides).
/// - If the denominator is `0` (e.g. a read consisting entirely of
///   soft-clips, or all `M` positions hit `N` bases), the function
///   returns `false` — there is no signal on which to base a drop
///   decision, and the conservative choice is to keep the read and
///   let downstream filters speak.
///
/// `seq`, `qual`, and `cigar` are the read's owned fields after
/// decoding; `ref_seq` is the reference slice covering
/// `[read.pos, read.pos + cigar_ref_span(cigar))`. Indexing into
/// these slices uses the standard CIGAR semantics — see
/// [`crate::pileup::walker::decompose`] for the
/// reference walk pattern this function mirrors.
fn read_exceeds_mismatch_fraction(
    cigar: &[CigarOp],
    seq: &[u8],
    qual: &[u8],
    ref_seq: &[u8],
    bq_floor: u8,
    threshold: f32,
) -> bool {
    let mut read_pos: usize = 0;
    let mut ref_pos: usize = 0;
    let mut mismatches: u32 = 0;
    let mut comparable_bases: u32 = 0;

    for op in cigar {
        match *op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                let n = n as usize;
                for k in 0..n {
                    let r = seq.get(read_pos + k).copied().unwrap_or(b'N');
                    let g = ref_seq.get(ref_pos + k).copied().unwrap_or(b'N');
                    let q = qual.get(read_pos + k).copied().unwrap_or(0);
                    let r_atgc = matches!(r, b'A' | b'C' | b'G' | b'T');
                    let g_atgc = matches!(g, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');
                    if r_atgc && g_atgc {
                        comparable_bases += 1;
                        if r != g.to_ascii_uppercase() && q >= bq_floor {
                            mismatches += 1;
                        }
                    }
                }
                read_pos += n;
                ref_pos += n;
            }
            CigarOp::Insertion(n) => {
                read_pos += n as usize;
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                ref_pos += n as usize;
            }
            CigarOp::SoftClip(n) => {
                read_pos += n as usize;
            }
            CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
        }
    }

    if comparable_bases == 0 {
        return false;
    }
    let fraction = (mismatches as f32) / (comparable_bases as f32);
    fraction > threshold
}

fn cigar_to_ops(cigar: &sam::alignment::record_buf::Cigar) -> Vec<CigarOp> {
    use sam::alignment::record::cigar::op::Kind;
    cigar
        .as_ref()
        .iter()
        .map(|op| {
            let len = op.len() as u32;
            match op.kind() {
                Kind::Match => CigarOp::Match(len),
                Kind::Insertion => CigarOp::Insertion(len),
                Kind::Deletion => CigarOp::Deletion(len),
                Kind::Skip => CigarOp::Skip(len),
                Kind::SoftClip => CigarOp::SoftClip(len),
                Kind::HardClip => CigarOp::HardClip(len),
                Kind::Pad => CigarOp::Padding(len),
                Kind::SequenceMatch => CigarOp::SeqMatch(len),
                Kind::SequenceMismatch => CigarOp::SeqMismatch(len),
            }
        })
        .collect()
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use noodles_cram as cram;

    use super::*;
    use crate::pileup::per_sample::record_specs::{
        RecordSpec, default_contigs, open_cram_from_records, record_spec,
    };

    fn default_seq(len: usize) -> Vec<u8> {
        b"A".repeat(len)
    }

    fn default_qual(len: usize) -> Vec<u8> {
        vec![30u8; len]
    }

    // --- Pure-type tests ---------------------------------------------

    #[test]
    fn p2_malformed_record_message_omits_empty_qname() {
        let with_qname = AlignmentInputError::MalformedRecord {
            path: PathBuf::from("/foo.cram"),
            qname: Some("R1".into()),
            source: io::Error::new(io::ErrorKind::InvalidData, "bad"),
        };
        assert!(
            with_qname.to_string().contains("qname='R1'"),
            "got {}",
            with_qname
        );
        let no_qname = AlignmentInputError::MalformedRecord {
            path: PathBuf::from("/foo.cram"),
            qname: None,
            source: io::Error::new(io::ErrorKind::InvalidData, "bad"),
        };
        assert!(
            !no_qname.to_string().contains("qname"),
            "qname clause leaked: {}",
            no_qname
        );
    }

    #[test]
    fn p1_contig_list_md5_wildcard_equality() {
        let with_md5 = ContigEntry {
            name: "chr1".into(),
            length: 100,
            md5: Some([1u8; 16]),
        };
        let same_md5 = with_md5.clone();
        let no_md5 = ContigEntry {
            name: "chr1".into(),
            length: 100,
            md5: None,
        };
        let different_md5 = ContigEntry {
            name: "chr1".into(),
            length: 100,
            md5: Some([2u8; 16]),
        };
        let different_name = ContigEntry {
            name: "chr2".into(),
            length: 100,
            md5: None,
        };
        let different_length = ContigEntry {
            name: "chr1".into(),
            length: 200,
            md5: None,
        };

        assert_eq!(with_md5, same_md5);
        assert_eq!(with_md5, no_md5);
        assert_eq!(no_md5, with_md5);
        assert_ne!(with_md5, different_md5);
        assert_ne!(with_md5, different_name);
        assert_ne!(with_md5, different_length);
    }

    // --- Group A: via from_open_alignment_files --------------------------------

    fn run_to_completion(
        mut reader: AlignmentMergedReader,
    ) -> (Vec<MappedRead>, FilterCounts, Option<AlignmentInputError>) {
        let mut out = Vec::new();
        let mut err = None;
        loop {
            match reader.next() {
                Some(Ok(r)) => out.push(r),
                Some(Err(e)) => {
                    err = Some(e);
                    break;
                }
                None => break,
            }
        }
        (out, *reader.filter_counts(), err)
    }

    fn pass_record(qname: &str, ref_id: usize, pos: u64) -> RecordSpec {
        // MAPQ above DEFAULT_MIN_MAPQ; clean flags; SEQ length above
        // DEFAULT_MIN_READ_LENGTH.
        let len = (DEFAULT_MIN_READ_LENGTH as usize) + 20;
        RecordSpec {
            qname: qname.into(),
            flag: 0,
            ref_id,
            pos,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(len as u32)],
            seq: default_seq(len),
            qual: default_qual(len),
            mate_ref_id: None,
            mate_pos: None,
        }
    }

    #[test]
    fn a1_single_stream_pass_through() {
        let cram_a = open_cram_from_records(
            "a.cram",
            vec![
                record_spec(pass_record("R1", 0, 100)),
                record_spec(pass_record("R2", 0, 200)),
                record_spec(pass_record("R3", 0, 300)),
            ],
        );
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none(), "unexpected error: {:?}", err);
        assert_eq!(out.len(), 3);
        assert_eq!(out[0].qname, b"R1");
        assert_eq!(out[1].qname, b"R2");
        assert_eq!(out[2].qname, b"R3");
        assert_eq!(out[0].pos, 100);
        assert_eq!(out[2].pos, 300);
        assert_eq!(out[0].source_file_index, 0);
        assert_eq!(counts, FilterCounts::default());
    }

    #[test]
    fn a2_multi_stream_merge_order() {
        let cram_a = open_cram_from_records(
            "a.cram",
            vec![
                record_spec(pass_record("A1", 0, 100)),
                record_spec(pass_record("A2", 0, 300)),
                record_spec(pass_record("A3", 0, 500)),
            ],
        );
        let cram_b = open_cram_from_records(
            "b.cram",
            vec![
                record_spec(pass_record("B1", 0, 150)),
                record_spec(pass_record("B2", 0, 200)),
                record_spec(pass_record("B3", 0, 400)),
            ],
        );
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a, cram_b],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, _, err) = run_to_completion(reader);
        assert!(err.is_none(), "unexpected error: {:?}", err);
        let positions: Vec<u64> = out.iter().map(|r| r.pos).collect();
        assert_eq!(positions, vec![100, 150, 200, 300, 400, 500]);
        let sources: Vec<usize> = out.iter().map(|r| r.source_file_index).collect();
        assert_eq!(sources, vec![0, 1, 1, 0, 1, 0]);
    }

    #[test]
    fn a3_tiebreaker_on_equal_coordinates() {
        let cram_a =
            open_cram_from_records("a.cram", vec![record_spec(pass_record("R_A", 0, 100))]);
        let cram_b =
            open_cram_from_records("b.cram", vec![record_spec(pass_record("R_B", 0, 100))]);
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a, cram_b],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, _, err) = run_to_completion(reader);
        assert!(err.is_none(), "unexpected error: {:?}", err);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].source_file_index, 0);
        assert_eq!(out[0].qname, b"R_A");
        assert_eq!(out[1].source_file_index, 1);
        assert_eq!(out[1].qname, b"R_B");
    }

    #[test]
    fn a4_out_of_order_within_a_single_stream() {
        let cram_a = open_cram_from_records(
            "a.cram",
            vec![
                record_spec(pass_record("R1", 0, 200)),
                record_spec(pass_record("R2", 0, 100)),
            ],
        );
        let mut reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let first = reader.next().expect("first").expect("ok");
        assert_eq!(first.pos, 200);
        let err = reader.next().expect("err item").expect_err("expected err");
        match err {
            AlignmentInputError::OutOfOrderRead {
                qname,
                prev_ref_id,
                prev_pos,
                this_ref_id,
                this_pos,
                ..
            } => {
                assert_eq!(qname, "R2");
                assert_eq!(prev_ref_id, 0);
                assert_eq!(prev_pos, 200);
                assert_eq!(this_ref_id, 0);
                assert_eq!(this_pos, 100);
            }
            other => panic!("unexpected error: {:?}", other),
        }
    }

    #[test]
    fn a4b_iterator_returns_none_after_first_error() {
        // Same shape as a4: out-of-order within a single stream.
        let cram_a = open_cram_from_records(
            "a.cram",
            vec![
                record_spec(pass_record("R1", 0, 200)),
                record_spec(pass_record("R2", 0, 100)),
            ],
        );
        let mut reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        assert!(reader.next().expect("ok item").is_ok()); // pos=200
        assert!(reader.next().expect("err item").is_err()); // OutOfOrderRead
        assert!(reader.next().is_none(), "iterator must fuse after Err");
        assert!(reader.next().is_none(), "fuse is sticky");
    }

    #[test]
    fn a5_duplicate_read_across_streams() {
        let dup = pass_record("R1", 0, 100);
        let cram_a = open_cram_from_records("a.cram", vec![record_spec(dup.clone())]);
        let cram_b = open_cram_from_records("b.cram", vec![record_spec(dup)]);
        let mut reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a, cram_b],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let first = reader.next().expect("first").expect("ok");
        assert_eq!(first.pos, 100);
        let err = reader.next().expect("err item").expect_err("expected err");
        match err {
            AlignmentInputError::DuplicateReadAcrossFiles {
                qname, ref_id, pos, ..
            } => {
                assert_eq!(qname, "R1");
                assert_eq!(ref_id, 0);
                assert_eq!(pos, 100);
            }
            other => panic!("unexpected error: {:?}", other),
        }
    }

    #[test]
    fn a6_dedup_buffer_clears_on_locus_advance() {
        let cram_a = open_cram_from_records(
            "a.cram",
            vec![
                record_spec(pass_record("R1", 0, 100)),
                record_spec(pass_record("R1", 0, 200)),
            ],
        );
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, _, err) = run_to_completion(reader);
        assert!(err.is_none(), "unexpected: {:?}", err);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].pos, 100);
        assert_eq!(out[1].pos, 200);
    }

    #[test]
    fn a7_min_mapq_filter() {
        assert_eq!(DEFAULT_MIN_MAPQ, 20);
        let mut below = pass_record("R_low", 0, 100);
        below.mapq = 5;
        let mut just_below = pass_record("R_just_low", 0, 200);
        just_below.mapq = DEFAULT_MIN_MAPQ - 1;
        let above = pass_record("R_high", 0, 300);
        let cram_a = open_cram_from_records(
            "a.cram",
            vec![
                record_spec(below),
                record_spec(just_below),
                record_spec(above),
            ],
        );
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].qname, b"R_high");
        assert_eq!(counts.low_mapq, 2);
    }

    #[test]
    fn a7b_missing_mapq_is_treated_as_zero_and_filtered_by_default() {
        // record_spec maps mapq=0 to "no MappingQuality set" because of
        // the `if spec.mapq > 0` guard, so this exercises the
        // mapping_quality()==None path. Pins the decision recorded in
        // the doc-comment on `DEFAULT_MIN_MAPQ`: a missing MAPQ is
        // treated as 0 and dropped under any non-zero minimum.
        let mut rec = pass_record("R", 0, 100);
        rec.mapq = 0;
        let cram = open_cram_from_records("a.cram", vec![record_spec(rec)]);
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert!(out.is_empty(), "MAPQ-missing read should be dropped");
        assert_eq!(counts.low_mapq, 1);
    }

    #[test]
    fn a8_min_mapq_none_disables() {
        let mut below = pass_record("R_low", 0, 100);
        below.mapq = 5;
        let mut just_below = pass_record("R_just_low", 0, 200);
        just_below.mapq = DEFAULT_MIN_MAPQ - 1;
        let above = pass_record("R_high", 0, 300);
        let cram_a = open_cram_from_records(
            "a.cram",
            vec![
                record_spec(below),
                record_spec(just_below),
                record_spec(above),
            ],
        );
        let cfg = AlignmentMergedReaderConfig {
            min_mapq: None,
            ..Default::default()
        };
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            cfg,
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert_eq!(out.len(), 3);
        assert_eq!(counts.low_mapq, 0);
    }

    fn six_flagged_records() -> Vec<RecordSpec> {
        let mut clean = pass_record("clean", 0, 100);
        let mut unmapped = pass_record("unmapped", 0, 200);
        unmapped.flag = FLAG_UNMAPPED;
        let mut secondary = pass_record("secondary", 0, 300);
        secondary.flag = FLAG_SECONDARY;
        let mut supplementary = pass_record("supplementary", 0, 400);
        supplementary.flag = FLAG_SUPPLEMENTARY;
        let mut qc_fail = pass_record("qc_fail", 0, 500);
        qc_fail.flag = FLAG_QC_FAIL;
        let mut duplicate = pass_record("duplicate", 0, 600);
        duplicate.flag = FLAG_DUPLICATE;
        // Note clean's MAPQ etc. are already valid; bump pos so it
        // sorts ahead of the others.
        clean.pos = 50;
        vec![
            clean,
            unmapped,
            secondary,
            supplementary,
            qc_fail,
            duplicate,
        ]
    }

    #[test]
    fn a9_each_flag_drop_one_at_a_time() {
        // Only `drop_qc_fail` and `drop_duplicate` are toggleable;
        // secondary/supplementary are unconditionally dropped (see
        // the policy note on `AlignmentMergedReaderConfig`) and so they
        // can't be isolated by toggling. Their behaviour is
        // covered by the all-defaults assertion at the bottom of
        // the test, which still expects each bucket count to be 1.
        type FlagCase = (
            fn(&mut AlignmentMergedReaderConfig),
            &'static str,
            fn(&FilterCounts) -> u64,
        );
        let cases: &[FlagCase] = &[
            (
                |c| {
                    c.drop_duplicate = false;
                },
                "qc_fail",
                |c| c.qc_fail,
            ),
            (
                |c| {
                    c.drop_qc_fail = false;
                },
                "duplicate",
                |c| c.duplicate,
            ),
        ];
        for (mutator, missing_qname, count_fn) in cases {
            let mut cfg = AlignmentMergedReaderConfig::default();
            mutator(&mut cfg);
            let recs: Vec<_> = six_flagged_records().into_iter().map(record_spec).collect();
            let cram_a = open_cram_from_records("a.cram", recs);
            let reader = AlignmentMergedReader::from_open_alignment_files(
                vec![cram_a],
                default_contigs(),
                "sample".into(),
                cfg,
                None,
            )
            .expect("reader");
            let (out, counts, err) = run_to_completion(reader);
            assert!(err.is_none(), "{}: {:?}", missing_qname, err);
            let qnames: Vec<&[u8]> = out.iter().map(|r| r.qname.as_slice()).collect();
            assert!(
                !qnames.contains(&missing_qname.as_bytes()),
                "{} should be filtered, got {:?}",
                missing_qname,
                qnames
            );
            assert_eq!(count_fn(&counts), 1, "filter {} count", missing_qname);
        }

        // All defaults: only `clean` survives, every bucket == 1.
        let recs: Vec<_> = six_flagged_records().into_iter().map(record_spec).collect();
        let cram_a = open_cram_from_records("a.cram", recs);
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none(), "{:?}", err);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].qname, b"clean");
        assert_eq!(counts.unmapped, 1);
        assert_eq!(counts.secondary, 1);
        assert_eq!(counts.supplementary, 1);
        assert_eq!(counts.qc_fail, 1);
        assert_eq!(counts.duplicate, 1);
    }

    #[test]
    fn a10_all_optional_flag_drops_disabled_keeps_only_unconditional_drops() {
        // Disabling every toggleable flag drop (`drop_qc_fail`,
        // `drop_duplicate`) leaves only the unconditional drops:
        // unmapped, secondary, and supplementary. Three of six
        // records survive (clean, qc_fail, duplicate); the other
        // three each tick exactly one filter bucket.
        let recs: Vec<_> = six_flagged_records().into_iter().map(record_spec).collect();
        let cram_a = open_cram_from_records("a.cram", recs);
        let cfg = AlignmentMergedReaderConfig {
            min_mapq: None,
            min_read_length: None,
            drop_qc_fail: false,
            drop_duplicate: false,
            ..Default::default()
        };
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            cfg,
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert_eq!(out.len(), 3);
        let expected = FilterCounts {
            unmapped: 1,
            secondary: 1,
            supplementary: 1,
            ..FilterCounts::default()
        };
        assert_eq!(counts, expected);
    }

    #[test]
    fn a11_min_read_length_drops_short_and_empty() {
        assert_eq!(DEFAULT_MIN_READ_LENGTH, 30);
        let long = pass_record("long", 0, 100);
        let mut short = pass_record("short", 0, 200);
        let short_len = (DEFAULT_MIN_READ_LENGTH as usize) - 5;
        short.cigar_ops = vec![CigarOp::Match(short_len as u32)];
        short.seq = default_seq(short_len);
        short.qual = default_qual(short_len);
        let mut empty = pass_record("empty", 0, 300);
        empty.cigar_ops = vec![CigarOp::Match(0)];
        empty.seq = Vec::new();
        empty.qual = Vec::new();

        let cram_a = open_cram_from_records(
            "a.cram",
            vec![record_spec(long), record_spec(short), record_spec(empty)],
        );
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none(), "{:?}", err);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].qname, b"long");
        assert_eq!(counts.too_short, 2);

        // None disables the filter.
        let long2 = pass_record("long", 0, 100);
        let mut short2 = pass_record("short", 0, 200);
        let short_len2 = (DEFAULT_MIN_READ_LENGTH as usize) - 5;
        short2.cigar_ops = vec![CigarOp::Match(short_len2 as u32)];
        short2.seq = default_seq(short_len2);
        short2.qual = default_qual(short_len2);
        let mut empty2 = pass_record("empty", 0, 300);
        empty2.cigar_ops = vec![CigarOp::Match(0)];
        empty2.seq = Vec::new();
        empty2.qual = Vec::new();
        let cram_a = open_cram_from_records(
            "a.cram",
            vec![record_spec(long2), record_spec(short2), record_spec(empty2)],
        );
        let cfg = AlignmentMergedReaderConfig {
            min_read_length: None,
            ..Default::default()
        };
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            cfg,
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert_eq!(out.len(), 3);
        assert_eq!(counts.too_short, 0);
    }

    #[test]
    fn a12_filter_precedence_is_hit_rate_ordered() {
        // Record marked duplicate AND unmapped AND MAPQ 0 → only
        // duplicate increments (highest-hit-rate filter runs first).
        let mut combo = pass_record("combo", 0, 100);
        combo.flag = FLAG_DUPLICATE | FLAG_UNMAPPED;
        combo.mapq = 0;
        // Record marked unmapped AND MAPQ 0 (not duplicate) → low_mapq.
        let mut mq_un = pass_record("mq_un", 0, 200);
        mq_un.flag = FLAG_UNMAPPED;
        mq_un.mapq = 0;
        // Record marked unmapped only → unmapped.
        let mut un_only = pass_record("un_only", 0, 300);
        un_only.flag = FLAG_UNMAPPED;
        // un_only.mapq stays at the pass-record default (60).

        let cram_a = open_cram_from_records(
            "a.cram",
            vec![record_spec(combo), record_spec(mq_un), record_spec(un_only)],
        );
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert!(out.is_empty(), "expected all dropped, got {:?}", out);
        assert_eq!(counts.duplicate, 1);
        assert_eq!(counts.low_mapq, 1);
        assert_eq!(counts.unmapped, 1);
        assert_eq!(counts.secondary, 0);
        assert_eq!(counts.supplementary, 0);
        assert_eq!(counts.qc_fail, 0);
    }

    #[test]
    fn a13_empty_stream() {
        let cram_a = open_cram_from_records("a.cram", vec![]);
        let mut reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        assert!(reader.next().is_none());
    }

    #[test]
    fn a14_mixed_empty_and_non_empty() {
        let cram_a = open_cram_from_records("a.cram", vec![]);
        let cram_b = open_cram_from_records(
            "b.cram",
            vec![
                record_spec(pass_record("R1", 0, 100)),
                record_spec(pass_record("R2", 0, 200)),
            ],
        );
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a, cram_b],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, _, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].qname, b"R1");
        assert_eq!(out[1].qname, b"R2");
    }

    #[test]
    fn a14b_truly_unmapped_record_is_filtered_not_errored() {
        // Build a record whose unmapped flag is set and whose pos is 0
        // (the record_spec `if spec.pos > 0` guard means no
        // alignment_start is set on the underlying RecordBuf — i.e. the
        // shape of a real-CRAM unmapped record).
        // `mapq` is set above `DEFAULT_MIN_MAPQ` so the low-MAPQ filter
        // (which sits earlier in the cascade) doesn't pre-empt the
        // unmapped bucket; the point of this test is the unmapped path.
        let unmapped = RecordSpec {
            qname: "U".into(),
            flag: FLAG_UNMAPPED,
            ref_id: 0,
            pos: 0,
            mapq: 60,
            cigar_ops: vec![],
            seq: vec![],
            qual: vec![],
            mate_ref_id: None,
            mate_pos: None,
        };
        let cram = open_cram_from_records("a.cram", vec![record_spec(unmapped)]);
        let mut reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        assert!(
            reader.next().is_none(),
            "unmapped record must be silently dropped"
        );
        assert_eq!(reader.filter_counts().unmapped, 1);
    }

    #[test]
    fn a15_position_above_u32_max_round_trips_through_mapped_read() {
        let big_pos: u64 = (u32::MAX as u64) + 100;
        let rec = pass_record("R", 0, big_pos);
        let cram_a = open_cram_from_records("a.cram", vec![record_spec(rec)]);
        let mut contigs = default_contigs();
        contigs.entries[0].length = big_pos + 1_000;
        let mut reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram_a],
            contigs,
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let read = reader.next().expect("ok").expect("ok");
        assert_eq!(read.pos, big_pos);
        assert!(reader.next().is_none());
    }

    #[test]
    fn g2_consecutive_id_record_is_dropped_with_counter() {
        // Build a record whose CIGAR has an adjacent I/D pair. The
        // record is otherwise pristine (passes MAPQ, length, flag
        // checks) so the only filter that can fire is G2.
        let len = (DEFAULT_MIN_READ_LENGTH as usize) + 20;
        let mut spec = pass_record("R", 0, 100);
        // Adjust seq/qual length to match the new CIGAR's read-consuming ops.
        // M(20) + I(2) + D(1) + M(rest) → read consumes 20 + 2 + (len - 20) = len + 2.
        let m_tail = (len - 20) as u32;
        spec.cigar_ops = vec![
            CigarOp::Match(20),
            CigarOp::Insertion(2),
            CigarOp::Deletion(1),
            CigarOp::Match(m_tail),
        ];
        spec.seq = default_seq(len + 2);
        spec.qual = default_qual(len + 2);
        let cram = open_cram_from_records("a.cram", vec![record_spec(spec)]);
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none(), "G2 drop should not surface an error");
        assert!(out.is_empty(), "G2 must drop the consecutive-I/D record");
        assert_eq!(counts.bad_cigar, 1);
    }

    #[test]
    fn g2_first_op_deletion_record_is_dropped_with_counter() {
        // Boundary deletion: D as the first non-clip op. We add a
        // leading soft-clip to verify the clip-stripping rule fires
        // on the *post-clip* first op.
        let len = (DEFAULT_MIN_READ_LENGTH as usize) + 20;
        let mut spec = pass_record("R", 0, 100);
        let m_tail = (len - 5) as u32;
        spec.cigar_ops = vec![
            CigarOp::SoftClip(3),
            CigarOp::Deletion(2),
            CigarOp::Match(m_tail),
        ];
        // Read consumes 3 (SoftClip) + (len - 5) (Match) = len - 2.
        spec.seq = default_seq(len - 2);
        spec.qual = default_qual(len - 2);
        let cram = open_cram_from_records("a.cram", vec![record_spec(spec)]);
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert!(out.is_empty());
        assert_eq!(counts.bad_cigar, 1);
    }

    #[test]
    fn g2_well_formed_record_passes_through() {
        // Negative control: a clean CIGAR (single interior insertion)
        // must NOT be filtered by G2. Also pins that we are not
        // accidentally rejecting non-G2 patterns when wiring the
        // cascade.
        let len = (DEFAULT_MIN_READ_LENGTH as usize) + 20;
        let mut spec = pass_record("R", 0, 100);
        let m_tail = (len - 20) as u32;
        spec.cigar_ops = vec![
            CigarOp::Match(20),
            CigarOp::Insertion(2),
            CigarOp::Match(m_tail),
        ];
        spec.seq = default_seq(len + 2);
        spec.qual = default_qual(len + 2);
        let cram = open_cram_from_records("a.cram", vec![record_spec(spec)]);
        let reader = AlignmentMergedReader::from_open_alignment_files(
            vec![cram],
            default_contigs(),
            "sample".into(),
            AlignmentMergedReaderConfig::default(),
            None,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert_eq!(out.len(), 1);
        assert_eq!(counts.bad_cigar, 0);
    }

    #[test]
    fn a16_empty_inputs_returns_no_inputs() {
        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("ref.fa");
        // The .fai is checked after the empty-input guard; touching a
        // file isn't required because we expect to fail before the
        // FASTA is even opened.
        let result =
            AlignmentMergedReader::new(&[], &fasta_path, AlignmentMergedReaderConfig::default());
        match result {
            Ok(_) => panic!("expected NoInputs, got Ok"),
            Err(AlignmentInputError::NoInputs) => {}
            Err(other) => panic!("expected NoInputs, got {:?}", other),
        }
    }

    // --- F1 mismatch-fraction helper: pure-logic tests ---------------
    //
    // These exercise `read_exceeds_mismatch_fraction` directly, with
    // synthetic inputs (no Repository, no CRAM). Integration of the
    // filter into the alignment_input pipeline is covered by the
    // `default_config_*` tests below.

    fn m_only(len: u32) -> Vec<CigarOp> {
        vec![CigarOp::Match(len)]
    }

    #[test]
    fn f1_zero_mismatches_passes() {
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b"ACGTACGTAC";
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_high_mismatch_drops() {
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b"AGGGAGGTAG"; // 4 mismatches in 10 → 40 %
        assert!(read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_threshold_boundary_is_strict_greater_than() {
        // Exactly at threshold: should NOT drop (filter is `>`, not `>=`).
        // 1 mismatch in 10 = 0.10 fraction; threshold = 0.10 → keep.
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b"AGGTACGTAC"; // 1 mismatch at pos 1
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
        // Just above threshold: 2/10 = 0.20.
        let ref_seq_2mm = b"AGGAACGTAC"; // 2 mismatches
        assert!(read_exceeds_mismatch_fraction(
            &cigar,
            seq,
            &qual,
            ref_seq_2mm,
            10,
            0.10
        ));
    }

    #[test]
    fn f1_bq_floor_excludes_low_quality_mismatches() {
        // 5 mismatches in 10 bases (= 50 %), but all with BQ=5 (< floor=10).
        // Filter should ignore them all → keep.
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![5u8; 10];
        let ref_seq = b"AGGGAGGGAC"; // 5 mismatches
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
        // Same input with BQ=20 → mismatches count → drop.
        let qual_high = vec![20u8; 10];
        assert!(read_exceeds_mismatch_fraction(
            &cigar, seq, &qual_high, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_bq_floor_zero_counts_every_mismatch() {
        // floor=0 → every mismatch counts regardless of BQ.
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![0u8; 10];
        let ref_seq = b"AGGGAGGGAC"; // 5 mismatches
        assert!(read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 0, 0.10
        ));
    }

    #[test]
    fn f1_n_in_ref_skipped_from_both_numerator_and_denominator() {
        // Read = ACGTACGTAC, Ref = NNNNNNNNNN. All M positions hit
        // N-in-ref → comparable_bases = 0 → keep (no signal).
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b"NNNNNNNNNN";
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_n_in_read_skipped_from_both_numerator_and_denominator() {
        // 5 read bases are N (skipped); the other 5 all match
        // → 0 mismatches in 5 comparable → keep.
        let cigar = m_only(10);
        let seq = b"NNNNNACGTA";
        let qual = vec![30u8; 10];
        let ref_seq = b"AAAAAACGTA";
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_soft_clips_excluded_from_denominator() {
        // 5S5M: only the 5 M-op bases count. Two of them mismatch
        // → 2/5 = 40 % → drop. The freebayes formula (with seqlen
        // denominator = 10) would give 2/10 = 20 % — still drop here
        // but this test pins the M-op-only behaviour we documented.
        let cigar = vec![CigarOp::SoftClip(5), CigarOp::Match(5)];
        let seq = b"AAAAAACGTA"; // first 5 = soft clip
        let qual = vec![30u8; 10];
        let ref_seq = b"GGGTA"; // ref covers only the 5 M bases
        assert!(read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_indels_advance_cursors_but_dont_count() {
        // 4M 2I 4M with all matches in M-ops and an insertion of
        // garbage in the middle. M-op bases = 8, mismatches = 0.
        // Insertion does not contribute to either count → keep.
        let cigar = vec![CigarOp::Match(4), CigarOp::Insertion(2), CigarOp::Match(4)];
        let seq = b"ACGTNNACGT"; // last 4 align after the 2-base ins
        let qual = vec![30u8; 10];
        let ref_seq = b"ACGTACGT";
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_no_m_op_returns_false() {
        // All-soft-clip read (shouldn't happen post-MAPQ-filter but
        // could). comparable_bases = 0 → keep (conservative).
        let cigar = vec![CigarOp::SoftClip(10)];
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b""; // ref_span = 0
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_default_config_has_filter_enabled_at_ten_percent() {
        // Pin the default config: filter on, threshold 0.10, floor 10.
        let cfg = AlignmentMergedReaderConfig::default();
        assert_eq!(cfg.max_read_mismatch_fraction, Some(0.10));
        assert_eq!(cfg.mismatch_bq_floor, 10);
        assert_eq!(DEFAULT_MAX_READ_MISMATCH_FRACTION, 0.10);
        assert_eq!(DEFAULT_MISMATCH_BQ_FLOOR, 10);
    }

    #[test]
    fn f1_filter_counts_default_high_mismatch_fraction_is_zero() {
        // The new field on FilterCounts must default to 0 — otherwise
        // previously-default-comparing tests would silently break.
        let counts = FilterCounts::default();
        assert_eq!(counts.high_mismatch_fraction, 0);
    }

    // --- G1 adaptor-boundary helper: pure-logic tests ----------------
    //
    // These exercise `compute_adaptor_boundary` directly, with synthetic
    // SAM-flag/TLEN/mate-pos inputs. Integration into MappedRead
    // construction is covered by the Group B tests via the `flag` /
    // `mate_pos` knobs already available on `RecordSpec`.

    /// Build the standard "properly paired, FR orientation" flag for a
    /// forward-strand read whose mate is reverse-strand. Tests pass
    /// this directly so the SAM bit fiddling is centralised.
    const FLAG_FR_FWD: u16 = FLAG_PAIRED | FLAG_MATE_REVERSE_STRAND;
    /// Same as `FLAG_FR_FWD` but for the reverse-strand mate.
    const FLAG_FR_REV: u16 = FLAG_PAIRED | FLAG_REVERSE_STRAND;

    #[test]
    fn g1_single_end_returns_none() {
        // No FLAG_PAIRED → boundary is undefined.
        assert_eq!(
            compute_adaptor_boundary(0, 0, 100, 150, Some(0), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_mate_unmapped_returns_none() {
        let flag = FLAG_FR_FWD | FLAG_MATE_UNMAPPED;
        assert_eq!(
            compute_adaptor_boundary(flag, 0, 100, 150, Some(0), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_self_unmapped_returns_none() {
        let flag = FLAG_FR_FWD | FLAG_UNMAPPED;
        assert_eq!(
            compute_adaptor_boundary(flag, 0, 100, 150, Some(0), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_zero_tlen_returns_none() {
        // TLEN=0 means the aligner could not determine the insert size,
        // typically because the mates are on different references or
        // the placement is ambiguous. Don't filter.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(200), 0),
            None
        );
    }

    #[test]
    fn g1_same_strand_pair_returns_none() {
        // Same-strand FF or RR: improperly oriented, geometry is not
        // the simple readthrough case.
        let same_strand_flag = FLAG_PAIRED; // both forward (no FLAG_MATE_REVERSE_STRAND)
        assert_eq!(
            compute_adaptor_boundary(same_strand_flag, 0, 100, 150, Some(0), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_mate_on_different_contig_returns_none() {
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(1), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_molecule_at_least_as_long_as_read_returns_none() {
        // 150bp read in a 200bp insert: fragment ≥ read, no readthrough.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(250), 200),
            None
        );
        // Edge case: |tlen| == seq_len, also no readthrough.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(250), 150),
            None
        );
    }

    #[test]
    fn g1_forward_strand_short_insert_boundary_is_start_plus_tlen() {
        // 150bp read, insert 50bp: read overruns by 100bp into the
        // 3′ adaptor. Boundary at read.start + |tlen| = 100 + 50 = 150.
        // Any base at ref_pos >= 150 is in adaptor.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(149), 50),
            Some(150)
        );
    }

    #[test]
    fn g1_reverse_strand_short_insert_boundary_is_mate_start_minus_one() {
        // The reverse-strand read in an FR pair: mate's start is the
        // forward-strand mate's leftmost ref pos (the molecule's 5′
        // end on the forward strand). Boundary = mate.start - 1.
        // ancient-DNA shape: 50bp molecule, mate starts at 100, our
        // reverse read end placed beyond.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_REV, 0, 200, 150, Some(0), Some(100), -50),
            Some(99)
        );
    }

    #[test]
    fn g1_negative_tlen_takes_absolute_value() {
        // Forward-strand mate that happens to be the second-in-pair
        // (TLEN convention is reverse-mate TLEN, forward-mate -TLEN).
        // |tlen| should be used so a negative sign does not break
        // the calculation.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(149), -50),
            Some(150)
        );
    }

    #[test]
    fn g1_ancient_dna_short_fragment_yields_boundary_inside_read_span() {
        // 70bp ancient-DNA molecule sequenced with 100bp reads:
        // mate starts at 1000, |tlen| = 70.
        // Forward mate (start 1000): boundary at 1000 + 70 = 1070.
        // Read covers ref [1000, 1100), so positions 1070..1100
        // are filtered (30 of 100 bases).
        let boundary_fwd =
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 1000, 100, Some(0), Some(1030), 70);
        assert_eq!(boundary_fwd, Some(1070));
        // Reverse mate (the same molecule's other end): boundary
        // at mate.start - 1 = 1000 - 1 = 999. The reverse mate
        // sits at ref [1030, 1130), but its leftward overrun is
        // anywhere ≤ 999 — outside its own placement, which means
        // its leftmost soft-clipped/`M` bases past 999 are flagged.
        let boundary_rev =
            compute_adaptor_boundary(FLAG_FR_REV, 0, 1030, 100, Some(0), Some(1000), -70);
        assert_eq!(boundary_rev, Some(999));
    }

    #[test]
    fn g1_default_record_has_no_boundary() {
        // A record without any of the paired flags / TLEN must yield
        // None — the no-op case the cursor relies on.
        assert_eq!(
            compute_adaptor_boundary(0, 0, 100, 150, None, None, 0),
            None
        );
    }

    // --- G2 cigar_is_bad helper: pure-logic tests --------------------
    //
    // These exercise `cigar_is_bad` directly. Integration with the
    // alignment_input filter cascade and the `FilterCounts.bad_cigar`
    // counter is covered by the Group B tests below.

    #[test]
    fn g2_well_formed_cigar_is_not_bad() {
        // Standard shapes that should pass: pure-match, single
        // interior indel, soft-clipped ends, multi-indel CIGARs
        // separated by M ops.
        assert!(!cigar_is_bad(&[CigarOp::Match(100)]));
        assert!(!cigar_is_bad(&[
            CigarOp::Match(50),
            CigarOp::Insertion(2),
            CigarOp::Match(50),
        ]));
        assert!(!cigar_is_bad(&[
            CigarOp::Match(50),
            CigarOp::Deletion(3),
            CigarOp::Match(50),
        ]));
        assert!(!cigar_is_bad(&[
            CigarOp::SoftClip(5),
            CigarOp::Match(90),
            CigarOp::SoftClip(5),
        ]));
        assert!(!cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Insertion(2),
            CigarOp::Match(40),
            CigarOp::Deletion(1),
            CigarOp::Match(20),
        ]));
    }

    #[test]
    fn g2_consecutive_id_is_bad() {
        // I followed by D.
        assert!(cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Insertion(2),
            CigarOp::Deletion(1),
            CigarOp::Match(20),
        ]));
        // D followed by I (the symmetric case).
        assert!(cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Deletion(1),
            CigarOp::Insertion(2),
            CigarOp::Match(20),
        ]));
    }

    #[test]
    fn g2_consecutive_same_kind_indels_are_not_bad() {
        // II or DD is non-canonical but not the GoodCigar pattern —
        // those should be merged or normalised by the aligner; we
        // don't reject on them here because the GATK rule is only
        // about *mixed* I/D pairs.
        assert!(!cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Insertion(2),
            CigarOp::Insertion(1),
            CigarOp::Match(20),
        ]));
        assert!(!cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Deletion(2),
            CigarOp::Deletion(1),
            CigarOp::Match(20),
        ]));
    }

    #[test]
    fn g2_first_op_deletion_is_bad() {
        assert!(cigar_is_bad(&[CigarOp::Deletion(2), CigarOp::Match(50),]));
    }

    #[test]
    fn g2_last_op_deletion_is_bad() {
        assert!(cigar_is_bad(&[CigarOp::Match(50), CigarOp::Deletion(2),]));
    }

    #[test]
    fn g2_first_op_deletion_after_clip_is_bad() {
        // GATK's rule explicitly says "with or without preceding
        // clips": a deletion after only soft/hard clips is still a
        // boundary deletion.
        assert!(cigar_is_bad(&[
            CigarOp::SoftClip(5),
            CigarOp::Deletion(2),
            CigarOp::Match(50),
        ]));
        assert!(cigar_is_bad(&[
            CigarOp::HardClip(5),
            CigarOp::SoftClip(3),
            CigarOp::Deletion(2),
            CigarOp::Match(50),
        ]));
    }

    #[test]
    fn g2_last_op_deletion_after_trailing_clip_is_bad() {
        assert!(cigar_is_bad(&[
            CigarOp::Match(50),
            CigarOp::Deletion(2),
            CigarOp::SoftClip(5),
        ]));
        assert!(cigar_is_bad(&[
            CigarOp::Match(50),
            CigarOp::Deletion(2),
            CigarOp::SoftClip(3),
            CigarOp::HardClip(2),
        ]));
    }

    #[test]
    fn g2_first_op_insertion_is_not_bad() {
        // Boundary insertions are not the GoodCigar rule. The walker
        // already drops first/last-op insertions as events (no
        // flanking evidence), but the read itself stays in the active
        // set and contributes Match events from the surrounding ops.
        assert!(!cigar_is_bad(&[CigarOp::Insertion(2), CigarOp::Match(50),]));
        assert!(!cigar_is_bad(&[CigarOp::Match(50), CigarOp::Insertion(2),]));
        assert!(!cigar_is_bad(&[
            CigarOp::SoftClip(3),
            CigarOp::Insertion(2),
            CigarOp::Match(50),
        ]));
    }

    #[test]
    fn g2_only_clips_is_not_bad() {
        // Pathological but not the GoodCigar pattern: a CIGAR of just
        // clips has no first/last *non-clip* op, so the boundary-
        // deletion check finds no candidate and returns false.
        // Reads like this are dropped earlier by the min-read-length
        // / unmapped-flag checks, so this is just a defensive case.
        assert!(!cigar_is_bad(&[CigarOp::SoftClip(50)]));
        assert!(!cigar_is_bad(&[
            CigarOp::HardClip(10),
            CigarOp::SoftClip(40),
        ]));
    }

    #[test]
    fn g2_empty_cigar_is_not_bad() {
        // An empty CIGAR has no ops at all — `windows(2)` yields
        // nothing, `find` returns None, neither rule fires. The
        // upstream record decode would normally reject empty CIGARs
        // before this point; this test pins the no-panic behaviour.
        assert!(!cigar_is_bad(&[]));
    }

    // --- Group B: via new (real CRAM + FASTA) ------------------------

    use crate::pileup::per_sample::cram_files::{
        ContigSpec, HeaderOverrides, build_cram, build_cram_with_major_version, build_fasta,
    };

    /// `expect_err` requires `Debug` on the success type;
    /// `AlignmentMergedReader` carries a boxed `Send` iterator and cannot be
    /// `Debug`. This helper unwraps the error variant and panics with
    /// the supplied message if the result was `Ok`.
    fn err_or_panic(
        result: Result<AlignmentMergedReader, AlignmentInputError>,
        message: &str,
    ) -> AlignmentInputError {
        match result {
            Ok(_) => panic!("{message}"),
            Err(e) => e,
        }
    }

    fn one_contig_chr1() -> Vec<ContigSpec> {
        vec![ContigSpec {
            name: "chr1".into(),
            length: 1_000,
        }]
    }

    fn pass_record_for_b(qname: &str, ref_id: usize, pos: u64) -> RecordBufForB {
        let len = (DEFAULT_MIN_READ_LENGTH as usize) + 20;
        record_spec(RecordSpec {
            qname: qname.into(),
            flag: 0,
            ref_id,
            pos,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(len as u32)],
            seq: default_seq(len),
            qual: default_qual(len),
            mate_ref_id: None,
            mate_pos: None,
        })
    }

    type RecordBufForB = noodles_sam::alignment::record_buf::RecordBuf;

    #[test]
    fn b1_header_parsing_on_known_good_cram() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("build_fasta");
        let records = vec![
            pass_record_for_b("R1", 0, 100),
            pass_record_for_b("R2", 0, 200),
        ];
        let (_cram_dir, input_path) = build_cram(
            &fasta_path,
            &contigs,
            &HeaderOverrides {
                read_groups: vec![("rg0".into(), Some("s1".into()))],
                ..Default::default()
            },
            &records,
        )
        .expect("build_cram");
        let reader = AlignmentMergedReader::new(
            &[input_path],
            &fasta_path,
            AlignmentMergedReaderConfig::default(),
        )
        .expect("AlignmentMergedReader::new");
        assert_eq!(reader.sample_name(), "s1");
        let contig_list = reader.contigs();
        assert_eq!(contig_list.entries.len(), 1);
        assert_eq!(contig_list.entries[0].name, "chr1");
        assert_eq!(contig_list.entries[0].length, 1_000);
    }

    #[test]
    fn b2_cram_4x_rejected_at_header_time() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("build_fasta");
        let (_cram_dir, input_path) = build_cram_with_major_version(4, 0).expect("build forced");
        let err = err_or_panic(
            AlignmentMergedReader::new(
                std::slice::from_ref(&input_path),
                &fasta_path,
                AlignmentMergedReaderConfig::default(),
            ),
            "expected UnsupportedCramVersion",
        );
        match err {
            AlignmentInputError::UnsupportedCramVersion { major, minor, path } => {
                assert_eq!(major, 4);
                assert_eq!(minor, 0);
                assert_eq!(path, input_path);
            }
            other => panic!("unexpected: {:?}", other),
        }
    }

    #[test]
    fn b3_non_coordinate_sort_rejected() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("build_fasta");
        let records = vec![pass_record_for_b("R1", 0, 100)];
        let (_cram_dir, input_path) = build_cram(
            &fasta_path,
            &contigs,
            &HeaderOverrides {
                sort_order: Some("queryname".into()),
                ..Default::default()
            },
            &records,
        )
        .expect("build_cram");
        let err = err_or_panic(
            AlignmentMergedReader::new(
                &[input_path],
                &fasta_path,
                AlignmentMergedReaderConfig::default(),
            ),
            "expected NotCoordinateSorted",
        );
        match err {
            AlignmentInputError::NotCoordinateSorted { sort_order, .. } => {
                assert_eq!(sort_order, "queryname");
            }
            other => panic!("unexpected: {:?}", other),
        }
    }

    #[test]
    fn b4_sample_tag_handling() {
        let contigs = one_contig_chr1();

        // 1) Single CRAM with one @RG SM:foo
        let (_fa1_dir, fa1) = build_fasta(&contigs).expect("fasta");
        let records = vec![pass_record_for_b("R1", 0, 100)];
        let (_c1_dir, c1) = build_cram(
            &fa1,
            &contigs,
            &HeaderOverrides {
                read_groups: vec![("rg0".into(), Some("foo".into()))],
                ..Default::default()
            },
            &records,
        )
        .expect("build");
        let reader =
            AlignmentMergedReader::new(&[c1], &fa1, AlignmentMergedReaderConfig::default())
                .expect("reader");
        assert_eq!(reader.sample_name(), "foo");

        // 2) Single CRAM with @RG that has no SM tag → MissingSampleTag
        let (_fa2_dir, fa2) = build_fasta(&contigs).expect("fasta");
        let (_c2_dir, c2) = build_cram(
            &fa2,
            &contigs,
            &HeaderOverrides {
                read_groups: vec![("rg0".into(), None)],
                ..Default::default()
            },
            &records,
        )
        .expect("build");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[c2], &fa2, AlignmentMergedReaderConfig::default()),
            "expected MissingSampleTag",
        );
        assert!(
            matches!(err, AlignmentInputError::MissingSampleTag { .. }),
            "got {:?}",
            err
        );

        // 3) Two CRAMs with disagreeing SMs → MultipleSampleNames.
        let (_fa3_dir, fa3) = build_fasta(&contigs).expect("fasta");
        let (_c3a_dir, c3a) = build_cram(
            &fa3,
            &contigs,
            &HeaderOverrides {
                read_groups: vec![("rg0".into(), Some("foo".into()))],
                ..Default::default()
            },
            &records,
        )
        .expect("build foo");
        let (_c3b_dir, c3b) = build_cram(
            &fa3,
            &contigs,
            &HeaderOverrides {
                read_groups: vec![("rg0".into(), Some("bar".into()))],
                ..Default::default()
            },
            &records,
        )
        .expect("build bar");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[c3a, c3b], &fa3, AlignmentMergedReaderConfig::default()),
            "expected MultipleSampleNames",
        );
        assert!(
            matches!(err, AlignmentInputError::MultipleSampleNames { .. }),
            "got {:?}",
            err
        );

        // 4) Single CRAM with two @RG entries that disagree on SM →
        //    MultipleSampleNamesInFile, NOT MultipleSampleNames.
        let (_fa4_dir, fa4) = build_fasta(&contigs).expect("fasta");
        let (_c4_dir, c4) = build_cram(
            &fa4,
            &contigs,
            &HeaderOverrides {
                read_groups: vec![
                    ("rg_foo".into(), Some("foo".into())),
                    ("rg_bar".into(), Some("bar".into())),
                ],
                ..Default::default()
            },
            &records,
        )
        .expect("build within-file");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[c4], &fa4, AlignmentMergedReaderConfig::default()),
            "expected MultipleSampleNamesInFile",
        );
        match err {
            AlignmentInputError::MultipleSampleNamesInFile {
                rg_a,
                sm_a,
                rg_b,
                sm_b,
                ..
            } => {
                assert_eq!(rg_a, "rg_foo");
                assert_eq!(sm_a, "foo");
                assert_eq!(rg_b, "rg_bar");
                assert_eq!(sm_b, "bar");
            }
            other => panic!("expected MultipleSampleNamesInFile, got {:?}", other),
        }
    }

    #[test]
    fn b4b_malformed_md5_rejected() {
        let contigs = one_contig_chr1();
        let records = vec![pass_record_for_b("R1", 0, 100)];

        // 1) M5 with the wrong length (7 chars).
        let (_fa1_dir, fa1) = build_fasta(&contigs).expect("fasta");
        let (_c1_dir, c1) = build_cram(
            &fa1,
            &contigs,
            &HeaderOverrides {
                md5_overrides: vec![("chr1".into(), "not_hex".into())],
                ..Default::default()
            },
            &records,
        )
        .expect("build short md5");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[c1], &fa1, AlignmentMergedReaderConfig::default()),
            "expected MalformedMd5",
        );
        match err {
            AlignmentInputError::MalformedMd5 { contig, detail, .. } => {
                assert_eq!(contig, "chr1");
                assert!(
                    detail.contains("expected 32 hex chars"),
                    "detail = {}",
                    detail
                );
            }
            other => panic!("expected MalformedMd5, got {:?}", other),
        }

        // 2) M5 with 30 hex chars (still wrong length, just barely).
        let (_fa2_dir, fa2) = build_fasta(&contigs).expect("fasta");
        let (_c2_dir, c2) = build_cram(
            &fa2,
            &contigs,
            &HeaderOverrides {
                md5_overrides: vec![("chr1".into(), "0".repeat(30))],
                ..Default::default()
            },
            &records,
        )
        .expect("build short md5 30");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[c2], &fa2, AlignmentMergedReaderConfig::default()),
            "expected MalformedMd5 short",
        );
        assert!(
            matches!(err, AlignmentInputError::MalformedMd5 { .. }),
            "got {:?}",
            err
        );

        // 3) M5 with 32 chars but a non-hex character.
        let (_fa3_dir, fa3) = build_fasta(&contigs).expect("fasta");
        let mut bad = "0".repeat(31);
        bad.push('Z');
        let (_c3_dir, c3) = build_cram(
            &fa3,
            &contigs,
            &HeaderOverrides {
                md5_overrides: vec![("chr1".into(), bad)],
                ..Default::default()
            },
            &records,
        )
        .expect("build non-hex md5");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[c3], &fa3, AlignmentMergedReaderConfig::default()),
            "expected MalformedMd5 non-hex",
        );
        match err {
            AlignmentInputError::MalformedMd5 { detail, .. } => {
                assert!(detail.contains("non-hex"), "detail = {}", detail);
            }
            other => panic!("expected MalformedMd5, got {:?}", other),
        }
    }

    #[test]
    fn b5_contig_list_mismatch_across_crams() {
        let contigs = one_contig_chr1();
        let (_fa_dir, fa) = build_fasta(&contigs).expect("fasta");

        // 1) Name/order disagreement: build a two-contig FASTA, then
        //    write two CRAMs whose @SQ entries are in opposite order.
        //    Both CRAMs match the FASTA (just with different orders),
        //    so the error fires from the cross-CRAM check rather than
        //    the FASTA agreement check. Order matters for coordinate
        //    sort, so disagreement is a hard error per the spec.
        let two_contigs = vec![
            ContigSpec {
                name: "chr1".into(),
                length: 1_000,
            },
            ContigSpec {
                name: "chr2".into(),
                length: 1_000,
            },
        ];
        let (_fa_two_dir, fa_two) = build_fasta(&two_contigs).expect("fasta two");
        let two_swapped = vec![
            ContigSpec {
                name: "chr2".into(),
                length: 1_000,
            },
            ContigSpec {
                name: "chr1".into(),
                length: 1_000,
            },
        ];
        let (_a1_dir, a1) = build_cram(&fa_two, &two_contigs, &HeaderOverrides::default(), &[])
            .expect("normal cram");
        let (_a2_dir, a2) = build_cram(&fa_two, &two_swapped, &HeaderOverrides::default(), &[])
            .expect("swapped cram");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[a1, a2], &fa_two, AlignmentMergedReaderConfig::default()),
            "expected ContigListMismatch (name)",
        );
        match err {
            AlignmentInputError::ContigListMismatch { detail, .. } => {
                assert!(detail.contains("name disagreement"), "detail = {}", detail);
            }
            other => panic!("unexpected: {:?}", other),
        }

        // 2) Length disagreement: a CRAM whose @SQ length is 2_000
        //    while the reference and the other CRAM say 1_000. We
        //    have to build the FASTA at the longer length so that
        //    the CRAM containing the long contig can be written, then
        //    write the short-length CRAM against the same FASTA but
        //    with a length override that disagrees with both the FASTA
        //    and the long CRAM. To make the test focus on the
        //    cross-CRAM check, we build the FASTA matching the longer
        //    CRAM, write the long CRAM, then write a second CRAM with
        //    a length override that disagrees with both.
        let long_contigs = vec![ContigSpec {
            name: "chr1".into(),
            length: 2_000,
        }];
        let (_fa2_dir, fa2) = build_fasta(&long_contigs).expect("fasta long");
        let (_b1_dir, b1) =
            build_cram(&fa2, &long_contigs, &HeaderOverrides::default(), &[]).expect("long cram");
        let (_b2_dir, b2) = build_cram(
            &fa2,
            &long_contigs,
            &HeaderOverrides {
                length_overrides: vec![("chr1".into(), 1_000)],
                ..Default::default()
            },
            &[],
        )
        .expect("short-length cram");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[b1, b2], &fa2, AlignmentMergedReaderConfig::default()),
            "expected ContigListMismatch (length)",
        );
        match err {
            AlignmentInputError::ContigListMismatch { detail, .. } => {
                assert!(
                    detail.contains("length disagreement"),
                    "detail = {}",
                    detail
                );
            }
            AlignmentInputError::FastaContigMismatch { detail, .. } => {
                // Acceptable alternate route: the second CRAM's @SQ
                // length disagrees with the FASTA before the cross-
                // CRAM check fires.
                assert!(
                    detail.contains("length disagreement"),
                    "detail = {}",
                    detail
                );
            }
            other => panic!("unexpected: {:?}", other),
        }

        // 3) Both have MD5 but they differ.
        let (_c1_dir, c1) = build_cram(
            &fa,
            &contigs,
            &HeaderOverrides {
                md5_overrides: vec![("chr1".into(), "0".repeat(32))],
                ..Default::default()
            },
            &[],
        )
        .expect("md5 0 cram");
        let (_c2_dir, c2) = build_cram(
            &fa,
            &contigs,
            &HeaderOverrides {
                md5_overrides: vec![("chr1".into(), "f".repeat(32))],
                ..Default::default()
            },
            &[],
        )
        .expect("md5 f cram");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[c1, c2], &fa, AlignmentMergedReaderConfig::default()),
            "expected ContigListMismatch (md5)",
        );
        match err {
            AlignmentInputError::ContigListMismatch { detail, .. } => {
                assert!(detail.contains("md5 disagreement"), "detail = {}", detail);
            }
            other => panic!("unexpected: {:?}", other),
        }
    }

    #[test]
    fn b6_fasta_agreement() {
        let contigs = one_contig_chr1();
        let records = vec![pass_record_for_b("R1", 0, 100)];

        // Happy path: matching FASTA.
        let (_fa_dir, fa) = build_fasta(&contigs).expect("fasta");
        let (_cram_dir, input_path) =
            build_cram(&fa, &contigs, &HeaderOverrides::default(), &records).expect("cram");
        AlignmentMergedReader::new(
            std::slice::from_ref(&input_path),
            &fa,
            AlignmentMergedReaderConfig::default(),
        )
        .expect("reader");

        // Missing .fai: build a real FASTA, then delete the .fai.
        let (_fa2_dir, fa2) = build_fasta(&contigs).expect("fasta2");
        let mut fai_path = fa2.as_os_str().to_owned();
        fai_path.push(".fai");
        std::fs::remove_file(PathBuf::from(fai_path)).expect("rm fai");
        // We can't open a CRAM without the FASTA repository being
        // indexable, but the missing-.fai check happens before any
        // CRAM is opened, so reuse the already-built CRAM path.
        let err = err_or_panic(
            AlignmentMergedReader::new(&[input_path], &fa2, AlignmentMergedReaderConfig::default()),
            "expected MissingFastaIndex",
        );
        assert!(
            matches!(err, AlignmentInputError::MissingFastaIndex { .. }),
            "got {:?}",
            err
        );

        // Length mismatch: FASTA at length 1_000, CRAM @SQ overridden
        // to 1_000 but a different FASTA at 500 — write a fresh FASTA
        // with a different length and reuse the original CRAM.
        let other_contigs = vec![ContigSpec {
            name: "chr1".into(),
            length: 500,
        }];
        let (_fa3_dir, fa3) = build_fasta(&other_contigs).expect("fasta3");
        let (_cram_b_dir, cram_b) =
            build_cram(&fa, &contigs, &HeaderOverrides::default(), &records).expect("cram b");
        let err = err_or_panic(
            AlignmentMergedReader::new(&[cram_b], &fa3, AlignmentMergedReaderConfig::default()),
            "expected FastaContigMismatch",
        );
        assert!(
            matches!(err, AlignmentInputError::FastaContigMismatch { .. }),
            "got {:?}",
            err
        );
    }

    #[test]
    fn b7_end_to_end_smoke() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");
        let records = vec![
            pass_record_for_b("R1", 0, 100),
            pass_record_for_b("R2", 0, 200),
            pass_record_for_b("R3", 0, 300),
        ];
        let (_cram_dir, input_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records)
                .expect("build_cram");
        let reader = AlignmentMergedReader::new(
            &[input_path],
            &fasta_path,
            AlignmentMergedReaderConfig::default(),
        )
        .expect("reader");
        let mut decoded: Vec<MappedRead> = Vec::new();
        for item in reader {
            decoded.push(item.expect("ok record"));
        }
        assert_eq!(decoded.len(), 3);
        assert_eq!(decoded[0].pos, 100);
        assert_eq!(decoded[1].pos, 200);
        assert_eq!(decoded[2].pos, 300);
        assert_eq!(decoded[0].qname, b"R1");
    }

    // --- Group C: AlignmentMergedReader::query (indexed, per-contig) ---------

    /// Two-contig fixture used by the query tests below.
    fn two_contigs_chr1_chr2() -> Vec<ContigSpec> {
        vec![
            ContigSpec {
                name: "chr1".into(),
                length: 1_000,
            },
            ContigSpec {
                name: "chr2".into(),
                length: 1_000,
            },
        ]
    }

    /// Build a `.crai` for `input_path` in memory and return both the
    /// pre-loaded `Arc<sam::Header>` (parsed off the same CRAM) and
    /// the `AlignmentIndex` wrapper. Mirror of what the per-chrom
    /// driver will do once at startup.
    fn load_header_and_index(input_path: &Path) -> (Arc<sam::Header>, AlignmentIndex) {
        // Header.
        let mut reader = cram::io::reader::Builder::default()
            .build_from_path(input_path)
            .expect("open cram for header");
        reader.read_file_definition().expect("read file definition");
        let header = reader.read_file_header().expect("read file header");
        // Index.
        let index = noodles_cram::fs::index(input_path).expect("build crai");
        (Arc::new(header), AlignmentIndex::Crai(Arc::new(index)))
    }

    /// Convert a `ContigSpec` list (test fixture vocabulary) into the
    /// project's `ContigList` shape that `AlignmentMergedReader::query`
    /// expects.
    fn contigs_for(specs: &[ContigSpec]) -> ContigList {
        ContigList {
            entries: specs
                .iter()
                .map(|s| ContigEntry {
                    name: s.name.clone(),
                    length: s.length,
                    md5: None,
                })
                .collect(),
        }
    }

    /// Multi-contig fixture, records placed only on the queried
    /// contig. Asserts the query yields those records with the right
    /// `ref_id` set. The complementary "records on contigs we did
    /// **not** ask for" path — i.e. noodles' multi-reference slice
    /// decode + post-decode ref_id filter inside
    /// [`OwnedIndexedCramRecords`] — is not exercised here: when a
    /// multi-contig CRAM is written with records on more than one
    /// contig, the writer can pack them into a single
    /// multi-reference slice whose decode path is sensitive to
    /// fixture details outside the in-module fixture builder's
    /// control. A realistic cross-contig fixture lands with the
    /// per-chrom driver's integration test (commit 4 of this plan).
    #[test]
    fn c1_query_in_multi_contig_fixture_yields_records_with_target_ref_id() {
        let contigs = two_contigs_chr1_chr2();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");
        let records = vec![
            pass_record_for_b("R1", 0, 100),
            pass_record_for_b("R2", 0, 800),
        ];
        let (_cram_dir, input_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records)
                .expect("build_cram");

        let (header, index) = load_header_and_index(&input_path);
        let reader = AlignmentMergedReader::query(
            std::slice::from_ref(&input_path),
            &fasta_path,
            contigs_for(&contigs),
            "s1".to_string(),
            std::slice::from_ref(&header),
            std::slice::from_ref(&index),
            "chr1",
            AlignmentMergedReaderConfig::default(),
        )
        .expect("query chr1");

        let decoded: Vec<MappedRead> = reader.map(|r| r.expect("ok record")).collect();
        assert_eq!(decoded.len(), 2);
        for r in &decoded {
            assert_eq!(r.ref_id, 0, "ref_id must be chr1 (0)");
        }
        assert_eq!(decoded[0].qname, b"R1");
        assert_eq!(decoded[1].qname, b"R2");
    }

    #[test]
    fn c2_query_preserves_coordinate_order_within_contig() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");
        let records = vec![
            pass_record_for_b("R1", 0, 100),
            pass_record_for_b("R2", 0, 300),
            pass_record_for_b("R3", 0, 500),
            pass_record_for_b("R4", 0, 700),
        ];
        let (_cram_dir, input_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records)
                .expect("build_cram");

        let (header, index) = load_header_and_index(&input_path);
        let reader = AlignmentMergedReader::query(
            std::slice::from_ref(&input_path),
            &fasta_path,
            contigs_for(&contigs),
            "s1".to_string(),
            std::slice::from_ref(&header),
            std::slice::from_ref(&index),
            "chr1",
            AlignmentMergedReaderConfig::default(),
        )
        .expect("query chr1");

        let positions: Vec<u64> = reader.map(|r| r.expect("ok record").pos).collect();
        assert_eq!(positions, vec![100, 300, 500, 700]);
    }

    /// Two CRAMs for the same sample, each with two chr1 reads at
    /// interleaving positions. The query path's k-way merge must
    /// emit them in coordinate order. We compare against the
    /// expected positions directly rather than against
    /// `AlignmentMergedReader::new` because `new` 's two-CRAM-with-records
    /// path is not currently exercised by any other test in this
    /// module and has surfaced as fragile during this commit's
    /// development; cross-validating the two readers is a follow-up
    /// once that path has its own coverage.
    #[test]
    fn c3_query_k_way_merge_yields_coordinate_sorted_records_across_inputs() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");

        let records_a = vec![
            pass_record_for_b("Ra1", 0, 100),
            pass_record_for_b("Ra2", 0, 500),
        ];
        let records_b = vec![
            pass_record_for_b("Rb1", 0, 200),
            pass_record_for_b("Rb2", 0, 600),
        ];
        let header_overrides = HeaderOverrides {
            read_groups: vec![("rg0".into(), Some("s1".into()))],
            ..Default::default()
        };
        let (_a_dir, cram_a) =
            build_cram(&fasta_path, &contigs, &header_overrides, &records_a).expect("build cram_a");
        let (_b_dir, cram_b) =
            build_cram(&fasta_path, &contigs, &header_overrides, &records_b).expect("build cram_b");

        let (header_a, index_a) = load_header_and_index(&cram_a);
        let (header_b, index_b) = load_header_and_index(&cram_b);
        let headers = vec![header_a, header_b];
        let indexes = vec![index_a, index_b];

        let query_positions: Vec<(Vec<u8>, u64)> = AlignmentMergedReader::query(
            &[cram_a, cram_b],
            &fasta_path,
            contigs_for(&contigs),
            "s1".to_string(),
            &headers,
            &indexes,
            "chr1",
            AlignmentMergedReaderConfig::default(),
        )
        .expect("query chr1")
        .map(|r| {
            let r = r.expect("ok record");
            (r.qname.clone(), r.pos)
        })
        .collect();

        assert_eq!(
            query_positions,
            vec![
                (b"Ra1".to_vec(), 100),
                (b"Rb1".to_vec(), 200),
                (b"Ra2".to_vec(), 500),
                (b"Rb2".to_vec(), 600),
            ]
        );
    }

    #[test]
    fn c4_query_empty_contig_yields_no_records() {
        let contigs = two_contigs_chr1_chr2();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");
        // All records on chr1; chr2 is empty.
        let records = vec![
            pass_record_for_b("R1", 0, 100),
            pass_record_for_b("R2", 0, 200),
        ];
        let (_cram_dir, input_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records)
                .expect("build_cram");

        let (header, index) = load_header_and_index(&input_path);
        let reader = AlignmentMergedReader::query(
            std::slice::from_ref(&input_path),
            &fasta_path,
            contigs_for(&contigs),
            "s1".to_string(),
            std::slice::from_ref(&header),
            std::slice::from_ref(&index),
            "chr2",
            AlignmentMergedReaderConfig::default(),
        )
        .expect("query chr2");

        let decoded: Vec<MappedRead> = reader.map(|r| r.expect("ok record")).collect();
        assert!(decoded.is_empty(), "chr2 has no reads in fixture");
    }

    #[test]
    fn c5_query_rejects_contig_not_in_canonical_list() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");
        let records = vec![pass_record_for_b("R1", 0, 100)];
        let (_cram_dir, input_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records)
                .expect("build_cram");

        let (header, index) = load_header_and_index(&input_path);
        let err = AlignmentMergedReader::query(
            std::slice::from_ref(&input_path),
            &fasta_path,
            contigs_for(&contigs),
            "s1".to_string(),
            std::slice::from_ref(&header),
            std::slice::from_ref(&index),
            "chrUnknown",
            AlignmentMergedReaderConfig::default(),
        );
        let err = err_or_panic(err, "query should reject unknown contig");
        assert!(matches!(
            err,
            AlignmentInputError::ContigNotInList { ref contig, .. } if contig == "chrUnknown"
        ));
    }

    #[test]
    fn c6_query_rejects_mismatched_per_input_handle_counts() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");
        let records = vec![pass_record_for_b("R1", 0, 100)];
        let (_cram_dir, input_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records)
                .expect("build_cram");

        let (header, index) = load_header_and_index(&input_path);

        // 2 crams, 1 header, 1 index → mismatch.
        let err = AlignmentMergedReader::query(
            &[input_path.clone(), input_path.clone()],
            &fasta_path,
            contigs_for(&contigs),
            "s1".to_string(),
            std::slice::from_ref(&header),
            std::slice::from_ref(&index),
            "chr1",
            AlignmentMergedReaderConfig::default(),
        );
        let err = err_or_panic(err, "query should reject mismatched handle counts");
        assert!(matches!(
            err,
            AlignmentInputError::PerInputHandleCountMismatch {
                inputs: 2,
                headers: 1,
                indexes: 1,
            }
        ));
    }

    // --- Group D: 2-CRAM-with-records diagnostics ---------------------
    //
    // Investigation surface for the streaming `AlignmentMergedReader::new`
    // failure on two CRAMs of the same sample with records on both,
    // surfaced during commit 2 of the var-calling-from-bam per-chrom
    // parallelism plan (see
    // doc/devel/reports/implementations/var_calling_from_bam_per_chromosome_2026-05-24.md
    // §"Open follow-ups" #4). Tests prefixed `d_` are diagnostic
    // probes; tests prefixed `d_regression_` will be locked in once
    // a fix or workaround ships.

    /// Helper: build two CRAMs with `records_a` and `records_b` for
    /// the same sample on chr1.
    fn build_two_crams_for_diagnostics(
        records_a: &[RecordBufForB],
        records_b: &[RecordBufForB],
    ) -> (
        tempfile::TempDir,
        PathBuf,
        tempfile::TempDir,
        tempfile::TempDir,
        PathBuf,
        PathBuf,
    ) {
        let contigs = one_contig_chr1();
        let (fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");
        let header_overrides = HeaderOverrides {
            read_groups: vec![("rg0".into(), Some("s1".into()))],
            ..Default::default()
        };
        let (a_dir, cram_a) =
            build_cram(&fasta_path, &contigs, &header_overrides, records_a).expect("build cram_a");
        let (b_dir, cram_b) =
            build_cram(&fasta_path, &contigs, &header_overrides, records_b).expect("build cram_b");
        (fasta_dir, fasta_path, a_dir, b_dir, cram_a, cram_b)
    }

    /// Regression for the 2-CRAM-with-records bug: two CRAMs, same
    /// sample, single contig, two records each at interleaving
    /// coordinates. Iteration must yield exactly four
    /// `MappedRead`s in coordinate order — no spurious trailing
    /// error from noodles' non-idempotent `read_container` at EOF.
    /// Fixed by the `eof_latched` field on [`OwnedCramRecords`].
    #[test]
    fn d1_two_crams_with_records_emits_all_records_in_coord_order() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");

        let header_overrides = HeaderOverrides {
            read_groups: vec![("rg0".into(), Some("s1".into()))],
            ..Default::default()
        };
        let (_a_dir, cram_a) = build_cram(
            &fasta_path,
            &contigs,
            &header_overrides,
            &[
                pass_record_for_b("Ra1", 0, 100),
                pass_record_for_b("Ra2", 0, 500),
            ],
        )
        .expect("build cram_a");
        let (_b_dir, cram_b) = build_cram(
            &fasta_path,
            &contigs,
            &header_overrides,
            &[
                pass_record_for_b("Rb1", 0, 200),
                pass_record_for_b("Rb2", 0, 600),
            ],
        )
        .expect("build cram_b");

        let positions: Vec<(Vec<u8>, u64)> = AlignmentMergedReader::new(
            &[cram_a, cram_b],
            &fasta_path,
            AlignmentMergedReaderConfig::default(),
        )
        .expect("merged reader open")
        .map(|r| {
            let r = r.expect("ok record");
            (r.qname.clone(), r.pos)
        })
        .collect();

        assert_eq!(
            positions,
            vec![
                (b"Ra1".to_vec(), 100),
                (b"Rb1".to_vec(), 200),
                (b"Ra2".to_vec(), 500),
                (b"Rb2".to_vec(), 600),
            ]
        );
    }

    /// Probe: open each of the two CRAMs alone via
    /// `AlignmentMergedReader::new` (passing a single-element slice), drain
    /// each to completion. Tests whether the issue is per-CRAM or
    /// specific to the merge of two.
    #[test]
    fn d2_each_cram_alone_via_new_works() {
        let (_fasta_dir, fasta_path, _a_dir, _b_dir, cram_a, cram_b) =
            build_two_crams_for_diagnostics(
                &[
                    pass_record_for_b("Ra1", 0, 100),
                    pass_record_for_b("Ra2", 0, 500),
                ],
                &[
                    pass_record_for_b("Rb1", 0, 200),
                    pass_record_for_b("Rb2", 0, 600),
                ],
            );

        let pos_a: Vec<u64> = AlignmentMergedReader::new(
            std::slice::from_ref(&cram_a),
            &fasta_path,
            AlignmentMergedReaderConfig::default(),
        )
        .expect("open a")
        .map(|r| r.expect("ok a").pos)
        .collect();
        assert_eq!(pos_a, vec![100, 500], "cram_a alone");

        let pos_b: Vec<u64> = AlignmentMergedReader::new(
            std::slice::from_ref(&cram_b),
            &fasta_path,
            AlignmentMergedReaderConfig::default(),
        )
        .expect("open b")
        .map(|r| r.expect("ok b").pos)
        .collect();
        assert_eq!(pos_b, vec![200, 600], "cram_b alone");
    }

    /// Probe: two CRAMs, but `cram_b` is empty (no records). Does the
    /// merge still trip the bug?
    #[test]
    fn d3_two_crams_one_empty_works() {
        let (_fasta_dir, fasta_path, _a_dir, _b_dir, cram_a, cram_b) =
            build_two_crams_for_diagnostics(
                &[
                    pass_record_for_b("Ra1", 0, 100),
                    pass_record_for_b("Ra2", 0, 500),
                ],
                &[],
            );

        let positions: Vec<u64> = AlignmentMergedReader::new(
            &[cram_a, cram_b],
            &fasta_path,
            AlignmentMergedReaderConfig::default(),
        )
        .expect("open")
        .map(|r| r.expect("ok").pos)
        .collect();
        assert_eq!(positions, vec![100, 500]);
    }

    /// Regression: the pre-fix bug emitted every real record
    /// correctly and then yielded a trailing `Err(MalformedRecord)`
    /// instead of `None`. This test collects all `Result`s without
    /// unwrapping and asserts exactly four `Ok` records, zero
    /// `Err` records — i.e. no trailing spurious error.
    #[test]
    fn d4_no_spurious_trailing_error_after_clean_eof() {
        let (_fasta_dir, fasta_path, _a_dir, _b_dir, cram_a, cram_b) =
            build_two_crams_for_diagnostics(
                &[
                    pass_record_for_b("Ra1", 0, 100),
                    pass_record_for_b("Ra2", 0, 500),
                ],
                &[
                    pass_record_for_b("Rb1", 0, 200),
                    pass_record_for_b("Rb2", 0, 600),
                ],
            );

        let results: Vec<Result<MappedRead, AlignmentInputError>> = AlignmentMergedReader::new(
            &[cram_a, cram_b],
            &fasta_path,
            AlignmentMergedReaderConfig::default(),
        )
        .expect("open")
        .collect();

        let (oks, errs): (Vec<_>, Vec<_>) = results.into_iter().partition(Result::is_ok);
        assert_eq!(
            oks.len(),
            4,
            "expected 4 Ok records; got {} (errs: {:?})",
            oks.len(),
            errs
        );
        assert!(
            errs.is_empty(),
            "expected zero trailing Err records; got: {errs:?}"
        );
    }

    /// Pins the upstream behaviour the `eof_latched` field on
    /// [`OwnedCramRecords`] / [`OwnedIndexedCramRecords`] is
    /// guarding against: `noodles_cram::io::Reader::read_container`
    /// is **not** idempotent at EOF. After the first `Ok(0)` (clean
    /// EOF), the next call returns an error (currently
    /// `Err(InvalidData, TryFromIntError)` on noodles-cram 0.93;
    /// the exact variant is upstream-internal and not guaranteed
    /// stable). If a noodles upgrade ever fixes this — i.e. the
    /// second call returns `Ok(0)` too — this test will fail
    /// loudly and the latch becomes redundant. Keep the test;
    /// remove the latch only when this test starts failing on the
    /// `Ok(0)` branch.
    #[test]
    fn d6_documents_noodles_read_container_non_idempotency_at_eof() {
        let (_fasta_dir, fasta_path, _a_dir, _b_dir, cram_a, _cram_b) =
            build_two_crams_for_diagnostics(
                &[
                    pass_record_for_b("Ra1", 0, 100),
                    pass_record_for_b("Ra2", 0, 500),
                ],
                &[],
            );

        let indexed_fasta_reader = noodles_fasta::io::indexed_reader::Builder::default()
            .build_from_path(&fasta_path)
            .expect("fasta reader");
        let adapter = fasta::repository::adapters::IndexedReader::new(indexed_fasta_reader);
        let repository = fasta::Repository::new(adapter);

        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(&cram_a)
            .expect("open cram_a");
        reader.read_file_definition().expect("file def");
        let _header = reader.read_file_header().expect("file header");

        // Drain by reading containers until Ok(0).
        let mut container = cram::io::reader::Container::default();
        loop {
            match reader.read_container(&mut container) {
                Ok(0) => break,
                Ok(_) => continue,
                Err(e) => panic!("unexpected error during initial drain: {e:?}"),
            }
        }

        let second = reader.read_container(&mut container);
        assert!(
            second.is_err(),
            "noodles-cram changed: second read_container after Ok(0) now \
             returned {second:?} (was Err on 0.93). If this is Ok(0), the \
             non-idempotency is fixed upstream and the `eof_latched` field \
             on OwnedCramRecords / OwnedIndexedCramRecords can be removed."
        );
    }

    /// Probe: two CRAMs with records on _different_ positions on chr1,
    /// where one CRAM's records all come BEFORE the other's. The
    /// merge therefore drains one reader fully before touching the
    /// other. If interleaving is the trigger, this should work.
    #[test]
    fn d5_two_crams_non_interleaving_positions() {
        let (_fasta_dir, fasta_path, _a_dir, _b_dir, cram_a, cram_b) =
            build_two_crams_for_diagnostics(
                &[
                    pass_record_for_b("Ra1", 0, 100),
                    pass_record_for_b("Ra2", 0, 200),
                ],
                &[
                    pass_record_for_b("Rb1", 0, 500),
                    pass_record_for_b("Rb2", 0, 600),
                ],
            );

        let positions: Vec<u64> = AlignmentMergedReader::new(
            &[cram_a, cram_b],
            &fasta_path,
            AlignmentMergedReaderConfig::default(),
        )
        .expect("open")
        .map(|r| r.expect("ok").pos)
        .collect();
        assert_eq!(positions, vec![100, 200, 500, 600]);
    }

    // --- M9: AlignmentIndexFormatMismatch ----------------------------
    //
    // The variant is the only safety net against a driver bug that
    // pairs a `.cram` path with a BAM-side index (or vice versa).
    // The catch-all branch was tightened by M6 — this test pins
    // the typed-error contract.

    #[test]
    fn query_returns_alignment_index_format_mismatch_on_cram_path_with_bam_csi_index() {
        use crate::bam::index_preflight::AlignmentIndex;
        use noodles_csi::binning_index::Indexer;
        use noodles_csi::binning_index::index::reference_sequence::index::BinnedIndex;

        // Build a real one-contig CRAM via the existing test
        // fixture (the merger needs a parseable header to reach
        // the dispatch match).
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("build_fasta");
        let (_cram_dir, cram_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &[])
                .expect("build_cram");
        let (header, _crai) = load_header_and_index(&cram_path);

        // Construct a BamCsi-typed AlignmentIndex with an empty
        // CSI Index — the match in `query` fires before any seek
        // / decode happens, so the payload is never consulted.
        let empty_csi: noodles_csi::Index = Indexer::<BinnedIndex>::new(14, 6).build(contigs.len());
        let bad_index = AlignmentIndex::BamCsi(Arc::new(empty_csi));

        let result = AlignmentMergedReader::query(
            std::slice::from_ref(&cram_path),
            &fasta_path,
            contigs_for(&contigs),
            "s1".to_string(),
            std::slice::from_ref(&header),
            std::slice::from_ref(&bad_index),
            "chr1",
            AlignmentMergedReaderConfig::default(),
        );
        let err = match result {
            Ok(_) => panic!("CRAM path + BamCsi index must error"),
            Err(e) => e,
        };

        match err {
            AlignmentInputError::AlignmentIndexFormatMismatch {
                path,
                file_format,
                index_format,
            } => {
                assert_eq!(path, cram_path);
                assert_eq!(file_format, "CRAM");
                assert_eq!(index_format, "CSI");
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }
}
