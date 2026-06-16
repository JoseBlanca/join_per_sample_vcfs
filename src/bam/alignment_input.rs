//! Alignment-file input: header validation, the per-read filter
//! cascade, and the shared `RecordBuf` → [`MappedRead`] conversion.
//! Format-agnostic at the filter seam; per-format header reading lives
//! in the sibling [`crate::bam::cram_input`] / [`crate::bam::bam_input`]
//! modules, and the pooled per-segment read fetch lives in
//! [`crate::bam::segment_reader`].
//!
//! The CRAM-specific design rationale lives in
//! `doc/devel/implementation_plans/per_sample_caller_cram_input.md`;
//! the BAM extension plan that motivated lifting the shared surface out
//! of the CRAM-named module lives in
//! `doc/devel/implementation_plans/bam_input_support.md`.

use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::bam::errors::AlignmentInputError;
use crate::bam::index_preflight::{
    AlignmentFileKind, AlignmentIndex, load_alignment_index, preflight_alignment_indexes,
};
use crate::fasta::{ContigEntry, ContigList};
use crate::pileup::walker::CigarOp;

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
    /// Index into the canonical `ContigList`.
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
    /// [`BaqSkipReason`](crate::pileup::per_sample::baq_engine::BaqSkipReason) currently does not
    /// exist — every BAQ skip reason rolls up into this counter. Split
    /// later if a deployment cares which reason dominates.
    pub baq_rejected: u64,
}

impl FilterCounts {
    /// Add `other`'s tallies into `self`, field by field. Used to total
    /// the per-region readers' filter counts into one run summary when
    /// the pileup is driven region by region.
    pub fn merge(&mut self, other: &FilterCounts) {
        self.unmapped += other.unmapped;
        self.secondary += other.secondary;
        self.supplementary += other.supplementary;
        self.qc_fail += other.qc_fail;
        self.duplicate += other.duplicate;
        self.low_mapq += other.low_mapq;
        self.too_short += other.too_short;
        self.high_mismatch_fraction += other.high_mismatch_fraction;
        self.bad_cigar += other.bad_cigar;
        self.baq_rejected += other.baq_rejected;
    }

    /// Increment the field a [`classify_pre_decode`] drop falls into.
    /// The single home for the `FilterBucket` → counter mapping, used by
    /// the segment reader that applies the pre-decode flag/MAPQ filter.
    pub(super) fn record_drop(&mut self, bucket: FilterBucket) {
        match bucket {
            FilterBucket::Unmapped => self.unmapped += 1,
            FilterBucket::Secondary => self.secondary += 1,
            FilterBucket::Supplementary => self.supplementary += 1,
            FilterBucket::QcFail => self.qc_fail += 1,
            FilterBucket::Duplicate => self.duplicate += 1,
            FilterBucket::LowMapq => self.low_mapq += 1,
            FilterBucket::BaqRejected => self.baq_rejected += 1,
        }
    }
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
///
/// Currently uncalled: the only former caller was the deleted whole-file
/// `AlignmentMergedReader::new`. Retained — and `#[allow(dead_code)]`'d
/// rather than deleted — because it is the canonical FASTA/`@SQ`
/// agreement check the indexed `load_pileup_inputs` path is expected to
/// adopt (it presently only checks that the `.fai` exists).
#[allow(dead_code)]
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

/// Build the noodles FASTA [`fasta::Repository`] used by the CRAM slice
/// decoder and the reader's F1/F3 read-reference checks.
///
/// The `Repository` is a whole-contig cache keyed by name: the first
/// fetch of any base on a contig loads (and caches) that contig's entire
/// sequence. It is `Clone` (an `Arc` bump) and `clear()`-able, so the
/// pileup region loop builds it **once** and shares it across every
/// region, clearing the cache only on contig transition to keep one
/// contig resident. Building it per region — as `query` used to —
/// reloads the whole contig for every region on it, the
/// `--regions` performance regression this centralisation fixes. See
/// `doc/devel/implementation_plans/fasta_reference_reading_unification.md`.
///
/// # Errors
///
/// - [`AlignmentInputError::MissingFastaIndex`] if the sibling `.fai`
///   does not exist.
/// - [`AlignmentInputError::Io`] if the indexed FASTA reader cannot be
///   built (unreadable FASTA / malformed `.fai`).
pub fn build_fasta_repository(fasta: &Path) -> Result<fasta::Repository, AlignmentInputError> {
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
    Ok(fasta::Repository::new(adapter))
}

// ---------------------------------------------------------------------
// ContigInterval
// ---------------------------------------------------------------------

/// A 1-based inclusive position range within a single contig, used to
/// narrow a segment read query (see
/// [`crate::bam::segment_reader::AlignmentFile::get_reads_from_segment`])
/// to a sub-contig region. `start <= end`, both `>= 1`.
///
/// The query yields every read whose reference footprint *overlaps*
/// the range — including reads that start before `start` and extend
/// into it — so a downstream pileup that clamps emitted columns to the
/// range still has the flanking reads it needs for BAQ. Passing `None`
/// to `query` (or a range that already spans the whole contig) reads
/// the entire contig with no per-record overlap check.
#[derive(Debug, Clone, Copy)]
pub struct ContigInterval {
    /// 1-based inclusive lower bound.
    pub start: u32,
    /// 1-based inclusive upper bound.
    pub end: u32,
}

impl ContigInterval {
    /// Whether `record`'s reference footprint `[alignment_start,
    /// alignment_end]` overlaps this interval. A record with no mapped
    /// position never overlaps.
    pub(crate) fn overlaps_record(&self, record: &sam::alignment::RecordBuf) -> bool {
        match (record.alignment_start(), record.alignment_end()) {
            (Some(first), Some(last)) => {
                usize::from(first) as u64 <= self.end as u64
                    && usize::from(last) as u64 >= self.start as u64
            }
            _ => false,
        }
    }
}

// ---------------------------------------------------------------------
// Indexed pileup inputs (metadata + per-input handles for segment reads)
// ---------------------------------------------------------------------

/// Per-sample alignment inputs prepared for the indexed (region-driven)
/// pileup path. Bundles the validated cross-file metadata (sample name,
/// canonical contig list) with the per-input handles the pooled segment
/// readers need (one `Arc<sam::Header>` and one [`AlignmentIndex`] per
/// input, parallel to the input-file slice) — `run_pileup` builds one
/// [`AlignmentFile`](crate::bam::segment_reader::AlignmentFile) per input
/// from them.
///
/// Built once at startup by [`load_pileup_inputs`]; the region loop
/// then fetches each region's segment reads per region, reusing these
/// handles.
#[derive(Debug)]
pub struct PileupInputs {
    /// Single sample name shared by every input (cross-validated).
    pub sample_name: String,
    /// Canonical contig list (identical across inputs), in `@SQ`
    /// order — the `chrom_id` space the [`crate::regions::RegionSet`]
    /// and the PSP writer use.
    pub contigs: ContigList,
    /// One parsed header per input, in input order.
    pub headers: Vec<Arc<sam::Header>>,
    /// One loaded alignment index per input, in input order.
    pub indexes: Vec<AlignmentIndex>,
}

/// Prepare [`PileupInputs`] for the indexed pileup path: pre-flight the
/// alignment indexes, then open each input once to validate it and
/// collect the handles `query` needs.
///
/// `build_index_if_missing` is forwarded to
/// [`preflight_alignment_indexes`]: `false` hard-errors on the first
/// input without a `.crai`/`.csi`/`.bai`; `true` builds the missing
/// index in place. The reference `.fai` is required (a missing one is
/// [`AlignmentInputError::MissingFastaIndex`]) — this loader does not
/// build it.
///
/// Cross-file invariants (identical `@SQ` lists, a single sample name,
/// coordinate sort order, per-contig `@SQ M5`) are enforced via
/// [`extract_header`]; mixed CRAM + BAM and unknown-extension inputs
/// are rejected by the pre-flight pass.
pub fn load_pileup_inputs(
    alignment_files: &[PathBuf],
    fasta: &Path,
    build_index_if_missing: bool,
) -> Result<PileupInputs, AlignmentInputError> {
    if alignment_files.is_empty() {
        return Err(AlignmentInputError::NoInputs);
    }

    // Index pre-flight first: a missing index (without the build
    // opt-in) should fail before any header parse. This pass also
    // rejects mixed CRAM + BAM and unknown extensions.
    preflight_alignment_indexes(alignment_files, build_index_if_missing)?;

    // FASTA repository — needed to open CRAM inputs. Requires the
    // sibling `.fai`.
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
    let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(
        indexed_fasta_reader,
    ));

    // One pass per input: open, read + validate the header, and
    // collect the header + index handles. Cross-file checks mirror
    // `new()` (contigs identical, single sample).
    let mut headers: Vec<Arc<sam::Header>> = Vec::with_capacity(alignment_files.len());
    let mut indexes: Vec<AlignmentIndex> = Vec::with_capacity(alignment_files.len());
    let mut canonical_contigs: Option<ContigList> = None;
    let mut canonical_sample: Option<String> = None;
    let mut reference_input_path: Option<PathBuf> = None;

    for path in alignment_files {
        let kind = AlignmentFileKind::from_path(path)
            .ok_or_else(|| AlignmentInputError::UnsupportedExtension { path: path.clone() })?;
        // The per-format header-only opener reads (and format-validates)
        // the header; we keep the header and discard the reader (no
        // records are decoded here).
        let sam_header = match kind {
            AlignmentFileKind::Cram => {
                let (_reader, header) =
                    crate::bam::cram_input::open_cram_reader_with_header(path, Some(&repository))?;
                header
            }
            AlignmentFileKind::Bam => {
                let (_reader, header) = crate::bam::bam_input::open_bam_reader_with_header(path)?;
                header
            }
        };
        let summary = extract_header(path, &sam_header)?;

        match &canonical_contigs {
            None => {
                canonical_contigs = Some(summary.contigs.clone());
                reference_input_path = Some(path.clone());
            }
            Some(existing) => {
                if let Err(detail) = existing.first_disagreement(&summary.contigs) {
                    return Err(AlignmentInputError::ContigListMismatch {
                        reference_path: reference_input_path
                            .clone()
                            .expect("reference_input_path set when canonical_contigs is Some"),
                        other_path: path.clone(),
                        detail,
                    });
                }
            }
        }
        match &canonical_sample {
            None => canonical_sample = Some(summary.sample_name.clone()),
            Some(existing) if existing != &summary.sample_name => {
                return Err(AlignmentInputError::MultipleSampleNames {
                    path_a: reference_input_path
                        .clone()
                        .expect("reference_input_path set when canonical_sample is Some"),
                    sm_a: existing.clone(),
                    path_b: path.clone(),
                    sm_b: summary.sample_name.clone(),
                });
            }
            Some(_) => {}
        }

        headers.push(Arc::new(sam_header));
        indexes.push(load_alignment_index(path)?);
    }

    Ok(PileupInputs {
        sample_name: canonical_sample.expect("at least one input validated above"),
        contigs: canonical_contigs.expect("at least one input validated above"),
        headers,
        indexes,
    })
}

pub(super) fn classify_pre_decode(
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

#[derive(Debug, Clone, Copy)]
pub(super) enum PreDecodeOutcome {
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

pub(super) fn record_buf_to_mapped_read(
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
pub(crate) fn cigar_ref_span(cigar: &[CigarOp]) -> u32 {
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
/// Rule 1 (consecutive `I`/`D`) targets the **aligner's** original
/// CIGAR — an adjacent I/D pair as emitted is the suspect signal. It
/// must run before left-alignment: the GATK port left-alignment can
/// *itself* merge two colliding indels into a canonical adjacent
/// `D`/`I` pair, which is a legitimate (parsimonious) representation,
/// not the aligner artefact this rule rejects. Running rule 1 first
/// checks only what the aligner produced.
pub(crate) fn cigar_is_bad(cigar: &[CigarOp]) -> bool {
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
pub(crate) fn read_exceeds_mismatch_fraction(
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
    use super::*;
    use crate::pileup::per_sample::record_specs::{RecordSpec, record_spec};

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

    // --- Real CRAM + FASTA fixtures (shared by Group C3) -------------

    use crate::pileup::per_sample::cram_files::{
        ContigSpec, HeaderOverrides, build_cram, build_fasta,
    };

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

    // --- Group C3: load_pileup_inputs (metadata + handles) ------------

    /// `load_pileup_inputs` with `build_if_missing = true` builds the
    /// missing `.crai` in place and returns the validated metadata plus
    /// one header + index handle per input.
    #[test]
    fn load_pileup_inputs_builds_missing_index_and_returns_metadata() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");
        let overrides = HeaderOverrides {
            read_groups: vec![("rg0".into(), Some("s1".into()))],
            ..Default::default()
        };
        let records = vec![pass_record_for_b("R1", 0, 100)];
        let (_cram_dir, input_path) =
            build_cram(&fasta_path, &contigs, &overrides, &records).expect("build_cram");

        let inputs = load_pileup_inputs(
            std::slice::from_ref(&input_path),
            &fasta_path,
            /* build = */ true,
        )
        .expect("load_pileup_inputs");

        assert_eq!(inputs.sample_name, "s1");
        assert_eq!(inputs.contigs.entries.len(), 1);
        assert_eq!(inputs.contigs.entries[0].name, "chr1");
        assert_eq!(inputs.headers.len(), 1);
        assert_eq!(inputs.indexes.len(), 1);

        // The `.crai` was created next to the input.
        let crai = {
            let mut s = input_path.clone().into_os_string();
            s.push(".crai");
            PathBuf::from(s)
        };
        assert!(
            crai.exists(),
            "build_if_missing should have written the .crai"
        );
    }

    /// With `build_if_missing = false`, a missing index hard-errors via
    /// the pre-flight pass (bridged into `AlignmentInputError`).
    #[test]
    fn load_pileup_inputs_errors_when_index_missing_and_build_disabled() {
        use crate::bam::errors::AlignmentIndexError;

        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("fasta");
        let overrides = HeaderOverrides {
            read_groups: vec![("rg0".into(), Some("s1".into()))],
            ..Default::default()
        };
        let records = vec![pass_record_for_b("R1", 0, 100)];
        let (_cram_dir, input_path) =
            build_cram(&fasta_path, &contigs, &overrides, &records).expect("build_cram");

        let err = load_pileup_inputs(
            std::slice::from_ref(&input_path),
            &fasta_path,
            /* build = */ false,
        )
        .expect_err("missing index must error when build is disabled");
        assert!(matches!(
            err,
            AlignmentInputError::AlignmentIndex(AlignmentIndexError::MissingAlignmentIndex { .. })
        ));
    }
}
