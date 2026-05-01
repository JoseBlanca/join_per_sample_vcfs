//! CRAM input slice: header validation, peek-and-scan merge, per-read
//! filter cascade.
//!
//! See `ia/feature_implementation_plans/per_sample_caller_cram_input.md`
//! for the design rationale.

use std::collections::VecDeque;
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::buffered_peekable::BufferedPeekable;
use crate::per_sample_caller::errors::CramInputError;

// ---------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------

/// Reads with MAPQ strictly below this are dropped. Matches bcftools'
/// default and the spec recommendation in
/// `ia/specs/per_sample_caller.md` §"Read filters".
pub const DEFAULT_MIN_MAPQ: u8 = 20;

/// Decoded SEQ length below this is dropped. Reads shorter than this
/// rarely contribute reliable alignments at the project's coverage
/// targets.
pub const DEFAULT_MIN_READ_LENGTH: u32 = 30;

/// `BufferedPeekable` look-ahead used per CRAM in the merge. The CRAM
/// decoder underneath already does its own slice-level batching; an
/// extra outer buffer would just delay records without saving work.
/// See `ia/specs/buffered_peekable.md` §"Why it exists in the project".
const PER_PEEKER_BUFFER_SIZE: usize = 1;

// ---------------------------------------------------------------------
// CIGAR
// ---------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match(u32),
    Insertion(u32),
    Deletion(u32),
    Skip(u32),
    SoftClip(u32),
    HardClip(u32),
    Padding(u32),
    SeqMatch(u32),
    SeqMismatch(u32),
}

// ---------------------------------------------------------------------
// MappedRead
// ---------------------------------------------------------------------

/// A read decoded out of a CRAM, owned and decoupled from noodles.
///
/// Every field is public — downstream stages (BAQ, the pileup walker)
/// need to read them directly. The struct is move-only on the hot
/// path; we avoid hidden clones by handing the decoded bytes through
/// the iterator rather than referencing them.
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
    pub source_file_index: usize,
}

// ---------------------------------------------------------------------
// ContigList
// ---------------------------------------------------------------------

/// One reference sequence (`@SQ`) entry: name, length, and optional MD5.
#[derive(Debug, Clone)]
pub struct ContigEntry {
    pub name: String,
    pub length: u64,
    pub md5: Option<[u8; 16]>,
}

impl PartialEq for ContigEntry {
    fn eq(&self, other: &Self) -> bool {
        if self.name != other.name || self.length != other.length {
            return false;
        }
        match (self.md5, other.md5) {
            (Some(a), Some(b)) => a == b,
            // Absent MD5 acts as a wildcard: a CRAM that omits M5 is
            // not contradicting a CRAM that carries one.
            (None, _) | (_, None) => true,
        }
    }
}

impl Eq for ContigEntry {}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ContigList {
    pub entries: Vec<ContigEntry>,
}

impl ContigList {
    /// First-difference report between two lists, intended for error
    /// messages. Returns `Ok(())` when the lists agree; otherwise an
    /// `Err` carrying a short string identifying which field
    /// disagreed and where.
    fn first_disagreement(&self, other: &Self) -> Result<(), String> {
        if self.entries.len() != other.entries.len() {
            return Err(format!(
                "@SQ list length differs ({} vs {})",
                self.entries.len(),
                other.entries.len()
            ));
        }
        for (i, (a, b)) in self.entries.iter().zip(other.entries.iter()).enumerate() {
            if a.name != b.name {
                return Err(format!(
                    "name disagreement at index {} ('{}' vs '{}')",
                    i, a.name, b.name
                ));
            }
            if a.length != b.length {
                return Err(format!(
                    "length disagreement at index {} (contig '{}': {} vs {})",
                    i, a.name, a.length, b.length
                ));
            }
            if let (Some(ma), Some(mb)) = (a.md5, b.md5)
                && ma != mb
            {
                return Err(format!(
                    "md5 disagreement at index {} (contig '{}')",
                    i, a.name
                ));
            }
        }
        Ok(())
    }
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
#[derive(Debug, Clone, Copy)]
pub struct CramMergedReaderConfig {
    /// `None` = no minimum; `Some(n)` = drop reads with MAPQ < n.
    pub min_mapq: Option<u8>,
    /// `None` = no minimum; `Some(n)` = drop reads with decoded SEQ
    /// length < n.
    pub min_read_length: Option<u32>,
    /// Drop reads with `flag & 0x4` set.
    pub drop_unmapped: bool,
    /// Drop reads with `flag & 0x100` set.
    pub drop_secondary: bool,
    /// Drop reads with `flag & 0x800` set.
    pub drop_supplementary: bool,
    /// Drop reads with `flag & 0x200` set.
    pub drop_qc_fail: bool,
    /// Drop reads with `flag & 0x400` set.
    pub drop_duplicate: bool,
}

// SAM/BAM flag bit constants — used both inside this module and by
// tests building synthetic records.
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_QC_FAIL: u16 = 0x200;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

impl Default for CramMergedReaderConfig {
    fn default() -> Self {
        Self {
            min_mapq: Some(DEFAULT_MIN_MAPQ),
            min_read_length: Some(DEFAULT_MIN_READ_LENGTH),
            drop_unmapped: true,
            drop_secondary: true,
            drop_supplementary: true,
            drop_qc_fail: true,
            drop_duplicate: true,
        }
    }
}

// ---------------------------------------------------------------------
// Per-CRAM lazy record stream
// ---------------------------------------------------------------------

/// One CRAM's worth of already-decoded records, plus the path string
/// used for error messages.
///
/// `pub(crate)` on purpose — `RecordBuf` is a noodles type, and
/// exposing it on a public type would re-couple our public API to
/// noodles. Tests live in the same crate, so `pub(crate)` is enough.
pub(crate) struct OpenCram {
    pub path_for_errors: PathBuf,
    pub records: Box<dyn Iterator<Item = io::Result<sam::alignment::RecordBuf>> + Send>,
}

// ---------------------------------------------------------------------
// Owned CRAM record iterator
// ---------------------------------------------------------------------

/// Owns a CRAM `Reader<File>` plus its `sam::Header` and a clone of the
/// `fasta::Repository`, and yields decoded `RecordBuf`s. Reproduces
/// the work `noodles_cram::io::reader::Records` does, but as an owned
/// iterator (`'static + Send`) so it can be packed into the merge's
/// `Box<dyn Iterator>`.
struct OwnedCramRecords {
    reader: cram::io::Reader<File>,
    header: sam::Header,
    repository: fasta::Repository,
    container: cram::io::reader::Container,
    pending: std::vec::IntoIter<sam::alignment::RecordBuf>,
}

impl OwnedCramRecords {
    fn refill(&mut self) -> io::Result<bool> {
        let n = self.reader.read_container(&mut self.container)?;
        if n == 0 {
            return Ok(true);
        }
        let compression_header = self.container.compression_header()?;
        let mut all_records: Vec<sam::alignment::RecordBuf> = Vec::new();
        for slice_result in self.container.slices() {
            let slice = slice_result?;
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            let cram_records = slice.records(
                self.repository.clone(),
                &self.header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;
            for record in &cram_records {
                let rb =
                    sam::alignment::RecordBuf::try_from_alignment_record(&self.header, record)?;
                all_records.push(rb);
            }
        }
        self.pending = all_records.into_iter();
        Ok(false)
    }
}

impl Iterator for OwnedCramRecords {
    type Item = io::Result<sam::alignment::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(rb) = self.pending.next() {
                return Some(Ok(rb));
            }
            match self.refill() {
                Ok(true) => return None,
                Ok(false) => continue,
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

// ---------------------------------------------------------------------
// Header validation
// ---------------------------------------------------------------------

/// Per-CRAM header summary, built once during pre-flight validation.
struct CramHeader {
    contigs: ContigList,
    sample_name: String,
}

/// Pull `@HD SO`, `@SQ` list, and the (single) `@RG SM` value out of a
/// `sam::Header`. Returns `CramInputError` if any required invariant
/// fails: SO != coordinate, missing SM, multiple distinct SMs in the
/// same file.
fn extract_header(path: &Path, sam_header: &sam::Header) -> Result<CramHeader, CramInputError> {
    // Sort order.
    let sort_order = sort_order_string(sam_header);
    if sort_order.as_deref() != Some("coordinate") {
        return Err(CramInputError::NotCoordinateSorted {
            path: path.to_path_buf(),
            sort_order: sort_order.unwrap_or_else(|| "<missing>".into()),
        });
    }

    // @SQ list.
    let mut entries: Vec<ContigEntry> = Vec::new();
    for (name_bstr, ref_seq_map) in sam_header.reference_sequences() {
        let name = String::from_utf8_lossy(name_bstr.as_ref()).into_owned();
        let length: u64 = usize::from(ref_seq_map.length()) as u64;
        let md5 = md5_from_reference_sequence(ref_seq_map);
        entries.push(ContigEntry { name, length, md5 });
    }
    let contigs = ContigList { entries };

    // @RG SM. All SMs in this single file must collapse to one value.
    let sample_name = extract_single_sample_name(path, sam_header)?;

    Ok(CramHeader {
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
    ref_seq_map: &sam::header::record::value::Map<
        sam::header::record::value::map::ReferenceSequence,
    >,
) -> Option<[u8; 16]> {
    use noodles_sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
    let raw = ref_seq_map.other_fields().get(&MD5_CHECKSUM)?;
    decode_md5_hex(raw.as_ref())
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
) -> Result<String, CramInputError> {
    use noodles_sam::header::record::value::map::read_group::tag::SAMPLE;
    let mut current: Option<String> = None;
    for (rg_id_bstr, rg_map) in sam_header.read_groups() {
        let sm_raw =
            rg_map
                .other_fields()
                .get(&SAMPLE)
                .ok_or_else(|| CramInputError::MissingSampleTag {
                    path: path.to_path_buf(),
                    read_group_id: String::from_utf8_lossy(rg_id_bstr.as_ref()).into_owned(),
                })?;
        let sm = String::from_utf8_lossy(sm_raw.as_ref()).into_owned();
        match &current {
            None => current = Some(sm),
            Some(existing) if existing != &sm => {
                return Err(CramInputError::MultipleSampleNames {
                    path_a: path.to_path_buf(),
                    sm_a: existing.clone(),
                    path_b: path.to_path_buf(),
                    sm_b: sm,
                });
            }
            Some(_) => {}
        }
    }
    current.ok_or_else(|| CramInputError::MissingSampleTag {
        path: path.to_path_buf(),
        read_group_id: "<no @RG entries>".into(),
    })
}

// ---------------------------------------------------------------------
// FASTA agreement
// ---------------------------------------------------------------------

/// Validate that the indexed FASTA's contigs match `contigs`
/// (name + length, in the same order). The CRAM `@SQ M5` is trusted
/// per `ia/specs/per_sample_caller.md` §"Pre-flight validation"; we do
/// not recompute it from the FASTA.
fn validate_fasta_agreement(
    fasta_path: &Path,
    canonical_contigs: &ContigList,
    cram_path: &Path,
) -> Result<(), CramInputError> {
    let fai_path = with_fai_extension(fasta_path);
    if !fai_path.exists() {
        return Err(CramInputError::MissingFastaIndex {
            fasta_path: fasta_path.to_path_buf(),
        });
    }
    let index = noodles_fasta::fai::fs::read(&fai_path).map_err(|source| CramInputError::Io {
        path: fai_path.clone(),
        source,
    })?;
    let fai_records: &[noodles_fasta::fai::Record] = index.as_ref();

    if fai_records.len() != canonical_contigs.entries.len() {
        return Err(CramInputError::FastaContigMismatch {
            fasta_path: fasta_path.to_path_buf(),
            cram_path: cram_path.to_path_buf(),
            detail: format!(
                "FASTA has {} contigs, CRAM has {}",
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
            return Err(CramInputError::FastaContigMismatch {
                fasta_path: fasta_path.to_path_buf(),
                cram_path: cram_path.to_path_buf(),
                detail: format!(
                    "name disagreement at index {} ('{}' in FASTA vs '{}' in CRAM)",
                    i, fai_name, contig.name
                ),
            });
        }
        let fai_length: u64 = fai_record.length();
        if fai_length != contig.length {
            return Err(CramInputError::FastaContigMismatch {
                fasta_path: fasta_path.to_path_buf(),
                cram_path: cram_path.to_path_buf(),
                detail: format!(
                    "length disagreement at index {} (contig '{}': {} in FASTA vs {} in CRAM)",
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
// CramMergedReader
// ---------------------------------------------------------------------

/// Per-file order tracker: previous `(ref_id, pos)` accepted for that
/// file, used to detect within-file regressions.
type PerFileOrder = Option<(usize, u64)>;

/// Key used by the duplicate-detection window: `(qname, flag, ref_id,
/// pos)` per `ia/specs/per_sample_caller.md` §"Duplicate-read detection
/// across CRAMs".
#[derive(Debug, Clone, PartialEq, Eq)]
struct WindowKey {
    qname: Vec<u8>,
    flag: u16,
    ref_id: usize,
    pos: u64,
}

/// Per-window entry: the key plus the source file index, so the
/// duplicate error can name both files.
#[derive(Debug, Clone)]
struct WindowEntry {
    key: WindowKey,
    source_file_index: usize,
}

type CramRecordsIter = Box<dyn Iterator<Item = io::Result<sam::alignment::RecordBuf>> + Send>;
type CramPeekable = BufferedPeekable<CramRecordsIter, sam::alignment::RecordBuf, io::Error>;

pub struct CramMergedReader {
    peekers: Vec<CramPeekable>,
    paths: Vec<PathBuf>,
    contigs: ContigList,
    sample_name: String,
    config: CramMergedReaderConfig,
    prev_per_file: Vec<PerFileOrder>,
    /// Reads accepted at the current `(ref_id, pos)` — cleared on
    /// advance.
    window: VecDeque<WindowEntry>,
    /// Anchor for the current window: `None` while empty, `Some` at
    /// the position of the entries currently in `window`.
    window_anchor: Option<(usize, u64)>,
    filter_counts: FilterCounts,
}

impl CramMergedReader {
    pub fn new(
        crams: &[PathBuf],
        fasta: &Path,
        config: CramMergedReaderConfig,
    ) -> Result<Self, CramInputError> {
        if crams.is_empty() {
            return Err(CramInputError::Io {
                path: PathBuf::new(),
                source: io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "at least one CRAM input is required",
                ),
            });
        }

        // FASTA repository — built once and shared across decoders.
        let fai_path = with_fai_extension(fasta);
        if !fai_path.exists() {
            return Err(CramInputError::MissingFastaIndex {
                fasta_path: fasta.to_path_buf(),
            });
        }
        let indexed_reader = noodles_fasta::io::indexed_reader::Builder::default()
            .build_from_path(fasta)
            .map_err(|source| CramInputError::Io {
                path: fasta.to_path_buf(),
                source,
            })?;
        let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
        let repository = fasta::Repository::new(adapter);

        // Open every CRAM, validate per-file invariants, and collect
        // per-file headers.
        let mut open_crams: Vec<OpenCram> = Vec::with_capacity(crams.len());
        let mut canonical_contigs: Option<ContigList> = None;
        let mut canonical_sample: Option<String> = None;
        let mut canonical_path: Option<PathBuf> = None;

        for cram_path in crams {
            let mut cram_reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .build_from_path(cram_path)
                .map_err(|source| CramInputError::OpenFailed {
                    path: cram_path.clone(),
                    source,
                })?;

            // Read file definition first to gate on CRAM major version
            // before any container is decoded.
            let file_definition = cram_reader.read_file_definition().map_err(|source| {
                CramInputError::OpenFailed {
                    path: cram_path.clone(),
                    source,
                }
            })?;
            let version = file_definition.version();
            if version.major() != 3 {
                return Err(CramInputError::UnsupportedCramVersion {
                    path: cram_path.clone(),
                    major: version.major(),
                    minor: version.minor(),
                });
            }

            let sam_header =
                cram_reader
                    .read_file_header()
                    .map_err(|source| CramInputError::OpenFailed {
                        path: cram_path.clone(),
                        source,
                    })?;
            let cram_header = extract_header(cram_path, &sam_header)?;

            // Cross-file checks: contigs and sample name must be
            // identical across every file in the input.
            match &canonical_contigs {
                None => {
                    canonical_contigs = Some(cram_header.contigs.clone());
                    canonical_path = Some(cram_path.clone());
                }
                Some(existing) => {
                    if let Err(detail) = existing.first_disagreement(&cram_header.contigs) {
                        return Err(CramInputError::ContigListMismatch {
                            reference_path: canonical_path.clone().unwrap_or_default(),
                            other_path: cram_path.clone(),
                            detail,
                        });
                    }
                }
            }
            match &canonical_sample {
                None => canonical_sample = Some(cram_header.sample_name.clone()),
                Some(existing) if existing != &cram_header.sample_name => {
                    return Err(CramInputError::MultipleSampleNames {
                        path_a: canonical_path.clone().unwrap_or_default(),
                        sm_a: existing.clone(),
                        path_b: cram_path.clone(),
                        sm_b: cram_header.sample_name,
                    });
                }
                _ => {}
            }

            // Build the owned record iterator and pack it into an
            // `OpenCram`. The reader and header are moved into the
            // owned iterator — neither escapes from this scope.
            let owned: CramRecordsIter = Box::new(OwnedCramRecords {
                reader: cram_reader,
                header: sam_header,
                repository: repository.clone(),
                container: cram::io::reader::Container::default(),
                pending: Vec::new().into_iter(),
            });
            open_crams.push(OpenCram {
                path_for_errors: cram_path.clone(),
                records: owned,
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
            canonical_path.as_deref().unwrap_or(fasta),
        )?;

        Self::from_open_crams(open_crams, canonical_contigs, canonical_sample, config)
    }

    /// Crate-internal core. Takes pre-built per-CRAM record streams
    /// plus the canonical contig list and sample name (the caller is
    /// responsible for cross-file header validation — `new` does this
    /// before calling here). All I/O for opening lives outside this
    /// function, which is why it is the single point of entry tests
    /// drive against.
    pub(crate) fn from_open_crams(
        open_crams: Vec<OpenCram>,
        contigs: ContigList,
        sample_name: String,
        config: CramMergedReaderConfig,
    ) -> Result<Self, CramInputError> {
        let n = open_crams.len();
        let mut peekers = Vec::with_capacity(n);
        let mut paths = Vec::with_capacity(n);
        for open in open_crams {
            paths.push(open.path_for_errors);
            peekers.push(BufferedPeekable::with_buffer_size(
                open.records,
                PER_PEEKER_BUFFER_SIZE,
            ));
        }
        Ok(Self {
            peekers,
            paths,
            contigs,
            sample_name,
            config,
            prev_per_file: vec![None; n],
            window: VecDeque::new(),
            window_anchor: None,
            filter_counts: FilterCounts::default(),
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
}

// ---------------------------------------------------------------------
// Iterator + merge logic
// ---------------------------------------------------------------------

impl Iterator for CramMergedReader {
    type Item = Result<MappedRead, CramInputError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Refill peekers' heads, dropping any head that fails the
            // pre-decode filter cascade. Each dropped head increments
            // its `FilterCounts` bucket.
            if let Err(e) = self.refill_heads() {
                return Some(Err(e));
            }

            // Pick the smallest surviving head; `None` means every
            // peeker is exhausted and we are done.
            let chosen_idx = self.argmin_head()?;

            // Validate per-file order against this file's previous
            // accepted `(ref_id, pos)`.
            let head = match self.peek_head_keys(chosen_idx) {
                Ok(Some(keys)) => keys,
                Ok(None) => continue,
                Err(e) => return Some(Err(e)),
            };
            if let Some((prev_ref, prev_pos)) = self.prev_per_file[chosen_idx]
                && (head.ref_id, head.pos) < (prev_ref, prev_pos)
            {
                return Some(Err(CramInputError::OutOfOrderRead {
                    path: self.paths[chosen_idx].clone(),
                    qname: String::from_utf8_lossy(&head.qname).into_owned(),
                    prev_pos: encode_order_key(prev_ref, prev_pos),
                    this_pos: encode_order_key(head.ref_id, head.pos),
                }));
            }

            // Maintain the duplicate-detection window: clear if we
            // advanced past the current anchor.
            self.advance_window_if_needed(head.ref_id, head.pos);

            // Detect duplicates against the current window.
            let new_key = WindowKey {
                qname: head.qname.clone(),
                flag: head.flag,
                ref_id: head.ref_id,
                pos: head.pos,
            };
            if let Some(other) = self.window.iter().find(|entry| entry.key == new_key) {
                let other_path = self.paths[other.source_file_index].clone();
                return Some(Err(CramInputError::DuplicateReadAcrossFiles {
                    qname: String::from_utf8_lossy(&new_key.qname).into_owned(),
                    path_a: other_path,
                    path_b: self.paths[chosen_idx].clone(),
                    ref_id: new_key.ref_id,
                    pos: new_key.pos,
                }));
            }

            // Consume the chosen peeker's head.
            let record = match self.peekers[chosen_idx].next() {
                Some(Ok(rb)) => rb,
                Some(Err(e)) => {
                    return Some(Err(CramInputError::MalformedRecord {
                        path: self.paths[chosen_idx].clone(),
                        qname: String::from_utf8_lossy(&head.qname).into_owned(),
                        source: e,
                    }));
                }
                None => continue,
            };

            // Convert RecordBuf -> MappedRead.
            let mapped = match record_buf_to_mapped_read(&record, chosen_idx) {
                Ok(m) => m,
                Err(e) => {
                    return Some(Err(CramInputError::MalformedRecord {
                        path: self.paths[chosen_idx].clone(),
                        qname: String::from_utf8_lossy(&head.qname).into_owned(),
                        source: e,
                    }));
                }
            };

            // Post-decode filter: minimum read length (depends on
            // SEQ length, only available after decode).
            if let Some(min) = self.config.min_read_length
                && (mapped.seq.len() as u32) < min
            {
                self.filter_counts.too_short += 1;
                continue;
            }

            // Accept: update window, prev_per_file.
            self.window.push_back(WindowEntry {
                key: new_key,
                source_file_index: chosen_idx,
            });
            self.prev_per_file[chosen_idx] = Some((mapped.ref_id, mapped.pos));
            return Some(Ok(mapped));
        }
    }
}

impl CramMergedReader {
    fn advance_window_if_needed(&mut self, ref_id: usize, pos: u64) {
        match self.window_anchor {
            Some(anchor) if anchor != (ref_id, pos) => {
                self.window.clear();
                self.window_anchor = Some((ref_id, pos));
            }
            None => {
                self.window_anchor = Some((ref_id, pos));
            }
            _ => {}
        }
    }

    fn refill_heads(&mut self) -> Result<(), CramInputError> {
        let config = self.config;
        for idx in 0..self.peekers.len() {
            loop {
                let drop_outcome = match self.peekers[idx].peek() {
                    Ok(Some(rb)) => classify_pre_decode(&config, rb),
                    Ok(None) => break,
                    Err(e) => {
                        return Err(CramInputError::MalformedRecord {
                            path: self.paths[idx].clone(),
                            qname: String::new(),
                            source: e,
                        });
                    }
                };
                match drop_outcome {
                    PreDecodeOutcome::Keep => break,
                    PreDecodeOutcome::Drop(bucket) => {
                        self.bump_filter_count(bucket);
                        // Consume the dropped head and re-peek.
                        let _ = self.peekers[idx].next();
                    }
                }
            }
        }
        Ok(())
    }

    /// Inspect a peekable head without consuming it. Returns the
    /// merge-order key plus the few fields we need to perform the
    /// out-of-order and duplicate checks. `Ok(None)` means the
    /// peeker is exhausted (peek is None at this point — should not
    /// happen between `refill_heads` and `argmin_head`, but treat
    /// defensively).
    fn peek_head_keys(&mut self, idx: usize) -> Result<Option<HeadKey>, CramInputError> {
        match self.peekers[idx].peek() {
            Ok(Some(rb)) => match head_key(rb) {
                Some(key) => Ok(Some(key)),
                None => {
                    // Unmapped or malformed head leaked through; treat
                    // it as a malformed record so the user sees a
                    // pointed error instead of a silent drop.
                    let qname = rb
                        .name()
                        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                        .unwrap_or_default();
                    Err(CramInputError::MalformedRecord {
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
            Err(e) => Err(CramInputError::MalformedRecord {
                path: self.paths[idx].clone(),
                qname: String::new(),
                source: e,
            }),
        }
    }

    fn argmin_head(&mut self) -> Option<usize> {
        let mut best: Option<(usize, (usize, u64))> = None;
        for idx in 0..self.peekers.len() {
            let Ok(Some(rb)) = self.peekers[idx].peek() else {
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
        }
    }
}

fn classify_pre_decode(
    config: &CramMergedReaderConfig,
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
    // 3. Supplementary (~1-5%).
    if config.drop_supplementary && (flag & FLAG_SUPPLEMENTARY) != 0 {
        return PreDecodeOutcome::Drop(FilterBucket::Supplementary);
    }
    // 4. Secondary (~1-5%).
    if config.drop_secondary && (flag & FLAG_SECONDARY) != 0 {
        return PreDecodeOutcome::Drop(FilterBucket::Secondary);
    }
    // 5. Unmapped (~0-5%).
    if config.drop_unmapped && (flag & FLAG_UNMAPPED) != 0 {
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

fn encode_order_key(ref_id: usize, pos: u64) -> u64 {
    // Pack (ref_id, pos) into a single u64 for the OutOfOrderRead error
    // message. `ref_id` lives in the high 32 bits.
    ((ref_id as u64) << 32) | (pos & 0xFFFF_FFFF)
}

#[derive(Debug, Clone, Copy)]
enum PreDecodeOutcome {
    Keep,
    Drop(FilterBucket),
}

#[derive(Debug, Clone, Copy)]
enum FilterBucket {
    Unmapped,
    Secondary,
    Supplementary,
    QcFail,
    Duplicate,
    LowMapq,
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
        source_file_index,
    })
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
    use crate::per_sample_caller::record_specs::{
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

    // --- Group A: via from_open_crams --------------------------------

    fn run_to_completion(
        mut reader: CramMergedReader,
    ) -> (Vec<MappedRead>, FilterCounts, Option<CramInputError>) {
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
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
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
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a, cram_b],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
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
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a, cram_b],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
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
        let mut reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
        )
        .expect("reader");
        let first = reader.next().expect("first").expect("ok");
        assert_eq!(first.pos, 200);
        let err = reader.next().expect("err item").expect_err("expected err");
        match err {
            CramInputError::OutOfOrderRead { qname, .. } => {
                assert_eq!(qname, "R2");
            }
            other => panic!("unexpected error: {:?}", other),
        }
    }

    #[test]
    fn a5_duplicate_read_across_streams() {
        let dup = pass_record("R1", 0, 100);
        let cram_a = open_cram_from_records("a.cram", vec![record_spec(dup.clone())]);
        let cram_b = open_cram_from_records("b.cram", vec![record_spec(dup)]);
        let mut reader = CramMergedReader::from_open_crams(
            vec![cram_a, cram_b],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
        )
        .expect("reader");
        let first = reader.next().expect("first").expect("ok");
        assert_eq!(first.pos, 100);
        let err = reader.next().expect("err item").expect_err("expected err");
        match err {
            CramInputError::DuplicateReadAcrossFiles {
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
    fn a6_duplicate_window_clears_on_advance() {
        let cram_a = open_cram_from_records(
            "a.cram",
            vec![
                record_spec(pass_record("R1", 0, 100)),
                record_spec(pass_record("R1", 0, 200)),
            ],
        );
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
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
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].qname, b"R_high");
        assert_eq!(counts.low_mapq, 2);
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
        let cfg = CramMergedReaderConfig {
            min_mapq: None,
            ..Default::default()
        };
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            cfg,
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
        type FlagCase = (
            fn(&mut CramMergedReaderConfig),
            &'static str,
            fn(&FilterCounts) -> u64,
        );
        let cases: &[FlagCase] = &[
            (
                |c| {
                    c.drop_secondary = false;
                    c.drop_supplementary = false;
                    c.drop_qc_fail = false;
                    c.drop_duplicate = false;
                },
                "unmapped",
                |c| c.unmapped,
            ),
            (
                |c| {
                    c.drop_unmapped = false;
                    c.drop_supplementary = false;
                    c.drop_qc_fail = false;
                    c.drop_duplicate = false;
                },
                "secondary",
                |c| c.secondary,
            ),
            (
                |c| {
                    c.drop_unmapped = false;
                    c.drop_secondary = false;
                    c.drop_qc_fail = false;
                    c.drop_duplicate = false;
                },
                "supplementary",
                |c| c.supplementary,
            ),
            (
                |c| {
                    c.drop_unmapped = false;
                    c.drop_secondary = false;
                    c.drop_supplementary = false;
                    c.drop_duplicate = false;
                },
                "qc_fail",
                |c| c.qc_fail,
            ),
            (
                |c| {
                    c.drop_unmapped = false;
                    c.drop_secondary = false;
                    c.drop_supplementary = false;
                    c.drop_qc_fail = false;
                },
                "duplicate",
                |c| c.duplicate,
            ),
        ];
        for (mutator, missing_qname, count_fn) in cases {
            let mut cfg = CramMergedReaderConfig::default();
            mutator(&mut cfg);
            let recs: Vec<_> = six_flagged_records().into_iter().map(record_spec).collect();
            let cram_a = open_cram_from_records("a.cram", recs);
            let reader = CramMergedReader::from_open_crams(
                vec![cram_a],
                default_contigs(),
                "sample".into(),
                cfg,
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
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
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
    fn a10_all_flag_drops_disabled_passes_everything() {
        let recs: Vec<_> = six_flagged_records().into_iter().map(record_spec).collect();
        let cram_a = open_cram_from_records("a.cram", recs);
        let cfg = CramMergedReaderConfig {
            min_mapq: None,
            min_read_length: None,
            drop_unmapped: false,
            drop_secondary: false,
            drop_supplementary: false,
            drop_qc_fail: false,
            drop_duplicate: false,
        };
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            cfg,
        )
        .expect("reader");
        let (out, counts, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert_eq!(out.len(), 6);
        assert_eq!(counts, FilterCounts::default());
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
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
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
        let cfg = CramMergedReaderConfig {
            min_read_length: None,
            ..Default::default()
        };
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            cfg,
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
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
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
        let mut reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
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
        let reader = CramMergedReader::from_open_crams(
            vec![cram_a, cram_b],
            default_contigs(),
            "sample".into(),
            CramMergedReaderConfig::default(),
        )
        .expect("reader");
        let (out, _, err) = run_to_completion(reader);
        assert!(err.is_none());
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].qname, b"R1");
        assert_eq!(out[1].qname, b"R2");
    }

    #[test]
    fn a15_position_above_u32_max_round_trips_through_mapped_read() {
        let big_pos: u64 = (u32::MAX as u64) + 100;
        let rec = pass_record("R", 0, big_pos);
        let cram_a = open_cram_from_records("a.cram", vec![record_spec(rec)]);
        let mut contigs = default_contigs();
        contigs.entries[0].length = big_pos + 1_000;
        let mut reader = CramMergedReader::from_open_crams(
            vec![cram_a],
            contigs,
            "sample".into(),
            CramMergedReaderConfig::default(),
        )
        .expect("reader");
        let read = reader.next().expect("ok").expect("ok");
        assert_eq!(read.pos, big_pos);
        assert!(reader.next().is_none());
    }

    // --- Group B: via new (real CRAM + FASTA) ------------------------

    use crate::per_sample_caller::cram_files::{
        ContigSpec, HeaderOverrides, build_cram, build_cram_with_major_version, build_fasta,
    };

    /// `expect_err` requires `Debug` on the success type;
    /// `CramMergedReader` carries a boxed `Send` iterator and cannot be
    /// `Debug`. This helper unwraps the error variant and panics with
    /// the supplied message if the result was `Ok`.
    fn err_or_panic(
        result: Result<CramMergedReader, CramInputError>,
        message: &str,
    ) -> CramInputError {
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
        let (_cram_dir, cram_path) = build_cram(
            &fasta_path,
            &contigs,
            &HeaderOverrides {
                read_groups: vec![("rg0".into(), Some("s1".into()))],
                ..Default::default()
            },
            &records,
        )
        .expect("build_cram");
        let reader =
            CramMergedReader::new(&[cram_path], &fasta_path, CramMergedReaderConfig::default())
                .expect("CramMergedReader::new");
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
        let (_cram_dir, cram_path) = build_cram_with_major_version(4, 0).expect("build forced");
        let err = err_or_panic(
            CramMergedReader::new(
                std::slice::from_ref(&cram_path),
                &fasta_path,
                CramMergedReaderConfig::default(),
            ),
            "expected UnsupportedCramVersion",
        );
        match err {
            CramInputError::UnsupportedCramVersion { major, minor, path } => {
                assert_eq!(major, 4);
                assert_eq!(minor, 0);
                assert_eq!(path, cram_path);
            }
            other => panic!("unexpected: {:?}", other),
        }
    }

    #[test]
    fn b3_non_coordinate_sort_rejected() {
        let contigs = one_contig_chr1();
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).expect("build_fasta");
        let records = vec![pass_record_for_b("R1", 0, 100)];
        let (_cram_dir, cram_path) = build_cram(
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
            CramMergedReader::new(&[cram_path], &fasta_path, CramMergedReaderConfig::default()),
            "expected NotCoordinateSorted",
        );
        match err {
            CramInputError::NotCoordinateSorted { sort_order, .. } => {
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
            CramMergedReader::new(&[c1], &fa1, CramMergedReaderConfig::default()).expect("reader");
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
            CramMergedReader::new(&[c2], &fa2, CramMergedReaderConfig::default()),
            "expected MissingSampleTag",
        );
        assert!(
            matches!(err, CramInputError::MissingSampleTag { .. }),
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
            CramMergedReader::new(&[c3a, c3b], &fa3, CramMergedReaderConfig::default()),
            "expected MultipleSampleNames",
        );
        assert!(
            matches!(err, CramInputError::MultipleSampleNames { .. }),
            "got {:?}",
            err
        );
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
            CramMergedReader::new(&[a1, a2], &fa_two, CramMergedReaderConfig::default()),
            "expected ContigListMismatch (name)",
        );
        match err {
            CramInputError::ContigListMismatch { detail, .. } => {
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
            CramMergedReader::new(&[b1, b2], &fa2, CramMergedReaderConfig::default()),
            "expected ContigListMismatch (length)",
        );
        match err {
            CramInputError::ContigListMismatch { detail, .. } => {
                assert!(
                    detail.contains("length disagreement"),
                    "detail = {}",
                    detail
                );
            }
            CramInputError::FastaContigMismatch { detail, .. } => {
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
            CramMergedReader::new(&[c1, c2], &fa, CramMergedReaderConfig::default()),
            "expected ContigListMismatch (md5)",
        );
        match err {
            CramInputError::ContigListMismatch { detail, .. } => {
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
        let (_cram_dir, cram_path) =
            build_cram(&fa, &contigs, &HeaderOverrides::default(), &records).expect("cram");
        CramMergedReader::new(
            std::slice::from_ref(&cram_path),
            &fa,
            CramMergedReaderConfig::default(),
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
            CramMergedReader::new(&[cram_path], &fa2, CramMergedReaderConfig::default()),
            "expected MissingFastaIndex",
        );
        assert!(
            matches!(err, CramInputError::MissingFastaIndex { .. }),
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
            CramMergedReader::new(&[cram_b], &fa3, CramMergedReaderConfig::default()),
            "expected FastaContigMismatch",
        );
        assert!(
            matches!(err, CramInputError::FastaContigMismatch { .. }),
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
        let (_cram_dir, cram_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records)
                .expect("build_cram");
        let reader =
            CramMergedReader::new(&[cram_path], &fasta_path, CramMergedReaderConfig::default())
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
}
