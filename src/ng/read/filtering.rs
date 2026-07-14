//! ng step 1 — read filtering: the whole-read keep/drop prelude. It takes the
//! reads of one sample's alignment file (BAM/CRAM) and yields the subset worth
//! carrying forward, plus a running tally of what was dropped and why. Every
//! decision is per-read, locus-independent, and content-preserving — filtering
//! *selects* reads, it never rewrites them.
//!
//! Design: `doc/devel/ng/spec/read_filtering.md` (the "why"),
//! `doc/devel/ng/arch/read_filtering.md` (types & interfaces).
//!
//! Read filtering is a **port** of the production filter stack in
//! [`crate::bam::alignment_input`]: it reuses that module's pure predicates
//! (`read_exceeds_mismatch_fraction`, `cigar_is_bad`), its `FLAG_*` /
//! `DEFAULT_*` constants, and its `RecordBuf → MappedRead` decode path as-is,
//! and supplies only its own driver and config.
//!
//! The whole step is here: the step-1-local types, the two-phase `verdict_*`
//! cascade, the `RawRecord`/`RecordSource` seam with the noodles BAM/CRAM
//! adapters, and the `ReadFilter` iterator that drives them into a stream of
//! kept reads with a running drop tally.

use crate::bam::alignment_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
    DEFAULT_MISMATCH_BQ_FLOOR, FLAG_DUPLICATE, FLAG_QC_FAIL, FLAG_SECONDARY, FLAG_SUPPLEMENTARY,
    FLAG_UNMAPPED, MappedRead, cigar_is_bad, cigar_ref_span, read_exceeds_mismatch_fraction,
    record_buf_to_mapped_read,
};
use crate::ng::ref_seq::{RawRefSeq, RefSeqError};
use crate::ng::types::{BaseQual, Bp, ContigId, MapQual, MismatchFraction};
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;
use std::io;

/// The filtering policy: which filters are active and their thresholds. Minimal
/// by design — one field per active filter, no dormant levers (downsampling,
/// read pooling enter only when they enter the pipeline). [`Default`] is the
/// production policy the lab runs with, its thresholds the reused `DEFAULT_*`
/// constants from [`crate::bam::alignment_input`]. Mirrors the filtering subset
/// of the existing `AlignmentMergedReaderConfig`.
///
/// `Option<T>` is the "no threshold" state, never a sentinel: `None` means *no
/// minimum*, `Some(q)` means *drop below `q`* — so `Some(0)` (drop nothing, but
/// a threshold is set) stays structurally distinct from `None` (no threshold).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReadFilterConfig {
    /// `None` = no minimum; `Some(q)` = drop reads with MAPQ `< q` (filter #2).
    pub min_mapq: Option<MapQual>,
    /// `None` = no minimum; `Some(n)` = drop reads shorter than `n` bp (#7).
    pub min_read_length: Option<Bp>,
    /// Drop reads flagged QC-fail (`FLAG_QC_FAIL`) — filter #6.
    pub drop_qc_fail: bool,
    /// Drop reads flagged as PCR/optical duplicates (`FLAG_DUPLICATE`) — #1.
    pub drop_duplicate: bool,
    /// `None` = filter #8 off (no reference access at all); `Some(x)` = drop a
    /// read whose quality-clearing `M`-op mismatch fraction exceeds `x`.
    pub max_read_mismatch_fraction: Option<MismatchFraction>,
    /// BQ floor below which a mismatch does not count toward filter #8. Only
    /// meaningful when `max_read_mismatch_fraction` is `Some`.
    pub mismatch_bq_floor: BaseQual,
}

impl Default for ReadFilterConfig {
    fn default() -> Self {
        Self {
            min_mapq: Some(MapQual(DEFAULT_MIN_MAPQ)),
            min_read_length: Some(Bp(DEFAULT_MIN_READ_LENGTH)),
            drop_qc_fail: true,
            drop_duplicate: true,
            // PANIC-FREE: the default fraction is a known-good in-range constant
            // (0.10, `DEFAULT_MAX_READ_MISMATCH_FRACTION`), so the checked
            // constructor cannot fail here. Guarded by the
            // `default_config_reproduces_the_production_filter_policy` test, which
            // exercises this exact path.
            max_read_mismatch_fraction: Some(
                MismatchFraction::try_new(DEFAULT_MAX_READ_MISMATCH_FRACTION)
                    .expect("DEFAULT_MAX_READ_MISMATCH_FRACTION is in [0, 1]"),
            ),
            mismatch_bq_floor: BaseQual(DEFAULT_MISMATCH_BQ_FLOOR),
        }
    }
}

/// The verdict for one read. `Keep` carries the read on; `Drop` records which
/// filter fired — the first one, per the hit-rate order — for the tally.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterVerdict {
    Keep,
    Drop(DropReason),
}

/// Which filter dropped a read. Variant names line up 1:1 with the
/// [`ReadFilterCounts`] fields, so a drop maps to exactly one counter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    Duplicate,
    LowMapq,
    Supplementary,
    Secondary,
    Unmapped,
    QcFail,
    TooShort,
    HighMismatchFraction,
    BadCigar,
}

/// A per-sample tally of the filtering pass — one counter per drop reason, plus
/// the kept count. The ng port of `FilterCounts`. Surfacing every drop is the
/// "no silent caps" discipline: a read that vanished must be accounted for. It
/// is a **running** tally — readable at any point, final once the input is
/// exhausted.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReadFilterCounts {
    pub kept: u64,
    pub duplicate: u64,
    pub low_mapq: u64,
    pub supplementary: u64,
    pub secondary: u64,
    pub unmapped: u64,
    pub qc_fail: u64,
    pub too_short: u64,
    pub high_mismatch_fraction: u64,
    pub bad_cigar: u64,
}

impl ReadFilterCounts {
    /// Tally one drop against its counter. The exhaustive `match` is the guard for
    /// the documented `DropReason` ↔ `ReadFilterCounts` 1:1 mapping: adding a
    /// `DropReason` variant without a counter here is a compile error, so the two
    /// cannot silently desync (mirrors production `FilterCounts::record_drop`).
    fn record_drop(&mut self, reason: DropReason) {
        match reason {
            DropReason::Duplicate => self.duplicate += 1,
            DropReason::LowMapq => self.low_mapq += 1,
            DropReason::Supplementary => self.supplementary += 1,
            DropReason::Secondary => self.secondary += 1,
            DropReason::Unmapped => self.unmapped += 1,
            DropReason::QcFail => self.qc_fail += 1,
            DropReason::TooShort => self.too_short += 1,
            DropReason::HighMismatchFraction => self.high_mismatch_fraction += 1,
            DropReason::BadCigar => self.bad_cigar += 1,
        }
    }
}

/// Phase one of the cascade — the flag/MAPQ filters (#1–#6), decided on an
/// undecoded record's `flag` and `mapq` alone. Reference-free and decode-free:
/// `Keep` means "decode and continue to phase two", `Drop` is charged to the
/// first filter that fires. Order is identical to production's
/// `classify_pre_decode` (hit-rate-ordered: duplicate, low-MAPQ, supplementary,
/// secondary, unmapped, QC-fail).
///
/// `mapq` is already resolved: SAM's "unavailable" (`0xFF`) is mapped to
/// `MapQual(0)` by the record source (Milestone C), so a non-zero `min_mapq`
/// drops it — matching production. `flag` is the raw SAM bitfield
/// (`MappedRead.flag`), tested against the reused `FLAG_*` constants.
fn verdict_pre_decode(flag: u16, mapq: MapQual, config: &ReadFilterConfig) -> FilterVerdict {
    // 1. Duplicate — a PCR/optical copy of another molecule (toggle).
    if config.drop_duplicate && (flag & FLAG_DUPLICATE) != 0 {
        return FilterVerdict::Drop(DropReason::Duplicate);
    }
    // 2. Low MAPQ — the aligner is unsure of the placement (toggle via threshold).
    if let Some(min) = config.min_mapq
        && mapq < min
    {
        return FilterVerdict::Drop(DropReason::LowMapq);
    }
    // 3. Supplementary — a chunk of a chimeric read (unconditional).
    if (flag & FLAG_SUPPLEMENTARY) != 0 {
        return FilterVerdict::Drop(DropReason::Supplementary);
    }
    // 4. Secondary — a duplicate projection of a primary alignment (unconditional).
    if (flag & FLAG_SECONDARY) != 0 {
        return FilterVerdict::Drop(DropReason::Secondary);
    }
    // 5. Unmapped — no placement, so no allele evidence (unconditional).
    if (flag & FLAG_UNMAPPED) != 0 {
        return FilterVerdict::Drop(DropReason::Unmapped);
    }
    // 6. QC fail — the sequencer/pipeline flagged the read (toggle).
    if config.drop_qc_fail && (flag & FLAG_QC_FAIL) != 0 {
        return FilterVerdict::Drop(DropReason::QcFail);
    }
    FilterVerdict::Keep
}

/// Phase two of the cascade — the decode-dependent filters, run on a decoded
/// [`MappedRead`] only after it clears phase one. Evaluated **cheapest-first**:
///
/// 1. **#7 too-short** — one length compare, no reference.
/// 2. **#9 bad-CIGAR** — a pure CIGAR scan, no reference.
/// 3. **#8 high-mismatch** — a reference fetch plus a per-base walk, and only
///    when `max_read_mismatch_fraction` is `Some`.
///
/// This is a deliberate reordering relative to the spec's #7/#8/#9 *table*: it
/// puts the one reference-touching filter (#8) last, so a read dropped for
/// being too short or for a malformed CIGAR never pays the reference fetch and
/// base walk. It honours the spec's stated "cheapest, most-often-firing first"
/// principle (spec §3), charges a both-failing read to the root cause
/// (`BadCigar`) rather than the symptom (`HighMismatchFraction`), and leaves the
/// keep/drop *set* unchanged (the filters are independent — a read failing any
/// is dropped regardless of order). The mismatch fraction is measured on the
/// read's own (un-left-aligned) CIGAR: left-alignment only shifts indels across
/// equal bases, so it does not change the match/mismatch tally, and it is
/// deferred to `pileup/` anyway (spec §6).
///
/// `ref_buf` is a caller-owned scratch buffer the reference bytes are read into
/// (reused across reads, so #8 costs no per-read allocation). It is written only
/// when #8 actually runs. A [`RefSeqError`] from the fetch is **fatal to the
/// run**, propagated rather than swallowed into a drop or a keep. This includes
/// an `OutOfBounds` window running past the contig end: a validly-aligned read
/// never covers reference positions the contig does not have, so an
/// out-of-bounds fetch signals a malformed record — and the fatal model treats
/// it, like a truncated file, as corrupt input to fail loudly on rather than
/// filter around (spec §7).
fn verdict_post_decode(
    read: &MappedRead,
    reference: &impl RawRefSeq,
    config: &ReadFilterConfig,
    ref_buf: &mut Vec<u8>,
) -> Result<FilterVerdict, RefSeqError> {
    // #7 — too short (decoded SEQ length). Cheapest: no CIGAR walk, no reference.
    if let Some(min) = config.min_read_length
        && (read.seq.len() as u32) < min.get()
    {
        return Ok(FilterVerdict::Drop(DropReason::TooShort));
    }

    // #9 — bad CIGAR (adjacent I/D, or a boundary deletion). Cheap: a pure scan
    // of the aligner's CIGAR, no reference. Run before #8 so a malformed read
    // never pays the reference fetch below.
    if cigar_is_bad(&read.cigar) {
        return Ok(FilterVerdict::Drop(DropReason::BadCigar));
    }

    // #8 — high mismatch fraction. The only reference-dependent filter, and the
    // only post-decode one that allocates work; skipped entirely when disabled.
    if let Some(max) = config.max_read_mismatch_fraction {
        let ref_span = cigar_ref_span(&read.cigar);
        // PANIC-FREE: reference coordinates are `u32` at the RefSeq boundary
        // (ref_seq.md, Decision 3). `ref_id` indexes the `u32` contig table and a
        // mapped read's 1-based position fits `u32` for any real contig
        // (< 4.29 Gbp), so neither conversion truncates on real input. A value
        // that did not fit would be a corrupt record — failing loudly is the
        // intended response under the fatal error model, not a silent `as`
        // truncation that would fetch the wrong window and mis-verdict the read.
        let contig = ContigId(u32::try_from(read.ref_id).expect("ref_id fits u32"));
        let pos = u32::try_from(read.pos).expect("read position fits u32");
        // Reads raw (un-canonicalised) bytes, matching production's
        // `RawContigRefCache` path so the ported filter behaves identically.
        // A zero span yields an empty slice → `read_exceeds_mismatch_fraction`
        // has no comparable bases and keeps the read (same as production).
        reference.fetch_raw_into(contig, pos, ref_span, ref_buf)?;
        if read_exceeds_mismatch_fraction(
            &read.cigar,
            &read.seq,
            &read.qual,
            ref_buf,
            config.mismatch_bq_floor.get(),
            max.get(),
        ) {
            return Ok(FilterVerdict::Drop(DropReason::HighMismatchFraction));
        }
    }

    Ok(FilterVerdict::Keep)
}

// ---------------------------------------------------------------------
// The record-source seam (input edge)
// ---------------------------------------------------------------------

/// A borrowed view of one alignment record — the seam that lets the flag/MAPQ
/// cascade (#1–#6) run *before* the record is decoded. [`Self::flag`] and
/// [`Self::mapq`] are cheap field reads on the still-packed record; [`Self::decode`]
/// is the expensive phase (base/quality decode + adaptor-boundary annotation →
/// [`MappedRead`]), run only on the reads that clear the pre-decode gate. It
/// borrows `&self`, not `self`, so the underlying buffer stays reusable for the
/// next read.
pub trait RawRecord {
    /// The SAM flag bitfield — the same `u16` [`MappedRead::flag`] carries; read
    /// by filters #1, #3–#6.
    fn flag(&self) -> u16;
    /// The SAM mapping quality; filter #2. "Unavailable" (SAM `0xFF`) resolves to
    /// [`MapQual`]`(0)`, matching production, so any non-zero minimum drops it.
    fn mapq(&self) -> MapQual;
    /// Decode the record into an owned [`MappedRead`] (filters #7–#9 read the
    /// result). Fallible because the reused decode path
    /// ([`record_buf_to_mapped_read`]) is: a record reaching here has passed the
    /// pre-decode cascade (so it is mapped), but a corrupt record — the unmapped
    /// flag clear yet no position — surfaces here as an `Err` that, like the #8
    /// fetch and `RecordSource::read_next`, is **fatal to the run** rather than a
    /// per-read drop or a panic (spec §7). (The spec's illustrative signature
    /// showed this infallible; it is fallible in practice because the reused
    /// production decoder is.)
    fn decode(&self) -> io::Result<MappedRead>;
}

/// The filter's input: fills a caller-owned buffer with the next record, reusing
/// its allocations. Modelled as a source that *fills* a buffer rather than an
/// `Iterator<Item = RawRecord>` precisely so one buffer survives the whole pass
/// (a std iterator's owned `Item` cannot borrow a reused buffer — the
/// lending-iterator problem, spec §5). `Ok(true)` = filled, `Ok(false)` = end of
/// input; an `Err` is **fatal to the run**.
pub trait RecordSource {
    /// The reused buffer type — a [`RawRecord`] the source refills in place. The
    /// `Default` bound is load-bearing: the [`ReadFilter`] iterator seeds *one*
    /// buffer with `Default::default()` and hands the same `&mut` to
    /// [`Self::read_next`] on every read, so the whole pass allocates one record.
    type Record: RawRecord + Default;
    /// The alignment file's SAM header. Its `@SQ` list defines the contigs the
    /// reads may reference; [`ReadFilter::new`] validates that every one resolves
    /// in the reference up front, so an in-loop reference error later means
    /// genuinely corrupt input rather than a mismatched reference.
    fn header(&self) -> &sam::Header;
    /// Fill `buf` with the next record, reusing its allocations. `Ok(true)` =
    /// filled; `Ok(false)` = end of input, after which `buf` holds an unspecified
    /// (stale) record the caller must not read; `Err` is fatal to the run.
    fn read_next(&mut self, buf: &mut Self::Record) -> io::Result<bool>;
}

/// An ng-owned adapter viewing a noodles [`RecordBuf`] as a [`RawRecord`]. The
/// production `RecordSource` refills the inner `record` in place;
/// [`RawRecord::flag`]/[`RawRecord::mapq`] are cheap field reads on the undecoded
/// record, and [`RawRecord::decode`] runs the reused
/// [`record_buf_to_mapped_read`] path (which itself reuses
/// `compute_adaptor_boundary`). `source_file_index` is the per-file tag stamped
/// onto the decoded [`MappedRead`]; the source sets it before each read.
///
/// ng-owned by design (spec §7, decision a): the dependency points ng → existing
/// code, so production never learns about ng.
///
/// Both fields are `pub(crate)` internal state, not a caller-set surface:
/// `source_file_index` is (re)stamped by the source on every [`RecordSource::read_next`],
/// and `record` is refilled in place — a caller-written value would just be
/// overwritten. No `Clone`: the whole point is to reuse one buffer, and cloning
/// would deep-copy the `RecordBuf` it exists to avoid copying.
#[derive(Debug, Default)]
pub struct NoodlesRawRecord {
    /// The reused undecoded record buffer.
    pub(crate) record: RecordBuf,
    /// The 0-based index of the file this record came from (stamped onto the
    /// decoded [`MappedRead`]).
    pub(crate) source_file_index: usize,
}

impl RawRecord for NoodlesRawRecord {
    fn flag(&self) -> u16 {
        self.record.flags().bits()
    }

    fn mapq(&self) -> MapQual {
        // `map(u8::from).unwrap_or(0)`: SAM `0xFF` (unavailable) decodes to `None`
        // → treated as MAPQ 0, matching production `classify_pre_decode`.
        MapQual(self.record.mapping_quality().map(u8::from).unwrap_or(0))
    }

    fn decode(&self) -> io::Result<MappedRead> {
        record_buf_to_mapped_read(&self.record, self.source_file_index)
    }
}

/// A [`RecordSource`] over a noodles BAM reader. `read_next` reads the next record
/// into the reused buffer via noodles' `read_record_buf` (true buffer reuse — no
/// per-read allocation), returning `Ok(false)` at end of input. It is
/// **unfiltered**: it hands over every record and lets [`ReadFilter`] (Milestone
/// D) decide, which is what keeps the whole filtering policy in one place
/// (spec §2.5).
///
/// BAM-specific because it leans on BAM's self-contained one-record-at-a-time
/// `read_record_buf`. CRAM does not fit that shape (it decodes at *container*
/// granularity into owned records and consults a reference at decode time), so it
/// has its own sibling source, [`CramRecordSource`], below.
pub struct BamRecordSource<R> {
    reader: bam::io::Reader<R>,
    /// BAM record decode ignores the header, but `read_record_buf` requires it.
    header: sam::Header,
    source_file_index: usize,
}

impl<R> BamRecordSource<R> {
    /// Wrap an already-opened BAM reader whose header has already been read.
    pub fn new(reader: bam::io::Reader<R>, header: sam::Header, source_file_index: usize) -> Self {
        Self {
            reader,
            header,
            source_file_index,
        }
    }
}

impl<R: io::Read> RecordSource for BamRecordSource<R> {
    type Record = NoodlesRawRecord;

    fn header(&self) -> &sam::Header {
        &self.header
    }

    fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        // Stamp the file tag before the read; the reader refills `buf.record` in
        // place (reusing its allocations).
        buf.source_file_index = self.source_file_index;
        let bytes_read = self.reader.read_record_buf(&self.header, &mut buf.record)?;
        Ok(bytes_read != 0)
    }
}

/// A [`RecordSource`] over a noodles CRAM reader. Unlike BAM's one-record-at-a-time
/// `read_record_buf`, CRAM is decoded at *container* granularity — a container
/// holds one or more slices whose records are decoded together against the
/// reference — so this source decodes one container into an owned buffer of
/// records and yields them one per [`RecordSource::read_next`], decoding the next
/// container when that buffer drains. `read_next` therefore *moves* an
/// already-decoded [`RecordBuf`] into the caller's buffer rather than refilling
/// it in place (the container buffer's allocations are what get reused across
/// containers). Still **unfiltered**: it hands over every record.
///
/// `reference_sequence_repository` is the FASTA reference CRAM slice decoding
/// consults; it must cover the contigs the CRAM references. This source drives
/// `read_container` + slice decode directly and passes this repository to the
/// slice decoder, so — unlike the reader's inherent `records()`/`query()` methods,
/// which this source does not use — a `cram::io::Reader` built *without* a
/// repository works fine here.
pub struct CramRecordSource<R> {
    reader: cram::io::Reader<R>,
    header: sam::Header,
    reference_sequence_repository: fasta::Repository,
    source_file_index: usize,
    /// The reused container buffer — its allocations are recycled across reads.
    container: cram::io::reader::Container,
    /// The current container's decoded records, yielded one per read.
    buffered_records: std::vec::IntoIter<RecordBuf>,
    /// Set once end of input is reached. noodles' `read_container` is not
    /// idempotent past the CRAM EOF marker (a second call reads past the end and
    /// errors), so this latch makes `read_next` return `Ok(false)` idempotently
    /// after the first EOF — matching the BAM source's behaviour.
    exhausted: bool,
}

/// Outcome of decoding one CRAM container into the record buffer.
enum ContainerRefill {
    /// A container was decoded; its records are now buffered.
    Decoded,
    /// No container left — end of input.
    EndOfInput,
}

impl<R> CramRecordSource<R> {
    /// Wrap an already-opened CRAM reader whose header has already been read.
    /// `reference_sequence_repository` is the FASTA reference used to decode CRAM
    /// slices; it must cover the CRAM's contigs. This repository is authoritative:
    /// the source decodes slices with it directly and never consults any
    /// repository baked into `reader` (the reader's own copy, if any, is used only
    /// by its inherent `records()`/`query()` methods, which this source bypasses).
    pub fn new(
        reader: cram::io::Reader<R>,
        header: sam::Header,
        reference_sequence_repository: fasta::Repository,
        source_file_index: usize,
    ) -> Self {
        Self {
            reader,
            header,
            reference_sequence_repository,
            source_file_index,
            container: cram::io::reader::Container::default(),
            buffered_records: Vec::new().into_iter(),
            exhausted: false,
        }
    }
}

impl<R: io::Read> CramRecordSource<R> {
    /// Decode the next container's records into `self.buffered_records`. This
    /// mirrors noodles' own internal `Records::read_container_records`, but passes
    /// our own repository clone — the reader's copy is private, and only the slice
    /// decoder needs it (we pass it explicitly).
    fn refill_from_next_container(&mut self) -> io::Result<ContainerRefill> {
        if self.reader.read_container(&mut self.container)? == 0 {
            return Ok(ContainerRefill::EndOfInput);
        }
        let compression_header = self.container.compression_header()?;
        let mut records: Vec<RecordBuf> = Vec::new();
        for slice in self.container.slices() {
            let slice = slice?;
            // The decoded block data and the borrowed `Record`s live only within
            // this loop body; converting to owned `RecordBuf` here keeps the
            // buffered result independent of those borrows.
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            for record in slice.records(
                self.reference_sequence_repository.clone(),
                &self.header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )? {
                records.push(RecordBuf::try_from_alignment_record(&self.header, &record)?);
            }
        }
        self.buffered_records = records.into_iter();
        Ok(ContainerRefill::Decoded)
    }
}

impl<R: io::Read> RecordSource for CramRecordSource<R> {
    type Record = NoodlesRawRecord;

    fn header(&self) -> &sam::Header {
        &self.header
    }

    fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        buf.source_file_index = self.source_file_index;
        loop {
            if let Some(record) = self.buffered_records.next() {
                buf.record = record;
                return Ok(true);
            }
            if self.exhausted {
                return Ok(false);
            }
            // Buffer drained → decode the next container (or latch end of input).
            match self.refill_from_next_container()? {
                ContainerRefill::Decoded => {}
                ContainerRefill::EndOfInput => {
                    self.exhausted = true;
                    return Ok(false);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------
// The ReadFilter iterator (the driver)
// ---------------------------------------------------------------------

/// A fatal, run-level error from a [`ReadFilter`] pass. It is yielded **in the
/// iterator's item stream** — a fatal condition (a failed record read, a failed
/// decode, or a reference-fetch failure) makes `next()` return `Some(Err(..))`
/// once and then `None` — so the caller cannot mistake it for a clean end of
/// input: `let read = read?;` propagates it. It is never folded into a per-read
/// drop or a silent EOF.
#[derive(Debug, thiserror::Error)]
pub enum ReadFilterError {
    /// The record source failed to read the next record (e.g. a truncated file).
    #[error("reading the next alignment record failed")]
    Source(#[source] io::Error),
    /// A record that cleared the pre-decode gate failed to decode — a corrupt
    /// record (the unmapped flag clear yet no position).
    #[error("decoding an alignment record failed")]
    Decode(#[source] io::Error),
    /// Filter #8's reference fetch failed. `ReadFilter::new` validates contigs up
    /// front, so this signals corrupt input or a read past a contig end.
    #[error("reference access failed during filtering")]
    Reference(#[source] RefSeqError),
}

/// Filters one sample's reads, lazily, as an
/// `Iterator<Item = Result<MappedRead, ReadFilterError>>`. Each `next()` reads
/// the next record into the single reused buffer, runs the pre-decode cascade
/// (#1–#6) on its flag/MAPQ, decodes only on survival, runs the decode-dependent
/// cascade (#7, #9, #8), tallies every drop, and returns the first read that
/// passes as `Ok(read)`. A read dropped pre-decode builds no `MappedRead` at all.
///
/// **The item is a `Result`** (unlike the spec's original `Item = MappedRead`,
/// revised so a fatal error cannot be silently lost): a fatal condition yields
/// `Some(Err(_))` once and then `None`, so `for read in &mut filter { let read =
/// read?; … }` surfaces it with no separate bookkeeping. This matches the noodles
/// readers this sits on (`records()` is `Item = io::Result<_>`). The iterator is
/// fused — after a clean `None` or a fatal `Err`, further `next()`s return `None`.
///
/// `counts()` is a **running** tally, readable at any point and final once the
/// iterator is exhausted (iterate by `&mut` so the filter — and its counts —
/// outlives the loop).
pub struct ReadFilter<S: RecordSource, R> {
    source: S,
    /// The single record buffer reused across every read.
    record_buf: S::Record,
    /// Raw reference bytes for filter #8; see `ref_seq.md`.
    reference: R,
    config: ReadFilterConfig,
    counts: ReadFilterCounts,
    /// Reused scratch for #8's reference fetch (touched only when #8 runs).
    ref_buf: Vec<u8>,
    /// Set on clean end of input or after a fatal error is yielded; makes the
    /// iterator fused (subsequent `next()`s return `None`).
    done: bool,
}

impl<S: RecordSource, R: RawRefSeq> ReadFilter<S, R> {
    /// Fail-fast setup: validate that every contig in the source header's `@SQ`
    /// list resolves in the reference, then seed the reused record buffer. A
    /// header/reference mismatch is surfaced here (`Err`) rather than mid-stream,
    /// which is what lets an in-loop reference error be treated as fatal-corrupt
    /// (spec §5).
    pub fn new(source: S, reference: R, config: ReadFilterConfig) -> Result<Self, RefSeqError> {
        // A read's `ref_id` indexes the `@SQ` list (ContigId(i) = the i-th `@SQ`),
        // so validating those indices covers every contig a read can reference. A
        // zero-length fetch at position 1 resolves iff the contig exists.
        let contig_count = source.header().reference_sequences().len();
        let mut probe = Vec::new();
        for i in 0..contig_count {
            // PANIC-FREE: a `@SQ` list long enough to overflow `u32` is not
            // representable by any real alignment file (`ref_id` is a 32-bit
            // field); the index fits by construction.
            let contig = ContigId(u32::try_from(i).expect("contig index fits u32"));
            reference.fetch_into(contig, 1, 0, &mut probe)?;
        }
        Ok(Self {
            source,
            record_buf: S::Record::default(),
            reference,
            config,
            counts: ReadFilterCounts::default(),
            ref_buf: Vec::new(),
            done: false,
        })
    }

    /// The running tally — current counts, final once iteration is exhausted.
    pub fn counts(&self) -> &ReadFilterCounts {
        &self.counts
    }

    /// Mark the iterator finished and yield a fatal error. Shared by the three
    /// fatal arms of `next()`.
    fn fail(&mut self, error: ReadFilterError) -> Option<Result<MappedRead, ReadFilterError>> {
        self.done = true;
        Some(Err(error))
    }
}

impl<S: RecordSource, R: RawRefSeq> Iterator for ReadFilter<S, R> {
    type Item = Result<MappedRead, ReadFilterError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Fused: once a clean EOF or a fatal error has been reported, stay stopped.
        if self.done {
            return None;
        }
        loop {
            match self.source.read_next(&mut self.record_buf) {
                Ok(true) => {}
                Ok(false) => {
                    self.done = true;
                    return None; // clean end of input
                }
                Err(error) => return self.fail(ReadFilterError::Source(error)),
            }

            // Phase one — flag/MAPQ, before decode. Exhaustive so a future
            // `FilterVerdict` variant cannot silently fall through to decode.
            match verdict_pre_decode(self.record_buf.flag(), self.record_buf.mapq(), &self.config) {
                FilterVerdict::Keep => {}
                FilterVerdict::Drop(reason) => {
                    self.counts.record_drop(reason);
                    continue;
                }
            }

            // Decode only the pre-decode survivors.
            let read = match self.record_buf.decode() {
                Ok(read) => read,
                Err(error) => return self.fail(ReadFilterError::Decode(error)),
            };

            // Phase two — length / CIGAR / mismatch, on the decoded read.
            match verdict_post_decode(&read, &self.reference, &self.config, &mut self.ref_buf) {
                Ok(FilterVerdict::Keep) => {
                    self.counts.kept += 1;
                    return Some(Ok(read));
                }
                Ok(FilterVerdict::Drop(reason)) => {
                    self.counts.record_drop(reason);
                    continue;
                }
                Err(error) => return self.fail(ReadFilterError::Reference(error)),
            }
        }
    }
}

impl<S: RecordSource, R: RawRefSeq> std::iter::FusedIterator for ReadFilter<S, R> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::FLAG_PAIRED;
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::pileup::walker::CigarOp;

    #[test]
    fn default_config_reproduces_the_production_filter_policy() {
        let config = ReadFilterConfig::default();
        assert_eq!(config.min_mapq, Some(MapQual(DEFAULT_MIN_MAPQ)));
        assert_eq!(config.min_read_length, Some(Bp(DEFAULT_MIN_READ_LENGTH)));
        assert!(config.drop_qc_fail);
        assert!(config.drop_duplicate);
        assert_eq!(
            config.max_read_mismatch_fraction.map(MismatchFraction::get),
            Some(DEFAULT_MAX_READ_MISMATCH_FRACTION)
        );
        assert_eq!(
            config.mismatch_bq_floor,
            BaseQual(DEFAULT_MISMATCH_BQ_FLOOR)
        );
    }

    #[test]
    fn counts_default_is_all_zero() {
        // Explicit all-zero literal (no `..`): pins every counter to 0 and forces
        // this test to be revisited if a counter field is ever added.
        assert_eq!(
            ReadFilterCounts::default(),
            ReadFilterCounts {
                kept: 0,
                duplicate: 0,
                low_mapq: 0,
                supplementary: 0,
                secondary: 0,
                unmapped: 0,
                qc_fail: 0,
                too_short: 0,
                high_mismatch_fraction: 0,
                bad_cigar: 0,
            }
        );
    }

    // ----- verdict_pre_decode (#1–#6) -------------------------------------

    /// A record's flag/MAPQ pair; the pre-decode cascade needs nothing else.
    fn pre(flag: u16, mapq: u8, config: &ReadFilterConfig) -> FilterVerdict {
        verdict_pre_decode(flag, MapQual(mapq), config)
    }

    #[test]
    fn pre_decode_keeps_a_clean_primary_read() {
        let cfg = ReadFilterConfig::default();
        assert_eq!(pre(0, 60, &cfg), FilterVerdict::Keep);
    }

    #[test]
    fn low_mapq_boundary_keeps_at_threshold_drops_one_below() {
        let cfg = ReadFilterConfig::default(); // min_mapq = Some(20)
        assert_eq!(pre(0, 20, &cfg), FilterVerdict::Keep);
        assert_eq!(pre(0, 19, &cfg), FilterVerdict::Drop(DropReason::LowMapq));
        // Unavailable MAPQ arrives as 0 (resolved by the record source) → dropped.
        assert_eq!(pre(0, 0, &cfg), FilterVerdict::Drop(DropReason::LowMapq));
    }

    #[test]
    fn no_mapq_minimum_keeps_any_quality() {
        let cfg = ReadFilterConfig {
            min_mapq: None,
            ..ReadFilterConfig::default()
        };
        assert_eq!(pre(0, 0, &cfg), FilterVerdict::Keep);
    }

    #[test]
    fn each_flag_bit_drops_to_its_own_bucket() {
        let cfg = ReadFilterConfig::default();
        assert_eq!(
            pre(FLAG_DUPLICATE, 60, &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );
        assert_eq!(
            pre(FLAG_SUPPLEMENTARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
        assert_eq!(
            pre(FLAG_SECONDARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Secondary)
        );
        assert_eq!(
            pre(FLAG_UNMAPPED, 60, &cfg),
            FilterVerdict::Drop(DropReason::Unmapped)
        );
        assert_eq!(
            pre(FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::QcFail)
        );
    }

    #[test]
    fn duplicate_and_qc_fail_toggles_off_keep_those_reads() {
        let cfg = ReadFilterConfig {
            drop_duplicate: false,
            drop_qc_fail: false,
            ..ReadFilterConfig::default()
        };
        assert_eq!(pre(FLAG_DUPLICATE, 60, &cfg), FilterVerdict::Keep);
        assert_eq!(pre(FLAG_QC_FAIL, 60, &cfg), FilterVerdict::Keep);
        // Supplementary/secondary/unmapped have no toggle — still dropped.
        assert_eq!(
            pre(FLAG_SUPPLEMENTARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
    }

    #[test]
    fn pre_decode_attribution_charges_the_first_firing_filter() {
        let cfg = ReadFilterConfig::default();
        // Duplicate (1) + unmapped (5) both set → charged to duplicate (earlier).
        assert_eq!(
            pre(FLAG_DUPLICATE | FLAG_UNMAPPED, 60, &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );
        // Supplementary (3) + secondary (4) → supplementary (earlier).
        assert_eq!(
            pre(FLAG_SUPPLEMENTARY | FLAG_SECONDARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
        // Low MAPQ (2) beats unmapped (5) when both apply.
        assert_eq!(
            pre(FLAG_UNMAPPED, 5, &cfg),
            FilterVerdict::Drop(DropReason::LowMapq)
        );
    }

    #[test]
    fn pre_decode_charges_filters_in_full_cascade_order() {
        // Every filter would fire; peel them off one at a time and confirm the
        // read is charged to each in the exact cascade order
        // (duplicate → low-MAPQ → supplementary → secondary → unmapped → QC-fail).
        let cfg = ReadFilterConfig::default();
        let all_flags =
            FLAG_DUPLICATE | FLAG_SUPPLEMENTARY | FLAG_SECONDARY | FLAG_UNMAPPED | FLAG_QC_FAIL;
        // mapq 0 also fails #2 throughout, so every stage after #1 has low-MAPQ
        // waiting behind it — proving each earlier filter really wins.
        assert_eq!(
            pre(all_flags, 0, &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );
        assert_eq!(
            pre(all_flags & !FLAG_DUPLICATE, 0, &cfg),
            FilterVerdict::Drop(DropReason::LowMapq)
        );
        // From here MAPQ is fine (60), so the next flag in order wins.
        assert_eq!(
            pre(all_flags & !FLAG_DUPLICATE, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
        assert_eq!(
            pre(FLAG_SECONDARY | FLAG_UNMAPPED | FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::Secondary)
        );
        assert_eq!(
            pre(FLAG_UNMAPPED | FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::Unmapped)
        );
        assert_eq!(
            pre(FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::QcFail)
        );
    }

    // ----- verdict_post_decode (#7, #9, #8) -------------------------------

    /// Build a mapped read at contig 0, pos 1, with the given decoded sequence,
    /// per-base qualities, and CIGAR. Flag/MAPQ are already-passed values.
    fn mapped(seq: &[u8], qual: &[u8], cigar: Vec<CigarOp>) -> MappedRead {
        MappedRead {
            qname: b"read".to_vec(),
            flag: FLAG_PAIRED,
            ref_id: 0,
            pos: 1,
            mapq: 60,
            cigar,
            seq: seq.to_vec(),
            qual: qual.to_vec(),
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            source_file_index: 0,
        }
    }

    /// A single contig of `n` adenines — the reference for the post-decode
    /// tests, where a `Match` read of `T`s mismatches every aligned base.
    fn poly_a_ref(n: usize) -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![vec![b'A'; n]])
    }

    /// Post-decode config isolating one filter at a time: #7 off, toggles off,
    /// #8 set to `max` (`None` disables it).
    fn post_config(max: Option<f32>) -> ReadFilterConfig {
        ReadFilterConfig {
            min_mapq: None,
            min_read_length: None,
            drop_qc_fail: false,
            drop_duplicate: false,
            max_read_mismatch_fraction: max.map(|x| MismatchFraction::try_new(x).unwrap()),
            mismatch_bq_floor: BaseQual(0),
        }
    }

    fn post(
        read: &MappedRead,
        reference: &impl RawRefSeq,
        config: &ReadFilterConfig,
    ) -> FilterVerdict {
        let mut buf = Vec::new();
        verdict_post_decode(read, reference, config, &mut buf).unwrap()
    }

    #[test]
    fn too_short_boundary_keeps_at_threshold_drops_one_below() {
        let reference = poly_a_ref(40);
        let cfg = ReadFilterConfig {
            min_read_length: Some(Bp(30)),
            ..post_config(None) // #8 already disabled by post_config(None)
        };
        let at = mapped(&vec![b'A'; 30], &vec![40; 30], vec![CigarOp::Match(30)]);
        let below = mapped(&vec![b'A'; 29], &vec![40; 29], vec![CigarOp::Match(29)]);
        assert_eq!(post(&at, &reference, &cfg), FilterVerdict::Keep);
        assert_eq!(
            post(&below, &reference, &cfg),
            FilterVerdict::Drop(DropReason::TooShort)
        );
    }

    #[test]
    fn bad_cigar_drops_the_two_ill_formed_shapes() {
        let reference = poly_a_ref(40);
        let cfg = post_config(None);
        // Adjacent insertion/deletion pair.
        let adjacent_indel = mapped(
            b"AAAAAAAA",
            &[40; 8],
            vec![
                CigarOp::Match(4),
                CigarOp::Insertion(1),
                CigarOp::Deletion(1),
                CigarOp::Match(4),
            ],
        );
        // Leading deletion (boundary deletion).
        let boundary_deletion = mapped(
            b"AAAAAAAA",
            &[40; 8],
            vec![CigarOp::Deletion(2), CigarOp::Match(8)],
        );
        assert_eq!(
            post(&adjacent_indel, &reference, &cfg),
            FilterVerdict::Drop(DropReason::BadCigar)
        );
        assert_eq!(
            post(&boundary_deletion, &reference, &cfg),
            FilterVerdict::Drop(DropReason::BadCigar)
        );
    }

    #[test]
    fn high_mismatch_boundary_keeps_at_threshold_drops_above() {
        // ref = 10×A; a Match(10) read with k mismatches has fraction k/10.
        let reference = poly_a_ref(10);
        let cfg = post_config(Some(0.10));
        // 1/10 = 0.10, not > 0.10 → kept (boundary is exclusive).
        let at = mapped(b"TAAAAAAAAA", &[40; 10], vec![CigarOp::Match(10)]);
        // 2/10 = 0.20 > 0.10 → dropped.
        let above = mapped(b"TTAAAAAAAA", &[40; 10], vec![CigarOp::Match(10)]);
        assert_eq!(post(&at, &reference, &cfg), FilterVerdict::Keep);
        assert_eq!(
            post(&above, &reference, &cfg),
            FilterVerdict::Drop(DropReason::HighMismatchFraction)
        );
    }

    #[test]
    fn low_quality_mismatches_do_not_count_toward_the_fraction() {
        let reference = poly_a_ref(10);
        let cfg = ReadFilterConfig {
            mismatch_bq_floor: BaseQual(10),
            ..post_config(Some(0.0))
        };
        // Two mismatches, both below the BQ floor → neither counts → kept even
        // at a zero threshold.
        let read = mapped(
            b"TTAAAAAAAA",
            &[5, 5, 40, 40, 40, 40, 40, 40, 40, 40],
            vec![CigarOp::Match(10)],
        );
        assert_eq!(post(&read, &reference, &cfg), FilterVerdict::Keep);
    }

    #[test]
    fn mismatch_filter_disabled_makes_no_reference_access() {
        // Empty reference: any fetch errors with UnknownContig. With #8 disabled
        // the fetch never happens, so a high-mismatch read is kept.
        let empty = InMemoryRefSeq::from_contigs(Vec::new());
        let read = mapped(b"TTTTTTTTTT", &[40; 10], vec![CigarOp::Match(10)]);

        let disabled = post_config(None);
        assert_eq!(post(&read, &empty, &disabled), FilterVerdict::Keep);

        // Enabling #8 on the same empty reference proves a fetch is attempted.
        let enabled = post_config(Some(0.10));
        let mut buf = Vec::new();
        assert!(matches!(
            verdict_post_decode(&read, &empty, &enabled, &mut buf),
            Err(RefSeqError::UnknownContig(_))
        ));
    }

    #[test]
    fn high_mismatch_fetch_past_contig_end_is_fatal() {
        // A read whose reference span runs past the contig end (Match(10) at
        // pos 1 on a 5-base contig). Under the fatal error model the OutOfBounds
        // fetch propagates as `Err`, not a per-read drop or keep — a
        // validly-aligned read cannot cover positions the contig lacks, so this
        // signals a malformed record.
        let reference = poly_a_ref(5);
        let read = mapped(b"TTTTTTTTTT", &[40; 10], vec![CigarOp::Match(10)]);
        let mut buf = Vec::new();
        assert!(matches!(
            verdict_post_decode(&read, &reference, &post_config(Some(0.10)), &mut buf),
            Err(RefSeqError::OutOfBounds { .. })
        ));
    }

    #[test]
    fn bad_cigar_is_charged_before_high_mismatch() {
        // A read that is BOTH a boundary deletion (#9) and, on its M positions,
        // fully mismatched against the reference (#8 at a zero threshold).
        let reference = poly_a_ref(10);
        let both = mapped(
            b"TTTTTTTT",
            &[40; 8],
            vec![CigarOp::Deletion(2), CigarOp::Match(8)],
        );
        // #9 fires first → charged BadCigar, and #8's reference walk is skipped.
        assert_eq!(
            post(&both, &reference, &post_config(Some(0.0))),
            FilterVerdict::Drop(DropReason::BadCigar)
        );
        // The same sequence with a well-formed CIGAR does reach #8 and drops there,
        // proving the read would have failed #8 too — attribution is the only
        // difference the ordering makes.
        let good_cigar = mapped(b"TTTTTTTT", &[40; 8], vec![CigarOp::Match(8)]);
        assert_eq!(
            post(&good_cigar, &reference, &post_config(Some(0.0))),
            FilterVerdict::Drop(DropReason::HighMismatchFraction)
        );
    }

    #[test]
    fn too_short_is_charged_before_bad_cigar() {
        let reference = poly_a_ref(10);
        let cfg = ReadFilterConfig {
            min_read_length: Some(Bp(30)),
            ..post_config(None)
        };
        // Short AND a boundary deletion → charged TooShort (#7 before #9).
        let read = mapped(
            b"AAAAA",
            &[40; 5],
            vec![CigarOp::Deletion(2), CigarOp::Match(5)],
        );
        assert_eq!(
            post(&read, &reference, &cfg),
            FilterVerdict::Drop(DropReason::TooShort)
        );
    }

    #[test]
    fn zero_reference_span_read_is_kept_by_the_mismatch_filter() {
        // An all-soft-clip read has ref_span 0 → empty slice → no comparable
        // bases → kept (matches production's skip-on-no-span behaviour).
        let reference = poly_a_ref(10);
        let read = mapped(b"TTTT", &[40; 4], vec![CigarOp::SoftClip(4)]);
        assert_eq!(
            post(&read, &reference, &post_config(Some(0.0))),
            FilterVerdict::Keep
        );
    }

    // ----- the record-source seam (C) -------------------------------------

    use noodles_core::Position;
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record::{Flags, MappingQuality};
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::ReferenceSequence;
    use std::num::NonZero;

    /// A trivial in-memory `RawRecord`: it carries a whole `MappedRead` and reads
    /// its `flag`/`mapq` back off it, so the cascade can be driven with no BAM.
    /// `decode_fails` makes `decode` return an error (to exercise the fatal
    /// decode-error path).
    #[derive(Clone)]
    struct FakeRecord {
        read: MappedRead,
        decode_fails: bool,
    }

    impl Default for FakeRecord {
        fn default() -> Self {
            Self {
                read: mapped(b"", &[], Vec::new()),
                decode_fails: false,
            }
        }
    }

    impl RawRecord for FakeRecord {
        fn flag(&self) -> u16 {
            self.read.flag
        }
        fn mapq(&self) -> MapQual {
            MapQual(self.read.mapq)
        }
        fn decode(&self) -> io::Result<MappedRead> {
            if self.decode_fails {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "fake decode failure",
                ));
            }
            Ok(self.read.clone())
        }
    }

    /// A `RecordSource` yielding a fixed slice of `FakeRecord`s, refilling the
    /// caller's buffer each call.
    struct FakeSource {
        records: Vec<FakeRecord>,
        header: sam::Header,
        next_index: usize,
    }

    impl FakeSource {
        fn new(records: Vec<FakeRecord>, header: sam::Header) -> Self {
            Self {
                records,
                header,
                next_index: 0,
            }
        }
    }

    impl RecordSource for FakeSource {
        type Record = FakeRecord;
        fn header(&self) -> &sam::Header {
            &self.header
        }
        fn read_next(&mut self, buf: &mut FakeRecord) -> io::Result<bool> {
            match self.records.get(self.next_index) {
                Some(record) => {
                    *buf = record.clone();
                    self.next_index += 1;
                    Ok(true)
                }
                None => Ok(false),
            }
        }
    }

    /// A 30 bp all-`A` mapped read carrying the given flag/MAPQ — long enough to
    /// clear the default `min_read_length`, and matching a poly-A reference.
    fn fake(flag: u16, mapq: u8) -> FakeRecord {
        let mut read = mapped(&vec![b'A'; 30], &vec![30u8; 30], vec![CigarOp::Match(30)]);
        read.flag = flag;
        read.mapq = mapq;
        FakeRecord {
            read,
            decode_fails: false,
        }
    }

    #[test]
    fn fake_source_drives_the_seam() {
        // Mimic what the Milestone-D iterator will do: read_next → pre-decode →
        // decode survivors → post-decode, over the fake source.
        let reference = poly_a_ref(30);
        let cfg = ReadFilterConfig::default();
        let mut source = FakeSource::new(
            vec![fake(FLAG_DUPLICATE, 60), fake(FLAG_PAIRED, 60)],
            one_contig_header(),
        );
        let mut buf = FakeRecord::default();
        let mut ref_buf = Vec::new();

        // Record 1 — duplicate: dropped pre-decode, never decoded.
        assert!(source.read_next(&mut buf).unwrap());
        assert_eq!(buf.flag(), FLAG_DUPLICATE);
        assert_eq!(buf.mapq(), MapQual(60));
        assert_eq!(
            verdict_pre_decode(buf.flag(), buf.mapq(), &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );

        // Record 2 — clean: pre-decode keep → decode → post-decode keep.
        assert!(source.read_next(&mut buf).unwrap());
        assert_eq!(
            verdict_pre_decode(buf.flag(), buf.mapq(), &cfg),
            FilterVerdict::Keep
        );
        let read = buf.decode().unwrap();
        assert_eq!(
            verdict_post_decode(&read, &reference, &cfg, &mut ref_buf).unwrap(),
            FilterVerdict::Keep
        );

        // End of input.
        assert!(!source.read_next(&mut buf).unwrap());
    }

    fn bam_record(
        qname: &str,
        ref_id: usize,
        start: usize,
        len: usize,
        mapq: u8,
        flags: Flags,
    ) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname)
            .set_reference_sequence_id(ref_id)
            .set_flags(flags)
            .set_mapping_quality(MappingQuality::new(mapq).expect("mapq in range"))
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, len)].into_iter().collect())
            .set_sequence(Sequence::from(vec![b'A'; len]))
            .set_quality_scores(QualityScores::from(vec![30u8; len]))
            .build()
    }

    /// A header with a single `chr1` reference sequence of the given length.
    fn contig_header(length: usize) -> sam::Header {
        sam::Header::builder()
            .set_header(Default::default())
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(length).unwrap()),
            )
            .build()
    }

    fn one_contig_header() -> sam::Header {
        contig_header(100)
    }

    /// Encode a header + records to an in-memory BAM (BGZF), returning the bytes.
    fn in_memory_bam(header: &sam::Header, records: &[RecordBuf]) -> Vec<u8> {
        let mut bytes = Vec::new();
        // Scope the writer so its `&mut bytes` borrow (held live by its `Drop`)
        // is released before `bytes` is returned.
        {
            let mut writer = bam::io::Writer::new(&mut bytes);
            writer.write_header(header).expect("write header");
            for record in records {
                writer
                    .write_alignment_record(header, record)
                    .expect("write record");
            }
            writer.try_finish().expect("finish bam");
        }
        bytes
    }

    #[test]
    fn noodles_raw_record_reads_flag_mapq_and_decodes() {
        let record = bam_record("r1", 0, 5, 4, 42, Flags::from(FLAG_DUPLICATE));
        let raw = NoodlesRawRecord {
            record,
            source_file_index: 7,
        };
        // Cheap pre-decode reads.
        assert_eq!(raw.flag(), FLAG_DUPLICATE);
        assert_eq!(raw.mapq(), MapQual(42));
        // Decode produces the MappedRead the cascade's phase two consumes.
        let read = raw.decode().unwrap();
        assert_eq!(read.ref_id, 0);
        assert_eq!(read.pos, 5);
        assert_eq!(read.mapq, 42);
        assert_eq!(read.seq, b"AAAA");
        assert_eq!(read.source_file_index, 7);
    }

    #[test]
    fn noodles_raw_record_maps_unavailable_mapq_to_zero() {
        // A record with no mapping quality (SAM 0xFF) → MapQual(0).
        let record = RecordBuf::builder()
            .set_reference_sequence_id(0)
            .set_flags(Flags::from(FLAG_PAIRED))
            .set_alignment_start(Position::try_from(1usize).unwrap())
            .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 4]))
            .build();
        let raw = NoodlesRawRecord {
            record,
            source_file_index: 0,
        };
        assert_eq!(raw.mapq(), MapQual(0));
    }

    #[test]
    fn noodles_raw_record_decode_errors_on_a_record_with_no_position() {
        // A default record has no reference_sequence_id / alignment_start, so the
        // reused decoder fails — the Err a corrupt record would surface (fatal).
        let raw = NoodlesRawRecord::default();
        let err = raw.decode().expect_err("missing position must fail decode");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn bam_record_source_reads_flag_mapq_and_decodes_through_a_real_bam() {
        let header = one_contig_header();
        let records = [
            bam_record("r1", 0, 10, 30, 60, Flags::from(FLAG_PAIRED)),
            bam_record(
                "r2",
                0,
                20,
                30,
                40,
                Flags::from(FLAG_PAIRED | FLAG_DUPLICATE),
            ),
        ];
        let bytes = in_memory_bam(&header, &records);

        let mut reader = bam::io::Reader::new(&bytes[..]);
        let read_header = reader.read_header().expect("read header");
        let mut source = BamRecordSource::new(reader, read_header, 3);

        let mut buf = NoodlesRawRecord::default();

        // Record 1 — pre-decode reads, then decode.
        assert!(source.read_next(&mut buf).unwrap());
        assert_eq!(buf.flag(), FLAG_PAIRED);
        assert_eq!(buf.mapq(), MapQual(60));
        let read = buf.decode().unwrap();
        assert_eq!(read.pos, 10);
        assert_eq!(read.seq.len(), 30);
        assert_eq!(read.source_file_index, 3);

        // Record 2 — the reused buffer is refilled in place.
        assert!(source.read_next(&mut buf).unwrap());
        assert_eq!(buf.flag(), FLAG_PAIRED | FLAG_DUPLICATE);
        assert_eq!(buf.mapq(), MapQual(40));

        // End of input.
        assert!(!source.read_next(&mut buf).unwrap());
    }

    #[test]
    fn bam_record_source_reuses_the_buffer_without_leaking_a_prior_record() {
        // A long record then a short one, both decoded through the SAME reused
        // buffer: if `read_record_buf` did not fully overwrite the buffer, the
        // 40-base record's tail would leak into the 10-base one.
        let header = one_contig_header();
        let records = [
            bam_record("long", 0, 10, 40, 60, Flags::from(FLAG_PAIRED)),
            bam_record("short", 0, 60, 10, 30, Flags::from(FLAG_PAIRED)),
        ];
        let bytes = in_memory_bam(&header, &records);

        let mut reader = bam::io::Reader::new(&bytes[..]);
        let read_header = reader.read_header().expect("read header");
        let mut source = BamRecordSource::new(reader, read_header, 0);
        let mut buf = NoodlesRawRecord::default();

        assert!(source.read_next(&mut buf).unwrap());
        assert_eq!(buf.decode().unwrap().seq.len(), 40);

        assert!(source.read_next(&mut buf).unwrap());
        let short = buf.decode().unwrap();
        assert_eq!(short.seq.len(), 10, "no stale tail from the 40-base record");
        assert_eq!(short.pos, 60);
        assert_eq!(short.mapq, 30);
    }

    /// A mapped record with an explicit sequence (so a test can put a
    /// non-reference base in it and check the reconstructed bytes).
    fn record_with_seq(qname: &str, start: usize, mapq: u8, flags: Flags, seq: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname)
            .set_reference_sequence_id(0usize)
            .set_flags(flags)
            .set_mapping_quality(MappingQuality::new(mapq).expect("mapq in range"))
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, seq.len())].into_iter().collect())
            .set_sequence(Sequence::from(seq.to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; seq.len()]))
            .build()
    }

    #[test]
    fn cram_record_source_reads_flag_mapq_and_decodes_through_a_real_cram() {
        use crate::bam::alignment_input::build_fasta_repository;
        use crate::pileup::per_sample::cram_files::{
            ContigSpec, HeaderOverrides, build_cram, build_fasta,
        };

        let contigs = vec![ContigSpec {
            name: "chr1".into(),
            length: 100,
        }];
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).unwrap();
        // The reference is all 'A'. Record 1 carries non-reference bases (`C`,
        // `G`), so CRAM must store and decode back substitution features —
        // exercising reference-diff reconstruction, not just a length.
        let r1_seq = b"AACAAAAAAG";
        let records = [
            record_with_seq("r1", 10, 60, Flags::from(FLAG_PAIRED), r1_seq),
            record_with_seq(
                "r2",
                40,
                40,
                Flags::from(FLAG_PAIRED | FLAG_DUPLICATE),
                b"AAAAA",
            ),
        ];
        let (_cram_dir, cram_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records).unwrap();
        let repository = build_fasta_repository(&fasta_path).unwrap();

        let file = std::fs::File::open(&cram_path).unwrap();
        let mut reader = cram::io::Reader::new(file);
        let header = reader.read_header().unwrap();
        let mut source = CramRecordSource::new(reader, header, repository, 5);
        let mut buf = NoodlesRawRecord::default();

        // Record 1 — pre-decode reads through the container-buffered source, then decode.
        assert!(source.read_next(&mut buf).unwrap());
        assert_eq!(buf.flag(), FLAG_PAIRED);
        assert_eq!(buf.mapq(), MapQual(60));
        let read = buf.decode().unwrap();
        assert_eq!(read.pos, 10);
        assert_eq!(
            read.seq, r1_seq,
            "the non-reference bases decode back exactly"
        );
        assert_eq!(read.source_file_index, 5);

        // Record 2.
        assert!(source.read_next(&mut buf).unwrap());
        assert_eq!(buf.flag(), FLAG_PAIRED | FLAG_DUPLICATE);
        assert_eq!(buf.mapq(), MapQual(40));

        // End of input (buffer drained + no further container).
        assert!(!source.read_next(&mut buf).unwrap());

        // `_fasta_dir`/`_cram_dir` keep the tempdirs alive through the read above.
    }

    #[test]
    fn cram_record_source_matches_noodles_records_across_multiple_containers() {
        use crate::bam::alignment_input::build_fasta_repository;
        use crate::pileup::per_sample::cram_files::{
            ContigSpec, HeaderOverrides, build_cram, build_fasta,
        };

        // noodles packs 10240 records per CRAM container, so > 10240 records force
        // at least two containers — exercising the drain-and-refill loop that is
        // the whole reason CramRecordSource differs from BamRecordSource.
        let record_count = 10_241usize;
        let contigs = vec![ContigSpec {
            name: "chr1".into(),
            length: 200,
        }];
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).unwrap();
        // Coordinate-sorted (non-decreasing positions), all within the contig.
        let records: Vec<RecordBuf> = (0..record_count)
            .map(|i| {
                let start = 1 + i / 100; // 1..=103, monotonic; footprint <= 113
                record_with_seq(
                    &format!("r{i}"),
                    start,
                    40,
                    Flags::from(FLAG_PAIRED),
                    b"AAAAAAAAAA",
                )
            })
            .collect();
        let (_cram_dir, cram_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records).unwrap();
        let repository = build_fasta_repository(&fasta_path).unwrap();

        // Ground truth: noodles' own whole-file record iterator (reader built WITH
        // the repository, which its inherent `records()` requires).
        let expected: Vec<RecordBuf> = {
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .build_from_path(&cram_path)
                .unwrap();
            let header = reader.read_header().unwrap();
            reader
                .records(&header)
                .collect::<io::Result<Vec<_>>>()
                .unwrap()
        };

        // Ours: a plain reader (no baked-in repository), repo passed explicitly.
        let mut reader = cram::io::Reader::new(std::fs::File::open(&cram_path).unwrap());
        let header = reader.read_header().unwrap();
        let mut source = CramRecordSource::new(reader, header, repository, 0);
        let mut buf = NoodlesRawRecord::default();
        let mut actual: Vec<RecordBuf> = Vec::new();
        while source.read_next(&mut buf).unwrap() {
            actual.push(buf.record.clone());
        }
        // Post-EOF is idempotently Ok(false).
        assert!(!source.read_next(&mut buf).unwrap());

        assert_eq!(
            actual.len(),
            record_count,
            "every record across both containers is yielded"
        );
        assert_eq!(
            actual, expected,
            "identical to noodles' own record iterator, in order"
        );
    }

    // ----- the ReadFilter iterator (D) ------------------------------------

    /// A read whose CIGAR is `Deletion(2), Match(seq.len())` (a boundary deletion
    /// → bad CIGAR), long enough to clear `min_read_length` so #9 — not #7 —
    /// fires. `seq` lets a test also make it high-mismatch (fails #8) to check the
    /// #9-before-#8 attribution.
    fn boundary_deletion_record(start: usize, seq: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name("bad")
            .set_reference_sequence_id(0usize)
            .set_flags(Flags::from(FLAG_PAIRED))
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(
                [Op::new(Kind::Deletion, 2), Op::new(Kind::Match, seq.len())]
                    .into_iter()
                    .collect(),
            )
            .set_sequence(Sequence::from(seq.to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; seq.len()]))
            .build()
    }

    /// The fixture: two kept reads plus one read per *mapped* drop reason (#1–#4,
    /// #6–#9), all on a single all-`A` contig. The #5 unmapped drop is covered
    /// separately (`read_filter_charges_an_unmapped_read_end_to_end`): a realistic
    /// unmapped read has MAPQ 0 and is charged to #2 first, and a fake MAPQ does
    /// not survive a CRAM round-trip. Returns the records and the
    /// `ReadFilterCounts` a correct pass must produce; the counts are asserted, not
    /// read order.
    fn drop_fixture() -> (Vec<RecordBuf>, ReadFilterCounts) {
        let clean = |name: &str, start: usize| {
            record_with_seq(
                name,
                start,
                60,
                Flags::from(FLAG_PAIRED),
                b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            )
        };
        let records = vec![
            clean("kept1", 10), // kept
            clean("kept2", 20), // kept
            record_with_seq(
                "dup",
                30,
                60,
                Flags::from(FLAG_PAIRED | FLAG_DUPLICATE),
                b"AAAAAAAAAA",
            ), // #1
            record_with_seq("lowmapq", 40, 5, Flags::from(FLAG_PAIRED), b"AAAAAAAAAA"), // #2 (mapq 5 < 20)
            record_with_seq(
                "supp",
                50,
                60,
                Flags::from(FLAG_PAIRED | FLAG_SUPPLEMENTARY),
                b"AAAAAAAAAA",
            ), // #3
            record_with_seq(
                "sec",
                60,
                60,
                Flags::from(FLAG_PAIRED | FLAG_SECONDARY),
                b"AAAAAAAAAA",
            ), // #4
            record_with_seq(
                "qcfail",
                70,
                60,
                Flags::from(FLAG_PAIRED | FLAG_QC_FAIL),
                b"AAAAAAAAAA",
            ), // #6
            record_with_seq("tooshort", 80, 60, Flags::from(FLAG_PAIRED), b"AAAAA"), // #7 (len 5 < 30)
            // #8: 5 non-reference bases out of 30 = 16.7% > 10%.
            record_with_seq(
                "highmismatch",
                90,
                60,
                Flags::from(FLAG_PAIRED),
                b"CCCCCAAAAAAAAAAAAAAAAAAAAAAAAA",
            ),
            // #9 bad CIGAR — AND high-mismatch (5 `C`s). Because it fails both #9
            // and #8, the exact counts below discriminate the ng #9-before-#8
            // order: it must land in `bad_cigar` (not `high_mismatch_fraction`).
            boundary_deletion_record(120, b"CCCCCAAAAAAAAAAAAAAAAAAAAAAAAA"), // #9
        ];
        let expected = ReadFilterCounts {
            kept: 2,
            duplicate: 1,
            low_mapq: 1,
            supplementary: 1,
            secondary: 1,
            unmapped: 0,
            qc_fail: 1,
            too_short: 1,
            high_mismatch_fraction: 1,
            bad_cigar: 1,
        };
        (records, expected)
    }

    /// The `@SQ`-length-200 all-`A` reference the drop fixture filters against.
    fn fixture_reference() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![vec![b'A'; 200]])
    }

    #[test]
    fn read_filter_bam_fixture_matches_hand_counted_drops() {
        let (records, expected) = drop_fixture();
        let header = one_contig_200_header();
        let bytes = in_memory_bam(&header, &records);

        let mut reader = bam::io::Reader::new(&bytes[..]);
        let read_header = reader.read_header().unwrap();
        let source = BamRecordSource::new(reader, read_header, 0);
        let mut filter =
            ReadFilter::new(source, fixture_reference(), ReadFilterConfig::default()).unwrap();

        let kept: Vec<MappedRead> = (&mut filter).collect::<Result<_, _>>().unwrap();
        assert_eq!(kept.len(), 2, "exactly the two clean reads survive");
        assert_eq!(*filter.counts(), expected);
    }

    #[test]
    fn read_filter_cram_fixture_matches_hand_counted_drops() {
        use crate::bam::alignment_input::build_fasta_repository;
        use crate::pileup::per_sample::cram_files::{
            ContigSpec, HeaderOverrides, build_cram, build_fasta,
        };

        let (records, expected) = drop_fixture();
        let contigs = vec![ContigSpec {
            name: "chr1".into(),
            length: 200,
        }];
        let (_fasta_dir, fasta_path) = build_fasta(&contigs).unwrap();
        let (_cram_dir, cram_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), &records).unwrap();
        let repository = build_fasta_repository(&fasta_path).unwrap();

        let mut reader = cram::io::Reader::new(std::fs::File::open(&cram_path).unwrap());
        let header = reader.read_header().unwrap();
        let source = CramRecordSource::new(reader, header, repository, 0);
        let mut filter =
            ReadFilter::new(source, fixture_reference(), ReadFilterConfig::default()).unwrap();

        let kept: Vec<MappedRead> = (&mut filter).collect::<Result<_, _>>().unwrap();
        assert_eq!(kept.len(), 2);
        assert_eq!(*filter.counts(), expected);
    }

    #[test]
    fn read_filter_new_rejects_a_contig_missing_from_the_reference() {
        // Header declares two contigs; the reference has only one.
        let header = sam::Header::builder()
            .set_header(Default::default())
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(100usize).unwrap()),
            )
            .add_reference_sequence(
                "chr2",
                Map::<ReferenceSequence>::new(NonZero::new(100usize).unwrap()),
            )
            .build();
        let source = FakeSource::new(Vec::new(), header);
        let result = ReadFilter::new(source, poly_a_ref(100), ReadFilterConfig::default());
        assert!(matches!(
            result,
            Err(RefSeqError::UnknownContig(ContigId(1)))
        ));
    }

    /// A `RecordSource` whose `read_next` always fails — to drive the fatal path.
    struct ErroringSource {
        header: sam::Header,
    }

    impl RecordSource for ErroringSource {
        type Record = FakeRecord;
        fn header(&self) -> &sam::Header {
            &self.header
        }
        fn read_next(&mut self, _buf: &mut FakeRecord) -> io::Result<bool> {
            Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated"))
        }
    }

    #[test]
    fn read_filter_source_read_error_is_fatal() {
        let source = ErroringSource {
            header: one_contig_header(),
        };
        let mut filter =
            ReadFilter::new(source, poly_a_ref(100), ReadFilterConfig::default()).unwrap();
        // The read error surfaces as the yielded item …
        assert!(matches!(
            (&mut filter).next(),
            Some(Err(ReadFilterError::Source(_)))
        ));
        // … and the iterator is fused (stays stopped) afterward.
        assert!((&mut filter).next().is_none());
    }

    #[test]
    fn read_filter_decode_error_is_fatal() {
        // A record that clears the pre-decode gate but fails to decode.
        let mut record = fake(FLAG_PAIRED, 60);
        record.decode_fails = true;
        let source = FakeSource::new(vec![record], one_contig_header());
        let mut filter =
            ReadFilter::new(source, poly_a_ref(100), ReadFilterConfig::default()).unwrap();
        assert!(matches!(
            (&mut filter).next(),
            Some(Err(ReadFilterError::Decode(_)))
        ));
        assert!((&mut filter).next().is_none());
    }

    #[test]
    fn read_filter_reference_error_mid_stream_is_fatal() {
        // A read whose #8 reference window runs past the (short) contig end. The
        // contig resolves in `new` (a zero-length probe), so this fails only
        // in-loop → fatal Reference error.
        let read = mapped(&vec![b'A'; 50], &vec![30u8; 50], vec![CigarOp::Match(50)]);
        let source = FakeSource::new(
            vec![FakeRecord {
                read,
                decode_fails: false,
            }],
            one_contig_header(),
        );
        // contig length 10 < the read's 50-base reference span.
        let mut filter =
            ReadFilter::new(source, poly_a_ref(10), ReadFilterConfig::default()).unwrap();
        assert!(matches!(
            (&mut filter).next(),
            Some(Err(ReadFilterError::Reference(
                RefSeqError::OutOfBounds { .. }
            )))
        ));
        assert!((&mut filter).next().is_none());
    }

    #[test]
    fn read_filter_charges_an_unmapped_read_end_to_end() {
        // The #5 counter the BAM/CRAM fixture omits: an unmapped read that clears
        // #2 (MAPQ 60) is charged to `unmapped` through the full iterator.
        let unmapped = fake(FLAG_PAIRED | FLAG_UNMAPPED, 60);
        let clean = fake(FLAG_PAIRED, 60);
        let source = FakeSource::new(vec![unmapped, clean], one_contig_header());
        let mut filter =
            ReadFilter::new(source, poly_a_ref(30), ReadFilterConfig::default()).unwrap();
        let kept: Vec<MappedRead> = (&mut filter).collect::<Result<_, _>>().unwrap();
        assert_eq!(kept.len(), 1, "only the clean read survives");
        assert_eq!(filter.counts().unmapped, 1);
        assert_eq!(filter.counts().kept, 1);
    }

    #[test]
    fn read_filter_over_an_empty_source_yields_nothing_and_zero_counts() {
        let source = FakeSource::new(Vec::new(), one_contig_header());
        let mut filter =
            ReadFilter::new(source, poly_a_ref(30), ReadFilterConfig::default()).unwrap();
        assert!((&mut filter).next().is_none());
        assert_eq!(*filter.counts(), ReadFilterCounts::default());
    }

    #[test]
    fn read_filter_counts_is_a_running_tally_before_exhaustion() {
        // A duplicate (dropped) then a clean read (kept): the first `next()`
        // returns the kept read, and `counts()` already reflects both.
        let source = FakeSource::new(
            vec![
                fake(FLAG_PAIRED | FLAG_DUPLICATE, 60),
                fake(FLAG_PAIRED, 60),
            ],
            one_contig_header(),
        );
        let mut filter =
            ReadFilter::new(source, poly_a_ref(30), ReadFilterConfig::default()).unwrap();
        assert!(matches!((&mut filter).next(), Some(Ok(_))));
        assert_eq!(filter.counts().duplicate, 1);
        assert_eq!(filter.counts().kept, 1);
        assert!((&mut filter).next().is_none());
        assert_eq!(filter.counts().kept, 1);
    }

    fn one_contig_200_header() -> sam::Header {
        contig_header(200)
    }
}
