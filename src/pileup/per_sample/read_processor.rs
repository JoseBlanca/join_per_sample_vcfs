//! The full per-read processing fold: one `MappedRead` in, one
//! `PreparedRead` (or a structured drop) out. This is the unit of
//! work the parallel read-processing stage runs, and it bundles
//! everything that happens to a read between the coordinate merge and
//! the walker:
//!
//! 1. **G2** — `cigar_is_bad` rejection (adjacent I/D, boundary
//!    deletion). Runs first, *before* left-alignment, so the check
//!    sees the aligner's own CIGAR.
//! 2. **F3** — `left_align_indels` canonicalisation (mutates the
//!    CIGAR in place). Needs the read's reference slice.
//! 3. **F1** — `read_exceeds_mismatch_fraction` rejection, computed
//!    against the left-aligned CIGAR.
//! 4. **BAQ** — `BaqEngine::process` caps each base quality and builds
//!    the `PreparedRead`.
//!
//! These steps used to be split: G2/F3/F1 ran serially inside
//! `AlignmentMergedReader::next`, only BAQ ran in parallel. Folding
//! them into one function lets the *whole* per-read cost run on the
//! worker threads, shrinking the serial floor to the coordinate
//! merge + walker. Every step here is a pure function of the read and
//! its reference window, so the fold is deterministic and order-
//! independent across reads — the property the parallel stage relies
//! on for byte-identical output.
//!
//! **Two reference sources, on purpose.** F3/F1 read *raw* reference
//! bytes (case-preserving, no `N` canonicalisation) from
//! [`RawContigRefCache`]; BAQ reads *uppercased* windows from its own
//! [`ManualEvictChromRefFetcher`]. This mirrors exactly what the old
//! serial path did — the reader's `fetch_ref_for_read` returned raw
//! bytes for F3/F1 while BAQ used the streaming fetcher — and the
//! distinction is load-bearing for byte-identity on soft-masked
//! (lowercase) references, where `left_align_indels`' case-sensitive
//! comparison would shift indels differently against uppercased bytes.

use std::sync::Arc;

use noodles_fasta as fasta;

use crate::bam::alignment_input::{
    MappedRead, cigar_is_bad, cigar_ref_span, read_exceeds_mismatch_fraction,
};
use crate::fasta::{ContigList, ManualEvictChromRefFetcher};
use crate::pileup::per_sample::baq_engine::{
    BaqEngine, BaqOutcome, BaqSkipReason, prepare_passthrough,
};
use crate::pileup::walker::{PreparedRead, indel_norm::left_align_indels};

/// Knobs for the reference-dependent read filters (F1). The F3
/// left-alignment and G2 CIGAR check take no configuration. BAQ is
/// configured separately on the [`BaqEngine`] itself.
#[derive(Debug, Clone, Copy)]
pub struct ReadProcessingConfig {
    /// F1 mismatch-fraction threshold; `None` disables the filter
    /// (mirrors `AlignmentMergedReaderConfig::max_read_mismatch_fraction`).
    pub max_read_mismatch_fraction: Option<f32>,
    /// BQ floor below which a mismatch does not count toward F1.
    pub mismatch_bq_floor: u8,
}

/// Per-worker cache of one resident contig's *raw* reference bytes,
/// used by F3/F1. Lifted verbatim from the old
/// `AlignmentMergedReader::fetch_ref_for_read`: a whole-contig
/// sequence pulled from the shared `Repository`, sliced raw (no
/// uppercasing, no `N` canonicalisation), cached by `ref_id` and
/// refreshed on a contig change.
///
/// The `Repository` is `Arc`-backed, so cloning one into each worker
/// is cheap and every worker shares the same resident contig bytes
/// (the cache only holds an `Arc` clone of the sequence).
pub struct RawContigRefCache {
    repository: fasta::Repository,
    contigs: ContigList,
    cached_contig: Option<(usize, Arc<fasta::record::Sequence>)>,
}

impl RawContigRefCache {
    pub fn new(repository: fasta::Repository, contigs: ContigList) -> Self {
        Self {
            repository,
            contigs,
            cached_contig: None,
        }
    }

    /// Raw reference bytes for the 1-based, half-open span
    /// `[pos, pos+ref_span)`, or `None` when there is nothing to fetch
    /// or the window is unusable — exactly the cases the old
    /// `fetch_ref_for_read` returned `None` for: zero span, unknown
    /// `ref_id`, repository miss, `pos == 0`, or a window past the
    /// contig end. Returns a borrowed slice into the cached sequence
    /// (the old code allocated a `Vec`; the bytes are identical).
    pub fn fetch_raw_slice(
        &mut self,
        ref_id: usize,
        pos_1based: u64,
        ref_span: u32,
    ) -> Option<&[u8]> {
        if ref_span == 0 {
            return None;
        }
        let contig = self.contigs.entries.get(ref_id)?;

        let needs_refresh = match &self.cached_contig {
            Some((id, _)) => *id != ref_id,
            None => true,
        };
        if needs_refresh {
            let seq = self.repository.get(contig.name.as_bytes())?.ok()?;
            self.cached_contig = Some((ref_id, seq));
        }
        let seq = &self.cached_contig.as_ref().expect("just set").1;

        let raw: &[u8] = AsRef::<[u8]>::as_ref(seq.as_ref());
        let start = (pos_1based as usize).checked_sub(1)?;
        let end = start.checked_add(ref_span as usize)?;
        if end > raw.len() {
            return None;
        }
        Some(&raw[start..end])
    }
}

/// Why [`process_read`] dropped a read instead of producing a
/// `PreparedRead`. Each variant maps onto exactly one of the
/// pipeline's existing counters, so the parallel stage can roll these
/// up into `FilterCounts` / `BaqSkipCounts` with the same totals the
/// old serial path produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    /// G2 — `cigar_is_bad`. Counts into `FilterCounts.bad_cigar`.
    BadCigar,
    /// F1 — `read_exceeds_mismatch_fraction`. Counts into
    /// `FilterCounts.high_mismatch_fraction`.
    HighMismatchFraction,
    /// BAQ skipped the read. Counts into `BaqSkipCounts` by reason.
    Baq(BaqSkipReason),
}

/// Result of [`process_read`]: a fully prepared read, or a structured
/// drop reason.
#[derive(Debug)]
pub enum ReadOutcome {
    Prepared(PreparedRead),
    Dropped(DropReason),
}

/// Run the full per-read fold (G2 → F3 → F1 → BAQ) on one read.
///
/// `baq` carries the per-worker BAQ engine and its uppercased-window
/// fetcher when BAQ is enabled, or `None` under `--no-baq` — in which
/// case G2/F3/F1 still run and the surviving read is built by
/// [`prepare_passthrough`] (raw base qualities, no capping). `raw_ref`
/// (and the BAQ fetcher, when present) must be bound to the read's
/// contig; the caller guarantees this by processing one single-contig
/// packet at a time, exactly as the BAQ stage already does. See the
/// module-level note for why F3/F1 and BAQ use different reference
/// sources.
///
/// **Reference-fetch failure is not an error.** When the read has no
/// reference span, or its window runs past the contig end, F3/F1 are
/// silently skipped and the read proceeds unchanged — matching the old
/// reader, whose `fetch_ref_for_read` returned `None` in those same
/// cases. BAQ then applies its own past-chrom-end skip.
pub fn process_read(
    mut read: MappedRead,
    baq: Option<(&mut BaqEngine, &mut ManualEvictChromRefFetcher)>,
    raw_ref: &mut RawContigRefCache,
    cfg: &ReadProcessingConfig,
) -> ReadOutcome {
    // G2 — reject before left-alignment so the check sees the
    // aligner's CIGAR, not a canonicalised one.
    if cigar_is_bad(&read.cigar) {
        return ReadOutcome::Dropped(DropReason::BadCigar);
    }

    // F3 + F1 both need the read's raw reference slice. Fetch once and
    // share. A zero ref span or an out-of-bounds window means "skip
    // F3/F1" (the old `fetch_ref_for_read` returned `None` here).
    let ref_span = cigar_ref_span(&read.cigar);
    if let Some(ref_seq) = raw_ref.fetch_raw_slice(read.ref_id, read.pos, ref_span) {
        // F3 — left-align indels in place. Read sequence/qualities are
        // untouched; only the CIGAR's indel anchoring moves. The
        // reference span is invariant under left-alignment, so the
        // slice stays valid for F1 below.
        left_align_indels(&mut read.cigar, &read.seq, ref_seq);

        // F1 — drop reads whose mismatch fraction (against the now
        // left-aligned CIGAR) exceeds the threshold.
        if let Some(threshold) = cfg.max_read_mismatch_fraction
            && read_exceeds_mismatch_fraction(
                &read.cigar,
                &read.seq,
                &read.qual,
                ref_seq,
                cfg.mismatch_bq_floor,
                threshold,
            )
        {
            return ReadOutcome::Dropped(DropReason::HighMismatchFraction);
        }
    }

    match baq {
        // BAQ on — caps base qualities and builds the PreparedRead.
        // Consumes the read by value so its buffers move into the
        // PreparedRead.
        Some((engine, baq_fetcher)) => match engine.process(read, baq_fetcher) {
            BaqOutcome::Capped(prepared) => ReadOutcome::Prepared(prepared),
            BaqOutcome::Skipped(reason) => ReadOutcome::Dropped(DropReason::Baq(reason)),
        },
        // --no-baq — passthrough build with raw base qualities. The
        // `ref_id` fits `u32` for every read the merged reader emits
        // (mapped, in-range); the old passthrough path asserted the
        // same.
        None => {
            let chrom_id = u32::try_from(read.ref_id).expect("ref_id fits u32");
            ReadOutcome::Prepared(prepare_passthrough(read, chrom_id))
        }
    }
}
