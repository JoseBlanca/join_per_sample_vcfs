//! Iterator adapter that sits between `AlignmentMergedReader` and
//! `pileup::walker::run`. Pulls a chunk of coordinate-sorted
//! `MappedRead`s, BAQ-caps each in parallel via rayon (one engine per
//! worker thread), and yields the resulting `PreparedRead`s in the
//! original order. The walker's
//! non-decreasing-coordinate invariant is preserved by rayon's
//! order-preserving `par_iter` plus our serial chunk-by-chunk drain.

use std::collections::VecDeque;
use std::path::PathBuf;

use noodles_fasta as fasta;
use rayon::prelude::*;

use crate::bam::alignment_input::MappedRead;
use crate::bam::errors::AlignmentInputError;
use crate::fasta::ContigList;
use crate::fasta::ManualEvictChromRefFetcher;
use crate::pileup::walker::PreparedRead;

use super::baq_engine::{BaqEngine, BaqSkipReason};
use super::read_processor::{
    DropReason, RawContigRefCache, ReadOutcome, ReadProcessingConfig, process_read,
};
use crate::baq::BaqConfig;

/// Default chunk size — reads per rayon batch. Picked at 1024 from the
/// `baq_stream_chunk_size` criterion bench
/// ([benches/baq_perf.rs](../../benches/baq_perf.rs)): on the
/// reference workstation (Intel i7-1260P, 16 logical CPUs) `/1024`
/// runs at ~88 ms vs `/4096` at ~85 ms — under 5 % wall-time gap for
/// 4× the memory footprint per chunk, so 1024 is the better default
/// for typical per-sample workloads. See
/// [doc/devel/reports/reviews/perf_baq_2026-05-12.md](../../doc/devel/reports/reviews/perf_baq_2026-05-12.md)
/// for the saved baseline.
pub const DEFAULT_BAQ_CHUNK_SIZE: usize = 1024;

/// Per-reason BAQ skip counters. `total` is what the pipeline rolls
/// up into `FilterCounts.baq_rejected`; the per-reason fields exist
/// for telemetry / debugging so a deployment can tell which skip
/// dominates without re-running with extra logging.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct BaqSkipCounts {
    pub total: u64,
    pub unmapped: u64,
    pub empty_query: u64,
    pub qual_absent: u64,
    pub no_match_in_cigar: u64,
    pub contains_ref_skip: u64,
    pub hmm_overflow: u64,
    pub ref_window_past_chrom_end: u64,
    pub pos_out_of_range: u64,
    pub read_too_long: u64,
    pub chrom_id_out_of_range: u64,
}

impl BaqSkipCounts {
    pub(super) fn bump(&mut self, reason: BaqSkipReason) {
        self.total += 1;
        match reason {
            BaqSkipReason::Unmapped => self.unmapped += 1,
            BaqSkipReason::EmptyQuery => self.empty_query += 1,
            BaqSkipReason::QualAbsent => self.qual_absent += 1,
            BaqSkipReason::NoMatchInCigar => self.no_match_in_cigar += 1,
            BaqSkipReason::ContainsRefSkip => self.contains_ref_skip += 1,
            BaqSkipReason::HmmOverflow => self.hmm_overflow += 1,
            BaqSkipReason::RefWindowPastChromEnd => self.ref_window_past_chrom_end += 1,
            BaqSkipReason::PosOutOfRange => self.pos_out_of_range += 1,
            BaqSkipReason::ReadTooLong => self.read_too_long += 1,
            BaqSkipReason::ChromIdOutOfRange => self.chrom_id_out_of_range += 1,
        }
    }

    /// Add `other`'s tallies into `self`, field by field. Totals the
    /// per-region BAQ skip counts when the pileup runs region by region.
    pub fn merge(&mut self, other: &BaqSkipCounts) {
        // Exhaustive destructure (no `..`): a new BaqSkipCounts field is a
        // compile error here until it is explicitly folded in (review M3).
        let BaqSkipCounts {
            total,
            unmapped,
            empty_query,
            qual_absent,
            no_match_in_cigar,
            contains_ref_skip,
            hmm_overflow,
            ref_window_past_chrom_end,
            pos_out_of_range,
            read_too_long,
            chrom_id_out_of_range,
        } = *other;
        self.total += total;
        self.unmapped += unmapped;
        self.empty_query += empty_query;
        self.qual_absent += qual_absent;
        self.no_match_in_cigar += no_match_in_cigar;
        self.contains_ref_skip += contains_ref_skip;
        self.hmm_overflow += hmm_overflow;
        self.ref_window_past_chrom_end += ref_window_past_chrom_end;
        self.pos_out_of_range += pos_out_of_range;
        self.read_too_long += read_too_long;
        self.chrom_id_out_of_range += chrom_id_out_of_range;
    }
}

/// Iterator adapter: `Iterator<Item = Result<MappedRead, AlignmentInputError>>`
/// → `Iterator<Item = Result<PreparedRead, AlignmentInputError>>` with a
/// rayon-parallel BAQ pass in between.
///
/// Errors from the upstream are propagated **after** any
/// already-batched successful reads, so the walker sees the survivors
/// first and the error becomes the final iterator item. The iterator
/// is fused: once exhausted (either by upstream returning `None` or
/// after returning an upstream `Err`), it returns `None` forever.
///
/// **Memory shape (after the `unified_chrom_ref_fetcher` migration).**
/// Each rayon worker constructs its own
/// [`ManualEvictChromRefFetcher`] in the `map_init` closure: per-worker
/// fetchers, no shared state, no `Arc`/`Mutex` on the fetch hot path.
/// After processing each read, the worker calls
/// `fetcher.evict_before(read.pos)` so the buffer stays at the
/// current-position window. Chunks are constrained not to span
/// chrom transitions (see `refill_batch`); each chunk's fetcher
/// binds to one contig.
pub struct BaqStream<R> {
    reads: R,
    cfg: BaqConfig,
    /// F1 mismatch-fraction config, applied per read inside the worker.
    proc_cfg: ReadProcessingConfig,
    /// Whether to BAQ-cap each read. `false` is the `--no-baq` path:
    /// G2/F3/F1 still run, but the surviving read is built by
    /// `prepare_passthrough` (raw base qualities). When `false` no
    /// per-worker BAQ engine or uppercased-window fetcher is built.
    apply_baq: bool,
    /// Path to the reference FASTA. Each rayon worker opens its own
    /// `File` from this path on first use; cheap (one open per worker).
    /// Drives the per-worker BAQ fetcher (uppercased windows).
    fasta_path: PathBuf,
    /// Shared FASTA repository for the per-worker F3/F1 raw-byte cache
    /// (case-preserving — see [`RawContigRefCache`]). `Arc`-backed, so
    /// cloning one per worker is cheap and contig bytes stay shared.
    repository: fasta::Repository,
    /// Contigs table for `chrom_id → name` resolution. Cloned cheaply
    /// into the rayon closure as a borrowed reference.
    contigs: ContigList,
    chunk_size: usize,
    /// Reusable input-chunk buffer drained from upstream each refill.
    /// Its capacity sticks at `chunk_size` after the first refill.
    chunk_buf: Vec<MappedRead>,
    /// A read picked from upstream that doesn't belong in the
    /// current chunk (different chrom). Stashed so the next batch
    /// starts with it. Preserves the chunk-per-chrom invariant.
    pending_next_read: Option<MappedRead>,
    /// Reusable parallel-output buffer; `collect_into_vec` clears and
    /// refills it in place.
    outcomes_buf: Vec<ReadOutcome>,
    /// Surviving `PreparedRead`s for the current chunk, drained
    /// front-to-back via `pop_front`. Lives on the stream so the
    /// backing buffer survives across refills.
    current_batch: VecDeque<PreparedRead>,
    pending_error: Option<AlignmentInputError>,
    upstream_done: bool,
    skip_counts: BaqSkipCounts,
    /// G2 rejections (`cigar_is_bad`) tallied in the worker. Rolled
    /// into `FilterCounts.bad_cigar` by the pipeline; counted here now
    /// that G2 runs in the read-processing stage rather than the reader.
    bad_cigar_count: u64,
    /// F1 rejections (`read_exceeds_mismatch_fraction`) tallied in the
    /// worker. Rolled into `FilterCounts.high_mismatch_fraction`.
    high_mismatch_count: u64,
}

impl<R> BaqStream<R>
where
    R: Iterator<Item = Result<MappedRead, AlignmentInputError>>,
{
    /// Construct a `BaqStream`. Panics if `chunk_size == 0` — a
    /// programmer error in every call site (silently coercing to 1
    /// would mask a misconfigured CLI flag as a perf bug).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reads: R,
        cfg: BaqConfig,
        proc_cfg: ReadProcessingConfig,
        apply_baq: bool,
        fasta_path: PathBuf,
        repository: fasta::Repository,
        contigs: ContigList,
        chunk_size: usize,
    ) -> Self {
        assert!(chunk_size > 0, "BaqStream chunk_size must be > 0");
        Self {
            reads,
            cfg,
            proc_cfg,
            apply_baq,
            fasta_path,
            repository,
            contigs,
            chunk_size,
            chunk_buf: Vec::new(),
            pending_next_read: None,
            outcomes_buf: Vec::new(),
            current_batch: VecDeque::new(),
            pending_error: None,
            upstream_done: false,
            skip_counts: BaqSkipCounts::default(),
            bad_cigar_count: 0,
            high_mismatch_count: 0,
        }
    }

    /// Skip counts accumulated so far. Useful for the per-sample run
    /// summary; the walker bumping into the upstream's `FilterCounts`
    /// reads `total` here.
    pub fn skip_counts(&self) -> &BaqSkipCounts {
        &self.skip_counts
    }

    /// G2 (`bad_cigar`) rejections tallied in the read-processing
    /// stage so far. The pipeline rolls this into
    /// `FilterCounts.bad_cigar`.
    pub fn bad_cigar_count(&self) -> u64 {
        self.bad_cigar_count
    }

    /// F1 (`high_mismatch_fraction`) rejections tallied in the
    /// read-processing stage so far. The pipeline rolls this into
    /// `FilterCounts.high_mismatch_fraction`.
    pub fn high_mismatch_count(&self) -> u64 {
        self.high_mismatch_count
    }

    fn refill_batch(&mut self) {
        if self.upstream_done && self.pending_next_read.is_none() {
            return;
        }
        self.chunk_buf.clear();

        // If a read from the previous batch was stashed because its
        // chrom_id didn't match, that read seeds this batch.
        if let Some(r) = self.pending_next_read.take() {
            self.chunk_buf.push(r);
        }

        // Fill until chunk_size or a chrom transition (whichever
        // comes first). The chunk-per-chrom invariant lets each
        // worker construct a single-contig fetcher for the chunk's
        // contig — the new `ManualEvictChromRefFetcher` is bound to
        // one contig at construction.
        while self.chunk_buf.len() < self.chunk_size {
            match self.reads.next() {
                Some(Ok(r)) => {
                    if let Some(first) = self.chunk_buf.first()
                        && r.ref_id != first.ref_id
                    {
                        // Different chrom — stash for the next batch.
                        self.pending_next_read = Some(r);
                        break;
                    }
                    self.chunk_buf.push(r);
                }
                Some(Err(e)) => {
                    self.pending_error = Some(e);
                    self.upstream_done = true;
                    break;
                }
                None => {
                    self.upstream_done = true;
                    break;
                }
            }
        }

        if self.chunk_buf.is_empty() {
            return;
        }

        // All reads in `chunk_buf` share `ref_id` by construction
        // (see the chrom-boundary check above).
        let chunk_ref_id = self.chunk_buf[0].ref_id;
        let contig_entry = match self.contigs.entries.get(chunk_ref_id) {
            Some(e) => e,
            None => {
                // Out-of-range ref_id — bail every read in this
                // chunk as ChromIdOutOfRange and continue.
                let n = self.chunk_buf.len() as u64;
                self.chunk_buf.clear();
                self.skip_counts.total += n;
                self.skip_counts.chrom_id_out_of_range += n;
                self.current_batch.clear();
                return;
            }
        };
        let contig_name: &str = &contig_entry.name;
        let fasta_path: &std::path::Path = &self.fasta_path;
        let cfg = self.cfg;
        let proc_cfg = &self.proc_cfg;
        let apply_baq = self.apply_baq;
        let repository = &self.repository;
        let contigs = &self.contigs;

        // par_drain consumes the chunk by-value (so the per-read fold
        // can move the read's cigar/seq/qname straight into the
        // resulting PreparedRead) while preserving both input order
        // and the chunk_buf's backing capacity. Each rayon worker
        // builds its own per-worker scratch via `map_init` — no shared
        // mutable state. `RawContigRefCache` serves F3/F1's
        // case-preserving bytes (always); the `(BaqEngine,
        // ManualEvictChromRefFetcher)` pair serves BAQ's uppercased
        // windows and is built only when BAQ is enabled. After
        // processing each read the worker calls
        // `fetcher.evict_before(read.pos)` so the per-worker BAQ buffer
        // stays at "current-position window plus forward growth";
        // see ManualEvictChromRefFetcher's module-level comment.
        self.chunk_buf
            .par_drain(..)
            .map_init(
                || {
                    let baq = apply_baq.then(|| {
                        let fetcher =
                            ManualEvictChromRefFetcher::for_contig(fasta_path, contig_name).expect(
                                "BAQ per-worker fetcher construction failed; FASTA path + \
                                 contig name were validated at BaqStream construction time",
                            );
                        (BaqEngine::new(cfg), fetcher)
                    });
                    let raw_ref = RawContigRefCache::new(repository.clone(), contigs.clone());
                    (baq, raw_ref)
                },
                |(baq, raw_ref), read| {
                    let pos = read.pos;
                    let baq_args = baq.as_mut().map(|(engine, fetcher)| (engine, fetcher));
                    let outcome = process_read(read, baq_args, raw_ref, proc_cfg);
                    if let Some((_, fetcher)) = baq.as_mut()
                        && let Ok(p) = u32::try_from(pos)
                    {
                        fetcher.evict_before(p);
                    }
                    outcome
                },
            )
            .collect_into_vec(&mut self.outcomes_buf);

        self.current_batch.clear();
        for outcome in self.outcomes_buf.drain(..) {
            match outcome {
                ReadOutcome::Prepared(p) => self.current_batch.push_back(p),
                ReadOutcome::Dropped(DropReason::Baq(reason)) => self.skip_counts.bump(reason),
                ReadOutcome::Dropped(DropReason::BadCigar) => self.bad_cigar_count += 1,
                ReadOutcome::Dropped(DropReason::HighMismatchFraction) => {
                    self.high_mismatch_count += 1
                }
            }
        }
    }
}

impl<R> Iterator for BaqStream<R>
where
    R: Iterator<Item = Result<MappedRead, AlignmentInputError>>,
{
    type Item = Result<PreparedRead, AlignmentInputError>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(p) = self.current_batch.pop_front() {
                return Some(Ok(p));
            }
            if let Some(e) = self.pending_error.take() {
                return Some(Err(e));
            }
            if self.upstream_done && self.pending_next_read.is_none() {
                return None;
            }
            self.refill_batch();
        }
    }
}
