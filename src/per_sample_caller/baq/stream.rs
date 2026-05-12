//! Iterator adapter that sits between `CramMergedReader` and
//! `pileup::walker::run`. Pulls a chunk of coordinate-sorted
//! `MappedRead`s, BAQ-caps each in parallel via rayon (one engine per
//! worker thread), and yields the resulting `PreparedRead`s in the
//! original order. The walker's
//! non-decreasing-coordinate invariant is preserved by rayon's
//! order-preserving `par_iter` plus our serial chunk-by-chunk drain.

use rayon::prelude::*;

use crate::per_sample_caller::cram_input::MappedRead;
use crate::per_sample_caller::errors::CramInputError;
use crate::per_sample_caller::pileup::{PreparedRead, RefSeqFetcher};

use super::BaqConfig;
use super::engine::{BaqEngine, BaqOutcome, BaqSkipReason};

/// Default chunk size — reads per rayon batch. 1024 is the plan's
/// starting point
/// ([baq.md commit 3](../../../ia/feature_implementation_plans/baq.md));
/// the benchmark in commit 4 will tune this.
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
}

impl BaqSkipCounts {
    fn bump(&mut self, reason: BaqSkipReason) {
        self.total += 1;
        match reason {
            BaqSkipReason::Unmapped => self.unmapped += 1,
            BaqSkipReason::EmptyQuery => self.empty_query += 1,
            BaqSkipReason::QualAbsent => self.qual_absent += 1,
            BaqSkipReason::NoMatchInCigar => self.no_match_in_cigar += 1,
            BaqSkipReason::ContainsRefSkip => self.contains_ref_skip += 1,
            BaqSkipReason::HmmOverflow => self.hmm_overflow += 1,
            BaqSkipReason::RefWindowPastChromEnd => self.ref_window_past_chrom_end += 1,
        }
    }
}

/// Iterator adapter: `Iterator<Item = Result<MappedRead, CramInputError>>`
/// → `Iterator<Item = Result<PreparedRead, CramInputError>>` with a
/// rayon-parallel BAQ pass in between.
///
/// Errors from the upstream are propagated **after** any
/// already-batched successful reads, so the walker sees the survivors
/// first and the error becomes the final iterator item. The iterator
/// is fused: once exhausted (either by upstream returning `None` or
/// after returning an upstream `Err`), it returns `None` forever.
pub struct BaqStream<'a, R, F>
where
    F: RefSeqFetcher + Sync,
{
    reads: R,
    cfg: BaqConfig,
    ref_fetcher: &'a F,
    chunk_size: usize,
    /// Reusable input-chunk buffer drained from upstream each refill.
    /// Its capacity sticks at `chunk_size` after the first refill.
    chunk_buf: Vec<MappedRead>,
    /// Reusable parallel-output buffer; `collect_into_vec` clears and
    /// refills it in place.
    outcomes_buf: Vec<BaqOutcome>,
    /// Surviving `PreparedRead`s for the current chunk, drained
    /// front-to-back via `current_idx`. Lives on the stream so the
    /// backing buffer survives across refills.
    current_batch: Vec<PreparedRead>,
    pending_error: Option<CramInputError>,
    upstream_done: bool,
    skip_counts: BaqSkipCounts,
}

impl<'a, R, F> BaqStream<'a, R, F>
where
    R: Iterator<Item = Result<MappedRead, CramInputError>>,
    F: RefSeqFetcher + Sync,
{
    pub fn new(reads: R, cfg: BaqConfig, ref_fetcher: &'a F, chunk_size: usize) -> Self {
        Self {
            reads,
            cfg,
            ref_fetcher,
            chunk_size: chunk_size.max(1),
            chunk_buf: Vec::new(),
            outcomes_buf: Vec::new(),
            current_batch: Vec::new(),
            pending_error: None,
            upstream_done: false,
            skip_counts: BaqSkipCounts::default(),
        }
    }

    /// Skip counts accumulated so far. Useful for the per-sample run
    /// summary; the walker bumping into the upstream's `FilterCounts`
    /// reads `total` here.
    pub fn skip_counts(&self) -> &BaqSkipCounts {
        &self.skip_counts
    }

    fn refill_batch(&mut self) {
        if self.upstream_done {
            return;
        }
        self.chunk_buf.clear();
        while self.chunk_buf.len() < self.chunk_size {
            match self.reads.next() {
                Some(Ok(r)) => self.chunk_buf.push(r),
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
        let cfg = self.cfg;
        let fetcher = self.ref_fetcher;
        // collect_into_vec clears the destination first and reuses its
        // backing buffer — no per-chunk allocation of `outcomes`.
        self.chunk_buf
            .par_iter()
            .map_init(
                || BaqEngine::new(cfg),
                |engine, read| engine.process(read, fetcher),
            )
            .collect_into_vec(&mut self.outcomes_buf);
        self.current_batch.clear();
        for outcome in self.outcomes_buf.drain(..) {
            match outcome {
                BaqOutcome::Capped(p) => self.current_batch.push(p),
                BaqOutcome::Skipped(reason) => self.skip_counts.bump(reason),
            }
        }
        // Reverse so `next()` can `pop()` from the back and still yield
        // in the original (coordinate-sorted) order.
        self.current_batch.reverse();
    }
}

impl<R, F> Iterator for BaqStream<'_, R, F>
where
    R: Iterator<Item = Result<MappedRead, CramInputError>>,
    F: RefSeqFetcher + Sync,
{
    type Item = Result<PreparedRead, CramInputError>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(p) = self.current_batch.pop() {
                return Some(Ok(p));
            }
            if let Some(e) = self.pending_error.take() {
                return Some(Err(e));
            }
            if self.upstream_done {
                return None;
            }
            self.refill_batch();
        }
    }
}
