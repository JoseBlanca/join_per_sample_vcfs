//! Span-addressable columnar reader over one sample's `.psp`.
//!
//! Wraps a [`BlockColumnReader`] and turns it into a *genomic-span*
//! interface: the caller asks for the
//! next servable span ([`peek_next_span`](ColumnSpanReader::peek_next_span),
//! free — from the block index) and then for the columns over a span
//! ([`read_span`](ColumnSpanReader::read_span), appended straight into a
//! [`SampleColumns`] with no row-shape `PileupRecord`). PSP block
//! boundaries are an implementation detail hidden here: a `read_span`
//! that ends inside a block decodes that block once and resumes from the
//! same cursor on the next call (the "tail buffer" is just the loaded
//! block plus the record cursor — no extra copy).
//!
//! This is the per-sample seam the streaming cohort producer drives:
//! peek every sample, pick the consensus span `W = min next end`, ask
//! every sample for `read_span(W)`, and the samples stay aligned with
//! no producer-side block bookkeeping. Alignment is a *convention* the
//! producer follows by asking for the peeked span — `read_span` accepts
//! any `end`.

use std::io::{Read, Seek};

use crate::psp::{BlockColumnReader, PspReadError, PspReader, ScalarDecodeError};
use crate::var_calling::columns::SampleColumns;

/// Per-sample columnar reader bound to one chromosome region. Serves
/// records in `[region_start, region_end]` (1-based inclusive) in
/// ascending position order across successive `read_span` calls.
pub struct ColumnSpanReader<R: Read + Seek> {
    blocks: BlockColumnReader<R>,
    chrom_id: u32,
    /// Inclusive upper bound of the region; records past it are never
    /// served (the next region / chromosome owns them).
    region_end: u32,
    /// Inclusive lower bound; records below it are skipped (only bites
    /// in the region's first block, whose start may precede it).
    floor: u32,
    /// Cursor within the currently-decoded block: index of the next
    /// unserved record, its cumulative allele offset, and its absolute
    /// 1-based position. `next_pos` is `None` when no block is decoded
    /// (a `read_span` will load the block at the cursor).
    next_record: usize,
    next_allele: usize,
    next_pos: Option<u32>,
    /// Reusable absolute-position buffer for the per-span append.
    abs_scratch: Vec<u32>,
}

impl<R: Read + Seek> ColumnSpanReader<R> {
    /// Open a reader over `[region_start, region_end]` of `chrom_id`,
    /// taking ownership of `reader`.
    pub fn new(reader: PspReader<R>, chrom_id: u32, region_start: u32, region_end: u32) -> Self {
        let mut this = Self::detached(reader);
        this.reset(chrom_id, region_start, region_end);
        this
    }

    /// Wrap `reader` without positioning it anywhere — `peek_next_span`
    /// is `None` and `read_span` serves nothing until [`reset`](Self::reset)
    /// names a region. The cohort producer builds one per sample up
    /// front and `reset`s it at every covered interval (reusing the
    /// decode buffers across the whole run).
    pub fn detached(reader: PspReader<R>) -> Self {
        Self {
            blocks: reader.into_column_blocks(),
            chrom_id: 0,
            region_end: 0,
            floor: 0,
            next_record: 0,
            next_allele: 0,
            next_pos: None,
            abs_scratch: Vec::new(),
        }
    }

    /// Recover the owned [`PspReader`].
    pub fn into_reader(self) -> PspReader<R> {
        self.blocks.into_reader()
    }

    /// The sample's on-disk block index — the producer unions these
    /// across samples to find the covered intervals (data-bearing
    /// segments) per chromosome.
    pub fn block_index(&self) -> &[crate::psp::BlockIndexEntry] {
        self.blocks.index()
    }

    /// Re-point the reader at a new region, reusing the decode buffers.
    /// The cohort producer calls this at every covered-interval /
    /// chromosome boundary so the zstd context + column slabs persist
    /// across the whole run.
    pub fn reset(&mut self, chrom_id: u32, region_start: u32, region_end: u32) {
        self.blocks.seek_to(chrom_id, region_start, region_end);
        self.chrom_id = chrom_id;
        self.region_end = region_end;
        self.floor = region_start;
        self.next_record = 0;
        self.next_allele = 0;
        self.next_pos = None;
    }

    /// Inclusive-end of the next span this reader can serve **without a
    /// new decode** (the rest of the loaded block, or — if none is
    /// loaded — the next block on this chromosome), clamped to
    /// `region_end`. `None` once the region is exhausted.
    ///
    /// The producer takes `min` of this across samples to pick the
    /// consensus span, so each subsequent `read_span` costs each sample
    /// at most one new block decode.
    pub fn peek_next_span(&self) -> Option<u32> {
        if let Some(p) = self.next_pos {
            if p > self.region_end {
                return None;
            }
            let last = self.blocks.index()[self.blocks.cur_block_idx()].last_pos;
            return Some(last.min(self.region_end));
        }
        let entry = self.blocks.peek_block()?;
        if entry.chrom_id != self.chrom_id || entry.first_pos > self.region_end {
            return None;
        }
        Some(entry.last_pos.min(self.region_end))
    }

    /// Append every not-yet-served record with position `< end` (and
    /// `>= region_start`, `<= region_end`) to `out`, in ascending
    /// order, decoding blocks as needed and leaving the cursor on the
    /// first record `>= end` for the next call.
    pub fn read_span(&mut self, end: u32, out: &mut SampleColumns) -> Result<(), PspReadError> {
        // Lift the scratch out of `self` so the per-block borrow of
        // `self.blocks` (which ties to all of `&self`) can coexist with
        // mutating the position buffer.
        let mut abs = std::mem::take(&mut self.abs_scratch);
        let result = self.read_span_inner(end, out, &mut abs);
        abs.clear();
        self.abs_scratch = abs;
        result
    }

    fn read_span_inner(
        &mut self,
        end: u32,
        out: &mut SampleColumns,
        abs: &mut Vec<u32>,
    ) -> Result<(), PspReadError> {
        // Serve `pos < end` and `pos <= region_end`.
        let upper = end.min(self.region_end.saturating_add(1));
        loop {
            // Decode the block at the cursor if none is live.
            if self.next_pos.is_none() {
                match self.blocks.peek_block() {
                    // `first_pos < upper` (with `upper <= region_end + 1`)
                    // means the block may carry a servable record.
                    Some(entry) if entry.chrom_id == self.chrom_id && entry.first_pos < upper => {
                        let loaded = self.blocks.load_current()?;
                        debug_assert!(loaded, "peek_block was Some, so the index is in range");
                        let first_pos = self.blocks.columns().expect("just loaded").first_pos;
                        self.next_record = 0;
                        self.next_allele = 0;
                        self.next_pos = Some(first_pos);
                    }
                    _ => return Ok(()),
                }
            }

            // Walk + append records of the live block, scoped so the
            // `cols` borrow ends before the cursor / block mutations.
            let (r_end, a_end, p_end, exhausted) = {
                let cols = self.blocks.columns().expect("block is loaded");
                let n_records = cols.n_records as usize;
                let floor = self.floor;
                let mut r = self.next_record;
                let mut a = self.next_allele;
                let mut p = self.next_pos.expect("block is loaded");

                // Skip records below the region floor (region's first
                // block only — later blocks start at or above it).
                while r < n_records && p < floor {
                    a += cols.n_alleles[r] as usize;
                    r += 1;
                    if r < n_records {
                        p = advance_pos(p, cols.delta_pos[r], r)?;
                    }
                }

                let serve_record_lo = r;
                let serve_allele_lo = a;
                abs.clear();
                while r < n_records && p < upper {
                    abs.push(p);
                    a += cols.n_alleles[r] as usize;
                    r += 1;
                    if r < n_records {
                        p = advance_pos(p, cols.delta_pos[r], r)?;
                    }
                }

                out.append_block_window(&cols, serve_record_lo, r, serve_allele_lo, abs);
                (r, a, p, r >= n_records)
            };

            self.next_record = r_end;
            self.next_allele = a_end;
            if exhausted {
                // Whole block served; move to the next and keep going —
                // the requested `end` may reach into it.
                self.blocks.advance();
                self.next_pos = None;
            } else {
                // Stopped on a record `>= end`; that's the resume point.
                self.next_pos = Some(p_end);
                return Ok(());
            }
        }
    }
}

impl<R: Read + Seek> crate::var_calling::loader::SpanColumnSource for ColumnSpanReader<R> {
    fn peek_next_span(&self) -> Option<u32> {
        ColumnSpanReader::peek_next_span(self)
    }

    fn read_span(&mut self, end: u32, out: &mut SampleColumns) -> Result<(), PspReadError> {
        ColumnSpanReader::read_span(self, end, out)
    }
}

/// Advance a delta-encoded position, mirroring the row reader's
/// overflow contract: a `delta` that does not fit `u32` is a corrupt /
/// truncated column, surfaced as the same typed error
/// `materialise_next_record` raises rather than a silent wrap.
fn advance_pos(prev: u32, delta: u64, entry: usize) -> Result<u32, PspReadError> {
    let delta_u32 = u32::try_from(delta).map_err(|_| PspReadError::ColumnElementDecode {
        column: "delta-pos".to_string(),
        entry,
        source: ScalarDecodeError::VarintOverflow,
    })?;
    Ok(prev.saturating_add(delta_u32))
}

#[cfg(test)]
mod tests {
    //! Equivalence tests: the columnar `read_span` path must build the
    //! exact same `SampleColumns` the row path (`region_records` →
    //! `push_record`) builds, for any window and any split into
    //! sub-spans — including splits that land inside a PSP block and
    //! across blocks that are misaligned between samples.

    use std::io::Cursor;

    use super::*;
    use crate::pileup_record::PileupRecord;
    use crate::psp::test_fixtures::writer_header;
    use crate::psp::writer::PspWriter;
    use crate::var_calling::test_helpers::{allele, record};

    /// A small fixture cohort of multi-allele records (varied allele
    /// counts, sequences, scalars, and chain-id lists) at ascending
    /// positions, dense enough that a small block target produces
    /// several blocks.
    fn fixture_records() -> Vec<PileupRecord> {
        (0..200u32)
            .map(|i| {
                let pos = 10 + i * 7;
                match i % 3 {
                    0 => record(pos, vec![allele(b"A", 12, -1.5, &[])]),
                    1 => record(
                        pos,
                        vec![
                            allele(b"A", 9, -1.0, &[]),
                            allele(b"T", 4, -2.0, &[(i as u64) + 1]),
                        ],
                    ),
                    _ => record(
                        pos,
                        vec![
                            allele(b"AC", 7, -0.5, &[]),
                            allele(b"AG", 3, -3.0, &[(i as u64) + 1, (i as u64) + 2]),
                            allele(b"A", 5, -1.2, &[]),
                        ],
                    ),
                }
            })
            .collect()
    }

    /// Serialize `records` to an in-memory `.psp` with the given block
    /// target (smaller → more, smaller blocks).
    fn psp_bytes(records: &[PileupRecord], block_target: usize) -> Vec<u8> {
        let mut writer = PspWriter::new_with_block_target(
            Cursor::new(Vec::new()),
            writer_header(1),
            block_target,
        )
        .expect("writer opens");
        for r in records {
            writer.write_record(r).expect("write_record");
        }
        writer.finish().expect("finish").into_inner()
    }

    /// Row path: `region_records` → `push_record` into `SampleColumns`.
    fn columns_via_records(bytes: &[u8], chrom: u32, start: u32, end: u32) -> SampleColumns {
        let mut reader = PspReader::new(Cursor::new(bytes.to_vec())).expect("reader opens");
        let mut sc = SampleColumns::empty();
        for r in reader.region_records(chrom, start, end) {
            sc.push_record(r.expect("record decodes"));
        }
        sc
    }

    /// Columnar path: `ColumnSpanReader::read_span` over `split_points`
    /// (then a final unbounded span to drain the rest).
    fn columns_via_span(
        bytes: &[u8],
        chrom: u32,
        start: u32,
        end: u32,
        split_points: &[u32],
    ) -> SampleColumns {
        let reader = PspReader::new(Cursor::new(bytes.to_vec())).expect("reader opens");
        let mut sc = SampleColumns::empty();
        let mut span = ColumnSpanReader::new(reader, chrom, start, end);
        for &sp in split_points {
            span.read_span(sp, &mut sc).expect("read_span");
        }
        span.read_span(u32::MAX, &mut sc).expect("read_span drain");
        sc
    }

    #[test]
    fn whole_region_one_span_matches_record_path() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        // Sanity: the fixture really is multi-block.
        let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
        assert!(
            reader.block_index().len() >= 3,
            "fixture should be multi-block"
        );

        let via_records = columns_via_records(&bytes, 0, 1, 1_000_000);
        let via_span = columns_via_span(&bytes, 0, 1, 1_000_000, &[]);
        assert_eq!(via_span, via_records);
    }

    #[test]
    fn split_spans_including_midblock_match_record_path() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        let via_records = columns_via_records(&bytes, 0, 1, 1_000_000);
        // Split points chosen to land between records and almost
        // certainly inside blocks (block boundaries are content-driven).
        let via_span = columns_via_span(&bytes, 0, 1, 1_000_000, &[100, 250, 251, 600, 999]);
        assert_eq!(via_span, via_records);
    }

    #[test]
    fn clamped_subwindow_matches_record_path() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        // A window whose ends fall mid-block on both sides.
        let (start, end) = (137, 941);
        let via_records = columns_via_records(&bytes, 0, start, end);
        let via_span = columns_via_span(&bytes, 0, start, end, &[400, 700]);
        assert_eq!(via_span, via_records);
        // The clamp really dropped records outside the window.
        let full = columns_via_records(&bytes, 0, 1, 1_000_000);
        assert!(via_records.n_records() < full.n_records());
    }

    #[test]
    fn empty_region_yields_no_records() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        // A gap with no records (positions are 10, 17, 24, ...; pick a
        // window between two of them).
        let via_span = columns_via_span(&bytes, 0, 11, 16, &[]);
        assert_eq!(via_span.n_records(), 0);
        assert_eq!(via_span, columns_via_records(&bytes, 0, 11, 16));
    }

    #[test]
    fn peek_reports_servable_ends_and_terminates() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut span = ColumnSpanReader::new(reader, 0, 1, 1_000_000);
        let mut sc = SampleColumns::empty();
        let mut peeks = Vec::new();
        // Drive by consensus-of-one: always ask for exactly the peeked
        // end, draining the file block by block.
        while let Some(end) = span.peek_next_span() {
            peeks.push(end);
            span.read_span(end.saturating_add(1), &mut sc).unwrap();
            if peeks.len() > 1000 {
                panic!("peek did not terminate");
            }
        }
        // Peeked ends are strictly ascending (one block's last_pos at a
        // time) and the block-by-block drain reconstructs the file.
        assert!(peeks.windows(2).all(|w| w[0] < w[1]), "peeks ascending");
        assert_eq!(sc, columns_via_records_owned(recs));
    }

    fn columns_via_records_owned(recs: Vec<PileupRecord>) -> SampleColumns {
        let mut sc = SampleColumns::empty();
        for r in recs {
            sc.push_record(r);
        }
        sc
    }

    /// The misalignment fixture: the *same* records written with two
    /// different block targets produce different block boundaries, yet
    /// driving each through `ColumnSpanReader` at a shared consensus
    /// span (`W = min(peek)`) reconstructs each sample exactly. This is
    /// the cross-sample scenario Stage 2's producer relies on.
    #[test]
    fn misaligned_blocks_served_at_consensus_span_match() {
        let recs = fixture_records();
        let bytes_a = psp_bytes(&recs, 256); // smaller blocks
        let bytes_b = psp_bytes(&recs, 1024); // larger blocks

        let reader_a = PspReader::new(Cursor::new(bytes_a.clone())).unwrap();
        let reader_b = PspReader::new(Cursor::new(bytes_b.clone())).unwrap();
        // Confirm the boundaries really differ.
        assert_ne!(
            reader_a.block_index().len(),
            reader_b.block_index().len(),
            "fixture should be misaligned across the two samples"
        );

        let mut span_a = ColumnSpanReader::new(reader_a, 0, 1, 1_000_000);
        let mut span_b = ColumnSpanReader::new(reader_b, 0, 1, 1_000_000);
        let mut sc_a = SampleColumns::empty();
        let mut sc_b = SampleColumns::empty();

        loop {
            let pa = span_a.peek_next_span();
            let pb = span_b.peek_next_span();
            let w = match (pa, pb) {
                (Some(a), Some(b)) => a.min(b),
                (None, None) => break,
                // One sample exhausted before the other: drain the rest
                // of whichever remains.
                (Some(a), None) => a,
                (None, Some(b)) => b,
            };
            span_a.read_span(w.saturating_add(1), &mut sc_a).unwrap();
            span_b.read_span(w.saturating_add(1), &mut sc_b).unwrap();
        }

        let expected = columns_via_records_owned(recs);
        assert_eq!(sc_a, expected, "small-block sample reconstructs");
        assert_eq!(sc_b, expected, "large-block sample reconstructs");
        assert_eq!(sc_a, sc_b, "both samples agree");
    }
}
