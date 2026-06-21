//! The per-sample cohort reader — [`SampleEvidenceCursor`] (arch
//! `doc/devel/architecture/ssr_call_reading.md` §3).
//!
//! One cursor per input `.ssr.psp`. It owns that sample's reader and yields the
//! sample's evidence at a queried locus, with the loci asked in **strictly ascending
//! catalog order**. It knows only its own sample; the cross-sample gather is the
//! merger's job (§4).
//!
//! The cursor "always holds the next locus" ([`held`](SampleEvidenceCursor::held)),
//! decoded one block at a time via an owning [`OwnedRecordsIter`], so the resident set
//! is one decoded block (not the whole file). Its monotonic-order contract — and why a
//! `last_query` guard makes the Absent case unambiguous — is documented on
//! [`evidence_at`](SampleEvidenceCursor::evidence_at).
//!
//! **Phase 1: synchronous inline decode.** The shared decode pool + prefetched futures
//! (arch §5) are a later, profiling-gated phase; here a block simply decodes inline
//! when the owning iterator crosses into it.

use std::cmp::Ordering;
use std::io::{Read, Seek};

use crate::psp::registry_ssr::{SsrKind, SsrLocusRecord};
use crate::psp::{OwnedRecordsIter, PspReadError, PspReader};
use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};

/// Errors raised while streaming a sample's `.ssr.psp` evidence.
#[derive(Debug, thiserror::Error)]
pub(crate) enum SsrCohortReadError {
    /// A decode error from the underlying `.ssr.psp` reader.
    #[error("decoding .ssr.psp")]
    Psp(#[from] PspReadError),
    /// A record names a per-file chromosome id with no entry in the cohort-global
    /// remap — the file and the catalog/cohort disagree on the chromosome table.
    #[error("record chromosome id {per_file_id} is outside the chromosome remap (len {remap_len})")]
    ChromIdOutOfRange { per_file_id: u32, remap_len: usize },
}

/// A per-sample coordinate cursor over one `.ssr.psp`: ask it for a locus (in
/// ascending catalog order) and it answers this sample's evidence there, or Absent.
pub(crate) struct SampleEvidenceCursor<R: Read + Seek> {
    /// The sample's records, streamed one block at a time. Owning, so the cursor can
    /// hold it across calls (§3 — a borrowing iterator would be self-referential).
    records: OwnedRecordsIter<R, SsrKind>,
    /// Per-file chromosome id → cohort-global id (indexed by the file's `chrom_id`).
    /// Built by the merger from the catalog's chromosome order; identity in tests.
    chrom_remap: Box<[u32]>,
    /// The cursor's next stored locus + its evidence, decoded ahead. `None` once the
    /// sample is exhausted. The cursor "always holds the next locus."
    held: Option<(LocusId, SampleEvidence)>,
    /// The last locus the merger asked for — the monotonic guard that separates a
    /// genuine Absent (`q` ahead of nothing this sample has) from a rewind (a bug).
    last_query: Option<LocusId>,
}

impl<R: Read + Seek> SampleEvidenceCursor<R> {
    /// Open a cursor over `reader`'s records, remapping per-file chromosome ids through
    /// `chrom_remap`. Preloads the first stored locus.
    ///
    /// The same-catalog precondition (every input declares the catalog the merger
    /// holds) is validated by the merger before it builds cursors, so it is not
    /// re-checked here.
    pub(crate) fn new(
        reader: PspReader<R>,
        chrom_remap: Box<[u32]>,
    ) -> Result<Self, SsrCohortReadError> {
        let mut cursor = Self {
            records: reader.into_records_of::<SsrKind>(),
            chrom_remap,
            held: None,
            last_query: None,
        };
        cursor.advance()?; // preload `held`
        Ok(cursor)
    }

    /// Ask for the sample's evidence at catalog locus `q`. **Precondition:** `q`
    /// strictly ascends across calls (catalog order).
    ///
    /// `Ok(Some(_))` = present; `Ok(None)` = Absent (this sample has no record at `q`,
    /// or it is exhausted); `Err` = a decode failure while refilling. (The arch doc
    /// writes the return as `Option<SampleEvidence>`; real decode can fail lazily on a
    /// block refill, so it is wrapped in a `Result` here.)
    ///
    /// Comparing `q` to the held front gives three outcomes: `q == held` → return it
    /// and advance; `q < held` → Absent (the sample lacks `q`); `q > held` → the merger
    /// walked past a locus this sample *has* (out-of-order, or catalog/sample drift) →
    /// panic. A rewind (`q <= last_query`) is caught by the guard *before* the match,
    /// so the `q < held` arm can only mean a forward Absent — never a backwards step.
    pub(crate) fn evidence_at(
        &mut self,
        q: LocusId,
    ) -> Result<Option<SampleEvidence>, SsrCohortReadError> {
        debug_assert!(
            self.last_query.is_none_or(|prev| q > prev),
            "evidence_at requires strictly ascending loci (out-of-order / rewind): \
             {q:?} after {:?}",
            self.last_query
        );
        self.last_query = Some(q);

        match self.held.as_ref().map(|(front, _)| q.cmp(front)) {
            None => Ok(None),                 // exhausted → Absent forever after
            Some(Ordering::Less) => Ok(None), // forward, sample lacks q → Absent
            Some(Ordering::Equal) => {
                let (_, evidence) = self.held.take().expect("held is Some on Equal");
                self.advance()?; // preload the next stored locus
                Ok(Some(evidence))
            }
            Some(Ordering::Greater) => panic!(
                "merger skipped a stored locus this sample has: asked {q:?}, holding {:?}",
                self.held.as_ref().map(|(front, _)| front)
            ),
        }
    }

    /// Refill `held` from the next stored record (crossing a block when the owning
    /// iterator needs to). `None` once no records remain.
    fn advance(&mut self) -> Result<(), SsrCohortReadError> {
        self.held = match self.records.next() {
            None => None,
            Some(Ok(record)) => Some(self.adapt(record)?),
            Some(Err(err)) => return Err(err.into()),
        };
        Ok(())
    }

    /// Convert a container record (per-file chromosome id, 1-based `[start, end)`) into
    /// the catalog frame (`LocusId`: cohort-global id, 0-based half-open) + the
    /// sample's evidence. Inverse of the Stage-1 writer's `+1` coordinate shift.
    fn adapt(
        &self,
        record: SsrLocusRecord,
    ) -> Result<(LocusId, SampleEvidence), SsrCohortReadError> {
        let chrom_id = *self.chrom_remap.get(record.chrom_id as usize).ok_or(
            SsrCohortReadError::ChromIdOutOfRange {
                per_file_id: record.chrom_id,
                remap_len: self.chrom_remap.len(),
            },
        )?;
        debug_assert!(record.start >= 1, "container coordinates are 1-based");
        let locus = LocusId {
            chrom_id,
            start: record.start - 1, // 1-based inclusive → 0-based inclusive
            end: record.end - 1,     // 1-based exclusive → 0-based exclusive
        };
        let evidence = SampleEvidence {
            seq_counts: record.observed,
            qc: SsrQc {
                depth: record.depth,
                n_filtered: record.n_filtered,
                mapped_reads: record.mapped_reads,
                n_low_quality: record.n_low_quality,
                n_border_off_end: record.n_border_off_end,
            },
        };
        Ok((locus, evidence))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::psp::test_fixtures::writer_header;
    use crate::psp::writer::PspWriter;
    use std::io::Cursor;

    fn obs(pairs: &[(&[u8], u32)]) -> Vec<(Box<[u8]>, u32)> {
        pairs
            .iter()
            .map(|(s, c)| (s.to_vec().into_boxed_slice(), *c))
            .collect()
    }

    /// A container record (1-based coords) with fixed QC scalars so the conversion is
    /// checkable.
    fn rec(chrom_id: u32, start: u32, end: u32, observed: Vec<(Box<[u8]>, u32)>) -> SsrLocusRecord {
        SsrLocusRecord {
            chrom_id,
            start,
            end,
            depth: 30,
            n_filtered: 2,
            mapped_reads: 33,
            n_low_quality: 1,
            n_border_off_end: 3,
            observed,
        }
    }

    /// Write a `.ssr.psp` in memory. `block_window_bp` small + spaced-out starts forces
    /// multiple blocks.
    fn ssr_psp(records: &[SsrLocusRecord], block_window_bp: u32) -> Vec<u8> {
        let mut writer = PspWriter::<_, SsrKind>::new_ssr_with_block_layout(
            Cursor::new(Vec::<u8>::new()),
            writer_header(2),
            16 * 1024 * 1024,
            block_window_bp,
        )
        .unwrap();
        for r in records {
            writer.write_locus(r).unwrap();
        }
        writer.finish().unwrap().into_inner()
    }

    fn cursor(
        records: &[SsrLocusRecord],
        block_window_bp: u32,
        remap: &[u32],
    ) -> SampleEvidenceCursor<Cursor<Vec<u8>>> {
        let bytes = ssr_psp(records, block_window_bp);
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        SampleEvidenceCursor::new(reader, remap.to_vec().into_boxed_slice()).unwrap()
    }

    fn locus(chrom_id: u32, start: u32, end: u32) -> LocusId {
        LocusId {
            chrom_id,
            start,
            end,
        }
    }

    #[test]
    fn streams_present_loci_with_coordinate_conversion() {
        // Container 1-based starts → expected 0-based LocusIds (start-1).
        let mut c = cursor(
            &[
                rec(0, 51, 57, obs(&[(b"CACACA", 12)])),
                rec(0, 61, 71, obs(&[(b"CACACA", 8), (b"CACAACA", 1)])),
                rec(1, 6, 12, obs(&[(b"AT", 4)])),
            ],
            64,
            &[0, 1],
        );

        let e0 = c.evidence_at(locus(0, 50, 56)).unwrap().expect("present");
        assert_eq!(e0.seq_counts, obs(&[(b"CACACA", 12)]));
        // QC scalars carried through from the container record.
        assert_eq!(e0.qc.depth, 30);
        assert_eq!(e0.qc.n_filtered, 2);
        assert_eq!(e0.qc.n_border_off_end, 3);

        let e1 = c.evidence_at(locus(0, 60, 70)).unwrap().expect("present");
        assert_eq!(e1.seq_counts.len(), 2);

        // chrom id 1, container [6,12) → 0-based [5,11).
        let e2 = c.evidence_at(locus(1, 5, 11)).unwrap().expect("present");
        assert_eq!(e2.seq_counts, obs(&[(b"AT", 4)]));
    }

    #[test]
    fn absent_for_a_locus_the_sample_lacks() {
        let mut c = cursor(
            &[
                rec(0, 51, 57, obs(&[(b"CA", 1)])),
                rec(0, 151, 163, obs(&[(b"CA", 2)])),
            ],
            64,
            &[0, 1],
        );
        assert!(c.evidence_at(locus(0, 50, 56)).unwrap().is_some());
        // q < held (next stored is 0-based 150) → Absent, does not consume the front.
        assert!(c.evidence_at(locus(0, 100, 106)).unwrap().is_none());
        // The later stored locus is still reachable.
        assert!(c.evidence_at(locus(0, 150, 162)).unwrap().is_some());
    }

    #[test]
    fn absent_before_the_first_stored_locus() {
        let mut c = cursor(&[rec(0, 51, 57, obs(&[(b"CA", 1)]))], 64, &[0, 1]);
        assert!(c.evidence_at(locus(0, 10, 16)).unwrap().is_none());
        assert!(c.evidence_at(locus(0, 50, 56)).unwrap().is_some());
    }

    #[test]
    fn absent_after_exhaustion_forever() {
        let mut c = cursor(&[rec(0, 51, 57, obs(&[(b"CA", 1)]))], 64, &[0, 1]);
        assert!(c.evidence_at(locus(0, 50, 56)).unwrap().is_some());
        assert!(c.evidence_at(locus(0, 60, 66)).unwrap().is_none());
        assert!(c.evidence_at(locus(0, 70, 76)).unwrap().is_none());
    }

    #[test]
    fn empty_observed_locus_is_present() {
        let mut c = cursor(&[rec(0, 151, 163, obs(&[]))], 64, &[0, 1]);
        let e = c
            .evidence_at(locus(0, 150, 162))
            .unwrap()
            .expect("present even with no observations");
        assert!(e.seq_counts.is_empty());
        assert_eq!(e.qc.depth, 30);
    }

    #[test]
    fn crosses_block_boundaries() {
        // Spaced-out starts + small window grid → distinct blocks.
        let records = [
            rec(0, 51, 57, obs(&[(b"CA", 1)])),
            rec(0, 101, 107, obs(&[(b"CA", 2)])),
            rec(0, 161, 167, obs(&[(b"CA", 3)])),
        ];
        let bytes = ssr_psp(&records, 16);
        // Confirm the fixture is genuinely multi-block.
        let probe = PspReader::new(Cursor::new(bytes.clone())).unwrap();
        assert!(
            probe.block_index().len() >= 2,
            "fixture must span >1 block, got {}",
            probe.block_index().len()
        );

        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut c = SampleEvidenceCursor::new(reader, vec![0, 1].into()).unwrap();
        assert_eq!(
            c.evidence_at(locus(0, 50, 56)).unwrap().unwrap().seq_counts,
            obs(&[(b"CA", 1)])
        );
        assert_eq!(
            c.evidence_at(locus(0, 100, 106))
                .unwrap()
                .unwrap()
                .seq_counts,
            obs(&[(b"CA", 2)])
        );
        assert_eq!(
            c.evidence_at(locus(0, 160, 166))
                .unwrap()
                .unwrap()
                .seq_counts,
            obs(&[(b"CA", 3)])
        );
    }

    #[test]
    fn chrom_remap_is_applied() {
        // per-file chrom id 1 → cohort-global id 7.
        let mut c = cursor(&[rec(1, 6, 12, obs(&[(b"AT", 4)]))], 64, &[3, 7]);
        assert!(c.evidence_at(locus(7, 5, 11)).unwrap().is_some());
    }

    #[test]
    #[should_panic(expected = "skipped a stored locus")]
    fn skipping_a_held_locus_panics() {
        let mut c = cursor(&[rec(0, 51, 57, obs(&[(b"CA", 1)]))], 64, &[0, 1]);
        // held = 0-based [50,56]; asking ahead of it (q > held) means the merger never
        // asked for a locus this sample has.
        let _ = c.evidence_at(locus(0, 60, 66));
    }

    #[test]
    #[should_panic(expected = "ascending")]
    fn rewind_panics() {
        let mut c = cursor(
            &[
                rec(0, 51, 57, obs(&[(b"CA", 1)])),
                rec(0, 61, 67, obs(&[(b"CA", 2)])),
            ],
            64,
            &[0, 1],
        );
        c.evidence_at(locus(0, 50, 56)).unwrap(); // consume the front; last_query=[50,56]
        let _ = c.evidence_at(locus(0, 40, 46)); // q < last_query → rewind
    }
}
