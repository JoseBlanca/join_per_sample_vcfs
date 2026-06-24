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
    /// The merger (walking the sorted catalog) asked for a locus *past* one this
    /// sample still holds — so the sample contains a locus the catalog does not
    /// list, i.e. it was genotyped against a different catalog.
    #[error(
        "sample contains locus {held:?} not present in the catalog \
         (asked for {asked:?}; built against a different catalog?)"
    )]
    LocusNotInCatalog { held: LocusId, asked: LocusId },
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
    /// and advance; `q < held` → Absent (the sample lacks `q`); `q > held` → the sample
    /// holds a locus the (sorted) catalog walked past, i.e. it was built against a
    /// different catalog → a hard [`SsrCohortReadError::LocusNotInCatalog`]. A rewind
    /// (`q <= last_query`) is caught by the `debug_assert` before the match — the merger
    /// validates catalog monotonicity, so a rewind cannot occur from on-disk data and
    /// this guard is belt-and-suspenders only.
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

        let front = self.held.as_ref().map(|(locus, _)| *locus);
        match front.map(|f| q.cmp(&f)) {
            None => Ok(None),                 // exhausted → Absent forever after
            Some(Ordering::Less) => Ok(None), // forward, sample lacks q → Absent
            Some(Ordering::Equal) => {
                // PANIC-FREE: `Equal` is only produced when `front` (a clone of
                // `self.held`) is `Some`, so `held` is `Some` here.
                let (_, evidence) = self.held.take().expect("held is Some on Equal");
                self.advance()?; // preload the next stored locus
                Ok(Some(evidence))
            }
            Some(Ordering::Greater) => Err(SsrCohortReadError::LocusNotInCatalog {
                // PANIC-FREE: any `Some(Ordering::_)` means `front` was `Some`.
                held: front.expect("front is Some on Greater"),
                asked: q,
            }),
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
        // The SSR decoder rejects `start == 0` / `span == 0` on read
        // (`registry_ssr` coordinate guards), so `start >= 1` and `end > start`
        // hold here and neither subtraction can underflow.
        debug_assert!(
            record.start >= 1 && record.end > record.start,
            "decoder guarantees well-formed 1-based coordinates"
        );
        let locus = LocusId {
            chrom_id,
            start: record.start - 1, // 1-based start → 0-based start
            end: record.end - 1,     // 1-based exclusive end → 0-based exclusive end
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
            n_widened: 0,
            n_window_truncated: 0,
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

    /// Pins the Stage-1 writer's `+1` coordinate shift against the cursor's `-1`:
    /// a 0-based tally written through the *real* `to_container_record` and read back
    /// through the cursor must recover its original 0-based coordinates. Catches a
    /// drift on either side that the hand-written fixtures cannot.
    #[test]
    fn stage1_to_cohort_coordinate_round_trip() {
        use crate::ssr::pileup::driver::to_container_record;
        use crate::ssr::pileup::locus_tally::SsrLocusObs;
        use std::collections::HashMap;

        let tally = SsrLocusObs {
            chrom: "chr1".into(),
            start: 16, // 0-based half-open [16, 22)
            end: 22,
            depth: 9,
            n_filtered: 1,
            mapped_reads: 10,
            n_low_quality: 0,
            n_border_off_end: 0,
            n_widened: 0,
            n_window_truncated: 0,
            observed: obs(&[(b"CACACA", 9)]),
        };
        let name_to_id = HashMap::from([("chr1", 0u32)]);
        let container = to_container_record(tally, &name_to_id).unwrap();

        let reader = PspReader::new(Cursor::new(ssr_psp(&[container], 64))).unwrap();
        let mut c = SampleEvidenceCursor::new(reader, vec![0, 1].into()).unwrap();

        let ev = c
            .evidence_at(locus(0, 16, 22)) // original 0-based coordinates
            .unwrap()
            .expect("the round-tripped locus is present at its original coordinates");
        assert_eq!(ev.seq_counts, obs(&[(b"CACACA", 9)]));
        assert_eq!(ev.qc.depth, 9);
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
    fn skipping_a_held_locus_is_a_hard_error() {
        let mut c = cursor(&[rec(0, 51, 57, obs(&[(b"CA", 1)]))], 64, &[0, 1]);
        // held = 0-based [50,56]; asking ahead (q > held) means the sample holds a
        // locus the catalog lacks → typed hard error, not a panic.
        let err = c.evidence_at(locus(0, 60, 66)).unwrap_err();
        assert!(matches!(err, SsrCohortReadError::LocusNotInCatalog { .. }));
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
