//! Pipelined read-processing stage.
//!
//! Replaces `BaqStream`'s bulk-synchronous batch barrier (pull a chunk →
//! `par_drain` → drain serially, with no overlap) with a three-role
//! pipeline over bounded `crossbeam` channels, so decode, the per-read
//! fold, and the walker all run concurrently:
//!
//! ```text
//!  SOURCE (1 thread)          WORKERS (M threads)        CONSUMER (walker thread)
//!  reader → packets    ─tx→   process_read per read  ─tx→  reorder by seq → walker
//!  (merge/dedup/order)        (G2·F3·F1·BAQ)               (serial, in order)
//! ```
//!
//! The source cuts packets on the chunk-size grid **and** at contig
//! boundaries, so every packet is single-contig — letting each worker
//! bind one [`ManualEvictChromRefFetcher`] per contig and reuse it (and
//! its raw-byte cache) across that contig's packets. Workers pull from a
//! shared MPMC queue, so a packet is processed by whichever worker is
//! free; each carries a monotonic sequence number so the consumer can
//! restore global coordinate order with a [`ReorderReadIter`] before
//! feeding the (serial) walker — preserving byte-identical output.
//!
//! Back-pressure is the bounded channel depth: the source blocks when
//! workers fall behind, workers block when the walker falls behind.

use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Mutex;

use crossbeam_channel::{Receiver, Sender};
use noodles_fasta as fasta;

use crate::bam::alignment_input::MappedRead;
use crate::bam::errors::AlignmentInputError;
use crate::baq::BaqConfig;
use crate::fasta::{ContigList, ManualEvictChromRefFetcher};
use crate::pileup::walker::PreparedRead;

use super::baq_engine::BaqEngine;
use super::baq_stream::BaqSkipCounts;
use super::read_processor::{
    DropReason, RawContigRefCache, ReadOutcome, ReadProcessingConfig, process_read,
};

/// Filter tallies the read-processing stage produces outside the
/// walker's own counters: BAQ skips (by reason) plus the G2 / F1
/// rejection counts that moved out of the reader. Accumulated per
/// packet by the workers and summed by the consumer.
#[derive(Debug, Default, Clone, Copy)]
pub struct StageCounts {
    pub baq_skips: BaqSkipCounts,
    pub bad_cigar: u64,
    pub high_mismatch: u64,
}

impl StageCounts {
    fn record_drop(&mut self, reason: DropReason) {
        match reason {
            DropReason::Baq(r) => self.baq_skips.bump(r),
            DropReason::BadCigar => self.bad_cigar += 1,
            DropReason::HighMismatchFraction => self.high_mismatch += 1,
        }
    }

    fn merge(&mut self, other: &StageCounts) {
        self.baq_skips.merge(&other.baq_skips);
        self.bad_cigar += other.bad_cigar;
        self.high_mismatch += other.high_mismatch;
    }
}

/// One single-contig chunk of reads handed from the source to a worker.
/// Non-empty and homogeneous in `ref_id` by construction (see
/// [`produce_packets`]).
pub struct RawPacket {
    seq: u64,
    reads: Vec<MappedRead>,
}

/// A worker's output for one [`RawPacket`]: the surviving prepared
/// reads (in input order) plus the filter tallies for that packet,
/// stamped with the packet's sequence number for reordering.
pub struct ResultPacket {
    seq: u64,
    survivors: Vec<PreparedRead>,
    counts: StageCounts,
}

/// Source role: pull coordinate-sorted reads from `reader`, pack them
/// into single-contig chunks of up to `chunk_size`, and send each with
/// a monotonic sequence number. An upstream read error is stashed into
/// `err_cell` (surfaced by the caller after the survivors drain, matching
/// the old error-shedding adapter) and ends production; reads already
/// buffered before the error are flushed first.
///
/// Returns early (dropping `tx`) if the consumer has gone away.
pub fn produce_packets<R>(
    reader: &mut R,
    tx: Sender<RawPacket>,
    err_cell: &Mutex<Option<AlignmentInputError>>,
    chunk_size: usize,
) where
    R: Iterator<Item = Result<MappedRead, AlignmentInputError>>,
{
    let mut seq = 0u64;
    let mut buf: Vec<MappedRead> = Vec::with_capacity(chunk_size);
    // A read pulled across a contig boundary; seeds the next packet.
    let mut pending: Option<MappedRead> = None;
    let mut upstream_done = false;

    while !upstream_done || pending.is_some() || !buf.is_empty() {
        if let Some(r) = pending.take() {
            buf.push(r);
        }
        while buf.len() < chunk_size && !upstream_done {
            match reader.next() {
                Some(Ok(r)) => {
                    if buf.first().is_some_and(|first| first.ref_id != r.ref_id) {
                        pending = Some(r);
                        break;
                    }
                    buf.push(r);
                }
                Some(Err(e)) => {
                    *err_cell.lock().expect("read-pipeline err_cell poisoned") = Some(e);
                    upstream_done = true;
                }
                None => upstream_done = true,
            }
        }
        if !buf.is_empty() {
            let packet = RawPacket {
                seq,
                reads: std::mem::take(&mut buf),
            };
            seq += 1;
            if tx.send(packet).is_err() {
                return; // consumer dropped; stop producing.
            }
            buf = Vec::with_capacity(chunk_size);
        }
    }
}

/// Worker role: process every read of each packet through the full
/// [`process_read`] fold and forward the survivors + tallies. Owns
/// per-contig BAQ scratch (an engine + uppercased-window fetcher, built
/// only when BAQ is enabled) and a persistent raw-byte reference cache,
/// all rebound when a packet's contig changes — which, because the
/// source emits packets in contig order, happens at most once per
/// contig per worker.
///
/// Returns when the packet channel disconnects (source done) or the
/// consumer has gone away.
#[allow(clippy::too_many_arguments)]
pub fn run_worker(
    rx: Receiver<RawPacket>,
    tx: Sender<ResultPacket>,
    baq_cfg: BaqConfig,
    proc_cfg: ReadProcessingConfig,
    apply_baq: bool,
    repository: &fasta::Repository,
    fasta_path: &Path,
    contigs: &ContigList,
) {
    let mut current_contig: Option<usize> = None;
    let mut baq: Option<(BaqEngine, ManualEvictChromRefFetcher)> = None;
    let mut raw_ref = RawContigRefCache::new(repository.clone(), contigs.clone());

    for packet in rx {
        // Packets are non-empty and single-contig by construction.
        let ref_id = packet.reads[0].ref_id;
        if current_contig != Some(ref_id) {
            current_contig = Some(ref_id);
            baq = apply_baq.then(|| {
                let name = &contigs.entries[ref_id].name;
                let fetcher = ManualEvictChromRefFetcher::for_contig(fasta_path, name).expect(
                    "read-pipeline per-worker fetcher construction failed; FASTA path + contig \
                     name were validated upstream",
                );
                (BaqEngine::new(baq_cfg), fetcher)
            });
            // `raw_ref` refreshes itself on the new `ref_id` lazily.
        }

        let mut survivors = Vec::with_capacity(packet.reads.len());
        let mut counts = StageCounts::default();
        for read in packet.reads {
            let pos = read.pos;
            let baq_args = baq.as_mut().map(|(engine, fetcher)| (engine, fetcher));
            match process_read(read, baq_args, &mut raw_ref, &proc_cfg) {
                ReadOutcome::Prepared(prepared) => survivors.push(prepared),
                ReadOutcome::Dropped(reason) => counts.record_drop(reason),
            }
            if let Some((_, fetcher)) = baq.as_mut()
                && let Ok(p) = u32::try_from(pos)
            {
                fetcher.evict_before(p);
            }
        }

        let result = ResultPacket {
            seq: packet.seq,
            survivors,
            counts,
        };
        if tx.send(result).is_err() {
            return; // consumer dropped.
        }
    }
}

/// Consumer-side adapter: receives [`ResultPacket`]s (arriving in
/// whatever order the workers finish) and yields their surviving
/// `PreparedRead`s in strict packet-sequence order, so the walker sees a
/// coordinate-sorted stream. Out-of-order packets are held in a
/// `BTreeMap` reorder buffer until the next expected sequence arrives.
/// Stage counts are accumulated as packets are consumed and are
/// available via [`Self::counts`] once the stream is drained.
pub struct ReorderReadIter {
    rx: Receiver<ResultPacket>,
    buffer: BTreeMap<u64, std::vec::IntoIter<PreparedRead>>,
    next_seq: u64,
    current: Option<std::vec::IntoIter<PreparedRead>>,
    counts: StageCounts,
    disconnected: bool,
}

impl ReorderReadIter {
    pub fn new(rx: Receiver<ResultPacket>) -> Self {
        Self {
            rx,
            buffer: BTreeMap::new(),
            next_seq: 0,
            current: None,
            counts: StageCounts::default(),
            disconnected: false,
        }
    }

    /// Stage counts summed across every packet consumed so far. Final
    /// once the iterator has been fully drained by the walker.
    pub fn counts(&self) -> StageCounts {
        self.counts
    }
}

impl Iterator for ReorderReadIter {
    type Item = PreparedRead;

    fn next(&mut self) -> Option<PreparedRead> {
        loop {
            if let Some(it) = self.current.as_mut() {
                if let Some(read) = it.next() {
                    return Some(read);
                }
                self.current = None;
                self.next_seq += 1;
            }
            if let Some(packet_reads) = self.buffer.remove(&self.next_seq) {
                self.current = Some(packet_reads);
                continue;
            }
            if self.disconnected {
                return None;
            }
            match self.rx.recv() {
                Ok(packet) => {
                    self.counts.merge(&packet.counts);
                    self.buffer.insert(packet.seq, packet.survivors.into_iter());
                }
                Err(_) => self.disconnected = true,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::walker::{CigarOp, MateRole};
    use std::sync::Arc;

    fn dummy_prepared(start: u32) -> PreparedRead {
        PreparedRead {
            chrom_id: 0,
            alignment_start: start,
            alignment_end: start,
            cigar: vec![CigarOp::Match(1)],
            seq: vec![b'A'],
            bq_baq: vec![30],
            mq_log_err: 0.0,
            mapq: 60,
            is_reverse_strand: false,
            qname: Arc::from("r"),
            mate_role: MateRole::Solo,
            adaptor_boundary: None,
        }
    }

    fn result_packet(seq: u64, starts: &[u32], bad_cigar: u64) -> ResultPacket {
        ResultPacket {
            seq,
            survivors: starts.iter().map(|&s| dummy_prepared(s)).collect(),
            counts: StageCounts {
                bad_cigar,
                ..Default::default()
            },
        }
    }

    #[test]
    fn reorders_out_of_order_packets_into_sequence() {
        let (tx, rx) = crossbeam_channel::unbounded();
        // Send packets out of order: 2, 0, 1.
        tx.send(result_packet(2, &[20, 21], 0)).unwrap();
        tx.send(result_packet(0, &[1, 2], 1)).unwrap();
        tx.send(result_packet(1, &[10], 0)).unwrap();
        drop(tx);

        let mut iter = ReorderReadIter::new(rx);
        let starts: Vec<u32> = (&mut iter).map(|p| p.alignment_start).collect();
        assert_eq!(starts, vec![1, 2, 10, 20, 21]);
        assert_eq!(iter.counts().bad_cigar, 1);
    }

    #[test]
    fn empty_survivor_packets_are_skipped_but_counted() {
        let (tx, rx) = crossbeam_channel::unbounded();
        tx.send(result_packet(0, &[], 3)).unwrap(); // all dropped
        tx.send(result_packet(1, &[5], 0)).unwrap();
        drop(tx);

        let mut iter = ReorderReadIter::new(rx);
        let starts: Vec<u32> = (&mut iter).map(|p| p.alignment_start).collect();
        assert_eq!(starts, vec![5]);
        assert_eq!(iter.counts().bad_cigar, 3);
    }

    #[test]
    fn empty_stream_yields_nothing() {
        let (tx, rx) = crossbeam_channel::unbounded::<ResultPacket>();
        drop(tx);
        let mut iter = ReorderReadIter::new(rx);
        assert!(iter.next().is_none());
        assert_eq!(iter.counts().bad_cigar, 0);
    }
}
