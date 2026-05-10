//! Heap-profile a single pileup walker run.
//!
//! Build and run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example dhat_pileup --features dhat-heap
//! ```
//!
//! Produces `dhat-heap.json` in the project working directory. Open it
//! at <https://nnethercote.github.io/dh_view/dh_view.html> to see
//! allocation sites ranked by total bytes / blocks / lifetimes.
//!
//! Fixture matches `pileup_walker_multi_op/L=5000` from
//! `benches/pileup_walker_scaling.rs` — the most allocation-pressured
//! configuration in the existing criterion baselines and the one called
//! out as the measurement target by findings L6–L12 in
//! `ia/reviews/perf_pileup_2026-05-10.md`.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::io;
use std::sync::Arc;
use std::sync::mpsc;
use std::thread;

use merge_per_sample_vcfs::per_sample_caller::pileup::{
    CigarOp, PileupRecord, PreparedRead, RefBaseFetcher, WalkerConfig, run,
};

struct ConstFasta {
    len: usize,
}

impl RefBaseFetcher for ConstFasta {
    fn fetch(&self, _chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        let lo = (start_1based - 1) as usize;
        let hi = lo + length as usize;
        if hi > self.len {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "off end"));
        }
        Ok(vec![b'A'; length as usize])
    }
}

fn build_multi_op_cigar(ref_span: u32, cycle_len: u32, ins_len: u32) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut consumed = 0u32;
    while consumed + cycle_len + 1 <= ref_span {
        ops.push(CigarOp::Match(cycle_len));
        ops.push(CigarOp::Insertion(ins_len));
        consumed += cycle_len;
    }
    let remainder = ref_span - consumed;
    if remainder > 0 {
        ops.push(CigarOp::Match(remainder));
    }
    ops
}

fn read_seq_len_for_cigar(cigar: &[CigarOp]) -> u32 {
    cigar
        .iter()
        .map(|op| match *op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Insertion(n)
            | CigarOp::SoftClip(n) => n,
            _ => 0,
        })
        .sum()
}

fn build_reads(read_len: u32, span: u32, coverage: u32) -> Vec<PreparedRead> {
    let num_reads = ((span as u64) * (coverage as u64) / (read_len as u64)) as u32;
    let last_start = span - read_len + 1;
    let cigar = build_multi_op_cigar(read_len, 50, 2);
    let read_consumed = read_seq_len_for_cigar(&cigar);
    let seq = vec![b'A'; read_consumed as usize];
    let bq_baq = vec![30u8; read_consumed as usize];

    let mut reads = Vec::with_capacity(num_reads as usize);
    for i in 0..num_reads {
        let start = 1 + ((i as u64) * (last_start as u64 - 1) / num_reads.max(1) as u64) as u32;
        reads.push(PreparedRead {
            chrom_id: 0,
            alignment_start: start,
            alignment_end: start + read_len - 1,
            cigar: cigar.clone(),
            seq: seq.clone(),
            bq_baq: bq_baq.clone(),
            mq_log_err: -3.0,
            is_reverse_strand: false,
            qname: Arc::from(format!("r{i}").as_str()),
            is_first_mate: true,
            has_mate: false,
            adaptor_boundary: None,
        });
    }
    reads.sort_by_key(|r| r.alignment_start);
    reads
}

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let span: u32 = 50_000;
    let coverage: u32 = 30;
    let read_len: u32 = 5_000;
    let fasta = ConstFasta {
        len: span as usize + 100,
    };

    let reads = build_reads(read_len, span, coverage);

    let (tx, rx) = mpsc::sync_channel::<PileupRecord>(64);
    let collector = thread::spawn(move || {
        let mut count = 0u64;
        while rx.recv().is_ok() {
            count += 1;
        }
        count
    });
    run(reads, &fasta, &tx, &WalkerConfig::default()).expect("walker run failed");
    drop(tx);
    let count = collector.join().expect("collector panicked");
    eprintln!("pileup records emitted: {count}");
}
