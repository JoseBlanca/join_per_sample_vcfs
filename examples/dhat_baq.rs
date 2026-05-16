//! Heap-profile a single `BaqEngine` pass over the bench's
//! `baq_engine_read_length/150` workload.
//!
//! Build and run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example dhat_baq --features dhat-heap
//! ```
//!
//! Produces `dhat-heap.json`. Open it at
//! <https://nnethercote.github.io/dh_view/dh_view.html> to see
//! allocation sites ranked by total bytes / blocks / lifetimes.
//!
//! Used by the section-3 measurement plan in
//! `ia/reviews/perf_baq_2026-05-12.md` to confirm or refute the
//! per-read allocation findings L1–L5.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::io;

use merge_per_sample_vcfs::per_sample_pileup::baq::{BaqConfig, BaqEngine, BaqOutcome};
use merge_per_sample_vcfs::per_sample_pileup::cram_input::{CigarOp, MappedRead};
use merge_per_sample_vcfs::per_sample_pileup::pileup::RefSeqFetcher;

/// Constant-base reference; mirrors the bench's `ConstRefFetcher`.
struct ConstRefFetcher {
    len: usize,
}

impl RefSeqFetcher for ConstRefFetcher {
    fn fetch(&self, _chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        let lo = (start_1based - 1) as usize;
        let hi = lo + length as usize;
        if hi > self.len {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "off end"));
        }
        Ok(vec![b'A'; length as usize])
    }
}

fn build_mapped_reads(read_len: u32, span: u32, coverage: u32) -> Vec<MappedRead> {
    let num_reads = ((span as u64) * (coverage as u64) / (read_len as u64)).max(1) as u32;
    let last_start = span - read_len + 1;
    let mut reads = Vec::with_capacity(num_reads as usize);
    for i in 0..num_reads {
        let start = 1 + ((i as u64) * (last_start as u64 - 1) / num_reads.max(1) as u64) as u32;
        reads.push(MappedRead {
            qname: format!("r{i}").into_bytes(),
            flag: 0,
            ref_id: 0,
            pos: start as u64,
            mapq: 60,
            cigar: vec![CigarOp::Match(read_len)],
            seq: vec![b'A'; read_len as usize],
            qual: vec![30u8; read_len as usize],
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            source_file_index: 0,
        });
    }
    reads.sort_by_key(|r| r.pos);
    reads
}

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let span: u32 = 50_000;
    let coverage: u32 = 30;
    let read_len: u32 = 150;
    let fetcher = ConstRefFetcher {
        len: span as usize + 200,
    };

    let reads = build_mapped_reads(read_len, span, coverage);

    let mut engine = BaqEngine::new(BaqConfig::default());
    let mut capped: u64 = 0;
    let n = reads.len();
    for r in reads {
        if let BaqOutcome::Capped(_) = engine.process(r, &fetcher) {
            capped += 1;
        }
    }
    eprintln!("reads processed: {n} (capped: {capped})");
}
