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

use pop_var_caller::baq::BaqConfig;
use pop_var_caller::per_sample_pileup::baq_engine::{BaqEngine, BaqOutcome};
use pop_var_caller::per_sample_pileup::cram_input::{CigarOp, MappedRead};
use pop_var_caller::per_sample_pileup::ref_fetcher::ManualEvictChromRefFetcher;

/// Build a tempfile FASTA with `length` 'A' bases and return a
/// [`ManualEvictChromRefFetcher`] over it. The BAQ engine now takes
/// the concrete fetcher type, so the heap profile uses the same
/// production code path that Stage 1 does.
fn build_fetcher(length: usize) -> (tempfile::TempDir, ManualEvictChromRefFetcher) {
    use std::io::Write as _;
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta_path = dir.path().join("ref.fa");
    let fai_path = dir.path().join("ref.fa.fai");
    let header = b">chr0\n";
    {
        let mut fa = std::fs::File::create(&fasta_path).expect("fa");
        fa.write_all(header).expect("hdr");
        for _ in 0..length {
            fa.write_all(b"A").expect("seq");
        }
        fa.write_all(b"\n").expect("nl");
    }
    {
        let mut fai = std::fs::File::create(&fai_path).expect("fai");
        writeln!(
            fai,
            "chr0\t{}\t{}\t{}\t{}",
            length,
            header.len(),
            length,
            length + 1,
        )
        .expect("fai");
    }
    let fetcher = ManualEvictChromRefFetcher::for_contig(&fasta_path, "chr0").expect("fetcher");
    (dir, fetcher)
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
    let (_dir, mut fetcher) = build_fetcher(span as usize + 200);

    let reads = build_mapped_reads(read_len, span, coverage);

    let mut engine = BaqEngine::new(BaqConfig::default());
    let mut capped: u64 = 0;
    let n = reads.len();
    for r in reads {
        let pos = r.pos;
        if let BaqOutcome::Capped(_) = engine.process(r, &mut fetcher) {
            capped += 1;
        }
        if let Ok(p) = u32::try_from(pos) {
            fetcher.evict_before(p);
        }
    }
    eprintln!("reads processed: {n} (capped: {capped})");
}
