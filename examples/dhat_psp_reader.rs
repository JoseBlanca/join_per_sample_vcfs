//! Heap-profile a single `PspReader` end-to-end drain.
//!
//! Build and run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example dhat_psp_reader --features dhat-heap
//! ```
//!
//! Produces `dhat-heap.json` in the project working directory. Open
//! it at <https://nnethercote.github.io/dh_view/dh_view.html> to see
//! allocation sites ranked by total bytes / blocks / lifetimes —
//! deferred since the 2026-05-13 review (L9 there) and asked for
//! again as L7 of the 2026-05-23 PSP reader review, because the H1
//! CSR-collapse finding's mechanism is "fewer allocations per block"
//! — and a sampling cycle profile can't show *how many* allocations
//! happen at a site, only how many *CPU samples* land there.
//!
//! Pre-H1, the dominant allocation sites should be:
//!   - `decode_bytes_split` per-allele `bytes[..].to_vec()`
//!     (block.rs ~521; review H1 target)
//!   - `decode_list_column` per-row `Vec::with_capacity(k)`
//!     (block.rs ~412; review L1 target)
//!   - `materialise_next_record` per-record `Vec::with_capacity(n_alleles)`
//!     (reader.rs ~665; will remain until/unless lending-iterator)
//!   - `materialise_next_record` per-emit `Vec<u8>::to_vec()` for
//!     `AlleleObservation.seq` and `chain_ids` (reader.rs ~671/678 —
//!     the same materialisation point, just moved from the `mem::take`
//!     shape to slice-and-to_vec after H1)
//!
//! Post-H1 expectation: the two block-decoder lines collapse to
//! "two allocations per column per block" (the slab + offsets,
//! reused across blocks once the scratch grows), and the
//! materialiser lines stay until/unless `AlleleObservation` migrates
//! off `Vec<u8>`.
//!
//! Workload: 500 000 records on one chromosome, **multi-allele**
//! shape (2–4 alleles per record, indel-shaped allele lengths
//! 1–10 bytes) — same shape as `multi_allele_500k` in
//! `benches/psp_reader_perf.rs`. Multi-allele is the H1-relevant
//! stress: the `allele-seqs` bytes column allocates one inner Vec
//! per allele, so allele-rich workloads make the
//! `decode_bytes_split` line dominate the dhat report.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, Write};

use pop_var_caller::per_sample_pileup::pileup::{
    AlleleObservation, AlleleSupportStats, PileupRecord,
};
use pop_var_caller::per_sample_pileup::psp::PspReader;
use pop_var_caller::per_sample_pileup::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use pop_var_caller::per_sample_pileup::psp::writer::PspWriter;

const NUM_RECORDS: usize = 500_000;

fn writer_header() -> WriterHeader {
    let mut params = BTreeMap::new();
    params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
    WriterHeader {
        format_version: (1, 0),
        sample: "dhat".to_string(),
        reference: "ref.fa".to_string(),
        created: "2026-05-23T11:00:00Z".parse().unwrap(),
        chromosomes: vec![ChromosomeEntry {
            name: "chr1".to_string(),
            length: NUM_RECORDS as u32 + 1024,
            md5: "0".repeat(32),
        }],
        writer: WriterProvenance {
            tool: "dhat".to_string(),
            version: "0.0.0".to_string(),
            subcommand: "psp-reader-dhat".to_string(),
            input_crams: vec!["a.cram".to_string()],
            input_fasta: "ref.fa".to_string(),
            parameters: params,
        },
    }
}

/// Multi-allele workload: 2-4 alleles per record, indel-shaped
/// allele lengths 1-10 bytes. Mirrors `build_multi_allele_records`
/// from `benches/psp_reader_perf.rs`.
fn build_records() -> Vec<PileupRecord> {
    let mut records = Vec::with_capacity(NUM_RECORDS);
    let bases = [b'A', b'C', b'G', b'T'];
    for i in 0..NUM_RECORDS {
        let pos = (i as u32) + 1;
        let n_alleles = 2 + (i % 3); // 2, 3, or 4 alleles per record
        let mut alleles = Vec::with_capacity(n_alleles);
        for a in 0..n_alleles {
            let seq_len = 1 + ((i + a) % 10); // 1..=10 bytes
            let seq: Vec<u8> = (0..seq_len).map(|k| bases[(i + a + k) & 3]).collect();
            alleles.push(AlleleObservation::new(
                seq,
                AlleleSupportStats::new(
                    10 + ((i as u32 + a as u32) % 20),
                    -20.0 - ((i + a) as f64).rem_euclid(10.0),
                    5 + ((i as u32 + a as u32) % 8),
                    2 + ((i as u32) % 3),
                    1 + ((i as u32 + a as u32) % 4),
                    0,
                    0,
                ),
                Vec::new(),
            ));
        }
        records.push(PileupRecord::new(0, pos, alleles));
    }
    records
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Build fixtures *outside* the dhat scope so the writer's +
    // fixture-construction allocator pressure does not pollute the
    // reader-side report.
    let records = build_records();
    let header = writer_header();

    // Serialise to a Vec<u8> outside the timed scope, then write to
    // a tempfile so the reader exercises the production
    // `BufReader<File>` shape (not `Cursor<&[u8]>`).
    let sink: Vec<u8> = Vec::with_capacity(8 * 1024 * 1024);
    let mut writer = PspWriter::new(sink, header)?;
    for r in &records {
        writer.write_record(r)?;
    }
    let serialised = writer.finish()?;
    eprintln!(
        "dhat psp reader: serialised {} records into {} bytes on disk",
        records.len(),
        serialised.len()
    );
    drop(records);

    let mut tempfile = tempfile::NamedTempFile::new()?;
    tempfile.write_all(&serialised)?;
    tempfile.flush()?;
    drop(serialised);
    let path = tempfile.path().to_path_buf();

    // ---------- start measurement ----------
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let file = File::open(&path)?;
    let buf = BufReader::with_capacity(64 * 1024, file);
    let mut reader = PspReader::new(buf)?;
    let mut count: u64 = 0;
    for r in reader.records() {
        let r = r?;
        // Force the materialised AlleleObservation Vecs to be alive
        // at dhat capture time so they show up in the live-bytes
        // total (they would otherwise be optimised out).
        std::hint::black_box(&r);
        count += 1;
    }
    // Profiler drop at end of main writes dhat-heap.json.
    eprintln!("dhat psp reader: drained {} records", count);
    assert_eq!(count, NUM_RECORDS as u64);
    Ok(())
}
