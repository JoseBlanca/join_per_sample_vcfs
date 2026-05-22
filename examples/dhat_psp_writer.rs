//! Heap-profile a single `PspWriter` run.
//!
//! Build and run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example dhat_psp_writer --features dhat-heap
//! ```
//!
//! Produces `dhat-heap.json` in the project working directory. Open
//! it at <https://nnethercote.github.io/dh_view/dh_view.html> to see
//! allocation sites ranked by total bytes / blocks / lifetimes — the
//! primary signal for the H1 / H2 / L1 / L2 wins in
//! `ia/reviews/perf_psp_writer_2026-05-13.md`.
//!
//! Workload: 1 000 000 records on one chromosome, SNP-typical
//! single-byte alleles. Sized so the writer triggers ≥ 2 block
//! flushes at the 16 MiB target (a full block + a partial), keeping
//! both `BlockAccumulator::append_record` and `flush_block` on the
//! measured path.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::io;

use pop_var_caller::per_sample_pileup::pileup::{
    AlleleObservation, AlleleSupportStats, PileupRecord,
};
use pop_var_caller::per_sample_pileup::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use pop_var_caller::per_sample_pileup::psp::writer::PspWriter;
use std::collections::BTreeMap;

const NUM_RECORDS: usize = 1_000_000;
const SECOND_ALLELE_PER_1000: u32 = 1;

fn writer_header() -> WriterHeader {
    let mut params = BTreeMap::new();
    params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
    WriterHeader {
        format_version: (1, 0),
        sample: "dhat".to_string(),
        reference: "ref.fa".to_string(),
        created: "2026-05-13T11:00:00Z".parse().unwrap(),
        chromosomes: vec![ChromosomeEntry {
            name: "chr1".to_string(),
            length: NUM_RECORDS as u32 + 1024,
            md5: "0".repeat(32),
        }],
        writer: WriterProvenance {
            tool: "dhat".to_string(),
            version: "0.0.0".to_string(),
            subcommand: "psp-writer-dhat".to_string(),
            input_crams: vec!["a.cram".to_string()],
            input_fasta: "ref.fa".to_string(),
            parameters: params,
        },
    }
}

fn build_records() -> Vec<PileupRecord> {
    let mut records = Vec::with_capacity(NUM_RECORDS);
    let bases = [b'A', b'C', b'G', b'T'];
    for i in 0..NUM_RECORDS {
        let pos = (i as u32) + 1;
        let ref_base = bases[i & 3];
        let two_alleles = (i as u32 % 1000) < SECOND_ALLELE_PER_1000;
        let mut alleles = Vec::with_capacity(if two_alleles { 2 } else { 1 });
        alleles.push(AlleleObservation::new(
            vec![ref_base],
            AlleleSupportStats::new(
                28 + ((i as u32) % 4),
                -42.0 - (i as f64).rem_euclid(7.0),
                14 + ((i as u32) % 3),
                5 + ((i as u32) % 3),
                1 + ((i as u32) % 3),
                            0,
                0,
),
            Vec::new(),
        ));
        if two_alleles {
            let alt_base = bases[(i + 1) & 3];
            alleles.push(AlleleObservation::new(
                vec![alt_base],
                AlleleSupportStats::new(
                    3 + ((i as u32) % 5),
                    -3.5 - (i as f64).rem_euclid(2.0),
                    1 + ((i as u32) % 2),
                    1,
                    0,
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

fn main() {
    // Build fixtures *outside* the dhat scope so the writer's
    // allocator pressure is isolated from per-record fixture
    // construction (1 M records × ≥ 2 small `Vec::new()`s per record
    // would otherwise dominate the report).
    let records = build_records();
    let header = writer_header();

    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let mut writer = PspWriter::new(io::sink(), header).expect("writer new");
    for r in &records {
        writer.write_record(r).expect("write_record");
    }
    writer.finish().expect("finish");

    eprintln!("psp writer records emitted: {}", records.len());
}
