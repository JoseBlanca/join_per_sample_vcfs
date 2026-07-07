//! Dump per-locus observed tract sequences from one or more `.ssr.psp` files for the
//! Phase-2 P2.0 slip-rate measurement: one TSV row per (sample, locus, observed sequence),
//! `sample \t chrom \t start \t seq \t count`. `start` is the record start stored in the psp
//! (== VCF POS for these loci), so it joins directly against the cohort VCF. The observed
//! sequence is the extracted repeat tract (ACGTN); joined with the VCF's called genotypes it
//! lets an offline analysis attribute each carrier's reads to its called allele and read out
//! the per-allele slip distribution (pure vs impure at fixed length).
//!
//!   cargo run --release --example ssr_slip_dump -- <a.ssr.psp> [b.ssr.psp ...] > reads.tsv

use pop_var_caller::psp::{PspReader, SsrKind};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: ssr_slip_dump <a.ssr.psp> [b.ssr.psp ...]");
        std::process::exit(2);
    }
    let stdout = std::io::stdout();
    let mut w = BufWriter::new(stdout.lock());
    writeln!(w, "sample\tchrom\tstart\tseq\tcount").unwrap();
    for path in &args[1..] {
        let sample = Path::new(path)
            .file_name()
            .and_then(|s| s.to_str())
            .and_then(|s| s.split('.').next())
            .unwrap_or("?");
        let reader = PspReader::new(File::open(path).expect("open psp")).expect("psp header");
        let chroms: Vec<String> = reader
            .header()
            .chromosomes
            .iter()
            .map(|c| c.name.clone())
            .collect();
        for rec in reader.into_records_of::<SsrKind>() {
            let r = rec.expect("decode ssr record");
            let chrom = &chroms[r.chrom_id as usize];
            for (seq, count) in &r.observed {
                let seq = std::str::from_utf8(seq).expect("ACGTN tract");
                writeln!(w, "{sample}\t{chrom}\t{}\t{seq}\t{count}", r.start).unwrap();
            }
        }
    }
}
