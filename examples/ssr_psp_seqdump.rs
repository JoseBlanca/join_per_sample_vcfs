//! Pool the distinct observed tract sequences at one catalog locus across many
//! `.ssr.psp` files, to see same-length interior variation that a length-keyed
//! ladder would collapse.
//!
//!   cargo run --release --example ssr_psp_seqdump -- <start> <psp1> [psp2 ...]
//!
//! `start` is the record `start` stored in the psp (== VCF POS for these loci).

use pop_var_caller::psp::{PspReader, SsrKind};
use std::collections::BTreeMap;
use std::fs::File;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: ssr_psp_seqdump <start> <psp1> [psp2 ...]");
        std::process::exit(2);
    }
    let target: u32 = args[1].parse().expect("start");
    // seq -> (total_count, n_samples_with_it)
    let mut pooled: BTreeMap<Vec<u8>, (u32, u32)> = BTreeMap::new();
    let mut n_samples_present = 0u32;
    for path in &args[2..] {
        let reader = PspReader::new(File::open(path).expect("open psp")).expect("psp header");
        for rec in reader.into_records_of::<SsrKind>() {
            let r = rec.expect("decode");
            if r.start != target {
                continue;
            }
            if !r.observed.is_empty() {
                n_samples_present += 1;
            }
            for (seq, c) in &r.observed {
                let e = pooled.entry(seq.to_vec()).or_insert((0, 0));
                e.0 += *c;
                e.1 += 1;
            }
        }
    }
    let mut rows: Vec<(&Vec<u8>, (u32, u32))> = pooled.iter().map(|(s, v)| (s, *v)).collect();
    rows.sort_by_key(|row| std::cmp::Reverse(row.1.0));
    println!(
        "# start={target} n_samples_with_reads={n_samples_present} n_distinct_seq={}",
        rows.len()
    );
    println!("len\ttotal_count\tn_samples\tseq");
    for (seq, (count, nsamp)) in rows {
        println!(
            "{}\t{}\t{}\t{}",
            seq.len(),
            count,
            nsamp,
            String::from_utf8_lossy(seq)
        );
    }
}
