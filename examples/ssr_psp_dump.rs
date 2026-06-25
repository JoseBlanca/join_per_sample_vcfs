//! Dump per-locus QC counts + observed sequences from a `.ssr.psp`.
//!
//! A diagnostic for the SSR benchmark: shows, per catalog locus, how many
//! reads survived each `ssr-pileup` gate (depth / n_filtered / n_low_quality
//! / n_border_off_end / n_window_truncated) and the distinct repeat-region
//! sequences tallied. Lets us tell apart "pileup dropped the reads" from
//! "genotyper called hom-ref" when our caller disagrees with HipSTR.
//!
//!   cargo run --release --example ssr_psp_dump -- <file.ssr.psp> [start_min] [start_max]
//!
//! With start_min/start_max (record `start`, the value stored in the psp),
//! only loci in that window are printed; otherwise every locus is dumped.

use pop_var_caller::psp::{PspReader, SsrKind};
use std::fs::File;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: ssr_psp_dump <file.ssr.psp> [start_min] [start_max]");
        std::process::exit(2);
    }
    let path = &args[1];
    let start_min: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);
    let start_max: u32 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(u32::MAX);

    let reader = PspReader::new(File::open(path).expect("open psp")).expect("psp header");
    println!(
        "chrom_id\tstart\tend\tdepth\tn_filtered\tmapped\tn_lowQ\tn_borderOff\tn_winTrunc\tn_obs\tobservations(len×count)"
    );
    for rec in reader.into_records_of::<SsrKind>() {
        let r = rec.expect("decode ssr record");
        if r.start < start_min || r.start > start_max {
            continue;
        }
        let mut obs: Vec<(usize, u32)> =
            r.observed.iter().map(|(seq, c)| (seq.len(), *c)).collect();
        obs.sort();
        let obs_str = obs
            .iter()
            .map(|(len, c)| format!("{len}bp×{c}"))
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.chrom_id,
            r.start,
            r.end,
            r.depth,
            r.n_filtered,
            r.mapped_reads,
            r.n_low_quality,
            r.n_border_off_end,
            r.n_window_truncated,
            r.observed.len(),
            obs_str,
        );
    }
}
