//! Print the per-sample summary (coverage-by-GC + observed het) stored in
//! a `.psp` metadata section — the hidden-paralog filter input.
//!
//! This is the canonical **consumer** of the summary section: it opens the
//! `.psp` with the schema-agnostic [`PspReader`], pulls the opaque metadata
//! bytes ([`PspReader::metadata`]), and parses them with the SNP-kind
//! [`SampleSummary::from_toml_bytes`]. No new reader API is needed — the
//! container stays schema-agnostic and the kind interprets the bytes.
//!
//! Usage: `cargo run --release --example dump_sample_summary -- <file.psp>`

use std::fs::File;
use std::io::BufReader;

use pop_var_caller::pileup::per_sample::pileup_to_psp::SampleSummaryAccumulators;
use pop_var_caller::psp::PspReader;
use pop_var_caller::sample_summary::coverage::CoverageBinScheme;
use pop_var_caller::sample_summary::het::HetClassifyParams;
use pop_var_caller::sample_summary::{
    DEFAULT_DEPTH_BIN_WIDTH, DEFAULT_DEPTH_BINS, DEFAULT_GC_BINS, DEFAULT_GC_WINDOW_BP,
    DEFAULT_HET_ERROR_RATE, DEFAULT_HET_LR_MARGIN, DEFAULT_HET_MIN_DEPTH, SampleSummary,
};

fn main() {
    let path = match std::env::args().nth(1) {
        Some(p) => p,
        None => {
            eprintln!("usage: dump_sample_summary <file.psp>");
            std::process::exit(2);
        }
    };

    let file = File::open(&path).unwrap_or_else(|e| {
        eprintln!("open {path}: {e}");
        std::process::exit(1);
    });
    let mut reader = PspReader::new(BufReader::new(file)).unwrap_or_else(|e| {
        eprintln!("open psp {path}: {e}");
        std::process::exit(1);
    });

    println!("sample: {}", reader.header().sample);
    // Prefer the stored summary section (new producer). For an older `.psp`
    // without one, re-derive it from the body — the coverage histogram is a
    // pure function of the per-position records (Premise 1), which is what
    // makes the section a cache rather than new information. This re-derive
    // path is also the empirical-parity (D2) tool.
    let (summary, source) = match reader.metadata() {
        Some(bytes) => {
            let s = SampleSummary::from_toml_bytes(bytes).unwrap_or_else(|e| {
                eprintln!("parse summary: {e}");
                std::process::exit(1);
            });
            (s, "stored section")
        }
        None => {
            let mut acc = SampleSummaryAccumulators::new(
                CoverageBinScheme {
                    window_bp: DEFAULT_GC_WINDOW_BP,
                    gc_bins: DEFAULT_GC_BINS,
                    depth_bin_width: DEFAULT_DEPTH_BIN_WIDTH,
                    depth_bins: DEFAULT_DEPTH_BINS,
                },
                HetClassifyParams {
                    min_depth: DEFAULT_HET_MIN_DEPTH,
                    error_rate: DEFAULT_HET_ERROR_RATE,
                    lr_margin: DEFAULT_HET_LR_MARGIN,
                },
            );
            for r in reader.records() {
                let record = r.unwrap_or_else(|e| {
                    eprintln!("decode record: {e}");
                    std::process::exit(1);
                });
                acc.observe_record(&record);
            }
            (acc.finish(), "re-derived from body")
        }
    };
    println!("summary source: {source}");

    let cov = &summary.coverage_by_gc;
    let het = &summary.heterozygosity;

    // Coverage: derive the overall covered-bases mean depth from the
    // histogram (Σ count·bin-midpoint ÷ Σ count) as a single-number sanity
    // figure against the prototype's per-window depth table.
    let depth_cols = cov.depth_bins as usize + 1;
    let mut total = 0.0_f64;
    let mut weighted = 0.0_f64;
    for gc in 0..cov.gc_bins as usize {
        for d in 0..depth_cols {
            let c = cov.counts[gc * depth_cols + d] as f64;
            // Bin midpoint; the overflow column uses its lower edge.
            let mid = if d < cov.depth_bins as usize {
                (d as f64 + 0.5) * cov.depth_bin_width
            } else {
                cov.depth_bins as f64 * cov.depth_bin_width
            };
            total += c;
            weighted += c * mid;
        }
    }
    let mean_depth = if total > 0.0 { weighted / total } else { 0.0 };

    println!("coverage-by-gc:");
    println!("  window_bp      : {}", cov.window_bp);
    println!("  gc_bins        : {}", cov.gc_bins);
    println!(
        "  depth bins     : {} regular (width {}) + 1 overflow",
        cov.depth_bins, cov.depth_bin_width
    );
    println!("  n_tiles        : {}", cov.n_tiles);
    println!("  n_skipped_tiles: {}", cov.n_skipped_tiles);
    println!("  est mean depth : {mean_depth:.3}  (histogram-weighted)");

    println!("heterozygosity:");
    println!("  n_variant_sites : {}", het.n_variant_sites);
    println!("  n_het_sites     : {}", het.n_het_sites);
    println!("  n_hom_alt_sites : {}", het.n_hom_alt_sites);
    println!("  n_ambiguous     : {}", het.n_ambiguous_sites);
    let confident = het.n_het_sites + het.n_hom_alt_sites;
    if confident > 0 {
        println!(
            "  Hobs (confident): {:.4}",
            het.n_het_sites as f64 / confident as f64
        );
    }
    println!(
        "  params          : min_depth {} eps {} margin {:.4}",
        het.min_depth, het.error_rate, het.lr_margin
    );
}
