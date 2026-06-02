//! Read-only PSP block-index stats — calibration for the
//! genomic-window block-cut (scaling work, 2026-06-02).
//!
//! For each `.psp` given on the command line, reads only the block
//! index (no payload decode) and reports, per file: block count, and
//! the distribution of genomic span per block (last minus first
//! position), records per block, and on-disk bytes per block. Then a
//! cross-file alignment check: for the most-populated chromosome,
//! how many block *boundaries* (`last_pos`) are shared across files
//! vs. unique — the direct measure of the misalignment that inflates
//! the cohort read's sync-step count.
//!
//! Run:
//!   cargo build --release --example psp_block_stats
//!   ./target/release/examples/psp_block_stats a.psp b.psp ...

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use pop_var_caller::psp::PspReader;

fn pct(sorted: &[u64], p: f64) -> u64 {
    if sorted.is_empty() {
        return 0;
    }
    let idx = ((sorted.len() as f64 - 1.0) * p).round() as usize;
    sorted[idx]
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<PathBuf> = std::env::args().skip(1).map(PathBuf::from).collect();
    if paths.is_empty() {
        eprintln!("usage: psp_block_stats <a.psp> [b.psp ...]");
        std::process::exit(2);
    }

    // Per-file boundary sets for the alignment check, keyed by chrom.
    // chrom -> file_idx -> set of last_pos boundaries.
    let mut boundaries_by_chrom: HashMap<u32, Vec<Vec<u32>>> = HashMap::new();

    for (file_idx, path) in paths.iter().enumerate() {
        let file = File::open(path)?;
        let reader = PspReader::new(BufReader::with_capacity(64 * 1024, file))?;
        let index = reader.block_index();

        let mut spans: Vec<u64> = Vec::with_capacity(index.len());
        let mut recs: Vec<u64> = Vec::with_capacity(index.len());
        let mut bytes: Vec<u64> = Vec::with_capacity(index.len());

        for (i, e) in index.iter().enumerate() {
            spans.push(u64::from(e.last_pos - e.first_pos + 1));
            recs.push(u64::from(e.n_records));
            // On-disk size = next block's offset - this one's. Skip the
            // last block (no successor offset in the index).
            if i + 1 < index.len() {
                bytes.push(index[i + 1].block_offset - e.block_offset);
            }
            boundaries_by_chrom
                .entry(e.chrom_id)
                .or_insert_with(|| vec![Vec::new(); paths.len()])[file_idx]
                .push(e.last_pos);
        }
        spans.sort_unstable();
        recs.sort_unstable();
        bytes.sort_unstable();

        let mean = |v: &[u64]| {
            if v.is_empty() {
                0
            } else {
                v.iter().sum::<u64>() / v.len() as u64
            }
        };
        println!(
            "{:<40} blocks={:>6} | span(bp) mean={:>7} p50={:>7} p10={:>7} p90={:>7} | recs/blk mean={:>6} | bytes/blk mean={:>7} (p50={})",
            path.file_name().unwrap().to_string_lossy(),
            index.len(),
            mean(&spans),
            pct(&spans, 0.5),
            pct(&spans, 0.10),
            pct(&spans, 0.90),
            mean(&recs),
            mean(&bytes),
            pct(&bytes, 0.5),
        );
    }

    // ── Cross-file alignment on the most-populated chromosome. ──
    if paths.len() > 1 {
        let (&chrom, per_file) = boundaries_by_chrom
            .iter()
            .max_by_key(|(_, v)| v.iter().map(Vec::len).sum::<usize>())
            .unwrap();
        // Count how many files share each distinct boundary.
        let mut share: HashMap<u32, usize> = HashMap::new();
        for f in per_file {
            for &b in f {
                *share.entry(b).or_insert(0) += 1;
            }
        }
        let n_files = paths.len();
        let distinct = share.len();
        let shared_by_all = share.values().filter(|&&c| c == n_files).count();
        let unique_to_one = share.values().filter(|&&c| c == 1).count();
        let avg_per_file: f64 =
            per_file.iter().map(Vec::len).sum::<usize>() as f64 / n_files as f64;
        println!("\n── alignment on chrom {chrom} across {n_files} files ──");
        println!("  avg block boundaries per file : {avg_per_file:.0}");
        println!("  distinct boundaries (union)   : {distinct}");
        println!("  shared by ALL files           : {shared_by_all}");
        println!("  unique to ONE file            : {unique_to_one}");
        println!(
            "  → union/per-file ratio        : {:.1}x  (1.0 = perfectly aligned; ~n_files = fully misaligned)",
            distinct as f64 / avg_per_file.max(1.0),
        );
    }

    Ok(())
}
