//! Probe the fitted single-copy coverage model for a `.psp` — the paralog
//! filter's per-sample input. Prints the fit outcome (accepted / rejected),
//! `single_copy_scale` (depth mode), σ₀, the depth median, the mode/median
//! guard ratio, and — the load-bearing number — the relative copy number the
//! model assigns to a window sitting at the *median* depth, plus how many σ₀
//! that lands from H1's mean of 1.0.
//!
//! Usage: `cargo run --release --example paralog_fit_probe -- <file.psp> ...`

use std::fs::File;
use std::io::BufReader;

use pop_var_caller::paralog::{CoverageFitConfig, SingleCopyCoverageModel};
use pop_var_caller::psp::PspReader;
use pop_var_caller::sample_summary::SampleSummary;

fn main() {
    let cfg = CoverageFitConfig::default();
    for path in std::env::args().skip(1) {
        let file = File::open(&path).expect("open psp");
        let reader = PspReader::new(BufReader::new(file)).expect("read psp");
        let sample = reader.header().sample.clone();
        let bytes = reader.metadata().expect("summary section");
        let summary = SampleSummary::from_toml_bytes(bytes).expect("parse summary");
        let cov = &summary.coverage_by_gc;

        // Depth marginal (sum over GC) and its weighted median depth.
        let depth_cols = cov.depth_bins as usize + 1;
        let mut marg = vec![0.0f64; depth_cols];
        let mut total = 0.0f64;
        for gc in 0..cov.gc_bins as usize {
            for (d, slot) in marg.iter_mut().enumerate() {
                let c = cov.counts[gc * depth_cols + d] as f64;
                *slot += c;
                total += c;
            }
        }
        let mid = |d: usize| -> f64 {
            if d < cov.depth_bins as usize {
                (d as f64 + 0.5) * cov.depth_bin_width
            } else {
                cov.depth_bins as f64 * cov.depth_bin_width
            }
        };
        let mut acc = 0.0f64;
        let mut median = 0.0f64;
        for (d, &m) in marg.iter().enumerate() {
            acc += m;
            if acc >= total / 2.0 {
                median = mid(d);
                break;
            }
        }

        println!("================ {sample}  ({path}) ================");
        println!("  median depth (marginal) : {median:.2}");
        match SingleCopyCoverageModel::fit(cov, &cfg) {
            Ok(model) => {
                let scale = model.single_copy_scale();
                let sigma0 = model.single_copy_depth_sd();
                // Relative copy number of a window at the median depth,
                // averaged over the two GC ends (0.35 / 0.55 are typical).
                let rcn_med = model.relative_copy_number(0.45, median);
                let z = (rcn_med - 1.0) / sigma0; // σ₀ from H1's mean (1.0)
                println!("  FIT ACCEPTED");
                println!("  single_copy_scale (mode): {scale:.3}");
                println!("  sigma0 (relative)       : {sigma0:.4}");
                println!("  mode/median ratio       : {:.3}", scale / median);
                println!("  rcn @ median depth      : {rcn_med:.3}");
                println!("  -> (rcn-1)/sigma0       : {z:.2}  σ₀ from H1 mean");
            }
            Err(e) => {
                println!("  FIT REJECTED: {e:?}");
            }
        }
    }
}
