//! Print the fitted single-copy depth SD (σ₀) for each `.psp` given on the
//! command line — the per-sample coverage-model yardstick the hidden-paralog
//! score uses. One line per sample: `<file>\t<sigma0 or "unfit">`.
//!
//! Usage: `cargo run --release --example tomato2_sigma0 -- a.psp b.psp ...`
//! Ad-hoc validation helper for the M8 tomato2 run (records the real σ₀ on the
//! regenerated sliding-window `.psp`).

use std::fs::File;
use std::io::BufReader;

use pop_var_caller::paralog::{CoverageFitConfig, SingleCopyCoverageModel};
use pop_var_caller::psp::PspReader;
use pop_var_caller::sample_summary::SampleSummary;

fn main() {
    let paths: Vec<String> = std::env::args().skip(1).collect();
    if paths.is_empty() {
        eprintln!("usage: tomato2_sigma0 <file.psp> [more.psp ...]");
        std::process::exit(2);
    }
    let cfg = CoverageFitConfig::default();
    let mut fit_sigmas: Vec<f64> = Vec::new();
    println!("file\tsigma0\tsingle_copy_scale");
    for path in &paths {
        let reader = match File::open(path)
            .map_err(|e| e.to_string())
            .and_then(|f| PspReader::new(BufReader::new(f)).map_err(|e| e.to_string()))
        {
            Ok(r) => r,
            Err(e) => {
                println!("{path}\topen-error: {e}\t");
                continue;
            }
        };
        let summary = reader
            .metadata()
            .and_then(|bytes| SampleSummary::from_toml_bytes(bytes).ok());
        let Some(summary) = summary else {
            println!("{path}\tno-summary\t");
            continue;
        };
        match SingleCopyCoverageModel::fit(&summary.coverage_by_gc, &cfg) {
            Ok(model) => {
                let sigma = model.single_copy_depth_sd();
                fit_sigmas.push(sigma);
                println!("{path}\t{sigma:.6}\t{:.4}", model.single_copy_scale());
            }
            Err(e) => println!("{path}\tunfit: {e:?}\t"),
        }
    }
    if !fit_sigmas.is_empty() {
        fit_sigmas.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = fit_sigmas.len();
        let mean = fit_sigmas.iter().sum::<f64>() / n as f64;
        let median = fit_sigmas[n / 2];
        eprintln!(
            "# fit {n} sample(s): sigma0 mean {mean:.4}, median {median:.4}, min {:.4}, max {:.4}",
            fit_sigmas[0],
            fit_sigmas[n - 1],
        );
    }
}
