//! Dump the rough-caller het baseline from a set of `.psp` files.
//!
//! For each `.psp` given on the command line it reads the `SampleSummary`
//! metadata (the Stage-1 rough-caller output) and prints the per-sample
//! het counts, the `obs_het` rate (`n_het_sites / callable_positions`), and
//! both inbreeding-coefficient flavours:
//!   - `F_ratio` — the het:hom-alt formula the SFS prior feeds on
//!     (`DiversityEstimate::inbreeding_coefficients`);
//!   - `F_rate`  — `clip(1 − obs_het / Hexp, 0, 0.99)`, the rate/Hexp formula
//!     from `paralog::inbreeding`, using the cohort `Hexp` estimated here.
//!
//! Usage: `cargo run --release --example het_baseline -- <a.psp> <b.psp> ...`
//! This is a throwaway measurement harness for the rough-SNP-heuristics work,
//! not a shipped tool.

use std::fs::File;
use std::path::Path;

use pop_var_caller::paralog::inbreeding::{inbreeding_coefficient, obs_het};
use pop_var_caller::psp::PspReader;
use pop_var_caller::sample_summary::SampleSummary;
use pop_var_caller::var_calling::diversity::DiversityEstimate;

// Species-range fallback θ; only used if the cohort has too few alt copies to
// self-estimate Hexp (it won't here). Value is inert for a real cohort.
const PRIOR_THETA: f64 = 0.01;

fn main() {
    let paths: Vec<String> = std::env::args().skip(1).collect();
    if paths.is_empty() {
        eprintln!("usage: het_baseline <a.psp> <b.psp> ...");
        std::process::exit(2);
    }

    let mut names = Vec::new();
    let mut summaries = Vec::new();
    for p in &paths {
        let f = File::open(p).unwrap_or_else(|e| panic!("open {p}: {e}"));
        let reader = PspReader::new(f).unwrap_or_else(|e| panic!("read {p}: {e}"));
        let meta = reader
            .metadata()
            .unwrap_or_else(|| panic!("{p}: no metadata section"));
        let summary = SampleSummary::from_toml_bytes(meta)
            .unwrap_or_else(|e| panic!("{p}: parse summary: {e}"));
        names.push(
            Path::new(p)
                .file_stem()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| p.clone()),
        );
        summaries.push(summary);
    }

    // Cohort Hexp (nucleotide diversity) + the ratio-based per-sample F.
    let diversity = DiversityEstimate::from_summaries(&summaries, PRIOR_THETA, None);
    let hexp = diversity.nucleotide_diversity;

    println!(
        "# {} samples | Hexp(theta)={:.6} src={:?}",
        summaries.len(),
        hexp,
        diversity.source
    );
    println!(
        "{:<28}{:>10}{:>9}{:>9}{:>9}{:>12}{:>11}{:>10}{:>9}",
        "sample", "callable", "n_het", "n_homalt", "n_ambig", "obs_het/kb", "obs_het", "F_ratio",
        "F_rate"
    );

    let mut tot_callable: u64 = 0;
    let mut tot_het: u64 = 0;
    let mut tot_homalt: u64 = 0;
    let mut tot_ambig: u64 = 0;
    for (i, s) in summaries.iter().enumerate() {
        let callable = s.coverage_by_gc.callable_positions;
        let h = &s.heterozygosity;
        let oh = obs_het(h, callable);
        let f_ratio = diversity.inbreeding_coefficients[i];
        let f_rate = inbreeding_coefficient(oh, hexp);
        println!(
            "{:<28}{:>10}{:>9}{:>9}{:>9}{:>12.3}{:>11.6}{:>10.4}{:>9.4}",
            names[i],
            callable,
            h.n_het_sites,
            h.n_hom_alt_sites,
            h.n_ambiguous_sites,
            oh * 1000.0,
            oh,
            f_ratio,
            f_rate,
        );
        tot_callable += callable;
        tot_het += h.n_het_sites;
        tot_homalt += h.n_hom_alt_sites;
        tot_ambig += h.n_ambiguous_sites;
    }

    let pooled_oh = if tot_callable > 0 {
        tot_het as f64 / tot_callable as f64
    } else {
        0.0
    };
    println!(
        "{:<28}{:>10}{:>9}{:>9}{:>9}{:>12.3}{:>11.6}{:>10}{:>9}",
        "POOLED",
        tot_callable,
        tot_het,
        tot_homalt,
        tot_ambig,
        pooled_oh * 1000.0,
        pooled_oh,
        "-",
        "-"
    );
}
