//! Data-first validation of the hidden-paralog filter (plan step R1).
//!
//! Runs the pure statistics pieces (Q2 coverage model, Q3 scorer, Q4
//! inbreeding, Q5 prior/FDR) end-to-end on the real tomato2 data — the
//! per-sample `.psp` files (coverage histogram + het counts, re-derived from
//! the body so `callable_positions` is present) and the cohort VCF (per-sample
//! `AD` and `GT`) — and reports the data-first checks the plan calls for:
//!
//! - coverage model: every sample fits (the mode/median guard passes); σ₀ ≈
//!   0.28; the single-copy peak is anchored at `relative_copy_number = 1.0`
//!   (by construction of the mode anchor);
//! - the prior π (expected ≈ 9 %) and the FDR-flagged set's LR profile;
//! - per-locus LRs + posteriors dumped to a TSV so a `uv` Python companion
//!   (`paralog_score_parity.py`) can do the loose Python↔Rust LR correlation
//!   against the prototype's `paralog_lr.parquet` (Rust has no parquet reader)
//!   and the F Spearman check.
//!
//! Mirrors `examples/dump_sample_summary.rs` (the summary consumer) but for the
//! whole model. Not a committed benchmark; it reads the operator's tomato2 data
//! off disk and writes TSVs.
//!
//! ```text
//! cargo run --release --example paralog_score_parity -- \
//!     --psp-dir benchmarks/tomato2/psp_files \
//!     --vcf     benchmarks/tomato2/results/cohort.vcf.gz \
//!     --out-lr      tmp/paralog_parity_lr.tsv \
//!     --out-samples tmp/paralog_parity_samples.tsv
//! ```

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use noodles_vcf::variant::record_buf::samples::sample::Value;
use noodles_vcf::variant::record_buf::samples::sample::value::Array;

use pop_var_caller::paralog::coverage_model::CoverageFitConfig;
use pop_var_caller::paralog::{
    EmConfig, LocusObservations, ParalogFdrCurve, ParalogLrHistogram, ParalogModelParams,
    ParalogPrior, ParalogScorePrecompute, SampleObservation, SingleCopyCoverageModel,
    inbreeding_coefficient, obs_het, score_locus_for_paralogy,
};
use pop_var_caller::psp::PspReader;
use pop_var_caller::sample_summary::coverage::{CoverageBinScheme, CoverageByGcAccumulator};
use pop_var_caller::sample_summary::het::{HetAccumulator, HetClassifyParams, SiteCounts};
use pop_var_caller::sample_summary::{
    DEFAULT_DEPTH_BIN_WIDTH, DEFAULT_DEPTH_BINS, DEFAULT_GC_BINS, DEFAULT_GC_WINDOW_BP,
    DEFAULT_HET_ERROR_RATE, DEFAULT_HET_LR_MARGIN, DEFAULT_HET_MIN_DEPTH,
    DEFAULT_HET_STRAND_BIAS_Z,
};

/// The GC / analysis window (matches the prototype `W = 500` and the
/// per-sample summary default).
const WINDOW_BP: u32 = DEFAULT_GC_WINDOW_BP;
/// FDR level for the flagged-set profile report.
const REPORT_FDR: f64 = 0.05;

struct Args {
    psp_dir: PathBuf,
    vcf: PathBuf,
    out_lr: PathBuf,
    out_samples: PathBuf,
}

fn parse_args() -> Args {
    let mut psp_dir = None;
    let mut vcf = None;
    let mut out_lr = PathBuf::from("tmp/paralog_parity_lr.tsv");
    let mut out_samples = PathBuf::from("tmp/paralog_parity_samples.tsv");
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut val = || it.next().expect("flag needs a value");
        match flag.as_str() {
            "--psp-dir" => psp_dir = Some(PathBuf::from(val())),
            "--vcf" => vcf = Some(PathBuf::from(val())),
            "--out-lr" => out_lr = PathBuf::from(val()),
            "--out-samples" => out_samples = PathBuf::from(val()),
            other => {
                eprintln!("unknown flag: {other}");
                std::process::exit(2);
            }
        }
    }
    Args {
        psp_dir: psp_dir.expect("--psp-dir required"),
        vcf: vcf.expect("--vcf required"),
        out_lr,
        out_samples,
    }
}

/// A window key: `(global chromosome id, tile index)` where
/// `tile index = (pos − 1) / WINDOW_BP`.
type WindowKey = (u32, u32);

/// Per-sample fitted state carried from the `.psp` pass into the scoring pass.
struct SampleState {
    model: SingleCopyCoverageModel,
    inbreeding: f64,
    obs_het: f64,
    callable_positions: u64,
    /// Dense per-variant-window `(gc fraction, mean depth)`, indexed by the
    /// variant-window dense id; `None` where the sample had no covered tile.
    window_depth: Vec<Option<(f32, f32)>>,
}

fn main() {
    let args = parse_args();

    // Canonical chromosome ids (assigned by first appearance in the VCF) and
    // the dense index of variant windows.
    let mut chrom_ids: HashMap<String, u32> = HashMap::new();
    let mut window_index: HashMap<WindowKey, usize> = HashMap::new();

    // ---- VCF pre-pass: collect variant windows + the cohort Σ 2pq ----
    // Hexp (the per-callable-position expected het) is `Σ2pq / callable`, so we
    // keep the *sum* here and divide by a representative callable count once the
    // per-sample `.psp` pass has produced it (arch Premise 3).
    let (twopq_sum, sample_names) = {
        let mut reader = noodles_vcf::io::reader::Builder::default()
            .build_from_path(&args.vcf)
            .expect("open vcf");
        let header = reader.read_header().expect("vcf header");
        let sample_names: Vec<String> = header.sample_names().iter().cloned().collect();

        let mut twopq_sum = 0.0f64;
        for result in reader.record_bufs(&header) {
            let record = result.expect("vcf record");
            let Some((chrom, pos)) = biallelic_snp_site(&record) else {
                continue;
            };
            let cid = intern_chrom(&mut chrom_ids, chrom);
            let key = (cid, (pos - 1) / WINDOW_BP);
            let next = window_index.len();
            window_index.entry(key).or_insert(next);

            // Cohort allele frequency from the GTs → 2pq.
            let (alt_alleles, called_alleles) = allele_counts(&record, sample_names.len());
            if called_alleles > 0 {
                let p = alt_alleles as f64 / called_alleles as f64;
                twopq_sum += 2.0 * p * (1.0 - p);
            }
        }
        (twopq_sum, sample_names)
    };
    let n_windows = window_index.len();
    eprintln!(
        "vcf pre-pass: {} samples, {n_windows} variant windows, Σ2pq = {twopq_sum:.1}",
        sample_names.len()
    );

    // The `.psp` files are named by SRR run id but their header carries the
    // SRS sample id that matches the VCF columns — map by the header sample.
    let psp_by_sample = map_psp_by_sample(&args.psp_dir);

    // ---- Per-sample `.psp` pass: fit the coverage model + per-window depth ----
    let params = ParalogModelParams::default();
    let mut states: Vec<Option<SampleState>> = Vec::with_capacity(sample_names.len());
    let mut sigma0s: Vec<f64> = Vec::new();
    let mut n_rejected = 0usize;
    for name in &sample_names {
        let Some(path) = psp_by_sample.get(name) else {
            eprintln!("  no .psp for sample {name}; skipping");
            states.push(None);
            continue;
        };
        match fit_sample(path, &mut chrom_ids, &window_index, n_windows) {
            Ok(state) => {
                sigma0s.push(state.model.single_copy_depth_sd());
                states.push(Some(state));
            }
            Err(msg) => {
                eprintln!("  sample {name}: coverage fit rejected: {msg}");
                n_rejected += 1;
                states.push(None);
            }
        }
    }

    // Hexp = Σ2pq / callable, on the per-callable-position scale that matches
    // the observed het rate. Use the median callable across fit samples as the
    // cohort reference (arch Premise 3: one cohort Hexp), then form each F.
    let mut callables: Vec<u64> = states
        .iter()
        .filter_map(|s| s.as_ref().map(|s| s.callable_positions))
        .collect();
    callables.sort_unstable();
    let callable_ref = callables
        .get(callables.len() / 2)
        .copied()
        .unwrap_or(1)
        .max(1);
    let hexp = twopq_sum / callable_ref as f64;
    for state in states.iter_mut().flatten() {
        state.inbreeding = inbreeding_coefficient(state.obs_het, hexp);
    }
    eprintln!("Hexp = Σ2pq/callable = {twopq_sum:.1}/{callable_ref} = {hexp:.6}",);

    let n_fit = states.iter().filter(|s| s.is_some()).count();
    sigma0s.sort_by(f64::total_cmp);
    let median_sigma0 = sigma0s.get(sigma0s.len() / 2).copied().unwrap_or(f64::NAN);
    eprintln!(
        "psp pass: {n_fit}/{} samples fit ({n_rejected} rejected); median σ₀ = {median_sigma0:.3} \
         (range {:.3}..{:.3})",
        sample_names.len(),
        sigma0s.first().copied().unwrap_or(f64::NAN),
        sigma0s.last().copied().unwrap_or(f64::NAN),
    );

    // Write the per-sample TSV (for the F-Spearman companion).
    write_samples_tsv(&args.out_samples, &sample_names, &states);

    // ---- VCF scoring pass: score every locus, fold LRs into the histogram ----
    let sigma0_by_sample: Vec<f64> = states
        .iter()
        .map(|s| {
            s.as_ref()
                .map(|s| s.model.single_copy_depth_sd())
                .unwrap_or(1.0)
        })
        .collect();
    // `F` is one cohort constant now (the caller's `--inbreeding-coefficient`),
    // not a per-individual `Hexp`-derived value — see the single-sample
    // reformulation. The per-sample `state.inbreeding`/`obs_het` above are kept
    // only as diagnostics in the output table, not fed to the score. The scorer's
    // log-prior tables are memoised once from (params, F) via the precompute.
    let cohort_f = pop_var_caller::var_calling::posterior_engine::DEFAULT_INBREEDING_COEFFICIENT;
    let inbreeding_by_sample: Vec<f64> = vec![cohort_f; sample_names.len()];
    let precompute = ParalogScorePrecompute::new(&params, &inbreeding_by_sample);

    let mut histogram = ParalogLrHistogram::with_defaults();
    let mut per_locus: Vec<(String, u32, f64)> = Vec::new();
    {
        let mut reader = noodles_vcf::io::reader::Builder::default()
            .build_from_path(&args.vcf)
            .expect("reopen vcf");
        let header = reader.read_header().expect("vcf header");
        let mut obs_buf: Vec<Option<SampleObservation>> = vec![None; sample_names.len()];
        for result in reader.record_bufs(&header) {
            let record = result.expect("vcf record");
            let Some((chrom, pos)) = biallelic_snp_site(&record) else {
                continue;
            };
            let Some(&cid) = chrom_ids.get(chrom) else {
                continue;
            };
            let Some(&win) = window_index.get(&(cid, (pos - 1) / WINDOW_BP)) else {
                continue;
            };

            for (col, slot) in obs_buf.iter_mut().enumerate() {
                *slot = build_observation(&record, col, &states, win);
            }
            // No minimum-sample gate — the LR self-gates (under-powered → LR≈0).
            // Only skip a locus with nothing to score (no usable sample).
            if obs_buf.iter().all(|o| o.is_none()) {
                continue;
            }
            let score = score_locus_for_paralogy(
                &LocusObservations { samples: &obs_buf },
                &sigma0_by_sample,
                &precompute,
            );
            let lr = score.paralog_log_likelihood_ratio;
            histogram.push(lr);
            per_locus.push((chrom.to_string(), pos, lr));
        }
    }
    eprintln!("scoring pass: {} loci scored", per_locus.len());

    // ---- Calibrate: EM prior + FDR curve, then stamp posteriors ----
    let prior = ParalogPrior::estimate(&histogram, &EmConfig::default());
    let curve = ParalogFdrCurve::from_histogram(&histogram, &prior);
    eprintln!(
        "prior: π = {:.4} ({:.2}% of loci; converged={}); 50% cut at LR = {:.2}",
        prior.prior_probability,
        100.0 * prior.prior_probability,
        prior.converged,
        prior.half_posterior_ratio(),
    );

    write_lr_tsv(&args.out_lr, &per_locus, &prior, &curve);
    report_flagged_profile(&per_locus, &curve);
}

/// Fit one sample: coverage histogram + het counts (→ σ₀, model, F) and the
/// per-variant-window `(gc, depth)` map, in a single pass over the `.psp` body.
fn fit_sample(
    path: &Path,
    chrom_ids: &mut HashMap<String, u32>,
    window_index: &HashMap<WindowKey, usize>,
    n_windows: usize,
) -> Result<SampleState, String> {
    let file = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let mut reader =
        PspReader::new(BufReader::with_capacity(64 * 1024, file)).map_err(|e| e.to_string())?;

    // Map this sample's chromosome ids to the canonical (VCF-assigned) ids.
    let psp_to_global: Vec<u32> = reader
        .header()
        .chromosomes
        .iter()
        .map(|c| intern_chrom(chrom_ids, &c.name))
        .collect();

    let mut coverage = CoverageByGcAccumulator::new(CoverageBinScheme {
        window_bp: WINDOW_BP,
        gc_bins: DEFAULT_GC_BINS,
        depth_bin_width: DEFAULT_DEPTH_BIN_WIDTH,
        depth_bins: DEFAULT_DEPTH_BINS,
    });
    let mut het = HetAccumulator::new(HetClassifyParams {
        min_depth: DEFAULT_HET_MIN_DEPTH,
        error_rate: DEFAULT_HET_ERROR_RATE,
        lr_margin: DEFAULT_HET_LR_MARGIN,
        strand_bias_z: DEFAULT_HET_STRAND_BIAS_Z,
    });
    let mut windows = WindowDepthCollector::new(window_index, n_windows);

    for r in reader.records() {
        let record = r.map_err(|e| e.to_string())?;
        let Some(ref_allele) = record.alleles.first() else {
            continue;
        };
        let ref_base = ref_allele.seq.first().copied().unwrap_or(b'N');
        let depth: u32 = record.alleles.iter().map(|a| a.support.num_obs).sum();
        coverage.observe(record.chrom_id, record.pos, ref_base, depth);
        let gcid = psp_to_global[record.chrom_id as usize];
        windows.observe(gcid, record.pos, ref_base, depth);
        if let Some(site) = SiteCounts::from_record(&record) {
            het.observe_site(site);
        }
    }

    let histogram = coverage.finish();
    let het_counts = het.finish();
    let window_depth = windows.finish();

    let model = SingleCopyCoverageModel::fit(&histogram, &CoverageFitConfig::default())
        .map_err(|e| e.to_string())?;
    let obs = obs_het(&het_counts, histogram.callable_positions);
    Ok(SampleState {
        model,
        inbreeding: 0.0, // set once the cohort Hexp is known (see caller)
        obs_het: obs,
        callable_positions: histogram.callable_positions,
        window_depth,
    })
}

/// Retains per-variant-window `(gc, mean depth)` while streaming a `.psp`
/// body, mirroring [`CoverageByGcAccumulator`]'s tiling but keeping the raw
/// per-window value for the windows that carry a variant.
struct WindowDepthCollector<'a> {
    index: &'a HashMap<WindowKey, usize>,
    out: Vec<Option<(f32, f32)>>,
    open: Option<(WindowKey, u64, u64, u64)>, // key, covered, gc, depth_sum
}

impl<'a> WindowDepthCollector<'a> {
    fn new(index: &'a HashMap<WindowKey, usize>, n_windows: usize) -> Self {
        Self {
            index,
            out: vec![None; n_windows],
            open: None,
        }
    }

    fn observe(&mut self, gcid: u32, pos: u32, ref_base: u8, depth: u32) {
        let key = (gcid, pos.saturating_sub(1) / WINDOW_BP);
        match self.open {
            Some((k, ..)) if k == key => {}
            _ => {
                self.finalize();
                self.open = Some((key, 0, 0, 0));
            }
        }
        if ref_base.eq_ignore_ascii_case(&b'N') {
            return;
        }
        if let Some((_, covered, gc, depth_sum)) = self.open.as_mut() {
            *covered += 1;
            if matches!(ref_base.to_ascii_uppercase(), b'G' | b'C') {
                *gc += 1;
            }
            *depth_sum = depth_sum.saturating_add(u64::from(depth));
        }
    }

    fn finalize(&mut self) {
        if let Some((key, covered, gc, depth_sum)) = self.open.take()
            && covered > 0
            && let Some(&idx) = self.index.get(&key)
        {
            let gc_frac = gc as f32 / covered as f32;
            let mean_depth = depth_sum as f32 / covered as f32;
            self.out[idx] = Some((gc_frac, mean_depth));
        }
    }

    fn finish(mut self) -> Vec<Option<(f32, f32)>> {
        self.finalize();
        self.out
    }
}

/// Build one sample's observation at a locus: its window's relative copy
/// number (via the fitted model) + the VCF `AD` counts. `None` if the sample
/// was not fit, has no covered window, or has no reads.
fn build_observation(
    record: &noodles_vcf::variant::RecordBuf,
    col: usize,
    states: &[Option<SampleState>],
    win: usize,
) -> Option<SampleObservation> {
    let state = states.get(col)?.as_ref()?;
    // `win < window_depth.len()`: the window index is frozen after the VCF
    // pre-pass and every collector is sized to that count — guard locally so a
    // future lazy-indexing change can't turn this into a data-dependent panic.
    let (gc, depth) = (*state.window_depth.get(win)?)?;
    let (ref_reads, alt_reads) = sample_ad(record, col)?;
    let total = ref_reads + alt_reads;
    if total == 0 {
        return None;
    }
    let rel = state
        .model
        .relative_copy_number(f64::from(gc), f64::from(depth));
    Some(SampleObservation {
        relative_copy_number: rel,
        alt_reads,
        total_reads: total,
        inbreeding_coefficient: state.inbreeding,
    })
}

/// `(chrom, 1-based pos)` if the record is a biallelic SNP, else `None`.
fn biallelic_snp_site(record: &noodles_vcf::variant::RecordBuf) -> Option<(&str, u32)> {
    if record.reference_bases().len() != 1 {
        return None;
    }
    let alts = record.alternate_bases().as_ref();
    if alts.len() != 1 || alts[0].len() != 1 {
        return None;
    }
    let pos = u32::try_from(usize::from(record.variant_start()?)).ok()?;
    Some((record.reference_sequence_name(), pos))
}

/// Assign / look up a canonical chromosome id for `name`.
fn intern_chrom(chrom_ids: &mut HashMap<String, u32>, name: &str) -> u32 {
    if let Some(&id) = chrom_ids.get(name) {
        return id;
    }
    let id = chrom_ids.len() as u32;
    chrom_ids.insert(name.to_string(), id);
    id
}

/// Cohort `(alt allele count, total called allele count)` from the GTs.
fn allele_counts(record: &noodles_vcf::variant::RecordBuf, n_samples: usize) -> (u64, u64) {
    let mut alt = 0u64;
    let mut called = 0u64;
    let samples = record.samples();
    for col in 0..n_samples {
        let Some(sample) = samples.get_index(col) else {
            continue;
        };
        if let Some(Some(Value::Genotype(gt))) = sample.get("GT") {
            for allele in gt.as_ref() {
                if let Some(idx) = allele.position() {
                    called += 1;
                    if idx > 0 {
                        alt += 1;
                    }
                }
            }
        }
    }
    (alt, called)
}

/// The `AD` `(ref, alt)` for sample column `col`, if present and biallelic.
fn sample_ad(record: &noodles_vcf::variant::RecordBuf, col: usize) -> Option<(u32, u32)> {
    let sample = record.samples().get_index(col)?;
    match sample.get("AD") {
        Some(Some(Value::Array(Array::Integer(values)))) => {
            let refc = (*values.first()?)?;
            let altc = (*values.get(1)?)?;
            Some((refc.max(0) as u32, altc.max(0) as u32))
        }
        _ => None,
    }
}

/// Map each `.psp` file in `dir` by its header sample id (the SRS accession
/// the VCF columns use), opening only the header of each.
fn map_psp_by_sample(dir: &Path) -> HashMap<String, PathBuf> {
    let mut map = HashMap::new();
    let Ok(entries) = std::fs::read_dir(dir) else {
        return map;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) != Some("psp") {
            continue;
        }
        let Ok(file) = File::open(&path) else {
            continue;
        };
        let Ok(reader) = PspReader::new(BufReader::with_capacity(64 * 1024, file)) else {
            continue;
        };
        let sample = reader.header().sample.clone();
        if let Some(prev) = map.insert(sample.clone(), path.clone()) {
            eprintln!(
                "  warning: duplicate .psp sample id {sample}: {} shadows {}",
                path.display(),
                prev.display()
            );
        }
    }
    map
}

fn write_samples_tsv(path: &Path, names: &[String], states: &[Option<SampleState>]) {
    let mut w = BufWriter::new(File::create(path).expect("create samples tsv"));
    writeln!(w, "sample\tsigma0\tobs_het\tF\tsingle_copy_scale").unwrap();
    for (name, state) in names.iter().zip(states.iter()) {
        if let Some(s) = state {
            writeln!(
                w,
                "{name}\t{:.5}\t{:.6}\t{:.4}\t{:.4}",
                s.model.single_copy_depth_sd(),
                s.obs_het,
                s.inbreeding,
                s.model.single_copy_scale(),
            )
            .unwrap();
        }
    }
    eprintln!("wrote {}", path.display());
}

fn write_lr_tsv(
    path: &Path,
    per_locus: &[(String, u32, f64)],
    prior: &ParalogPrior,
    curve: &ParalogFdrCurve,
) {
    let mut w = BufWriter::new(File::create(path).expect("create lr tsv"));
    writeln!(w, "chrom\tpos\tlr\tpost\tqval").unwrap();
    for (chrom, pos, lr) in per_locus {
        writeln!(
            w,
            "{chrom}\t{pos}\t{lr:.6}\t{:.6}\t{:.6}",
            prior.paralog_posterior(*lr),
            curve.q_of_lr(*lr),
        )
        .unwrap();
    }
    eprintln!("wrote {}", path.display());
}

fn report_flagged_profile(per_locus: &[(String, u32, f64)], curve: &ParalogFdrCurve) {
    let mut flagged = Vec::new();
    let mut kept = Vec::new();
    for (_, _, lr) in per_locus {
        if curve.q_of_lr(*lr) <= REPORT_FDR {
            flagged.push(*lr);
        } else {
            kept.push(*lr);
        }
    }
    eprintln!(
        "FDR {:.0}% flagged: {} loci (median LR {:.2}); kept: {} loci (median LR {:.2})",
        100.0 * REPORT_FDR,
        flagged.len(),
        median(&mut flagged),
        kept.len(),
        median(&mut kept),
    );
}

fn median(xs: &mut [f64]) -> f64 {
    if xs.is_empty() {
        return f64::NAN;
    }
    xs.sort_by(f64::total_cmp);
    xs[xs.len() / 2]
}
