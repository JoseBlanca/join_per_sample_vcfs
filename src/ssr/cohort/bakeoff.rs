//! The Stage-2 read-likelihood bake-off harness (plan
//! `ssr_stutter_scoring_model_bakeoff.md` §4): simulate a cohort with **known**
//! genotypes, score + genotype its reads under a chosen [`ReadLikelihoodModel`], and
//! compare the calls to truth.
//!
//! This is the simulate→score→genotype→compare driver. It is test-only (it pulls in
//! the [`sim`](super::sim) generator) and emits a reproducible model × metric table so
//! the production model choice rests on measured accuracy, not priors. Models A/B/C
//! plug in through the trait as they land; the generative side is varied separately
//! (the G1/G2/G3 fairness axis — §5) by feeding different [`SimCohortSpec`]s.
//!
//! **Fairness (Step 2 cut):** the genotyper runs on a [`ParamSet`] built **directly
//! from the simulation's truth chemistry** ([`param_set_from_truth`]) rather than from
//! the pre-pass estimator, so a model is judged on the *shape* of its `Qᵣ`, not on how
//! well the estimator recovered its parameters (plan §9 — "compare models, not
//! estimators"). The pre-pass (and Model B's out-of-frame leak fix) is scored
//! separately once the models exist.

use std::collections::{BTreeMap, HashMap};
use std::time::{Duration, Instant};

use crate::ssr::cohort::candidate_set::{CandidateCfg, assemble_candidates};
use crate::ssr::cohort::em::{EmCfg, run_locus_em_with};
use crate::ssr::cohort::em_init::seed_locus;
use crate::ssr::cohort::param_estimation::{
    DEFAULT_G0_FALLBACK_P, G0PseudocountDecay, ParamSet, PurityLevel, StutterShape,
};
use crate::ssr::cohort::read_model::ReadLikelihoodModel;
use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
use crate::ssr::cohort::sim::{GenerativeNoise, SimCohortSpec, TruthTable, simulate_with};

/// Diploid is the only supported ploidy (matches the simulator and the EM contract).
const PLOIDY: u8 = 2;

/// The outcome of scoring one generative dataset with one model: how the called
/// genotypes compared to the known truth, plus the wall time spent in the per-locus
/// genotyping EM (the `Qᵣ` hot path dominates it).
///
/// Counts (not ratios) are stored so several runs can be merged and so the
/// determinism check can compare them exactly; the ratios are derived on demand. The
/// runtime is advisory — wall time is not deterministic, so it is never part of an
/// equality assertion.
#[derive(Debug, Clone, Default)]
pub(crate) struct BakeoffMetrics {
    /// Present (sample, locus) genotypes compared to truth.
    pub(crate) genotypes_compared: usize,
    /// Of those, the count whose called allele-length multiset equals truth exactly.
    pub(crate) genotypes_concordant: usize,
    /// `|called − true|` in repeat units, one entry per aligned allele copy (a het
    /// contributes two). A no-call contributes nothing here (it is counted as a
    /// non-concordant comparison instead).
    pub(crate) allele_length_errors: Vec<u32>,
    /// Per emitted (non-no-call) genotype: its reported GQ and whether it matched truth.
    /// The reliability of GQ against this empirical correctness is the calibration metric
    /// (a model can be accurate yet over-confident — plan §6).
    pub(crate) gq_correct: Vec<(u8, bool)>,
    /// Wall time spent in the genotyping EM across all loci (advisory; non-deterministic).
    pub(crate) scoring_wall: Duration,
}

impl BakeoffMetrics {
    /// Fraction of compared genotypes that matched truth exactly (the primary metric).
    /// `1.0` when nothing was compared (vacuously concordant).
    pub(crate) fn concordance(&self) -> f64 {
        if self.genotypes_compared == 0 {
            return 1.0;
        }
        self.genotypes_concordant as f64 / self.genotypes_compared as f64
    }

    /// Mean per-allele-copy length error in repeat units (`0.0` if no copies aligned).
    pub(crate) fn mean_allele_length_error(&self) -> f64 {
        if self.allele_length_errors.is_empty() {
            return 0.0;
        }
        let total: u64 = self.allele_length_errors.iter().map(|&e| e as u64).sum();
        total as f64 / self.allele_length_errors.len() as f64
    }

    /// The largest single-allele length error seen (`0` if none).
    pub(crate) fn max_allele_length_error(&self) -> u32 {
        self.allele_length_errors.iter().copied().max().unwrap_or(0)
    }

    /// The allele-length-error distribution as `error → count` (the near-miss profile).
    pub(crate) fn allele_error_histogram(&self) -> BTreeMap<u32, usize> {
        let mut hist = BTreeMap::new();
        for &error in &self.allele_length_errors {
            *hist.entry(error).or_insert(0) += 1;
        }
        hist
    }

    /// Expected Calibration Error: the GQ-weighted gap between a model's stated
    /// confidence and its empirical correctness (`0.0` = perfectly calibrated, higher =
    /// more over/under-confident). Genotypes are bucketed by GQ in 10-Phred bins; each
    /// populated bin contributes `(n_bin/N)·|accuracy − confidence|`, where confidence is
    /// the bin's mean `1 − 10^(−GQ/10)`. `0.0` when nothing was emitted.
    pub(crate) fn expected_calibration_error(&self) -> f64 {
        if self.gq_correct.is_empty() {
            return 0.0;
        }
        // Bucket index = GQ/10 (so 0..=9, with GQ 90..=99 sharing the top bucket).
        let mut bucket_conf = [0.0f64; 10];
        let mut bucket_correct = [0u32; 10];
        let mut bucket_count = [0u32; 10];
        for &(gq, correct) in &self.gq_correct {
            let bucket = (gq as usize / 10).min(9);
            bucket_conf[bucket] += 1.0 - 10f64.powf(-(gq as f64) / 10.0);
            bucket_correct[bucket] += u32::from(correct);
            bucket_count[bucket] += 1;
        }
        let total = self.gq_correct.len() as f64;
        let mut ece = 0.0;
        for bucket in 0..10 {
            let n = bucket_count[bucket];
            if n == 0 {
                continue;
            }
            let accuracy = bucket_correct[bucket] as f64 / n as f64;
            let confidence = bucket_conf[bucket] / n as f64;
            ece += (n as f64 / total) * (accuracy - confidence).abs();
        }
        ece
    }
}

/// Build the frozen [`ParamSet`] the genotyping EM consumes **directly from the
/// simulation truth** — `ε`, stutter level, and sample groups are taken verbatim from
/// the chemistry the reads were drawn under (the Step-2 fairness choice).
///
/// The per-period parent shape `θ_period` is a cohort-shared *seed* the EM refines per
/// locus, but truth carries a shape per sample *group*, which may diverge. We seed each
/// period from the **first group's** shape and let the per-locus `θ_locus` refit adapt;
/// this is the one place the harness is not a pure mirror of truth, so it is called out
/// (a single-group cohort — the common G1 case — has no divergence at all). The `G₀`
/// decay is the coded fallback (the bake-off measures genotyping, not the allele-prior
/// fit). The per-`(group, period)` cell shapes are left empty: the per-locus EM scores
/// every sample at a locus under one `θ_locus`, so it never reads them.
pub(crate) fn param_set_from_truth(spec: &SimCohortSpec, truth: &TruthTable) -> ParamSet {
    let periods: Vec<u8> = {
        let mut ps: Vec<u8> = spec.loci.iter().map(|l| l.motif.period() as u8).collect();
        ps.sort_unstable();
        ps.dedup();
        ps
    };
    let seed_shape: StutterShape = truth
        .chemistry
        .first()
        .map(|c| c.shape)
        .expect("a simulated cohort has at least one sample group");

    let stutter_shape_parent: HashMap<u8, StutterShape> =
        periods.iter().map(|&p| (p, seed_shape)).collect();
    let pseudocount_decay_per_loci_group: HashMap<u8, G0PseudocountDecay> = periods
        .iter()
        .map(|&p| {
            (
                p,
                G0PseudocountDecay {
                    p: DEFAULT_G0_FALLBACK_P,
                },
            )
        })
        .collect();

    ParamSet {
        error_per_sample_group: truth.chemistry.iter().map(|c| c.error).collect(),
        stutter_shape_parent,
        stutter_shape_by_cell: HashMap::new(),
        level_seed: truth.chemistry.iter().map(|c| c.level).collect(),
        pseudocount_decay_per_loci_group,
        group_of_sample: truth.group_of_sample.clone(),
        f0_seed: 0.0,
        purity_level: PurityLevel::none(),
    }
}

/// Simulate `spec` under the clean G1 generator and score it under `model`. Thin
/// wrapper over [`run_bakeoff_with`].
pub(crate) fn run_bakeoff<M: ReadLikelihoodModel>(
    spec: &SimCohortSpec,
    model: &M,
) -> BakeoffMetrics {
    run_bakeoff_with(spec, &GenerativeNoise::none(), model)
}

/// Simulate the cohort described by `spec` under `noise` (the generative axis G1/G2/G3),
/// genotype every locus under `model` on truth-derived parameters, and score the calls
/// against the known genotypes.
///
/// The truth-derived parameters describe the *clean* chemistry only (whole-unit shape,
/// level, ε) — they do **not** encode the noise. That is deliberate: scoring messy reads
/// (G2/G3) on clean parameters is exactly the model-mismatch robustness the bake-off
/// weighs (plan §5/§9). The loci are matched back to the spec by start coordinate
/// (synthetic specs use distinct starts), so a partially-absent locus is handled
/// correctly rather than relying on a positional 1:1 with the merged stream.
pub(crate) fn run_bakeoff_with<M: ReadLikelihoodModel>(
    spec: &SimCohortSpec,
    noise: &GenerativeNoise,
    model: &M,
) -> BakeoffMetrics {
    let cohort = simulate_with(spec, noise);
    let params = param_set_from_truth(spec, &cohort.truth);

    // Spec-locus index keyed by start coordinate, so a merged locus maps to its truth
    // column without assuming the stream is positionally 1:1 with the spec.
    let locus_idx_by_start: HashMap<u32, usize> = spec
        .loci
        .iter()
        .enumerate()
        .map(|(idx, l)| (l.start, idx))
        .collect();
    debug_assert_eq!(
        locus_idx_by_start.len(),
        spec.loci.len(),
        "bake-off specs must use distinct locus start coordinates"
    );

    let rung_cfg = RungCfg::dev_default();
    // Force locus admission for the bake-off. The G2/G3 generative cases deliberately
    // dial the out-of-frame fraction *above* the production admission threshold to stress
    // the read-likelihood models, so the real gate would filter those loci to a no-call
    // *before any model runs* — identically for every model, leaving nothing to compare.
    // The bake-off isolates the Qᵣ model (downstream of admission), so it admits every
    // locus (`max_out_of_frame_frac = 1.0`) and lets the model genotype the messy reads.
    let cand_cfg = CandidateCfg {
        max_out_of_frame_frac: 1.0,
        ..CandidateCfg::dev_default()
    };
    let em_cfg = EmCfg::dev_default();

    let mut metrics = BakeoffMetrics::default();
    for item in cohort.merger() {
        let (_, locus) = item.expect("simulated cohort merges cleanly");
        let locus_idx = locus_idx_by_start[&locus.locus.start];

        let rungs = build_rungs(&locus, &rung_cfg);
        let candidates = assemble_candidates(&locus, &rungs, PLOIDY, &cand_cfg);
        let seed = seed_locus(
            &locus,
            &rungs,
            &candidates,
            &params,
            PLOIDY,
            rung_cfg.prominence,
        );
        let f_present = vec![0.0; locus.samples.len()];

        let started = Instant::now();
        let call = run_locus_em_with(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            PLOIDY,
            &em_cfg,
            &f_present,
            &params.level_seed,
            model,
        )
        .0;
        metrics.scoring_wall += started.elapsed();

        for (k, &global) in locus.present.iter().enumerate() {
            let Some(truth_genotype) = cohort.truth.genotype(global as usize, locus_idx) else {
                continue; // present in the stream but absent in truth — cannot happen, skip safely
            };
            score_genotype(
                &call.calls[k].genotype_units,
                &truth_genotype.allele_units,
                call.calls[k].gq,
                &mut metrics,
            );
        }
    }
    metrics
}

/// Compare one called allele-length multiset to truth, folding the result into
/// `metrics`: an exact multiset match is concordant, and each aligned copy's
/// `|called − true|` enters the error distribution. A no-call (empty `called`) counts
/// as a non-concordant comparison with no error entries.
fn score_genotype(called: &[u16], truth: &[u16], gq: u8, metrics: &mut BakeoffMetrics) {
    metrics.genotypes_compared += 1;
    if called.is_empty() {
        return; // no-call: a miss, but no |Δ| or GQ to attribute
    }
    let mut called_sorted = called.to_vec();
    let mut truth_sorted = truth.to_vec();
    called_sorted.sort_unstable();
    truth_sorted.sort_unstable();

    let correct = called_sorted == truth_sorted;
    if correct {
        metrics.genotypes_concordant += 1;
    }
    metrics.gq_correct.push((gq, correct));
    // Per-copy absolute error on the sorted alignment (well-defined when ploidy matches,
    // which it does for the diploid contract).
    if called_sorted.len() == truth_sorted.len() {
        for (&c, &t) in called_sorted.iter().zip(&truth_sorted) {
            metrics
                .allele_length_errors
                .push((c as i32 - t as i32).unsigned_abs());
        }
    }
}

/// Render a labelled set of bake-off results as a fixed-width table — the reproducible
/// artefact the model decision is read off (one row per model, or per
/// scoring × generative cell once the 3×3 sweep lands).
pub(crate) fn format_results_table(rows: &[(&str, &BakeoffMetrics)]) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    writeln!(
        out,
        "{:<24} {:>11} {:>10} {:>8} {:>8} {:>7} {:>10}",
        "model × generative", "concordance", "mean |Δ|", "max |Δ|", "ECE", "n", "scoring ms"
    )
    .unwrap();
    for (label, m) in rows {
        writeln!(
            out,
            "{:<24} {:>11.4} {:>10.4} {:>8} {:>8.4} {:>7} {:>10.1}",
            label,
            m.concordance(),
            m.mean_allele_length_error(),
            m.max_allele_length_error(),
            m.expected_calibration_error(),
            m.genotypes_compared,
            m.scoring_wall.as_secs_f64() * 1e3,
        )
        .unwrap();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::param_estimation::{PerBaseError, SampleGroupId, StutterLevel};
    use crate::ssr::cohort::read_model::{ClassicStutterModel, HipstrModel};
    use crate::ssr::cohort::sim::{SimChemistry, SimCohortSpec, SimGenotype, SimLocus, SimSample};
    use crate::ssr::types::Motif;

    fn motif(bytes: &[u8]) -> Motif {
        Motif::new(bytes).unwrap()
    }

    /// A clean single-group G1 cohort: low stutter, low ε, two CA loci, a mix of
    /// homozygotes and well-separated hets at good depth — the case every model should
    /// nail, so it pins the harness wiring (a break shows up as a concordance collapse).
    fn clean_g1_spec(seed: u64) -> SimCohortSpec {
        let chem = SimChemistry {
            error: PerBaseError(0.001),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.04,
                slope: 0.0,
            },
        };
        // Homozygotes and separated hets (≥2 units apart so the het is resolvable).
        let genos = [
            SimGenotype::homozygous(8, 2),
            SimGenotype::homozygous(6, 2),
            SimGenotype::diploid(6, 10),
            SimGenotype::homozygous(10, 2),
            SimGenotype::diploid(7, 12),
            SimGenotype::homozygous(9, 2),
        ];
        SimCohortSpec {
            seed,
            loci: vec![
                SimLocus {
                    chrom: "chr1".into(),
                    start: 40,
                    motif: motif(b"CA"),
                    ref_units: 8,
                },
                SimLocus {
                    chrom: "chr1".into(),
                    start: 400,
                    motif: motif(b"CA"),
                    ref_units: 9,
                },
            ],
            groups: vec![chem],
            samples: genos
                .iter()
                .enumerate()
                .map(|(i, g)| SimSample {
                    name: format!("S{i}"),
                    group: SampleGroupId(0),
                    genotypes: vec![Some(g.clone()), Some(g.clone())],
                })
                .collect(),
            depth: 100,
        }
    }

    #[test]
    fn param_set_from_truth_mirrors_the_simulation_chemistry() {
        let spec = clean_g1_spec(1);
        let cohort = simulate_with(&spec, &GenerativeNoise::none());
        let params = param_set_from_truth(&spec, &cohort.truth);

        assert_eq!(params.error_per_sample_group, vec![PerBaseError(0.001)]);
        assert_eq!(params.level_seed[0].baseline, 0.04);
        assert_eq!(params.group_of_sample, vec![SampleGroupId(0); 6]);
        // The CA period (2) seeds its parent shape from the (only) group's shape.
        assert_eq!(params.stutter_shape_parent[&2].decay, 0.1);
        assert!(params.pseudocount_decay_per_loci_group.contains_key(&2));
    }

    #[test]
    fn model_b_is_highly_concordant_on_a_clean_g1_cohort() {
        let metrics = run_bakeoff(&clean_g1_spec(2026), &ClassicStutterModel);
        // 6 samples × 2 loci = 12 genotypes compared.
        assert_eq!(metrics.genotypes_compared, 12);
        assert!(
            metrics.concordance() >= 0.9,
            "Model B should nail a clean low-stutter cohort, got {:.3}\n{}",
            metrics.concordance(),
            format_results_table(&[("classic (B)", &metrics)]),
        );
        // Near-misses, when they happen, are small — the delimiter's old hard collapse
        // would show large |Δ|.
        assert!(
            metrics.max_allele_length_error() <= 2,
            "allele-length errors should be near-misses, got max {}",
            metrics.max_allele_length_error(),
        );
    }

    #[test]
    fn model_a_runs_end_to_end_and_is_reasonable_on_g1() {
        // Model A (HipSTR in/out-of-frame) plugs into the same harness; G1 has no
        // out-of-frame draws, so this exercises its in-frame branch and pins the wiring.
        let metrics = run_bakeoff(&clean_g1_spec(2026), &HipstrModel);
        assert_eq!(metrics.genotypes_compared, 12);
        assert!(
            metrics.concordance() >= 0.8,
            "Model A should be reasonable on a clean cohort, got {:.3}",
            metrics.concordance(),
        );
    }

    /// A richer cohort for the sweep: 16 samples (homozygotes + separated hets across a
    /// range of unit lengths, incl. a long allele) over three CA loci, so concordance is
    /// stable enough to compare models. One sample group (clean chemistry).
    fn sweep_spec(seed: u64) -> SimCohortSpec {
        let chem = SimChemistry {
            error: PerBaseError(0.002),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.15,
            },
            level: StutterLevel {
                baseline: 0.06,
                slope: 0.004, // longer alleles stutter a touch more
            },
        };
        let genos = [
            SimGenotype::homozygous(8, 2),
            SimGenotype::homozygous(12, 2),
            SimGenotype::diploid(6, 10),
            SimGenotype::homozygous(10, 2),
            SimGenotype::diploid(7, 14),
            SimGenotype::homozygous(16, 2),
            SimGenotype::diploid(9, 13),
            SimGenotype::homozygous(11, 2),
        ];
        let samples = (0..16)
            .map(|i| SimSample {
                name: format!("S{i}"),
                group: SampleGroupId(0),
                genotypes: vec![
                    Some(genos[i % genos.len()].clone()),
                    Some(genos[(i + 3) % genos.len()].clone()),
                    Some(genos[(i + 5) % genos.len()].clone()),
                ],
            })
            .collect();
        SimCohortSpec {
            seed,
            loci: vec![
                SimLocus {
                    chrom: "chr1".into(),
                    start: 60,
                    motif: motif(b"CA"),
                    ref_units: 10,
                },
                SimLocus {
                    chrom: "chr1".into(),
                    start: 600,
                    motif: motif(b"CA"),
                    ref_units: 12,
                },
                SimLocus {
                    chrom: "chr1".into(),
                    start: 1200,
                    motif: motif(b"CA"),
                    ref_units: 8,
                },
            ],
            groups: vec![chem],
            samples,
            depth: 60,
        }
    }

    /// THE BAKE-OFF (post-decision form): the surviving models — A (production, HipSTR)
    /// and B (the retained reference) — across the G1/G2/G3 generative axis plus a depth
    /// sweep, emitting the reproducible table (run with `--nocapture`). Model C was
    /// eliminated and removed (it collapsed on out-of-frame reads); the full three-model
    /// numbers that drove the decision live in the implementation report. This now guards
    /// that both survivors stay concordant across the generative cases.
    #[test]
    fn bakeoff_sweep_and_depth_robustness_emit_the_table() {
        let spec = sweep_spec(2026);
        let generatives = [
            ("G1-clean", GenerativeNoise::none()),
            ("G2-oof", GenerativeNoise::hipstr_like()),
            ("G3-messy", GenerativeNoise::messy()),
        ];

        // Run every (model, generative) cell; store owned metrics so the table can borrow.
        let mut cells: Vec<(String, BakeoffMetrics)> = Vec::new();
        for (gname, noise) in &generatives {
            cells.push((
                format!("A-hipstr  / {gname}"),
                run_bakeoff_with(&spec, noise, &HipstrModel),
            ));
            cells.push((
                format!("B-classic / {gname}"),
                run_bakeoff_with(&spec, noise, &ClassicStutterModel),
            ));
        }
        let rows: Vec<(&str, &BakeoffMetrics)> =
            cells.iter().map(|(l, m)| (l.as_str(), m)).collect();
        println!(
            "\n=== Stage-2 read-likelihood bake-off (survivors A vs B) ===\n{}",
            format_results_table(&rows)
        );

        // Depth-robustness sweep on G3-messy (the case the decision weighed most).
        println!("--- depth robustness on G3-messy (concordance) ---");
        for depth in [10u32, 20, 40, 80] {
            let mut s = spec.clone();
            s.depth = depth;
            let a = run_bakeoff_with(&s, &GenerativeNoise::messy(), &HipstrModel);
            let b = run_bakeoff_with(&s, &GenerativeNoise::messy(), &ClassicStutterModel);
            println!(
                "depth {depth:>3}:  A {:.3}   B {:.3}",
                a.concordance(),
                b.concordance()
            );
        }

        // Both survivors stay highly concordant across ALL three generative cases
        // (incl. G3-messy) — the property the production choice rests on.
        for (label, m) in &cells {
            assert_eq!(
                m.genotypes_compared, 48,
                "{label} should compare 48 genotypes"
            );
            assert!(
                m.concordance() >= 0.9,
                "{label} should stay concordant, got {:.3}",
                m.concordance()
            );
        }
    }

    #[test]
    fn production_admission_gate_admits_realistic_but_rejects_excessive_out_of_frame() {
        use crate::ssr::cohort::candidate_set::{Admission, assemble_candidates};
        use crate::ssr::cohort::rung_ladder::build_rungs;

        // Confirm the support-aware periodicity gate on *real simulated* cohorts, through
        // the production `dev_default` config (threshold 10% out-of-frame) — NOT the
        // bake-off's admit-everything override.
        let spec = sweep_spec(2026);
        let rc = RungCfg::dev_default();
        let cc = CandidateCfg::dev_default();
        let admits = |noise: &GenerativeNoise| -> Vec<Admission> {
            simulate_with(&spec, noise)
                .merger()
                .map(|item| {
                    let (_, locus) = item.unwrap();
                    let rungs = build_rungs(&locus, &rc);
                    assemble_candidates(&locus, &rungs, PLOIDY, &cc).admit
                })
                .collect()
        };

        // Realistic 5% out-of-frame: every locus is still admitted (the bug fix — the old
        // presence-based gate filtered these on the first stray read).
        let realistic = GenerativeNoise {
            out_of_frame_rate: 0.05,
            impurity_rate: 0.0,
            extra_substitution: 0.0,
        };
        assert!(
            admits(&realistic).iter().all(|a| *a == Admission::Pass),
            "5% out-of-frame loci must admit, got {:?}",
            admits(&realistic)
        );

        // Excessive 30% out-of-frame: the gate still does its job and rejects them.
        let excessive = GenerativeNoise {
            out_of_frame_rate: 0.30,
            impurity_rate: 0.0,
            extra_substitution: 0.0,
        };
        assert!(
            admits(&excessive)
                .iter()
                .all(|a| *a == Admission::NotPeriodic),
            "30% out-of-frame loci must be rejected, got {:?}",
            admits(&excessive)
        );
    }

    #[test]
    fn results_table_renders_multiple_models() {
        // The reproducible artefact: the production model A against the retained B baseline.
        let a = run_bakeoff(&clean_g1_spec(2026), &HipstrModel);
        let b = run_bakeoff(&clean_g1_spec(2026), &ClassicStutterModel);
        let table = format_results_table(&[("hipstr (A)", &a), ("classic (B)", &b)]);
        assert!(table.contains("concordance"));
        assert!(table.contains("hipstr (A)"));
        assert!(table.contains("classic (B)"));
        assert_eq!(table.lines().count(), 3); // header + two rows
    }

    #[test]
    fn bakeoff_accuracy_metrics_are_deterministic() {
        // The accuracy metrics are a pure function of the spec + model (only wall time
        // varies), so two runs must agree exactly — the byte-identity contract surfaced
        // at the harness level.
        let a = run_bakeoff(&clean_g1_spec(7), &ClassicStutterModel);
        let b = run_bakeoff(&clean_g1_spec(7), &ClassicStutterModel);
        assert_eq!(a.genotypes_compared, b.genotypes_compared);
        assert_eq!(a.genotypes_concordant, b.genotypes_concordant);
        assert_eq!(a.allele_length_errors, b.allele_length_errors);
        assert_eq!(a.gq_correct, b.gq_correct); // calibration data deterministic too
    }

    #[test]
    fn score_genotype_counts_match_miss_no_call_and_calibration() {
        let mut m = BakeoffMetrics::default();
        score_genotype(&[6, 10], &[6, 10], 60, &mut m); // exact (order-insensitive)
        score_genotype(&[10, 6], &[6, 10], 40, &mut m); // exact, reversed
        score_genotype(&[6, 9], &[6, 10], 20, &mut m); // one copy off by 1 (a miss)
        score_genotype(&[], &[6, 10], 0, &mut m); // no-call

        assert_eq!(m.genotypes_compared, 4);
        assert_eq!(m.genotypes_concordant, 2);
        // Three aligned genotypes (the no-call contributes none) → 6 copies; one error of 1.
        assert_eq!(m.allele_length_errors.len(), 6);
        assert_eq!(m.max_allele_length_error(), 1);
        assert_eq!(m.allele_error_histogram(), BTreeMap::from([(0, 5), (1, 1)]));
        // Calibration records the three emitted genotypes (no-call excluded), with the
        // GQ-20 one being the only incorrect call.
        assert_eq!(m.gq_correct, vec![(60, true), (40, true), (20, false)]);
        assert!(m.expected_calibration_error() >= 0.0);
    }

    #[test]
    fn calibration_error_is_zero_for_a_perfectly_calibrated_set() {
        // GQ 0 ⇒ confidence 0; an all-wrong GQ-0 bucket is perfectly calibrated (ECE 0).
        let mut m = BakeoffMetrics::default();
        for _ in 0..5 {
            score_genotype(&[6, 9], &[6, 10], 0, &mut m); // wrong, GQ 0
        }
        assert!(
            m.expected_calibration_error() < 1e-9,
            "all-wrong GQ-0 calls are perfectly calibrated, got {}",
            m.expected_calibration_error()
        );
    }
}
