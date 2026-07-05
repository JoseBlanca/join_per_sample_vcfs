//! Integration tests for the Stage 6 posterior engine.
//!
//! These build synthetic `MergedRecord` streams directly — Stage 5
//! round-trip is exercised by the per-group merger's own tests, and
//! Stage 6's contract is independent of Stage 5's internal layout
//! per the implementation plan.

use pop_var_caller::pileup_record::AlleleSupportStats;
use pop_var_caller::var_calling::per_group_merger::{
    MergedAllele, MergedRecord, PerGroupMergerError, genotype_order,
};
use pop_var_caller::var_calling::posterior_engine::backends::ExactMath;
use pop_var_caller::var_calling::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig, PosteriorEngineError, PosteriorRecord,
};

fn simple_alleles(seqs: &[&[u8]]) -> Vec<MergedAllele> {
    seqs.iter()
        .map(|s| MergedAllele {
            seq: s.to_vec(),
            is_compound: false,
            constituents: Vec::new(),
        })
        .collect()
}

fn merged_record_simple(
    chrom_id: u32,
    pos: u32,
    alleles: Vec<&[u8]>,
    ploidy: u8,
    likelihoods: Vec<Vec<f64>>,
) -> MergedRecord {
    let alleles = simple_alleles(&alleles);
    let n_samples = likelihoods.len();
    let n_alleles = alleles.len();
    let n_genotypes = genotype_order(ploidy, n_alleles).len();
    let ref_len = alleles[0].seq.len() as u32;
    let scalars = vec![AlleleSupportStats::default(); n_samples * n_alleles];
    let other_scalars = vec![AlleleSupportStats::default(); n_samples];
    let chain_anchor_flags = vec![false; n_samples * n_alleles];
    let log_likelihoods: Vec<f64> = likelihoods.into_iter().flatten().collect();
    MergedRecord {
        chrom_id,
        start: pos,
        end: pos + ref_len - 1,
        alleles,
        ploidy,
        n_samples,
        n_genotypes,
        scalars,
        other_scalars,
        chain_anchor_flags,
        log_likelihoods,
    }
}

fn collect_ok(records: Vec<MergedRecord>) -> Vec<PosteriorRecord> {
    // Pin to `ExactMath` so the strict-tolerance assertions below
    // (e.g. row-sum within `1e-9`) keep the bit-identical-baseline
    // semantics they were written against. Backend equivalence under
    // the looser approximate-parity budget is enforced separately by
    // the `tests_math_backend_accuracy` harness in the engine module.
    let upstream = records.into_iter().map(Ok::<_, PerGroupMergerError>);
    PosteriorEngine::with_math_backend(
        upstream,
        PosteriorEngineConfig::with_project_defaults(),
        ExactMath,
    )
    .collect::<Result<Vec<_>, PosteriorEngineError>>()
    .expect("posterior engine emitted error")
}

#[test]
fn streams_three_synthetic_records_in_genomic_order() {
    let records = vec![
        merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            2,
            vec![
                vec![0.0, -50.0, -50.0],
                vec![-50.0, -50.0, 0.0],
                vec![-50.0, 0.0, -50.0],
            ],
        ),
        merged_record_simple(
            1,
            200,
            vec![b"A", b"G"],
            2,
            vec![
                vec![0.0, -50.0, -50.0],
                vec![0.0, -50.0, -50.0],
                vec![0.0, -50.0, -50.0],
            ],
        ),
        merged_record_simple(
            1,
            300,
            vec![b"AC", b"A"],
            2,
            vec![
                vec![-50.0, 0.0, -50.0],
                vec![-50.0, 0.0, -50.0],
                vec![-50.0, -50.0, 0.0],
            ],
        ),
    ];

    let out = collect_ok(records);
    assert_eq!(out.len(), 3);
    let starts: Vec<u32> = out.iter().map(|r| r.locus.start).collect();
    assert_eq!(starts, vec![100, 200, 300]);

    assert_eq!(out[0].best_genotype, vec![0, 2, 1]);
    assert!(out[1].qual_phred < 1.0, "qual = {}", out[1].qual_phred);
    assert_eq!(out[2].best_genotype[2], 2);
}

#[test]
fn cohort_prior_pulls_lone_weak_alt_back_to_ref() {
    let mut likelihoods: Vec<Vec<f64>> = (0..4).map(|_| vec![0.0, -50.0, -50.0]).collect();
    likelihoods.push(vec![-1.0, 0.0, -1.0]);
    let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, likelihoods);
    let pr = collect_ok(vec![record]).into_iter().next().unwrap();
    assert_eq!(pr.best_genotype, vec![0, 0, 0, 0, 0]);
    assert!(pr.allele_frequencies[1] < 0.1);
}

/// The **empirical-Bayes large-cohort regime** (arch §10, Milestone 4): five
/// samples with weak-but-consistent het evidence let the cohort overrule the
/// rare-allele species prior and call het.
///
/// This is the unit-level proof of cohort strength-borrowing. Each sample alone
/// (the single-sample species prior) would call hom-ref against the rare-allele
/// prior; pooled via the leave-one-out frequency-posterior update
/// (`α'_s = α_species + expected cohort allele copies of the *other* samples`),
/// the accumulated evidence sharpens each sample's prior toward het. Landed with
/// the flat-first-step + leave-one-out EM (arch §10.2).
#[test]
fn cohort_evidence_overcomes_rare_allele_prior_when_all_samples_agree() {
    let likelihoods: Vec<Vec<f64>> = (0..5).map(|_| vec![-2.0, 0.0, -2.0]).collect();
    let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, likelihoods);
    let pr = collect_ok(vec![record]).into_iter().next().unwrap();
    for &g in &pr.best_genotype {
        assert_eq!(g, 1, "best_genotype = {:?}", pr.best_genotype);
    }
    assert!(
        pr.allele_frequencies[1] > 0.15,
        "p̂[alt] = {}",
        pr.allele_frequencies[1]
    );
}

/// Cohort sharpening generalises past biallelic diploid: a **triallelic**
/// diploid site where five samples all weakly favour the same REF/ALT1
/// heterozygote is called het by the cohort, while the identical evidence from a
/// lone sample stays hom-ref (the leave-one-out single-sample guarantee, arch
/// §10.2a). Genotype order for 3 alleles is `AA, AC, CC, AG, CG, GG`, so the
/// A/C het is index 1 and the all-ref hom is index 0.
#[test]
fn cohort_sharpens_a_multiallelic_het_a_lone_sample_leaves_hom_ref() {
    // Weakly favour the A/C het (index 1) over everything else.
    let weak_ac_het = || {
        let mut ll = vec![-2.0_f64; 6];
        ll[1] = 0.0;
        ll
    };

    let cohort: Vec<Vec<f64>> = (0..5).map(|_| weak_ac_het()).collect();
    let record = merged_record_simple(1, 100, vec![b"A", b"C", b"G"], 2, cohort);
    let pr = collect_ok(vec![record]).into_iter().next().unwrap();
    // ALT2 (G) is unused → pruned to the biallelic [A, C]; the A/C het remains
    // index 1 in the pruned genotype order (AA, AC, CC).
    for &g in &pr.best_genotype {
        assert_eq!(g, 1, "cohort best_genotype = {:?}", pr.best_genotype);
    }
    assert!(
        pr.allele_frequencies[1] > 0.15,
        "p̂[C] = {}",
        pr.allele_frequencies[1]
    );

    // The same weak evidence, alone, must not be talked into a het.
    let lone = merged_record_simple(1, 100, vec![b"A", b"C", b"G"], 2, vec![weak_ac_het()]);
    let pr_lone = collect_ok(vec![lone]).into_iter().next().unwrap();
    assert_eq!(
        pr_lone.best_genotype,
        vec![0],
        "lone best_genotype = {:?}",
        pr_lone.best_genotype
    );
}

/// Cohort sharpening at higher ploidy: a **tetraploid** biallelic site where
/// five samples all weakly favour the balanced (2-ALT-copy) genotype is called
/// heterozygous by the cohort, while a lone sample with the same evidence stays
/// hom-ref. Genotype order for ploidy 4 / 2 alleles is by ALT-copy count
/// `0,1,2,3,4`, so the 2-copy genotype is index 2 and the all-ref hom is 0.
#[test]
fn cohort_sharpens_a_polyploid_het_a_lone_sample_leaves_hom_ref() {
    let weak_two_copy = || {
        let mut ll = vec![-2.0_f64; 5];
        ll[2] = 0.0;
        ll
    };

    let cohort: Vec<Vec<f64>> = (0..5).map(|_| weak_two_copy()).collect();
    let record = merged_record_simple(1, 100, vec![b"A", b"C"], 4, cohort);
    let pr = collect_ok(vec![record]).into_iter().next().unwrap();
    for &g in &pr.best_genotype {
        assert_eq!(g, 2, "cohort best_genotype = {:?}", pr.best_genotype);
    }
    assert!(
        pr.allele_frequencies[1] > 0.15,
        "p̂[alt] = {}",
        pr.allele_frequencies[1]
    );

    let lone = merged_record_simple(1, 100, vec![b"A", b"C"], 4, vec![weak_two_copy()]);
    let pr_lone = collect_ok(vec![lone]).into_iter().next().unwrap();
    assert_eq!(
        pr_lone.best_genotype,
        vec![0],
        "lone best_genotype = {:?}",
        pr_lone.best_genotype
    );
}

/// The leave-one-out cohort prior couples samples through `expected_counts`, but
/// with unambiguous per-sample evidence (no genotype-probability ties) the
/// per-sample calls stay exactly equivariant under a sample-order permutation.
/// (The order-dependent floating sum in `expected_counts` can only flip a call
/// at an exact tie, which this fixture avoids.)
#[test]
fn cohort_loo_per_sample_calls_follow_sample_permutation() {
    // Distinct, sharp evidence: hom-ref / het / hom-alt / hom-ref / het.
    let rows: Vec<Vec<f64>> = vec![
        vec![0.0, -30.0, -60.0],
        vec![-30.0, 0.0, -30.0],
        vec![-60.0, -30.0, 0.0],
        vec![0.0, -30.0, -60.0],
        vec![-30.0, 0.0, -30.0],
    ];
    let base = collect_ok(vec![merged_record_simple(
        1,
        100,
        vec![b"A", b"C"],
        2,
        rows.clone(),
    )])
    .into_iter()
    .next()
    .unwrap();

    let mut reversed = rows.clone();
    reversed.reverse();
    let permuted = collect_ok(vec![merged_record_simple(
        1,
        100,
        vec![b"A", b"C"],
        2,
        reversed,
    )])
    .into_iter()
    .next()
    .unwrap();

    let n = base.best_genotype.len();
    assert_eq!(n, permuted.best_genotype.len());
    for i in 0..n {
        assert_eq!(
            base.best_genotype[i],
            permuted.best_genotype[n - 1 - i],
            "sample {i} call not equivariant under sample-order reversal: {:?} vs reversed {:?}",
            base.best_genotype,
            permuted.best_genotype
        );
    }
}

/// A cohort never converges on its flat first iteration: the emitted genotypes
/// come from the leave-one-out prior, so the EM runs iteration 1 (flat) plus at
/// least one steady-state (LOO) iteration, and `converged` reflects LOO
/// stability rather than the likelihood-only seed (arch §10.2b).
#[test]
fn cohort_runs_at_least_one_leave_one_out_iteration_before_converging() {
    let likelihoods: Vec<Vec<f64>> = (0..5).map(|_| vec![-2.0, 0.0, -2.0]).collect();
    let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, likelihoods);
    let pr = collect_ok(vec![record]).into_iter().next().unwrap();
    assert!(
        pr.diagnostics.iterations >= 2,
        "cohort converged on the flat step alone: iterations = {}",
        pr.diagnostics.iterations
    );
    assert!(pr.diagnostics.converged, "expected convergence");
}

#[test]
fn engine_output_is_bit_identical_across_runs() {
    let likelihoods: Vec<Vec<f64>> = (0..3)
        .map(|s| {
            let base = -(s as f64);
            vec![base, base - 1.0, base - 2.0]
        })
        .collect();
    let record_a = merged_record_simple(1, 100, vec![b"A", b"C"], 2, likelihoods.clone());
    let record_b = merged_record_simple(1, 100, vec![b"A", b"C"], 2, likelihoods);
    let out_a = collect_ok(vec![record_a]);
    let out_b = collect_ok(vec![record_b]);
    assert_eq!(out_a, out_b);
}

#[test]
fn config_override_propagates_through_engine() {
    let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![-2.0, 0.0, -2.0]]);
    // Crank the inbreeding to 1 ⇒ heterozygote prior mass is 0 ⇒
    // best genotype is forced to a homozygote even though the
    // likelihood prefers the heterozygote. `#[non_exhaustive]` on
    // `PosteriorEngineConfig` disallows struct-update syntax across
    // crates, so set the field after construction.
    #[allow(clippy::field_reassign_with_default)]
    let config = {
        let mut c = PosteriorEngineConfig::default();
        c.fixation_index_default = 1.0;
        c
    };
    let upstream = std::iter::once(Ok::<_, PerGroupMergerError>(record));
    let pr = PosteriorEngine::with_config(upstream, config)
        .next()
        .unwrap()
        .expect("posterior");
    assert_ne!(pr.best_genotype[0], 1);
}

#[test]
fn upstream_error_propagates_as_posterior_engine_upstream_error() {
    let rec = merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
    let err = PerGroupMergerError::RefFetch {
        chrom_id: 1,
        start: 200,
        end: 200,
        source: std::io::Error::other("synthetic"),
    };
    let upstream = vec![Ok(rec), Err(err)].into_iter();
    let out: Vec<_> = PosteriorEngine::new(upstream).collect();
    assert!(out[0].is_ok());
    assert!(matches!(out[1], Err(PosteriorEngineError::Upstream(_))));
}

#[test]
fn tetraploid_record_emits_posteriors_with_correct_n_genotypes() {
    // Tetraploid biallelic: genotype_order(4, 2).len() == 5.
    let n_genotypes = genotype_order(4, 2).len();
    let rec = merged_record_simple(1, 100, vec![b"A", b"C"], 4, vec![vec![0.0; n_genotypes]]);
    let pr = collect_ok(vec![rec]).into_iter().next().unwrap();
    assert_eq!(pr.n_genotypes, n_genotypes);
    let sum: f64 = pr.posteriors_row(0).iter().sum();
    assert!((sum - 1.0).abs() < 1e-9, "row sum = {sum}");
}
