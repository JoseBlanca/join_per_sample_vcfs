//! Integration tests for the contamination side-pass.
//!
//! These exercise the public crate-level API surface
//! (`var_calling::contamination_estimation::estimate_contamination`)
//! end-to-end with the [`PosteriorEngine`] consumer side. Both stages
//! are iterator-shaped, so the tests build synthetic upstream streams
//! directly rather than rounding through `.psp` on-disk artefacts —
//! the on-disk path is exercised by `psp_to_pileup_integration.rs`.

use merge_per_sample_vcfs::per_sample_pileup::pileup::{
    AlleleObservation, AlleleSupportStats, PileupRecord,
};
use merge_per_sample_vcfs::var_calling::contamination_estimation::{
    ContaminationEstimateSource, ContaminationEstimates, ContaminationEstimationConfig,
    ContaminationEstimationError, StoppingMode, estimate_contamination,
};
use merge_per_sample_vcfs::var_calling::per_group_merger::{
    MergedAllele, MergedRecord, PerGroupMergerError, genotype_order,
};
use merge_per_sample_vcfs::var_calling::per_position_merger::{
    PerPositionMergerError, PerPositionPileups,
};
use merge_per_sample_vcfs::var_calling::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig, PosteriorRecord,
};

// ----- Fixture builders ------------------------------------------------

fn obs(seq: &[u8], num_obs: u32, mean_bq_err: f64) -> AlleleObservation {
    let q_sum = if num_obs == 0 {
        0.0
    } else {
        mean_bq_err.ln() * f64::from(num_obs)
    };
    let support = AlleleSupportStats::new(num_obs, q_sum, num_obs / 2, num_obs / 2, 0);
    AlleleObservation::new(seq.to_vec(), support, Vec::new())
}

fn rec(pos: u32, alleles: Vec<AlleleObservation>) -> PileupRecord {
    PileupRecord::new(0, pos, alleles)
}

fn pos_pileups(pos: u32, per_sample: Vec<Option<PileupRecord>>) -> PerPositionPileups {
    PerPositionPileups {
        chrom_id: 0,
        pos,
        per_sample,
    }
}

/// Three-sample synthetic stream: sample 0 contaminated at `c_true`,
/// sample 1 hom-REF (clean), sample 2 hom-ALT (provides cohort
/// polymorphism so Step 1a accepts).
fn synth_three_sample_stream(
    n_sites: u32,
    n_reads: u32,
    c_true: f64,
) -> Vec<Result<PerPositionPileups, PerPositionMergerError>> {
    let n_contam = (f64::from(n_reads) * c_true).round() as u32;
    let n_own = n_reads - n_contam;
    let mean_err = 0.001;
    (0..n_sites)
        .map(|i| {
            let s0 = rec(
                i + 1,
                vec![obs(b"A", n_own, mean_err), obs(b"C", n_contam, mean_err)],
            );
            let s1 = rec(i + 1, vec![obs(b"A", n_reads, mean_err)]);
            let s2 = rec(
                i + 1,
                vec![obs(b"A", 0, mean_err), obs(b"C", n_reads, mean_err)],
            );
            Ok(pos_pileups(i + 1, vec![Some(s0), Some(s1), Some(s2)]))
        })
        .collect()
}

fn merged_record(
    chrom_id: u32,
    pos: u32,
    allele_seqs: Vec<&[u8]>,
    ploidy: u8,
    per_sample_scalars: Vec<Vec<AlleleSupportStats>>,
    per_sample_log_likelihoods: Vec<Vec<f64>>,
) -> MergedRecord {
    let alleles: Vec<MergedAllele> = allele_seqs
        .iter()
        .map(|s| MergedAllele {
            seq: s.to_vec(),
            is_compound: false,
            constituents: Vec::new(),
        })
        .collect();
    let n_samples = per_sample_scalars.len();
    let n_alleles = alleles.len();
    let n_genotypes = genotype_order(ploidy, n_alleles).len();
    let ref_len = alleles[0].seq.len() as u32;
    let scalars: Vec<AlleleSupportStats> = per_sample_scalars.into_iter().flatten().collect();
    assert_eq!(scalars.len(), n_samples * n_alleles);
    let other_scalars = vec![AlleleSupportStats::default(); n_samples];
    let chain_anchor_flags = vec![false; n_samples * n_alleles];
    let log_likelihoods: Vec<f64> = per_sample_log_likelihoods.into_iter().flatten().collect();
    assert_eq!(log_likelihoods.len(), n_samples * n_genotypes);
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

// ----- Tests -----------------------------------------------------------

/// Build a side-pass config in fixed-N mode with the project defaults
/// adjusted for the small synthetic cohorts the integration tests use
/// (the production floor of 5 samples/batch would suppress every sample).
fn small_cohort_fixed_n_config(num_sites: u32) -> ContaminationEstimationConfig {
    let mut cfg = ContaminationEstimationConfig::with_project_defaults();
    cfg.stopping_mode = StoppingMode::FixedSites { num_sites };
    cfg.block_size = 500;
    cfg.min_batch_size_for_contamination = 2;
    cfg
}

#[test]
fn side_pass_recovers_three_percent_contamination_end_to_end() {
    // Inject 3% contamination on sample 0; recover within tolerance.
    let stream = synth_three_sample_stream(4000, 50, 0.03);
    let estimates = estimate_contamination(
        stream.into_iter(),
        3,
        vec![0, 0, 0],
        1,
        small_cohort_fixed_n_config(3000),
    )
    .expect("side-pass should succeed");

    assert!(
        matches!(
            estimates.source,
            ContaminationEstimateSource::SidePass {
                sites_processed: 3000,
                ..
            }
        ),
        "unexpected source: {:?}",
        estimates.source
    );

    let c_s_0 = estimates.effective_c_s(0);
    assert!(
        (c_s_0 - 0.03).abs() < 0.02,
        "c_s for sample 0 = {}; expected ~0.03",
        c_s_0
    );
}

#[test]
fn side_pass_output_feeds_stage_6_without_error() {
    // The original "no error" smoke test, kept as-is: run the side-
    // pass, plug the ContaminationEstimates into a PosteriorEngine,
    // confirm emission succeeds and the posteriors are well-shaped.
    let stream = synth_three_sample_stream(2000, 50, 0.02);
    let estimates = estimate_contamination(
        stream.into_iter(),
        3,
        vec![0, 0, 0],
        1,
        small_cohort_fixed_n_config(1500),
    )
    .expect("side-pass should succeed");

    let alleles_seqs = vec![&b"A"[..], &b"C"[..]];
    let scalars_per_sample = vec![
        vec![
            AlleleSupportStats::new(49, 0.001_f64.ln() * 49.0, 24, 24, 0),
            AlleleSupportStats::new(1, 0.001_f64.ln(), 0, 0, 0),
        ];
        3
    ];
    let n_genotypes = genotype_order(2, 2).len();
    let lls_per_sample = vec![vec![0.0_f64; n_genotypes]; 3];
    let record = merged_record(0, 100, alleles_seqs, 2, scalars_per_sample, lls_per_sample);

    let mut config = PosteriorEngineConfig::default();
    config.contamination = Some(estimates);
    let upstream = std::iter::once(Ok::<MergedRecord, PerGroupMergerError>(record));
    let outputs: Vec<_> = PosteriorEngine::with_config(upstream, config).collect();
    assert_eq!(outputs.len(), 1);
    let pr: &PosteriorRecord = outputs[0].as_ref().expect("engine should emit");
    for s in 0..pr.n_samples {
        let row_sum: f64 = pr.posteriors_row(s).iter().sum();
        assert!(
            (row_sum - 1.0).abs() < 1e-9,
            "posterior row {} sums to {}",
            s,
            row_sum
        );
    }
}

#[test]
fn mixture_branch_shifts_posteriors_relative_to_no_contamination_baseline() {
    // B2 main assertion: feeding the engine a `ContaminationEstimates`
    // with a non-trivial `c_s` produces *different* posteriors than
    // the no-contamination path on the same record. A regression that
    // silently bypassed `compute_mixture_log_likelihoods` (e.g. a
    // wiring change making `effective_c_s` return 0) would make
    // these byte-identical. Uses `from_user_supplied` directly to
    // decouple the test from the side-pass's accuracy and pick a
    // c_s large enough to make the EM shift unambiguous.
    let estimates =
        ContaminationEstimates::from_user_supplied(vec![Some(0.5)], vec![[0.5, 0.5, 0.0]], vec![0])
            .expect("valid inputs");

    // Sample with very ambiguous data: 5 REF + 5 ALT reads. Without
    // contamination this looks like a clean heterozygote; with c_s=0.5
    // (every read 50/50 own-vs-contam) the mixture posterior is
    // visibly different. Small read count + lowered REF pseudocount
    // ensures the prior doesn't drown out the likelihood shift.
    let alleles_seqs = vec![&b"A"[..], &b"C"[..]];
    let scalars_per_sample = vec![vec![
        AlleleSupportStats::new(5, 0.001_f64.ln() * 5.0, 2, 2, 0),
        AlleleSupportStats::new(5, 0.001_f64.ln() * 5.0, 2, 2, 0),
    ]];
    // Heterozygote-leaning baseline lls (only used in the `c_s=0`
    // run; the mixture path recomputes from scalars).
    let lls_per_sample = vec![vec![-3.0_f64, 0.0, -3.0]];

    let record_contam = merged_record(
        0,
        100,
        alleles_seqs.clone(),
        2,
        scalars_per_sample.clone(),
        lls_per_sample.clone(),
    );
    let record_baseline =
        merged_record(0, 100, alleles_seqs, 2, scalars_per_sample, lls_per_sample);

    let mut config_contam = PosteriorEngineConfig::default();
    config_contam.contamination = Some(estimates);
    config_contam.ref_pseudocount = 0.1; // don't let pseudocount swamp 10-read signal
    let pr_contam: PosteriorRecord = PosteriorEngine::with_config(
        std::iter::once(Ok::<MergedRecord, PerGroupMergerError>(record_contam)),
        config_contam,
    )
    .next()
    .unwrap()
    .expect("engine should emit (contam)");

    let mut config_baseline = PosteriorEngineConfig::default();
    config_baseline.ref_pseudocount = 0.1;
    let pr_baseline: PosteriorRecord = PosteriorEngine::with_config(
        std::iter::once(Ok::<MergedRecord, PerGroupMergerError>(record_baseline)),
        config_baseline,
    )
    .next()
    .unwrap()
    .expect("engine should emit (baseline)");

    let mut any_diff = false;
    for s in 0..pr_contam.n_samples {
        for (a, b) in pr_contam
            .posteriors_row(s)
            .iter()
            .zip(pr_baseline.posteriors_row(s).iter())
        {
            if (a - b).abs() > 1e-6 {
                any_diff = true;
            }
        }
    }
    assert!(
        any_diff,
        "mixture branch produced output identical to c_s=0 path — \
         did `compute_mixture_log_likelihoods` get silently bypassed? \
         contam = {:?}, baseline = {:?}",
        pr_contam.posteriors_row(0),
        pr_baseline.posteriors_row(0),
    );
}

#[test]
fn fixed_n_mode_returns_insufficient_sites_with_recommendation() {
    let stream = synth_three_sample_stream(100, 50, 0.02);
    let result = estimate_contamination(
        stream.into_iter(),
        3,
        vec![0, 0, 0],
        1,
        small_cohort_fixed_n_config(10_000),
    );
    match result {
        Err(ContaminationEstimationError::InsufficientSites {
            requested,
            found,
            recommendation,
        }) => {
            assert_eq!(requested, 10_000);
            assert!(found <= 100);
            assert!(
                recommendation.suggested_num_sites.is_some() || !recommendation.message.is_empty()
            );
        }
        other => panic!("expected InsufficientSites, got {other:?}"),
    }
}

#[test]
fn user_supplied_estimates_round_trip_via_zero_constructor() {
    // ContaminationEstimates::zero(...) is the constructor the CLI
    // takes when --contamination-batches is not supplied; the engine
    // must accept it without erroring and produce posteriors
    // identical to the contamination=None path. This locks the
    // backwards-compatibility property.
    let alleles_seqs = vec![&b"A"[..], &b"C"[..]];
    let scalars_per_sample = vec![
        vec![
            AlleleSupportStats::new(50, 0.001_f64.ln() * 50.0, 25, 25, 0),
            AlleleSupportStats::default(),
        ];
        2
    ];
    let n_genotypes = genotype_order(2, 2).len();
    let lls_per_sample = vec![vec![0.0_f64, -50.0, -50.0]; 2];
    assert_eq!(lls_per_sample[0].len(), n_genotypes);

    let record_baseline = merged_record(
        0,
        100,
        alleles_seqs.clone(),
        2,
        scalars_per_sample.clone(),
        lls_per_sample.clone(),
    );
    let record_zero = merged_record(0, 100, alleles_seqs, 2, scalars_per_sample, lls_per_sample);

    let baseline: Vec<_> = PosteriorEngine::new(std::iter::once(Ok::<
        MergedRecord,
        PerGroupMergerError,
    >(record_baseline)))
    .collect();
    let mut config = PosteriorEngineConfig::default();
    config.contamination =
        Some(ContaminationEstimates::zero(2, vec![0, 0], 1).expect("valid inputs"));
    let mirrored: Vec<_> = PosteriorEngine::with_config(
        std::iter::once(Ok::<MergedRecord, PerGroupMergerError>(record_zero)),
        config,
    )
    .collect();

    match (&baseline[0], &mirrored[0]) {
        (Ok(b), Ok(m)) => assert_eq!(b, m),
        other => panic!("expected matching Ok records, got {other:?}"),
    }
}
