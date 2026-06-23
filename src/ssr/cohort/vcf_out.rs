//! VCF record formatting + the FP-control / output semantics (implementation plan
//! §3, C4 + E2): the allele-balance defence, emit-iff-variable, site
//! `QUAL = −10·log₁₀ P(monomorphic)`, the SSR FILTER vocabulary, and the apparent-
//! `F_IS` warning.
//!
//! `format_vcf_record` turns one [`LocusCall`] into a VCF data line
//! (`CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <samples>`), with `GT:GQ:REPCN`
//! per sample. `apply_fp_control` runs *first* — it penalises a depth-inflated false
//! het (a wildly imbalanced "het" is a homozygote + stutter/error, the SNP-path
//! blind spot) down to a low GQ / no-call. The site QUAL here is a posterior-based
//! proxy; the exact-AF convolution kernel (verify-fix #7b) is an F2 refinement.

use crate::ssr::cohort::candidate_set::{Admission, CandidateSet};
use crate::ssr::cohort::em::{LocusCall, SampleCall};
use crate::ssr::cohort::types::{CohortLocus, SampleEvidence};

/// FP-control + output thresholds (pinned in F2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FpControlCfg {
    /// A het whose minor-allele read fraction is below this is suspect (a
    /// depth-inflated false het); its GQ is scaled down toward zero.
    pub(crate) min_allele_balance: f64,
    /// A sample whose (penalised) GQ falls below this is converted to a no-call.
    pub(crate) no_call_gq: u8,
    /// Cap on the site QUAL.
    pub(crate) qual_cap: f64,
    /// Mean `F` above which the cohort gets an apparent-`F_IS` warning.
    pub(crate) f_is_warn_threshold: f64,
}

impl FpControlCfg {
    pub(crate) fn dev_default() -> Self {
        Self {
            min_allele_balance: 0.20,
            no_call_gq: 15,
            qual_cap: 200.0,
            f_is_warn_threshold: 0.50,
        }
    }
}

/// The FILTER column value for an admission verdict.
fn filter_text(admit: Admission) -> &'static str {
    match admit {
        Admission::Pass => "PASS",
        Admission::NotPeriodic => "notPeriodic",
        Admission::TooManyAlleles => "tooManyAlleles",
        Admission::LowDepth => "lowDepth",
    }
}

/// The minor-allele read fraction of a het call (`0` … `0.5`), by hard-attributing
/// each read to the nearer called allele length. A homozygote returns `1.0`.
fn allele_balance(evidence: &SampleEvidence, call: &SampleCall, period: usize) -> f64 {
    if call.genotype_units.len() < 2 || call.genotype_units[0] == call.genotype_units[1] {
        return 1.0;
    }
    let (a, b) = (call.genotype_units[0], call.genotype_units[1]);
    let (mut sa, mut sb) = (0u64, 0u64);
    for (obs, count) in &evidence.seq_counts {
        let units = (obs.len() / period) as i32;
        if (units - a as i32).abs() <= (units - b as i32).abs() {
            sa += *count as u64;
        } else {
            sb += *count as u64;
        }
    }
    let total = sa + sb;
    if total == 0 {
        return 1.0;
    }
    sa.min(sb) as f64 / total as f64
}

/// Apply the allele-balance FP defence in place: scale down an imbalanced het's GQ
/// and convert sub-threshold calls to no-calls.
pub(crate) fn apply_fp_control(locus: &CohortLocus, call: &mut LocusCall, cfg: &FpControlCfg) {
    let period = locus.motif.period();
    for (k, sample_call) in call.calls.iter_mut().enumerate() {
        if sample_call.allele_indices.is_empty() {
            continue; // already a no-call
        }
        let balance = allele_balance(&locus.samples[k], sample_call, period);
        if balance < cfg.min_allele_balance {
            // A depth-inflated false het: scale GQ by how far below the floor it is.
            let scale = (balance / cfg.min_allele_balance).clamp(0.0, 1.0);
            sample_call.gq = (sample_call.gq as f64 * scale).round() as u8;
        }
        if sample_call.gq < cfg.no_call_gq {
            *sample_call = SampleCall::no_call();
        }
    }
}

/// Whether the site is polymorphic: at least one present sample carries a
/// non-reference allele in its (post-FP-control) call.
pub(crate) fn is_variable(call: &LocusCall, candidates: &CandidateSet) -> bool {
    call.calls
        .iter()
        .flat_map(|c| c.allele_indices.iter())
        .any(|&idx| idx != candidates.ref_idx)
}

/// Site `QUAL = −10·log₁₀ P(monomorphic)` — a posterior proxy (the exact-AF kernel
/// is F2). `P(monomorphic)` is the product over samples of each being homozygous for
/// the cohort modal allele; a confident non-modal sample drives it toward 0 (high
/// QUAL).
pub(crate) fn site_qual(call: &LocusCall, candidates: &CandidateSet, cfg: &FpControlCfg) -> f64 {
    let modal = call
        .pi
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(candidates.ref_idx);

    let mut log_p_mono = 0.0_f64;
    for (k, sample_call) in call.calls.iter().enumerate() {
        if sample_call.allele_indices.is_empty() {
            continue; // a no-call contributes no evidence either way
        }
        let is_hom_modal = sample_call.allele_indices.iter().all(|&i| i == modal);
        let p_mono_s = if is_hom_modal {
            call.posterior_hom[k].clamp(1e-6, 1.0)
        } else {
            (1.0 - sample_call.posterior).clamp(1e-6, 1.0)
        };
        log_p_mono += p_mono_s.log10();
    }
    (-10.0 * log_p_mono).clamp(0.0, cfg.qual_cap)
}

/// An apparent-`F_IS` warning when the cohort mean `F` is implausibly high (almost
/// always a data artifact — relatedness, contamination, or mis-estimated stutter —
/// not real inbreeding). `None` if within range.
pub(crate) fn f_is_warning(f_per_sample: &[f64], cfg: &FpControlCfg) -> Option<String> {
    if f_per_sample.is_empty() {
        return None;
    }
    let mean = f_per_sample.iter().sum::<f64>() / f_per_sample.len() as f64;
    (mean > cfg.f_is_warn_threshold).then(|| {
        format!(
            "apparent F_IS unusually high (mean {mean:.3} > {:.2}); check for relatedness, \
             contamination, or mis-estimated stutter before reading it as inbreeding",
            cfg.f_is_warn_threshold
        )
    })
}

/// Format one locus call as a VCF data line.
///
/// `chrom` is the contig name (the work-item carries only a cohort-global id, so the
/// caller supplies the name). Allele numbering is VCF-standard: the reference
/// candidate is `0`, the remaining candidates become `ALT` alleles `1..` in
/// candidate order.
pub(crate) fn format_vcf_record(
    chrom: &str,
    locus: &CohortLocus,
    candidates: &CandidateSet,
    call: &LocusCall,
    qual: f64,
) -> String {
    let k = candidates.alleles.len();

    // Map each candidate index to its VCF allele number (ref → 0, others → 1..).
    let mut vcf_number = vec![0usize; k];
    let mut next = 1usize;
    for (idx, slot) in vcf_number.iter_mut().enumerate() {
        if idx == candidates.ref_idx {
            *slot = 0;
        } else {
            *slot = next;
            next += 1;
        }
    }

    let ref_allele = String::from_utf8_lossy(&candidates.alleles[candidates.ref_idx]);
    let alt: Vec<String> = (0..k)
        .filter(|&idx| idx != candidates.ref_idx)
        .map(|idx| String::from_utf8_lossy(&candidates.alleles[idx]).into_owned())
        .collect();
    let alt_field = if alt.is_empty() {
        ".".to_string()
    } else {
        alt.join(",")
    };

    let pos = locus.locus.start + 1; // VCF is 1-based
    let period = locus.motif.period();

    let mut sample_fields = Vec::with_capacity(call.calls.len());
    for sample in &call.calls {
        if sample.allele_indices.is_empty() {
            sample_fields.push("./.:.:.".to_string());
            continue;
        }
        let mut numbers: Vec<usize> = sample
            .allele_indices
            .iter()
            .map(|&c| vcf_number[c])
            .collect();
        numbers.sort_unstable();
        let gt = numbers
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join("/");
        let repcn = sample
            .genotype_units
            .iter()
            .map(u16::to_string)
            .collect::<Vec<_>>()
            .join(",");
        sample_fields.push(format!("{gt}:{}:{repcn}", sample.gq));
    }

    format!(
        "{chrom}\t{pos}\t.\t{ref_allele}\t{alt_field}\t{qual_field}\t{}\tPERIOD={period}\tGT:GQ:REPCN\t{}",
        filter_text(call.admit),
        sample_fields.join("\t"),
        qual_field = if qual > 0.0 {
            format!("{qual:.1}")
        } else {
            ".".to_string()
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::em::SampleCall;
    use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};
    use crate::ssr::types::Motif;

    fn ca_seq(units: u16) -> Box<[u8]> {
        std::iter::repeat_n(*b"CA", units as usize)
            .flatten()
            .collect()
    }

    fn locus() -> CohortLocus {
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 40,
                end: 56,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGG".as_slice()),
            ca_seq(8),
        );
        locus.push(
            0,
            SampleEvidence {
                seq_counts: vec![(ca_seq(8), 100)],
                qc: SsrQc::default(),
            },
        );
        locus.push(
            1,
            SampleEvidence {
                seq_counts: vec![(ca_seq(6), 50), (ca_seq(10), 50)],
                qc: SsrQc::default(),
            },
        );
        locus
    }

    fn candidates() -> CandidateSet {
        CandidateSet {
            // ref (8) at index 0; ALTs 6 and 10.
            alleles: vec![ca_seq(8), ca_seq(6), ca_seq(10)],
            ref_idx: 0,
            admit: Admission::Pass,
        }
    }

    #[test]
    fn formats_a_passing_record_with_gt_gq_repcn() {
        let call = LocusCall {
            calls: vec![
                SampleCall {
                    allele_indices: vec![0, 0], // hom ref (8/8)
                    genotype_units: vec![8, 8],
                    posterior: 0.99,
                    gq: 40,
                },
                SampleCall {
                    allele_indices: vec![1, 2], // het 6/10
                    genotype_units: vec![6, 10],
                    posterior: 0.95,
                    gq: 30,
                },
            ],
            pi: vec![0.5, 0.25, 0.25],
            posterior_hom: vec![0.99, 0.05],
            admit: Admission::Pass,
        };
        let line = format_vcf_record("chr1", &locus(), &candidates(), &call, 50.0);
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[0], "chr1");
        assert_eq!(cols[5], "50.0"); // QUAL
        assert_eq!(cols[1], "41"); // 40 + 1
        assert_eq!(cols[3], "CACACACACACACACA"); // REF = (CA)x8
        assert_eq!(cols[4], "CACACACACACA,CACACACACACACACACACA"); // ALT 6,10
        assert_eq!(cols[6], "PASS");
        assert_eq!(cols[7], "PERIOD=2");
        assert_eq!(cols[8], "GT:GQ:REPCN");
        assert_eq!(cols[9], "0/0:40:8,8");
        assert_eq!(cols[10], "1/2:30:6,10");
    }

    #[test]
    fn formats_a_no_call_filtered_record() {
        let call = LocusCall {
            calls: vec![SampleCall {
                allele_indices: vec![],
                genotype_units: vec![],
                posterior: 0.0,
                gq: 0,
            }],
            pi: vec![1.0],
            posterior_hom: vec![0.0],
            admit: Admission::LowDepth,
        };
        let cands = CandidateSet {
            alleles: vec![ca_seq(8)],
            ref_idx: 0,
            admit: Admission::LowDepth,
        };
        let line = format_vcf_record("chr1", &locus(), &cands, &call, 0.0);
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[5], "."); // QUAL absent
        assert_eq!(cols[4], "."); // no ALT
        assert_eq!(cols[6], "lowDepth");
        assert_eq!(cols[9], "./.:.:.");
    }

    // ── E2: FP control + output semantics ──

    fn het_call(gq: u8, posterior: f64) -> SampleCall {
        SampleCall {
            allele_indices: vec![1, 2], // 6/10 het
            genotype_units: vec![6, 10],
            posterior,
            gq,
        }
    }

    /// A sample whose reads are 95% at length 6 (so the 6/10 "het" is a depth-
    /// inflated false het) and a balanced 6/10 sample.
    fn imbalance_locus() -> CohortLocus {
        let mut l = locus(); // ref tract 8; samples replaced below
        l.present.clear();
        l.samples.clear();
        // Sample 0: 190 reads at 6, 10 at 10 → balance ≈ 0.05.
        l.push(
            0,
            SampleEvidence {
                seq_counts: vec![(ca_seq(6), 190), (ca_seq(10), 10)],
                qc: SsrQc::default(),
            },
        );
        // Sample 1: balanced 100/100 → balance ≈ 0.5.
        l.push(
            1,
            SampleEvidence {
                seq_counts: vec![(ca_seq(6), 100), (ca_seq(10), 100)],
                qc: SsrQc::default(),
            },
        );
        l
    }

    #[test]
    fn allele_balance_flags_an_imbalanced_het_and_keeps_a_balanced_one() {
        let l = imbalance_locus();
        let call = het_call(40, 0.99);
        assert!(allele_balance(&l.samples[0], &call, 2) < 0.1);
        assert!(allele_balance(&l.samples[1], &call, 2) > 0.45);
    }

    #[test]
    fn fp_control_no_calls_a_depth_inflated_false_het() {
        let l = imbalance_locus();
        let mut call = LocusCall {
            calls: vec![het_call(40, 0.99), het_call(40, 0.99)],
            pi: vec![0.2, 0.4, 0.4],
            posterior_hom: vec![0.0, 0.0],
            admit: Admission::Pass,
        };
        apply_fp_control(&l, &mut call, &FpControlCfg::dev_default());
        // The imbalanced het collapsed to a no-call; the balanced het survived.
        assert!(
            call.calls[0].allele_indices.is_empty(),
            "imbalanced het → no-call"
        );
        assert_eq!(call.calls[1].gq, 40, "balanced het keeps its GQ");
    }

    #[test]
    fn is_variable_true_with_an_alt_call_false_when_all_ref() {
        let cands = candidates();
        let variant = LocusCall {
            calls: vec![het_call(40, 0.9)],
            pi: vec![0.5, 0.25, 0.25],
            posterior_hom: vec![0.1],
            admit: Admission::Pass,
        };
        assert!(is_variable(&variant, &cands));

        let hom_ref = LocusCall {
            calls: vec![SampleCall {
                allele_indices: vec![0, 0],
                genotype_units: vec![8, 8],
                posterior: 0.99,
                gq: 40,
            }],
            pi: vec![1.0, 0.0, 0.0],
            posterior_hom: vec![0.99],
            admit: Admission::Pass,
        };
        assert!(!is_variable(&hom_ref, &cands));
    }

    #[test]
    fn site_qual_is_high_for_a_confident_variant_and_low_for_monomorphic() {
        let cfg = FpControlCfg::dev_default();
        let cands = candidates();

        let variant = LocusCall {
            calls: vec![het_call(40, 0.99), het_call(40, 0.99)],
            pi: vec![0.2, 0.4, 0.4],
            posterior_hom: vec![0.0, 0.0],
            admit: Admission::Pass,
        };
        let mono = LocusCall {
            calls: vec![
                SampleCall {
                    allele_indices: vec![0, 0],
                    genotype_units: vec![8, 8],
                    posterior: 0.99,
                    gq: 40,
                },
                SampleCall {
                    allele_indices: vec![0, 0],
                    genotype_units: vec![8, 8],
                    posterior: 0.99,
                    gq: 40,
                },
            ],
            pi: vec![1.0, 0.0, 0.0],
            posterior_hom: vec![0.99, 0.99],
            admit: Admission::Pass,
        };
        assert!(site_qual(&variant, &cands, &cfg) > site_qual(&mono, &cands, &cfg));
        assert!(site_qual(&variant, &cands, &cfg) > 20.0);
    }

    #[test]
    fn f_is_warning_fires_only_when_mean_f_is_implausible() {
        let cfg = FpControlCfg::dev_default();
        assert!(f_is_warning(&[0.05, 0.0, 0.1], &cfg).is_none());
        assert!(f_is_warning(&[0.8, 0.9, 0.7], &cfg).is_some());
    }
}
