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

use std::io::{self, Write};

use crate::psp::header::ParsedChromosome;
use crate::ssr::cohort::attribution::nearest_parent;
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

/// The non-PASS admission verdicts declared in the header, with their descriptions, in
/// declaration order. The IDs are taken from [`filter_text`] so the header declarations
/// and the per-record FILTER values stay one source of truth (PASS is implicit and not
/// declared, per VCF convention).
const FILTER_DESCRIPTIONS: &[(Admission, &str)] = &[
    (
        Admission::NotPeriodic,
        "locus allele-length distribution is inconsistent with the motif period",
    ),
    (
        Admission::TooManyAlleles,
        "more candidate alleles segregate than the caller admits",
    ),
    (Admission::LowDepth, "insufficient depth to call the locus"),
];

/// Write the SSR VCF header (arch `ssr_call_driver.md` §5): `##fileformat`, any cohort
/// warnings, one `##contig` per chromosome (the lengths come from the `.ssr.psp`
/// headers via [`CohortMerger::chromosomes`](crate::ssr::cohort::merge::CohortMerger::chromosomes)),
/// the `PERIOD` INFO, the `GT`/`GQ`/`REPCN` FORMATs, the SSR FILTER vocabulary, then the
/// `#CHROM … <samples>` column line in cohort order.
///
/// A dedicated SSR writer rather than the SNP `src/vcf/writer.rs`: the field vocabulary
/// differs and the data lines are the plain text [`format_vcf_record`] emits, so the
/// whole SSR VCF stays in one serialization style. `warnings` are emitted as
/// `##ssrCallWarning=…` meta lines (e.g. the apparent-`F_IS` notice, [`f_is_warning`]).
///
/// **Input contract:** each `warning` must be a single line (a newline would split it
/// into a bogus meta line — `f_is_warning` honours this); the caller is responsible for
/// VCF-clean sample/contig names (no tabs), and for **sample-name uniqueness** — the
/// driver validates that before calling here, since duplicate `#CHROM` columns are an
/// invalid VCF. `##contig` deliberately carries `ID` + `length` only; the `md5` on
/// `ParsedChromosome` is available but omitted (arch §5).
pub(crate) fn write_vcf_header<W: Write>(
    out: &mut W,
    chromosomes: &[ParsedChromosome],
    sample_names: &[String],
    warnings: &[String],
) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.4")?;
    for warning in warnings {
        writeln!(out, "##ssrCallWarning={warning}")?;
    }
    for chrom in chromosomes {
        writeln!(out, "##contig=<ID={},length={}>", chrom.name, chrom.length)?;
    }
    writeln!(
        out,
        "##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"Repeat unit length in bases\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scaled genotype quality\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=REPCN,Number=.,Type=Integer,Description=\"Repeat copy number of each called allele\">"
    )?;
    for (admit, description) in FILTER_DESCRIPTIONS {
        writeln!(
            out,
            "##FILTER=<ID={},Description=\"{description}\">",
            filter_text(*admit)
        )?;
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for name in sample_names {
        write!(out, "\t{name}")?;
    }
    writeln!(out)?;
    Ok(())
}

/// The minor-allele read fraction of a het call (`0` … `0.5`), by hard-attributing
/// each read to the nearer called allele length. A homozygote returns `1.0`.
fn allele_balance(evidence: &SampleEvidence, call: &SampleCall, period: usize) -> f64 {
    if call.genotype_units.len() < 2 || call.genotype_units[0] == call.genotype_units[1] {
        return 1.0;
    }
    let parents = [call.genotype_units[0], call.genotype_units[1]];
    let (mut sa, mut sb) = (0u64, 0u64);
    for (obs, count) in &evidence.seq_counts {
        let units = (obs.len() / period) as i32;
        // The shared nearest-parent rule breaks an equidistant tie toward the first allele,
        // matching the previous `<=` (review Mi1).
        let (parent_idx, _delta) =
            nearest_parent(units, &parents).expect("a het has two parent alleles");
        if parent_idx == 0 {
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
    // `total_cmp` is a total order over f64 (NaN-safe): `pi` is a normalized, finite
    // frequency vector from `run_pi_em`, so this is also panic-free, but `total_cmp` removes
    // the `partial_cmp().unwrap()` outright (review M4). Empty `pi` → the ref candidate.
    let modal = call
        .pi
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
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
    // `-10·log₁₀ P(mono)` is `-0.0` when `log_p_mono` is exactly 0 (a fully monomorphic
    // or all-no-call site); the trailing `+ 0.0` normalizes that negative zero so QUAL
    // prints `0.0`, never `-0.0` (which some VCF consumers reject).
    (-10.0 * log_p_mono).clamp(0.0, cfg.qual_cap) + 0.0
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
///
/// **Dense over the cohort (B1):** `call.calls` covers only the *present* samples
/// (`call.calls[k]` belongs to cohort sample `locus.present[k]`), but a VCF data line
/// must carry one column per cohort sample in cohort order. The row is built at width
/// `n_samples` (= the `#CHROM` column count), each present call placed at its
/// `locus.present[k]` position and every absent sample left as the `./.:.:.` no-call
/// placeholder.
pub(crate) fn format_vcf_record(
    chrom: &str,
    locus: &CohortLocus,
    candidates: &CandidateSet,
    call: &LocusCall,
    qual: f64,
    n_samples: usize,
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

    // PANIC-FREE: allele bytes are A/C/G/T/N-validated at the merge boundary (review Mi13),
    // so they are always valid ASCII/UTF-8 — `from_utf8` cannot error here, and a broken
    // upstream guarantee surfaces loudly rather than as a silent U+FFFD in REF/ALT.
    const ACGTN_GUARANTEE: &str = "alleles are A/C/G/T/N-validated at the merge boundary";
    let ref_bytes = candidates.alleles[candidates.ref_idx].as_ref();
    let alt_bytes: Vec<&[u8]> = (0..k)
        .filter(|&idx| idx != candidates.ref_idx)
        .map(|idx| candidates.alleles[idx].as_ref())
        .collect();

    // A length-zero allele (a complete tract deletion) cannot be written as an empty
    // REF/ALT field. When one is present, anchor the whole record on the reference base
    // just left of the tract (the VCF indel convention) so no field is ever empty; the
    // common case (no empty allele) still prints the tract sequences directly.
    let needs_anchor = ref_bytes.is_empty() || alt_bytes.iter().any(|a| a.is_empty());
    let mut pos = locus.locus.start + 1; // VCF is 1-based (tract start)
    let plain = |b: &[u8]| std::str::from_utf8(b).expect(ACGTN_GUARANTEE).to_string();
    let (ref_allele, alt): (String, Vec<String>) = if needs_anchor {
        // PANIC-FREE: the merger sets `left_anchor` for every real locus; it is `None`
        // only at a flankless contig-start locus, where a length-zero spanning-read
        // allele cannot arise. Fall back to `N` (a valid base) rather than emit an
        // invalid empty field, and only shift POS when a real left base was consumed.
        let anchor = locus.left_anchor.unwrap_or(b'N');
        if locus.left_anchor.is_some() {
            pos -= 1; // the record now starts at the anchor base
        }
        let anchored = |b: &[u8]| {
            let mut s = String::with_capacity(b.len() + 1);
            s.push(anchor as char);
            s.push_str(std::str::from_utf8(b).expect(ACGTN_GUARANTEE));
            s
        };
        (
            anchored(ref_bytes),
            alt_bytes.iter().map(|a| anchored(a)).collect(),
        )
    } else {
        (
            plain(ref_bytes),
            alt_bytes.iter().map(|a| plain(a)).collect(),
        )
    };
    let alt_field = if alt.is_empty() {
        ".".to_string()
    } else {
        alt.join(",")
    };

    let period = locus.motif.period();

    // Dense over the cohort: every absent sample stays a `./.:.:.` no-call; each present
    // call lands in its own cohort column via `locus.present` (B1).
    let mut sample_fields = vec!["./.:.:.".to_string(); n_samples];
    for (k, sample) in call.calls.iter().enumerate() {
        if sample.allele_indices.is_empty() {
            continue; // a present-but-no-call sample keeps the placeholder
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
        sample_fields[locus.present[k] as usize] = format!("{gt}:{}:{repcn}", sample.gq);
    }

    // Every emitted record carries a numeric QUAL — `site_qual` always returns a finite
    // value in `[0, qual_cap]`, so a clamped-to-zero variable locus prints `0.0`, not `.`
    // (downstream QUAL filters then read "lowest", not "missing"; review Mi14). `.` is
    // reserved for a genuinely unscored site, which this path never produces.
    format!(
        "{chrom}\t{pos}\t.\t{ref_allele}\t{alt_field}\t{qual:.1}\t{}\tPERIOD={period}\tGT:GQ:REPCN\t{}",
        filter_text(call.admit),
        sample_fields.join("\t"),
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
                    allele_support: Default::default(),
                },
                SampleCall {
                    allele_indices: vec![1, 2], // het 6/10
                    genotype_units: vec![6, 10],
                    posterior: 0.95,
                    gq: 30,
                    allele_support: Default::default(),
                },
            ],
            pi: vec![0.5, 0.25, 0.25],
            posterior_hom: vec![0.99, 0.05],
            admit: Admission::Pass,
        };
        let line = format_vcf_record("chr1", &locus(), &candidates(), &call, 50.0, 2);
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
    fn zero_unit_allele_is_anchored_so_no_field_is_empty() {
        // A complete tract deletion (a 0-unit allele) would otherwise emit an empty ALT —
        // an invalid VCF line. With the merger's left-anchor base set, the record anchors
        // on the reference base just left of the tract: REF and ALT both carry it and POS
        // shifts left by one. (Real-data regression: SL4.0ch01:32961838.)
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
        locus.left_anchor = Some(b'G'); // the merger sets this from the catalog frame
        locus.push(
            0,
            SampleEvidence {
                seq_counts: vec![(ca_seq(8), 50), (Box::from(&b""[..]), 50)],
                qc: SsrQc::default(),
            },
        );
        let candidates = CandidateSet {
            alleles: vec![ca_seq(8), Box::from(&b""[..])], // ref (8 units) + a 0-unit deletion
            ref_idx: 0,
            admit: Admission::Pass,
        };
        let call = LocusCall {
            calls: vec![SampleCall {
                allele_indices: vec![0, 1], // het: reference vs full deletion
                genotype_units: vec![8, 0],
                posterior: 0.9,
                gq: 30,
                allele_support: Default::default(),
            }],
            pi: vec![0.5, 0.5],
            posterior_hom: vec![0.1],
            admit: Admission::Pass,
        };
        let line = format_vcf_record("chr1", &locus, &candidates, &call, 42.0, 1);
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[1], "40"); // POS shifted left to the anchor (41 - 1)
        assert_eq!(cols[3], "GCACACACACACACACA"); // REF = anchor + (CA)x8
        assert_eq!(cols[4], "G"); // ALT = anchor only (the deletion) — NOT empty
        assert!(!cols[4].is_empty(), "ALT must never be an empty field");
        assert_eq!(cols[9], "0/1:30:8,0");
    }

    #[test]
    fn site_qual_normalizes_negative_zero() {
        // An all-no-call (or fully monomorphic) site gives log_p_mono == 0, so
        // -10*log_p_mono is -0.0; QUAL must print 0.0, never -0.0 (real-data regression:
        // every lowDepth filtered locus hit this).
        let cands = CandidateSet {
            alleles: vec![ca_seq(8)],
            ref_idx: 0,
            admit: Admission::LowDepth,
        };
        let call = LocusCall {
            calls: vec![SampleCall {
                allele_indices: vec![],
                genotype_units: vec![],
                posterior: 0.0,
                gq: 0,
                allele_support: Default::default(),
            }],
            pi: vec![1.0],
            posterior_hom: vec![0.0],
            admit: Admission::LowDepth,
        };
        let q = site_qual(&call, &cands, &FpControlCfg::dev_default());
        assert_eq!(q, 0.0);
        assert!(q.is_sign_positive(), "QUAL must be +0.0, not -0.0");
        assert_eq!(format!("{q:.1}"), "0.0");
    }

    #[test]
    fn formats_a_no_call_filtered_record() {
        let call = LocusCall {
            calls: vec![SampleCall {
                allele_indices: vec![],
                genotype_units: vec![],
                posterior: 0.0,
                gq: 0,
                allele_support: Default::default(),
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
        let line = format_vcf_record("chr1", &locus(), &cands, &call, 0.0, 1);
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[5], "0.0"); // QUAL clamped to zero is printed numerically (Mi14)
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
            allele_support: Default::default(),
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
                allele_support: Default::default(),
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
                    allele_support: Default::default(),
                },
                SampleCall {
                    allele_indices: vec![0, 0],
                    genotype_units: vec![8, 8],
                    posterior: 0.99,
                    gq: 40,
                    allele_support: Default::default(),
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

    // ── H3: the VCF header writer ──

    fn contig(name: &str, length: u32) -> ParsedChromosome {
        ParsedChromosome {
            name: name.to_string(),
            length,
            md5: "0".repeat(32),
        }
    }

    fn header_to_string(
        chroms: &[ParsedChromosome],
        samples: &[String],
        warnings: &[String],
    ) -> String {
        let mut out = Vec::new();
        write_vcf_header(&mut out, chroms, samples, warnings).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn write_vcf_header_emits_all_sections() {
        let text = header_to_string(
            &[contig("chr1", 248_956_422), contig("chr2", 242_193_529)],
            &["sampleA".to_string(), "sampleB".to_string()],
            &["apparent F_IS unusually high".to_string()],
        );
        let lines: Vec<&str> = text.lines().collect();

        assert_eq!(lines[0], "##fileformat=VCFv4.4");
        assert!(text.contains("##ssrCallWarning=apparent F_IS unusually high"));
        assert!(text.contains("##contig=<ID=chr1,length=248956422>"));
        assert!(text.contains("##contig=<ID=chr2,length=242193529>"));
        assert!(text.contains("##INFO=<ID=PERIOD,"));
        assert!(text.contains("##FORMAT=<ID=GT,"));
        assert!(text.contains("##FORMAT=<ID=GQ,"));
        assert!(text.contains("##FORMAT=<ID=REPCN,"));
        // The column line is last, tab-separated, with the samples in order.
        assert_eq!(
            lines.last().unwrap(),
            &"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB"
        );
    }

    #[test]
    fn write_vcf_header_filter_ids_match_filter_text_and_omit_pass() {
        let text = header_to_string(&[contig("chr1", 1000)], &["s".to_string()], &[]);
        // Every non-PASS admission verdict is declared, with the same ID filter_text emits.
        for admit in [
            Admission::NotPeriodic,
            Admission::TooManyAlleles,
            Admission::LowDepth,
        ] {
            assert!(
                text.contains(&format!("##FILTER=<ID={},", filter_text(admit))),
                "missing FILTER for {admit:?}"
            );
        }
        // PASS is implicit, never declared.
        assert!(!text.contains("##FILTER=<ID=PASS"));
    }

    #[test]
    fn write_vcf_header_without_warnings_emits_no_warning_line() {
        let text = header_to_string(&[contig("chr1", 1000)], &["s".to_string()], &[]);
        assert!(!text.contains("##ssrCallWarning="));
    }
}
