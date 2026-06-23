//! Minimal VCF record formatting for `ssr-call` (implementation plan §3, C4; full
//! semantics — site QUAL, the SSR FILTER vocabulary, no-call refinement — are E2).
//!
//! For the walking skeleton this turns one [`LocusCall`] into a single VCF data line
//! (`CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <samples>`), with `GT:GQ:REPCN`
//! per sample. Site `QUAL` is `.` until E2 computes `−10·log₁₀ P(monomorphic)`.

use crate::ssr::cohort::candidate_set::{Admission, CandidateSet};
use crate::ssr::cohort::em::LocusCall;
use crate::ssr::cohort::types::CohortLocus;

/// The FILTER column value for an admission verdict.
fn filter_text(admit: Admission) -> &'static str {
    match admit {
        Admission::Pass => "PASS",
        Admission::NotPeriodic => "notPeriodic",
        Admission::TooManyAlleles => "tooManyAlleles",
        Admission::LowDepth => "lowDepth",
    }
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
        "{chrom}\t{pos}\t.\t{ref_allele}\t{alt_field}\t.\t{}\tPERIOD={period}\tGT:GQ:REPCN\t{}",
        filter_text(call.admit),
        sample_fields.join("\t")
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
        let line = format_vcf_record("chr1", &locus(), &candidates(), &call);
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[0], "chr1");
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
        let line = format_vcf_record("chr1", &locus(), &cands, &call);
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[4], "."); // no ALT
        assert_eq!(cols[6], "lowDepth");
        assert_eq!(cols[9], "./.:.:.");
    }
}
