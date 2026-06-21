//! Core Stage-2 (`ssr-call`) work-item types — the Phase-1 spine shared with the
//! parameter pre-pass and the genotyping EM (arch
//! `doc/devel/architecture/ssr_call_reading.md` §2).
//!
//! The cohort reader/merger turns N per-sample `.ssr.psp` files into a stream of
//! [`CohortLocus`] — one per catalog locus that at least one sample covers. The
//! types here are deliberately tiny per locus (a handful of distinct sequences per
//! sample); the cost of Stage 2 is decompression, not assembly.

use crate::ssr::types::Motif;

/// A cohort-wide locus coordinate key — the merge key shared by every sample file
/// and the catalog.
///
/// Coordinates are the **catalog's frame**: a cohort-global chromosome id and a
/// **0-based half-open** `[start, end)` tract. Each per-sample `.ssr.psp` record
/// (which is id-keyed in its own per-file dictionary and 1-based) is converted into
/// this frame at the cursor boundary, so the merger only ever compares `LocusId`s in
/// one frame.
///
/// `Ord` is lexicographic `(chrom_id, start, end)` — the ascending order the merger
/// walks the catalog in, and the order the per-sample cursor's monotonic guard
/// relies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct LocusId {
    /// Cohort-global chromosome id (assigned by the merger from catalog chromosome
    /// order — *not* any single file's per-file id).
    pub(crate) chrom_id: u32,
    /// 0-based start of the repeat tract (catalog frame).
    pub(crate) start: u32,
    /// Exclusive end of the repeat tract (catalog frame); the tract is `[start, end)`.
    pub(crate) end: u32,
}

/// The per-locus, per-sample QC scalars carried alongside the evidence. Mirrors the
/// lean QC set the Stage-1 tally records (`SsrLocusObs`'s QC fields); nothing here is
/// derivable from `seq_counts`, so it travels separately.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct SsrQc {
    /// Usable primary reads considered at the locus.
    pub(crate) depth: u32,
    /// Reads dropped by the admission gate (low-MAPQ / duplicate / qc-fail / short).
    pub(crate) n_filtered: u32,
    /// Dup-free primary coverage denominator.
    pub(crate) mapped_reads: u32,
    /// Reads dropped by the first-quartile quality gate.
    pub(crate) n_low_quality: u32,
    /// Reads whose flank ran off the read end (allele ≥ read length).
    pub(crate) n_border_off_end: u32,
}

/// One sample's evidence at one locus — the distinct repeat-tract sequences it
/// observed with their supporting read counts, plus the QC scalars.
///
/// This is the cohort-side mirror of the Stage-1 `SsrLocusObs` for a single present
/// sample: the cursor's decode adapter maps a container `SsrLocusRecord` into this
/// (its `observed` column becomes [`Self::seq_counts`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SampleEvidence {
    /// Distinct repeat-region sequences → supporting read count, **sorted by bytes**
    /// (the deterministic order the Stage-1 writer emits). Not alleles yet — raw
    /// candidate observations the EM will genotype.
    pub(crate) seq_counts: Vec<(Box<[u8]>, u32)>,
    /// Per-sample QC at this locus.
    pub(crate) qc: SsrQc,
}

/// The Stage-2 analysis unit: one catalog locus plus the evidence of every sample
/// that covered it.
///
/// The catalog supplies the **frame** ([`Self::locus`], [`Self::motif`],
/// [`Self::ref_frame`]) — a coordinate/alignment frame, **not** an allele claim. The
/// cohort evidence is stored **sparse, struct-of-arrays**: [`Self::samples`] holds
/// only the present samples, and [`Self::present`] is the parallel vector of their
/// cohort sample indices (`samples[k]` belongs to sample `present[k]`). Absent
/// samples are simply omitted; a `CohortLocus` is built only when ≥1 sample is
/// present (all-absent loci are dropped by the merger — sparse-omit).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CohortLocus {
    /// The cohort-match key (catalog frame).
    pub(crate) locus: LocusId,
    /// The repeat unit (the stutter unit the genotyper's stutter kernel uses).
    pub(crate) motif: Motif,
    /// Reference **tract + flanks** (not just the tract): the Δ frame and the
    /// alignment frame the genotyper aligns reads against.
    pub(crate) ref_frame: Box<[u8]>,
    /// Cohort sample indices of the present samples, ascending; parallel to
    /// [`Self::samples`].
    pub(crate) present: Vec<u32>,
    /// Evidence for each present sample; `samples[k]` belongs to sample `present[k]`.
    pub(crate) samples: Vec<SampleEvidence>,
}

impl CohortLocus {
    /// Start an empty work-item for `locus` with its catalog frame; the merger then
    /// [`push`](Self::push)es each present sample's evidence in ascending sample-index
    /// order.
    pub(crate) fn new(locus: LocusId, motif: Motif, ref_frame: Box<[u8]>) -> Self {
        Self {
            locus,
            motif,
            ref_frame,
            present: Vec::new(),
            samples: Vec::new(),
        }
    }

    /// Append one present sample's evidence. Callers must push in **ascending**
    /// `sample_idx` order so [`present`](Self::present) stays sorted (the merger asks
    /// cursors in index order, which already satisfies this).
    pub(crate) fn push(&mut self, sample_idx: u32, evidence: SampleEvidence) {
        debug_assert!(
            self.present.last().is_none_or(|&prev| prev < sample_idx),
            "CohortLocus::push requires ascending, distinct sample indices"
        );
        self.present.push(sample_idx);
        self.samples.push(evidence);
    }

    /// Number of present samples.
    pub(crate) fn present_count(&self) -> usize {
        self.present.len()
    }

    /// Whether no sample is present (the merger drops these — sparse-omit).
    pub(crate) fn is_empty(&self) -> bool {
        self.present.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn motif(bytes: &[u8]) -> Motif {
        Motif::new(bytes).expect("valid motif")
    }

    fn evidence(seqs: &[(&[u8], u32)]) -> SampleEvidence {
        SampleEvidence {
            seq_counts: seqs.iter().map(|(s, c)| (Box::from(*s), *c)).collect(),
            qc: SsrQc::default(),
        }
    }

    #[test]
    fn locus_id_orders_by_chrom_then_start_then_end() {
        let a = LocusId {
            chrom_id: 0,
            start: 100,
            end: 110,
        };
        let b = LocusId {
            chrom_id: 0,
            start: 100,
            end: 120,
        }; // same start, later end
        let c = LocusId {
            chrom_id: 0,
            start: 200,
            end: 205,
        }; // later start
        let d = LocusId {
            chrom_id: 1,
            start: 0,
            end: 5,
        }; // later chrom dominates

        assert!(a < b, "end breaks ties after chrom+start");
        assert!(b < c, "start dominates end");
        assert!(c < d, "chrom dominates start");
        assert!(a < d);
    }

    #[test]
    fn cohort_locus_starts_empty() {
        let cl = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 10,
                end: 16,
            },
            motif(b"AT"),
            Box::from(b"ATATAT".as_slice()),
        );
        assert!(cl.is_empty());
        assert_eq!(cl.present_count(), 0);
    }

    #[test]
    fn push_keeps_present_and_samples_parallel() {
        let mut cl = CohortLocus::new(
            LocusId {
                chrom_id: 2,
                start: 50,
                end: 58,
            },
            motif(b"AT"),
            Box::from(b"ATATATAT".as_slice()),
        );
        cl.push(0, evidence(&[(b"ATATAT", 7)]));
        cl.push(3, evidence(&[(b"ATATATAT", 5), (b"ATATAT", 2)]));

        assert!(!cl.is_empty());
        assert_eq!(cl.present_count(), 2);
        assert_eq!(cl.present, vec![0, 3]);
        assert_eq!(cl.samples.len(), 2);
        // samples[k] belongs to present[k]
        assert_eq!(cl.samples[0].seq_counts[0].1, 7);
        assert_eq!(cl.samples[1].seq_counts.len(), 2);
    }

    #[test]
    #[should_panic(expected = "ascending")]
    fn push_rejects_non_ascending_sample_index() {
        let mut cl = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 0,
                end: 4,
            },
            motif(b"AT"),
            Box::from(b"ATAT".as_slice()),
        );
        cl.push(5, evidence(&[(b"AT", 1)]));
        cl.push(5, evidence(&[(b"AT", 1)])); // not strictly greater → debug panic
    }
}
