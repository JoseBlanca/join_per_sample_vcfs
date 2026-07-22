//! ng locus generation — typed regions → a sample's loci (the shared shape).
//!
//! The middle arrow of the caller spine: region typing says *what the reference is*
//! at every position, read ingestion says *what reads a sample has*; this step joins
//! them into **loci**, a locus being a stretch of genome with one sample's evidence
//! attached. The evidence is the same kind of thing for every kind of locus — the
//! distinct sequences the reads showed, each with its support — so one type
//! ([`SampleLocusObservations`]) serves a candidate SNP, a microsatellite, and a
//! repeat cluster alike; [`LocusKind`] tags which.
//!
//! This module owns the **shared shape**: the locus type, the `LocusGenerator`
//! contract, the dispatcher, and the count-only `NoLoci` generator (landing across
//! this plan's milestones). Each real generator plugs in from its own module —
//! `ssr.rs` (STR), `pileup/` (generic). See `doc/devel/ng/spec/locus_generation.md`
//! (design) and `doc/devel/ng/arch/locus_generation.md` (types & interfaces).

use crate::ng::read::input::{IngestError, SampleReads};
use crate::ng::ref_seq::RefSeqError;
use crate::ng::region_typing::TypedRegionError;
use crate::ng::region_typing::segment_criteria::Motif;
use crate::ng::types::GenomeRegion;
use crate::pileup_record::ChainId;

/// One sample's locus: the stretch of genome it covers, and what that sample's reads
/// showed there.
///
/// **Owned, no lifetimes** — a cohort stage merges these across samples and a future
/// artifact writes them, so it must outlive every buffer it was built from
/// (spec §3). The evidence is uniform across kinds; [`kind`](Self::kind) is what names
/// the locus, not which fields are populated.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleLocusObservations {
    /// The stretch this locus covers — one base for a candidate SNP, several for an
    /// indel, the tract for a microsatellite.
    pub region: GenomeRegion,
    /// The reference (REF) bases over `region` — what a wider-span projection needs
    /// when samples merge.
    pub reference_bases: Box<[u8]>,
    /// The distinct sequences the reads showed, each with its support. **Observations,
    /// not alleles** — they become alleles when something calls them.
    pub observed_sequences: Vec<ObservedSequence>,
    /// Reads that covered this locus but produced no observation at all. A scalar with
    /// no positions: *no coverage* and *coverage that said nothing* are different
    /// states, and only one means "look at the mapping" (spec §3).
    pub reads_without_observation: u32,
    /// Reads a depth cap discarded. **Non-zero means the support counts are a
    /// subsample, not the depth** (spec §3).
    pub reads_discarded_by_cap: u32,
    /// What kind of locus this is, and the extras that kind carries — read off the
    /// type, never inferred from which fields are populated.
    pub kind: LocusKind,
}

impl SampleLocusObservations {
    /// Read depth at each position of `region`, in order — **derived, not stored**.
    ///
    /// A [`Complete`](ReadCoverage::Complete) observation counts its `num_obs` at every
    /// position; a [`PartialLeft(n)`](ReadCoverage::PartialLeft) at the leftmost `n`, a
    /// [`PartialRight(n)`](ReadCoverage::PartialRight) at the rightmost `n`. The returned
    /// vector has exactly `region.len()` entries.
    ///
    /// This is *observation* depth and only exact per locus: it omits reads that
    /// covered the tract but anchored no border (they are in
    /// [`reads_without_observation`](Self::reads_without_observation), a scalar with no
    /// positions), and overlapping loci on the generic path double-count if summed. The
    /// paralog filter owns both caveats (spec §3, §11).
    pub fn num_obs_along_locus(&self) -> Vec<u32> {
        let len = self.region.len() as usize;
        let mut depth = vec![0u32; len];
        for obs in &self.observed_sequences {
            // A partial's covered extent cannot exceed the locus length — that is a
            // producer invariant, enforced where `ReadCoverage` is minted (the STR
            // generator). Here we *clamp* rather than `debug_assert`: this is a
            // consumer-side derivation run over whole cohorts, and a debug-only guard
            // compiles out of the release build this repo actually runs (a trap it has
            // recorded hitting twice). Clamping keeps the derivation panic-free on any
            // input; a bad extent is caught at the producer, not here.
            let (from, to) = match obs.read_coverage {
                ReadCoverage::Complete => (0, len),
                ReadCoverage::PartialLeft(n) => (0, (n as usize).min(len)),
                ReadCoverage::PartialRight(n) => (len - (n as usize).min(len), len),
            };
            for slot in &mut depth[from..to] {
                *slot = slot.saturating_add(obs.num_obs);
            }
        }
        depth
    }

    /// The observations a likelihood may score directly — the
    /// [`Complete`](ReadCoverage::Complete) ones.
    ///
    /// A partial is a lower bound that mis-scores as a *short* allele until a censored
    /// likelihood models it (step 7), so reaching the partials is a deliberate act:
    /// this iterator is the guard (spec §3).
    pub fn complete_observations(&self) -> impl Iterator<Item = &ObservedSequence> + '_ {
        self.observed_sequences
            .iter()
            .filter(|obs| obs.read_coverage == ReadCoverage::Complete)
    }
}

/// One distinct sequence the reads showed at a locus, with its support.
///
/// The fields between `num_obs` and `chain_ids` are the per-read moments the SNP
/// filters read (strand bias, base-quality error, the MAPQ multi-mapper test); an STR
/// model reduces to `num_obs` alone. Modelled on production's per-allele shape
/// (`AlleleObservation` + `AlleleSupportStats`), minus the anchor-relative
/// read-position-bias fields, which are the generic path's (spec §3).
#[derive(Debug, Clone, PartialEq)]
pub struct ObservedSequence {
    /// The observed bases — allele content, in **read** coordinates.
    pub bases: Box<[u8]>,
    /// How much of the locus a read of this sequence spanned. **Part of the
    /// identity**: a [`Complete`](ReadCoverage::Complete) and a
    /// [`PartialLeft`](ReadCoverage::PartialLeft) of the same `bases` are different
    /// evidence and stay separate entries (spec §3).
    pub read_coverage: ReadCoverage,
    /// How many reads showed this sequence. The whole support on the STR path, and the
    /// one field every model on both paths reduces to.
    pub num_obs: u32,
    /// Reads on the forward strand — strand bias.
    pub num_fwd: u32,
    /// Σ per-read log-error over the supporting reads — the freebayes per-read error
    /// term (production's `q_sum`).
    pub q_sum: f64,
    /// Σ MAPQ over the supporting reads. With `mapq_sum_sq` and `num_obs` it recovers
    /// the mean and variance the MAPQ Welch's-t multi-mapper filter reads.
    pub mapq_sum: u32,
    /// Σ MAPQ² over the supporting reads (see `mapq_sum`).
    pub mapq_sum_sq: u64,
    /// Phase-chain ids of the reads folded here — what lets a later step chain
    /// observations at neighbouring loci into a haplotype.
    pub chain_ids: Vec<ChainId>,
}

/// How much of a locus a single read spanned — **one read's span, not depth**.
///
/// `Complete` means the read reached **both** borders of the locus; a partial ran off
/// its own end partway, carrying how many of the locus's positions it reached from
/// that border, in **locus** coordinates (spec §3). A partial is a *censored*
/// observation: the sequence is at least this long, but not how long.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReadCoverage {
    /// The read reached both borders of the locus.
    Complete,
    /// The read was flush with the left border and ran off its right end after
    /// covering the leftmost `n` positions of the locus.
    PartialLeft(u16),
    /// The read was flush with the right border and ran off its left end after
    /// covering the rightmost `n` positions of the locus.
    PartialRight(u16),
}

/// The kind of locus, plus whatever that kind adds to the shared evidence fields.
///
/// `#[non_exhaustive]` because a kind's payload fills in as the generator that mints
/// it is written. Shared vocabulary (`ng_step_interfaces.md` §2); authoritative in
/// spec §3.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum LocusKind {
    /// A SNP/indel candidate site. Its evidence is the observed alleles; no extras.
    Generic,
    /// A microsatellite tract — carries the motif and the flanks the read model needs.
    Ssr(SsrDetail),
    /// A repeat cluster with no clean flanks, coarser than a single tract. Its payload
    /// is the bundle generator's to decide (deferred, spec §11).
    SsrBundle,
}

/// What an [`LocusKind::Ssr`] locus carries — grouped so a repeat's motif and flanks
/// are present or absent together, never half.
///
/// The STR generator mints these, splitting the flanks out of the reference bases it
/// fetches (spec §3; `locus_generation_ssr.md` §1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SsrDetail {
    /// The repeat unit.
    pub motif: Motif,
    /// The reference flank left of the tract — the read model's alignment anchor.
    /// Clamped at the contig end, so it may be shorter than the right flank.
    pub left_flank: Box<[u8]>,
    /// The reference flank right of the tract (see `left_flank`).
    pub right_flank: Box<[u8]>,
}

/// Generates a sample's loci from one segment of kind `S`, streaming **one locus at a
/// time**.
///
/// `S` is the segment payload the generator consumes — `SsrSegment` for the STR generator.
/// It is a parameter on the *contract*, not an associated type inside each implementation,
/// so two generators for the same kind stay interchangeable behind `Box<dyn
/// LocusGenerator<S>>` (spec §4). A generator holds its own accessors (reference, aligner,
/// scratch) as fields, so the only per-call context is the segment and the sample's reads.
pub trait LocusGenerator<S> {
    /// Start a new segment: reset progress. Does no gathering and cannot fail.
    fn begin_segment(&mut self, region: GenomeRegion);

    /// The next locus of the segment begun, or `None` once it has no more.
    ///
    /// Called repeatedly with the same `segment` until it returns `None`; returning `None`
    /// immediately is a normal outcome, not a failure. The `segment` must be the one whose
    /// region was passed to the preceding [`begin_segment`](Self::begin_segment) — the two
    /// calls are paired, and nothing in the types enforces it. `&mut self` because a
    /// generator owns reusable scratch (alignment matrices, sampling buffers) that must not
    /// be reallocated per segment.
    fn next_locus(
        &mut self,
        segment: &S,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError>;
}

/// A generator that produces no loci and reports why.
///
/// One implementation covers every kind, because it ignores the segment entirely — the
/// count-only fallback every region kind with no real generator resolves to, so that "we
/// produce nothing here" is a configuration with a reason attached rather than a silent
/// gap (spec §5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NoLoci {
    /// Why this kind produces no loci — the fact the dispatcher accounts by.
    pub reason: UnhandledReason,
}

/// **Why** a kind produces no loci — a boundary we chose vs. a gap not yet filled.
///
/// Not cosmetic: the two answer different questions ("what will this caller never cover?"
/// vs "how much does it not cover *yet*?") and must not be added together (spec §5, §7).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnhandledReason {
    /// Deliberately outside the caller's scope — e.g. satellite arrays. Permanent.
    OutOfScope,
    /// No generator written yet. Temporary by construction.
    NotImplemented,
}

impl<S> LocusGenerator<S> for NoLoci {
    fn begin_segment(&mut self, _region: GenomeRegion) {}

    fn next_locus(
        &mut self,
        _segment: &S,
        _reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        Ok(None)
    }
}

/// A fatal, run-level failure of locus generation.
///
/// `#[non_exhaustive]`; every variant wraps an upstream error — this step mints none of its
/// own. A read that yields no observation is a tallied per-read outcome, never an error; an
/// error means the run is broken (spec §6). A reference fetch can surface two ways — through
/// the upstream walk ([`TypedRegion`](Self::TypedRegion)) or a generator's own fetch
/// ([`Reference`](Self::Reference)) — and they stay distinct because they fail in different
/// places.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum LocusGenerationError {
    /// The upstream typed-region walk failed.
    #[error("typed-region walk failed during locus generation")]
    TypedRegion(#[from] TypedRegionError),
    /// A read query failed mid-stream, or the alignment input was malformed (the open
    /// already succeeded upstream; this fires while a generator pulls reads).
    #[error("read access failed during locus generation")]
    Reads(#[from] IngestError),
    /// A reference fetch failed — a broken reference, or a region past a contig end.
    #[error("reference fetch failed during locus generation")]
    Reference(#[from] RefSeqError),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::{ContigId, Position};

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// An observation of `bases` with `num_obs` reads at a given coverage — the moment
    /// fields are irrelevant to the depth derivation, so they are fixed.
    fn obs(bases: &[u8], read_coverage: ReadCoverage, num_obs: u32) -> ObservedSequence {
        ObservedSequence {
            bases: Box::from(bases),
            read_coverage,
            num_obs,
            num_fwd: 0,
            q_sum: 0.0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            chain_ids: Vec::new(),
        }
    }

    fn locus(region: GenomeRegion, observed: Vec<ObservedSequence>) -> SampleLocusObservations {
        SampleLocusObservations {
            region,
            reference_bases: Box::from(&b""[..]),
            observed_sequences: observed,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// A minimal real `SampleReads` over the read-ingestion test fixtures — one indexed BAM
    /// naming one sample, opened against the fixture reference. Constructing `SampleReads`
    /// needs alignment files, so this is the cheapest honest handle; the `NoLoci` path never
    /// reads it, but the signature requires one. Returns the temp dirs so they outlive the
    /// handle.
    fn sample_reads_over_fixture() -> (tempfile::TempDir, tempfile::TempDir, SampleReads) {
        use crate::ng::read::filtering::ReadFilterConfig;
        use crate::ng::read::input::test_fixtures::{
            fixture_reference, header, indexed_bam, matching_contigs, read_named_with_length,
        };

        let (reference_dir, reference) = fixture_reference(false);
        let records = vec![read_named_with_length("r0", 0, 1, 30)];
        let (bam_dir, bam_path) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            &records,
        );
        let reads = SampleReads::open(&[bam_path], &reference, ReadFilterConfig::default(), false)
            .expect("the fixture sample opens");
        (reference_dir, bam_dir, reads)
    }

    /// `NoLoci` is a `LocusGenerator` for any segment type, emits no locus, and carries its
    /// reason for the dispatcher to account by (spec §5).
    #[test]
    fn no_loci_emits_nothing_and_carries_its_reason() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let mut generator = NoLoci {
            reason: UnhandledReason::OutOfScope,
        };
        // Driven over `()` as the segment — NoLoci ignores it, as it does every kind. The
        // segment type must be named because NoLoci implements the trait for *every* `S`.
        LocusGenerator::<()>::begin_segment(&mut generator, region(1, 5));
        let out = LocusGenerator::<()>::next_locus(&mut generator, &(), &reads).unwrap();
        assert!(out.is_none(), "NoLoci produces no locus");
        assert_eq!(generator.reason, UnhandledReason::OutOfScope);
    }

    /// The types compose into a locus of each kind — a smoke test that the shared
    /// shape holds together before the contract and dispatcher land on it.
    #[test]
    fn a_locus_of_each_kind_can_be_built() {
        let generic = SampleLocusObservations {
            region: region(100, 100),
            reference_bases: Box::from(&b"A"[..]),
            observed_sequences: vec![ObservedSequence {
                bases: Box::from(&b"T"[..]),
                read_coverage: ReadCoverage::Complete,
                num_obs: 9,
                num_fwd: 5,
                q_sum: -12.0,
                mapq_sum: 540,
                mapq_sum_sq: 32_400,
                chain_ids: vec![1, 2],
            }],
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        };
        assert_eq!(generic.region.len(), 1);
        assert_eq!(generic.observed_sequences[0].num_obs, 9);

        let ssr = SampleLocusObservations {
            region: region(10_442, 10_461),
            reference_bases: Box::from(&b"ATATATATATATATATATAT"[..]),
            observed_sequences: Vec::new(),
            reads_without_observation: 3,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Ssr(SsrDetail {
                motif: Motif::new(b"AT").unwrap(),
                left_flank: Box::from(&b"CCCGGG"[..]),
                right_flank: Box::from(&b"TTTAAA"[..]),
            }),
        };
        // Zero coverage is a real observation: the locus exists with an empty table.
        assert!(ssr.observed_sequences.is_empty());
        assert!(matches!(ssr.kind, LocusKind::Ssr(_)));

        let bundle = SampleLocusObservations {
            region: region(200, 260),
            reference_bases: Box::from(&b"N"[..]),
            observed_sequences: Vec::new(),
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::SsrBundle,
        };
        assert_eq!(bundle.kind, LocusKind::SsrBundle);
    }

    /// A complete and a partial observation of the *same* bases are distinct evidence,
    /// so they must not compare equal — the property the `(bases, read_coverage)`
    /// dedup key rests on (spec §3).
    #[test]
    fn same_bases_differ_by_read_coverage() {
        let complete = ObservedSequence {
            bases: Box::from(&b"ATATAT"[..]),
            read_coverage: ReadCoverage::Complete,
            num_obs: 1,
            num_fwd: 1,
            q_sum: 0.0,
            mapq_sum: 60,
            mapq_sum_sq: 3_600,
            chain_ids: Vec::new(),
        };
        let partial = ObservedSequence {
            read_coverage: ReadCoverage::PartialLeft(6),
            ..complete.clone()
        };
        assert_ne!(complete, partial);
        assert_ne!(complete.read_coverage, partial.read_coverage);
    }

    /// Depth derives correctly from read-coverage (spec §13.5): the vector has
    /// `region.len()` entries; a `Complete` raises every position, a `PartialLeft(n)`
    /// only the leftmost `n`, a `PartialRight(n)` only the rightmost `n`. A 10-position
    /// locus with one complete (×3), one left-partial reaching 4 (×2), one right-partial
    /// reaching 3 (×5).
    #[test]
    fn depth_derives_from_read_coverage() {
        let l = locus(
            region(1, 10),
            vec![
                obs(b"AAAAAAAAAA", ReadCoverage::Complete, 3),
                obs(b"AAAA", ReadCoverage::PartialLeft(4), 2),
                obs(b"AAA", ReadCoverage::PartialRight(3), 5),
            ],
        );
        // positions:            1  2  3  4  5  6  7  8  9 10
        //   complete ×3:        3  3  3  3  3  3  3  3  3  3
        //   left(4)  ×2:       +2 +2 +2 +2  .  .  .  .  .  .
        //   right(3) ×5:        .  .  .  .  .  .  . +5 +5 +5
        assert_eq!(l.num_obs_along_locus(), vec![5, 5, 5, 5, 3, 3, 3, 8, 8, 8],);
    }

    /// A single-base locus (a candidate SNP) has one depth position, raised by every
    /// complete observation over it.
    #[test]
    fn single_base_locus_has_one_depth_position() {
        let l = locus(
            region(42, 42),
            vec![
                obs(b"A", ReadCoverage::Complete, 7),
                obs(b"T", ReadCoverage::Complete, 2),
            ],
        );
        assert_eq!(l.num_obs_along_locus(), vec![9]);
    }

    /// No observations → depth is zero at every position, still `region.len()` long
    /// (the zero-coverage locus is a real one, not an absent one).
    #[test]
    fn no_observations_is_all_zero_full_length() {
        assert_eq!(
            locus(region(1, 4), Vec::new()).num_obs_along_locus(),
            vec![0; 4]
        );
    }

    /// A partial claiming to reach further than the locus is long is clamped, not an
    /// out-of-bounds index — the defensive guard, on **both** ends. The right arm's
    /// clamp is what keeps `len - n` from underflowing, so it is exercised too.
    #[test]
    fn partial_reach_beyond_locus_is_clamped() {
        let left = locus(
            region(1, 3),
            vec![obs(b"AAA", ReadCoverage::PartialLeft(9), 4)],
        );
        assert_eq!(left.num_obs_along_locus(), vec![4, 4, 4]);

        let right = locus(
            region(1, 3),
            vec![obs(b"AAA", ReadCoverage::PartialRight(9), 4)],
        );
        assert_eq!(right.num_obs_along_locus(), vec![4, 4, 4]);
    }

    /// `complete_observations()` yields only the complete entries — the guard that a
    /// partial (a lower bound) is never scored as an exact allele.
    #[test]
    fn complete_observations_excludes_partials() {
        let l = locus(
            region(1, 6),
            vec![
                obs(b"ATATAT", ReadCoverage::Complete, 4),
                obs(b"ATATAT", ReadCoverage::PartialLeft(6), 2),
                obs(b"ATGTAT", ReadCoverage::Complete, 3),
                obs(b"ATAT", ReadCoverage::PartialRight(4), 1),
            ],
        );
        // Both completes, and only the completes — a partial is never scored as exact.
        let complete: Vec<&[u8]> = l
            .complete_observations()
            .map(|o| o.bases.as_ref())
            .collect();
        assert_eq!(complete, vec![&b"ATATAT"[..], &b"ATGTAT"[..]]);
        assert_eq!(l.complete_observations().count(), 2);
    }
}
