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
}
