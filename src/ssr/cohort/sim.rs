//! A2 — the cohort simulator: synthetic cohorts with **known** genotypes, stutter
//! (shape × level), `ε`, and sample groups, emitted as the `.ssr.psp` bytes the
//! Stage-2 reader consumes (implementation plan §A2; arch `ssr_call_parameters.md`
//! §1 — the generative side of the same model the kernel scores).
//!
//! This is the ground truth every later step recovers: a generated cohort
//! round-trips through [`CohortMerger`](super::merge::CohortMerger) to
//! [`CohortLocus`](super::types::CohortLocus), and its [`TruthTable`] answers "what
//! did we put in?" for the recover-what-we-put-in asserts of milestones C–E.
//!
//! It supports per-**group** stutter *shape* divergence (not just level), the
//! mononucleotide / C1 case, and well-separated-het cohorts (the CG-seed) — the
//! cases D3/M3 and the genotyper need exercised.
//!
//! **Forward model** (one read of a parent allele of `units` repeat units):
//! 1. **slip** — with probability `level(units) = baseline + slope·units` (clamped
//!    to `[0,1]`) the read slips; otherwise it is faithful (`Δ = 0`). On a slip the
//!    sign is `+` with probability `up_rate / (up_rate + down_rate)`, and the
//!    magnitude `m ≥ 1` is geometric with continuation probability `decay` (capped
//!    at [`MAX_SLIP`]). The read length is `clamp(units ± m, 1, …)` units.
//! 2. **substitution** — each base of the realized tract is, independently with
//!    probability `ε`, replaced by a different uniform base (within-tract
//!    substitutions only — C1; length is set by the slip alone).
//!
//! Determinism: every (sample, locus) cell draws from a [`SplitMix64`] seeded from
//! `(global_seed, sample_idx, locus_idx)`, so the whole cohort is a pure function
//! of the spec — same spec ⇒ byte-identical `.ssr.psp`.

use std::collections::BTreeMap;
use std::io::Cursor;

use crate::psp::PspReader;
use crate::psp::registry_ssr::SsrLocusRecord;
use crate::ssr::catalog::io::CatalogReader;
use crate::ssr::cohort::merge::CohortMerger;
use crate::ssr::cohort::param_estimation::{
    MAX_SLIP, PerBaseError, SampleGroupId, StutterLevel, StutterShape,
};
use crate::ssr::cohort::test_support::{REF_MD5, catalog_reader, reader, ssr_header, ssr_psp};
use crate::ssr::types::{Locus, Motif};

/// A minimal deterministic PRNG (SplitMix64) — dependency-free and reproducible, so
/// the simulator's determinism does not rest on an external crate's stability.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A float in `[0, 1)` with 53 bits of mantissa.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A Bernoulli draw — true with probability `p`.
    fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }

    /// A uniform index in `0..n` (`n > 0`).
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// The seed for one (sample, locus) cell — mixes the global seed with both indices
/// so cells are independent yet reproducible.
fn cell_seed(global: u64, sample_idx: usize, locus_idx: usize) -> u64 {
    let mixed = global
        ^ (sample_idx as u64).wrapping_mul(0xD1B5_4A32_D192_ED03)
        ^ (locus_idx as u64).wrapping_mul(0xA076_1D64_78BD_642F);
    // One round of SplitMix64 to decorrelate the mixed bits.
    SplitMix64::new(mixed).next_u64()
}

/// The generative chemistry of one sample group: `ε`, stutter shape, stutter level.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct SimChemistry {
    pub(crate) error: PerBaseError,
    pub(crate) shape: StutterShape,
    pub(crate) level: StutterLevel,
}

/// The true genotype at one (sample, locus): each allele copy's length in **repeat
/// units** (`allele_units.len()` = ploidy; a homozygote repeats a value).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SimGenotype {
    pub(crate) allele_units: Vec<u16>,
}

impl SimGenotype {
    pub(crate) fn homozygous(units: u16, ploidy: usize) -> Self {
        Self {
            allele_units: vec![units; ploidy],
        }
    }

    pub(crate) fn diploid(a: u16, b: u16) -> Self {
        Self {
            allele_units: vec![a, b],
        }
    }
}

/// One locus in the simulated catalog. The realized tract is `motif` tiled
/// `ref_units` times; the reference frame embeds a fixed 6 bp flank each side, so
/// `start` must be `>= 6`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SimLocus {
    pub(crate) chrom: String,
    pub(crate) start: u32,
    pub(crate) motif: Motif,
    pub(crate) ref_units: u16,
}

const LEFT_FLANK: &[u8] = b"GGGGGG";
const RIGHT_FLANK: &[u8] = b"TTTTTT";

/// One simulated sample: its sample group, and a per-locus genotype (`None` = the
/// sample is absent at that locus). `genotypes` is parallel to the spec's loci.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SimSample {
    pub(crate) name: String,
    pub(crate) group: SampleGroupId,
    pub(crate) genotypes: Vec<Option<SimGenotype>>,
}

/// The full description of a synthetic cohort to generate.
#[derive(Debug, Clone)]
pub(crate) struct SimCohortSpec {
    /// The global RNG seed — the whole cohort is a pure function of it + the spec.
    pub(crate) seed: u64,
    /// Loci, **in catalog order** (sorted by chromosome then start).
    pub(crate) loci: Vec<SimLocus>,
    /// Per-group chemistry, indexed by [`SampleGroupId`].
    pub(crate) groups: Vec<SimChemistry>,
    /// The samples.
    pub(crate) samples: Vec<SimSample>,
    /// Reads generated per present (sample, locus).
    pub(crate) depth: u32,
}

/// What was put in — queried by the recovery asserts of later milestones.
#[derive(Debug, Clone)]
pub(crate) struct TruthTable {
    /// Sample index → its sample group.
    pub(crate) group_of_sample: Vec<SampleGroupId>,
    /// Per-group chemistry (the `ε` / shape / level the reads were drawn under).
    pub(crate) chemistry: Vec<SimChemistry>,
    /// `genotypes[sample][locus]` — the true genotype, or `None` if absent.
    pub(crate) genotypes: Vec<Vec<Option<SimGenotype>>>,
}

impl TruthTable {
    /// The true genotype at `(sample_idx, locus_idx)`, if the sample is present.
    pub(crate) fn genotype(&self, sample_idx: usize, locus_idx: usize) -> Option<&SimGenotype> {
        self.genotypes[sample_idx][locus_idx].as_ref()
    }

    /// The chemistry the given sample's reads were generated under.
    pub(crate) fn chemistry_of(&self, sample_idx: usize) -> &SimChemistry {
        &self.chemistry[self.group_of_sample[sample_idx].0 as usize]
    }
}

/// A generated cohort: the realized catalog loci, each sample's `.ssr.psp` bytes,
/// and the truth table.
#[derive(Debug, Clone)]
pub(crate) struct SimCohort {
    /// The realized catalog loci (catalog order).
    pub(crate) catalog_loci: Vec<Locus>,
    /// Per-sample `(name, .ssr.psp bytes)`, parallel to the spec's samples.
    pub(crate) sample_psp: Vec<(String, Vec<u8>)>,
    /// Distinct chromosome names, sorted — the per-file chromosome frame.
    pub(crate) chroms: Vec<String>,
    /// What was put in.
    pub(crate) truth: TruthTable,
}

impl SimCohort {
    /// Open a [`CohortMerger`] over freshly-cloned in-memory readers — the real
    /// reader path, so a round-trip test exercises decode + merge, not a shortcut.
    pub(crate) fn merger(&self) -> CohortMerger<Cursor<Vec<u8>>, Cursor<Vec<u8>>> {
        let catalog: CatalogReader<Cursor<Vec<u8>>> = catalog_reader(REF_MD5, &self.catalog_loci);
        let samples: Vec<(String, PspReader<Cursor<Vec<u8>>>)> = self
            .sample_psp
            .iter()
            .map(|(name, bytes)| (name.clone(), reader(bytes.clone())))
            .collect();
        CohortMerger::from_parts(catalog, samples).expect("simulated cohort opens")
    }
}

/// Tile `motif` `units` times into a tract sequence.
fn build_tract(motif: &Motif, units: u16) -> Vec<u8> {
    let m = motif.as_bytes();
    let mut tract = Vec::with_capacity(m.len() * units as usize);
    for _ in 0..units {
        tract.extend_from_slice(m);
    }
    tract
}

/// Realize one [`SimLocus`] into a catalog [`Locus`] (fixed 6 bp flanks).
fn realized_locus(sl: &SimLocus) -> Locus {
    let tract = build_tract(&sl.motif, sl.ref_units);
    let mut ref_bytes = Vec::with_capacity(LEFT_FLANK.len() + tract.len() + RIGHT_FLANK.len());
    ref_bytes.extend_from_slice(LEFT_FLANK);
    ref_bytes.extend_from_slice(&tract);
    ref_bytes.extend_from_slice(RIGHT_FLANK);
    let start = sl.start;
    let end = start + tract.len() as u32;
    let ref_bytes_start = start - LEFT_FLANK.len() as u32;
    Locus::new(
        sl.chrom.clone().into(),
        start,
        end,
        sl.motif,
        1.0,
        ref_bytes.into_boxed_slice(),
        ref_bytes_start,
    )
    .expect("simulated locus is well-formed")
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Draw a read's realized repeat-unit length given its parent allele's `units`.
fn slip_length(rng: &mut SplitMix64, chem: &SimChemistry, units: u16) -> u16 {
    let level = (chem.level.baseline + chem.level.slope * units as f64).clamp(0.0, 1.0);
    if !rng.chance(level) {
        return units; // faithful
    }
    let up_share = {
        let denom = chem.shape.up_rate + chem.shape.down_rate;
        if denom > 0.0 {
            chem.shape.up_rate / denom
        } else {
            0.5
        }
    };
    let sign: i32 = if rng.chance(up_share) { 1 } else { -1 };
    let mut magnitude: i32 = 1;
    while (magnitude as usize) < MAX_SLIP && rng.chance(chem.shape.decay) {
        magnitude += 1;
    }
    let slipped = units as i32 + sign * magnitude;
    slipped.clamp(1, u16::MAX as i32) as u16
}

/// Apply within-tract substitutions at rate `ε`.
fn apply_substitutions(rng: &mut SplitMix64, tract: &mut [u8], eps: f64) {
    if eps <= 0.0 {
        return;
    }
    for base in tract.iter_mut() {
        if rng.chance(eps) {
            let mut replacement = BASES[rng.below(BASES.len())];
            while replacement == *base {
                replacement = BASES[rng.below(BASES.len())];
            }
            *base = replacement;
        }
    }
}

/// Generate one (sample, locus) cell's distinct observed sequences, ascending by
/// bytes (the Stage-1 writer's contract). A `BTreeMap` keyed by the sequence keeps
/// that order for free.
fn simulate_cell(
    seed: u64,
    chem: &SimChemistry,
    genotype: &SimGenotype,
    motif: &Motif,
    depth: u32,
) -> Vec<(Box<[u8]>, u32)> {
    let mut rng = SplitMix64::new(seed);
    let ploidy = genotype.allele_units.len();
    let mut counts: BTreeMap<Box<[u8]>, u32> = BTreeMap::new();
    for _ in 0..depth {
        let parent = genotype.allele_units[rng.below(ploidy)];
        let len = slip_length(&mut rng, chem, parent);
        let mut tract = build_tract(motif, len);
        apply_substitutions(&mut rng, &mut tract, chem.error.0);
        *counts.entry(tract.into_boxed_slice()).or_insert(0) += 1;
    }
    counts.into_iter().collect()
}

/// Generate a cohort from its spec.
pub(crate) fn simulate(spec: &SimCohortSpec) -> SimCohort {
    // The per-file chromosome frame: distinct catalog chromosome names, sorted.
    let mut chroms: Vec<String> = spec.loci.iter().map(|l| l.chrom.clone()).collect();
    chroms.sort();
    chroms.dedup();
    let chrom_id = |name: &str| chroms.iter().position(|c| c == name).unwrap() as u32;

    // Loci must arrive in catalog order so the per-sample records stay monotone.
    debug_assert!(
        spec.loci
            .windows(2)
            .all(|w| (chrom_id(&w[0].chrom), w[0].start) < (chrom_id(&w[1].chrom), w[1].start)),
        "SimCohortSpec.loci must be sorted by (chromosome, start)"
    );

    let catalog_loci: Vec<Locus> = spec.loci.iter().map(realized_locus).collect();

    let chrom_refs: Vec<&str> = chroms.iter().map(String::as_str).collect();
    let sample_psp: Vec<(String, Vec<u8>)> = spec
        .samples
        .iter()
        .enumerate()
        .map(|(s_idx, sample)| {
            let chem = &spec.groups[sample.group.0 as usize];
            let records: Vec<SsrLocusRecord> = spec
                .loci
                .iter()
                .enumerate()
                .filter_map(|(l_idx, locus)| {
                    let genotype = sample.genotypes[l_idx].as_ref()?;
                    let observed = simulate_cell(
                        cell_seed(spec.seed, s_idx, l_idx),
                        chem,
                        genotype,
                        &locus.motif,
                        spec.depth,
                    );
                    let depth: u32 = observed.iter().map(|(_, c)| c).sum();
                    let tract_units_len = build_tract(&locus.motif, locus.ref_units).len() as u32;
                    Some(SsrLocusRecord {
                        chrom_id: chrom_id(&locus.chrom),
                        start: locus.start + 1, // container coords are 1-based
                        end: locus.start + 1 + tract_units_len,
                        depth,
                        n_filtered: 0,
                        mapped_reads: depth,
                        n_low_quality: 0,
                        n_border_off_end: 0,
                        observed,
                    })
                })
                .collect();
            let bytes = ssr_psp(ssr_header(&chrom_refs, REF_MD5), &records);
            (sample.name.clone(), bytes)
        })
        .collect();

    let truth = TruthTable {
        group_of_sample: spec.samples.iter().map(|s| s.group).collect(),
        chemistry: spec.groups.clone(),
        genotypes: spec.samples.iter().map(|s| s.genotypes.clone()).collect(),
    };

    SimCohort {
        catalog_loci,
        sample_psp,
        chroms,
        truth,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::types::CohortLocus;

    fn motif(bytes: &[u8]) -> Motif {
        Motif::new(bytes).unwrap()
    }

    /// A clean two-group cohort: group 0 low stutter, group 1 higher stutter; two
    /// di-nucleotide loci; sample A a homozygote, sample B a separated het.
    fn clean_spec(seed: u64) -> SimCohortSpec {
        let low = SimChemistry {
            error: PerBaseError(0.001),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.03,
                slope: 0.0,
            },
        };
        let high = SimChemistry {
            error: PerBaseError(0.002),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.2,
            },
            level: StutterLevel {
                baseline: 0.10,
                slope: 0.0,
            },
        };
        SimCohortSpec {
            seed,
            loci: vec![
                SimLocus {
                    chrom: "chr1".into(),
                    start: 20,
                    motif: motif(b"CA"),
                    ref_units: 6,
                },
                SimLocus {
                    chrom: "chr1".into(),
                    start: 100,
                    motif: motif(b"CA"),
                    ref_units: 6,
                },
            ],
            groups: vec![low, high],
            samples: vec![
                SimSample {
                    name: "A".into(),
                    group: SampleGroupId(0),
                    genotypes: vec![
                        Some(SimGenotype::homozygous(6, 2)),
                        Some(SimGenotype::homozygous(8, 2)),
                    ],
                },
                SimSample {
                    name: "B".into(),
                    group: SampleGroupId(1),
                    genotypes: vec![
                        // separated het: 4 vs 9 units (≥2 apart)
                        Some(SimGenotype::diploid(4, 9)),
                        None, // absent at the second locus
                    ],
                },
            ],
            depth: 80,
        }
    }

    fn collect(cohort: &SimCohort) -> Vec<(u64, CohortLocus)> {
        cohort
            .merger()
            .map(|item| item.expect("merge yields a locus"))
            .collect()
    }

    #[test]
    fn generated_cohort_round_trips_through_the_merger() {
        let cohort = simulate(&clean_spec(42));
        let loci = collect(&cohort);

        // Two loci emitted; locus 0 has both samples, locus 1 only sample A.
        assert_eq!(loci.len(), 2);
        let (_, l0) = &loci[0];
        assert_eq!(l0.locus.start, 20);
        assert_eq!(l0.motif, motif(b"CA"));
        assert_eq!(l0.present, vec![0, 1]);

        let (_, l1) = &loci[1];
        assert_eq!(l1.locus.start, 100);
        assert_eq!(l1.present, vec![0]); // B is absent here

        // Every present sample has at least one observed sequence.
        for (_, cl) in &loci {
            for ev in &cl.samples {
                assert!(!ev.seq_counts.is_empty());
            }
        }
    }

    #[test]
    fn modal_observation_matches_the_homozygote_allele() {
        let cohort = simulate(&clean_spec(7));
        let (_, l0) = &collect(&cohort)[0];

        // Sample A (index 0 in `present`) is homozygous 6 units of CA → its most
        // supported observed sequence is the faithful tract (CA × 6).
        let a = &l0.samples[0];
        let (modal_seq, _) = a
            .seq_counts
            .iter()
            .max_by_key(|(_, count)| *count)
            .expect("A has evidence");
        assert_eq!(&**modal_seq, &build_tract(&motif(b"CA"), 6)[..]);
    }

    #[test]
    fn truth_table_is_queryable() {
        let cohort = simulate(&clean_spec(1));
        let truth = &cohort.truth;

        assert_eq!(truth.genotype(0, 0).unwrap().allele_units, vec![6, 6]);
        assert_eq!(truth.genotype(1, 0).unwrap().allele_units, vec![4, 9]);
        assert!(truth.genotype(1, 1).is_none()); // B absent at locus 1
        assert_eq!(truth.chemistry_of(1).level.baseline, 0.10);
        assert_eq!(
            truth.group_of_sample,
            vec![SampleGroupId(0), SampleGroupId(1)]
        );
    }

    #[test]
    fn same_seed_is_byte_identical_and_different_seed_diverges() {
        let a = simulate(&clean_spec(123));
        let b = simulate(&clean_spec(123));
        assert_eq!(a.sample_psp, b.sample_psp, "same seed ⇒ identical bytes");

        let c = simulate(&clean_spec(124));
        assert_ne!(
            a.sample_psp, c.sample_psp,
            "a different seed should perturb the reads"
        );
    }

    #[test]
    fn mononucleotide_locus_simulates_and_round_trips() {
        // The C1 / period-1 case: a single-base motif.
        let chem = SimChemistry {
            error: PerBaseError(0.001),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 1.0,
                decay: 0.15,
            },
            level: StutterLevel {
                baseline: 0.08,
                slope: 0.0,
            },
        };
        let spec = SimCohortSpec {
            seed: 99,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 30,
                motif: motif(b"A"),
                ref_units: 10,
            }],
            groups: vec![chem],
            samples: vec![SimSample {
                name: "S".into(),
                group: SampleGroupId(0),
                genotypes: vec![Some(SimGenotype::homozygous(10, 2))],
            }],
            depth: 60,
        };
        let cohort = simulate(&spec);
        let loci = collect(&cohort);
        assert_eq!(loci.len(), 1);
        assert_eq!(loci[0].1.motif, motif(b"A"));
        assert!(!loci[0].1.samples[0].seq_counts.is_empty());
    }

    #[test]
    fn higher_level_group_slips_more_often() {
        // Same allele, two chemistries differing only in stutter level: the higher
        // level produces a smaller faithful fraction. Guards the level knob.
        let allele = SimGenotype::homozygous(10, 2);
        let motif = motif(b"CA");
        let faithful = build_tract(&motif, 10).into_boxed_slice();

        let faithful_fraction = |level: f64| {
            let chem = SimChemistry {
                error: PerBaseError(0.0),
                shape: StutterShape {
                    up_rate: 1.0,
                    down_rate: 1.0,
                    decay: 0.2,
                },
                level: StutterLevel {
                    baseline: level,
                    slope: 0.0,
                },
            };
            let counts = simulate_cell(555, &chem, &allele, &motif, 4000);
            let total: u32 = counts.iter().map(|(_, c)| c).sum();
            let faithful_count = counts
                .iter()
                .find(|(s, _)| *s == faithful)
                .map(|(_, c)| *c)
                .unwrap_or(0);
            faithful_count as f64 / total as f64
        };

        assert!(
            faithful_fraction(0.05) > faithful_fraction(0.40),
            "a higher stutter level must lower the faithful fraction"
        );
    }
}
