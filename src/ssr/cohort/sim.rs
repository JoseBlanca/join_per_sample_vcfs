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
    fn index_below(&mut self, n: usize) -> usize {
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

/// Read-level generative noise *on top of* the whole-unit stutter + ε forward model —
/// the knobs that turn the clean in-frame generator (G1) into the bake-off's harder
/// generative cases (plan §5). Each is independent and length-/composition-only, so the
/// **truth genotype (allele lengths in units) is unchanged** — these stress the read
/// likelihood's robustness, not the labels it is scored against.
///
/// - **G1** = [`none`](Self::none) (whole-unit slips only — the original model).
/// - **G2** = [`hipstr_like`](Self::hipstr_like) (adds out-of-frame draws).
/// - **G3** = [`messy`](Self::messy) (out-of-frame + interruptions + a heavier ε tail).
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct GenerativeNoise {
    /// Probability a read's tract length is knocked **off the unit grid** by ±1 bp (an
    /// out-of-frame event), independent of the whole-unit slip.
    pub(crate) out_of_frame_rate: f64,
    /// Probability a read's tract carries one motif-breaking substitution (a sequence
    /// interruption), length-preserving.
    pub(crate) impurity_rate: f64,
    /// Extra per-base substitution added on top of the chemistry's ε (the heavier tail).
    pub(crate) extra_substitution: f64,
}

impl GenerativeNoise {
    /// G1 — the clean whole-unit generator (no extra noise).
    pub(crate) fn none() -> Self {
        Self::default()
    }

    /// G2 — HipSTR-style: whole-unit slips plus a moderate out-of-frame fraction.
    pub(crate) fn hipstr_like() -> Self {
        Self {
            out_of_frame_rate: 0.10,
            impurity_rate: 0.0,
            extra_substitution: 0.0,
        }
    }

    /// G3 — messy realistic: out-of-frame draws, interruptions, and a heavier ε tail.
    /// The case production robustness is weighed on (plan §5/§7).
    pub(crate) fn messy() -> Self {
        Self {
            out_of_frame_rate: 0.12,
            impurity_rate: 0.12,
            extra_substitution: 0.015,
        }
    }
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

/// Tile `motif` to exactly `len_bp` bases (the final unit may be partial — an
/// out-of-frame length). Equals [`build_tract`] when `len_bp` is a unit multiple.
fn build_tract_bp(motif: &Motif, len_bp: usize) -> Vec<u8> {
    let m = motif.as_bytes();
    (0..len_bp).map(|i| m[i % m.len()]).collect()
}

/// Force one motif-breaking substitution at a random tract position (a sequence
/// interruption); length-preserving, so the allele's unit length is unchanged.
fn force_interruption(rng: &mut SplitMix64, tract: &mut [u8]) {
    if tract.is_empty() {
        return;
    }
    let pos = rng.index_below(tract.len());
    let mut replacement = BASES[rng.index_below(BASES.len())];
    while replacement == tract[pos] {
        replacement = BASES[rng.index_below(BASES.len())];
    }
    tract[pos] = replacement;
}

/// Draw one read's realized tract from a parent allele of `parent_units`: a whole-unit
/// slip ([`slip_length`]), then the optional out-of-frame / interruption / heavier-ε
/// noise of `GenerativeNoise`. With [`GenerativeNoise::none`] this is exactly the G1
/// forward model (whole-unit slip + ε substitutions).
fn draw_read(
    rng: &mut SplitMix64,
    chem: &SimChemistry,
    parent_units: u16,
    motif: &Motif,
    noise: &GenerativeNoise,
) -> Box<[u8]> {
    let units = slip_length(rng, chem, parent_units);
    let period = motif.period() as i32;
    let mut len_bp = units as i32 * period;
    if noise.out_of_frame_rate > 0.0 && rng.chance(noise.out_of_frame_rate) {
        len_bp += if rng.chance(0.5) { 1 } else { -1 };
    }
    let len_bp = len_bp.max(1) as usize;

    let mut tract = build_tract_bp(motif, len_bp);
    if noise.impurity_rate > 0.0 && rng.chance(noise.impurity_rate) {
        force_interruption(rng, &mut tract);
    }
    apply_substitutions(
        rng,
        &mut tract,
        (chem.error.0 + noise.extra_substitution).min(1.0),
    );
    tract.into_boxed_slice()
}

/// Apply within-tract substitutions at rate `ε`.
fn apply_substitutions(rng: &mut SplitMix64, tract: &mut [u8], eps: f64) {
    if eps <= 0.0 {
        return;
    }
    for base in tract.iter_mut() {
        if rng.chance(eps) {
            let mut replacement = BASES[rng.index_below(BASES.len())];
            while replacement == *base {
                replacement = BASES[rng.index_below(BASES.len())];
            }
            *base = replacement;
        }
    }
}

/// Generate one (sample, locus) cell's distinct observed sequences under the clean G1
/// forward model (whole-unit slip + ε). Thin wrapper over [`simulate_cell_with`].
fn simulate_cell(
    seed: u64,
    chem: &SimChemistry,
    genotype: &SimGenotype,
    motif: &Motif,
    depth: u32,
) -> Vec<(Box<[u8]>, u32)> {
    simulate_cell_with(seed, chem, genotype, motif, depth, &GenerativeNoise::none())
}

/// Generate one (sample, locus) cell's distinct observed sequences, ascending by bytes
/// (the Stage-1 writer's contract — a `BTreeMap` keyed by the sequence keeps that order
/// for free), under the supplied [`GenerativeNoise`].
fn simulate_cell_with(
    seed: u64,
    chem: &SimChemistry,
    genotype: &SimGenotype,
    motif: &Motif,
    depth: u32,
    noise: &GenerativeNoise,
) -> Vec<(Box<[u8]>, u32)> {
    let mut rng = SplitMix64::new(seed);
    let ploidy = genotype.allele_units.len();
    let mut counts: BTreeMap<Box<[u8]>, u32> = BTreeMap::new();
    for _ in 0..depth {
        let parent = genotype.allele_units[rng.index_below(ploidy)];
        let tract = draw_read(&mut rng, chem, parent, motif, noise);
        *counts.entry(tract).or_insert(0) += 1;
    }
    counts.into_iter().collect()
}

/// Generate a cohort from its spec under the clean G1 forward model. Thin wrapper over
/// [`simulate_with`].
pub(crate) fn simulate(spec: &SimCohortSpec) -> SimCohort {
    simulate_with(spec, &GenerativeNoise::none())
}

/// Generate a cohort from its spec under the supplied [`GenerativeNoise`] (the bake-off
/// generative axis G1/G2/G3). The truth table is the same regardless of noise — the
/// noise perturbs reads, not the genotypes they are scored against.
pub(crate) fn simulate_with(spec: &SimCohortSpec, noise: &GenerativeNoise) -> SimCohort {
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
                    let observed = simulate_cell_with(
                        cell_seed(spec.seed, s_idx, l_idx),
                        chem,
                        genotype,
                        &locus.motif,
                        spec.depth,
                        noise,
                    );
                    let depth: u32 = observed.iter().map(|(_, c)| c).sum();
                    let tract_units_len = (locus.motif.period() * locus.ref_units as usize) as u32;
                    Some(SsrLocusRecord {
                        chrom_id: chrom_id(&locus.chrom),
                        start: locus.start + 1, // container coords are 1-based
                        end: locus.start + 1 + tract_units_len,
                        depth,
                        n_filtered: 0,
                        mapped_reads: depth,
                        n_low_quality: 0,
                        n_border_off_end: 0,
                        n_widened: 0,
                        n_window_truncated: 0,
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

    #[test]
    fn separated_het_deposits_support_at_both_allele_lengths() {
        // A well-separated het (4 vs 9 units, CG-seed): both allele tracts must
        // appear as supported observations — a generator that drew from one copy,
        // or mis-tiled, would show a single mode.
        let chem = SimChemistry {
            error: PerBaseError(0.0),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.05,
                slope: 0.0,
            },
        };
        let motif = motif(b"CA");
        let counts = simulate_cell(2024, &chem, &SimGenotype::diploid(4, 9), &motif, 400);

        let support = |units: u16| {
            let tract = build_tract(&motif, units).into_boxed_slice();
            counts
                .iter()
                .find(|(s, _)| *s == tract)
                .map(|(_, c)| *c)
                .unwrap_or(0)
        };
        assert!(
            support(4) > 50,
            "low allele under-supported: {}",
            support(4)
        );
        assert!(
            support(9) > 50,
            "high allele under-supported: {}",
            support(9)
        );
    }

    #[test]
    fn g2_out_of_frame_noise_produces_non_unit_lengths() {
        // G1 produces only unit-multiple tract lengths; G2 (out-of-frame) must produce
        // some reads whose length is off the period grid — the in/out-of-frame axis.
        let chem = SimChemistry {
            error: PerBaseError(0.0),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 1.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.05,
                slope: 0.0,
            },
        };
        let motif = motif(b"CA");
        let allele = SimGenotype::homozygous(10, 2);
        let period = motif.period();

        let g1 = simulate_cell_with(7, &chem, &allele, &motif, 500, &GenerativeNoise::none());
        assert!(
            g1.iter().all(|(s, _)| s.len() % period == 0),
            "G1 must stay on the unit grid"
        );

        let g2 = simulate_cell_with(
            7,
            &chem,
            &allele,
            &motif,
            500,
            &GenerativeNoise::hipstr_like(),
        );
        assert!(
            g2.iter().any(|(s, _)| s.len() % period != 0),
            "G2 must produce some out-of-frame (non-unit-length) reads"
        );
    }

    #[test]
    fn g3_impurity_breaks_the_motif_tiling() {
        // G3 interruptions must yield reads that are NOT a clean motif tiling even at an
        // in-frame length (a motif-breaking substitution), with ε = 0 so only the
        // interruption knob can cause it.
        let chem = SimChemistry {
            error: PerBaseError(0.0),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 1.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.0, // no slips: isolate the interruption knob
                slope: 0.0,
            },
        };
        let motif = motif(b"CA");
        let allele = SimGenotype::homozygous(8, 2);
        let noise = GenerativeNoise {
            out_of_frame_rate: 0.0,
            impurity_rate: 0.5,
            extra_substitution: 0.0,
        };
        let counts = simulate_cell_with(3, &chem, &allele, &motif, 400, &noise);
        let clean = build_tract(&motif, 8).into_boxed_slice();
        let impure_reads: u32 = counts
            .iter()
            .filter(|(s, _)| **s != *clean && s.len() == clean.len())
            .map(|(_, c)| *c)
            .sum();
        assert!(
            impure_reads > 0,
            "G3 impurity must produce interrupted same-length reads"
        );
    }

    #[test]
    fn differing_group_shape_changes_the_slip_magnitude_distribution() {
        // Two chemistries identical except `shape.decay` (per-group SHAPE — the M3
        // axis): the steeper-tailed shape must produce larger average slip
        // magnitudes. Guards that the forward model actually reads `decay`.
        let allele = SimGenotype::homozygous(12, 2);
        let motif = motif(b"CA");
        let faithful = build_tract(&motif, 12).into_boxed_slice();

        let mean_abs_slip = |decay: f64| {
            let chem = SimChemistry {
                error: PerBaseError(0.0),
                shape: StutterShape {
                    up_rate: 1.0,
                    down_rate: 1.0,
                    decay,
                },
                level: StutterLevel {
                    baseline: 1.0, // every read slips, isolating the magnitude
                    slope: 0.0,
                },
            };
            let counts = simulate_cell(31, &chem, &allele, &motif, 6000);
            let mut slipped_reads = 0u64;
            let mut total_abs_units = 0i64;
            for (seq, count) in &counts {
                if *seq == faithful {
                    continue; // Δ = 0 contributes nothing to mean magnitude
                }
                let units = (seq.len() / motif.period()) as i64;
                total_abs_units += (units - 12).abs() * i64::from(*count);
                slipped_reads += u64::from(*count);
            }
            total_abs_units as f64 / slipped_reads as f64
        };

        assert!(
            mean_abs_slip(0.5) > mean_abs_slip(0.05),
            "a heavier-tailed shape (larger decay) must slip by larger magnitudes"
        );
    }
}
