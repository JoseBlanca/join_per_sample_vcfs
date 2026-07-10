//! Rung-ladder & confident-genotype resolution — the shared peak primitive both
//! Stage-2 halves call (arch `ssr_call_parameters.md` §2, spec §5,
//! implementation plan §3.4, B1).
//!
//! Two things live here: [`build_rungs`] pools the cohort's observed sequences into
//! the **length-keyed rung ladder** (the cross-sample coordinate system §5 reads), and
//! [`resolve_confident_genotype`] runs the **model-based** confident-genotype gate on
//! one sample — the 1-vs-2-allele **BIC resolution test** (D1b, spec
//! [`ssr_bic_confident_genotype.md`]): it scores the reads under the best one-allele
//! (homozygote) vs best two-allele (het) model with the C2 read likelihood `Qᵣ` and
//! admits the second allele only when it earns a purity-tuned complexity penalty,
//! guarded by a depth floor, a length-separation clean-seed rule, and cohort
//! recurrence. A *confident genotype* (a homozygote, or a resolved het — including a
//! **same-length** interruption het) is the labelled seed the pre-pass estimates
//! chemistry from.
//!
//! The gate replaced the earlier length-histogram heuristic (which could not see a
//! same-length het and needed separate separation/dosage-balance rules `Qᵣ` now
//! subsumes). Pure-allele lengths are exact multiples of the period; B1 keys rungs by
//! `seq.len() / period` units (substitution variants share a length, the
//! in-tract-no-indel invariant of §6). [`sample_clear_peaks`] still exposes the
//! length-histogram peaks for candidate assembly (C1).

use std::collections::BTreeMap;

use crate::ssr::cohort::likelihood::{LikelihoodScratch, read_given_genotype, read_likelihood};
use crate::ssr::cohort::param_estimation::StutterShape;
use crate::ssr::cohort::types::{CohortLocus, SampleEvidence};
use crate::ssr::types::Motif;

/// A distinct observed sequence at a rung with its cohort support — the element a rung
/// holds.
///
/// `reads` is the total cohort read count; `samples` is the number of **distinct samples**
/// that observed it. The sample count is the recurrence signal the same-length admission
/// bar keys on (spec §5.2): a real interruption haplotype recurs across its carriers, a
/// per-base substitution error is sporadic and lands at a different base each time. It is
/// an integer count → order-free, so it does not disturb cross-thread determinism (arch §4).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RungSeq {
    /// The distinct tract sequence.
    pub(crate) seq: Box<[u8]>,
    /// Total cohort reads supporting this sequence at its length.
    pub(crate) reads: u32,
    /// Number of distinct samples that observed this sequence (spec §5.2 recurrence).
    pub(crate) samples: u32,
}

/// One resolved peak with its labelled parent allele.
///
/// (Shape, not final — B1 may extend it.) Carries what the dosage-consistency and
/// cohort-recurrence checks need: the parent allele's tract sequence, its length in
/// repeat units, and the read support at the peak.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PeakAllele {
    /// The peak's parent allele (tract sequence).
    pub(crate) allele: Box<[u8]>,
    /// The allele's length in **repeat units**.
    pub(crate) repeat_len: u16,
    /// Read support at the peak.
    pub(crate) support: u32,
}

/// A confidently resolved genotype: 1..=ploidy clear peaks, each a labelled parent
/// allele (homozygote = the 1-peak case; a separated het = 2 peaks).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum ResolvedGenotype {
    /// The resolved peaks (length `1..=ploidy`).
    Peaks(Vec<PeakAllele>),
}

/// Why a (sample, locus) failed confident resolution and is left to the soft EM
/// rather than forced into a label (arch §2).
///
/// The model-based gate (D1b) replaces the old length-histogram rejections: the BIC
/// 1-vs-2-allele test decides hom vs het by likelihood, so a homozygote-plus-heavy-
/// stutter simply resolves to a homozygote (the second allele earns no BIC gain — the
/// model subsumes the old dosage-balance heuristic, so there is no `DosageInconsistent`
/// reason any more, spec §3). `NonRecurrent` and `Thin` are unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum UnresolvedReason {
    /// The gate did not produce a clean seed het: the BIC test admitted a second allele
    /// but the two alleles are **length-adjacent** (closer than `separation_min`, yet not
    /// same-length), so their stutter skirts overlap and neither a het seed (skirts
    /// entangled) nor a homozygote seed (would bury a real allele in stutter) is safe —
    /// so the sample is left to the soft EM. Named for the historical merged-het case.
    Merged,
    /// A resolved allele does not recur across enough samples (a length-separated allele
    /// that is not a recurrent peak, or a same-length minority sequence seen too rarely
    /// to be told from a per-read artefact — spec §2.4/§3).
    NonRecurrent,
    /// Too little depth to resolve the peaks at all.
    Thin,
}

/// The outcome of confident-genotype resolution for one (sample, locus).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Resolution {
    /// Resolved at full confidence — a chemistry seed.
    Confident(ResolvedGenotype),
    /// Not resolved — left to the soft EM, with the reason recorded.
    Unresolved(UnresolvedReason),
}

/// Thresholds for rung admission + confident-genotype resolution.
///
/// All are deliberately exposed (no hidden defaults); the numbers are pinned on the
/// simulator in F2. [`dev_default`](Self::dev_default) gives the working values the
/// tests + the burn-in seed start from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct RungCfg {
    /// Prominence floor: a length is a *clear* local maximum only when its support
    /// stands **more than** this many reads above each adjacent (±1 unit) length
    /// (spec §5: "> 3 reads above each adjacent rung").
    pub(crate) prominence: u32,
    /// A resolved allele must recur as a *clear peak* in at least this many samples
    /// (the cohort-recurrence guard against a stutter artifact masquerading as an
    /// allele). Counts the sample itself. The driver **clamps this to the cohort size**
    /// (`min(2, n_samples)`) so a single-sample run degrades gracefully: with no cohort
    /// to corroborate against, the allele need only be a clear peak in that one sample,
    /// and the `min_depth` / `prominence` / het-BIC gates carry the confidence instead.
    /// Multi-sample runs (`n_samples ≥ 2`) are unaffected — it stays 2.
    pub(crate) recurrence_k: u32,
    /// A **length-separated** confident het must have its two alleles at least this
    /// many repeat units apart to be a clean seed: closer than this, the alleles' stutter
    /// skirts overlap and cannot be cleanly labelled (spec §4.3 requires ≥ 2). This does
    /// **not** apply to a **same-length** het — its two alleles are told apart by
    /// composition, not length, so there is no skirt-overlap ambiguity (spec §2.3/§3).
    pub(crate) separation_min: u16,
    /// A sample needs at least this total depth to attempt resolution.
    pub(crate) min_depth: u32,
}

impl RungCfg {
    /// Working defaults (recalibrated in F2).
    pub(crate) fn dev_default() -> Self {
        Self {
            prominence: 3,
            recurrence_k: 2,
            separation_min: 2,
            min_depth: 10,
        }
    }
}

/// The coded **seed** parameters the model-based confident-genotype gate (D1) scores
/// the read likelihood `Qᵣ` with — the "burn-in start" of Mark-2 §4.3/§4.4
/// (spec [`ssr_bic_confident_genotype.md`] §4).
///
/// The gate needs `ε`, a stutter shape, and a stutter level to score reads, but the
/// pre-pass is what *measures* those — a bootstrap. D1 breaks it with these coded
/// seeds; the gate is prior hygiene (the EM can overrule a biased prior), and the case
/// it targets (the same-length het) is robust to a rough seed because its signal is
/// composition, not length. The data-driven co-evolution that replaces this seed with
/// the pre-pass's own fitted params — re-running the gate inside the burn-in loop — is
/// roadmap D2, out of D1 scope. Every field is an exposed dev value (no hidden
/// defaults); the numbers are pinned on the simulator in F2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct GateParams {
    /// Seed per-base error `ε` for `Qᵣ`.
    pub(crate) eps: f64,
    /// Seed stutter shape (direction split + geometric decay of multi-unit slips).
    pub(crate) shape: StutterShape,
    /// Seed stutter level (≈ `P(Δ ≠ 0)`); constant in the seed (no length slope).
    pub(crate) level: f64,
    /// Seed uniform-outlier weight `λ` for the genotype-likelihood floor.
    pub(crate) lambda: f64,
    /// The purity-tuned BIC coefficient: admit the second allele (call a het) iff
    /// `lnL̂₂ − lnL̂₁ > het_admission_cost · ln(n)` (`n` = reads scored). BIC-neutral
    /// for one extra parameter is `½`; the seed runs it **well above** ½ so almost no
    /// false het sneaks into the chemistry seed (a hidden het poisons `ε`/`θ`; a
    /// discarded homozygote only costs a little data — spec §2.1). Pinned in F2.
    pub(crate) het_admission_cost: f64,
    /// Cap on how many of the sample's distinct observed sequences are searched as
    /// candidate alleles (top-M by read support, then bytes) — bounds the pair search
    /// and keeps scattered low-count error variants out of it (spec §2.4).
    pub(crate) max_candidates: usize,
}

impl GateParams {
    /// Working seed values (recalibrated in F2). The shape is symmetric with a mild
    /// decay — a deliberately neutral seed. (Its literals coincide with the pre-pass's
    /// `FALLBACK_SHAPE` today, but the two are *independent* seeds, each pinned
    /// separately in F2, so they are not coupled and may diverge.) `het_admission_cost`
    /// is deliberately conservative (§2.1).
    pub(crate) fn dev_default() -> Self {
        Self {
            eps: 0.005,
            shape: StutterShape {
                up_rate: 0.5,
                down_rate: 0.5,
                decay: 0.1,
            },
            level: 0.05,
            lambda: 0.01,
            het_admission_cost: 3.0,
            max_candidates: 6,
        }
    }
}

/// The pooled cohort rung ladder at one locus: the occupied repeat-unit lengths and,
/// per length, the distinct sequences (with cohort counts) and how many samples
/// peak there.
///
/// "Length" is `seq.len() / period` repeat units. The ladder is the cross-sample
/// coordinate system S2 reads and candidate assembly (C1) unions over; B1 uses the
/// per-length **peak recurrence** for the resolution guard.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Rungs {
    /// The locus motif period (bytes per repeat unit).
    period: usize,
    /// Ascending occupied lengths (repeat units) — the rung coordinate system.
    lengths: Vec<u16>,
    /// Per length: total reads observed across the cohort.
    cohort_support: BTreeMap<u16, u32>,
    /// Per length: number of samples in which that length is a clear local maximum.
    peak_recurrence: BTreeMap<u16, u32>,
    /// Per length: distinct sequences at that length with their cohort support
    /// (a rung holds a *set* of sequences — substitution / interruption variants).
    seqs_by_length: BTreeMap<u16, Vec<RungSeq>>,
}

impl Rungs {
    /// The ascending occupied rung lengths (repeat units).
    pub(crate) fn lengths(&self) -> &[u16] {
        &self.lengths
    }

    /// Total cohort reads observed at `length` units.
    pub(crate) fn cohort_support(&self, length: u16) -> u32 {
        self.cohort_support.get(&length).copied().unwrap_or(0)
    }

    /// How many samples have a clear local maximum at `length` units.
    pub(crate) fn peak_recurrence(&self, length: u16) -> u32 {
        self.peak_recurrence.get(&length).copied().unwrap_or(0)
    }

    /// The distinct sequences (with cohort support) at `length` units, byte-sorted, if any.
    pub(crate) fn seqs_at(&self, length: u16) -> Option<&[RungSeq]> {
        self.seqs_by_length.get(&length).map(Vec::as_slice)
    }

    /// The locus motif period.
    pub(crate) fn period(&self) -> usize {
        self.period
    }

    /// The cohort modal allele length (the most-supported rung; ties → the shorter
    /// length). `None` only for an empty locus. This is the centre of the `G₀` prior
    /// (spec §5.5 — the mode, *not* the reference).
    pub(crate) fn modal_length(&self) -> Option<u16> {
        self.cohort_support
            .iter()
            .max_by(|(la, sa), (lb, sb)| sa.cmp(sb).then_with(|| lb.cmp(la)))
            .map(|(length, _)| *length)
    }
}

/// One sample's length histogram (repeat units → total supporting reads).
fn sample_histogram(evidence: &SampleEvidence, period: usize) -> BTreeMap<u16, u32> {
    let mut histogram = BTreeMap::new();
    for (seq, count) in &evidence.seq_counts {
        let units = (seq.len() / period) as u16;
        *histogram.entry(units).or_insert(0) += count;
    }
    histogram
}

/// Whether `length` is a *clear* local maximum in `histogram`: its support exceeds
/// each adjacent (±1 unit) length's support by more than `prominence` reads (a
/// missing neighbour counts as zero support).
fn is_clear_peak(histogram: &BTreeMap<u16, u32>, length: u16, prominence: u32) -> bool {
    let support = histogram.get(&length).copied().unwrap_or(0);
    if support == 0 {
        return false;
    }
    let lower = length
        .checked_sub(1)
        .and_then(|l| histogram.get(&l).copied())
        .unwrap_or(0);
    let upper = length
        .checked_add(1)
        .and_then(|l| histogram.get(&l).copied())
        .unwrap_or(0);
    support > lower + prominence && support > upper + prominence
}

/// Pool the cohort's observed sequences into the rung ladder (spec §5 levels 1–2).
pub(crate) fn build_rungs(locus: &CohortLocus, cfg: &RungCfg) -> Rungs {
    let period = locus.motif.period().max(1);
    let mut cohort_support: BTreeMap<u16, u32> = BTreeMap::new();
    let mut peak_recurrence: BTreeMap<u16, u32> = BTreeMap::new();
    let mut seqs_by_length: BTreeMap<u16, Vec<RungSeq>> = BTreeMap::new();

    for evidence in &locus.samples {
        // Each `(seq, count)` in a sample's `seq_counts` is one distinct sample observing
        // that sequence (the Stage-1 contract keys `seq_counts` on distinct sequences), so
        // every match/insert here increments the per-sequence distinct-sample tally by one.
        for (seq, count) in &evidence.seq_counts {
            let units = (seq.len() / period) as u16;
            *cohort_support.entry(units).or_insert(0) += count;
            let bucket = seqs_by_length.entry(units).or_default();
            match bucket.iter_mut().find(|rs| rs.seq.as_ref() == seq.as_ref()) {
                Some(rs) => {
                    rs.reads += count;
                    rs.samples += 1;
                }
                None => bucket.push(RungSeq {
                    seq: seq.clone(),
                    reads: *count,
                    samples: 1,
                }),
            }
        }
        let histogram = sample_histogram(evidence, period);
        for &length in histogram.keys() {
            if is_clear_peak(&histogram, length, cfg.prominence) {
                *peak_recurrence.entry(length).or_insert(0) += 1;
            }
        }
    }

    // Deterministic, distinct-sequence order within each rung (bytes ascending).
    for bucket in seqs_by_length.values_mut() {
        bucket.sort_by(|a, b| a.seq.cmp(&b.seq));
    }

    let lengths = cohort_support.keys().copied().collect();
    Rungs {
        period,
        lengths,
        cohort_support,
        peak_recurrence,
        seqs_by_length,
    }
}

/// One sample's clear local maxima as labelled peaks (shared by the confident-
/// genotype gate and candidate assembly, C1). Ascending by length.
pub(crate) fn sample_clear_peaks(
    evidence: &SampleEvidence,
    period: usize,
    prominence: u32,
) -> Vec<PeakAllele> {
    let histogram = sample_histogram(evidence, period);
    histogram
        .keys()
        .copied()
        .filter(|&length| is_clear_peak(&histogram, length, prominence))
        .map(|length| PeakAllele {
            allele: representative_sequence(evidence, length, period),
            repeat_len: length,
            support: histogram[&length],
        })
        .collect()
}

/// The sample's most-supported sequence at `length` units (its representative
/// parent allele for a resolved peak). `length` is known to be occupied.
fn representative_sequence(evidence: &SampleEvidence, length: u16, period: usize) -> Box<[u8]> {
    evidence
        .seq_counts
        .iter()
        .filter(|(seq, _)| (seq.len() / period) as u16 == length)
        // Pick the most-supported sequence; break count ties on the lexicographically
        // smallest sequence (so `max_by` returns it, the tie-break order is reversed).
        .max_by(|(a_seq, a_c), (b_seq, b_c)| a_c.cmp(b_c).then_with(|| b_seq.cmp(a_seq)))
        .map(|(seq, _)| seq.clone())
        .expect("a resolved peak length is occupied by the sample")
}

/// One candidate allele the gate scores: a distinct sequence the sample observed, with
/// its repeat-unit length and read support. Borrows the sample's bytes (no copy until a
/// candidate is chosen as a peak).
struct GateCandidate<'a> {
    seq: &'a [u8],
    /// The candidate's length in **repeat units** (the same quantity `PeakAllele` and the
    /// `Rungs` API call `repeat_len` / `length`).
    repeat_len: u16,
    support: u32,
}

/// Run the **model-based** confident-genotype gate on one sample: the 1-vs-2-allele
/// BIC resolution test (Mark-2 §4.3 Q-P7, spec [`ssr_bic_confident_genotype.md`] §2).
///
/// Scores the sample's reads under the best one-allele (homozygote) model vs the best
/// two-allele (het) model with the read likelihood `Qᵣ` (`seed` supplies the coded
/// scoring params; `scratch` is a reusable `Qᵣ` buffer), and admits the second allele
/// only when its log-likelihood gain beats the purity-tuned BIC penalty
/// `het_admission_cost · ln(n)`. The candidate alleles are the sample's own distinct
/// observed sequences (so a same-length interruption het can be posed — the two alleles
/// share a length but differ in composition). Guards kept: a minimum-depth skip; a
/// **length-separation** clean-seed guard (a length-adjacent het is not seedable, but a
/// *same-length* het is — §2.3); and cohort recurrence, sequence-keyed for same-length
/// alleles (§2.4). The old dosage-balance heuristic is subsumed by `Qᵣ` predicting the
/// stutter skirt.
///
/// Determinism: the whole decision is computed within this one locus (never across a
/// thread boundary), summing read log-likelihoods in `seq_counts` (byte-sorted) order,
/// with every argmax tie broken on the smaller candidate sequence / lower index — so it
/// is byte-identical across thread counts with plain fixed-order f64 (arch §4).
pub(crate) fn resolve_confident_genotype(
    sample: &SampleEvidence,
    rungs: &Rungs,
    ploidy: u8,
    cfg: &RungCfg,
    motif: &Motif,
    seed: &GateParams,
    scratch: &mut LikelihoodScratch,
) -> Resolution {
    // The gate's coded seed must keep the outlier weight positive (it floors
    // `P(read | genotype) > 0`, so `genotype_ln_likelihood` never takes `ln(0)`) and
    // search at least one candidate.
    debug_assert!(seed.lambda > 0.0, "the outlier weight must floor P > 0");
    debug_assert!(
        seed.max_candidates >= 1,
        "must search at least one candidate"
    );

    let depth: u32 = sample.seq_counts.iter().map(|(_, c)| c).sum();
    if depth < cfg.min_depth {
        return Resolution::Unresolved(UnresolvedReason::Thin);
    }

    let candidates = top_candidates(sample, rungs.period, seed.max_candidates);
    if candidates.is_empty() {
        // Adequate depth ⇒ a non-empty `seq_counts`, so this is reached only under a
        // misconfigured `max_candidates == 0` (guarded above in debug). Fall back safely.
        return Resolution::Unresolved(UnresolvedReason::Thin);
    }

    // Qᵣ(read | candidate) once per (read, candidate); the matrix is shared by both the
    // one-allele and two-allele searches. `distinct_seq_count` is the sample's distinct
    // sequence count — the local proxy for the outlier floor `λ/D` (it enters both
    // hypotheses identically, so its exact value does not steer the comparison).
    let qr = gate_score_matrix(sample, &candidates, motif, seed, scratch);
    let counts: Vec<u32> = sample.seq_counts.iter().map(|(_, c)| *c).collect();
    let distinct_seq_count = sample.seq_counts.len().max(1);

    // Best one-allele (homozygote) model.
    let (best_a, ln_l1) = best_one_allele(
        &qr,
        &counts,
        candidates.len(),
        seed.lambda,
        distinct_seq_count,
    );
    let mut chosen_indices = vec![best_a];

    // Best two-allele (het) model, diploid in v1. Admit the second allele only when its
    // gain earns the purity-tuned penalty, and only when the pair is a clean seed.
    if ploidy >= 2 && candidates.len() >= 2 {
        let (a, b, ln_l2) = best_two_alleles(
            &qr,
            &counts,
            candidates.len(),
            seed.lambda,
            distinct_seq_count,
        );
        let earns_penalty = ln_l2 - ln_l1 > seed.het_admission_cost * f64::from(depth).ln();
        if earns_penalty {
            let (len_a, len_b) = (candidates[a].repeat_len, candidates[b].repeat_len);
            let length_gap = i32::from(len_a).abs_diff(i32::from(len_b));
            // Same-length het (gap 0): clean — the two alleles are told apart by
            // composition. Length-separated het: clean only when ≥ separation_min apart
            // (skirts do not overlap). A length-adjacent pair is not a seedable het.
            if len_a == len_b || length_gap >= u32::from(cfg.separation_min) {
                chosen_indices = vec![a, b];
            } else {
                return Resolution::Unresolved(UnresolvedReason::Merged);
            }
        }
    }

    // Build the labelled peaks (a copy of the chosen candidate bytes).
    let mut peaks: Vec<PeakAllele> = chosen_indices
        .iter()
        .map(|&i| PeakAllele {
            allele: Box::from(candidates[i].seq),
            repeat_len: candidates[i].repeat_len,
            support: candidates[i].support,
        })
        .collect();

    // Cohort recurrence (spec §2.4/§3): every chosen allele must recur across ≥ k
    // samples. Sequence-keyed for a same-length pair (a private per-read substitution
    // artefact is rejected here); length-keyed otherwise. A non-recurrent admitted
    // allele fails the whole genotype — we neither seed it as a het (unsure it is real)
    // nor collapse to a homozygote (that would re-inflate ε with its reads).
    if !all_alleles_recurrent(&peaks, rungs, cfg.recurrence_k) {
        return Resolution::Unresolved(UnresolvedReason::NonRecurrent);
    }

    // Order for a stable genotype: by length, then composition (a same-length pair).
    peaks.sort_by(|x, y| {
        x.repeat_len
            .cmp(&y.repeat_len)
            .then_with(|| x.allele.cmp(&y.allele))
    });
    Resolution::Confident(ResolvedGenotype::Peaks(peaks))
}

/// The sample's distinct observed sequences as gate candidates, capped to the top
/// `max` by read support (ties on the smaller bytes), then returned in **byte order**
/// (the canonical order the argmax tie-breaks rely on). Keeping only the best-supported
/// sequences bounds the pair search and keeps scattered low-count error variants out of
/// it (spec §2.4).
fn top_candidates(sample: &SampleEvidence, period: usize, max: usize) -> Vec<GateCandidate<'_>> {
    let mut cands: Vec<GateCandidate> = sample
        .seq_counts
        .iter()
        .map(|(seq, count)| GateCandidate {
            seq,
            repeat_len: (seq.len() / period) as u16,
            support: *count,
        })
        .collect();
    cands.sort_by(|a, b| b.support.cmp(&a.support).then_with(|| a.seq.cmp(b.seq)));
    cands.truncate(max);
    cands.sort_by(|a, b| a.seq.cmp(b.seq));
    cands
}

/// `Qᵣ(read | candidate)` for every (read, candidate) pair, row-major over the sample's
/// reads (`seq_counts` order) and column-major over the byte-ordered candidates.
fn gate_score_matrix(
    sample: &SampleEvidence,
    candidates: &[GateCandidate<'_>],
    motif: &Motif,
    seed: &GateParams,
    scratch: &mut LikelihoodScratch,
) -> Vec<Vec<f64>> {
    sample
        .seq_counts
        .iter()
        .map(|(obs, _)| {
            candidates
                .iter()
                .map(|cand| {
                    read_likelihood(
                        obs,
                        cand.seq,
                        motif,
                        &seed.shape,
                        seed.level,
                        seed.eps,
                        scratch,
                    )
                })
                .collect()
        })
        .collect()
}

/// Total read log-likelihood of a genotype: `Σ_read count · ln P(read | genotype)`,
/// summed in fixed `seq_counts` order (deterministic within the locus). `genotype` is a
/// multiset of candidate indices (a homozygote repeats one; a het holds two).
fn genotype_ln_likelihood(
    qr: &[Vec<f64>],
    counts: &[u32],
    genotype: &[usize],
    lambda: f64,
    distinct_seq_count: usize,
) -> f64 {
    qr.iter()
        .zip(counts)
        .map(|(row, &count)| {
            f64::from(count) * read_given_genotype(row, genotype, lambda, distinct_seq_count).ln()
        })
        .sum()
}

/// The best single allele (`(A, A)` homozygote) and its total log-likelihood. Ties break
/// on the lower candidate index (= smaller bytes, since candidates are byte-ordered) —
/// `reduce` keeps the earlier element, so it must not become `max_by` (which keeps the
/// last), or the determinism contract breaks.
fn best_one_allele(
    qr: &[Vec<f64>],
    counts: &[u32],
    n_candidates: usize,
    lambda: f64,
    distinct_seq_count: usize,
) -> (usize, f64) {
    (0..n_candidates)
        .map(|a| {
            (
                a,
                genotype_ln_likelihood(qr, counts, &[a, a], lambda, distinct_seq_count),
            )
        })
        .reduce(|best, cur| if cur.1 > best.1 { cur } else { best })
        .expect("at least one candidate")
}

/// The best allele pair (`(A, B)` het) and its total log-likelihood. Requires
/// `n_candidates ≥ 2`. Ties break on the lexicographically smaller `(a, b)` index pair
/// (= smaller composition, since candidates are byte-ordered) — same `reduce`-keeps-the-
/// earlier idiom as [`best_one_allele`], over the ascending `(a, b)` pair enumeration.
fn best_two_alleles(
    qr: &[Vec<f64>],
    counts: &[u32],
    n_candidates: usize,
    lambda: f64,
    distinct_seq_count: usize,
) -> (usize, usize, f64) {
    (0..n_candidates)
        .flat_map(|a| ((a + 1)..n_candidates).map(move |b| (a, b)))
        .map(|(a, b)| {
            (
                a,
                b,
                genotype_ln_likelihood(qr, counts, &[a, b], lambda, distinct_seq_count),
            )
        })
        .reduce(|best, cur| if cur.2 > best.2 { cur } else { best })
        .expect("n_candidates >= 2 guaranteed by the caller")
}

/// Whether every chosen allele recurs across ≥ `k` samples (spec §2.4/§3). A peak that
/// **shares its length** with another chosen peak (a same-length het) is checked by its
/// distinct sequence's per-sample tally (`RungSeq.samples`) — the sporadic-vs-systematic
/// recurrence test; any other peak is checked by its length's clear-peak recurrence.
fn all_alleles_recurrent(peaks: &[PeakAllele], rungs: &Rungs, k: u32) -> bool {
    peaks.iter().all(|peak| {
        let shares_length = peaks
            .iter()
            .filter(|other| other.repeat_len == peak.repeat_len)
            .count()
            > 1;
        let recurrence = if shares_length {
            sequence_recurrence(rungs, peak.repeat_len, &peak.allele)
        } else {
            rungs.peak_recurrence(peak.repeat_len)
        };
        recurrence >= k
    })
}

/// How many distinct samples observed `seq` at `length` units (`0` if absent) — the
/// per-sequence recurrence signal a same-length allele is checked on.
fn sequence_recurrence(rungs: &Rungs, length: u16, seq: &[u8]) -> u32 {
    rungs
        .seqs_at(length)
        .and_then(|seqs| seqs.iter().find(|rs| rs.seq.as_ref() == seq))
        .map_or(0, |rs| rs.samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peak(allele: &[u8], repeat_len: u16, support: u32) -> PeakAllele {
        PeakAllele {
            allele: Box::from(allele),
            repeat_len,
            support,
        }
    }

    #[test]
    fn confident_homozygote_is_one_peak() {
        let res = Resolution::Confident(ResolvedGenotype::Peaks(vec![peak(b"ATATAT", 3, 40)]));
        match res {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 1);
                assert_eq!(peaks[0].repeat_len, 3);
            }
            other => panic!("expected a confident 1-peak genotype, got {other:?}"),
        }
    }

    #[test]
    fn confident_separated_het_carries_two_labelled_peaks() {
        let res = Resolution::Confident(ResolvedGenotype::Peaks(vec![
            peak(b"ATAT", 2, 22),
            peak(b"ATATATATAT", 5, 19),
        ]));
        match res {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 2);
                assert_eq!(peaks[0].repeat_len, 2);
                assert_eq!(peaks[1].repeat_len, 5);
            }
            other => panic!("expected a confident 2-peak genotype, got {other:?}"),
        }
    }

    #[test]
    fn unresolved_reasons_are_distinct() {
        assert_ne!(
            Resolution::Unresolved(UnresolvedReason::Merged),
            Resolution::Unresolved(UnresolvedReason::Thin)
        );
    }

    // ── B1: build_rungs + resolve_confident_genotype ──

    use crate::ssr::cohort::types::{CohortLocus, LocusId, SampleEvidence, SsrQc};
    use crate::ssr::types::Motif;

    /// A `CA`-tiled tract of `units` repeat units (12 bp for `units = 6`).
    fn ca_seq(units: u16) -> Box<[u8]> {
        std::iter::repeat_n(*b"CA", units as usize)
            .flatten()
            .collect()
    }

    /// A `SampleEvidence` over a CA (period-2) motif from `(units, count)` bins —
    /// each bin becomes the tiled sequence `CA × units`. Sorted ascending by bytes
    /// (the Stage-1 contract).
    fn ca_evidence(bins: &[(u16, u32)]) -> SampleEvidence {
        let mut seq_counts: Vec<(Box<[u8]>, u32)> = bins
            .iter()
            .map(|&(units, count)| {
                let mut seq = Vec::new();
                for _ in 0..units {
                    seq.extend_from_slice(b"CA");
                }
                (seq.into_boxed_slice(), count)
            })
            .collect();
        seq_counts.sort_by(|(a, _), (b, _)| a.cmp(b));
        SampleEvidence {
            seq_counts,
            qc: SsrQc::default(),
        }
    }

    fn ca_cohort(samples: Vec<SampleEvidence>) -> CohortLocus {
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 20,
                end: 32,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGGCACACATTTTTT".as_slice()),
            Box::from(b"CACACA".as_slice()),
        );
        for (idx, evidence) in samples.into_iter().enumerate() {
            locus.push(idx as u32, evidence);
        }
        locus
    }

    fn resolve_sample0(locus: &CohortLocus, cfg: &RungCfg) -> Resolution {
        let rungs = build_rungs(locus, cfg);
        let mut scratch = LikelihoodScratch::new();
        resolve_confident_genotype(
            &locus.samples[0],
            &rungs,
            2,
            cfg,
            &locus.motif,
            &GateParams::dev_default(),
            &mut scratch,
        )
    }

    #[test]
    fn build_rungs_pools_lengths_support_and_peak_recurrence() {
        // Two homozygous-6 samples (faithful peak at 6, ±1 stutter skirt).
        let sample = ca_evidence(&[(5, 5), (6, 50), (7, 5)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());

        assert_eq!(rungs.lengths(), &[5, 6, 7]);
        assert_eq!(rungs.cohort_support(6), 100); // 50 + 50
        assert_eq!(rungs.peak_recurrence(6), 2); // both samples peak at 6
        assert_eq!(rungs.peak_recurrence(5), 0); // stutter band, never a peak
        let modal = rungs.seqs_at(6).unwrap();
        assert_eq!(modal.len(), 1);
        assert_eq!(modal[0].reads, 100); // 50 + 50 cohort reads
        assert_eq!(modal[0].samples, 2); // both samples observed the length-6 sequence
    }

    #[test]
    fn build_rungs_tallies_distinct_samples_per_same_length_sequence() {
        // Two same-length (6-unit / 12 bp) sequences: a pure `CA×6` and an interrupted
        // variant of equal length. The pure allele recurs in 3 samples, the interrupted in
        // 2, and a singleton substitution noise variant in 1 — the per-sequence sample
        // tally must separate them (the §5.2 recurrence signal).
        let pure = ca_seq(6); // CA×6, 12 bp
        let interrupted: Box<[u8]> = {
            let mut s = pure.to_vec();
            s[5] = b'T'; // flip one interior base → same length, distinct sequence
            s.into_boxed_slice()
        };
        let noise: Box<[u8]> = {
            let mut s = pure.to_vec();
            s[7] = b'G'; // a different sporadic substitution
            s.into_boxed_slice()
        };
        let sample_of = |seqs: &[&[u8]]| {
            let mut seq_counts: Vec<(Box<[u8]>, u32)> =
                seqs.iter().map(|s| (Box::<[u8]>::from(*s), 20)).collect();
            seq_counts.sort_by(|(a, _), (b, _)| a.cmp(b));
            SampleEvidence {
                seq_counts,
                qc: SsrQc::default(),
            }
        };
        let samples = vec![
            sample_of(&[pure.as_ref(), interrupted.as_ref()]),
            sample_of(&[pure.as_ref(), interrupted.as_ref()]),
            sample_of(&[pure.as_ref(), noise.as_ref()]),
        ];
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 20,
                end: 32,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGGCACACACACACATTTTTT".as_slice()),
            pure.clone(),
        );
        for (idx, ev) in samples.into_iter().enumerate() {
            locus.push(idx as u32, ev);
        }
        let rungs = build_rungs(&locus, &RungCfg::dev_default());

        let at6 = rungs.seqs_at(6).unwrap();
        let tally = |seq: &[u8]| at6.iter().find(|rs| rs.seq.as_ref() == seq).unwrap();
        assert_eq!(tally(&pure).samples, 3);
        assert_eq!(tally(&interrupted).samples, 2);
        assert_eq!(tally(&noise).samples, 1);
        // Bytes-ascending order within the rung (determinism contract).
        assert!(at6.windows(2).all(|w| w[0].seq <= w[1].seq));
    }

    #[test]
    fn resolve_confident_genotype_calls_a_clean_homozygote() {
        let sample = ca_evidence(&[(4, 1), (5, 5), (6, 50), (7, 5), (8, 1)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        match resolve_sample0(&locus, &RungCfg::dev_default()) {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 1);
                assert_eq!(peaks[0].repeat_len, 6);
                assert_eq!(&*peaks[0].allele, b"CACACACACACA"); // CA × 6
            }
            other => panic!("expected a confident homozygote, got {other:?}"),
        }
    }

    #[test]
    fn resolve_confident_genotype_calls_a_separated_het_with_two_labelled_peaks() {
        let sample = ca_evidence(&[(3, 1), (4, 30), (5, 4), (8, 4), (9, 28), (10, 1)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        match resolve_sample0(&locus, &RungCfg::dev_default()) {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 2);
                assert_eq!(peaks[0].repeat_len, 4);
                assert_eq!(peaks[1].repeat_len, 9);
            }
            other => panic!("expected a confident separated het, got {other:?}"),
        }
    }

    #[test]
    fn resolve_confident_genotype_rejects_a_one_apart_merged_het() {
        // Alleles 5 and 6 units, balanced (30:28): the BIC test admits a two-allele
        // model, but the pair is length-adjacent (1 apart, not same-length), so the
        // clean-seed guard rejects it — the skirts overlap and it cannot be seeded.
        let sample = ca_evidence(&[(4, 5), (5, 30), (6, 28), (7, 5)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        assert_eq!(
            resolve_sample0(&locus, &RungCfg::dev_default()),
            Resolution::Unresolved(UnresolvedReason::Merged)
        );
    }

    #[test]
    fn resolve_confident_genotype_subsumes_homozygote_plus_heavy_stutter_as_a_hom() {
        // A tall peak at 8 with a −2 stutter satellite at 6: `Qᵣ` predicts the satellite
        // as 8's stutter, so the second allele earns no BIC gain and the sample resolves
        // to a homozygote — no dosage-balance heuristic needed (the model subsumes it).
        let sample = ca_evidence(&[(6, 6), (7, 3), (8, 50), (9, 3)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        match resolve_sample0(&locus, &RungCfg::dev_default()) {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 1);
                assert_eq!(peaks[0].repeat_len, 8);
            }
            other => panic!("expected a confident homozygote at 8, got {other:?}"),
        }
    }

    #[test]
    fn resolve_confident_genotype_rejects_a_non_recurrent_allele() {
        // A clean homozygote, but the allele peaks in only this one sample, so it
        // fails the recurrence_k = 2 guard.
        let sample = ca_evidence(&[(5, 5), (6, 50), (7, 5)]);
        let locus = ca_cohort(vec![sample]);
        assert_eq!(
            resolve_sample0(&locus, &RungCfg::dev_default()),
            Resolution::Unresolved(UnresolvedReason::NonRecurrent)
        );
    }

    #[test]
    fn resolve_confident_genotype_rejects_a_thin_sample() {
        let sample = ca_evidence(&[(6, 3)]); // depth 3 < min_depth 10
        let locus = ca_cohort(vec![sample.clone(), sample]);
        assert_eq!(
            resolve_sample0(&locus, &RungCfg::dev_default()),
            Resolution::Unresolved(UnresolvedReason::Thin)
        );
    }

    /// `CA × units` with the interior base at `flip` changed to `T` — a same-length,
    /// composition-distinct interruption variant of the pure tract.
    fn interrupted_ca(units: u16, flip: usize) -> Box<[u8]> {
        let mut seq = ca_seq(units).to_vec();
        seq[flip] = b'T';
        seq.into_boxed_slice()
    }

    /// A cohort of `n` identical samples, each with the given `(sequence, count)` reads.
    fn seq_cohort(reads: &[(Box<[u8]>, u32)], n: usize) -> CohortLocus {
        let sample = || {
            let mut seq_counts: Vec<(Box<[u8]>, u32)> = reads.to_vec();
            seq_counts.sort_by(|(a, _), (b, _)| a.cmp(b));
            SampleEvidence {
                seq_counts,
                qc: SsrQc::default(),
            }
        };
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 20,
                end: 38,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGGCACACACACACACACACATTTTTT".as_slice()),
            ca_seq(9),
        );
        for idx in 0..n {
            locus.push(idx as u32, sample());
        }
        locus
    }

    #[test]
    fn resolve_confident_genotype_calls_a_same_length_interruption_het() {
        // The capability D1 exists for: a pure-9 / interrupted-9 het (both 18 bp, one
        // interior base apart). The length histogram sees one peak at 9 and would call a
        // homozygote; the BIC test, scoring composition via Qᵣ, admits both alleles at
        // the SAME length. Recurrence is satisfied — both sequences recur in all samples.
        let pure = ca_seq(9);
        let interrupted = interrupted_ca(9, 7);
        let reads = [(pure.clone(), 30), (interrupted.clone(), 30)];
        let locus = seq_cohort(&reads, 3);
        match resolve_sample0(&locus, &RungCfg::dev_default()) {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 2, "a same-length het carries two alleles");
                assert_eq!(peaks[0].repeat_len, peaks[1].repeat_len, "same length");
                assert_eq!(peaks[0].repeat_len, 9);
                let mut alleles: Vec<&[u8]> = peaks.iter().map(|p| p.allele.as_ref()).collect();
                alleles.sort_unstable();
                let mut expected = [pure.as_ref(), interrupted.as_ref()];
                expected.sort_unstable();
                assert_eq!(alleles, expected, "the two distinct compositions, by bytes");
            }
            other => panic!("expected a confident same-length het, got {other:?}"),
        }
    }

    #[test]
    fn resolve_confident_genotype_does_not_invent_an_allele_from_an_error_halo() {
        // The mirror bias (spec §2.4): a homozygote whose non-faithful reads are a FLAT
        // SPREAD of low-count single-base error variants (not one concentrated allele).
        // The BIC test must NOT admit a second allele — the diploid mixing penalty on the
        // 90 faithful reads dwarfs the gain from a 1–2-read error variant — so `ε` is not
        // deflated by reclassifying error as allele. Resolves to the homozygote.
        let pure = ca_seq(9);
        let reads = [
            (pure.clone(), 90),
            (interrupted_ca(9, 3), 2),
            (interrupted_ca(9, 5), 2),
            (interrupted_ca(9, 7), 1),
            (interrupted_ca(9, 9), 1),
            (interrupted_ca(9, 11), 1),
        ];
        let locus = seq_cohort(&reads, 2);
        match resolve_sample0(&locus, &RungCfg::dev_default()) {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(
                    peaks.len(),
                    1,
                    "a dispersed error halo is not a second allele"
                );
                assert_eq!(&*peaks[0].allele, pure.as_ref());
            }
            other => panic!("expected a confident homozygote, got {other:?}"),
        }
    }

    /// A `SampleEvidence` from arbitrary `(sequence, count)` reads, byte-sorted (the
    /// Stage-1 contract) — for cohorts whose samples differ in composition.
    fn evidence_of(reads: &[(Box<[u8]>, u32)]) -> SampleEvidence {
        let mut seq_counts = reads.to_vec();
        seq_counts.sort_by(|(a, _), (b, _)| a.cmp(b));
        SampleEvidence {
            seq_counts,
            qc: SsrQc::default(),
        }
    }

    #[test]
    fn genotype_ln_likelihood_is_finite_for_a_no_match_read() {
        // A read matching no candidate (an all-zero Qᵣ row) stays finite via the λ/D
        // outlier floor — never ln(0) = −∞ (the invariant the gate's `lambda > 0` seed
        // guarantees; without it the argmax would be poisoned by −∞/NaN).
        let qr = vec![vec![0.0, 0.0]];
        let counts = [10u32];
        let ln = genotype_ln_likelihood(&qr, &counts, &[0, 1], 0.01, 4);
        assert!(ln.is_finite(), "outlier floor keeps ln finite, got {ln}");
    }

    #[test]
    fn resolve_confident_genotype_rejects_a_same_length_het_with_a_private_minority() {
        // Sample 0 is a pure-9 / interrupted-9 het; the other two samples carry only the
        // pure allele. The interrupted sequence therefore recurs in ONE sample
        // (< recurrence_k = 2) and is rejected NonRecurrent — the sequence-keyed guard
        // (spec §2.4) that stops a private per-read artefact from deflating ε.
        let pure = ca_seq(9);
        let interrupted = interrupted_ca(9, 7);
        let het = evidence_of(&[(pure.clone(), 30), (interrupted.clone(), 30)]);
        let pure_only = evidence_of(&[(pure.clone(), 60)]);
        let locus = ca_cohort(vec![het, pure_only.clone(), pure_only]);
        assert_eq!(
            resolve_sample0(&locus, &RungCfg::dev_default()),
            Resolution::Unresolved(UnresolvedReason::NonRecurrent)
        );
    }

    #[test]
    fn resolve_confident_genotype_admits_a_het_exactly_at_separation_min() {
        // Alleles at 4 and 6 units: the length gap equals separation_min (2), the
        // minimum clean-seed separation → admitted (boundary of the `>=` guard).
        let sample = ca_evidence(&[(3, 1), (4, 30), (5, 2), (6, 30), (7, 1)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        match resolve_sample0(&locus, &RungCfg::dev_default()) {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 2);
                assert_eq!((peaks[0].repeat_len, peaks[1].repeat_len), (4, 6));
            }
            other => panic!("gap == separation_min should admit a het, got {other:?}"),
        }
    }

    #[test]
    fn resolve_confident_genotype_resolves_a_single_distinct_sequence_as_a_hom() {
        // One distinct sequence at adequate depth: the two-allele search is skipped
        // (candidates.len() < 2), so `best_two_alleles` is never indexed out of bounds.
        let sample = ca_evidence(&[(6, 50)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        match resolve_sample0(&locus, &RungCfg::dev_default()) {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => assert_eq!(peaks.len(), 1),
            other => panic!("a single distinct sequence is a homozygote, got {other:?}"),
        }
    }

    #[test]
    fn resolve_confident_genotype_never_admits_a_het_at_ploidy_one() {
        // A clean separated-het shape, but ploidy 1 → the two-allele search is skipped
        // (the `ploidy >= 2` guard), so a haploid locus resolves to a single allele.
        let sample = ca_evidence(&[(3, 1), (4, 30), (5, 4), (8, 4), (9, 28), (10, 1)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let mut scratch = LikelihoodScratch::new();
        let res = resolve_confident_genotype(
            &locus.samples[0],
            &rungs,
            1,
            &RungCfg::dev_default(),
            &locus.motif,
            &GateParams::dev_default(),
            &mut scratch,
        );
        match res {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => assert_eq!(peaks.len(), 1),
            other => panic!("a haploid locus resolves to one allele, got {other:?}"),
        }
    }

    #[test]
    fn resolve_confident_genotype_resolves_at_depth_equal_to_min_depth() {
        // depth == min_depth (10) is adequate — the skip guard is strict `<`, not `<=`.
        let sample = ca_evidence(&[(6, 10)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        assert!(matches!(
            resolve_sample0(&locus, &RungCfg::dev_default()),
            Resolution::Confident(_)
        ));
    }
}
