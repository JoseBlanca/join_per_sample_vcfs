# SSR pre-pass D1 — model-based confident-genotype gate (architecture)

**Status:** draft, 2026-07-06, branch `ssr-interruptions`. The *how-the-code-is-wired*
companion to spec
[ssr_bic_confident_genotype.md](../specs/ssr_bic_confident_genotype.md) and to Mark-2
[ssr_cohort_mark2.md](../specs/ssr_cohort_mark2.md) §4.3/§4.4. Where they disagree on
intent, the spec wins; on module/struct/signature layout, this doc wins. This is a
**localized change inside the existing pre-pass** — the confident-genotype gate and
the read-attribution step — not a new stage.

Companions:
- [ssr_call_parameters.md](ssr_call_parameters.md) — the pre-pass this edits.
- [ssr_call_genotyping.md](ssr_call_genotyping.md) — the downstream per-locus pipeline (untouched).
- spec [ssr_bic_confident_genotype.md](../specs/ssr_bic_confident_genotype.md) — intent, model form, guards, determinism.

> **The one invariant.** *The gate's verdict is a pure function of `(the sample's
> reads, the cohort rungs, the seed params)`, computed entirely within one locus. It
> never reads shared mutable state and never crosses a thread boundary, so the pre-pass
> stays byte-identical across `--threads`.*

---

## 1. The two edited call sites (and nothing else)

```text
run_prepass_stats(loci, ploidy, cfg)                     prepass.rs   [UNCHANGED shape]
  └ loci.par_iter().fold(..)                             one locus → one thread
      └ accumulate_locus(stats, locus, rungs, ploidy, cfg, params, scratch)   [EDITED]
          ├ resolve_confident_genotype(sample, rungs, ploidy, cfg, params, scratch)  [REPLACED body]
          │     └ the BIC 1-vs-2-allele test (this doc §2)
          └ per-read attribution: nearest_parent → nearest_called_by_sequence  [SWAPPED, §3]
```

Two functions change. Everything below `accumulate_locus` in the call graph
(`SlipProfile`/`SampleStutterStats` accumulation, `purity_slip`, the `G₀`
allele-spread) keeps its integer-count shape; only *which read lands in which bin*
changes. Everything the EM/candidate-set/`G₀`/burn-in own is untouched (spec §8).

---

## 2. The gate — `resolve_confident_genotype`, re-bodied

### 2.1 Signature change

```rust
// BEFORE
pub(crate) fn resolve_confident_genotype(
    sample: &SampleEvidence, rungs: &Rungs, ploidy: u8, cfg: &RungCfg,
) -> Resolution

// AFTER — add the seed params it scores with and a reusable scratch
pub(crate) fn resolve_confident_genotype(
    sample: &SampleEvidence,
    rungs: &Rungs,
    ploidy: u8,
    cfg: &RungCfg,
    seed: &GateParams,           // NEW: what Qᵣ scores with (spec §4)
    scratch: &mut LikelihoodScratch,   // NEW: reused per locus, allocates once
) -> Resolution
```

`Resolution` / `ResolvedGenotype::Peaks(Vec<PeakAllele>)` are **reused as-is** — the
return shape is unchanged (1 or 2 `PeakAllele`). A same-length het is simply two
`PeakAllele` sharing a `repeat_len` but carrying different `allele` bytes — the struct
already supports it (`allele` is the sequence, `repeat_len` the length). The
downstream reader (`accumulate_locus`) is what must stop keying on `repeat_len` alone
(§3).

### 2.2 `GateParams` — the seed the gate scores with (new, small)

```rust
/// The coded seed parameters the D1 gate scores Qᵣ with (spec §4 — the "burn-in
/// start"). Every field is an exposed dev default, no hidden constants; pinned on the
/// simulator in F2. The data-driven co-evolution that replaces this seed with the
/// pre-pass's own fitted params is roadmap D2 (the burn-in loop), out of D1 scope.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct GateParams {
    /// Seed per-base error ε for Qᵣ.
    pub(crate) eps: f64,
    /// Seed stutter shape (direction split + geometric decay) — the existing FALLBACK_SHAPE.
    pub(crate) shape: StutterShape,
    /// Seed stutter level (≈ P(Δ ≠ 0)); constant in D1 (no length slope in the seed).
    pub(crate) level: f64,
    /// Seed outlier weight λ for the genotype-likelihood floor.
    pub(crate) lambda: f64,
    /// The purity-tuned BIC coefficient: admit the 2nd allele iff
    /// lnL̂₂ − lnL̂₁ > het_admission_cost · ln(n). Default well above ½ (spec §2.1).
    pub(crate) het_admission_cost: f64,
    /// Cap on the number of the sample's distinct sequences searched as candidates
    /// (top-M by read support, then bytes) — bounds the O(M²) pair search.
    pub(crate) max_candidates: usize,
}

impl GateParams {
    pub(crate) fn dev_default() -> Self { /* FALLBACK_SHAPE, eps≈0.005, level≈0.05,
        lambda≈0.01, het_admission_cost≈<conservative, F2>, max_candidates=6 */ }
}
```

`StutterShape`, `FALLBACK_SHAPE`, `LikelihoodScratch`, `read_likelihood`,
`read_given_genotype` all already exist — the gate *composes* them.

### 2.3 The algorithm (deterministic throughout)

```rust
fn resolve_confident_genotype(sample, rungs, ploidy, cfg, seed, scratch) -> Resolution {
    let depth = sample.total_reads();
    if depth < cfg.min_depth { return Unresolved(Thin); }              // pre-BIC skip (kept)

    // Candidates = the sample's OWN distinct observed sequences (byte-sorted already),
    // capped to the top-M by read support then bytes (deterministic).
    let cands: Vec<(&[u8], u16 /*units*/, u32 /*support*/)> = top_m_candidates(sample, period, seed.max_candidates);
    if cands.is_empty() { return Unresolved(Thin); }

    // Precompute Qᵣ(read | cand) once per (read, candidate) — reused by both models.
    // reads = sample.seq_counts in fixed byte order; qr[r][c].
    let qr = score_matrix(&cands, &sample.seq_counts, rungs, period, seed, scratch);

    // M₁: best single allele.  ln L̂₁ = max_c Σ_r count·ln[(1−λ)·qr[r][c] + λ/D]
    let (best1, ln_l1) = argmax_single(&qr, &sample.seq_counts, D, seed.lambda);   // ties → smaller bytes / lower idx

    let mut peaks = vec![ peak_from(cands[best1]) ];

    // M₂ (diploid only in v1): best pair.  ln L̂₂ = max_(a<b) Σ_r count·ln[(1−λ)·½(qr[r][a]+qr[r][b]) + λ/D]
    if ploidy >= 2 && cands.len() >= 2 {
        let (a, b, ln_l2) = argmax_pair(&qr, &sample.seq_counts, D, seed.lambda);   // ties → smaller (bytes,bytes)
        let n = depth as f64;
        if ln_l2 - ln_l1 > seed.het_admission_cost * n.ln() {          // the BIC admission (spec §2.1)
            peaks = vec![ peak_from(cands[a]), peak_from(cands[b]) ];  // may share repeat_len (same-length het!)
        }
    }

    // Cohort recurrence (kept, generalised): every admitted allele recurs across samples.
    if !all_recurrent(&peaks, rungs, cfg.recurrence_k) {              // §2.4
        return Unresolved(NonRecurrent);
    }

    // Order peaks for a stable, dosage-readable genotype: by length, then bytes (same-length pair).
    peaks.sort_by(|x, y| x.repeat_len.cmp(&y.repeat_len).then_with(|| x.allele.cmp(&y.allele)));
    Confident(Peaks(peaks))
}
```

Notes:
- **`peak_from`** fills `PeakAllele { allele, repeat_len, support }` from the candidate
  (`support` = the sample's read count at that sequence; used downstream for the
  allele-copies tally, unchanged semantics).
- **`argmax_single` / `argmax_pair`** reduce Σ over reads in `seq_counts` byte order
  (fixed), break every tie on the smaller sequence bytes then lower candidate index.
  Plain f64 in fixed order — deterministic because it never leaves the locus (spec §5).
- **`M₁ ⊂ M₂`** so `ln_l2 ≥ ln_l1`; the admission is purely the penalty test.

### 2.4 Cohort recurrence, sequence-keyed (`all_recurrent`)

For each admitted peak, require it recurs as a real allele across ≥ `recurrence_k`
samples. Two cases, both already available on `Rungs`:
- **length-separated allele** (its length is not shared by another admitted peak): use
  the existing per-length `rungs.peak_recurrence(len)` (unchanged).
- **same-length allele** (Phase 1's case): use the per-sequence distinct-sample tally
  `RungSeq.samples` from `rungs.seqs_at(len)` — a real interruption allele recurs
  across carriers, a per-read substitution artefact does not. This is the cohort-level
  false-allele defence the single-sample likelihood cannot give (spec §3), and the
  reason `RungSeq.samples` was added in Phase 1.

`UnresolvedReason` shape: keep `Thin`, `NonRecurrent`; `Merged` now means "the BIC test
kept one allele" (the model's verdict); `DosageInconsistent` is **removed** (the model
subsumes it — spec §3). Removing a variant is a small ripple in the enum's `match`
sites and its tests — the arch chooses removal over keeping a dead variant.

---

## 3. The coupled attribution swap in `accumulate_locus`

Resolving a same-length het is inert unless the reads are split between the two
same-length alleles by composition. `accumulate_locus` currently attributes each read
with **length-only** `nearest_parent(read_units, &peak_units)`, which ties on two
same-length peaks and dumps both alleles' reads onto whichever wins the tie —
re-contaminating `ε`.

Swap that one call to the **existing** sequence-aware
`nearest_called_by_sequence(obs, read_units, &called, eps, scratch)`
([attribution.rs](../../../src/ssr/cohort/attribution.rs)), where
`called: &[(&[u8], u16)]` is the admitted peaks' `(allele bytes, repeat_len)`. Its
contract (already documented in `attribution.rs`): length distance first, composition
**only** to break a same-length tie — exactly the same-length het case — and the slip
`Δ` is identical whichever same-length allele wins, so the `θ`/slip statistics are
effect-neutral for length-separated genotypes (no behaviour change there) and *only*
the same-length split is new. The `eps` it uses is `seed.eps` (the same seed the gate
scored with — attribution and gate agree, mirroring the EM's "attribution uses the same
metric as Qᵣ" rule).

Downstream of attribution, the existing per-peak accounting (`peak_faithful`,
`peak_slipped`, `peak_interruptions`, `compare_bases`, `purity_slip`) is **already
keyed per-peak by index** — it was written Phase-1-ready. So once `peak_idx` comes from
`nearest_called_by_sequence` instead of `nearest_parent`, the same-length het's reads
flow into the correct peak's bins with no further change:
- an `interrupted-18` read → `peak_idx` = the interrupted allele → `compare_bases`
  against *its own* sequence → all-match → `base_match` (not `base_mismatch`). **ε
  de-contaminated.**
- its `purity_slip` cell keys on that peak's own `interruption_count`, so the purity
  contrast (P2.2b) also sharpens.

The `BIAS NOTE` at `prepass.rs:223` (het `ε` mildly biased high) narrows: with
same-length hets now resolved and split, the residual is only the genuine inner-valley
soft-split approximation, unchanged.

`accumulate_locus` grows two params (`seed: &GateParams`, `scratch: &mut …`) threaded
from `run_prepass_stats`; the scratch is a per-fold (per-thread) reusable buffer, so no
cross-thread state.

---

## 4. Wiring — where the seed and scratch come from

```text
run_prepass(loci, ploidy, cfg, g0_cfg)                    [+ seed: GateParams]
run_prepass_stats(loci, ploidy, cfg)                      [+ seed: GateParams]
  └ par_iter().fold(|| (PrepassStats::default(), scratch_default()), |(acc, scr), locus| {
        let rungs = build_rungs(locus, cfg);
        accumulate_locus(&mut acc, locus, &rungs, ploidy, cfg, &seed, &mut scr);
        (acc, scr)
    })
    .map(|(acc, _scr)| acc)
    .reduce(PrepassStats::default, |a, b| a.merge(&b))     [reduce UNCHANGED — integer merge]
```

The `fold` accumulator gains a per-thread `LikelihoodScratch` alongside the
`PrepassStats` (allocate once per thread, reused across that thread's loci — the
scratch-buffer idiom). The reduce merges only `PrepassStats` (integer counts), so it
stays order-independent and byte-identical. The `GateParams` seed is a `Copy` value
passed by shared ref — read-only, no contention.

Where `run_prepass` is called (the driver / D2 entry) passes `GateParams::dev_default()`
for now; D2 will pass its burn-in's current fitted params here instead (spec §4). That
call-site is the clean D2 seam.

---

## 5. Cost and scratch

- `Qᵣ` is the expensive primitive. The gate scores `M candidates × R distinct reads`
  per confident-eligible sample, `M ≤ max_candidates` (dev 6). The score matrix is
  computed once and shared by both `M₁` and the `M₂` pair search, so `Qᵣ` is evaluated
  `M·R` times, not `M²·R`. `R` = the sample's distinct sequence count (small — a tract
  histogram, not raw reads). This is comparable to one EM E-step over the locus; the
  pre-pass runs it once per (sample, locus), so it is well within budget.
- `LikelihoodScratch` (the `Qᵣ` placement-variant + DP buffers) is reused across the
  whole matrix and across the thread's loci — one allocation per thread.

---

## 6. Test surface (unit-level; the plan §… owns the sim/benchmark gates)

- `rung_ladder.rs`: a clean homozygote → 1 peak; a length-separated het → 2 peaks (as
  before); a **same-length het** (pure-18 / interrupted-18, both cohort-recurrent) →
  2 same-length peaks with distinct `allele` bytes (**new capability**); a 1-apart
  balanced pair → 1 peak (BIC rejects, "contributes nothing"); a hom+heavy-stutter →
  1 peak (model subsumes dosage); a non-recurrent same-length minority → `NonRecurrent`;
  thin → `Thin`. Determinism: gate verdict identical under `--threads 1` vs `K`.
- `prepass.rs`: on an injected same-length-het simulator with known `ε`, the recovered
  `ε` matches truth (where the pre-swap heuristic over-estimates it); `--threads 1` vs
  `K` byte-identical params; the existing recovery/checkpoint-2 tests still pass
  (length-separated behaviour effect-neutral).

---

## 7. What this doc deliberately does not change

`em.rs`, `candidate_set.rs`, `allele_freq_prior.rs` (`G₀`), `inbreeding.rs`,
`vcf_out.rs`, the reduce in `run_prepass_stats`, the sufficient-statistic structs, and
Stage 1. The gate's hypotheses are the sample's own sequences, not the cohort candidate
set; the params it scores with are a seed, not the EM's output. Any pull on these is a
scope-creep smell — stop and ask (spec §8).
