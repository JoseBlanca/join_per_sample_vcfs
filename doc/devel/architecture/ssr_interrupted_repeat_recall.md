# SSR interrupted-repeat recall — sequence-keyed alleles (architecture)

**Status:** draft, 2026-07-06, branch `ssr-interruptions`. The *how-the-code-is-wired*
companion to spec [ssr_interrupted_repeat_recall.md](../specs/ssr_interrupted_repeat_recall.md)
(Phase 1 structural, Phase 2 statistical) and to Mark-2
[ssr_cohort_mark2.md](../specs/ssr_cohort_mark2.md) §5/§6/§7. Where they disagree on
intent, the spec wins; on module/struct/signature layout, this doc wins. Refines the
genotyping arch sketch [ssr_call_genotyping.md](ssr_call_genotyping.md) — this is a
small, localized change *inside* that pipeline, not a new stage.

Companions:
- [ssr_call_genotyping.md](ssr_call_genotyping.md) — the S1→S3→EM→CALL per-locus pipeline this edits.
- [ssr_call_parameters.md](ssr_call_parameters.md) — the pre-pass (Phase 2 fits here).
- spec [ssr_interrupted_repeat_recall.md](../specs/ssr_interrupted_repeat_recall.md) — intent, risks, validation.

> **The one structural invariant.** *A candidate allele is a **sequence**; a rung
> holds a **set** of them; an allele's identity is its **candidate index**, never its
> tract length.* The EM, genotype enumeration, `G₀`, `is_variable`, `site_qual`, and
> ALT emission already obey this (spec §4). Phase 1 is the finite job of converting the
> **four remaining length-keyed holdouts** — nomination, the seed, the slip attribution,
> and the allele-balance FP term — to key on the candidate/sequence. Nothing else moves.

---

## 1. Where the change sits

Per `CohortLocus`, inside the existing pipeline ([ssr_call_genotyping.md §1](ssr_call_genotyping.md)),
the touched boxes are marked `◆` (Phase 1) and `◇` (Phase 2):

```
 CohortLocus ─► S1 ASSEMBLE ─►  S3 SCORE        ─► per-locus EM      ─► CALL
                ◆ set-valued     (unchanged —       ◆ slip attrib.       ◆ allele-balance
                  nomination +    already scores      by sequence         on deconvolved
                  §5.2 admission   sequences)         ◇ per-allele        per-allele resp.
                ◆ seq-aware seed                       stutter rate       (REPCN unchanged)
```

- **S3 (`Qᵣ`) is untouched** — it already scores each observed sequence against each
  *candidate sequence* (`compute_data_ll` builds `obs_qr[obs][c] = q_r(obs, alleles[c], …)`),
  so two same-length candidates already get different likelihoods. This is the linchpin
  the whole change rests on (spec §4, verified).
- The four `◆` edits are structural (no new statistics). The one `◇` edit (Phase 2) is
  the only new statistics, and it is measurement-gated (spec §6.5/§8).

---

## 2. Phase 1 — module / struct / signature changes

### 2.1 `rung_ladder.rs` — expose the admitted-sequence set

The rung ladder **already stores** the set of distinct sequences per length —
`seqs_by_length: BTreeMap<u16, Vec<SeqCount>>`, sorted by sequence bytes at
`build_rungs` (the `bucket.sort_by(|(a,_),(b,_)| a.cmp(b))` pass). `seqs_at(length)`
already returns the sorted slice. **No storage change.** Add only a thin accessor that
yields the sequences at a length in deterministic (byte-sorted) order for nomination:

```rust
/// The distinct sequences observed at `length`, byte-sorted (for §5.2 admission).
pub(crate) fn seqs_at(&self, length: u16) -> Option<&[SeqCount]>   // exists
```

`cohort_support`, `peak_recurrence`, and `modal_length` are unchanged — `G₀` still keys
off `cohort_support` / `modal_length` (length, not sequence), which is correct: `G₀` is a
length-offset prior and has nothing to say about composition (spec §4).

### 2.2 `candidate_set.rs` — set-valued nomination + the §5.2 admission bar

Replace the single-representative step with a set:

```rust
// was: fn cohort_representative(rungs, length) -> Option<Box<[u8]>>   (most-supported one)
fn cohort_alleles(rungs: &Rungs, length: u16, total_at_len: u32, cfg: &CandidateCfg)
    -> impl Iterator<Item = Box<[u8]>>            // every seq clearing §5.2, byte-sorted
```

Admission predicate (spec §5.2) — **recurrence-based, all three joined** (`&&`), computed
from `seqs_at(length)` plus the per-sequence sample count:

```
promote(seq) ⇔  cohort_reads(seq)     ≥ cfg.min_same_length_reads
             && distinct_samples(seq) ≥ cfg.min_same_length_samples
             && cohort_reads(seq) as f64 ≥ cfg.min_same_length_fraction * total_at_len as f64
```

`SeqCount` today carries only a cohort read count, not a per-sequence **sample** count.
The distinct-sample count is needed for the middle condition, so **`build_rungs` must also
tally, per (length, sequence), the number of distinct samples that observed it** — a
second small accumulator alongside `seqs_by_length` (a `BTreeMap<(u16, Box<[u8]>), u32>`,
or widen `SeqCount` to `(seq, reads, samples)`). This is the one new sufficient statistic
Phase 1 adds; it is an integer count → order-free, determinism-safe.

New `CandidateCfg` fields (defaults from spec §5.2, pinned in the sweep):

```rust
pub(crate) struct CandidateCfg {
    // … existing: prominence, min_cohort_depth, max_candidate_alleles, max_out_of_frame_frac
    pub(crate) min_same_length_reads: u32,      // dev 8
    pub(crate) min_same_length_samples: u32,    // dev 3
    pub(crate) min_same_length_fraction: f64,   // dev 0.10
}
```

Wiring in `assemble_candidates`:
- Reference allele stays candidate 0, seeded unconditionally (unchanged).
- Per nominated length, `push_unique` **every** `cohort_alleles(…)` result (was: the one
  representative). `push_unique` already dedups by exact bytes, so a same-length ALT equal
  to REF is not duplicated.
- `MAX_CANDIDATE_ALLELES` still caps the union → `TooManyAlleles`. The §5.2 bar filters the
  substitution cloud **before** the cap sees it (only recurrent structure survives), so the
  cap counts real alleles. (Regression watch — spec §10 — a hypervariable locus with genuine
  same-length structure at many lengths could cross the cap; the benchmark checks no
  currently-`Pass` locus flips.)
- **Ordering:** candidates within a length enter in byte-sorted order (from `seqs_at`), so
  ALT column order and GT indices are deterministic across threads (spec §5.1).

`is_periodic` / locus admission is unchanged — it is length-only and modal-anchored, so
set-valued rungs don't change which lengths exist (spec §6; existing
`odd_reference_interrupted_repeat_is_admitted` test still holds).

### 2.3 `em_init.rs` — sequence-aware seed

`candidate_of_length` currently returns the *first* candidate at a length. Change it to
match by the sample's representative sequence:

```rust
fn candidate_for(candidates: &CandidateSet, sample: &SampleEvidence, period, length)
    -> Option<usize>   // the candidate whose sequence == the sample's repr. at `length`,
                       // else the closest by align_subst (the same metric §2.4 uses)
```

**Honest limitation (spec §5.3):** a same-length *het* sample shows one length peak with
one representative, so the seed still assigns **both** allele copies to the majority
composition (a hom seed); the minority same-length allele starts at its `G₀` floor.
Recovery is the EM's job — the composition-aware `Qᵣ` (§1) lets the E-step move onto the
het, and the pseudocount floor keeps `π_B > 0`. This is documented, not fixed, in Phase 1;
a seed-level het proposal (converging with the Mark-2 BIC gate) is a follow-up.

### 2.4 `attribution.rs` / `em.rs` — the sequence-aware attribution primitive

Today `nearest_parent(read_units, parent_units)` attributes by length. Add the sequence
primitives (one hard, one soft), both scoring with the read model so *attribution and the
likelihood never disagree* (spec §5.3):

```rust
/// Hard: the called allele a read best matches, + its length slip (for the θ_locus stats).
fn nearest_called_by_sequence(obs, called: &[&[u8]], ctx: &ReadScoringContext, scratch)
    -> (usize, i32);       // argmax q_r (or align_subst on the Δ=0 tie); tie → lower index

/// Soft: a read's responsibility to each called allele (for allele balance, §2.5).
fn allele_responsibilities(obs, called: &[&[u8]], ctx, scratch) -> SmallVec<[f64; 2]>;
                           // q_r(obs|a) normalized over the called alleles
```

- **`attribute_locus`** (the `θ_locus` slip refit) swaps `nearest_parent` for
  `nearest_called_by_sequence`. Integer slip counts are preserved → order-free reduce,
  determinism unchanged. **Effect-neutral for same-length alleles** (both give the same Δ),
  but the code stops being length-keyed and the primitive is shared with §2.5.
- Both reuse `align_subst` (exact match, else the `(ε/3)^mismatch` closed form) — the exact
  metric inside `q_r`, so a read's parent by attribution is its parent by likelihood.

### 2.5 `vcf_out.rs` — sequence-aware allele balance (the issue-2 fix)

`allele_balance` currently (a) short-circuits to `1.0` when the two `genotype_units` are
equal — reading a same-length het as a homozygote — and (b) attributes reads by length.
Both disable the FP defence for exactly the same-length het Phase 1 adds, and Mark-2 §6
(2026-07-06 amendment) mandates the sequence-keyed form. Rewire:

- **Trigger on allele *identity*** — skip only when the two `allele_indices` are equal, not
  when the two `genotype_units` are equal.
- **Deconvolve by composition** — the per-allele support is the soft responsibility sum
  `n_a = Σ_reads count · resp_a` from §2.4's `allele_responsibilities`, not a length bin.
  The binomial-tail balance test then runs on `(n_A, n_B)` vs the ½ expectation, exactly as
  Mark-2 §6 specifies ("deconvolved per-allele responsibilities, not raw lengths").

Where the responsibilities come from is the one new data wire — see §3. `REPCN` emission,
`is_variable`, `site_qual`, GT/ALT formatting are **unchanged** (all already candidate-index
keyed; REPCN stays as the HipSTR-style length annotation, spec §5.4).

### 2.6 What is deliberately **not** touched (spec §4)

`allele_freq_prior.rs` (`g0_pseudocounts` — a `Vec` parallel to candidates, same-length
candidates get equal value, no collision), `em.rs` `enumerate_diploid_genotypes` /
`genotype_prior` / `run_pi_em` / `final_calls` (all index-keyed), and `vcf_out.rs`
`is_variable` / `site_qual` / ALT emission. Confirming these stay untouched is part of the
plan's review, not a code change.

---

## 3. The one new data wire — deconvolved responsibilities for allele balance

The balance test (§2.5) needs, per sample, the deconvolved per-allele counts for its
**called** genotype `(A, B)`. Those come from `q_r(obs|A)` and `q_r(obs|B)`, which the EM
already computes inside `compute_data_ll` (`obs_qr`), but which is local and discarded.

**Decision (Q-I1): compute the deconvolved counts in `em.rs`, at the final E-step, and
attach them to the call — do *not* re-attribute in `vcf_out.rs`.** `final_calls` already
holds `data_ll` and can reach the read model + candidates; after it picks `(A, B)` it
recomputes just the two `q_r` columns for the sample's observed sequences (cheap — a handful
of distinct sequences × 2 alleles), forms `allele_responsibilities`, and stores the summed
`(n_A, n_B)` (or the finished balance penalty) on `SampleCall`:

```rust
struct SampleCall {
    allele_indices: Vec<usize>,
    genotype_units: Vec<u16>,          // kept for REPCN (length annotation)
    posterior: f64,
    gq: u16,
    allele_support: SmallVec<[f64; 2]>, // NEW: deconvolved per-called-allele responsibility
}
```

`vcf_out.rs` then consumes `allele_support` directly — no read model, no `ctx`, no
length attribution at the VCF stage. This keeps the "same `Qᵣ` everywhere" invariant and
localizes the sequence-awareness where the read model already lives.

*Fallback (Q-I1-alt):* recompute `q_r` in `vcf_out.rs` (needs `ctx` threaded through). Only
if attaching to `SampleCall` proves awkward; measured cost is the same (a re-score of the
called two alleles).

---

## 4. Determinism

Unchanged from the Mark-2 contract (byte-identical across and within thread counts):

- `seqs_by_length` and the new per-sequence sample tally are `BTreeMap`s, byte-sorted →
  set-valued nomination emits candidates in a fixed order regardless of thread.
- Candidate assembly, the seed, the EM, attribution, and the balance deconvolution are all
  **per-locus** pure functions of that locus's data — computed identically on whichever
  worker owns the locus (the atomic-per-locus property, Mark-2 §4.4).
- Slip sufficient statistics stay **integer counts** (order-free reduce). The soft
  responsibility sums `(n_A, n_B)` are summed in **fixed per-locus read order** (the sample's
  `seq_counts` is byte-sorted by the Stage-1 contract), so the balance penalty is
  reproducible.
- All same-length tie-breaks resolve on sequence bytes / lower candidate index — never on a
  `HashMap` iteration order or a float race.

---

## 5. Phase 2 — per-allele stutter (module changes, measurement-gated)

Only after Phase 1 lands and the §6.5 measurement settles the model form (spec §8 step 5).
Changes are additive and collapse to Phase 1 when no impure alleles exist.

- **`param_estimation.rs`** — a per-allele **purity measure** (interruption count /
  impure-base fraction — form settled by the §6.5 measurement, Q-I2); the covariate default
  `f(length, period, motif, purity)` fit on the data-rich pooled skirts (purity kept only if
  it survives the fixed-length test, §6.5a); a per-allele shrinkage of the chosen stutter
  parameter toward `f` by read count.
- **`read_model/`** — apply the per-allele term at the existing per-candidate scoring site.
  **Which parameter it attaches to — `level`, the `decay`/shape, or both — is gated by the
  §6.5b measurement (Q-I3), not assumed.** A pure allele at an average locus lands on today's
  arithmetic.
- **`em.rs`** — feed each called allele's slip counts (from §2.4's sequence attribution) into
  the per-allele estimator; shrink toward the covariate default on thin data.

---

## 6. Out of scope here (tracked elsewhere)

- **Same-length paralog FP guard** — a systematic same-length variant (collapsed paralog /
  CNV / mismap) recurs across carriers and clears §5.2; neither `G₀` nor sequence-aware
  allele balance rejects a clean ~½ balance. Needs a coverage-excess / cohort-heterozygosity
  guard adapted from `src/var_calling/paralog_filter/` — **`doc/devel/TODO.txt`**, gated by a
  benchmark FP audit (spec §10).
- **The Mark-2 BIC confident-genotype gate (D1)** — makes the pre-pass gate sequence-aware
  and removes the `ε` contamination Phase 1 leaves; a separate follow-up (spec §5.5, TODO).

---

## 7. Open decisions

| Q | Decision | Where settled |
|---|---|---|
| **Q-I1** allele-balance responsibilities | compute deconvolved `(n_A,n_B)` in `em.rs` final E-step, attach to `SampleCall` (recommended); recompute-in-`vcf_out` is the fallback | this doc §3 |
| **Q-I2** purity measure | interruption count vs impure-base fraction | Phase 2 measurement, spec §6.5 / plan P2.0 |
| **Q-I3** purity attaches to level / decay / both | data decision | spec §6.5b / plan P2.0 |
| **Q-I4** per-sequence sample tally storage | widen `SeqCount` vs a side `BTreeMap<(len,seq),u32>` | P1.1 types |
