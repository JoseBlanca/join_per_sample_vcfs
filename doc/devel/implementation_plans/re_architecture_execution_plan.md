# Re-architecture — phased execution plan

Companion to the two design docs — read them first:
- [re_architecture_streaming_pipeline.md](re_architecture_streaming_pipeline.md) — the architecture, constraints, build strategy.
- [re_architecture_module_outline.md](re_architecture_module_outline.md) — the modules, types, and functions (the map).

This document turns the plan's §7 build strategy + rollout sketch into
**ordered, byte-identity-gated, independently-reviewable phases**. It says
*when* to build *what*, and *how each phase is verified* — not the design
(that's the docs above).

---

## Principles (the spine of every phase)

1. **Parallel package.** Build in a fresh **`var_calling_new::`**; the old
   `var_calling::` keeps shipping (CLI stays on it) until the final swap.
2. **Rebuild structure, copy kernels.** Build the *structure* from scratch
   (readers, producer, wiring, record-based grouping); **copy the
   byte-identity-sensitive numeric kernels verbatim** (the SIMD EM —
   `run_em_columnar` / `compute_log_likelihoods_columnar` / `posterior_engine`
   / `kernels` — and the `per_group_merger` numerics). Rewriting proven
   floating-point math is how byte-identity is lost.
3. **The old package is the live oracle.** The byte-identity gate is an
   in-process test that runs **both** pipelines on the same `.psp` cohort and
   diffs the VCFs (strip `^##`, md5, QUAL excluded — QUAL is non-deterministic
   since `03e2221`). The retained `var_calling::` is the branch base (= `main`
   @ `a71a078`), so "new vs old" *is* "new vs `main`." **Stand the harness up
   in P0; it must be green from P4 onward.**
4. **Each phase leaves the tree green** (`cargo fmt --check`, `clippy
   --all-targets --all-features -D warnings`, `cargo test`) and is reviewable
   on its own. Within a phase: **types first, then implementation**, with a
   pause to review before the next phase.
5. **Correctness before optimisation.** Reach byte-identity with a simple
   decode first (P1–P4); the column-selective read lever (P5) and the deferred
   memory levers (post-P6) come *after*, each re-gated.
6. **Measure, then decide.** The memory/perf levers (byte-budget queue #8,
   low-memory mode #9) are **not built on spec** — only after P6 measures the
   new architecture against `main` and the numbers call for them.

**Invariant to hold from P2 on** (plan §2.2 / appendix §B): the cohort-wide
in-memory data is **only the consolidated records — variable positions, AC /
`min_alt_obs` already applied**. Never materialise a full-coverage cohort
structure.

---

## Phase 0 — Scaffold, kernels, oracle, re-profile

**Goal:** the empty new package, the copied kernels, and the byte-identity
harness — before any new production logic.

**Build:**
- Create the `var_calling_new::` package (flat module skeleton per appendix
  §A–§G: `sample_reader`, `cohort_integration`, `types`, `pileup_overlaps`,
  `em_posterior_calc`, `vcf_writer`, `pipeline`, `dust`).
- **Copy the verbatim kernels** into it: `per_group_merger` numerics,
  `posterior_engine` + `kernels` (the SIMD EM), and the `#[cfg(test)]`
  row-shape grouper (`build_overlapping_variant_group`) + `OverlappingVariantGroup`
  / `MergedRecord` — to be promoted to production in P3.
- **Oracle harness**: an integration test (or example) that runs old
  `run_var_calling` and the new entry on the same `.psp` fixtures
  (3-tomato + the multichrom fixture) and diffs the VCFs. It *fails* until P4
  — that's expected; it's the scaffold the later phases turn green.
- **Re-profile at representative N** (the plan's measurement-first item): pull
  the `[profile.profiling]` build + `tmp/lowmem_measure/dhat_attrib.py`; record
  the `main` baseline (RSS + wall, N=50, T=4/T=8, default) for P6 to compare
  against. (No production code change.)

**Gate:** package compiles; copied kernels pass their copied unit tests; the
oracle harness builds and runs (red, by design).

---

## Phase 1 — Data types & per-sample reader

**Goal:** the data model and the per-sample read path.

**Build (types first):**
- Types (appendix §C): `CohortPileupRecord`, `PileupCohortChunk`
  (`{ chunk_order, records, ref_span }`), `RefSpan`, `Variant`, `CalledChunk`.
- `SamplePspReader` (appendix §A) — per sample, over `psp::PspReader`;
  `peek_next_span` (segment boundary) + segment read. `!Send`.
- `SamplePspChunk` — one segment; light accessors (`positions` / `nonref_obs`
  / `ref_spans`, decoded at `new`) + the typed getters
  (`take_seq`/`take_chain_ids`/`take_fixed`). **Decode may be simple here**
  (decode the needed columns correctly; the *skip-the-rest* optimisation is
  P5).

**Gate:** unit tests — decode a known `.psp` segment; assert the light
accessors and `take_*` results equal a full decode of the same block (reuse
the psp test fixtures). Tree green.

---

## Phase 2 — Producer (the record stream)

**Goal:** `CohortChunkIntegrator` emits `PileupCohortChunk`s (record streams).

**Build:**
- `CohortPerPositionMerge` — **revive `per_position_merger`** as the cohort
  join: walk all samples' light columns by position → variable positions
  (AC / `min_alt_obs`, reusing `derive_is_kept`'s rule), then build one
  `CohortPileupRecord` per variable position (heavy via `take_*`).
- Dust: apply at step 2 (drop masked positions, after the obs decision).
- Safe-gap cut (`find_block_cut`, reused) → chunk boundaries between whole
  records; `RefSpan` fetch (producer-side, monotonic-forward); `chunk_order`
  stamping.
- The `run` loop: lockstep segment read → merge → records → cut → ship.

**Gate:** unit tests against the oracle's intermediate — the set of variable
positions (and their per-sample records) the producer emits must match what
old `var_calling::` keeps (its keep-mask + per-position view) for the
fixtures. Confirm the cohort-wide-variable-only invariant holds (no
full-coverage structure). Tree green.

---

## Phase 3 — Caller (grouping + EM)

**Goal:** `VariantCaller` turns one chunk into `Variant`s.

**Build:**
- `pileup_overlaps` — **promote** `build_overlapping_variant_group` to
  production: the chunk's `CohortPileupRecord`s → `OverlappingPileupRecords`
  (bridging gaps left by dropped dust positions).
- `em_posterior_calc` — wire the **copied** `per_group_merger` + SIMD EM
  verbatim: `OverlappingPileupRecords` → `Variant`s (per-record EM, SoA +
  SIMD kept). REF sliced per group from the chunk's `RefSpan`.
- `VariantCaller::call_chunk` composing the two; `run` loop emitting one
  `CalledChunk` per chunk (**even when empty**, for the gapless order).

**Gate:** unit tests — for a fixed chunk, the `Variant`s equal the old code's
output for the same positions (kernels are copied, so this should hold by
construction). Tree green.

---

## Phase 4 — Writer + wiring → **byte-identity milestone**

**Goal:** first end-to-end run; the oracle turns green.

**Build:**
- `VcfWriter` — `CalledChunk` stream → VCF, reorder by `chunk_order`
  (`BTreeMap`), with the contiguity debug-assert.
- `pipeline` — spawn producer + W callers + writer; two bounded
  `crossbeam-channel` hand-offs (**count cap**, not byte-budget); first-error
  propagation; the new entry point.

**Gate (the headline A/B point):** the **oracle goes green** — new vs old
VCFs byte-identical (QUAL excluded) on the 3-tomato + multichrom fixtures, at
T=1 and T=8, default mode. This is where byte-identity is *proven*. Tree
green; full new-package test suite passes.

---

## Phase 5 — Column-selective decode (the read lever)

**Goal:** decode only what's needed — light to decide variable, heavy only
for variable positions; skip the rest in the psp block.

**Build:** extend the psp decode (`decode_block_payload`) to decompress/decode
only the requested columns and seek past the others (manifest `compressed_len`);
wire `SamplePspChunk`'s `take_*` / light accessors onto it (plan §4.1).

**Gate:** **oracle still green** (byte-identical — column-selective must
produce the same `CohortPileupRecord`s for variable positions). Record the
read-side reduction (decode volume / allocator) vs P4.

---

## Phase 6 — Measure vs `main` → decide deferred levers

**Goal:** the plan's measurement step (plan §6) — the scorecard.

**Do:** wall-time + peak RSS at N=50 whole-genome, T=4 and T=8, default mode,
new vs the P0 `main` baseline. Note: the cohort-wide-variable-only saving
already exists on `main`, so it's the floor, not a new win; the genuine new
headroom is column-selective not unpacking heavy columns for dropped
positions.

**Decide (data-driven):**
- byte-budget queue admission (#8, plan §4.2) — build only if the count-cap's
  `N × threads` RAM is a problem;
- low-memory `SamplePspChunk` re-read mode (#9) — build only if RAM still too
  high;
- the DUST-coupling confirm (plan §8) — verify column-selective didn't change
  the mask coupling (expected: no).

Anything not demanded by the numbers stays unbuilt.

---

## Phase 7 — Swap (the cutover)

**Goal:** `var_calling_new::` becomes `var_calling::`.

**Do, in one commit:**
- delete old `var_calling::` + its now-redundant **unit** tests;
- **port** the correctness / byte-identity **integration** tests to the new
  entry point (retarget, don't delete);
- rename `var_calling_new::` → `var_calling::`; repoint the CLI;
- keep the oracle test only if a frozen golden VCF is retained (the old
  in-process oracle is gone once old code is deleted — snapshot a golden VCF
  in P4/P6 to guard the swap).

**Gate:** full suite green; one final VCF diff against the P6 golden snapshot.

---

## Dependency / ordering notes

- **P1 → P2 → P3** can overlap once their types (P1) exist, but P4 needs all
  three. P5 needs P4 (optimise on a byte-identical baseline). P6 needs P5
  (measure the final read path). P7 is last.
- **Byte-identity is provable only at P4** (full pipeline → VCF). Earlier
  phases gate on component tests + the copied-kernel guarantee.
- **The deferred memory levers (#8, #9) are post-P6**, gated on measurement —
  they are not phases of the core build.

## Out of scope (unchanged by this plan)

- The calling math, EM/posterior engine, VCF schema (plan §9 non-goals).
- `contamination_estimation` and the `psp::` format layer (used, not
  rebuilt — except P5 extends `decode_block_payload`).
