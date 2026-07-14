# Tandem-repeat scanner — implementation plan

**Status:** draft, 2026-07-14. The build order for the **tandem-repeat scanner**: the
`src/ng/tandem_repeat.rs` module (its types, the `find_tandem_repeats` interval finder, and the
`RegionScanner` region seam), plus its **validation against the current `trf-mod` catalog**
(parity). Design is settled in [`../spec/ssr_repeat_scanner.md`](../spec/ssr_repeat_scanner.md)
(spec) and [`../arch/ssr_repeat_scanner.md`](../arch/ssr_repeat_scanner.md) (types & interfaces).
This roadmap turns that design into build order; it is **not** a place for new design — every open
question is resolved in the spec §9.

**Scope note (2026-07-14): production stays on `trf-mod` for now.** This plan builds the scanner
and *proves* it reproduces the `trf-mod` catalog, but it **does not** replace the production
dependency — no `postprocess` retype, no `catalog::run` swap, no `trf.rs`/`Containerfile` removal.
Whether and when to make that swap is a **future decision**; the target integration is designed in
spec §6–§7 and stands ready, unexecuted. The scanner is a standalone sequence primitive (spec §1,
§4), a shared primitive **in** the ng tree (not one of the ng pipeline steps).

---

## Scope

**In:**

- `src/ng/tandem_repeat.rs` — all the types (`PeriodRange` + `PeriodRangeError`, `ScanParams`,
  `SegmentOptions`, `RepeatInterval`, `RegionSpan`, `RepeatRegion`, `Region`, `ScanError`); the
  `find_tandem_repeats` interval finder; the `RegionScanner` region seam (resident-slice **and**
  windowed-`ChromRefFetcher` constructors).
- **Validation against `trf-mod` (no production change):** the committed golden catalog captured
  from the current `trf-mod` path, and a parity test that runs the scanner through the *unchanged*
  post-filter (via a test-only bridge) and compares the resulting `Locus` set to the golden.

**Out (future decision / later plans / deferred):**

- **Replacing `trf-mod` in production** — retyping `postprocess` (`TrfRecord` → `RepeatInterval`),
  pointing `catalog::run` at `find_tandem_repeats`, the header `detector` field, and deleting
  `trf.rs` / the `--trf-mod-path`·`--temp-dir` knobs / the `Containerfile` install. **A future
  decision**, designed in spec §6–§7; not executed here.
- **The end-to-end `ssr_tomato1` genotyping check** — it needs a scanner-*built* production catalog,
  which is the deferred swap above; revisit it when that swap is decided. (Locus-set parity, D2,
  already proves catalog equivalence.)
- **The ng snp/str caller consuming the region seam** — the seam's *live* consumer and its
  end-to-end verification; its own plan (spec §1). This plan builds and unit-tests the seam; its
  only downstream use here is self-contained.
- **Stitch pass** (pathological large/multi-indel split), **per-period scoring weights**,
  **phase-aware interrupted-tract retention** (post-filter side) — all deferred until measured
  (spec §10).
- **Removing the vendored `TRF-mod/` tree** — repo-hygiene follow-up (spec §10).

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The algorithmic heart before the plumbing.** The `find_tandem_repeats` finder (B) is built
  and unit-tested before the seam that resolves its output (C); within C, the resident-slice
  coverage merge (C1) before the windowed streaming (C2). The detection/tiling logic is the
  thing that must be right; prove it first, plumb it second.
- **Simplest impl first, as the oracle for the next.** C1 (single-pass tiling over a resident
  slice) is the parity oracle for C2 (windowed streaming): the seam must be **window-count
  invariant**, so C2's output over any `window_bp` must equal C1's single-window output — the
  boundary-halo correctness gate, mirroring the cohort DUST split-invariance test.
- **Reuse over rewrite.** The parity test runs the **unchanged** `postprocess::build_loci` (via a
  test-only `RepeatInterval` → `TrfRecord` bridge), so production is untouched; the seam streams
  over the existing `ChromRefFetcher`; the parity test follows the `src/baq/` golden-fixture
  pattern. No post-filter logic is re-derived.
- **Verify against ground truth.** The north-star test is **catalog parity with the current
  `trf-mod` path** — the scanner→post-filter `Locus` set reproduces the `trf-mod`→post-filter
  `Locus` set (≥ 99 % recall, reviewed diff; spec §6, §9.2), not self-consistency.
- **Incremental, with pauses.** Land one milestone, stop for review, then the next.
- **Ungated / container builds.** `tandem_repeat` is a plain module in `src/ng/`; all `cargo` via
  `./scripts/dev.sh` (CLAUDE.md), a native host build at completion.

## Preconditions (already in place)

- **The current `trf-mod` → `postprocess` path works in the dev container** — `trf-mod` is
  installed (`Containerfile` ~101–114), and `catalog::run` builds a catalog from it. This is the
  **parity oracle source** for the golden fixture (D1). It stays in place — production is not
  changed by this plan.
- **`postprocess::build_loci`** and the `Locus`/`Motif`/`TrfRecord` types
  ([`src/ssr/catalog/postprocess.rs`](../../../../src/ssr/catalog/postprocess.rs),
  `src/ssr/catalog/trf.rs`, `src/ssr/types.rs`) — the unchanged post-filter, the golden's element
  type, and the `TrfRecord::for_test` bridge the parity test feeds scanner intervals through.
- **`ChromRefFetcher` / `StreamingChromRefFetcher`**
  ([`src/fasta/fetcher.rs`](../../../../src/fasta/fetcher.rs)) — the sliding-buffer reader the
  seam (C2) streams over; already used by the cohort/DUST path.
- **The `src/baq/` + `baq_tests.rs` fixture pattern** — the template for the D2 parity test.

---

## The steps

### Milestone A — module scaffold + types (types, no logic)

**A1. Scaffold `src/ng/tandem_repeat.rs`.**  ✅
The file with its `#[cfg(test)]` block; wire `pub(crate) mod tandem_repeat;` into `src/ng/mod.rs`.
One file, no folder — one algorithm, no bake-off; a shared primitive in the ng tree, not a
pipeline step. *Source:* arch §Module home.

**A2. Input types.**  ✅
`PeriodRange` (private `u8` fields, `new(min, max) -> Result<Self, PeriodRangeError>` rejecting
`min == 0` and `min > max`, `min()`/`max()` accessors); `PeriodRangeError` (`#[non_exhaustive]`
`thiserror`, `ZeroMin` / `MinExceedsMax`); `ScanParams` + `SegmentOptions` (plain config, visible
`Default`s — `{2,7,2}` and `{1000, 100_000, 0, 0}`); `DEFAULT_PERIODS = (1, 6)`. Unit tests:
`PeriodRange::new` accepts `(1,6)`/`(6,6)`, rejects `(0,_)` and `(3,2)`; `Default`s match the
named values. *Source:* spec §5, arch §1.1.

**A3. Output types.**  ✅
`RepeatInterval` (`start`/`end` `u32`, `period` `u8`, `score` `i32`); `RegionSpan { start, end }`;
`RepeatRegion { span, intervals }`; the three-kind `Region` enum; `ScanError`
(`#[non_exhaustive]`, `Fetch { source: ChromRefFetchError }`). No logic. *Source:* spec §4, §3.6,
arch §1.2–§1.3, §2.2.

> **Checkpoint A:** the module compiles; `PeriodRange` validation and the `Default`s are tested.
> Pause for review.

### Milestone B — the interval finder (the algorithmic heart, pure)

**B1. Scoring signal + Ruzzo–Tompa maximal segments.**  ☐
The per-period lag-*p* score walk (`+match_reward` on a case-insensitive ACGT match to `seq[j-p]`,
else `-mismatch_penalty`; non-ACGT never matches) and the O(n) Ruzzo–Tompa maximal-scoring-segment
pass over it, as an internal helper returning `[j0, j1]` segments for one period. Unit tests on the
helper: a clean tract yields one segment at the exact bounds; an `N`-run yields none. *Source:*
spec §3.1–§3.2, §3.5, arch §2.1.

**B2. `find_tandem_repeats` — emission over the period range.**  ☐
Loop periods `min..=max`, map each segment to `RepeatInterval { start: j0-p, end: j1+1, period: p,
score }`, apply the `min_copies` emission floor, concatenate. Unit tests (consumer-agnostic, spec
§9): clean `(CAG)k` (exact `start/end/period/score`); a substitution-interrupted tract (**one**
segment); a **single-indel**-interrupted tract (**one** span — the impure-tolerance property); a
large/multi-indel tract (the pathological split into two); a soft-masked (lowercase) tract (found);
a period-2 tract that also matches at p = 4/6 (all emitted — no period de-dup). *Depends:* B1.
*Source:* spec §3.2–§3.5, arch §2.1.

> **Checkpoint B:** the finder is decided and unit-tested in isolation — clean, impure
> (substitution + single indel merged), pathological-split, `N`/soft-mask, multi-period. Pause for
> review.

### Milestone C — the region seam

**C1. Coverage merge → `Region` tiling over a resident slice.**  ☐
`RegionScanner::over_slice(seq, periods, params, opts)` and the `Iterator<Item = Result<Region,
ScanError>>` impl for the resident case: run `find_tandem_repeats`, sort intervals by start, sweep
to union overlapping/abutting intervals into merged spans, classify each `Repeat` vs `Satellite`
by `max_repeat_len`, emit `Unique` for the gaps, apply the `merge_gap` / `min_repeat_len`
smoothing. Unit tests: the tiling contract (ordered, non-overlapping, union == `[0, n)`);
overlapping different-period intervals merge into one `Repeat` carrying both; `Unique, Repeat,
Unique` around a lone repeat; a `> max_repeat_len` run → `Satellite`; `merge_gap`/`min_repeat_len`
smoothing; empty input → nothing. *Depends:* B2. *Source:* spec §3.6, arch §1.3, §2.2.

**C2. Windowed streaming over `ChromRefFetcher`.**  ☐
`RegionScanner::new(fetcher, contig_len, periods, params, opts)`: walk in `window_bp` cores with a
right halo of `max_repeat_len` and a left margin of `periods.max`, run the C1 merge per window,
emit regions whose start falls in the core (halo/attribution dedup), coalesce satellites spanning
windows, and fuse the iterator on a `ChromRefFetchError` (`ScanError::Fetch`). Unit tests: a
`> max_repeat_len` satellite coalesces across windows; **window-count invariance** — identical
regions across two `window_bp` settings and vs `over_slice` (C1 is the oracle); a mid-scan fetch
error surfaces as a terminal `Err`. *Depends:* C1; the `ChromRefFetcher` precondition. *Source:*
spec §3.6, arch §2.2, §3 (error model).

> **Checkpoint C:** the module is complete — both interfaces built and unit-tested; the seam is
> window-count invariant against the resident oracle. Pause for review.

### Milestone D — validate against `trf-mod` (parity, no production change)

**D1. Capture & commit the golden catalog.**  ☐
Run the current `trf-mod` → `postprocess` path in the container on a small tomato-reference subset
(spec §9.5) and commit the resulting catalog TSV as a test fixture (under `tests/data/…` or the
module's fixtures, per the `baq` pattern). Production is untouched — this only *reads* the existing
path to snapshot the parity **oracle**. *Depends:* none (existing code). *Source:* spec §6, §9.

**D2. Golden-fixture parity test — the port anchor.**  ☐
In a test (living in `src/ssr/catalog`, where it can see `build_loci` + `TrfRecord`), run
`find_tandem_repeats` on the D1 reference, bridge each `RepeatInterval` to the existing post-filter
input via `TrfRecord::for_test(start, end, period as u16, score, b"")` (a **test-only** adapter —
no production retype), run the **unchanged** `build_loci`, and assert the resulting `Locus` set
reproduces the golden set at **≥ 99 % recall** with a reviewed, understood diff (any miss traced to
a concrete cause; extra loci inspected as genuine STRs). Follows the `baq_tests.rs` fixture pattern.
*Depends:* D1, B2. *Source:* spec §6, §9.2.

> **Checkpoint D:** the scanner is *proven* to reproduce the `trf-mod` catalog at the agreed bar,
> with **zero change to the production `trf-mod` path**. The module is built, tested, and validated;
> whether to swap it into production is a separate future decision (spec §6–§7). Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | type-level tests — `PeriodRange::new` validation; `Default`s = the named values |
| B | consumer-agnostic unit tests: clean / substitution / single-indel (all one interval), multi-indel split, `N`, soft-mask, multi-period emission |
| C | seam unit tests: tiling contract, multi-period merge, satellite (+ coalescing), smoothing; **window-count invariance vs the `over_slice` oracle**; terminal fetch-error |
| D | **catalog parity** vs the committed `trf-mod` golden (`Locus`-set recall ≥ 99 %, reviewed diff) through the *unchanged* post-filter — production left on `trf-mod` |

## Out of scope (future decision / next plans / deferred)

- **Replacing `trf-mod` in production** — the `postprocess` retype, `catalog::run` swap, header
  `detector` field, and `trf.rs` / `Containerfile` removal — a **future decision**, designed and
  ready in spec §6–§7, not executed here.
- **The end-to-end `ssr_tomato1` genotyping check** — needs a scanner-built production catalog (the
  deferred swap); revisit with that decision.
- **The ng snp/str caller** consuming the region seam — its own plan (spec §1).
- **Stitch pass** (pathological multi-indel), **per-period weights**, **phase-aware interrupted-
  tract retention** (post-filter side) — deferred until measured (spec §10).
- **Removing the vendored `TRF-mod/` tree** — repo-hygiene follow-up (spec §10).
- **Scoring-weight tuning** — start `2/7`; the empirical sweep against the D2 fixture refines it
  (spec §3.3, §9.3). A tuning follow-up, not a correctness gate.
