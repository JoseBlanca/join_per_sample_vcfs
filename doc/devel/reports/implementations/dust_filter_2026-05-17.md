# Stage 3 — sdust low-complexity filter (implementation report)

Date: 2026-05-17.

Implements the Stage 3 low-complexity filter described in
[doc/devel/specs/calling_pipeline_architecture.md §"Stage 3"](../../specs/calling_pipeline_architecture.md#L875),
following the plan at
[doc/devel/implementation_plans/dust_filter.md](../../implementation_plans/dust_filter.md).

## Plan

Two layers in a single new module
[src/var_calling/dust_filter.rs](../../../src/var_calling/dust_filter.rs):

1. **`sdust_mask`** — pure function over a `&[u8]` reference slice
   that returns the sorted, non-overlapping list of low-complexity
   intervals. Direct port of `sdust_core` from `lh3/sdust` (vendored
   at [sdust/](../../../sdust/), gitignored). Triplet-pair scoring,
   suffix-trim invariant, perfect-interval candidate list, N-break
   convention — all transcribed line-by-line from
   [sdust/sdust.c](../../../sdust/sdust.c).
2. **`DustFilter<I, F>`** — iterator adaptor that wraps an upstream
   `Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>`
   plus a `RefSeqFetcher`. On the first item for each chromosome it
   loads the whole reference and runs `sdust_mask` once, then answers
   pass/skip per upstream position by sweeping the resulting interval
   list (`O(1)` amortised per position).

## Assumptions

All silent choices were already documented in the plan; in the
implementation they came out as:

- **Window is left-anchored.** `window_left == p` for the decision
  on `p`; any-of-`w`-windows-containing-`p` semantics are emergent
  from how `find_perfect` accumulates per-position decisions as the
  window slides past.
- **Strict `>` threshold.** Mask iff `10·s > T·L`, matching
  [sdust.c:111](../../../sdust/sdust.c#L111). Pinned by a dedicated
  off-by-one unit test (the high-threshold sensitivity smoke at
  `T = 5000`).
- **Lowercase ACGT = uppercase ACGT.** Direct port of
  `seq_nt4_table` from [sdust.c:22-39](../../../sdust/sdust.c#L22).
  Pinned by `sdust_lowercase_treated_as_uppercase`.
- **N (or any non-ACGT byte) breaks the triplet stream, but does
  *not* reset the window/count/score state.** Direct port of
  [sdust.c:151-155](../../../sdust/sdust.c#L151). The C source
  resets only `l` (run length) and `t` (rolling word) at the N
  branch; stale window state shifts out as new triplets are pushed.
  The plan's "Assumptions" section flagged this as a faithful-port
  choice and the golden-vector test confirms it produces the same
  output as the binary.
- **Config is validated at construction.** `DustFilterConfig::new`
  is the fallible constructor; `DustFilter::new` is infallible.
  Mirrors the direction the Stage 6 review picked up for
  `PosteriorEngineConfig` (Mi12).
- **No bypass mode inside `DustFilter`.** `--no-complexity-filter`
  is the caller's concern; if a `DustFilter` exists, it filters.
- **`T = 20` at the API surface, not the equivalent density `2.0`.**
  Matches the universal sdust convention.
- **Per-chromosome batch + sweep**, not a streaming sdust integrated
  with the upstream `(chrom_id, pos)` iterator. Substantially
  simpler implementation and matches the canonical reference
  exactly. Memory cost is one chromosome's reference bytes
  transient plus the masked-interval list resident; both bounded.

## Changes made

- **New module**:
  [src/var_calling/dust_filter.rs](../../../src/var_calling/dust_filter.rs)
  (≈1180 lines including tests). Public surface:
  - `sdust_mask(seq: &[u8], window: u32, threshold: u32) -> SdustIntervals`
  - `DustFilter<I, F>` + `DustFilter::new` + `DustFilter::config`
  - `DustFilterConfig` (private fields, `new`/`window`/`threshold`/`Default`)
  - `DustFilterError` (`Upstream`, `RefFetch`, `UnknownChromId`, `InvalidWindow`; `#[non_exhaustive]`)
  - `SdustIntervals = Vec<(u32, u32)>` type alias
  - `DEFAULT_DUST_WINDOW = 64`, `DEFAULT_DUST_THRESHOLD = 20` constants
- **Module wiring**: [src/var_calling/mod.rs](../../../src/var_calling/mod.rs)
  gained `pub mod dust_filter;` between `contamination_estimation`
  and `per_group_merger`.
- **Golden-vector fixture**: a `GOLDEN_SNIPPETS` const inside
  [dust_filter.rs](../../../src/var_calling/dust_filter.rs)
  carrying six hand-picked snippets and their expected
  interval lists. The `expected` values were generated once at
  implementation time by running `lh3/sdust -w 64 -t 20` against
  each snippet and copying its BED-style output. After that
  one-time check, the committed values are the source of truth —
  `lh3/sdust` itself is *not* a build- or test-time dependency
  of this project.

## Tests added

25 new tests in `var_calling::dust_filter::tests`:

### `sdust_mask` (algorithmic core)

1. `sdust_empty_input` — empty `&[u8]` → empty interval list.
2. `sdust_high_complexity_passes` — deterministic xorshift-random
   ACGT 256-bp string at default config → no masked intervals.
3. `sdust_homopolymer_is_masked` — 128-bp poly-A → single interval
   `(0, 128)`.
4. `sdust_dinucleotide_repeat_is_masked` — 64×"AT" → masked.
5. `sdust_trinucleotide_repeat_is_masked` — 42×"ATG" → masked.
6. `sdust_lowercase_treated_as_uppercase` — uppercase and
   lowercase inputs produce identical masks.
7. `sdust_n_breaks_the_run` — two A-runs separated by an N → at
   least two distinct intervals, neither spanning the N.
8. `sdust_other_non_acgt_also_breaks_run` — confirms any non-ACGT
   byte (here whitespace) behaves like N.
9. `sdust_high_threshold_disables_masking` — `T = 5000` on a
   homopolymer → empty interval list. Pins that the threshold knob
   is wired through and that `>` is strict (a window scoring just
   at-threshold would still mask if the comparison were `≥`).
10. `sdust_intervals_are_sorted_and_disjoint` — multiple
    low-complexity regions in one slice come out in
    increasing-start order with no overlaps.
11. `sdust_mask_panics_on_tiny_window` — `window = 2` panics with
    "window must be" message.

### Golden-vector

12. `sdust_matches_vendored_binary_on_golden_snippets` — runs
    `sdust_mask` on six hand-picked snippets (homopolymer with
    flanks, dinucleotide and trinucleotide repeats, short
    homopolymer at a boundary, N inside a low-complexity tract,
    long CAG-style trinucleotide tract) and asserts byte-identical
    output to the vendored `lh3/sdust` binary at `W=64 T=20`. This
    is the load-bearing test for port fidelity.

### `DustFilterConfig`

13. `config_new_accepts_defaults` — `new(64, 20)` succeeds with
    the right getters.
14. `config_new_accepts_boundaries` — `new(3, _)` and
    `new(4096, _)` both succeed.
15. `config_new_rejects_too_small_window` — `new(2, _)` →
    `Err(InvalidWindow { window: 2 })`.
16. `config_new_rejects_too_large_window` — `new(4097, _)` →
    `Err(InvalidWindow { window: 4097 })`.
17. `config_default_passes_validation` — `Default::default()`
    yields the documented defaults.

### `DustFilter` (iterator plumbing)

Uses a `StubFetcher` (per-chromosome `Vec<u8>`, optional
fail-for-chrom, call-count tracking) and a `std::iter`-backed
upstream.

18. `filter_passes_through_high_complexity` — every upstream
    position survives on a random ACGT chromosome.
19. `filter_drops_homopolymer_positions` — every upstream position
    is dropped on a chromosome that's entirely poly-A.
20. `filter_fetcher_call_count_pinned` — 50 positions on chr0 +
    50 on chr1 trigger exactly one fetch per chromosome.
21. `filter_resets_mask_across_chromosomes` — same `pos` is
    dropped on a low-complexity chromosome and passes on a
    high-complexity one.
22. `filter_surfaces_upstream_error_and_latches` — upstream
    `Err(...)` propagates as `DustFilterError::Upstream(_)`; every
    subsequent `next()` returns `None`.
23. `filter_surfaces_ref_fetch_error_and_latches` — fetcher
    `io::Error` propagates as
    `DustFilterError::RefFetch { chrom_id, .. }` and latches.
24. `filter_surfaces_unknown_chrom_id_and_latches` — upstream
    chrom_id not in the chromosomes table → `UnknownChromId` and
    latches.
25. `filter_empty_upstream_does_no_work` — empty upstream yields
    nothing and never calls the fetcher.

## Validation results

Run inside the dev container (`./scripts/dev.sh`):

| Command | Result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean (6 warnings fixed: 1× `collapsible_if` in lib via let-chain, 3× `repeat_n` and 1× `io::Error::other` in tests) |
| `cargo test --all-targets --all-features` | 693 lib + 109 integration tests pass (was 668 lib + 109 integration, +25 from this module) |

## Tradeoffs and follow-ups

- **Per-chromosome batch is not streaming-pure.** The plan documents
  this trade: simpler implementation, exact match against
  `lh3/sdust`, transient memory of one chromosome's reference bytes
  during the `sdust_mask` call. The walker already holds one
  chromosome's reference in memory, so this is the established
  footprint, not new ground. Streaming sdust integrated with the
  upstream iterator is listed as an out-of-scope follow-up to
  revisit only if profiling shows the transient cost matters.
- **Golden-vector fixture is the committed source of truth.**
  The six snippets' `expected` intervals were generated once at
  implementation time against `lh3/sdust`; from then on, the
  committed values are what the test asserts against. This
  intentionally severs the project's dev tooling from any external
  C dependency — adding new snippets later is an off-repo task
  (set up `lh3/sdust` independently, run it, paste the output).
  Trade: a new-snippet addition with wrong `expected` values would
  not be caught automatically; code review is the safety net.
- **CLI parser bindings deferred.** The library API ships; the
  `--complexity-window` / `--complexity-threshold` /
  `--no-complexity-filter` flags will land with the cohort
  subcommand, same convention as Stages 5 and 6.
- **End-to-end PspReader → … → DustFilter integration test
  deferred.** Same dependency on the cohort CLI fixture as Stages
  5 and 6 — listed in the Stage 3 block's Open items for that
  later PR.
- **N-break does not reset window/count state.** This matches
  `lh3/sdust` exactly, and the golden-vector test pins it. If the
  semantics ever turn out to matter for cohort-level FPR, revisit
  by adding a per-N hard-reset variant guarded by a config knob —
  but only with evidence.
- **No benchmark added.** Plan called out that benching the filter
  in isolation against a synthetic input would over-fit; a real
  bench is justified once the cohort CLI runs end-to-end.
