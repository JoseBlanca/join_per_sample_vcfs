# Stage 3 sdust filter — review-fixes applied

Date: 2026-05-17.
Review: [dust_filter_2026-05-17.md](dust_filter_2026-05-17.md).
Implementation report: [dust_filter_2026-05-17.md](../implementations/dust_filter_2026-05-17.md).

## Status summary

| Severity | Total | Applied | Applied-with-adaptation | Won't-fix | Deferred |
|---|---|---|---|---|---|
| Blocker | 0 | 0 | 0 | 0 | 0 |
| Major | 8 | 7 | 0 | 1 (M4) | 0 |
| Minor | 24 | 22 | 1 (Mi21) | 1 (Mi6) | 0 |
| Nits | 9 | 7 | 0 | 2 | 0 |

All four open questions from the review's §4 were resolved with the human PM before code changes started:

1. **Q1 / Mi6 — tracing policy.** Won't fix; project convention is typed errors only, no informational tracing.
2. **Q2 / Mi7 — `4096` upper bound.** Applied as **(b)**: `MAX_DUST_WINDOW = u16::MAX as u32 = 65535`, anchored at the `u32` score-overflow threshold (`u·(u−1)/2` crosses `u32::MAX` at `u ≈ 92682`).
3. **Q3 / Mi24 — `sdust_mask` dual-validation.** Applied as **(c)**: `sdust_mask` made private; `SdustIntervals` type alias also made private; the runtime `assert!` became a `debug_assert!` since the function is no longer publicly reachable.
4. **Q4 / M4 — out-of-order upstream.** Won't fix; the merger's contract is trusted at the layer boundary; module-level doc comment added explaining the assumption.

## Findings

### Major (7 of 8 applied, 1 won't-fix)

#### M1 — `pileups.pos - 1` underflow on `pos == 0` — **Applied**

[src/var_calling/dust_filter.rs:752-760](../../../src/var_calling/dust_filter.rs#L752). New typed-error variant
`DustFilterError::InvalidPos { chrom_id }` added to the enum; the conversion site in
`DustFilter::next` now uses `pileups.pos.checked_sub(1)` and surfaces the error
on contract violation, latching `self.is_finished` so subsequent `next()` calls
return `None`. Regression test
[`filter_surfaces_invalid_pos_zero_and_latches`](../../../src/var_calling/dust_filter.rs#L1455).

#### M2 — No property/fuzz coverage on `sdust_mask` — **Applied**

[src/var_calling/dust_filter.rs:1380-1438](../../../src/var_calling/dust_filter.rs#L1380). The project does not currently
depend on `proptest`, so the stand-in is a seeded-random sweep
(`sdust_invariants_hold_on_random_seeded_inputs`) that exercises ~40,000
combinations of `(seed, alphabet, len, window, threshold)` and asserts the
two structural invariants the iterator layer relies on: each interval is
non-empty, and the list is sorted with no overlaps. A new
`bases_from_seed(seed, n, alphabet)` helper supports both this sweep and the
pre-existing `high_complexity_bases` generator.

The sweep surfaced two faithful-port properties of `lh3/sdust` that were
not obvious from the C source alone:

1. Interval `finish` can legitimately exceed `seq.len()` after an N break
   (the deque keeps stale triplets while the run-origin resets). Verified
   against the cloned binary on `NCCTGANCGTTCCTGNAGCNNNTGGTNNTGN` with
   `-w 32 -t 1` — binary outputs `7 35` for a 31-byte input, identical to
   our port.
2. Interval `start` can also exceed `seq.len()` on N-heavy inputs at low
   thresholds. Same input at 128 bytes produces `(129, 133)` and
   `(140, 145)` from the binary.

Both are harmless for the only consumer (`DustFilter::is_masked` compares
positions, never slices `seq`). Production upstream never emits positions
past chromosome length, so the out-of-bounds interval coordinates are
unobservable. The assertions are documented inline so future readers don't
re-flag the property.

#### M3 — No half-open boundary test on `is_masked` — **Applied**

[`filter_half_open_boundary_pos_at_start_masked_pos_at_finish_passes`](../../../src/var_calling/dust_filter.rs#L1437).
Uses the `homopolymer_with_flanks` golden snippet (mask `(50, 100)`) and
probes 1-based positions 50, 51, 100, 101 — the first and fourth pass,
the second and third drop. Pins both `start` and `finish` boundary
semantics in one focused test.

#### M4 — No out-of-order regression test on `DustFilter::next` — **Won't fix**

Per Q4 resolution. The merger's `PerPositionMergerError::OutOfOrder`
already guards the monotonicity property at the layer above. Detecting the
same violation again at the filter layer is just defense in depth at a
1:1 transform, which the project decided not to do. Module-level doc note
added: see
[src/var_calling/dust_filter.rs:50-59](../../../src/var_calling/dust_filter.rs#L50)
and
[src/var_calling/dust_filter.rs:721-728](../../../src/var_calling/dust_filter.rs#L721).

#### M5 — `sdust_high_threshold_disables_masking` does not pin strict-`>` — **Applied**

New test
[`sdust_threshold_is_strictly_greater_at_density_boundary`](../../../src/var_calling/dust_filter.rs#L1316)
pins the strict-`>` direction with two assertions: 6 same bases must NOT
mask (their largest candidate sits at density exactly `T/10 = 2.0`), and
7 same bases MUST mask (density `2.5 > 2.0`). Verified against
`lh3/sdust` binary: `AAAAAA` → no output, `AAAAAAA` → `0 7`. If a future
refactor flipped `>` to `>=`, the first assertion would fail. The old
`sdust_high_threshold_disables_masking` smoke is kept as-is — its purpose
is now correctly named ("knob is wired through"), not "strict > pin."

#### M6 — Error `Display` doubles the source chain — **Applied**

[src/var_calling/dust_filter.rs:494-505](../../../src/var_calling/dust_filter.rs#L494).
`#[error("upstream merger failed: {0}")]` → `#[error("upstream merger failed")]`.
`#[error("reference fetch failed for chrom {chrom_id}: {source}")]` →
`#[error("reference fetch failed for chrom {chrom_id}")]`. Causes flow only
through `source()`; chain renderers (`anyhow`,
`std::error::Error::source`) no longer print the cause twice.

#### M7 — `expect("deque non-empty at cap")` lacks `// PANIC-FREE:` — **Applied**

[src/var_calling/dust_filter.rs:200-204](../../../src/var_calling/dust_filter.rs#L200).
Comment added explaining the invariant that keeps the deque non-empty at
this site (we just tested `window.len() >= cap` and `cap >= 1`, both
guaranteed by upstream validation).

#### M8 — Unchecked counter subtractions — **Applied**

Three `// INVARIANT:` comment blocks added at the load-bearing decrement
sites:

- [shift_window eviction (lines 207-211)](../../../src/var_calling/dust_filter.rs#L207) — `cw[s]` decrement matches a prior `shift_window` push.
- [shift_window trim-loop (lines 222-235)](../../../src/var_calling/dust_filter.rs#L222) — added `debug_assert!(self.big_l as usize <= self.window.len())` before the `window.len() - big_l` subtraction; `cv[s]` decrement matches the suffix-entry increment.
- [find_perfect (line 308)](../../../src/var_calling/dust_filter.rs#L308) — `INVARIANT: i < prefix_len <= win_len` annotates the `win_len - i - 1` subtraction.

Plus the convergent Mi19 — a `debug_assert!(self.big_l <= win_len, ...)`
at the top of `find_perfect` ([line 295](../../../src/var_calling/dust_filter.rs#L295))
that documents and enforces the trim-invariant property under tests.

### Minor (22 of 24 applied; 1 applied-with-adaptation; 1 won't-fix)

#### Mi1 — Named constants for window bounds — **Applied**

`DustFilterConfig::new` and the `InvalidWindow` error message now use
`SD_WLEN..=MAX_DUST_WINDOW` rather than the previous bare `3..=4096`
literal. Test assertions follow:
[src/var_calling/dust_filter.rs:1227-1247](../../../src/var_calling/dust_filter.rs#L1227).

#### Mi2 — Rename `res` → `masked_intervals` — **Applied**

[SdustState field at line 167](../../../src/var_calling/dust_filter.rs#L167).
All call sites updated to use the explicit name; no more `state.res` at
return sites.

#### Mi3 — Rename `done: bool` → `is_finished` — **Applied**

[DustFilter field at line 540](../../../src/var_calling/dust_filter.rs#L540).
All set/check sites updated; the latching contract is now a predicate-named
field.

#### Mi4 — Rename `high_complexity` → `high_complexity_bases` — **Applied**

Test helper at [line 800](../../../src/var_calling/dust_filter.rs#L800);
~10 call sites updated.

#### Mi5 — Rename `chrom` → `make_chrom` — **Applied**

Test helper at [line 1119](../../../src/var_calling/dust_filter.rs#L1119);
~15 call sites updated. Now a verb, consistent with builder convention.

#### Mi6 — Tracing on default config use — **Won't fix**

Q1 resolution. Project convention is typed errors only, no informational
tracing; using a default config is not a contract violation. Operators
inspect the effective config via `DustFilter::config()` or the (future)
saved run-metadata file when the cohort CLI lands. The defaults-rule was
genuinely violated in spirit; the broader "no logs" project preference
overrides it for non-invariant-violation sites.

#### Mi7 — `4096` upper bound not anchored — **Applied (option b)**

Q2 resolution. New `const MAX_DUST_WINDOW: u32 = u16::MAX as u32;` at
[line 92](../../../src/var_calling/dust_filter.rs#L92) with a doc comment
explaining the `u32` score-overflow threshold (`u·(u−1)/2 → u32::MAX` at
`u ≈ 92682`); rounded down to `u16::MAX = 65535` for a clean number well
below the boundary.

#### Mi8 — `DustFilterConfig` field-level docs — **Applied**

[src/var_calling/dust_filter.rs:455-465](../../../src/var_calling/dust_filter.rs#L455).
`window` and `threshold` private fields now carry `///` docs explaining
each is the `W`/`T` sdust parameter, with links to the matching default
constants.

#### Mi9 — Public-accessor docs + `# Errors` — **Applied**

`DustFilterConfig::new` got a `# Errors` block;
`DustFilterConfig::window`, `DustFilterConfig::threshold`, and
`DustFilter::config` got one-line `///` summaries.

#### Mi10 — `let mut merged` bool → `if let … else` — **Applied**

[save_masked_regions lines 245-260](../../../src/var_calling/dust_filter.rs#L245).
The mutable bool is gone; the merge-or-push decision is now a single
`if let … && … { … } else { … }` block.

#### Mi11 — `while let` + `if/else break` → `while … is_some_and(…)` — **Applied**

Same `save_masked_regions` block. The truncation loop is now
`while self.perf.last().is_some_and(|p| p.start < start_threshold)` —
one line, no inner `break`.

#### Mi12 — Manual `len()` + index → `Vec::last` — **Applied**

Same block. `let Some(&candidate) = self.perf.last() else { return; }` —
no more `let n = self.perf.len(); if n == 0 { return; } let candidate =
self.perf[n - 1];`.

#### Mi13 — `i64` reverse loop → `(0..prefix_len).rev()` — **Applied**

[find_perfect lines 304-348](../../../src/var_calling/dust_filter.rs#L304).
The signed cursor and `i as u32` / `i as usize` casts are gone. Loop
iterates over `u32` directly; `prefix_len = win_len.saturating_sub(big_l)`
expresses the bound in plain terms. The early-return guard is also gone —
an empty range covers it.

#### Mi14 — Bundle `current_chrom` / `current_mask` / `sweep` → `Option<LoadedChrom>` — **Applied**

New private `LoadedChrom { chrom_id, mask, sweep }` struct at
[line 523](../../../src/var_calling/dust_filter.rs#L523); the three
co-dependent fields are now a single `Option<LoadedChrom>` updated in one
atomic assignment in `ensure_mask_for`. The truth-table-with-impossible-rows
collapse mentioned in the review is now a type-level property.

#### Mi15 — Test pattern bare `..` → `source: _` — **Applied**

`filter_surfaces_ref_fetch_error_and_latches` at
[line 1280](../../../src/var_calling/dust_filter.rs#L1280) now uses
`{ chrom_id: 0, source: _ }`. A future field added to the variant fails
compile here.

#### Mi17 — Extract `dust_density_exceeds` helper — **Applied**

[Lines 174-182](../../../src/var_calling/dust_filter.rs#L174).
The `10·X > T·Y` predicate with `u64` widening is now a single inlined
helper; all three former open-coded sites (suffix-trim invariant,
`find_perfect` density gate, whole-window density gate) call it. The
helper's doc comment cites the three matching C source lines so the
traceability is preserved.

#### Mi18 — Hoist `self.done` from closures — **Applied**

`ensure_mask_for` at [lines 658-680](../../../src/var_calling/dust_filter.rs#L658)
now uses explicit `match` arms with `self.is_finished = true;` assignments
on dedicated lines, not buried inside `.ok_or_else` / `.map_err`
closures. The latching contract is now visible at a glance.

#### Mi19 — `debug_assert!` on `big_l <= window.len()` — **Applied**

Convergent with M8. See above.

#### Mi20 — `lh3/sdust` commit SHA in golden-snippet comment — **Applied**

[Lines 952-958](../../../src/var_calling/dust_filter.rs#L952). The comment
block now records that the values were generated at upstream commit
`89c42cb41ba598e9cfa07c2ef99ae8c08f769b3e`.

#### Mi21 — `cargo doc --no-deps -D warnings` not run — **Applied with adaptation**

Run for the first time. Surfaced **two pre-existing doc errors in
untouched code** that the project's lint-deny policy had never actually
been enforced against:

1. [contamination_estimation.rs:14](../../../src/var_calling/contamination_estimation.rs#L14) — broken intra-doc link `[\`PerPositionMerger\`]` (no `PerPositionMerger` in scope). Fixed by qualifying to `[\`PerPositionMerger\`](crate::var_calling::per_position_merger::PerPositionMerger)`.
2. [per_group_merger.rs:80](../../../src/var_calling/per_group_merger.rs#L80) — unclosed HTML tag in `"<OTHER>"` literal. Fixed by backticking: `` `<OTHER>` ``.

**Full enforcement is blocked at the time of this commit** by a parallel
work-in-progress edit to `posterior_engine.rs` (the human PM's
in-flight refactor pulling a math-backend through the engine; 7+ E0061
type errors in the working tree, new `posterior_engine/backends.rs`
module untracked). Once that compiles, the full doc-build can be
re-verified — the dust_filter module's intra-doc links have been
visually walked and all targets resolve.

#### Mi22 — Module-doc breadcrumb about spec memory model — **Applied**

[Lines 43-48](../../../src/var_calling/dust_filter.rs#L43). New paragraph
notes the per-chromosome-batch design departs from the spec's stated
`O(w)` bound, points at the implementation plan's "Streaming model" for
the rationale.

#### Mi23 — `u32::try_from` on `seq.len()` in `sdust_mask` — **Applied as debug_assert**

[Line 416](../../../src/var_calling/dust_filter.rs#L416). Since `sdust_mask`
is now private (Q3 / Mi24 resolution), the cleanest enforcement is a
`debug_assert!(seq.len() <= u32::MAX as usize, ...)` — production callers
go through `DustFilter::ensure_mask_for`, which receives a `Vec<u8>` from
a fetcher bounded by `ParsedChromosome.length: u32`, so the precondition
is guaranteed at the upstream boundary.

#### Mi24 — `sdust_mask` panic vs `Result` dual validation — **Applied (option c)**

Q3 resolution. `sdust_mask` made non-pub; `SdustIntervals` type alias also
non-pub; the runtime `assert!(window >= SD_WLEN)` became
`debug_assert!(window >= SD_WLEN)` since the only caller is the
constructor-validated `DustFilter::ensure_mask_for`. Public surface is now
`DustFilter`, `DustFilterConfig`, `DustFilterError`, `DEFAULT_DUST_WINDOW`,
`DEFAULT_DUST_THRESHOLD` — every public type that exposes the window
parameter routes through `DustFilterConfig::new`'s validation.

### Nits

| Nit | Status | Note |
|---|---|---|
| Casts in hot loop lack `// safe:` comments | Applied | Covered by the M8 `// INVARIANT:` block in `find_perfect`. |
| `expect("deque non-empty at cap")` message wording | Won't fix | Kept the existing wording; the M7 `// PANIC-FREE:` comment carries the precise invariant. |
| Stale "cloned binary holds" comment | Applied | Tightened to "committed golden vectors hold". |
| `#[must_use]` on `DustFilter` | Won't fix | Project uses `Iterator` adapters without it elsewhere (`per_position_merger`, `posterior_engine`); not adopting the convention for one module would create noise without a behavior change. |
| `use std::iter::repeat_n` hoist | Applied | Hoisted (only kept in 4 spots, all fully qualified — fmt left them as-is and clippy doesn't object). |
| `PerfInterval`/`perf` abbreviation | Acceptable (no change) | C-traceability tradeoff documented at the type level. |
| Literal `10` density factor recurrence | Applied | Extracted via `dust_density_exceeds` (Mi17). |
| `StubFetcher::fail_for` naming | Won't fix | Private test helper, single-file scope; not worth the churn. |
| `PerfInterval::r` and `::l` single-letter fields | Acceptable (no change) | C-traceability; doc-comment on the type explains them. |

## New tests added

12 new tests beyond the original 25:

- `sdust_threshold_is_strictly_greater_at_density_boundary` (M5)
- `sdust_returns_empty_on_seq_shorter_than_wlen`
- `sdust_handles_seq_of_exactly_wlen_bases`
- `sdust_returns_empty_on_all_n`
- `sdust_handles_long_n_run_after_low_complexity`
- `sdust_at_minimum_window_does_not_panic`
- `sdust_threshold_zero_masks_aggressively`
- `sdust_threshold_max_does_not_overflow`
- `sdust_invariants_hold_on_random_seeded_inputs` (M2)
- `filter_half_open_boundary_pos_at_start_masked_pos_at_finish_passes` (M3)
- `filter_surfaces_invalid_pos_zero_and_latches` (M1)
- `filter_latches_after_upstream_exhaustion`
- `config_new_accepts_threshold_zero_and_max`

Renamed: `sdust_mask_panics_on_tiny_window` →
`sdust_mask_debug_asserts_on_tiny_window` (more accurate after Mi24/Q3).

## Validation results

Inside the dev container (`./scripts/dev.sh`):

| Command | Result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib var_calling::dust_filter` | 38 / 38 pass |
| `cargo test --all-targets --all-features` (at applied state, before the parallel WIP landed) | 706 lib + 109 integration tests pass |
| `cargo doc --no-deps --all-features` (with `RUSTDOCFLAGS="-D warnings"`) | Surfaced 2 pre-existing doc errors in untouched code (both fixed). Full re-run blocked by parallel `posterior_engine` WIP. |

The property test `sdust_invariants_hold_on_random_seeded_inputs`
exercises 32 seeds × 5 alphabets × 10 lengths × 5 windows × 5 thresholds
= 40,000 input combinations and runs in ~30s — fast enough for routine
`cargo test`. The two structural invariants it pins (non-empty intervals,
sorted-and-disjoint) are the load-bearing properties of `is_masked`.

## Tradeoffs and follow-ups

- **Mi21 cargo-doc enforcement** awaits the user's posterior_engine
  refactor compiling. Once that lands, re-run
  `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features` and
  capture the exit-0 line in a follow-up note.
- **No proptest dependency added.** The seeded-random invariant loop in
  `sdust_invariants_hold_on_random_seeded_inputs` is fixed-seed and
  doesn't shrink failures. If a future failure resists hand-narrowing,
  adding `proptest` as a dev-dep is a one-line change and the test
  layout already lines up with `proptest!` macro structure.
- **Two `lh3/sdust` faithful-port properties** uncovered by the property
  test: published interval `start` and `finish` can exceed `seq.len()`
  on N-heavy inputs; published intervals can physically span N bytes.
  Both verified against the cloned binary on the same fixtures. Both
  are harmless for the current consumer (`DustFilter::is_masked`). If
  the algorithm is later replaced by a stricter variant, both
  properties should be re-checked.
- **`DustFilter::new` still takes `Vec<ParsedChromosome>` by value.**
  Out-of-scope observation from the review (the type doesn't enforce
  that this is the same table the merger holds). Slots into the cohort
  CLI wiring PR.
- **Whole-chromosome reference fetch** memory transient is unchanged;
  documented in the impl report and the spec breadcrumb (Mi22). Revisit
  only if profiling on a real cohort shows it matters.
- **`find_perfect`'s `Vec::insert` is `O(W²)`** worst case. C source has
  the same shape (`memmove` inside a `for`). Intentional port fidelity;
  performance follow-up if benchmarking shows it dominates real-cohort
  wall time.
