# ng tandem-repeat scanner — Milestone A (types + scaffold)

*Implementation report, 2026-07-14. Plan:
[`doc/devel/ng/impl_plan/ssr_repeat_scanner.md`](../../ng/impl_plan/ssr_repeat_scanner.md)
(Milestone A, steps A1–A3). Design: [spec](../../ng/spec/ssr_repeat_scanner.md) +
[arch](../../ng/arch/ssr_repeat_scanner.md).*

## What landed

A new module `src/ng/tandem_repeat.rs` carrying the scanner's **type vocabulary only** (no
detection logic — that is Milestone B/C), wired into `src/ng/mod.rs` as `pub mod tandem_repeat;`.

- **A1 — scaffold.** The file + its `#[cfg(test)]` block; module declaration in `ng/mod.rs`
  (one line + a doc-comment mention). One file, no folder — one algorithm, no bake-off.
- **A2 — input types.** `PeriodRange` (constrained: private `u8` fields, checked
  `new(min, max) -> Result<_, PeriodRangeError>` rejecting `min == 0` then `min > max`,
  `min()`/`max()` accessors); `PeriodRangeError` (`#[non_exhaustive]` `thiserror`:
  `ZeroMin` / `MinExceedsMax`); `ScanParams` + `SegmentOptions` (plain config, manual `Default`);
  and the named default consts (`DEFAULT_PERIODS = (1,6)`, `DEFAULT_MATCH_REWARD = 2`,
  `DEFAULT_MISMATCH_PENALTY = 7`, `DEFAULT_MIN_COPIES = 2`, `DEFAULT_MAX_REPEAT_LEN = 1000`,
  `DEFAULT_WINDOW_BP = 100_000`, `DEFAULT_MERGE_GAP = 0`, `DEFAULT_MIN_REPEAT_LEN = 0`).
- **A3 — output types.** `RepeatInterval` (`start`/`end` `u32`, `period` `u8`, `score` `i32`);
  `RegionSpan { start, end }`; `RepeatRegion { span, intervals: Box<[RepeatInterval]> }`; the
  three-kind `Region` enum (`Repeat` / `Satellite` / `Unique`); `ScanError` (`#[non_exhaustive]`,
  `Fetch { source: ChromRefFetchError }` — reusing `crate::fasta::ChromRefFetchError`).

## Tests

Six type-level unit tests (`cargo test --lib ng::tandem_repeat` — 6 pass): `PeriodRange::new`
accepts `(1,6)`/`(6,6)`/`(2,6)` and rejects `(0,_)` (→ `ZeroMin`, checked before the empty-range
rule) and `(3,2)` (→ `MinExceedsMax`); `DEFAULT_PERIODS == (1,6)` and is a valid range; the two
`Default`s equal their named-const values.

## Deviations (recorded)

- **Visibility `pub`, not the arch's illustrative `pub(crate)`.** Matches the sibling ng modules
  (`types.rs`, `read/filtering.rs`), which expose their vocabulary as reachable API — this also
  lets the Milestone-A types be defined ahead of their B/C consumers without tripping `dead_code`
  under `clippy -D warnings`. A documented note sits at the top of the module.
- **Extra named default consts** beyond the arch's `DEFAULT_PERIODS` — the individual scoring /
  segmentation defaults are pulled out as named `pub const`s (code-shape: no magic numbers), and
  the `Default` impls reference them.

## Validation

- `cargo fmt -- --check`: clean (after `cargo fmt`).
- `cargo test --lib ng::tandem_repeat`: 6 pass.
- `clippy --all-targets -- -D warnings`: **zero findings in `tandem_repeat.rs`.** The container's
  Rust 1.95.0 auto-bump surfaced 21 *pre-existing* `useless_vec`-style lints in unrelated files
  (`ssr/cohort/freebayes_emit.rs`, `ng/read/filtering.rs` tests, `examples/ssr_psp_seqdump.rs`);
  those are repo-baseline debt, not from this step, addressed in a separate `chore` commit
  (the WIP `examples/ssr_psp_seqdump.rs` left to its owner).

## Next

Milestone B — the `find_tandem_repeats` interval finder (lag-`p` scoring + Ruzzo–Tompa maximal
segments), the algorithmic heart, with its consumer-agnostic unit tests.
