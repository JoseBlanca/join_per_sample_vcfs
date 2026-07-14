# ng tandem-repeat scanner — Milestone A review

*Review report, 2026-07-14. Scope: the Milestone-A working-tree diff (`src/ng/tandem_repeat.rs`
new + `src/ng/mod.rs` one-line wire). A **focused self-review** proportionate to a types-only
milestone — no algorithmic surface (the detection logic is Milestone B/C), so the review
categories that bite are naming, defaults, errors, and module structure.*

## Verdict: approve, 0 Blocker / 0 Major

### Naming (clear)
Domain nouns throughout (`PeriodRange`, `ScanParams`, `SegmentOptions`, `RepeatInterval`,
`RegionSpan`, `RepeatRegion`, `Region`, `ScanError`); `PeriodRange::new` a verb; the default
`pub const`s are self-describing with units in their doc comments. No `utils`/`common`/`data`
placeholder names.

### Defaults (visible)
Every default is a named `pub const` with a doc comment giving its meaning and source (spec §3.3
for the 7/2 ratio, spec §3.6 for the 1 kb cap); the two `Default` impls reference them. No magic
numbers. `Option`-vs-sentinel is not in play here (the "off" smoothing knobs are genuinely `0`, a
documented no-op, not a disable sentinel standing in for absence).

### Errors (explicit, future-proof)
`PeriodRangeError` and `ScanError` are both `#[non_exhaustive]` `thiserror` enums with a per-variant
doc comment stating when they fire. The fallible surface is decided and minimal: only
`PeriodRange::new` is construction-fallible (a caller bug); detection stays infallible and
`ScanError` is reserved for the Milestone-C reference-read path. `ScanError` reuses
`crate::fasta::ChromRefFetchError` rather than minting a duplicate.

### Module structure (correct)
One file, no folder (no bake-off); `pub` items matching the ng-sibling convention (deviation from
the arch's `pub(crate)` recorded in the module doc and the impl report). `PeriodRange` is the one
constrained newtype — private fields behind a checked constructor; the unconstrained config structs
keep `pub` fields. `Eq`/`Hash`-safety: no floating-point fields anywhere (scores/weights are `i32`),
so the derived `PartialEq, Eq` are sound.

### Correctness (types only)
`PeriodRange::new` orders its checks so `(0, 0)` is `ZeroMin` (not an empty-range error) — asserted
in a test. `RepeatRegion.intervals` is `Box<[RepeatInterval]>` (owned, immutable slice — right for a
finished region). `RepeatInterval` is `Copy`; `RepeatRegion`/`Region` are `Clone` (not `Copy`, due
to the boxed slice) — appropriate.

## Nits (non-blocking, not applied)
- `DEFAULT_PERIODS` is a `(u8, u8)` tuple rather than a `const PeriodRange`, because `PeriodRange::new`
  is fallible/non-`const`. Acceptable; a `const`-constructor is not worth the ceremony for a
  two-field type. Left as-is.

## Findings applied
None — no actionable findings surfaced. Recorded as a no-op fixes stage.
