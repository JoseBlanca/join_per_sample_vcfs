# ng tandem-repeat scanner — Milestone B review

*Review report, 2026-07-14. Scope: the Milestone-B diff in `src/ng/tandem_repeat.rs`
(`canonical_base`, `maximal_scoring_subsequences`, `find_tandem_repeats`, and their tests). A
focused review of the algorithmic heart; the north-star correctness check is the brute-force
property test, which is stronger than inspection for the Ruzzo–Tompa pass.*

## Verdict: approve, 0 Blocker / 0 Major

### Correctness (the heart)
- **Ruzzo–Tompa pass** is validated against an independent O(n⁴) definitional oracle over crafted
  + 400 pseudo-random inputs — the strongest available guarantee. The online rule-2 flush (the one
  non-textbook part) is covered by the same test, so its equivalence to offline RT is proven, not
  asserted.
- **Coordinate mapping** `[k0, k1] → [k0, k1 + p + 1)` is pinned by an *exact* interval assertion on
  a clean `(CAG)*8` (`{0, 24, 3, 42}`), the case where an off-by-one would show immediately.
- **Match predicate** correctly excludes `N`/non-ACGT (an `N`-run emits nothing) and is
  case-insensitive (soft-masked tract found).
- **Impure-repeat behaviour** matches the spec §3.4 claims by test: substitution and single indel
  each stay one interval.

### Overflow / widths
`i64` running score (no contig overflows it); saturating `i32` cast on emission (documented, and a
segment that large is satellite-capped downstream); `u32` coordinates per the project convention.
The `k1 + p + 1` and `(p..n)` arithmetic cannot underflow (`j >= p` guarantees `j - p >= 0`; the
`p >= n` guard skips periods with no comparable position).

### Complexity / memory
O(n) per period common-case; the `rposition` stack search is O(k) per positive element, so
worst-case O(n·k) on adversarial input. The rule-2 flush keeps `k` ≈ open segments on real DNA, and
the Milestone-C window bounds it regardless. Acceptable and documented; not a blocker.

### Naming / structure
`maximal_scoring_subsequences`, `canonical_base`, `RtSeg` are clear; the RT rules are commented by
number. `find_tandem_repeats` is `pub`, the helpers private. No magic numbers (weights come from
`ScanParams`).

## Nits (non-blocking, not applied)
- The O(n·k) worst case could be made strict O(n) with a monotone-stack invariant, but the flush +
  window make it a non-issue in practice; deferred rather than complicate the proven code.

## Findings applied
None actionable. Recorded as a no-op fixes stage. (One bug *was* found and fixed during
implementation — in the test oracle, not the production code — see the impl report.)
