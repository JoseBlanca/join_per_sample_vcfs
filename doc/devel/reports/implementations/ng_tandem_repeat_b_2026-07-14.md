# ng tandem-repeat scanner — Milestone B (the interval finder)

*Implementation report, 2026-07-14. Plan:
[`doc/devel/ng/impl_plan/ssr_repeat_scanner.md`](../../ng/impl_plan/ssr_repeat_scanner.md)
(Milestone B, steps B1–B2). Design: [spec](../../ng/spec/ssr_repeat_scanner.md) §3.1–§3.5 +
[arch](../../ng/arch/ssr_repeat_scanner.md) §2.1.*

## What landed

The **algorithmic heart** — `find_tandem_repeats(seq, PeriodRange, &ScanParams) ->
Vec<RepeatInterval>` — in `src/ng/tandem_repeat.rs`, plus two internal helpers:

- **B1 — `canonical_base` + `maximal_scoring_subsequences`.** `canonical_base` upper-cases an
  ACGT byte or returns `None` (so `N == N` never matches). `maximal_scoring_subsequences` is the
  Ruzzo–Tompa (1999) linear pass over a score iterator, emitting every maximal scoring segment
  `(start, end, score)` in start order. It runs **online**: on the Ruzzo–Tompa "rule 2" case (a
  new positive element with no left-smaller stacked segment), every stacked segment is provably
  unabsorbable by any future element, so the stack is flushed and reset — keeping the working
  stack ≈ open-segment-count (not O(n)) on the sparse-positive signal the scanner produces.
- **B2 — `find_tandem_repeats`.** Per period `p in min..=max`: score position `j` (`+match_reward`
  on a canonical ACGT match to `seq[j-p]`, else `-mismatch_penalty`), run the RT pass, map each
  segment `[k0, k1]` to the tract `[k0, k1 + p + 1)`, and emit it when its implied copy count
  `(end - start) / p ≥ min_copies`. Raw overlapping intervals, no period de-dup (a consumer's job).

## Correctness anchor — a brute-force property test

Ruzzo–Tompa is subtle, so its correctness does **not** rest on inspection. A second, independent
`brute_segments` oracle computes maximal scoring subsequences straight from the definition
(O(n⁴): positive score, every proper sub-range strictly smaller, and maximal under containment),
and two tests assert the RT pass matches it: on 10 crafted cases (all-positive, all-negative,
bridged vs unbridged dips, the classic `[2,-7,2]` split) and on **400 deterministic-LCG
pseudo-random** ±-score arrays. *(Writing this test immediately caught a bug — in my first
`brute_segments`, not the RT: the prefix-sum characterisation I first used missed the
maximality-under-containment condition and over-reported sub-intervals. The RT pass was correct
from the start; the oracle was fixed.)*

## The remaining `find_tandem_repeats` tests (consumer-agnostic, spec §9)

Clean `(CAG)*8` → one exact interval `{0, 24, 3, 42}`; a substitution-interrupted tract → **one**
span; a **single-indel**-interrupted tract → **one** span (the lag-`p` locality property, §3.4);
an `N`-run → nothing; a soft-masked (lowercase) tract → found; an `(AT)*10` tract → emitted at
periods 2, 4, **and** 6 (no de-dup) but not 3/5; and the `min_copies` floor drops a lone
coincidental match.

## Deviations / notes (recorded)

- **Online rule-2 flush** is an addition beyond a textbook offline Ruzzo–Tompa — a memory
  optimisation (bounded working stack), proven-equivalent and pinned by the property test.
- **Integer widths.** The running score/`l`/`r` are `i64` (no overflow on any contig); the emitted
  `RepeatInterval.score` is `i32::try_from(...).unwrap_or(i32::MAX)` — saturating, since a segment
  that large is a giant tract the region seam caps as satellite anyway. Coordinates are `u32`
  (project genomic-coordinate convention; contigs < 4 Gbp).
- **Time/memory.** O(n) per period in the common case; the RT search is `rposition` from the stack
  top, so worst-case (adversarial, not real DNA) is O(n·k). The region seam (Milestone C) windows
  the walk, bounding both regardless.

## Validation

`cargo fmt --check` clean; `cargo clippy --lib --tests -- -D warnings` clean on the module;
`cargo test --lib ng::tandem_repeat` → **15 pass**.

## Next

Milestone C — the `RegionScanner` region seam (coverage merge → repeat/satellite/unique tiling
over a resident slice, then windowed streaming over `ChromRefFetcher`).
