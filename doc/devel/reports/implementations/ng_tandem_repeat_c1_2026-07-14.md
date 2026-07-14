# ng tandem-repeat scanner — Milestone C1 (region tiling, resident)

*Implementation report, 2026-07-14. Plan:
[`doc/devel/ng/impl_plan/ssr_repeat_scanner.md`](../../ng/impl_plan/ssr_repeat_scanner.md)
(Milestone C, step C1). Design: [spec](../../ng/spec/ssr_repeat_scanner.md) §3.6 +
[arch](../../ng/arch/ssr_repeat_scanner.md) §1.3, §2.2.*

## What landed

The **region seam's tiling core**, resident path, in `src/ng/tandem_repeat.rs`:

- **`tile(seq, periods, params, opts) -> Vec<Region>`** — the coverage-merge core: run
  `find_tandem_repeats`, sort by start, union overlapping/abutting intervals into merged repeat
  spans (keeping each span's constituent intervals), apply the `SegmentOptions` smoothing
  (`merge_gap` bridges short unique gaps; `min_repeat_len` drops short repeat blips back to
  unique), classify each surviving span `Repeat` vs `Satellite` by `max_repeat_len`, and emit
  `Unique` for every gap. The output tiles `[0, seq.len())` exactly, coordinate-ordered, with no
  two same-kind tiles adjacent.
- **`RegionScanner`** — the seam iterator (`Item = Result<Region, ScanError>`), with the
  **`over_slice`** constructor that precomputes `tile()` for a resident sequence. The resident
  path is infallible (always `Ok`); the windowed streaming constructor (C2) is a separate
  increment.

## Tests (6, spec §9)

A `triple` + `assert_tiles` helper checks the tiling contract (ordered, gap-free over `[0, n)`,
non-overlapping, no adjacent same-kind) on every case: empty input → nothing; a lone repeat →
`Unique, Repeat, Unique` with the Repeat covering the exact tract; overlapping periods 2/4/6 over
`(AT)*10` → **one** `Repeat` carrying all three intervals; the same tract with a 10 bp cap →
`Satellite`; `merge_gap` bridges two `(CAG)*5` tracts across a 12 bp gap (2 repeats → 1); and
`min_repeat_len` reclassifies a `(CAG)*3` blip to unique (the contig becomes one `Unique`). 21
module tests pass total.

## Deviation (recorded) — the sealed `ChromRefFetcher`

The arch (§2.2) sketched `RegionScanner<F: ChromRefFetcher>` with `over_slice` wrapping a slice as
a fetcher. `ChromRefFetcher` is a **sealed** trait ([`fasta/fetcher.rs`](../../../../src/fasta/fetcher.rs)),
so no in-memory implementation is possible — the sketch is infeasible. `over_slice` therefore
operates **directly on the resident slice** (no fetcher), and `RegionScanner` is a plain struct
(not generic over `F`). The **contract is unchanged** (repeat/satellite/unique tiling); only the
type shape adapts to the real API. The windowed streaming path (C2) is what actually needs a
fetcher — see the C2 decision surfaced at this checkpoint.

## Validation

`cargo fmt --check` clean; `cargo clippy --lib --tests -- -D warnings` clean on the module;
`cargo test --lib ng::tandem_repeat` → **21 pass**.

## Next

Milestone C2 — windowed streaming over a `ChromRefFetcher` (memory-bounded, window-count
invariant against `tile()` as the oracle). **Surfaced for a decision** given the sealed trait
(streaming is file-backed only) and that it has no live consumer yet.
