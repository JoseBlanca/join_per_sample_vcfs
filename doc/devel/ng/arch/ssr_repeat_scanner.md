# Tandem-repeat scanner: types & interfaces

*Status: architecture draft (2026-07-14), companion to the spec
[`../spec/ssr_repeat_scanner.md`](../spec/ssr_repeat_scanner.md) (the design and its rationale).
Sits alongside the shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary +
step traits) and [`module_layout.md`](module_layout.md), but the scanner is **not an ng step** —
it is a standalone sequence primitive the ng caller (and the STR catalog) consume. Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): domain nouns for
types, verbs for functions, newtypes for constrained scalars; **STR** in prose ↔ `ssr`/repeat
terms in code. Signatures are illustrative; the **contract** is the deliverable. See the spec for
the "why" behind every decision.*

*Doc-vs-code location: this doc lives in `ng/arch/` (its motivating consumer is the ng caller and
its spec is in `ng/spec/`), and the code lives in the same tree at `src/ng/tandem_repeat.rs`.*

## Module home

One file, `src/ng/tandem_repeat.rs`, with its `#[cfg(test)]` block beside it; promoted to a
`src/ng/tandem_repeat/` directory only if it outgrows one file. It is a **shared sequence
primitive in the ng tree** — the home of the caller generation that primarily drives it — **not**
one of the ng pipeline steps and **not** under `src/ssr/` (the STR caller). It has **no
swappable-trait bake-off**: there is one detection algorithm, so it is plain functions/structs,
not a `dyn` trait family (code-shape: a step with no competitors is a file, not a trait). It stays
policy-free; the catalog consumes it cross-module from `src/ssr/catalog/`.

## 1. Types

### 1.1 Inputs — the scope and scoring knobs

The only configuration surface: which periods to scan, how to score, and (for the region seam)
how to shape the partition and walk. `PeriodRange` is the one **constrained** newtype (hidden
field, checked constructor — period 0 is meaningless and `min > max` is empty); `ScanParams` and
`SegmentOptions` are unconstrained config structs with visible `Default`s (defaults are the
scanner's general starting values, tuned empirically — spec §3.3, §9).

```rust
/// The inclusive period range (motif lengths, in bp) to scan. Constrained: `1 <= min <= max`.
/// The algorithm imposes no ceiling of its own — a consumer passes the range it wants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PeriodRange { min: u8, max: u8 }   // field private; construct via `new`

impl PeriodRange {
    /// The only constructor. Rejects `min == 0` and `min > max` loudly (a consumer bug).
    pub(crate) fn new(min: u8, max: u8) -> Result<Self, PeriodRangeError>;
    pub(crate) fn min(self) -> u8;
    pub(crate) fn max(self) -> u8;
}

/// The scanner's general default period range: 1..=6. Each consumer overrides as it sees fit
/// (the catalog / ng snp-str caller uses 2..=6). A named const, not a magic literal.
pub(crate) const DEFAULT_PERIODS: (u8, u8) = (1, 6);

/// Lag-p self-comparison scoring plus the minimum-copies emission floor. The match/mismatch
/// ratio is the mismatch tolerance (a ratio `r` holds tracts to purity `r/(1+r)`; spec §3.3).
/// `Default` = { 2, 7, 2 } — a starting point, tuned from experiments (spec §9), not a commitment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ScanParams {
    pub match_reward: i32,       // > 0
    pub mismatch_penalty: i32,   // > 0
    pub min_copies: u32,         // emission floor; a segment shorter than min_copies·period is dropped
}

/// Region-seam shaping + walk knobs — separate from `ScanParams` because they shape the
/// partition and the streaming, not the detection. `Default` = { 1000, 100_000, 0, 0 }.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SegmentOptions {
    pub max_repeat_len: u32,   // above this a repeat region is Satellite, not a genotypeable STR (1 kb, spec §3.6)
    pub window_bp: u32,        // streaming window core; a memory knob, region-count-invariant
    pub merge_gap: u32,        // smoothing: bridge unique gaps shorter than this (0 = off)
    pub min_repeat_len: u32,   // smoothing: reclassify sub-threshold repeat blips as unique (0 = off)
}

/// Construction-time error for `PeriodRange`. Fires only on a caller bug, never on data.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub(crate) enum PeriodRangeError {
    #[error("period range min must be >= 1 (period 0 is meaningless)")]
    ZeroMin,
    #[error("period range min ({min}) exceeds max ({max})")]
    MinExceedsMax { min: u8, max: u8 },
}
```

### 1.2 Output — the interval finder's natural shape

```rust
/// One detected tandem-repeat interval. Coordinates are 0-based half-open (`[start, end)`);
/// `period` is the lag it was found at; `score` is the Ruzzo–Tompa segment total. Nothing
/// consumer-specific here — it happens to carry exactly the fields the STR post-filter reads
/// (spec §4), which is convenience, not coupling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct RepeatInterval {
    pub start: u32,
    pub end: u32,
    pub period: u8,
    pub score: i32,
}
```

### 1.3 Output — the region seam's tiling

The seam yields an ordered, gap-free, three-kind tiling of the scanned span. `Repeat` carries its
constituent intervals (so the caller reads the repeat structure without re-scanning); the other
two carry only a `RegionSpan`.

```rust
/// A half-open coordinate pair, 0-based. The plain span `Satellite`/`Unique` carry and the
/// union `RepeatRegion.span` holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct RegionSpan { pub start: u32, pub end: u32 }

/// A genotypeable tandem repeat: merged repeat coverage no longer than `max_repeat_len`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RepeatRegion {
    pub span: RegionSpan,                       // union of `intervals`' spans
    pub intervals: Box<[RepeatInterval]>, // >= 1, coordinate-ordered; the overlaps merged here
}

/// One tile of the reference. Consecutive tiles never share a kind (each is maximal), are
/// coordinate-ordered and non-overlapping, and together cover the scanned span exactly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Region {
    Repeat(RepeatRegion),   // -> the STR path
    Satellite(RegionSpan),        // repeat coverage > max_repeat_len (satellite DNA): mask/skip, not STR, not SNP
    Unique(RegionSpan),           // no repeat coverage -> the SNP path
}
```

## 2. Interfaces

Two entry points over one core: the **interval finder** (a pure function) and the **region seam**
(a streaming iterator built on it).

### 2.1 The interval finder

```rust
/// Find every tandem-repeat interval in `seq` whose period lies in `periods`. Pure and total:
/// comparison is case-insensitive (bases upper-cased) and any non-ACGT base never matches, so it
/// is defined on arbitrary bytes. Returns raw, possibly-overlapping intervals (one region can
/// match at several periods) — the finder does NOT de-duplicate periods or resolve overlaps;
/// that is a consumer's job.
pub(crate) fn find_tandem_repeats(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
) -> Vec<RepeatInterval>;
```

**Contract.** For each `p in periods.min..=periods.max`, independently: score position `j` (`j >=
p`) `+match_reward` when `upper(seq[j]) == upper(seq[j-p])` and both are ACGT, else
`-mismatch_penalty`; take every Ruzzo–Tompa maximal scoring segment `[j0, j1]`; emit
`RepeatInterval { start: j0 - p, end: j1 + 1, period: p, score }` when the implied copy count
`(end - start) / p >= params.min_copies`. Substitutions and single indels stay inside one segment
(the lag-p comparison is local; spec §3.4); results across periods are concatenated. Deterministic
— a pure function of the inputs.

### 2.2 The region seam

A struct implementing `Iterator<Item = …>`, matching the repo's iterator-as-seam convention (the
ng `ReadFilter`). It streams a contig in windows through the existing sliding-buffer reference
reader, so peak memory is `~window_bp + max_repeat_len`, not the contig length (spec §3.6).

```rust
/// Streams a repeat/satellite/unique tiling of one contig, windowed and memory-bounded. Built on
/// `find_tandem_repeats` + a coverage merge, over the sliding-buffer `ChromRefFetcher`
/// (`src/fasta/fetcher.rs`) — the same ~1 MB-resident, monotonic-forward reader the cohort/DUST
/// path uses. A repeat spanning a window cut is emitted exactly once, by the window whose core
/// holds its start (the per-chunk-DUST halo/attribution pattern); the `max_repeat_len` cap bounds
/// the halo, so an STR-sized repeat is never split.
pub(crate) struct RegionScanner<F: ChromRefFetcher> { /* fetcher, contig_len, periods, params, opts, cursor */ }

impl<F: ChromRefFetcher> RegionScanner<F> {
    /// Stream `contig_len` bases from `fetcher` (a single-contig reader positioned at its start).
    pub(crate) fn new(
        fetcher: F,
        contig_len: u32,
        periods: PeriodRange,
        params: ScanParams,
        opts: SegmentOptions,
    ) -> Self;

    /// Convenience for a small, already-resident sequence (unit tests, in-memory callers):
    /// wraps `seq` in a trivial in-memory reader. Same iteration contract.
    pub(crate) fn over_slice(
        seq: &[u8],
        periods: PeriodRange,
        params: ScanParams,
        opts: SegmentOptions,
    ) -> RegionScanner<impl ChromRefFetcher>;
}

impl<F: ChromRefFetcher> Iterator for RegionScanner<F> {
    /// `Err` is terminal — a mid-scan reference-read failure (corrupt/truncated file) fuses the
    /// iterator (spec §8: the reference read is the seam's core, so an I/O error is fatal to it).
    type Item = Result<Region, ScanError>;
    fn next(&mut self) -> Option<Result<Region, ScanError>>;
}

/// Region-seam iteration error. Only a reference-read failure surfaces here; detection itself is
/// infallible.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub(crate) enum ScanError {
    #[error("reading the reference failed mid-scan")]
    Fetch { #[source] source: ChromRefFetchError },
}
```

**Contract.** The `Ok` items are coordinate-ordered, pairwise non-overlapping, of alternating
kind (two of the same kind never abut), and their union is exactly `[0, contig_len)`; every
`Repeat.intervals` is non-empty and its `span` is their union; a merged repeat span longer than
`opts.max_repeat_len` is `Satellite` (coalesced across windows into one span); `merge_gap` /
`min_repeat_len` apply the optional smoothing; an empty input yields nothing. The tiling is
**window-count-invariant** — a pure function of `(contig, periods, params, {max_repeat_len,
merge_gap, min_repeat_len})`; `window_bp` changes only memory and I/O granularity, never the
regions (the regression gate, §Test shape). Access to `fetcher` is monotonically non-decreasing,
satisfying its sliding-buffer contract by construction.

## 3. Design decisions — decided

Distilled from the spec; see it for the reasoning. Open items carry `OPEN:`.

- **Lives in the ng tree (`src/ng/tandem_repeat.rs`), not `src/ssr` — decided.** A policy-free
  shared primitive placed in the ng tree (the caller generation that primarily drives it), not
  under the STR caller. Rejected `src/ssr/tandem_repeat.rs` (misrepresents scope). Supersedes an
  earlier top-level `src/tandem_repeat.rs` lean (spec §4).
- **Two interfaces over one core — decided.** A low-level `find_tandem_repeats` (raw overlapping
  intervals) for consumers that resolve overlaps themselves (the catalog), and a high-level
  `RegionScanner` (resolved tiling) for the caller's router. Rejected a single interface: the two
  consumers want opposite shapes (spec §1, §3, §3.6).
- **Detection = lag-p self-comparison + Ruzzo–Tompa maximal segments — decided.** Clean-room from
  Benson's idea (compare to a p-shifted copy; runs of matches) built on the independent Ruzzo–Tompa
  algorithm; brute-forcing each period drops TRF's k-tuple heuristic and significance model. Not
  derived from the AGPL TRF-mod source (spec §2, clean-room note).
- **Impure-repeat tolerance is intrinsic, not a pass — decided.** The local lag-p comparison keeps
  substitutions (2-mismatch blip) and single indels (~p-mismatch burst) inside one interval; only
  large/clustered indels split. No wraparound DP, no stitch pass in v1 (spec §3.4).
- **Satellite cap, default 1 kb — decided.** Merged repeat coverage longer than `max_repeat_len` is
  `Satellite`, excluded from STR analysis (we do not genotype satellite DNA); the cap also bounds
  the streaming halo. Rejected: no cap (would treat satellites as STRs and unbound the halo)
  (spec §3.6).
- **Three-kind `Region` (Repeat / Satellite / Unique) — decided.** The router needs a three-way
  route; folding satellite into `Unique` would make the caller SNP-call satellite DNA, folding it
  into `Repeat` would feed satellites to STR genotyping. Both wrong (spec §3.6).
- **Windowed streaming via `ChromRefFetcher`, in scope — decided.** The seam reuses the existing
  sliding-buffer reader (not a new reader) and the per-chunk-DUST halo/attribution pattern, so it is
  memory-bounded over a whole contig and window-count-invariant. Rejected: hold the contig resident
  (spec §3.6; the memory-efficiency thesis).
- **Region-seam error model = terminal `Result` item — decided.** `Iterator::Item =
  Result<Region, ScanError>`; a mid-scan reference-read error fuses the iterator. This diverges
  from read filtering's `Item = MappedRead` (fatal-not-per-item) *because here the fallible
  reference read is the iteration's core*, not a side filter — so an I/O error must reach the
  caller, not be swallowed (spec §8). `find_tandem_repeats` and detection stay infallible; only
  `PeriodRange::new` is construction-fallible.
- **`PeriodRange` is the one constrained newtype — decided.** `1 <= min <= max`, checked
  constructor, `PeriodRangeError` on a caller bug. `ScanParams`/`SegmentOptions` are plain config
  with `Default`s and documented positivity invariants (code-shape: constrained scalars hide the
  field; unconstrained config does not) (spec §5). OPEN(impl): whether to validate
  `ScanParams`/`SegmentOptions` positivity in a constructor or document-only — a boy-scout call.

### Catalog-consumer decisions (spec §6)

- **Post-filter consumes `RepeatInterval` directly; `TrfRecord` deleted — decided.** `build_loci`
  reads only `start/end/period/score`, exactly `RepeatInterval`'s fields, so the change is a
  mechanical type substitution with the post-filter *logic* untouched (parity by construction). The
  never-read `frac_match`/`pattern` vanish (spec §6). One knock-on from the `period: u8` choice
  (§1.1): `TrfRecord.period` was `u16`, so `postprocess`'s `MIN_PERIOD`/`MAX_PERIOD` constants
  narrow `u16 → u8` too — still mechanical, values unchanged, no cast left at the comparison.
- **Catalog uses the raw finder, no satellite cap — decided.** It needs per-interval periods for
  its compound/bundle logic and must reproduce the current `trf-mod`→postprocess loci, so it calls
  `find_tandem_repeats` (not `RegionScanner`) and applies **no** `max_repeat_len` cap — adding one
  would diverge from the golden (spec §6).
- **Header `trf_mod_version` → `detector` — decided.** A built-in-detector provenance field
  (name + scan weights) replaces the external-binary version; small, localised `io.rs` change,
  outside `postprocess.rs` (spec §6).

## 4. Reconciliation with existing code

Verify against the code when implementing; consolidation points, not new duplicates.

| name | existing code | action |
|---|---|---|
| `RepeatInterval` | `TrfRecord` ([`src/ssr/catalog/trf.rs`](../../../../src/ssr/catalog/trf.rs)) | **replaces** it — same read fields, minus never-read `frac_match`/`pattern`; `trf.rs` deleted (spec §7) |
| `find_tandem_repeats` | `trf::run_on_contig` (subprocess + BED parse) | **new**; slots into `catalog::run` where `run_on_contig` was |
| `RegionScanner` / `Region` / `RepeatRegion` / `RegionSpan` | — | **new**; the router seam, no prior analogue |
| `ChromRefFetcher` / `StreamingChromRefFetcher` | [`src/fasta/fetcher.rs`](../../../../src/fasta/fetcher.rs) | **reuse as-is** — the sliding-buffer windowed reader `RegionScanner` streams over |
| window halo + split-invariance | cohort per-chunk DUST (`src/var_calling/dust_filter.rs` + cohort integration) | **pattern reuse** — halo/attribution + the window-count-invariance gate |
| `postprocess::build_loci` (+ `drop_bundles`/`is_close`) | [`src/ssr/catalog/postprocess.rs`](../../../../src/ssr/catalog/postprocess.rs) | **retype** input `TrfRecord` → `RepeatInterval`; logic unchanged |
| `CatalogHeader.detector` | `CatalogHeader.trf_mod_version` ([`src/ssr/catalog/io.rs`](../../../../src/ssr/catalog/io.rs)) | rename/repurpose the header field + its test constructors |
| `CatalogError::Trf*`, `--trf-mod-path`/`--temp-dir`, Containerfile install | `catalog/mod.rs`, `ssr_catalog.rs`, `Containerfile` | **delete** after parity banks (spec §7) |
| parity-fixture pattern | `src/baq/` + `baq_tests.rs` | template for the golden-catalog parity test |

## 5. Open items

- **Scoring weights (`ScanParams::Default`)** — start `2 / 7`; the final values are tuned from the
  parity fixture and downstream benchmarks (spec §3.3, §9). An empirical confirmation, not a design
  fork.
- **Per-consumer scan range** — module default `1..=6`; catalog / ng caller passes `2..=6`. A
  consumer-side choice that may move after experiments (spec §9). Not a scanner decision.
- **`window_bp` default** — `~100_000` as a memory/I/O granularity knob; the exact value is a
  perf-tuning confirmation (region-count-invariant, so it cannot affect correctness).
- **OPEN(impl): `ScanParams`/`SegmentOptions` validation** — constructor-checked vs document-only
  positivity (boy-scout, §3).
- **`ChromRefFetcher` concrete surface** — the exact `fetch`/`iter_bases` calls and the in-memory
  fetcher for `over_slice` resolve against `fetcher.rs` at implementation time. Not a design
  decision.

## 6. Test & bench shape (spec §9)

- **Algorithm unit tests** beside `tandem_repeat.rs`, consumer-agnostic. Interval finder: a clean
  tract (exact `start/end/period/score`); substitution-interrupted (one segment); single-indel
  (one span — the impure-tolerance property); large/multi-indel (the pathological split); `N`-run
  (nothing); soft-masked tract (found); a period-2 tract that also matches at p = 4/6 (all emitted —
  no period de-dup); `PeriodRange` validation. Region seam: the tiling contract (ordered,
  non-overlapping, union == `[0, n)`); overlapping different-period intervals merge into one
  `Repeat` carrying both; `Unique, Repeat, Unique` around a lone repeat; a > `max_repeat_len` run →
  `Satellite` (coalescing across windows); **window-count invariance** (identical regions across two
  `window_bp` settings — the boundary-halo correctness gate, mirroring the cohort DUST
  split-invariance test); `merge_gap`/`min_repeat_len` smoothing; empty input.
- **Catalog consumer** gets the golden-catalog **parity** test (snapshot the current
  `trf-mod`→postprocess `Locus` set on a small tomato-reference subset; assert
  `find_tandem_repeats`→`build_loci` reproduces it at ≥ 99 % recall with a reviewed diff) and the
  **end-to-end** `benchmarks/ssr_tomato1/` genotyping/HipSTR-concordance check (spec §6).
- **No `bench/`**: one detection algorithm, no bake-off frontier to plot; the regression anchors are
  the parity + window-invariance tests.
