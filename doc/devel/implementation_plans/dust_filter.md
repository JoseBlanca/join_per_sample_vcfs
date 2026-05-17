# Stage 3 — sdust low-complexity filter

Proposal date: 2026-05-17.

## Domain intent

A streaming reference-only filter that sits between the multi-way
per-position iterator (Stage 2 → Stage 4 seam) and the variant grouper.
It silently drops `PerPositionPileups` items whose reference position is
flagged low-complexity by the symmetric DUST algorithm (sdust), so no
`OverlappingVarGroup` is ever constructed for those positions and no
merged record is emitted.

```
.psp_i ─┐
.psp_j ─┼─► PerPositionMerger ─► DustFilter ─► VariantGrouper ─► ...
.psp_k ─┘                            ▲
                                     │
                          reference FASTA ──┘
```

Spec — algorithm and rationale:
[calling_pipeline_architecture.md §"Stage 3 — low-complexity filter"](../specs/calling_pipeline_architecture.md#L875).
Upstream iterator that the filter wraps:
[src/var_calling/per_position_merger.rs](../../src/var_calling/per_position_merger.rs).
Existing reference fetcher we will reuse:
[src/per_sample_pileup/ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs).
**Reference implementation we port from:** [`lh3/sdust`](https://github.com/lh3/sdust),
vendored locally at [sdust/](../../sdust/) (gitignored). The Rust
algorithm in this plan is a near-direct port of `sdust_core` from
[sdust/sdust.c](../../sdust/sdust.c); when in doubt, the C source is
authoritative.

## Why now

Stages 1, 2, 4, 5 and the Stage 6 engine + contamination side-pass are
all shipped (see [PROJECT_STATUS.md](../../PROJECT_STATUS.md)). Stage 3
is the last missing piece between the multi-way merger and the variant
grouper — without it the per-position stream feeding Stage 4 lacks the
low-complexity drop the spec mandates as the *single most effective
simple filter against pipeline false positives*. Once it lands, the
remaining work to a usable multi-sample SNP caller CLI is wiring, not
new components.

## Algorithm — symmetric DUST (sdust), as ported from `lh3/sdust`

The algorithm walks the reference one base at a time, encoding three
consecutive ACGT bases as a 6-bit triplet word, and maintains a list
of **perfect intervals** — currently-tracked maximal subintervals
(inside a window of up to `W` bases) whose score density `s / L`
exceeds `T / 10`, where `L` is the subinterval's triplet count and
`s = Σ_t f_t · (f_t − 1) / 2` is the standard triplet-pair score.

Default `W = 64`, `T = 20` (lh3/sdust + minimap2 defaults). With
`T = 20`, the density threshold is `> 2.0` triplet-pair counts per
triplet.

### State (direct port of [sdust.c:41-46](../../sdust/sdust.c#L41))

```text
window           triplet deque of length ≤ W − SD_WLEN + 1 = W − 2
cw : [u32; 64]   triplet counts across the whole window
cv : [u32; 64]   triplet counts in the active suffix of length L
rw : u32         score `Σ cw[t]·(cw[t]−1)/2` over the whole window, kept incrementally
rv : u32         score over the active suffix, kept incrementally
L  : u32         length (in triplets) of the active suffix
P  : Vec<PerfInterval>   candidate maximal subintervals, sorted by
                          descending start then ascending finish
res: Vec<(u32, u32)>     emitted masked intervals (BED-style 0-based half-open)
```

`PerfInterval` carries `{ start: u32, finish: u32, r: u32, l: u32 }`
([sdust.c:12-15](../../sdust/sdust.c#L12)).

### Per-base step (port of `sdust_core` loop body — [sdust.c:140-156](../../sdust/sdust.c#L140))

For each input byte at index `i`, map to `b ∈ {0,1,2,3,4}` via the
ACGT-or-N table ([sdust.c:22-39](../../sdust/sdust.c#L22)):

1. **`b < 4` (an A/C/G/T base).** Shift `b` into the rolling 6-bit word
   `t = ((t << 2) | b) & 0x3F`. Increment the "consecutive ACGT" counter
   `l`. Once `l ≥ 3` we have a complete triplet, and:
   - Compute the new window's left edge `start = max(l − W, 0) + (i + 1 − l)`
     (in *base* coordinates relative to the start of the current run of
     ACGT bases).
   - `save_masked_regions(start)`: drop any `PerfInterval` whose `start
     < start_of_new_window` from the tail of `P` into `res`, merging
     adjacent or overlapping intervals on the way
     ([sdust.c:87-101](../../sdust/sdust.c#L87)).
   - `shift_window(t)`: push `t` into the deque; if the deque was full,
     pop the leftmost triplet `s` and update `cw[s]`, `rw`, optionally
     `cv[s]`/`rv`/`L`. Then `cw[t]++; rw += cw[t]; cv[t]++; rv += cv[t]`.
     Finally — and this is the *suffix-trim invariant* — while
     `10 · cv[t] > 2T`, repeatedly drop the leftmost triplet `s` of the
     suffix, decrement `cv[s]`, subtract it from `rv`, and decrement `L`
     ([sdust.c:65-85](../../sdust/sdust.c#L65)).
   - If the whole-window density exceeds the threshold
     (`10 · rw > L · T`), call `find_perfect` to insert any new
     maximal subintervals discovered by the latest base extension
     ([sdust.c:148](../../sdust/sdust.c#L148)).

2. **`b = 4` (N, or end-of-sequence sentinel).** Flush every remaining
   `PerfInterval` from `P` into `res`, then reset `l = 0; t = 0`. N
   effectively breaks the reference into independent pieces — perfect
   intervals never span an N
   ([sdust.c:151-155](../../sdust/sdust.c#L151)).

### `find_perfect` (port of [sdust.c:103-127](../../sdust/sdust.c#L103))

Scans backwards through the window starting from the suffix boundary,
maintaining a *local* copy of `cv` as it grows the candidate
subinterval to the left. For each candidate length, if its density
exceeds threshold *and* it dominates any existing perfect interval
overlapping the same region (per the `new_r * max_l ≥ max_r * new_l`
cross-multiplied comparison at [sdust.c:117](../../sdust/sdust.c#L117)),
inserts it at the right rank in `P`. The suffix-trim invariant
bounds the scan distance, which is why `shift_window` enforces it
eagerly.

### Output

A sorted, non-overlapping list of 0-based, half-open base intervals.
Coordinates are *base* coordinates, not triplet coordinates: an
interval `(start, finish)` masks bases `seq[start..finish]`. The
"+ (SD_WLEN - 1)" adjustment at [sdust.c:122](../../sdust/sdust.c#L122)
converts the triplet-count length back to a base-count finish.

### N and lowercase

`seq_nt4_table` ([sdust.c:22-39](../../sdust/sdust.c#L22)) maps:
- `A/a → 0`, `C/c → 1`, `G/g → 2`, `T/t → 3`
- everything else → `4`

So lowercase ACGT is treated identically to uppercase, and N (along
with any other non-ACGT character) breaks the triplet stream. We port
this table verbatim; the spec already documents the N-handling
convention.

## Streaming model: per-chromosome batch + sweep lookup

We do **not** integrate sdust's per-base state machine with the
upstream `(chrom_id, pos)` stream. The upstream's "skip positions no
sample covered" access pattern makes a streaming-merged implementation
substantially more complex than the canonical batch sdust, and the
benefit is tiny: the masked-interval list for a single chromosome is
empirically a small fraction of its length (low-complexity regions are
typically < 1% of a WGS reference).

Instead:

1. On the first emission of a new `chrom_id`, fetch the whole
   chromosome's reference via
   `RefSeqFetcher::fetch(chrom_id, 1, chrom_length)`.
2. Run `sdust_mask(reference, window, threshold)` — a pure Rust
   function that ports `sdust_core` — to produce the
   `Vec<(u32, u32)>` masked-interval list for that chromosome.
3. Maintain a sweep pointer `sweep: usize` into the masked-interval
   list. Upstream positions are monotonic per chromosome, so the
   pointer advances monotonically: for each upstream `pos`, advance
   `sweep` past every interval whose `finish ≤ pos0` (0-based), then
   check whether `current_mask[sweep]` covers `pos0`. O(1) amortised
   per position.
4. When `chrom_id` changes, drop the previous chromosome's interval
   list, reset `sweep = 0`, and repeat.

Memory: bounded by the *current* chromosome's masked-interval list,
which in practice is tens of MB for human chr1 (~250 Mbp, ~1% masked,
~2.5 million bases, ~30 K merged intervals × 8 bytes each ≈ 240 KB).
This is the worst real-world case; for plant / animal references with
less repeat content the list is smaller. The previous chromosome's
list is dropped before the next is loaded, so the steady-state
footprint is one chromosome's worth.

The reference-base buffer (one full chromosome of ACGT bytes) is
also held during sdust execution, then dropped — for human chr1 that
is ~250 MB transient. This is the same shape `ChromBoundaryRefFetcher`
already holds for its noodles-cache; per the existing fetcher's
docstring, steady-state one-chrom-resident is the design contract.

## Out of scope

- **Modifying `.psp` files.** Stage 1 output is untouched; the filter
  is a cohort-level policy applied at iterator time, not a property
  of the per-sample artefact. Spec line 922-926.
- **Soft masking / a "low-confidence" tier.** A position is either
  dropped or passed through. No third state.
- **Classical (asymmetric, BLAST) DUST.** The spec moved to sdust as
  of 2026-05-17; the classical algorithm is referenced for historical
  context only.
- **Region-specific overrides** (per-chromosome thresholds, BED-
  region exemptions). The standing-item BED-region skip on the Stage
  1 CLI ([PROJECT_STATUS.md](../../PROJECT_STATUS.md)) is a separate
  feature.
- **Persisting a mask file.** The filter is purely a transformation
  on the iterator stream. No intermediate artefact is written.
- **Reference loading.** The filter does not open the FASTA; the
  caller hands it a fetcher.
- **A streaming sdust (state machine integrated with upstream
  positions).** Per-chromosome batch is simpler and matches the
  canonical reference; see "Streaming model" above for the trade.

## API shape

New module: `src/var_calling/dust_filter.rs`, declared from
[src/var_calling/mod.rs](../../src/var_calling/mod.rs) alongside
`per_position_merger`, `variant_grouping`, `per_group_merger`,
`posterior_engine`, `contamination_estimation`.

### Layer 1 — pure `sdust_mask`

The algorithmic core, with no I/O and no upstream-iterator
entanglement. Tested directly against the cloned `lh3/sdust` binary.

```rust
// src/var_calling/dust_filter.rs

/// Sorted, non-overlapping masked intervals on a single chromosome.
/// `(start, finish)` is BED-style 0-based half-open: bases at indices
/// `start..finish` of the input slice are low-complexity.
pub type SdustIntervals = Vec<(u32, u32)>;

/// Run sdust over `seq` (an ACGT/N byte slice) and return the
/// sorted, non-overlapping list of low-complexity intervals.
///
/// Direct port of `sdust_core` from [`lh3/sdust`](../../sdust/sdust.c).
/// `window` maps to that implementation's `W`; `threshold` maps to
/// `T`. The masking condition is `10·score > T·triplet_count`
/// (strictly greater), as at [sdust.c:111](../../sdust/sdust.c#L111).
///
/// `window` must be ≥ 3 (the minimum to fit one triplet). Validate
/// at the [`DustFilterConfig::new`] call site; this function takes
/// the validated config as a precondition.
pub fn sdust_mask(seq: &[u8], window: u32, threshold: u32) -> SdustIntervals;
```

### Layer 2 — `DustFilter` iterator adaptor

The thin upstream-wiring wrapper. Holds the current chromosome's
mask and a sweep pointer; delegates all algorithmic work to
`sdust_mask`.

```rust
use crate::per_sample_pileup::pileup::RefSeqFetcher;
use crate::per_sample_pileup::psp::header::ParsedChromosome;
use crate::var_calling::per_position_merger::{
    PerPositionMergerError, PerPositionPileups,
};

pub const DEFAULT_DUST_WINDOW: u32 = 64;
pub const DEFAULT_DUST_THRESHOLD: u32 = 20;

/// Validated DUST filter configuration. Fields are private so that
/// the only way to obtain a `DustFilterConfig` is via [`new`] or
/// [`Default`], both of which guarantee the window/threshold pair
/// is in range. Mirrors the direction the Stage 6 review picked up
/// for `PosteriorEngineConfig` (Mi12 —
/// [PROJECT_STATUS.md:244](../../PROJECT_STATUS.md#L244)).
#[derive(Debug, Clone, Copy)]
pub struct DustFilterConfig {
    window: u32,
    threshold: u32,
}

impl DustFilterConfig {
    /// `window` must be in `3..=4096`. `threshold` is the sdust `T`
    /// parameter — for example `20` for the lh3/sdust default
    /// density-threshold of `> 2.0` triplet-pairs per triplet.
    /// Threshold is unbounded above (a very large `T` is the
    /// effective "do not mask anything" config).
    pub fn new(window: u32, threshold: u32) -> Result<Self, DustFilterError> {
        if !(3..=4096).contains(&window) {
            return Err(DustFilterError::InvalidWindow { window });
        }
        Ok(Self { window, threshold })
    }

    pub fn window(&self) -> u32 { self.window }
    pub fn threshold(&self) -> u32 { self.threshold }
}

impl Default for DustFilterConfig {
    fn default() -> Self {
        Self {
            window: DEFAULT_DUST_WINDOW,
            threshold: DEFAULT_DUST_THRESHOLD,
        }
    }
}

pub struct DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: RefSeqFetcher,
{ /* fields private */ }

impl<I, F> DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: RefSeqFetcher,
{
    /// Wrap `upstream` so that every yielded item whose reference
    /// position is flagged low-complexity by sdust is silently
    /// dropped.
    ///
    /// Infallible — the config has already been validated by
    /// [`DustFilterConfig::new`] and `chromosomes` is taken
    /// verbatim from the upstream merger's
    /// [`PerPositionMerger::chromosomes`].
    pub fn new(
        upstream: I,
        fetcher: F,
        chromosomes: Vec<ParsedChromosome>,
        config: DustFilterConfig,
    ) -> Self;

    pub fn config(&self) -> DustFilterConfig;
}

impl<I, F> Iterator for DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: RefSeqFetcher,
{
    type Item = Result<PerPositionPileups, DustFilterError>;
}

#[non_exhaustive]
#[derive(thiserror::Error, Debug)]
pub enum DustFilterError {
    /// Upstream merger surfaced an error. Boxed for the same
    /// `Result`-size reason as `PerPositionMergerError::PerSampleReader`.
    #[error("upstream merger failed: {0}")]
    Upstream(#[source] Box<PerPositionMergerError>),

    /// Reference fetch failed. Carries the chromosome we were
    /// trying to load.
    #[error("reference fetch failed for chrom {chrom_id}: {source}")]
    RefFetch {
        chrom_id: u32,
        #[source]
        source: std::io::Error,
    },

    /// Upstream yielded a `chrom_id` not present in the
    /// `chromosomes` table the filter was constructed with.
    /// Defensive: in production the table comes from the same
    /// merger that emits the records, so this fires only on
    /// pathological mocks.
    #[error("upstream emitted unknown chrom_id {chrom_id}")]
    UnknownChromId { chrom_id: u32 },

    /// `DustFilterConfig::window` is outside the supported range.
    #[error("invalid DUST window {window}: must be in 3..=4096")]
    InvalidWindow { window: u32 },
}
```

### Wiring

```rust
let merger = PerPositionMerger::new(iters, sample_names, chromosomes)?;
let filter_chromosomes = merger.chromosomes().to_vec();
let fetcher = ChromBoundaryRefFetcher::new(&fasta_path, contigs)?;
let config = DustFilterConfig::new(64, 20)?;
let filter = DustFilter::new(merger, fetcher, filter_chromosomes, config);
for pileups in filter {
    let pileups = pileups?;
    // hand to VariantGrouper / ...
}
```

The filter doesn't own the merger — it consumes it as an iterator,
same shape as `Iterator::filter`.

## Internal state and `Iterator::next`

```rust
struct DustFilter<I, F> {
    upstream: I,
    fetcher: F,
    chromosomes: Vec<ParsedChromosome>,
    config: DustFilterConfig,

    /// `chrom_id` whose mask is currently loaded. `None` before the
    /// first emission and after a reset.
    current_chrom: Option<u32>,
    /// Masked intervals on `current_chrom`. Empty when between
    /// chromosomes.
    current_mask: SdustIntervals,
    /// Sweep pointer into `current_mask`. Indexes the next interval
    /// whose `finish > last_seen_pos0`. Advances monotonically.
    sweep: usize,

    /// Latched on first upstream error or internal error.
    done: bool,
}
```

`Iterator::next`:

1. If `done`, return `None`.
2. Pull `pileups = upstream.next()?`. On `Err`, latch `done` and
   surface as `DustFilterError::Upstream`.
3. **Chromosome switch.** If `Some(pileups.chrom_id) !=
   self.current_chrom`, load the new mask:
   a. Look up `chrom_length` in `self.chromosomes[pileups.chrom_id
      as usize].length`. Out-of-range → latch and surface
      `DustFilterError::UnknownChromId`.
   b. `let seq = self.fetcher.fetch(pileups.chrom_id, 1,
      chrom_length).map_err(|e| DustFilterError::RefFetch {
      chrom_id, source: e })?;`. Latches on error.
   c. `self.current_mask = sdust_mask(&seq, self.config.window,
      self.config.threshold);`
   d. `self.sweep = 0; self.current_chrom = Some(pileups.chrom_id);`
4. **Sweep advance.** Convert `pileups.pos` (1-based per
   [`PileupRecord::pos`](../../src/var_calling/per_position_merger.rs#L40))
   to 0-based: `let pos0 = pileups.pos - 1;`. Advance `self.sweep`
   while `self.sweep < self.current_mask.len() &&
   self.current_mask[self.sweep].1 <= pos0` (interval's `finish` is
   half-open, so `finish <= pos0` means the interval ended before
   `pos0`).
5. **Pass/skip decision.**
   - If `self.sweep < self.current_mask.len()` and
     `self.current_mask[self.sweep].0 <= pos0` (interval's `start`
     ≤ `pos0`), then `pos0` is inside the current interval — drop
     `pileups` and **loop back to step 2** without yielding.
   - Otherwise yield `Ok(pileups)`.

Errors latch — same shape as the merger.

### Coordinate conversion

- `PileupRecord::pos` is **1-based** per the `.psp` spec.
- `sdust_mask` returns **0-based, half-open** intervals.
- The filter converts at the lookup boundary (step 4) and never
  again. Document the convention in the module-level doc comment.

## Threshold semantics

The condition that triggers masking is **strictly greater than**:
`10 · score > T · triplet_count` ([sdust.c:111](../../sdust/sdust.c#L111)).
An interval whose density exactly equals `T / 10` is **not** masked.
A dedicated off-by-one test pins this.

The user-visible threshold is `T = 20` — the lh3/sdust + minimap2
convention — even though the internal math uses `T / 10`. We do not
re-scale to "density `> 2.0`" at the API surface, because everyone
who has ever invoked sdust at the command line knows `T = 20`.
Document the relationship in `DustFilterConfig::threshold()`'s
doc comment.

## Configuration and CLI

Per the spec's parameter table
([calling_pipeline_architecture.md:1017-1023](../specs/calling_pipeline_architecture.md#L1017)):

| Flag | Default | Effect |
|---|---|---|
| *(default)* | on | Apply sdust with `w = 64`, `T = 20`. |
| `--no-complexity-filter` | — | Bypass the filter entirely. |
| `--complexity-window N` | 64 | Override `W` (maximum subinterval length). |
| `--complexity-threshold T` | 20 | Override `T` (score-density threshold; mask iff `10·s > T·L`). |

CLI parser bindings land with the cohort subcommand (same convention
as `posterior_engine` and `contamination_estimation` —
[PROJECT_STATUS.md:282](../../PROJECT_STATUS.md#L282)). This plan
ships only the library API.

The bypass case (`--no-complexity-filter`) is **not** modelled inside
`DustFilter`. The cohort wiring chooses at construction time between
the bare merger and `DustFilter::new(merger, …)`; the filter has no
"identity mode". This keeps the hot path branch-free and the
library API honest: if a `DustFilter` exists, it filters.

## Test strategy

Tests split across two layers, matching the API split.

### `sdust_mask` unit tests

Pin algorithmic behaviour against in-memory ACGT byte slices. No
upstream iterator, no fetcher.

- **Empty input** → empty interval list.
- **Pure-ACGT random-ish high-complexity 256-bp slice, default
  config** → empty interval list. Pins that the filter does not
  over-mask realistic sequence.
- **Pure homopolymer (`"AAAA…A"`, 128 bp), default config** → one
  interval spanning the whole homopolymer. Test on a length
  comfortably larger than the default window to defeat any
  off-by-one at the window edge.
- **Perfect dinucleotide repeat (`"ATAT…"`, 128 bp)** → one
  interval spanning the repeat.
- **Perfect trinucleotide repeat (`"ATGATG…"`, 128 bp)** → one
  interval spanning the repeat.
- **High-complexity flank + homopolymer + high-complexity flank**
  → exactly one interval, covering only the homopolymer (or
  whatever sdust's reference flag for that exact input is — see
  golden-vector below).
- **N-handling.** `"AAAA…N…AAAA…"` should produce two intervals,
  one per ACGT run, because N breaks the triplet stream
  ([sdust.c:151-155](../../sdust/sdust.c#L151)).
- **Lowercase ACGT** behaves like uppercase (per `seq_nt4_table`
  — [sdust.c:22-39](../../sdust/sdust.c#L22)).
- **Off-by-one at the threshold.** Construct (or compute against
  the `lh3/sdust` binary on a small input) a slice whose top
  candidate subinterval has density **exactly** `T / 10`. It
  must **not** be masked (`>`, not `≥`). Pair with a slice that
  scores `T / 10 + ε` and **must** be masked.
- **Threshold sensitivity smoke.** Same homopolymer, run with
  `T = 5000` — now no candidate density exceeds `T / 10 = 500`,
  so the interval list is empty.

### `DustFilter` unit tests

Synthetic upstream (a `std::iter`-based
`Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>`)
plus a stub `RefSeqFetcher` over an in-memory `Vec<u8>` per
chromosome.

- **Empty upstream.** Filter yields nothing immediately; no fetch
  happens.
- **Single chromosome, all positions pass.** All-ACGT-random
  reference; every upstream position is yielded unchanged.
- **Single chromosome, mid-chromosome homopolymer.** Reference has
  a low-complexity stretch in the middle; upstream emits
  consecutive positions across the whole chromosome; filter drops
  exactly the positions inside the masked region. Pin both the
  start and the end of the masked range against the `sdust_mask`
  output (the filter's job is just a sweep over that list).
- **Two chromosomes.** Two reference contigs; upstream alternates
  some positions in chr1, then some in chr2. Confirms the filter
  reloads the mask on chromosome change and drops the previous
  mask. Stub fetcher's call count is `2`, not `(2 * n_positions)`.
- **Cross-chromosome reset.** chr1 ends in a low-complexity tract,
  chr2 starts with a high-complexity region at the same `pos`.
  Position `(chr2, pos)` is **not** masked even though
  `(chr1, pos)` was. Pins that the sweep pointer and mask reset
  on chromosome change.
- **Upstream error latches.** Synthetic upstream yields
  `Err(PerPositionMergerError::OutOfOrder { … })` after a few
  records. Filter surfaces it as `DustFilterError::Upstream`, and
  every subsequent `next()` is `None`.
- **Reference fetch error.** Stub fetcher returns `io::Error` for
  a specific `chrom_id`. Filter surfaces
  `DustFilterError::RefFetch { chrom_id, … }` and latches.
- **Unknown chrom_id.** Upstream yields a `chrom_id` that doesn't
  exist in the filter's `chromosomes` table. Filter surfaces
  `DustFilterError::UnknownChromId` and latches.
- **Invalid config.** `DustFilterConfig::new(2, 20)` →
  `Err(DustFilterError::InvalidWindow { window: 2 })`. Same for
  `window = 4097`. Boundary cases `window = 3` and `window = 4096`
  return `Ok`.
- **Sweep pointer monotonicity smoke.** Stub fetcher records call
  count; with N upstream positions on one chromosome, exactly one
  fetch call happens (not N). Pins that the per-position work
  inside the filter is O(1) amortised.

### Golden-vector test

Goal: lock `sdust_mask`'s per-position decisions to the canonical
`lh3/sdust` implementation. Catches any drift in the port — an
off-by-one in the suffix trim, a wrong threshold direction, an
encoding bug in the triplet word, an N-handling regression.

**Reference tool.** The cloned `lh3/sdust` binary at
[sdust/sdust](../../sdust/sdust) (or wherever the `make` in
[sdust/Makefile](../../sdust/Makefile) drops it). The same C source
we are porting.

**Fixture format.** JSON under `tests/golden/dust_filter/`:

```json
{
  "tool": "lh3/sdust @ <git-sha>",
  "tool_args": ["-w", "64", "-t", "20"],
  "window": 64,
  "threshold": 20,
  "snippets": [
    {
      "name": "homopolymer_with_flanks",
      "bases": "ACGTACGT...AAAAAAAAAAAAA...ACGTACGT",
      "intervals": [[16, 80]]
    },
    ...
  ]
}
```

`intervals` is the list `lh3/sdust` printed for that snippet at the
given `-w` / `-t`. Each entry is `[start, finish]`, BED-style 0-based
half-open — matches sdust's stdout exactly.

**Snippet selection.** Six snippets, each 200-500 bp, chosen to
cover:
1. A pure homopolymer flanked by high-complexity sequence.
2. A perfect dinucleotide repeat (e.g. `"ATAT…"`) flanked by
   high-complexity sequence.
3. A perfect trinucleotide repeat.
4. Mixed-complexity tract with low-complexity stretches abutting
   high-complexity sequence — catches off-by-one at the window
   edge.
5. Snippet containing N inside an otherwise low-complexity stretch
   — catches N-handling drift.
6. Snippet engineered to score close to threshold throughout —
   catches `>` vs `≥` drift in the threshold comparison.

**Test mechanics.** For each snippet:
1. Call `sdust_mask(snippet.bases.as_bytes(), fixture.window,
   fixture.threshold)`.
2. Assert the returned `Vec<(u32, u32)>` equals
   `snippet.intervals` exactly.

**Regeneration.** A short shell script at
`tests/golden/dust_filter/REGENERATE.sh` records the exact
`sdust -w <W> -t <T>` invocations that produced each snippet's
interval list. Re-runnable if the cloned `lh3/sdust` HEAD changes
or if the fixture needs extending. Document the `sdust` build
step (`cd sdust && make`) in `REGENERATE.sh` so anyone can run it
from scratch.

## Validation

Inside the dev container (`./scripts/dev.sh`):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo build --examples`
- `cargo build --benches`

No new benchmark in this plan. The `sdust_mask` call is `O(chrom_len)`
per chromosome (amortised — bounded by the suffix-trim invariant);
the per-position sweep is `O(1)` amortised. Benching against a
synthetic input would over-fit. A benchmark is justified once the
cohort CLI runs end-to-end and we can measure DUST's share of real
cohort wall time.

## Assumptions / silent choices

- **Per-chromosome batch, not streaming.** sdust runs once per
  chromosome on the full reference; the iterator adaptor just
  sweeps the resulting interval list. Simpler than a streaming
  port, matches the canonical reference exactly, and the memory
  cost (tens of MB for the interval list + ~chrom_len for the
  reference buffer) is bounded and small relative to the rest of
  the pipeline.
- **Config is validated at construction.** `DustFilterConfig::new`
  returns `Result`; `DustFilter::new` is infallible. Mirrors the
  posterior-engine review's Mi12 direction.
- **`>` not `≥` at the threshold.** `10 · score > T · L`, matching
  [sdust.c:111](../../sdust/sdust.c#L111). Pinned by a dedicated
  test.
- **`T = 20` at the API surface, not the equivalent density `2.0`.**
  Matches the universal sdust convention.
- **N (and any other non-ACGT byte) breaks the triplet stream.**
  Direct port of [sdust.c:22-39 + sdust.c:151-155](../../sdust/sdust.c#L22).
  Lowercase ACGT is treated as uppercase (same table entry).
- **No bypass mode inside `DustFilter`.** `--no-complexity-filter`
  is handled at the call site (skip filter construction entirely).
- **Filter generic on `Iterator<…>` + `RefSeqFetcher`.** Tests use
  `std::iter`-based upstreams and `Vec<u8>`-backed stub fetchers.
- **Errors latch.** Once any upstream error, fetch error, or
  unknown-chrom_id error surfaces, every subsequent `next()`
  returns `None`. Same shape as the walker and merger.
- **Coordinate conversion happens at the sweep boundary only.**
  1-based `PileupRecord::pos` → 0-based for sdust lookup; the
  filter never operates on mixed coordinate systems internally.

## Risks

- **Port fidelity to `lh3/sdust`.** The C source uses a small
  number of subtle invariants (the suffix-trim bound, the
  "dominated by existing perfect interval" cross-multiplied
  comparison at [sdust.c:117](../../sdust/sdust.c#L117), the
  flush-on-N convention). The required golden-vector test against
  the cloned binary is the load-bearing mitigation.
- **Reference availability for fetch.** A chromosome whose recorded
  `length` exceeds what the FASTA actually contains surfaces as
  `DustFilterError::RefFetch`. Correct shape, but the user-facing
  error is on a chrom_id, not a file path — the fetcher's own
  error message carries the path detail via the boxed
  `io::Error`.
- **Memory transient during sdust execution.** ~250 MB for human
  chr1 reference bytes + masked-interval list during the
  `sdust_mask` call. Dropped on chromosome change. The walker
  already holds one chrom's worth of reference at a time, so this
  is the established footprint; no surprise.
- **`ChromBoundaryRefFetcher` is non-`Sync`.** Stage 3 is
  intentionally sequential (per the spec's "streaming filter"
  wording and per the merger's single-threaded consumer pattern),
  so this is fine for v1. If a future per-chromosome
  parallelisation of Stages 3-5 emerges, switch to
  `SyncRefFetcher`
  ([ref_fetcher.rs:109-138](../../src/per_sample_pileup/ref_fetcher.rs#L109)).

## Out-of-scope follow-ups

- **Streaming sdust integrated with the upstream iterator.** Only
  worth pursuing if the per-chromosome-batch's transient memory
  becomes a measured problem on a real run.
- **Per-region threshold overrides.** A BED-driven map of
  `(region → threshold)` for selectively tightening or relaxing
  the filter in problem regions.
- **Cross-stage fusion** (e.g. fuse DUST inside the merger to save
  one iterator hop). Don't pre-optimise; profile first.
- **Per-chromosome parallelism.** Stage 3 can in principle run one
  worker per chromosome, paired with `SyncRefFetcher`. Only
  worth doing if cohort wall time becomes Stage-3-bound, which is
  unlikely given the per-chromosome sdust pass is fast.

## File touch list

- `src/var_calling/dust_filter.rs` — new file: `sdust_mask`,
  `DustFilter`, `DustFilterConfig`, `DustFilterError`,
  `SdustIntervals`, `DEFAULT_DUST_WINDOW`, `DEFAULT_DUST_THRESHOLD`,
  full `#[cfg(test)]` module.
- `src/var_calling/mod.rs` — `pub mod dust_filter;` added alongside
  the existing stage modules.
- `tests/golden/dust_filter/` — JSON fixture + `REGENERATE.sh`.
- `PROJECT_STATUS.md` — flip Stage 3's status block from "not yet
  planned" to "planned" with a link to this plan; update the
  "Next task" standing-candidates list.
