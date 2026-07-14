# A tandem-repeat scanner (a reusable sequence primitive)

*Status: design spec (2026-07-14; decoupled from the catalog, then given a streaming region
seam + satellite cap, 2026-07-14 — all open questions resolved, ready to implement). **No code
yet — this settles the design.** The scanner is a general tandem-repeat detector with two
interfaces: a low-level interval finder (`sequence + period range → repeat intervals`) and a
memory-bounded, windowed **region seam** that streams a repeat/satellite/unique tiling of the
reference (`RegionScanner`, §3.6). It lives at `src/ng/tandem_repeat.rs`. Its first consumer
replaces the external `trf-mod` shell-out in Stage 0
(`ssr-catalog`, [`../architecture/ssr_catalog.md`](../../architecture/ssr_catalog.md) §2–§3); the
ng snp/str caller is the second, consuming the region seam. The algorithm knows nothing about
either. Its code-facing companion is [`../arch/ssr_repeat_scanner.md`](../arch/ssr_repeat_scanner.md)
(types & interfaces). Naming follows the project convention: **STR** in prose, `ssr`/repeat terms
in code.*

*Clean-room note (legal, not stylistic). This design is derived **only** from (a) the
published Benson (1999) tandem-repeats-finder paper — the idea of comparing a window to a
period-shifted copy and growing runs of matches — and (b) our own code. It is **not** derived
from the AGPL-v3 `TRF-mod` source, which was not read. The scanner is a fresh, much smaller
construction (short periods only, no statistical-significance model, no wraparound alignment),
built on the Ruzzo–Tompa maximal-scoring-subsequence algorithm (Ruzzo & Tompa 1999 —
independent of TRF).*

---

## 1. What the scanner *is* — a decoupled algorithm, not a catalog part

The scanner is a **pure function over a sequence**: give it a byte slice and a period range,
and it returns the tandem-repeat intervals it finds, each as `(start, end, period, score)`. It
holds **no policy** about *why* you want repeats — no catalog, no purity floor, no
homopolymer rule, no motif convention. Those belong to whoever calls it. This decoupling is the
central design decision (it was previously entangled with the catalog, and is now deliberately
separated): **the algorithm accepts the range of repeat periods to explore, and that is the
whole of its configuration surface beyond the two scoring weights.**

It exposes **two interfaces** over one core (§3): a low-level **interval finder**
(`find_tandem_repeats`) that returns raw, possibly-overlapping candidates, and a high-level
**region seam** (`RegionScanner`, §3.6) that streams a resolved repeat-vs-unique tiling of the
reference. Two consumers are anticipated, each suited to one interface, and the design serves
both without favouring either:

- **The STR catalog (Stage 0)** — its first user, replacing `trf-mod`. It uses the low-level
  **interval finder** per contig (its post-filter needs per-interval periods for its
  compound/bundle logic) and feeds the intervals to the existing post-filter (§6). All
  catalog-specific choices live there, not in the scanner.
- **The SNP/SSR caller** — expected to use tandem-repeat context "as we see fit" (routing loci,
  annotating or down-weighting indels in repetitive regions). It consumes the **region seam**:
  as it walks the genome, `RegionScanner` tells it which spans are repetitive and which are
  unique. The scanner must not bake in any assumption that would make this use awkward — which
  is exactly why period range is an argument, not a constant, and why the seam is a streaming
  iterator that drops into the caller's existing walk.

**Goals of the algorithm:**

- **Use-agnostic.** The only inputs are the sequence, the period range, and the scoring
  weights. No consumer's policy leaks in; the output is raw intervals a caller can filter,
  trim, or annotate however it likes.
- **Impure-repeat tolerant** for the cases that occur in practice — substitutions *and* single
  indels stay one interval (§3.4), an intrinsic property of the lag-*p* method, not a bolt-on.
- **Simple and fast.** A linear scan per period, O(n · P) for P periods over a sequence of
  length n, with an O(n) Ruzzo–Tompa pass. Deterministic (a pure function of its inputs).
- **Self-contained, MIT-clean.** No external process, no AGPL surface, one small module.

**Non-goals — of the *algorithm* (deliberately excluded, a fuller detector could do them):**

- **No full wraparound / indel-aware DP alignment.** The scanner does not run TRF's wraparound
  dynamic program. It is still impure-tolerant for substitutions and single indels via the
  local lag-*p* comparison (§3.4); what it forgoes is only the pathological tail — a large
  indel, or several clustered indels, whose combined mismatch burst outruns the flanking copies
  can still split a tract into two intervals (§3.4). Rare in practice.
- **No statistical-significance model.** TRF scores a repeat against a random-DNA null and
  reports a significance. The scanner reports a raw alignment-style score only; deciding what
  score/length/purity is "significant" is the caller's, not the algorithm's.
- **No motif, copy-number, or purity output.** The scanner reports coordinates, period, and a
  score. Anything richer (the motif bytes, copy number, purity, canonical class) a caller
  derives from the sequence itself.

**The scanner does *not* know about, and never references:** the catalog, `postprocess.rs`,
period-1/homopolymer policy, a period ≤ 6 ceiling, `ref_seq`, the CLI, or any file format.
Those appear only in the **consumer** sections (§6, §7), never in the algorithm (§3–§5).

---

## 2. Algorithmic grounding (from the paper, not the source)

Benson's insight, reduced to what a short-period finder needs: **a tandem repeat of period
*p* is a run of positions where base *i* matches base *i − p*.** Comparing the sequence to a
copy of itself shifted by *p* turns "is this region periodic with period *p*?" into "is there
a long run of matches in the shifted comparison?". Substitutions inside the tract are isolated
mismatches in that run; the run survives them when the surrounding matches outweigh them. That
is the whole idea we take: **compare to a *p*-shifted copy; grow maximal runs of matches with
a mismatch-tolerant cutoff; a run *is* a repeat interval of period *p*.**

Everything statistical in the paper (the probabilistic match model, the significance criteria,
the *k*-tuple detection heuristic that makes the *general*-period case tractable) we **do not
need**, because a caller asks for a bounded set of periods and we brute-force each one directly
— one linear pass per period is cheap and removes the heuristic entirely. What remains is a
scoring-and-segmentation problem, which the Ruzzo–Tompa maximal-scoring-subsequence algorithm
solves exactly in one linear pass (§3.2).

---

## 3. The algorithm

The scanner exposes **two** interfaces over the same core (types in §5), for the two shapes of
consumer:

```rust
/// Low-level: find every tandem-repeat interval in `seq` whose period lies in
/// `periods`. Returns raw, possibly-overlapping intervals (one region can match
/// at several periods). Pure and total — case-insensitive, non-ACGT never matches.
/// Consumers that do their own overlap resolution use this (e.g. the STR catalog,
/// whose post-filter needs per-interval periods for its compound/bundle logic, §6).
pub(crate) fn find_tandem_repeats(
    seq: &[u8],
    periods: PeriodRange,     // the caller's requested period range — the only scope knob
    params: &ScanParams,      // match/mismatch weights + the minimum-copies emission floor
) -> Vec<RepeatInterval>;

/// High-level: walk `seq` and yield an ordered, gap-free tiling of it into
/// tandem-repeat regions and unique (non-repetitive) regions — the **routing seam**
/// the snp/ssr caller consumes to decide, position by position, which path a locus
/// takes (§3.6). Built on `find_tandem_repeats` + a coverage merge.
pub(crate) struct RegionScanner<'a> { /* … */ }   // impl Iterator<Item = Region>
```

The low-level finder runs the following for each period `p` in `periods.min..=periods.max`,
independently, and concatenates the results (§3.1–§3.5); the region seam is §3.6.

### 3.1 One score signal per period

Plain-English: for a candidate period *p*, walk the sequence and, at each position *j ≥ p*,
ask "does base *j* equal base *j − p*?". A yes earns a **match reward**; a no costs a
**mismatch penalty**. A tandem repeat of period *p* shows up as a stretch where the running
total climbs — a high-scoring segment of this signal.

Precisely, for `seq[0..n]` and period `p`, the integer score at each position `j ∈ [p, n)`:

```
score[j] =  +params.match_reward       if is_base(seq[j]) && upper(seq[j]) == upper(seq[j-p])
            -params.mismatch_penalty    otherwise           // includes any non-ACGT base (§3.5)
```

`is_base(b)` is true only for `A/C/G/T` (case-insensitive). The two weights are positive
integers; their ratio is the mismatch tolerance (§3.3).

### 3.2 Maximal scoring segments = repeat intervals (Ruzzo–Tompa)

A repeat interval is a **maximal scoring subsequence** of `score[]`: a contiguous run whose
total is positive and that cannot be extended, nor trimmed at either end, without lowering its
score. The Ruzzo–Tompa algorithm finds **all** of them in a single O(n) left-to-right pass with
a small stack — no seeding heuristic, no arbitrary window, no per-period parameter beyond the
two weights. It is a standard published algorithm (Ruzzo & Tompa 1999), unrelated to TRF, and
it is exactly the "grow maximal runs with a drop cutoff" behaviour, made deterministic and
parameter-light.

A maximal segment spanning score-positions `[j0, j1]` (inclusive) maps back to the tract it
certifies. Because `score[j]` compares `seq[j]` to `seq[j-p]`, the earliest base involved is
`j0 - p` and the latest is `j1`, so the interval as **0-based half-open** is:

```
start = j0 - p
end   = j1 + 1
```

Worked example — a perfect `(CAG)k` tract of length `L` beginning at position `a`: matches at
period 3 occur at `j = a+3 … a+L-1`, so `j0 = a+3`, `j1 = a+L-1`, giving `start = a`,
`end = a+L` — the exact tract. The `score` is the segment total (`≈ match_reward · (L − 3)`
when perfect), an alignment-style length-and-purity proxy the caller may rank or gate on.

**Emission floor.** A segment is emitted only if its implied copy count `(end − start) / p ≥
params.min_copies` (default 2 — a repeat needs at least two copies to exist at all). This is a
general "what counts as a repeat" knob, defaulted sanely and caller-overridable; it is *not* a
consumer's copy-number policy (the catalog applies its own, stricter, per-period floors *after*
the scanner, §6).

### 3.3 Mismatch tolerance — a scoring ratio

The scoring ratio governs sensitivity. A region of purity ρ (fraction of matching positions)
has expected per-position score `ρ · match_reward − (1 − ρ) · mismatch_penalty`; for that to
stay positive — so a tract of purity ρ stays one merged segment rather than fragmenting — we
need

```
mismatch_penalty / match_reward  <  ρ / (1 − ρ).
```

So the ratio *is* the "lowest purity that survives as one interval" knob: a ratio of `r` holds
tracts down to purity `r / (1 + r)`. **Defaults: `match_reward = 2`, `mismatch_penalty = 7`**
(ratio 3.5 → holds tracts to purity ≈ 0.78). The defaults are the scanner's, chosen to be a
reasonable general tolerance; a consumer with a specific purity target (the catalog's floor is
0.8) can pass weights that match it. *(That the default ratio lands near TRF's long-documented
+2 / −7 alignment weights is corroboration from published parameters, not a port — those
weights are numeric constants in the paper and manual, not source we read.)*

A single substitution flips **two** comparisons to mismatches (position `k` vs `k−p`, and `k+p`
vs `k`), a local cost of `2 · mismatch_penalty`, recovered by the surrounding matches — so
isolated substitutions never split a segment. A single indel is a *bounded* burst, not a
permanent break (§3.4).

### 3.4 Impure repeats — substitutions and single indels stay one interval

This is the property that makes the scanner useful on real, imperfect sequence — **and it is
intrinsic to the method, not an added pass** — because the lag-*p* comparison `seq[j]` vs
`seq[j−p]` is *local*. It asks "is there a copy of me *p* positions back?", which for a
period-*p* tract stays true **regardless of an upstream indel**. An indel therefore does *not*
desync the rest of the tract; it only disturbs the ~*p* comparisons whose window straddles it.

- **Substitution** — a 2-position mismatch blip (§3.3); the segment stays merged.
- **Single indel** (insertion or deletion, any length) — a *bounded* burst of ~*p*+1 mismatches
  where the comparison window straddles the indel, after which the downstream copies' own
  period-*p* structure is picked straight back up and matches resume. The burst costs
  ~`(p+1) · mismatch_penalty`; any tract of a few copies each side has far more flanking match
  score than that, so the Ruzzo–Tompa segment does **not** split — the whole interrupted tract
  is emitted as **one** interval `[start, end)` spanning the indel. (Contrast a *fixed-anchor*
  alignment, which an indel desyncs permanently; the lag-*p* comparison is immune to that.)
- **Whole-motif-multiple indel** — not even a burst: a pure copy-number change that preserves
  phase, so every comparison still matches. One clean segment.

So the scanner emits the **full** interrupted tract as a single interval — it never throws the
downstream half away. A consumer that wants only the clean in-phase part (as today's catalog
post-filter does, by trimming) can do that itself; a consumer that later wants to keep the
whole interrupted tract (phase-aware) already has the full span to work from. Keeping the full
segment rather than pre-trimming is deliberate: **the algorithm should not pre-decide how much
of an interrupted repeat its caller wants** (§7 records the catalog and future homes).

Only the pathological tail escapes: a *large* indel, or several indels clustered within one
tract, can produce a burst that outruns the flanking score and splits the segment. Rare, and a
caller that cares can stitch same-period neighbours (deferred, §8). It is not built in: the
local-comparison property already covers the cases that occur in practice.

### 3.5 Non-ACGT and case

The scanner **upper-cases** each base before comparing (so soft-masked, i.e. lowercase,
repeats are still found) and treats any non-`ACGT` base as a **non-match** (so an `NNNN…` run,
where `N == N` would otherwise read as a perfect repeat, scores negatively and yields nothing).
Both are general, use-agnostic rules — a caller never has to pre-clean its sequence.

### 3.6 The region-segmentation seam — repeat vs unique, streamed

The raw intervals are the right input for a consumer that resolves overlaps itself (the
catalog). But the **caller's router** wants the opposite: a single, resolved answer to "walking
the reference, where am I inside a tandem repeat and where am I in unique sequence?" — so it can
send repeat regions down the STR/indel-aware path and unique regions down the plain SNP path.
That is a *partition* of the sequence, not a bag of overlapping intervals. `RegionScanner`
provides it as an **iterator**, so it composes with the caller's existing genome walk and
streams rather than materialising a whole-contig `Vec`.

**What it yields.** An ordered, non-overlapping, **gap-free** tiling of the scanned span into
three region kinds, each maximal (two of the same kind never abut):

```rust
pub(crate) enum Region {
    /// A genotypeable tandem repeat: repeat coverage no longer than `max_repeat_len`
    /// (§ satellite cap). The union of one or more overlapping `find_tandem_repeats`
    /// intervals (periods may differ). This is the STR the caller routes to the STR path.
    Repeat(RepeatRegion),
    /// Repeat coverage *longer* than `max_repeat_len` — satellite DNA. Repetitive, so
    /// NOT unique sequence, but too long to be a genotypeable STR locus, so excluded
    /// from STR analysis. The caller neither SNP-calls nor STR-genotypes it (mask/skip).
    Satellite(RegionSpan),
    /// No tandem-repeat coverage — unique sequence, the SNP path's domain.
    Unique(RegionSpan),
}
```

A `Repeat` region carries its merged span **and** the constituent intervals, so the caller can
act on the repeat structure (dominant period/motif for indel context, or just "this is
repetitive") without re-scanning; `Satellite`/`Unique` carry only a `RegionSpan`:

```rust
pub(crate) struct RepeatRegion {
    pub span: RegionSpan,                       // union of the intervals' spans
    pub intervals: Box<[RepeatInterval]>, // ≥1, coordinate-ordered; the overlaps merged here
}
```

**The satellite cap — a max STR length (decided: default 1 kb).** We do not want to analyse
satellite DNA as if it were a microsatellite. So a merged repeat region longer than
`max_repeat_len` (default **1000 bp** — comfortably above any real STR locus) is emitted as
`Satellite`, not `Repeat`, and never enters STR analysis. This cap does double duty: it also
**bounds the window halo** that makes streaming correct (below).

**Streaming over a whole contig, memory-bounded (in scope, not deferred).** The seam does *not*
require the contig in memory. It walks the reference in windows using the existing
[`ChromRefFetcher`](../../../../src/fasta/fetcher.rs) sliding-buffer reader
(`StreamingChromRefFetcher`, ~1 MB resident regardless of contig size — the same primitive the
cohort per-chromosome workers and DUST already use for monotonic forward scans). For each
window it fetches a core of `window_bp` (default ~100 kb) plus a right **halo** of
`max_repeat_len` and a small left margin of `periods.max`, runs `find_tandem_repeats` on that
slice, merges coverage, and **yields the regions whose start falls in the core** — so a
repeat spanning a window cut is emitted exactly once, by the window whose core contains its
start (the per-chunk-DUST halo/attribution pattern). Because an STR is ≤ `max_repeat_len`, the
halo guarantees any core-starting STR is fully contained in the fetched slice and never split;
a run that reaches `max_repeat_len` is a `Satellite`, coalesced across windows into one span as
the scan streams through it. Peak memory is therefore ~`window_bp + max_repeat_len`, not the
contig length — the memory-efficiency thesis, applied to detection.

**How one window is resolved (coverage merge).** Run `find_tandem_repeats`, sort the intervals
by start, sweep once to union overlapping/abutting intervals into merged spans (grouping the
intervals in each), classify each merged span `Repeat` vs `Satellite` by the cap, and emit
`Unique` for every gap — so the output tiles the core exactly. The sweep is O(m) after an
O(m log m) sort (m = interval count in the window).

**Contract.** Segments are coordinate-ordered, pairwise non-overlapping, and their union is
exactly the scanned span; every `Repeat.intervals` is non-empty and its `span` is their union;
an empty input yields nothing; the output is **window-count-invariant** (a pure function of the
contig, `periods`, `params`, and options — `window_bp` changes memory, not the regions).
`RegionScanner` follows the repo's iterator-as-seam convention (the ng `ReadFilter` in
[`read_filtering.md`](read_filtering.md) §5 is the shape to match — a struct
implementing `Iterator`).

**Optional smoothing (a knob, defaulted off).** A raw coverage partition can fragment — a 1 bp
unique gap between two repeats, or a 4 bp lone repeat blip. Two optional parameters tame that
for the router: `merge_gap` (bridge unique gaps shorter than this into the flanking repeat) and
`min_repeat_len` (reclassify sub-threshold repeat regions as unique). Both default to 0 (pure
coverage partition); a router sets them to taste. These plus `window_bp` and `max_repeat_len`
live in a small `SegmentOptions` passed to `RegionScanner::new`, kept separate from `ScanParams`
because they shape the *partition and the walk*, not the *detection*.

---

## 4. Where the module lives

Because the scanner is shared (STR catalog now, ng snp/str caller later) and holds no
consumer-specific policy, it does **not** live under `src/ssr/` — that tree is the STR *caller*,
and burying a general primitive there would signal ownership it does not have and complicate
reuse. It lives in the **ng module tree**, the home of the new caller generation that primarily
drives it (it is a shared primitive *in* that tree, not one of the pipeline steps):

```
src/ng/tandem_repeat.rs     // the scanner: find_tandem_repeats + its types + #[cfg(test)] tests
```

**Decided (2026-07-14):** `src/ng/tandem_repeat.rs`, promoted to a `src/ng/tandem_repeat/`
directory only if it outgrows one file. It stays **policy-free**; the catalog consumes it
cross-module from `src/ssr/catalog/`, a normal dependency (§6). Rejected `src/ssr/tandem_repeat.rs`
(misrepresents scope — under the STR caller — and complicates reuse). *Supersedes an earlier draft
that leaned a top-level `src/tandem_repeat.rs` sibling to `src/baq/`; the ng tree was chosen as the
home of the caller generation the scanner primarily serves, and it remains use-agnostic there.*

---

## 5. The types the algorithm introduces

All small, use-agnostic, and owned by `tandem_repeat` (a consumer imports them, never
redefines them):

```rust
/// The inclusive range of periods (motif lengths, in bp) to scan for. The caller
/// chooses it; the algorithm imposes no ceiling of its own. `min >= 1` and
/// `min <= max` (validated on construction — period 0 is meaningless).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PeriodRange { min: u8, max: u8 }   // checked constructor; e.g. new(1, 6)

/// Scoring for the lag-p self-comparison plus the minimum-copies emission floor.
/// The match/mismatch ratio is the mismatch tolerance (§3.3). `Default` is the
/// scanner's general default (2 / 7 / 2), which a consumer may override.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ScanParams {
    pub match_reward: i32,       // > 0; default 2
    pub mismatch_penalty: i32,   // > 0; default 7
    pub min_copies: u32,         // emission floor, default 2
}

/// One detected tandem-repeat interval. Coordinates are 0-based half-open
/// (`[start, end)`); `period` is the lag it was found at; `score` is the
/// Ruzzo–Tompa segment total. Nothing catalog- or caller-specific here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct RepeatInterval {
    pub start: u32,
    pub end: u32,
    pub period: u8,
    pub score: i32,
}
```

`RepeatInterval` is deliberately the algorithm's *natural* output shape, not a shape borrowed
from any consumer. It happens to carry exactly the four fields the catalog post-filter needs
(§6), which is why the catalog can consume it directly — but that is a convenience, not a
coupling: the SNP caller reads the same fields for its own purposes.

The **region-seam** types — the three-kind `Region`, `RepeatRegion`, and the half-open
`RegionSpan { start, end }` — are defined with the seam in §3.6. `RegionScanner` streams over a
[`ChromRefFetcher`](../../../../src/fasta/fetcher.rs) (windowed, memory-bounded), with a `&[u8]`
convenience constructor for small/test inputs. Its `SegmentOptions` — separate from `ScanParams`
because they shape the partition and the walk, not the detection — carry the knobs:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SegmentOptions {
    pub max_repeat_len: u32,   // above this a repeat is Satellite, not a genotypeable STR; default 1000
    pub window_bp: u32,        // streaming window core; memory knob, region-invariant; default ~100_000
    pub merge_gap: u32,        // bridge unique gaps shorter than this; smoothing, default 0 (off)
    pub min_repeat_len: u32,   // reclassify sub-threshold repeat blips as unique; smoothing, default 0 (off)
}
```

---

## 6. Consumer #1 — the STR catalog (replacing `trf-mod`)

This section is **catalog policy**, kept out of §3–§5 on purpose. The catalog's Stage-0
orchestrator (`catalog::run`) stops spawning `trf-mod` per contig and instead calls
`tandem_repeat::find_tandem_repeats` on each contig, then feeds the intervals to the existing,
**unchanged** post-filter.

- **Period range — the consumer's choice, passed in (§9.4).** The module default is
  `PeriodRange::new(1, 6)`; the catalog / ng snp-str caller passes `2, 6` (period-1 homopolymers
  are the post-filter's to drop anyway). The *scanner* has no opinion — the ceiling of 6 and the
  homopolymer question are consumer concerns, and may change after the experiments.
- **The catalog uses the raw interval finder, no satellite cap.** It consumes
  `find_tandem_repeats` directly (flat `RepeatInterval`s), *not* the `RegionScanner` seam — its
  post-filter needs per-interval periods for compound/bundle logic, and it must reproduce the
  current `trf-mod`→postprocess loci for parity, so it applies **no** `max_repeat_len` cap (that
  cap is a live-caller routing feature, not catalog policy; adding it would diverge from the
  golden). The satellite cap, windowing, and region tiling are the *other* consumer's seam.
- **The catalog must pre-filter before `build_loci` (validated finding, Milestone D).** `trf-mod`
  hands the post-filter a *clean* candidate set (significant repeats, redundancy eliminated). The
  raw scanner is deliberately permissive (`min_copies = 2`), so it also emits low-copy noise (in
  aperiodic sequence) and every period-multiple of a real tract; fed straight to `build_loci` that
  noise trips `drop_bundles` — which runs *before* the copy-number floor — and cascades the real
  loci away. So `catalog::run` must apply, before `build_loci`, the two cleanups `trf-mod` bakes
  in: the **per-period copy floor** and **period-multiple redundancy elimination** (`IsRedundant`).
  This is *catalog policy* — it stays out of the use-agnostic scanner and the unchanged
  post-filter. With it, the scanner reproduces the golden catalog at **16/16 recall** (the parity
  test `src/ssr/catalog/scanner_parity.rs`); it is required whenever the production swap lands.
- **The post-filter is unchanged, consuming `RepeatInterval` directly.** `postprocess::build_loci`
  currently takes `Vec<TrfRecord>` and reads only `start/end/period/score` — exactly
  `RepeatInterval`'s fields. So `TrfRecord` is **deleted** and `build_loci` (and its internal
  `drop_bundles`/`is_close`) take `Vec<RepeatInterval>` / `&RepeatInterval` instead. This is a
  mechanical type substitution — the post-filter *logic* (period-2–6 scope, compound-drop,
  bundle-drop, minimal-trim, purity recompute, `ref_seq` embed, homopolymer drop) is untouched,
  so its behaviour is preserved by construction (which is what makes the parity oracle valid,
  below). The `frac_match`/`pattern` fields `TrfRecord` carried were never read; they vanish.
- **All catalog policy stays in the catalog.** Period-1 homopolymers are dropped by
  `postprocess::MIN_PERIOD` as today; the per-period copy-number floors, purity floor, compound
  and bundle drops, and `ref_seq` embedding are all unchanged post-filter steps. The scanner
  contributes none of this.
- **Catalog header provenance.** `CatalogHeader::trf_mod_version` (an external binary's version)
  is replaced by a **built-in-detector** field — `detector` (e.g. `"tandem-repeat-scanner
  <tool_version>"`) plus the scan weights — so the header stays self-describing. A small,
  localised change to [`io.rs`](../../../../src/ssr/catalog/io.rs) and the few test constructors
  that build a header; **outside** `postprocess.rs`, and the only format-visible change.

**Validation of this consumer — golden-fixture *catalog* parity** (the real risk; mirrors the
BAQ port precedent, `src/baq/` + `baq_tests.rs`). Before the swap lands, run the **current**
`trf-mod` → `postprocess` path in the dev container on a small real multi-contig reference (a
few hundred kb of genuine STR diversity) and snapshot the resulting **`Locus` set** as a
committed fixture. The parity test then runs `find_tandem_repeats` → `build_loci` on the same
reference and compares its `Locus` set to the golden one. Comparing the *final* loci (not raw
intervals) is deliberate: the post-filter normalises away most detector-boundary differences,
so this measures what we actually ship. Because single indels stay one interval (§3.4), today's
single-phase post-filter trims a phase-shifted tract to its clean in-phase portion exactly as
it does for `trf-mod`'s own spanning call — so no parity gap arises there. The parity **bar**
is a decision for review (§9): *leaning ≥ 99 % recall of golden loci with a reviewed,
understood diff*, not byte-exact set-equality (different detectors; post-filtering absorbs
boundary noise). Extra scanner-only loci are inspected too (should be genuine STRs TRF's
significance model rejected, not post-filter escapes).

**End-to-end:** drive the scanner-built catalog through the STR benchmark harness
([`benchmarks/ssr_tomato1/`](../../../../benchmarks/ssr_tomato1/)) and confirm genotyping accuracy
and HipSTR-concordance are unchanged from the `trf-mod`-built catalog.

---

## 7. What becomes removable (after catalog parity is banked)

Once `catalog::run` calls `find_tandem_repeats` and parity + e2e pass, the `trf-mod` plumbing
is dead and deleted:

| removed | file |
|---|---|
| `run_on_contig`, `parse_bed_line`, `locate_trf_mod`, `version`, `TrfRecord`, `TRF_MOD_BIN`, `BED_COLUMNS` | `src/ssr/catalog/trf.rs` (whole file) |
| `--trf-mod-path`, `--temp-dir` args + their `CatalogConfig` fields | `ssr_catalog.rs`, `catalog/mod.rs` |
| `CatalogError::{TrfModNotFound, TrfSpawn, TrfRun, TrfVersion, TrfParse}` | `catalog/mod.rs` |
| `trf_mod_version` header field (→ `detector`, §6) | `catalog/io.rs` + test constructors |
| the `trf-mod` clone/build/install step | `Containerfile` (lines ~101–114) |

Removing `--temp-dir` also drops Stage 0's disk-scratch requirement entirely — the scanner
needs no temp files.

---

## 8. Cross-cutting concerns

- **Performance & memory.** O(n · P) with an O(1) per-position step and an O(n) Ruzzo–Tompa
  stack; the caller already holds the sequence resident. For the catalog this is strictly
  *less* work than the subprocess + FASTA-write + BED-read it replaces, and the existing
  per-contig rayon fan-out parallelises it unchanged; determinism (a pure function of inputs)
  is preserved across thread counts.
- **Error model.** The scanner is total over any byte input (non-`ACGT` → non-match), so
  `find_tandem_repeats` returns `Vec<RepeatInterval>` infallibly — no error type of its own.
  `PeriodRange::new` is the one fallible surface (rejects `min < 1` or `min > max`). The
  catalog's `Trf*` error variants are deleted (§7).
- **Concurrency.** No shared state, no subprocess, no temp files — a pure function, safe to call
  from any thread; a *simpler* concurrency story than the temp-dir-per-task spawn it replaces.

---

## 9. Reuse map & open questions

**Reuse map**

| what | existing code | how reused |
|---|---|---|
| the STR post-filter (period scope, compound/bundle drop, trim, purity, `ref_seq`, homopolymer drop) | `postprocess::build_loci` | **unchanged** logic; input type `TrfRecord` → `tandem_repeat::RepeatInterval` (§6) |
| catalog header / writer / reader | `catalog/io.rs` | reused; one field swap `trf_mod_version` → `detector` (§6) |
| per-contig fan-out + ordered collector | `catalog::run`, `ssr_catalog.md` §8 | unchanged; `find_tandem_repeats` slots in where `run_on_contig` was |
| sliding-buffer reference reader | `ChromRefFetcher` / `StreamingChromRefFetcher` (`src/fasta/fetcher.rs`) | `RegionScanner` streams a contig in windows through it (§3.6) — ~1 MB resident, monotonic forward, as the cohort/DUST path uses it |
| per-chunk halo / split-invariance | cohort per-chunk DUST (`project_cohort_perf_dust_per_chunk`) | the window halo + attribution + window-count-invariance gate (§3.6, §9 tests) |
| parity-fixture pattern | `src/baq/`, `baq_tests.rs` | template for the §6 golden-catalog parity test |
| **parity oracle** | committed golden catalog from the current `trf-mod` path (§6) | the `Locus` set the scanner+post-filter must reproduce |

**Tests.** The *algorithm* gets its own use-agnostic unit tests beside `tandem_repeat.rs`,
independent of any consumer. For the **interval finder**: a clean tract (exact
`start/end/period/score`); a substitution-interrupted tract (one segment); a
**single-indel**-interrupted tract (one span — the impure-tolerance property); a
large/multi-indel tract (the pathological split); an `N`-run (nothing); a soft-masked tract
(found); a period-2 tract that also matches at p = 4/6 (all emitted — the finder does *not*
de-duplicate periods; that is a consumer's job); and `PeriodRange` validation. For the
**region seam** (`RegionScanner`): the tiling contract (segments ordered, non-overlapping, union
== the scanned span); two overlapping different-period intervals merge into one `Repeat` region
carrying both; a repeat flanked by unique sequence yields `Unique, Repeat, Unique`; a repeat
longer than `max_repeat_len` yields `Satellite` (and a satellite spanning several windows
coalesces into one span); **window-count invariance** (identical regions across two `window_bp`
settings — the boundary-halo/attribution correctness test, the analogue of the cohort DUST
split-invariance gate); `merge_gap` bridges a small gap and `min_repeat_len` reclassifies a blip;
empty input yields nothing. The *catalog* consumer gets the golden-catalog parity + e2e tests
(§6).

**Open questions**

1. **Module placement (§4) — resolved:** `src/ng/tandem_repeat.rs` (a `src/ng/tandem_repeat/`
   directory only if it grows). Rejected `src/ssr/tandem_repeat.rs` (misrepresents scope,
   complicates reuse); supersedes an earlier top-level `src/tandem_repeat.rs` lean.
2. **Parity bar (§6) — resolved:** ≥ 99 % recall of golden loci with a reviewed, understood
   diff, not exact set-equality.
3. **Scoring defaults (§3.3) — resolved (starting value):** begin at `2 / 7` and **tune from the
   empirical evidence** gathered in our experiments. The default is a starting point, not a
   commitment; the parity fixture and downstream benchmarks drive the final weights.
4. **Scan range (§6) — resolved:** *not the module's decision.* The module default is
   `PeriodRange::new(1, 6)`; each consumer chooses its own — the catalog / ng snp-str caller
   uses `2, 6`. Both are consumer-side and may change after the experiments.
5. **Golden reference (§6) — resolved:** a subset of the tomato reference used downstream, for
   continuity with the e2e benchmark.

*All open questions are now resolved; the remaining empirical items (scoring weights, exact
scan ranges) are explicitly "tune during implementation from experiments," not open design
forks.*

---

## 10. Deferred, with a recommended home

- **Same-period stitch pass for the pathological split (§3.4).** Single indels already stay one
  interval; only a large/multi-indel tract splits. If a consumer ever shows material loss to
  that, merge nearby same-period intervals. Home: an opt-in pass in `tandem_repeat.rs` (a
  `stitch_gap` param), off by default so emission stays simple. Not built until measured.
- **Per-period scoring weights.** A single global ratio is the starting point; if a consumer
  needs a different tolerance per period, make the weights period-indexed. Home: `ScanParams`.
  Defer until measured.
- **Retaining phase-shifted interrupted tracts as one locus (catalog-side).** The scanner
  already emits the full interrupted span (§3.4); *keeping* the whole thing (both phases) as a
  single STR locus is a post-filter change — phase-aware purity / allele modelling, the
  direction of [`ssr_interrupted_repeat_recall.md`](../../specs/ssr_interrupted_repeat_recall.md). Home:
  `postprocess.rs` (or its successor), **not** the scanner — which is exactly why the scanner
  emits the full span rather than pre-trimming.
- **Removing the vendored `TRF-mod/` tree** once the scanner ships and the golden fixture is
  self-contained. Home: a repo-hygiene follow-up, not this change.
