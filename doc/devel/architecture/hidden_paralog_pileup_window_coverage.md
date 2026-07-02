# Precomputing the paralog window coverage in the pileup

**Status:** design, 2026-07-02, branch `tomato2-paralog-filter`.

**Supersedes** the window-coverage machinery of
[hidden_paralog_varcalling_wiring.md](hidden_paralog_varcalling_wiring.md)
(Approach A / the S3 + S6c window spill) and largely dissolves
[hidden_paralog_inline_scoring.md](hidden_paralog_inline_scoring.md) (the reason
that document is complicated goes away here). The score itself
([hidden_paralog_locus_statistic.md](hidden_paralog_locus_statistic.md),
[hidden_paralog_single_sample_scoring.md](hidden_paralog_single_sample_scoring.md))
is unchanged.

---

## 0. The idea in one paragraph

The paralog score needs, at each variant locus and for each sample, the **mean
read coverage in a window around the locus** and the window's **GC fraction** —
the "is this region covered about twice as much as normal?" signal. Today the
var-calling phase recomputes that window coverage on the fly, which is genuinely
hard there because it needs to look *ahead* of the current position while running
inside a chunked, memory-optimised, multi-threaded pipeline. But that number is
trivial to compute in the **pileup** phase, where each sample is a single clean
position-by-position stream. So: compute the windowed coverage and GC **once, in
the pileup**, and store them as two per-position columns in the `.psp`. The
var-calling phase then does no window work at all — it just **looks up** the two
values sitting next to each variant.

---

## 1. Why the pileup, not var-calling

To get "coverage in a window centred on position *p*" you must average over
positions on *both sides* of *p* — you have to see positions after *p*.

- **In var-calling that look-ahead is expensive.** The genome is split into
  chunks handed to parallel worker threads; the per-position depth is decoded by
  a front stage; windows straddle chunk boundaries; a window isn't final until
  the stream has passed its far edge. Every one of the awkward options we
  explored (a scratch file joined later; a live cross-thread channel; scoring in
  the workers with a "carry" for boundary windows) exists only to work around
  this look-ahead.
- **In the pileup that look-ahead is trivial.** A sample's pileup is one
  position-by-position stream from start to finish. "Centred on *p*" just means
  keeping a small buffer of the next ~½-window of positions — a few hundred
  entries in a straight loop, entirely inside one sample. No chunks, no threads,
  no boundaries.

The pileup is where this computation belongs. Moving it there removes the root
complexity instead of engineering around it.

---

## 2. What is computed and stored

For every **covered** reference position *p* of a sample, the pileup computes a
**centred sliding window** of width `W` (the existing `--gc-window-bp`, default
500 bp → *p* ± 250) over the covered positions in `[p − W/2, p + W/2]`:

- **windowed mean coverage** = Σ(per-position depth) / (covered positions in the
  window);
- **windowed GC fraction** = (G/C reference bases in the window) / (covered
  positions in the window).

Both are exactly the quantities the current per-tile accumulator already produces
([coverage.rs:218-219](../../../src/sample_summary/coverage.rs#L218-L219)) — the
only change is that the window is a *sliding, centred* window per position rather
than a *fixed, non-overlapping* tile, and the value is **kept** rather than only
histogrammed.

These two values become **two new per-position columns** in the `.psp` body,
written next to each existing per-position row. (Pre-alpha: a plain format
extension, no version negotiation.)

GC is stored as a column too — even though it is a reference property and nearly
shared across samples — precisely because storing it means the var-calling phase
needs to know *nothing* about windows or the reference for scoring: it reads two
numbers per locus and is done. Recomputing GC from the reference in var-calling
would drag a reference-streaming step (and a subtly different GC definition) back
into the phase we are trying to simplify. The small redundancy buys a clean
downstream.

---

## 3. One computation, two consumers

A single **sliding-window accumulator** runs over the sample's per-position
stream. At each position it produces one `(GC, mean coverage)` value, which is
used in two places:

1. **The coverage-by-GC histogram** (unchanged in shape and role) — the value is
   binned exactly as today, and the histogram is still what
   [`SingleCopyCoverageModel::fit`](../../../src/paralog/coverage_model.rs) reads
   to estimate the sample's typical single-copy coverage, its scatter **σ₀**, and
   the GC-bias curve. The only difference: the samples fed in are now
   sliding-window values, not fixed-tile values.
2. **The per-position `.psp` columns** — the same value is written out, so a
   variant at that position can read it back later.

Because the **yardstick** (σ₀, fit from the histogram) and the **observation**
(the stored column, looked up at a locus) come from the *identical* window
definition, they are automatically on the same scale. That is the consistency the
score needs: "is this coverage surprising?" is judged with the same ruler it was
measured with. (Fitting σ₀ on 500 bp tiles while scoring on sliding windows would
mis-calibrate it — sliding windows are smoother, so their scatter is smaller.)

Note on counting: the histogram now receives one sample per covered **position**
rather than per **tile** (≈ W× more, and neighbouring windows overlap so the
samples are correlated). This does not bias the fit — a density estimate's mode
and scatter are unaffected by correlated repeats — but the bookkeeping field
`n_tiles` becomes, in effect, `n_positions`.

---

## 4. Pileup-side mechanics

The per-position stream is the pileup→psp seam
([pileup_to_psp.rs:55-76](../../../src/pileup/per_sample/pileup_to_psp.rs#L55-L76)):
`observe_record` already computes each covered position's depth and its reference
base. The sliding-window accumulator is fed that same `(pos, ref_base, depth)`
stream.

The one new mechanic is the **centred-window delay**. To finalise position *p*'s
value we need positions up to *p* + W/2, but the seam wants to write *p*'s row
when the walker emits it. So the seam holds a bounded buffer of pending rows and
emits row *p* (now carrying its finished windowed values) only once the stream
reaches *p* + W/2. Concretely, a two-pointer sliding window over the
covered-position stream: a deque of `(pos, depth, is_gc)` for positions within the
current window span, plus the pending rows awaiting their far edge. Memory is
`O(W)` positions; output order is unchanged (just lagged by ≤ W/2).

Edge cases (all pre-existing in spirit, none new to the model):

- **Contig boundaries** — windows do not span chromosomes; the accumulator resets
  per contig (as the tile accumulator does today).
- **Analysis-region (`--regions`) edges** — a window near a region boundary
  averages only the covered positions that were actually pileup'd, so it is
  truncated there. The current tile approach truncates identically; this is not a
  regression. Document it where the column is defined.
- **Uncovered gaps** — only covered positions enter the sum and the divisor, so a
  window straddling an uncovered gap behaves exactly like the tile accumulator's
  "GC-defined covered positions" rule.

---

## 5. Var-calling-side simplification

With the columns present, the var-calling paralog path loses all of its window
machinery. **Deleted:**

- `WindowMeanDepthAccumulator`
  ([window_coverage.rs](../../../src/var_calling/paralog_filter/window_coverage.rs))
  — the per-sample window-depth re-derivation;
- `ReferenceWindowGc`
  ([window_gc.rs](../../../src/var_calling/paralog_filter/window_gc.rs)) — the
  reference GC walk (for scoring; the reference-base consistency guard can stay if
  wanted);
- `WindowSpillBuilder`
  ([window_producer.rs](../../../src/var_calling/paralog_filter/window_producer.rs))
  and the whole **window spill** — `WindowSpillRecord`, `WindowSpillWriter/Reader`,
  and `WindowJoin`
  ([spill.rs](../../../src/var_calling/paralog_filter/spill.rs));
- the fold-loop wiring that feeds them (`feed_window_builder`,
  `begin_window_interval`), and the tile-key join in the calibrate/write passes.

**Replaced by:** a per-locus lookup. Each sample's `.psp` already streams through
the producer; the two new columns decode alongside the existing per-position data,
so at a variant locus each sample's `(GC, mean coverage)` is simply *there*, in
that position's decoded row. `build_observation` reads them directly instead of
joining a window record.

---

## 6. What this unblocks: scoring once, in parallel, the easy way

The whole reason the "score once, in parallel" question was hard is that one of
the score's two inputs (the window coverage) was not ready when and where the
record was. Once it is a column sitting next to the record, that difficulty
evaporates:

- **The caller worker already holds everything it needs.** The compacted chunk a
  worker calls carries each variant row's per-sample data, now including the
  `(GC, coverage)` columns. After the worker computes the allele counts, it has
  both score inputs — so it can compute the paralog LR right there, riding the
  workers' existing parallelism, with **no** cross-thread join and **no**
  straddling-window carry (the window mean is precomputed per position). This is
  the clean form of "Option 1" from the inline-scoring note, now with its one
  caveat removed.
- The LR is stored once (the `paralog_lr` field already added to the spill
  record) and folded into the histogram; after the main pass, EM → cut; a single
  write pass applies it. No calibrate re-scoring, no window join.

So the sink-channel gymnastics of
[hidden_paralog_inline_scoring.md](hidden_paralog_inline_scoring.md) are no longer
needed. That document's *goal* (score once, keep the value) stands and is even
easier here; its *mechanism* is superseded.

---

## 7. Costs and trade-offs

- **`.psp` size:** two `f32` per covered position per sample. This is the same
  order as the per-position data the `.psp` already stores, and the windowed
  values are smooth, so they compress well (delta / downsample if it ever
  matters). Judged acceptable.
- **Pileup work:** one extra streaming accumulator + a bounded reorder buffer in
  the seam. Cheap relative to BAQ / the walker.
- **`.psp` regeneration:** existing artefacts must be regenerated to gain the
  columns — but they already must be regenerated for the coverage summary
  section, so this rides along. (See the summary-missing hard-error already in
  place.)
- **A subtle correctness gain, not a cost:** the stored mean excludes `N` and is
  computed on exactly the window the model is fit on, removing the documented
  "var-calling depth doesn't N-exclude" approximation and the fit/score GC
  mismatch.

---

## 8. Open design details (settled in the plan)

1. **Exact window convention** — `W` = `--gc-window-bp`; centred half-open
   `[p − W/2, p + W/2)` vs inclusive `± W/2`. Pin one; keep fit and score
   identical by construction (same accumulator).
2. **Column placement / encoding in the `.psp`** — light vs heavy column, `f32`
   vs a scaled integer, and how the reader surfaces them to the var-calling
   decode.
3. **Whether the reference-base consistency guard** (write pass) stays; it no
   longer depends on the GC walk but is a cheap coordinate-drift check.
4. **Histogram field rename** (`n_tiles` → `n_positions`) and any downstream that
   reads it.
