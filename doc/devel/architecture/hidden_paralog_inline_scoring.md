# Scoring the hidden-paralog LR once, inside the main pass

**Status:** design, 2026-07-02, branch `tomato2-paralog-filter`. Follows
[hidden_paralog_single_sample_scoring.md](hidden_paralog_single_sample_scoring.md)
(which removed the last cohort-wide quantity from the score) and supersedes the
two-scoring-pass structure described in
[hidden_paralog_varcalling_wiring.md](hidden_paralog_varcalling_wiring.md) §S4/§S5.

**One-line change:** because the likelihood ratio no longer depends on anything
genome-wide, compute it **once**, during the main SNP-calling pass, and store it
in the spill — instead of computing it in the calibrate pass and *recomputing*
it in the write pass. This deletes a whole spill read, the separate window spill
+ join, and the "recompute-not-cache" machinery, and it moves the expensive work
onto the pass that is already parallel.

---

## 0. Why this is possible now (plain English)

The paralog score at a locus weighs two stories — one real variant (H1) vs. a
hidden collapsed paralog (H2) — using that locus's read depth and allele counts,
plus two things fixed for the whole run: each sample's coverage model (σ₀, GC
curve) and one cohort inbreeding number `F`. After the single-individual
reformulation, **that is the complete list** — there is no cohort expected
heterozygosity, no `Hexp`, no per-sample quantity accumulated across the genome.
The only genome-wide number left, π, is used to place the *cut*, not to compute
a locus's LR.

So a locus's LR can be computed the moment its own data exists. It does **not**
have to wait until the whole genome has streamed. The two-pass structure was
built when the LR still needed `Hexp` (a genome-wide accumulation); that
constraint is gone, so the structure can collapse.

---

## 1. What we do today (three passes, two scorings)

1. **Main pass** (parallel: producer pool + caller/EM workers). Calls every
   locus and, instead of writing the VCF, **spills** two streams: the full
   `PosteriorRecord` (record spill) and, separately, the per-window coverage
   (window spill). It cannot decide keep/drop yet — the cut needs π, which needs
   all loci.
2. **Calibrate pass** (serial). Reads the record spill, **joins** each locus to
   its window (by tile key), **computes the LR**, folds it into a fixed-size
   histogram; then EM over the histogram → π → FDR cut.
3. **Write pass** (serial). Reads the record spill again, **re-joins**,
   **recomputes the identical LR**, applies the cut, writes survivors.

The scoring — the expensive transcendental H1/H2 marginalisation — runs **twice**
(steps 2 and 3), **single-threaded**, after the 16-way main pass has finished.
That is the +11.6× wall regression (perf review 2026-07-02). The `-33 %`
memoisation already landed; this document removes the ×2 and the serialism at the
source.

---

## 2. The target: score once, in the main pass; histogram inline; one write pass

1. **Main pass.** As each locus is called, compute its LR once, **fold it into
   the histogram inline**, and spill `(record, LR)` in one stream. No separate
   window spill. The scoring is parallel (§3).
2. **After the main pass** (no spill read): the histogram is already complete →
   EM → π → FDR cut.
3. **Single write pass.** Read the `(record, LR)` spill once; apply the cut
   (`calibration.flags(LR)` — a pure curve lookup, no scoring); route survivors
   through the unchanged `VcfWriter`.

Two reads of the temp file → **one** (the write pass; still needed because the
cut isn't known until every LR is folded, and the records must be revisited to
emit them). Two scorings → **one**. And that one runs during the pass that
already owns all the cores.

---

## 3. Where the scoring runs (the load-bearing detail)

The LR needs two inputs that are produced in **different stages** of the main
pass:

- **window coverage** (each sample's mean depth over the locus's ~500 bp window,
  + the window's reference GC) — folded on the **main thread**, in the
  fold/plan stage, from `ChunkPlan.per_sample_segments` (per-position depth via
  `for_each_window_observation`). The caller workers never see per-position
  depth.
- **per-sample allele depth** (REF/ALT counts) — produced by the **caller
  worker's** per-group merger, inside `call_chunk`.

They **meet only after calling.** So the LR cannot be computed in the fold stage
(no AD yet) nor purely inside a caller worker (no window coverage). It is
computed where the two streams join, which is at/after the sink (the point that
already receives called records in genomic order). Two viable placements:

- **B — batch-and-score at the sink stage (recommended default).** The sink
  receives `CalledChunk`s (records + AD) and reorders them to genomic order. The
  main-thread window fold runs *ahead* of the callers, so by the time a record
  reaches the sink its window is closed and its coverage is available. Join
  record ⋈ window in memory (both monotone in coordinate → a bounded lockstep
  merge near the frontier — the window builder already keeps a bounded `pending`
  map), collect a batch, and **`rayon`-score the batch** (the pure, per-locus
  scorer), folding each LR into the histogram and spilling `(record, LR)`. The
  join is in memory, so **the window spill file disappears**. Scoring is parallel
  (a rayon batch), the sink thread only orchestrates.
- **A — score inside the caller worker (further optimisation).** Attach each
  chunk's *closed* variant-window coverages to the `RawCohortChunk` before
  dispatch (the fold stage has them), so the worker computes the LR right after
  it calls the chunk — riding the existing caller parallelism with no extra
  stage. Caveat: the window straddling the chunk's `cut` boundary is not closed
  at dispatch, so its variants' scoring must be **deferred** to the next chunk (a
  small carry between consecutive chunks). More plumbing; revisit only if B's
  sink-batch scoring shows up as a bottleneck.

Both compute the LR **once**, in parallel. B is simpler (all windows closed at
the join, no straddle deferral) and is the recommended first implementation; A is
a later refinement.

**Is it "parallel for free"?** Nearly. The window fold stays a cheap serial
main-thread step (summing depths, not transcendentals). The expensive LR runs in
parallel. It adds CPU work *during* the main pass rather than after it — but the
producer is decode-bound and often has spare cores, so scoring overlaps decode
stalls instead of running on an otherwise-idle machine afterwards. Even in the
worst case it is one parallel scoring instead of two serial ones.

---

## 4. What is removed / simplified

- **The window spill + `WindowJoin`** — gone. The record ⋈ window join happens in
  memory near the frontier (B) or by attaching coverage to the chunk (A); the LR
  it produces is what gets stored, so the coverage inputs need not be re-read.
- **The calibrate spill read** — gone. The histogram is folded inline during the
  main pass; after it, EM runs on the finished histogram with no read.
- **The write pass's LR recompute** — gone. The write pass reads the stored LR
  and does a pure curve lookup.
- **The "recompute-not-cache" invariant** (write_pass.rs §doc) — gone, and
  *strengthened*: the cut is applied to the exact LR the histogram was built
  from, byte-for-byte, instead of relying on the scorer being bit-reproducible
  across two independent evaluations.
- **`min_samples` / two-pass `hexp` threading** — already removed upstream.

The spill record grows by one `f64` (the LR; `NaN` = unscored → kept, matching
today's "non-finite LR is never flagged"). 8 bytes × loci ≈ 3 MB — negligible
against the multi-GB record spill.

---

## 5. Memory and byte-identity

- **Memory stays flat.** The histogram is a few KB. The record ⋈ window join
  buffer is bounded near the frontier (the window builder's existing `pending`
  bound). The rayon batch holds `N` records at once — the one new cost, the same
  RSS-vs-batch-size knob the perf review flagged; bound it by batch size. No
  genome-wide structure is retained.
- **Byte-identity holds.** Survivors still route through the unchanged
  `VcfWriter` in genomic order (the write pass feeds them post-reorder exactly as
  today). The histogram fold is order-independent integer bin counts, so π and
  the cut are identical regardless of scoring order. Storing the LR (vs.
  recomputing) makes calibrate↔write consistency *exact* by construction.
- **The reference-consistency guard** (write pass) is retained unchanged.

---

## 6. Open design points (for the implementation plan)

1. **Placement A vs B** — start with B (sink-batch + rayon), measure, consider A
   only if the sink stage bottlenecks.
2. **Batch size / thread budget** — the score batch competes with decode + EM for
   cores during the main pass. Pick a batch size that overlaps decode stalls
   without starving the callers; reuse / coordinate with the existing
   `resolve_split` budget rather than opening a third pool. RSS-vs-batch sweep is
   the gate.
3. **Spill record layout** — add the `LR: f64` field to `ParalogSpillRecord`;
   retire `WindowSpillRecord` / `WindowSpillBuilder` / `WindowJoin` once B lands.
4. **Histogram concurrency** — per-worker histograms merged at batch end (integer
   bins, associative), or one mutex-guarded histogram (contention is trivial at
   batch granularity). Prefer per-worker + merge.
5. **The un-scored classes** (not a SNP, no usable samples, no window) still spill
   with `LR = NaN` and are kept — one code path, no special "kept" spill flag.

---

## 7. Incremental implementation (byte-identity gated)

Each step keeps `--no-paralog-filter` byte-identical and the two-pass tests green
until the last step flips the structure.

1. **Add the `LR` field to the spill record** (written `NaN` for now; unused on
   read). Byte-identity of the off-path unaffected.
2. **Fold the histogram inline + write the real LR** in the main pass (placement
   B: join at the sink, rayon-score the batch). Keep the calibrate/write passes
   as-is but have calibrate read the *stored* LR instead of rescoring — proves the
   stored LR reproduces the old histogram (π / cut identical) before deleting
   anything.
3. **Delete the calibrate scoring** (histogram now comes from the main pass) and
   **the write-pass recompute** (read the stored LR). Retire the window spill +
   `WindowJoin`.
4. **Measure** end-to-end on tomato2 (wall + RSS) against the current two-pass +
   L1 baseline; confirm the T1 drop-profile is unchanged.

Step 2 is the substance; steps 1/3 are mechanical wrappers around it, and the
"calibrate reads the stored LR and must reproduce the old π/cut" check in step 2
is the correctness gate for the whole change.
