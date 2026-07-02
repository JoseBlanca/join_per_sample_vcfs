# Implementation plan — precomputed paralog window coverage in the pileup

Architecture:
[hidden_paralog_pileup_window_coverage.md](../architecture/hidden_paralog_pileup_window_coverage.md).
Branch `tomato2-paralog-filter`.

**Goal.** Compute the paralog window coverage + GC once in the pileup (a centred
sliding window per covered position), store them as two per-position `.psp`
columns, feed the coverage-by-GC histogram from the *same* sliding-window value,
and collapse the var-calling side to a per-locus lookup — deleting the window
accumulator, the window spill, and the join. Then score each locus once, in the
caller workers, and store the LR.

**Pacing.** Each milestone is independently testable and left green; pause
between milestones. Byte-identity notes are called out per step.

---

## Invariants to hold throughout

- **`--no-paralog-filter` stays byte-identical.** The new `.psp` columns are pure
  additions; the existing per-position record stream must round-trip unchanged,
  and the non-paralog VCF must not move.
- **The paralog drop set *may* shift** once scoring uses the new coverage (it is
  now `N`-excluded and computed on the same window the model is fit on). That is
  expected and is validated by profile coherence in M8, not by byte-identity.
- **Fit and score use the identical window** by construction (one accumulator).

---

## M1 — the sliding-window coverage accumulator (Stage-1 math)

Pure module work in [src/sample_summary/coverage.rs](../../../src/sample_summary/coverage.rs).

1. Settle the window convention: width `W = --gc-window-bp` (default 500),
   centred on `p`, over covered positions in `[p − W/2, p + W/2]` (pin the exact
   half-open/inclusive rule + odd/even `W` in the doc-comment). Fit and score
   share it because they share this accumulator.
2. Replace the fixed-tile fold with a **centred sliding window** that, per covered
   position, yields `(GC fraction, windowed mean coverage)` over the window's
   covered positions. Internally: a two-pointer deque of `(pos, depth, is_gc)`
   for positions within the current span; a position `p` is *finalised* when the
   stream reaches `p + W/2`.
3. Interface: the accumulator emits finalised `(pos, gc, mean_coverage)` values in
   order as the stream advances (plus a `drain()` for the tail at contig
   end/`finish`), **and** bins each finalised value into the existing
   `CoverageByGcHistogram`. So it serves both consumers (M3 writes the emitted
   values to columns; the histogram is unchanged in shape).
4. Reset per contig; only covered (non-`N`) positions enter the sum and divisor
   (same rule the tile accumulator used).
5. Rename `n_tiles` → `n_positions` on the histogram (now one sample per covered
   position, ≈ `callable_positions`); leave `callable_positions` (still the
   covered-position count = the het denominator).

**Tests:** the windowed mean/GC over hand-built streams (uniform, ramp, an
`N`-gap, a contig boundary, a short region edge, odd/even `W`); the finalise-lag
ordering; determinism. Histogram parity on a uniform stream (a flat sample still
lands in one GC/depth cell).

**Gate:** unit tests green. No `.psp`/pipeline change yet — this module is not
wired in.

---

## M2 — two per-position `.psp` columns

Format extension in the `.psp` writer + reader
([src/psp/](../../../src/psp/)). Pre-alpha: plain addition, no version
negotiation.

1. Add two per-position columns — `windowed_gc: f32` and `windowed_coverage: f32`
   — to the `.psp` body, decided as light vs heavy columns and `f32` vs a scaled
   integer (measure size; `f32` is the simple default).
2. Writer accepts the two values per record; reader surfaces them.
3. A record with no windowed value yet (shouldn't happen once M3 lags correctly,
   but define it) uses a `NaN` sentinel → treated as "no coverage evidence" →
   the sample contributes nothing at that locus (matches today's absent-window
   behaviour).

**Tests:** `.psp` round-trip of the two columns (incl. `NaN`); existing
record-stream round-trip **byte-identical** for the untouched fields.

**Gate:** psp round-trip tests green; the columns are written `NaN` by any
not-yet-updated caller (so nothing else breaks).

---

## M3 — wire the seam (pileup → columns + histogram)

[src/pileup/per_sample/pileup_to_psp.rs](../../../src/pileup/per_sample/pileup_to_psp.rs).

1. Replace `SampleSummaryAccumulators`' fixed-tile coverage accumulator with the
   M1 sliding-window one.
2. Add the **centred-window delay buffer**: hold pending `.psp` rows; as the
   accumulator finalises position `p`'s `(gc, coverage)`, attach them to the
   buffered row for `p` and write it (M2 columns). Output order unchanged, lagged
   by ≤ W/2.
3. Drain at region/contig/stream end so the tail rows get their (possibly
   truncated) windowed values.
4. Keep feeding the histogram from the same finalised values (M1 does the fold).

**Tests:** the seam writes the correct windowed values for a small stream
(cross-check against a direct recompute); the tail drains; region-clamp still
writes exactly the in-range rows with their windowed values; the existing seam
round-trip test still passes for the untouched fields.

**Gate:** pileup emits the columns end-to-end; existing pileup + seam tests green;
a freshly-generated `.psp` carries sane windowed values (spot-check on a fixture).

---

## M4 — coverage-model fit on the new histogram

[src/paralog/coverage_model.rs](../../../src/paralog/coverage_model.rs) —
expected to need **no code change** (it fits from the histogram, whose shape is
unchanged), but it must be validated on the new (sliding-window, per-position,
correlated, W×-denser) counts.

1. Confirm `fit` is robust to the count scale (mode = arg-max density, σ₀ =
   1.4826·MAD over the density) — correlated repeats don't bias mode/MAD, but
   pin it with a test on a synthetic sliding-window histogram.
2. Re-run the R1-style coverage-model validation on regenerated tomato2 `.psp`:
   σ₀ should land near the prototype value **but is expected to be smaller than
   the old tile-fit σ₀** (sliding windows are smoother) — record the new value;
   this is the yardstick the score now uses.

**Gate:** fit tests green; tomato2 σ₀ recorded and sane (a short note, not a full
report).

---

## M5 — var-calling reads the columns (parallel path in place, old path still live)

[src/var_calling/paralog_filter/](../../../src/var_calling/paralog_filter/).

1. Decode the two new columns in the producer alongside the existing per-position
   data, so at a variant locus each sample's `(windowed_gc, windowed_coverage)`
   is available in that row.
2. Point `build_observation`
   ([calibrate.rs](../../../src/var_calling/paralog_filter/calibrate.rs)) at the
   decoded columns instead of the joined window record: `relative_copy_number =
   coverage_model.relative_copy_number(windowed_gc, windowed_coverage)`.
3. **Do not delete the old window machinery yet.** Keep it running and add a
   **cross-check** (debug): the column-sourced `(gc, coverage)` must match the
   window-spill-sourced values within tolerance on regenerated `.psp` (they won't
   be bit-identical — different window definition — so assert *closeness* and log
   the distribution, or gate the check behind a flag). This proves the lookup is
   wired correctly before removing the fallback.

**Gate:** the paralog two-pass runs sourcing coverage from the columns; the
cross-check shows the lookup matches the intended window; all tests green.

---

## M6 — delete the window machinery

Once M5's lookup is trusted, remove the now-dead code and the window spill.

1. Delete `WindowMeanDepthAccumulator`
   ([window_coverage.rs](../../../src/var_calling/paralog_filter/window_coverage.rs)),
   `ReferenceWindowGc`
   ([window_gc.rs](../../../src/var_calling/paralog_filter/window_gc.rs)) — keep
   `reference_base_matches` / the coordinate guard if we want it,
   `WindowSpillBuilder`
   ([window_producer.rs](../../../src/var_calling/paralog_filter/window_producer.rs)),
   and `WindowSpillRecord` / `WindowSpillWriter` / `WindowSpillReader` /
   `WindowJoin` ([spill.rs](../../../src/var_calling/paralog_filter/spill.rs)).
2. Remove the fold-loop wiring (`feed_window_builder`, `begin_window_interval`,
   the `window_builder` state and its `finish()` in
   [pipeline.rs](../../../src/var_calling/pipeline.rs)) and the second spill file.
3. `score_spilled_locus` / `calibrate` / `run_write_pass` drop the `windows`
   join argument; the coverage comes from the record's columns.

**Gate:** clippy clean (no dead code); the two-pass filter runs coverage-from-columns
only; tests green; `--no-paralog-filter` byte-identical.

---

## M7 — score once, in the caller workers, store the LR

The payoff the whole line of work was aiming at (see architecture §6).

1. In `call_chunk` (the caller worker), after the allele counts exist, compute the
   paralog LR for each variant row using the per-sample `(gc, coverage)` columns +
   the up-front coverage model + `F` (all available to the worker). No join, no
   look-ahead — the window value is precomputed.
2. Store the LR in the spill record's `paralog_lr` field (already added, commit
   `10ee90e`) and fold it into a histogram (per-worker histograms merged, or a
   mutex-guarded one — integer bins, order-independent).
3. After the main pass: EM over the histogram → π → cut (no calibrate read).
4. Single write pass: read `(record, paralog_lr)`, apply the cut, write. Delete the
   calibrate scoring pass.

**Byte-identity:** survivors still route through the unchanged `VcfWriter` in
genomic order; the histogram fold is order-independent; storing the LR makes the
cut apply to the exact value it was built from.

**Gate:** two-pass filter runs with a single scoring, in parallel; tests green;
π/cut identical to an M6-era recompute on the same input (cross-check once, then
drop the recompute).

---

## M8 — validate on tomato2

1. Regenerate the 59-sample `.psp` with the columns.
2. Run the filter on/off; confirm the **T1 drop-profile** still shows the
   coverage-excess + het-excess signature, π is sane, and freebayes still emits
   the dropped loci (profile coherence — there is no truth set).
3. Measure wall + peak RSS vs the current two-pass + L1 baseline; the window
   spill and the double scoring are gone, so expect a substantial wall drop and
   the disk spill shrinking to the record spill alone.
4. Short impl report under `doc/devel/reports/implementations/`.

---

## Dependency order

```
M1 ─┬─> M3 ─┬─> M4
    │       └─> M5 ──> M6 ──> M7 ──> M8
M2 ─┘
```

M1 + M2 are independent and can land in either order; M3 needs both. M4 (fit
validation) and M5 (var-calling read) both follow M3; M5→M6→M7 is the collapse
sequence; M8 validates the whole.

---

## Notes / risks

- **Storage.** Two `f32`/covered-position/sample is the one real cost; measure the
  `.psp` size delta on a tomato sample in M2 and decide `f32` vs scaled-int then.
- **Region edges** truncate the window (fewer covered positions) — no worse than
  today's tiles; document at the column definition.
- **The old `hidden_paralog_inline_scoring.md`** and the window-spill parts of
  `hidden_paralog_varcalling_wiring.md` are superseded by this line of work;
  update their status once M6/M7 land.
- **`paralog_lr` field** (already committed) is reused by M7; until M7 the sink
  keeps writing `NaN` there, harmless.
