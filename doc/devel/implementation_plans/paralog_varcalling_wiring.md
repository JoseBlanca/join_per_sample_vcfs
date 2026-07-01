# Hidden-paralog filter — var-calling wiring — implementation plan

**Status:** draft, 2026-07-01, branch `tomato2-paralog-filter`. Turns the settled
wiring architecture
([hidden_paralog_varcalling_wiring.md](../architecture/hidden_paralog_varcalling_wiring.md))
into build order. This is the **detailed** form of Milestone S in the parent plan
[paralog_filter_model.md](paralog_filter_model.md) (which sketched S1–S5 before
the pipeline was mapped); it supersedes that sketch.

The pure statistics core (Q1–Q5) and its data-first validation (R1) are done and
committed. This plan is the production integration: read the per-sample summaries,
score every candidate variant, calibrate the prior + FDR from a streaming
histogram, and drop the paralog-suspect calls from the VCF — all with memory flat
in both variant count and sample count.

## Domain intent

`var-calling` currently streams `.psp` → VCF in one pass. The paralog filter can't
live in that pass, because a locus's verdict depends on every other locus (`Hexp`
needs all allele frequencies; π + the FDR cut need all LRs). So the run becomes:
score every candidate variant, spill it to an ephemeral binary file while
accumulating the global quantities, calibrate once at the end, then read the spill
back and write the surviving calls. The design decisions are settled in the
architecture doc; this plan builds them.

## The settled decisions this plan implements (do not reopen)

- **One `.psp` body read.** Per-window `(gc, mean depth)` is aggregated **inline in
  the existing fold** — the per-allele observation counts are already a phase-1
  "light" column, so total depth (`ref_obs + nonref_obs`) is free for every
  position; GC comes from the REF span the fold already fetches. No second body
  pass, no schema change (arch §2, §4).
- **`Hexp` is accumulated during the main caller pass** from
  `PosteriorRecord.allele_frequencies`, on the **per-callable-position scale**
  (`Σ2pq / callable`, *not* mean over variant sites — the R1 fix; arch §6).
- **Two spill reads, never cached.** Both spill passes stream from disk holding one
  locus at a time; the write pass **recomputes** each LR (bit-identical, since the
  scorer is pure) rather than carrying the calibrate pass's LRs (arch §3).
- **EM on the histogram, not the data.** The streaming pass folds each LR into the
  fixed-size histogram (no π needed); batch EM runs on the ~2000 bins afterward.
  A fixed-prior **fallback** is wired to Q5's `converged` flag (arch §6).
- **Hard removal, on by default.** Flagged loci are dropped (the `emit_or_drop`
  convention, a `records_dropped_paralog` stat); `--paralog-fdr` defaults to ≈ 1%,
  `--no-paralog-filter` disables and restores the single-pass path (arch §7, §8).
- **Spill = reuse the `.psp` block-writer** infrastructure; ephemeral, deleted on
  success and failure (arch §5).

## Principles (how the order was chosen)

- **Types first, then implementation, within every step; pause between steps.** The
  standing project + incremental rules.
- **Standalone pieces before the integration.** S1–S5 each build and unit-test a
  piece against fixtures with no change to `run_var_calling`'s control flow; S6 is
  the single integration step that rewires the pipeline. This keeps the risky
  hot-path and control-flow change isolated and late.
- **Byte-identity is the safety net.** The existing streaming path must stay
  byte-identical when the filter is off; every step that touches shared code
  proves that before adding the new behaviour.
- **Memory flatness is a per-step obligation, not just T2.** Every step that
  handles per-locus data holds one locus (or one window) at a time; a step that
  would retain a genome-wide structure is wrong even if it passes its test.

---

## Milestone S — the wiring

### S1. `ParalogPrePass` — per-sample models + `obs_het` + the `Hexp` accumulator ☐

The up-front, per-sample state, plus the running accumulator for the one global
quantity that comes from the caller pass.

- **Types first:** `ParalogPrePass` holding, per sample, the fitted
  `SingleCopyCoverageModel`, `obs_het`, and `callable_positions`; and a
  `HexpAccumulator` that folds a locus's cohort allele frequencies into `Σ2pq` and
  finalises to `Hexp = Σ2pq / callable_ref` (a representative callable count, e.g.
  the cohort median). `F` is *not* set yet — it needs the finished `Hexp`.
- **Build:** a pure function that, given the open `PspReader`s' metadata sections
  (coverage histogram + het counts), fits each model and forms `obs_het`
  (Q2 + Q4); and the `HexpAccumulator::observe(&allele_frequencies)` /
  `finish(callable_ref) -> Hexp`.
- **Deliverable:** the pre-pass built for the tomato2 metadata; a test that
  `HexpAccumulator` matches an independent `Σ2pq` sum on a fixture, and that a
  degenerate sample (coverage fit rejected) is carried as absent, not fatal.
- *Depends:* Q2, Q4. *Source:* arch §3 step 1, §6.

### S2. The spill file — streaming writer + reader ☐

The ephemeral per-locus store, written once and read twice.

- **Types first:** `ParalogSpillRecord` (position, alleles, per-sample GT/GQ/AD,
  cohort AF, per-sample window `(gc, mean depth)`, existing FILTER, QUAL — enough
  to both score the locus and write its VCF record); `ParalogSpillWriter`
  (streaming append, block-buffered) and `ParalogSpillReader` (one-pass
  read-back), built on the `.psp` block-writer infrastructure (`psp/writer.rs`).
- **Build:** append + block-flush + a read-back iterator; the file lives in scratch
  and is deleted on drop / on run end (success and failure).
- **Deliverable:** a round-trip test (write N records, read them back identical);
  a test that the reader holds one record at a time (no full materialisation); the
  temp file is registered for cleanup and gone after the reader is dropped.
- *Depends:* none (self-contained). *Source:* arch §5.

**Checkpoint (after S2):** the two standalone containers — the pre-pass and the
spill — exist and are tested in isolation. Pause.

### S3. Inline per-window coverage in the fold ☐

Aggregate each sample's 500 bp-window `(gc, mean depth)` inline in the existing
fold, without a second body read and without breaking the streaming path.

- **Types first:** a `WindowCoverageAccumulator` — the O(1) running state
  (`covered count`, `G/C count`, `depth sum`, open-window key) that mirrors the
  producer's `CoverageByGcAccumulator` tiling but retains per-window values;
  emits a window's `(gc, mean depth)` when the coordinate-ordered stream crosses
  the boundary.
- **Build:** feed it from the fold's light columns (`depth = ref_obs + nonref_obs`,
  both light; GC from the fetched REF span), one sample at a time; attach a closed
  window's `(gc, mean depth)` to the variant loci inside it (a small per-window
  buffer while the window is open).
- **Deliverable:** on a fixture, the inline per-window `(gc, mean depth)` matches an
  independent per-window computation and reproduces the producer's tiling (so it is
  consistent with the model's fit); **the existing var-calling output is
  byte-identical** when the accumulator runs but its result is unused (proves the
  hot-path change is side-effect-free).
- *Depends:* S2 (spill is where the window value lands). *Source:* arch §2, §4.

### S4. Calibrate pass — spill → LR → histogram → π + FDR cut ☐

The first spill read: score every locus and derive the global calibration.

- **Types first:** `ParalogCalibration { prior: ParalogPrior, curve:
  ParalogFdrCurve, lr_threshold: f64 }`, and a `calibrate(prepass, hexp, spill,
  params, fdr_target) -> ParalogCalibration`.
- **Build:** stream the spill; per locus turn window depth → relative copy number
  (Q2) and build `LocusObservations` (Q3 inputs) using the pre-pass's `F` (formed
  here from the finished `Hexp`), compute the LR, **fold it into the fixed-size
  histogram and discard the locus**. Then `ParalogPrior::estimate` (EM on the
  histogram) → π, with the **fixed-prior fallback** on `!converged`; build the
  `q_of_lr` curve and resolve the LR threshold for `fdr_target`.
- **Deliverable:** on a fixture spill, produces π + threshold with **no per-locus
  retention** (asserted); a test that π + threshold match a full-vector reference
  within bin tolerance; a test that the non-convergent path takes the fallback and
  warns.
- *Depends:* Q3, Q5, S1, S2, S3. *Source:* arch §3 step 5, §6.

**Checkpoint (after S4):** given a spill, the whole calibration (π, FDR cut) is
produced RAM-flat. The statistics are wired to the spill. Pause.

### S5. Write pass — spill → recompute LR → apply cut → VCF ☐

The second spill read: emit the surviving calls.

- **Types first:** the write-pass driver + a `records_dropped_paralog` counter in
  the writer stats; the new FILTER-free header provenance (an `##INFO`/`##command`
  or header line recording the target FDR, π, and the resolved LR cut).
- **Build:** stream the spill a second time; per locus **recompute the LR** from the
  spilled inputs (bit-identical to S4), look up `q_of_lr(LR)`; if `≤ threshold`,
  **drop** the record (increment the counter) — otherwise write it (mirroring the
  allele-balance drop in `vcf_writer.rs::emit_or_drop`). Write to the temp VCF path;
  on finish, atomic-rename and delete the spill (also on failure).
- **Deliverable:** on a fixture spill + a fixed calibration, the surviving VCF has
  exactly the expected records dropped; header records the provenance; the spill is
  gone afterward; a determinism test (same inputs → identical VCF, since the LR
  recompute is pure).
- *Depends:* S4; the VCF writer (`vcf_writer.rs`, `vcf/header.rs`). *Source:* arch
  §3 steps 6–7, §7.

### S6. Orchestration + CLI — wire the two-pass flow into `run_var_calling` ☐

The single integration step: make `run_var_calling` run the two-pass flow when the
filter is on, and the current single pass when it is off.

- **Types first:** the `--paralog-fdr <q>` (default ≈ 0.01) and
  `--no-paralog-filter` flags on `CohortPipelineArgs` (`cli/shared_args.rs`,
  following the allele-balance flags), threaded into the pipeline config.
- **Build:** when the filter is on, the main caller pass writes to the spill
  (S2) instead of the VCF while accumulating `Hexp` (S1) and the inline window
  coverage (S3); then form `F`, run the calibrate pass (S4) and the write pass (S5).
  When off (`--no-paralog-filter`, or `--paralog-fdr 0`), the pipeline keeps its
  current single-pass, direct-to-VCF behaviour untouched.
- **Deliverable:** an end-to-end `var-calling` run on tomato2 emits a filtered VCF;
  the temp spill is gone afterward; **the existing byte-identity / determinism
  integration tests pass with the filter off** (re-baselined to pass
  `--no-paralog-filter`), and **new tests** cover the filtered path (a locus is
  dropped, π + cut recorded in the header, a non-default `--paralog-fdr` changes the
  drop count). Flag parsed/validated/defaulted; help text present.
- *Depends:* S1–S5. *Source:* arch §1, §7, §8; `pipeline.rs`, `cli/shared_args.rs`.

**Checkpoint (after S6):** a real `var-calling` run applies the filter end-to-end
with bounded RAM and a cleaned-up spill — the first production artefact. Pause.

---

## Milestone T — end-to-end validation + cost

### T1. End-to-end behaviour on tomato2 ☐

Confirm the production wiring reproduces the isolated-statistics result (R1). The
FILTER'd callset's **dropped** set should carry the paralog profile (het excess +
coverage excess), planted/known introgression-like loci should survive, and π
should land ≈ 9%, matching the R1 numbers.

- **Deliverable:** a short report in `doc/devel/reports/implementations/`, compared
  against the R1 validation.
- *Depends:* S6. *Source:* spec §7; the R1 report.

### T2. Cost / RSS check — the gate ☐

Because the filter is **on by default**, its cost is paid on every run, so this is
a **gate, not a footnote**. Confirm memory stays flat in variant and sample count
(no genome-wide LR vector; the spill on disk, the histogram a few KB, the window
accumulator O(1) per sample), and record the wall/RSS delta and the spill size on a
representative run vs the filter-off path.

- **Deliverable:** before/after numbers in the T1 report; a fix only if RAM is not
  flat or the delta is unacceptable for a default-on path.
- *Depends:* S6. *Source:* arch §9; [[low_memory_mode_branch]].

---

## Checkpoints

- After **S2**: the pre-pass and the spill exist, tested in isolation. Pause.
- After **S4**: given a spill, the full calibration (π + FDR cut) is produced
  RAM-flat. Pause.
- After **S6**: a real `var-calling` run applies the filter end-to-end, bounded
  RAM, spill cleaned up. Pause.
- After **T2**: cost confirmed acceptable for a default-on filter.

## Risk register (what to watch)

- **S3 hot-path change** — the fold is the pipeline's throughput-critical stage;
  the byte-identity-when-unused test is the guard. Keep the accumulator's cost off
  the critical path (it is O(1) per position).
- **S6 control-flow change** — the two-pass path is a materially different shape
  from the current single pass; keep the off path a literal no-op over today's
  code, so byte-identity is trivially true.
- **On-by-default fallout** — existing integration tests re-baseline; T2 must clear
  before this lands as the default (the owner can still choose to ship it
  off-by-default if T2 disappoints — that is a one-line default change).

## Open items carried to implementation

- Spill framing detail (columns / block layout) — settled in S2 against the
  block-writer API.
- LR-histogram bin resolution — confirm/tune against the tomato2 LR range in S4.
- Exact default `--paralog-fdr` value — pinned in S6 against the T1 flagged-set
  profile.
- Optional: emit the paralog LR/posterior as an `INFO` on *surviving* records for
  downstream inspection (does not change the callset) — an S5 option.
