# Wiring the hidden-paralog filter into var-calling

**Status:** draft, 2026-07-01, branch `tomato2-paralog-filter`. Settles the
data-flow architecture for **Milestone S** of the implementation plan
([paralog_filter_model.md](../implementation_plans/paralog_filter_model.md)) —
how the pure statistics pieces (built and validated in Q1–Q5, R1) plug into the
cohort caller `var_calling::pipeline::run_var_calling` and turn into a `FILTER`
verdict in the VCF.

This is the piece [`hidden_paralog_locus_statistic.md`](hidden_paralog_locus_statistic.md)
Premise 6 deferred as "sketched, not this slice." The statistics functions'
input types are exactly what the caller can produce (that much was settled), but
the *data flow* — when each global quantity becomes available, where the
per-window coverage comes from, how we keep memory flat — was not, and that is
what this doc resolves against the pipeline as it actually exists.

---

## 0. Glossary

Terms used throughout, in plain language (formulas come in the sections that
use them):

- **Locus** — one candidate variant site (a position with at least one
  non-reference allele in the cohort).
- **`Hobs`** (observed heterozygosity) — how often a sample is actually
  heterozygous, as a *rate*: heterozygous sites divided by all callable
  positions. A per-sample number.
- **`Hexp`** (expected heterozygosity) — how often a sample *would* be
  heterozygous if the population were mating at random, given the cohort's
  allele frequencies. One number shared by the whole cohort.
- **`F`** (inbreeding coefficient) — how much less heterozygous a sample is than
  random mating predicts: `F = 1 − Hobs/Hexp`. Near 0 for an outbred sample,
  near 1 for a selfer. Feeds the genotype prior in the score.
- **Relative copy number** — a window's read depth expressed in copies, `1.0 =
  one copy`, after dividing out the sample's sequencing depth and GC bias (the
  coverage model of Q2).
- **LR** (likelihood ratio) — the per-locus paralog score: `log P(data | hidden
  paralog) − log P(data | real variant)`. Positive favours a hidden paralog.
- **π** (prior paralog fraction) — the estimated fraction of loci genome-wide
  that are hidden paralogs. *One number for the whole run.* Says how *common*
  paralogs are, not *which* loci are paralogs.
- **Posterior** — the probability that *this* locus is a hidden paralog:
  `sigmoid(LR + logit π)`.
- **q-value / FDR** — the false-discovery rate: if we flag every locus with a
  score at least this high, the q-value is the fraction of the flagged set we
  expect to be real variants wrongly flagged. The **cut is chosen on the
  q-value**, so we control how many real variants we sacrifice.
- **Spill file** — a temporary scratch file we write once and read back, used to
  hold per-locus data on disk so it does not have to sit in RAM. Deleted when
  the run ends.

---

## 1. Why the verdict cannot be decided in one streaming pass

Today `run_var_calling` is a single streaming pass: it reads the per-sample
`.psp` files, calls genotypes, and writes each VCF record as it goes, holding
only a small window of the genome in memory at a time. That is what keeps its
memory flat as the number of samples grows.

The paralog filter cannot fit inside that single pass, because **its verdict for
one locus depends on every other locus**. Two of the quantities it needs are
*global* — they are not known until the whole genome has been seen:

- **`Hexp`** (needed to form each sample's `F`, which the score uses) is built
  from the cohort's allele frequencies across *all* loci.
- **π** and the **FDR cut** (needed to turn a locus's LR into a keep/flag
  decision) are built from the distribution of *all* loci's LRs.

So a locus's LR needs `Hexp` (all loci seen once), and a locus's verdict needs
π (all LRs computed). There is no way to decide the first locus's verdict before
the last locus has been read. This forces the classic shape: **stream the
per-locus data out to a spill file, compute the global quantities once at the
end, then read the spill back to write the final VCF.** The rest of this doc is
the concrete form of that shape.

When the filter is switched **off** (the default until it is proven in
production, and always available via the CLI), none of this applies: the
pipeline keeps its current single-pass, direct-to-VCF behaviour untouched.

---

## 2. What the pipeline already hands us

The caller was mapped against the paralog filter's needs. Three facts make the
wiring a construction step rather than a redesign, and one is the hard part.

**Good: the per-sample summaries are available up front.** All `PspReader`s are
opened together near the top of `run_var_calling` (`pipeline.rs:235`), before
any records stream. Each sample's `.psp` metadata section — the coverage-by-GC
histogram and the het counts (`n_het_sites`, `callable_positions` from P1) — is
readable there. So the per-sample **coverage model (Q2)** and **`obs_het`
(Q4)** are fit once, up front, at no streaming cost.

**Good: the cohort allele frequency is already computed per locus.** The record
that reaches the VCF writer, `PosteriorRecord`, carries `allele_frequencies`
(the EM's `p̂` per allele) plus each sample's read counts (`scalars`,
i.e. per-allele `num_obs`) and best genotype. So **`Hexp = Σ 2p(1−p)` can be
accumulated for free** as records stream past — exactly what arch Premise 3
wanted ("inside the existing per-locus genotype pass, not a separate pass"). No
new allele-frequency computation is needed.

**Better than feared: the depth is already in the light columns.** The score
needs, per locus per sample, the sample's **500 bp-window relative copy
number**, which needs the window's mean read depth — an average over *all*
covered positions in the tile, overwhelmingly non-variant. It first looked as
though the pipeline discards that (it decides which positions are variant early
and defers the per-allele detail). But the per-allele **observation counts are a
phase-1 "light" column** (`is_light_tag` = `{0x01, 0x02, 0x03, 0x10}`, where
`0x10` is the obs count) — they have to be, because finding variants needs them.
Only the per-allele *statistics* (quality, strand, MAPQ, sequence bytes, chain
ids) are heavy. So **total depth = `ref_obs + nonref_obs` is a single addition
from data the fold already decodes, for every position including monomorphic
ones.** The window's *mean* still has to be aggregated (§4) — but the
ingredients are already streaming past, at no extra decode and no schema change.

---

## 3. The pass structure (mapping onto the eight-step flow)

The recommended structure, annotated with which of the operator's eight steps
each part covers:

1. **Set-up, at reader-open (steps 1 partial).** Read each sample's `.psp`
   metadata → fit the coverage model, form `obs_het`. Cheap, `O(1)` per sample.
   `Hexp` is *not* known yet (it needs the whole genome) — only `obs_het` and
   the models are ready here.

2. **Main pass — the existing caller, redirected to the spill (steps 2 partial,
   3).** Run the pipeline as today (produce `PosteriorRecord`s), but instead of
   writing the VCF, **spill each locus** to the temp file: its position, alleles,
   per-sample genotype/GQ/AD, and cohort allele frequencies — everything the
   final VCF will need. While spilling, accumulate `Hexp = Σ 2p(1−p)`. The LR is
   **not** computed here: it needs `F`, which needs the finished `Hexp`.

3. **Per-window depth — folded into the main pass, not a separate pass (step 2
   partial; §4).** As the main pass streams each sample's light columns in
   coordinate order, an O(1) running accumulator tiles them into 500 bp windows
   (`depth = ref_obs + nonref_obs`, both light; GC from the reference bases the
   fold already fetches) and emits each window's `(gc, mean depth)`, attached to
   the variant loci in it and spilled. No second read of the `.psp` bodies.

4. **Form `F` (steps 1 partial).** With `Hexp` complete: `F_s = 1 −
   obs_het_s / Hexp` for every sample. **`Hexp` is on the per-callable-position
   scale** — `Σ2pq` divided by the callable-position count, *not* a mean over
   variant sites (the bug R1 caught; see §6). One-time, tiny.

5. **Score + calibrate pass — stream the spill once (steps 2 finish, 4).** Read
   the spill locus by locus, streaming from disk. For each locus turn window
   depth into relative copy number (the coverage model), build the observations,
   compute the **LR** (Q3, now that `F` exists), **fold that LR into the
   fixed-size LR histogram (Q5), and discard the locus** — nothing per-locus is
   retained across iterations. The histogram is a few kilobytes; there is never a
   genome-wide vector of LRs or of spilled records in RAM. From the finished
   histogram, estimate **π** and build the **`q_of_lr` curve**, then resolve the
   LR threshold for the CLI's target FDR (step 4's "how many / which" — see §6).

6. **Write pass — stream the spill again (steps 5, 6, 7).** Read the spill a
   second time, again streaming from disk. For each locus, **recompute its LR
   from the same spilled inputs** (`score_locus_for_paralogy` is pure and
   deterministic, so this is bit-identical to the calibrate pass — the applied
   cut matches the histogram it came from), look up its q-value on the curve, and
   either write the record or drop it (§7). Write the VCF to a temp path, one
   record at a time.

7. **Finish (step 8).** Atomic-rename the finished VCF to the destination,
   delete the spill (on success *and* on failure), and report π, the cut, and
   the flagged count to the operator.

The spill is read **twice** (once to score, once to write) because π sits between
them: no locus's verdict exists until every LR is in the histogram. That second
read is the price of a global verdict, and it is a disk read, not RAM.

> **The non-negotiable invariant (do not "optimize" it away).** Neither pass
> caches the spill in memory. The calibrate pass streams the spill, folds each LR
> into the fixed-size histogram, and drops the locus; the write pass **re-reads
> the spill from disk** and **recomputes** each LR rather than carrying the
> calibrate pass's LRs forward. Holding the LRs (or the records) from the first
> pass to avoid the second read would make RAM grow with the variant count —
> exactly what the spill exists to prevent. The spill read is fast; the second
> read is cheap; a genome-wide in-RAM vector is not. Recompute is safe because
> `score_locus_for_paralogy` is a pure function of its spilled inputs, so the
> write pass's LR is bit-identical to the calibrate pass's.

> **A note on step ordering.** Steps 2 and 3 are the *same* pass: the window
> depth (step 3) is aggregated inline in the main pass's fold (§4), so the `.psp`
> bodies are read once. Steps 5 and 6 are the two reads of the *spill*, not the
> bodies.

---

## 4. The per-window coverage — computed inline in the fold

The score needs each locus's per-sample window depth. §2 established the good
news: the per-position depth is already in the phase-1 light columns, so it
streams past the fold at no extra decode cost. The window *mean* still has to be
aggregated; the memory-scaling thesis (RAM flat in *sample* count) decides how.

**Recommended — aggregate the window inline in the main-pass fold.** The fold
already visits every covered position of every sample in coordinate order. We add
a small per-sample running accumulator — `(covered count, G/C count, depth sum)`
for the currently-open 500 bp window — and:

- **depth** per position is `ref_obs + nonref_obs` (both light);
- **GC** per position is the reference base, which the fold already fetches (the
  REF span), so no new data;
- when the stream crosses a window boundary, the window's `(gc, mean depth)` is
  final; it is attached to the variant loci in that window and spilled.

This is **memory-flat** (one open window per sample, `O(samples)` running state,
freed as the stream advances — no per-window map, no genome-wide structure),
**single-pass** (no re-reading the `.psp` bodies), needs **no schema change**
(the obs counts are already light), and is **consistent with the fitted model by
construction** — it reproduces the exact `Σ num_obs per covered position, tiled
by 500 bp, GC from the ref base` that the producer's coverage histogram (which
the model was fit on) was built from. The one wrinkle: a window's mean is not
final until its last position streams by, so a variant record buffers briefly
until its window closes — a tiny per-window buffer (windows are 500 bp, variants
sparse), local to the fold. This is the path.

**Fallback — a dedicated second pass over the `.psp` bodies, one sample at a
time** (what the R1 harness did). Re-open each sample's `.psp` after the main
pass, tile the depth into windows, keep only variant windows. Also memory-flat
(`O(windows for one sample)`, freed per sample), but it reads the bodies a second
time and duplicates the tiling logic. Keep it in reserve only if folding the
window accumulation into the hot caller path proves to entangle it badly; the
inline approach is preferred.

**Rejected — hold all windows for all samples in memory** (`O(samples ×
windows)` RAM, ≈ 13 GB at 1000 samples): breaks the memory-scaling thesis.

---

## 5. The spill file

**What it is and is not.** An **ephemeral scratch file**: written once, read
twice, deleted when the run ends (success or failure). It lives in scratch
space, never a durable output. It is *not* a `.psp` — the `.psp` is a first-class
per-sample artefact built once and reused across many callings; the spill exists
only to bound memory within a single run.

**Format.** The crate has **no parquet/arrow dependency**, so parquet is out
without adding one. Two in-house options: (a) reuse the `.psp` block-writer
pattern (`psp/writer.rs` — buffer records into a block, flush on a size
threshold, index the blocks), which the SSR work already generalised toward a
kind-agnostic container; or (b) a simple length-framed binary of per-locus
records. **Recommendation: reuse the block-writer/container infrastructure** —
it already gives us streaming writes, bounded buffering, and a read-back
iterator, and it is the codebase's one columnar abstraction. The exact framing
is an implementation choice for S2; the architecture only requires: streaming
append, one-pass read-back, and self-contained cleanup.

**What it holds per locus.** Everything the final VCF needs (position, alleles,
per-sample GT/GQ/AD, QUAL, existing FILTER, cohort AF) plus the paralog scoring
inputs (per-sample relative copy number, or the raw window `(gc, depth)` + a
model reference). Because it holds the whole VCF payload, its size is ~VCF-sized
— the resource we deliberately spend to keep RAM flat.

---

## 6. The decision criterion — answering "which is the criteria?"

Two different numbers, easy to conflate:

- **π answers "how many."** It is the estimated genome-wide *fraction* of loci
  that are hidden paralogs — one number, found by the EM loop of Q5 (repeatedly
  re-score every locus with the current π and re-average). On tomato2 it is ≈ 9%.
  π describes how common paralogs are; it does **not**, by itself, decide any
  single locus.

- **The FDR cut answers "which ones."** Each locus gets a posterior
  `sigmoid(LR + logit π)` and, from ranking all loci by that posterior, a
  **q-value** (the false-discovery rate at that cut). A locus is flagged when its
  `q_of_lr(LR) ≤ target FDR`. The target FDR is the **operating knob** — set low
  (default high-confidence ≈ 1%) so we flag only loci we are confident about and
  sacrifice very few real variants. This is the introgression-safety knob: a real
  introgressed variant has normal coverage, so its LR is negative and it never
  approaches the cut.

So the criterion for flagging a locus is **not** "π says paralogs exist" and
**not** "LR > 0"; it is **"this locus's q-value is at or below the CLI's target
FDR."** π is an input to that; the FDR target is the control.

**Estimating π does not fight the streaming design (the apparent chicken-and-egg).**
Counting paralog sites needs a threshold, which needs π — a circularity. The EM
loop resolves it without a hard count: each locus contributes a *soft* count (its
probability of being a paralog, `sigmoid(LR + logit π)`), and π is the average;
start from a guess, average, repeat. Crucially, **the EM iterates over the
fixed-size LR histogram, not over the data**: the streaming pass just drops each
LR into a bin (that needs no π at all), and the EM then loops over the ~2000 bins
in memory — a few dozen iterations, microseconds, no re-read of the spill or the
`.psp`. So the histogram is what decouples the single-pass streaming from the
iterative estimate; we get the exact, data-adaptive π at negligible cost.

Two alternatives to full EM were weighed and set aside as the *primary* path,
kept only as a **fallback**: a fixed "reasonable" π (rejected — π varies by
species/cohort/divergence/depth, so a constant mis-calibrates every posterior),
and a fixed LR cut-off (rejected — the LR scale is dataset-dependent and a raw
cut does *not* control the false-discovery rate, losing the introgression-safety
guarantee). The one place a fixed π earns its keep is when the EM fails to
converge: `ParalogPrior` carries a `converged` flag (Q5), so on a pathological
dataset we fall back to a documented default π and warn, rather than trusting a
junk estimate.

**`Hexp` scale (the R1 correction, load-bearing here).** `F = 1 − Hobs/Hexp`
compares two per-callable-position rates. `Hobs = n_het / callable_positions` is
already such a rate. `Hexp` must be too: the expected heterozygosity **per
callable position**, `Σ 2p(1−p) / callable`, *not* the mean of `2pq` over variant
sites. Getting this wrong (R1's first run did) saturates every sample's `F` at
0.99 and inflates π. S1 accumulates `Σ2pq` in the main pass and divides by a
representative callable count (e.g. the cohort median) once the pass is done.

---

## 7. The verdict — hard removal *(settled: owner chose removal)*

A locus whose q-value is at or below the target FDR is **dropped from the
output VCF entirely** — not written. This matches the crate's existing filter
convention exactly: the allele-balance and min-QUAL filters already drop records
in the `emit_or_drop` seam (`vcf_writer.rs`), and the paralog filter is one more
such drop, counted in the writer's stats (`records_dropped_paralog`, alongside
`records_dropped_allele_balance`). No new `FILTER` value is declared.

Consequences to keep in mind:

- **The cut must be conservative.** Because a dropped call is gone with no record
  of it, the introgression-safety of the FDR target does the load-bearing work: a
  real introgressed variant has normal coverage → negative LR → never near the
  cut. The default FDR (§8) is set low precisely so removal is safe.
- **Provenance in the header, not per-record.** Since flagged records are gone,
  the run's parameters — the target FDR, the estimated π, and the resolved LR
  cut — are recorded in the **VCF header** (S5) so a reader knows the filter ran
  and with what settings.
- **Optional: an INFO score on the survivors.** We *may* still stamp the paralog
  LR/posterior as an INFO field on the *kept* records, for downstream inspection,
  without changing which records are written. Left as an S4 implementation option;
  it does not affect the callset.

The mechanism is the `emit_or_drop` drop-path in `vcf_writer.rs`, mirroring the
allele-balance filter (`record_fails_allele_balance` → increment a dropped
counter → return without writing).

---

## 8. The CLI knob — on by default *(settled: owner chose on-by-default)*

Two flags on `CohortPipelineArgs` (`cli/shared_args.rs`), following the
allele-balance flags' pattern:

- `--paralog-fdr <q>` — the target false-discovery rate; the filter's operating
  point. **Defaults to ≈ 0.01** (introgression-safe), so the filter is **on by
  default**.
- `--no-paralog-filter` — turns the filter off, restoring the current
  single-pass, direct-to-VCF behaviour.

The chosen FDR, the estimated π, and the resolved LR cut are recorded in the VCF
header for provenance (S5).

**On-by-default is a real change to the pipeline's default path**, worth stating
plainly: every cohort `var-calling` run now (a) spills to disk and reads it back
twice instead of writing the VCF in one streaming pass, and (b) reads the `.psp`
bodies a second time for the coverage pass (§4). So:

- **T2 (cost/RSS) becomes a gate, not a footnote** — the default path's wall-time
  and peak-RSS delta must be acceptable, since every run pays it.
- **Existing var-calling integration tests will change output** (some variants are
  now dropped) and must be re-baselined as part of S4 — expected, not a
  regression, but it means the byte-identity tests need the filter turned *off*
  to keep asserting the streaming path, plus new tests for the filtered path.

---

## 9. Keeping memory flat (the T2 obligation)

The design is memory-flat in both variant count and sample count, by
construction:

- **Coverage models + `obs_het`:** `O(samples)` small structs, held for the run.
- **`Hexp`:** one running sum.
- **LR histogram:** fixed size (a few KB), independent of variant count.
- **Per-locus data:** on the **spill (disk)**, never a genome-wide RAM vector.
  Both spill passes stream from disk and hold **one locus at a time**; the write
  pass recomputes each LR rather than caching the calibrate pass's — so RAM does
  not grow with the variant count in either pass (see the invariant note in §3).
- **Per-window depth (§4 inline fold):** an O(1) running window accumulator per
  sample, freed as the coordinate-ordered stream advances — no per-window map.

Nothing in the added path grows with the number of loci in RAM, and nothing
grows with the number of samples beyond the `O(samples)` small per-sample state
the caller already keeps. T2 (the cost/RSS check) verifies this holds on a real
run.

---

## 10. Settled decisions + what remains for implementation

**Settled by the owner (2026-07-01):**

- **Verdict = hard removal** (§7): flagged loci are dropped from the VCF, like the
  existing filters; provenance goes in the header.
- **On by default** (§8): `--paralog-fdr` defaults to ≈ 1%; `--no-paralog-filter`
  disables. The spill + two-pass + coverage-pass path is now the default path.

**Recommendations settled here (revisit only if measurement forces it):**

- **Coverage sourcing = inline fold** (§4): aggregate the 500 bp window inline in
  the main-pass fold (the obs counts are already phase-1 light, so depth is free;
  GC from the already-fetched REF span). Single-pass, memory-flat, no schema
  change. The dedicated second body pass is kept only as a fallback.
- **Spill framing = reuse the `.psp`/container block-writer** (§5); exact framing
  finalised in S2.

**Left to confirm during implementation:**

- **LR-histogram bin resolution** — Q5 defaults to 2000 bins over `[−100, 100]`;
  confirm or tune against the tomato2 LR range in S3.
- **Default FDR value** — "≈ 1%" is the intent; the exact default is pinned in S5
  against the T1 flagged-set profile.

S1–S5 now implement directly against this structure: S1 the set-up + `Hexp`, S2
the main pass + spill + coverage pass + histogram, S3 the calibration, S4 the
write pass + removal, S5 the CLI. **The on-by-default choice makes T2 a gate on
S4**, and requires re-baselining the existing var-calling tests (with the filter
off for the byte-identity path).
