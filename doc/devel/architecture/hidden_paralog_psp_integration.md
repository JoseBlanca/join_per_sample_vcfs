# Storing the hidden-paralog per-sample summaries in the `.psp`

**Status:** draft, 2026-06-29. Architecture working doc for the
`tomato2-paralog-filter` branch. Decides *how* the two per-sample
summary statistics the paralog filter needs (spec
[`hidden_paralog_filter.md`](../specs/hidden_paralog_filter.md) §4, §5,
§8) get produced and carried so var-calling consumes them without a
separate pre-pass. Grown one agreed premise at a time; the model
*intent* lives in the spec, the *shape* lives here.

The two statistics:

- **Per-sample observed heterozygosity** — the individual's overall het
  rate (tracked directly, not the derived F = 1 − Hobs/Hexp).
- **Per-sample coverage-by-GC** — the depth∼GC relationship + the
  single-copy coverage scale, used to GC-correct the windowed coverage
  that is the introgression-safe core of the filter (§4).

---

## Grounding facts (what the code already does)

These four facts decide the premises; established by reading the pileup
walker, the Stage-1 pipeline, and the `.psp` writer/header/trailer.

1. **The `.psp` is a full per-position pileup, not a variant-sites
   file.** The walker advances one reference position at a time and
   folds a Match for every covered read at every covered position;
   [`open_record::process_position`](../../../src/pileup/walker/open_record.rs)
   opens a record wherever any contributor has an event — i.e. every
   covered base. (Roundtrip test: two pure-REF reads over `ACGTA` emit
   **5 records**,
   [`pileup_to_psp.rs:143`](../../../src/pileup/per_sample/pileup_to_psp.rs#L143).)
   REF is always `alleles[0]`; per-position fragment depth = Σ
   `allele-obs-count` over the record's alleles — post-filter, dedup,
   mate-overlap-resolved, exactly the depth §4 wants. So **windowed
   depth is derivable from the body**, but only by scanning the whole
   body plus the reference for GC.

2. **The per-sample walk is serial on one thread.** Stage 1 is producer
   → BAQ worker pool → **a single consumer thread that runs the walker
   and writes records in coordinate order**
   ([`stage1_pipeline.rs:351`](../../../src/pop_var_caller/stage1_pipeline.rs#L351)).
   Parallelism is *across* samples (separate processes) and inside BAQ;
   the record stream itself is single-threaded. An accumulator hung off
   the record stream is free of synchronization and per-sample-parallel
   for free.

3. **The header is committed before the body.** `PspWriter::new` writes
   the framed head TOML immediately (fixed `u64` length prefix); blocks
   stream; the block index + fixed 32-byte trailer are written at
   `finish()`. A **body-derived** per-sample summary therefore cannot
   live in the head TOML without a second pass or seek-back — it has to
   be a trailing section, the same shape as the block index. (Detail
   for Premise 4, but load-bearing, so noted here.)

4. **Observed het is a single-sample observable.** It is the fraction of
   *this individual's* sites that are heterozygous (spec §5/§7 track Hobs
   directly, not the cohort-dependent F = 1 − Hobs/Hexp). The pileup
   calls no joint genotypes, but a *rough* per-site genotype from the
   allele fractions is enough for this coarse prior — so het, too, is
   estimable in the Stage-1 walk without the cohort EM.

---

## Premise 1 — where each statistic is computed *(SETTLED)*

**Both statistics are computed in the Stage-1 pileup and stored in the
`.psp`.** Var-calling does **no** pre-pass — it reads both per-sample
summaries from each `.psp` it already opens. This is the symmetric spec
§8 "fold the coverage model into Stage 1" vision, extended to het.

### 1a. Coverage-by-GC

A pure function of the per-sample record stream + reference GC; the walk
already visits every position single-threaded with the `ref_fetcher` in
hand (facts 1, 2); per-sample-parallel for free.

### 1b. Observed heterozygosity — a single-sample observable, estimated rough in the walk

**Decisive point:** the spec (§5, §7) tracks the **observed** per-sample
het rate *directly*, **not** the derived F = 1 − Hobs/Hexp. Only Hexp
needs cohort allele frequencies; **Hobs is inherently a single-sample
quantity** — the fraction of that individual's sites that are
heterozygous. So no cohort genotypes are needed, and an earlier draft
that routed het through a Stage-2 cohort pre-pass was over-engineered.

Hobs is used as a **prior** that places a sample on the selfer↔outbred
axis (tomato2 F ranged 0.00–0.97, §7). Separating F≈0 from F≈0.9 is
coarse, so a rough estimate suffices. In the same walk we call a rough
per-site genotype for each variant site from its `(k_alt, n_total)`
fragment counts — **not** a depth-blind VAF threshold but a **binomial
het-vs-hom likelihood ratio** (`logLR = n·ln½ − k·ln(1−ε) − (n−k)·ln ε`;
the binomial coefficient cancels, so it is a couple of multiplies). A
confidence margin `M` splits variant sites three ways: **confident het**
(`logLR > +M`), **confident hom-alt** (`logLR < −M`), and **ambiguous**
(`|logLR| ≤ M`). This is depth-aware by construction — at high depth the
LR is sharp and almost nothing is ambiguous; at ~6× many sites land in
ambiguous, exactly as they should. `Hobs = n_het / (n_het + n_hom_alt)`
(the confident ratio); `n_ambiguous` is kept as the **uncertainty
signal** (large in low-coverage samples → trust their Hobs less).

**Residual bias, bounded:** collapsed paralogs inflate the het count
(their capped-below-1 VAF often reads ~0.5), biasing Hobs up / F down.
Tolerable because (1) it's a coarse prior, and (2) the coverage signal
is built in the same walk, so the het tally can optionally be **gated to
exclude excess-depth positions** (rough running single-copy depth as the
cutoff) to de-bias it for free. Build the plain estimator first; add the
gate only if the bias shows up.

### The coverage-by-GC computation, concretely (the "where / how / threading" answer)

- **Computation site:** the consumer/writer side of Stage 1 — the
  `drive_*_to_psp` loop that already touches every record in coordinate
  order and holds the `ref_fetcher`. One place, already single-threaded,
  no new synchronization.
- **GC is a reference property, not a read property.** The GC covariate
  for a window is computed from the FASTA — independent of the reads.
  The reads only supply depth. Window scale is fixed at the analysis
  window (Premise 3, 500 bp default), so no separate fragment-scale
  machinery is needed.
- **Accumulation grain = tiled, non-overlapping windows, over covered
  positions only.** The walker emits one record per *covered* position
  in coordinate order. We keep running per-tile accumulators —
  `covered_count`, `gc_count` (REF base ∈ {G,C}), `depth_sum` (Σ
  `allele-obs-count` over the record's alleles = total fragment depth at
  that position). At each tile boundary we emit one sample
  `(gc_count/covered_count, depth_sum/covered_count)` — GC fraction and
  **covered-bases mean depth**, both over exactly the covered positions
  — then reset. Uncovered positions never emit a record, so they
  self-exclude; a position whose REF base is `N` is excluded from both
  `gc_count` and `covered_count` (GC undefined). Both quantities are
  measured over the *same* positions the scoring stage sees in the body
  (REF base is `alleles[0].seq`, depth is Σ obs), so training and
  application units coincide exactly.
- **Accumulator = a 2D count histogram** `[GC bin][depth bin]`, one cell
  bumped per tile. No steady-state allocation, no locks.
- **Lifetime:** the accumulator must **persist across regions** and be
  finalized **once**. With `--regions` / multi-region input, several
  walkers feed one shared `PspWriter`
  ([`drive_region_into_writer`](../../../src/pileup/per_sample/pileup_to_psp.rs#L84));
  the accumulator lives next to the writer (not per-region) and is
  reduced at `finish()`.
- **At finish:** the raw histogram is serialized into the metadata
  section (Premise 2 stores it raw; the curve / single-copy-scale fit is
  downstream, in var-calling).
- **Het rides the same walk** (Premise 1b): a running tally of
  `(n_het, n_hom_alt, n_ambiguous)` over variant sites from the per-site
  binomial-LR genotype calls, also persisted across regions and finalized
  once. It reduces to **four counts** (the three + `n_variant` = their
  sum), not a histogram (see Premise 2's coverage-vs-het asymmetry).
- **Cost:** a few-KB array + a handful of counters, O(1) per record, on a
  thread that is already the serial bottleneck's *consumer*, not the
  bottleneck. Expected to be in the noise; to be measured, not assumed.

Open sub-questions deferred: depth-axis bounding for the histogram; the
exact rough-genotype thresholds and whether the het tally is
coverage-gated (Premise 1b).

---

## Premise 2 — what exactly each summary stores *(SETTLED)*

The two summaries store at **different fidelities, on purpose**:

- **Coverage-by-GC → raw 2D histogram.** Its consumer is a *model still
  being calibrated* (binned-median curve, single-copy-scale estimator,
  maybe quantiles), so keeping raw sufficient statistics lets var-calling
  re-fit without re-running the pileup. Mirrors the `.psp` carrying
  `mapq-sum` / `mapq-sum-sq` raw and computing Welch's-t downstream.
- **Observed het → four counts** over variant sites: `n_het_sites`,
  `n_hom_alt_sites`, `n_ambiguous_sites` (and `n_variant_sites` = their
  sum). `Hobs = n_het / (n_het + n_hom_alt)` and its support (and the
  ambiguous fraction = low-coverage uncertainty) are all recoverable.
  There is **no downstream fit** to preserve raw stats for — the
  binomial-LR three-way classification at a *fixed* margin `M` is a
  settled computation — so a logLR histogram would buy nothing the four
  counts don't. This is the asymmetry: raw distributions are worth it only
  when a calibrating model consumes them (coverage), not when the
  statistic is settled (het). The price of counts-not-histogram: `M` is
  frozen at pileup time (recorded), so retuning the het/hom boundary means
  re-running the pileup — acceptable for a settled threshold.

### Coverage histogram details

**Store the raw 2D histogram** `[GC bin][depth bin]` of tiles — *not* a
pre-reduced curve. The model fit (binned-median depth∼GC curve +
single-copy scale) happens in the var-calling consumer.

What a tile contributes (Premise 1b): one `(GC fraction, covered-bases
mean depth)` sample, both over the tile's covered positions. **Covered-
bases mean** (no covered-fraction floor): mean depth over only the
positions that have a record, matching what the scoring stage computes
from the body.

Stored fields (the trailing section, Premise 4):
- the `[GC bin][depth bin]` count matrix (`u32` cells);
- **bin schemes**, so the consumer can interpret cells: GC bin count /
  edges (default: uniform 1% bins, 0..1), depth bin width + cap +
  overflow bin (defaults TBD by tuning — fine resolution near the ~6×
  single-copy mode, a catch-all tail bin for high-copy paralogs);
- window width (= 500, Premise 3) and tile-anchor convention;
- support counts: total tiles, tiles excluded for all-`N`, etc.

Bin-resolution defaults are tuning, not architecture — they are
recorded in the section so they can change without breaking readers.
The curve and single-copy scale are **derived downstream**, not stored.

### Het counts

`n_het_sites`, `n_hom_alt_sites`, `n_ambiguous_sites`, `n_variant_sites`
(four integers, the last = sum of the first three) + the rough-genotype
parameters recorded for reproducibility: `min_depth`, the error rate `ε`,
and the confidence margin `M`. `Hobs = n_het / (n_het + n_hom_alt)` is
formed downstream; `n_ambiguous` is the per-sample uncertainty weight.

## Premise 3 — the GC window scale *(SETTLED)*

**One fixed window scale, CLI-configurable, default 500 bp.** It serves
as *both* the GC covariate window (Stage 1, building the curve) and the
analysis window at scoring time (Stage 2) — a single width so the
training and application units coincide. We do **not** pursue the spec
§9 fragment-scale refinement now; a fixed window is simpler and 500 bp
already sits at a good operating point (§4.3, single-copy mode robust-SD
≈ 0.23).

Decisions:
- **CLI:** a window-size flag on the `pileup` subcommand (e.g.
  `--gc-window-bp`), default 500.
- **Tiled, non-overlapping** windows anchored at genome coordinates
  (0, 500, 1000, …), so a scoring-time locus maps to exactly one tile.
- **The width is recorded in the `.psp`** alongside the stored curve, so
  the var-calling consumer computes per-window `gc_rel` over the *same*
  width without its own flag (it inherits the producer's choice).

Owed within this premise (carried into Premise 2's detail): how a tile
with N / reference-gap bases (GC undefined) or partial coverage is
handled — exclude, or down-weight.

## Premise 4 — where in the `.psp` format the summary lives *(SETTLED)*

**Both** per-sample summaries are `.psp` concerns now (Premise 1b): the
coverage-by-GC histogram **and** the het counts. Both are body-derived
(known only after the walk), so neither can go in the head (fact 3);
both land in a new metadata section. The core `.psp` writer stays
schema-agnostic — the SNP-specific accumulator lives in the **Stage-1
driver** (the `drive_*_to_psp` seam that already loops every record),
*not* in the core `PspWriter`; the driver serializes the summaries and
hands the bytes to the writer at `finish()` through a generic
"attach metadata section" API, opaque to the core, interpreted by the
`snp` kind.

**Decision — one general metadata section, located by arithmetic, trailer untouched:**

File layout becomes `head · blocks · index · metadata · trailer`. The
metadata section is:

- **General-purpose, not summary-specific.** Future-proofs the §8
  pre-pass family (SSR params, contamination) — each is a new key in the
  doc, not a new format change. The `snp` kind populates it with the
  coverage histogram + het counts; the core defines only "an optional
  compressed TOML metadata section."
- **Written forward-only at `finish()`** with real values in hand — no
  space reservation, no seek-back-and-patch, so the writer stays
  `W: Write` (no `Seek` bound). Placed between the block index and the
  trailer.
- **Located by arithmetic, so the 32-byte trailer is byte-identical.**
  The reader reads the section as the bytes between index-end
  (`index_offset + index_byte_length`, both already in the trailer) and
  trailer-start. **A zero-length gap = no metadata**, so files/kinds
  that don't write one stay valid. The only format change is "the reader
  stops requiring index-end == trailer-start," gated by a format-version
  minor bump — no trailer struct change, no section table.
- **zstd-framed**, so it carries its own frame checksum: the index is
  covered by the trailer's `index_checksum`, the metadata by its frame
  checksum, the trailer by its magic — no unchecked region.
- **TOML**, matching the head and the contamination artefact
  ([`contamination_artefact.rs`](../../../src/pop_var_caller/contamination_artefact.rs#L184)),
  one serialization stack. The coverage histogram (Premise 2) is a flat
  row-major integer array + the bin-scheme / window scalars; the het
  summary is the four counts + recorded parameters (`min_depth`, `ε`,
  `M`) — both plain TOML.

**Why not header-reserve-and-patch** (considered): the size *is* fixed
by config so reservation is feasible, but it forces `W: Write + Seek` +
a two-phase patch in `finish()` (losing the forward-only writer), and a
reserved region can't be zstd-compressed. The tail section keeps the
writer forward-only and compresses the sparse histogram.

**Wrinkle:** TOML is config-shaped, not ideal for a dense numeric
matrix; compression makes the size moot at modest bin resolution
(~3 K cells), and the doc can carry a binary sub-blob later if
resolution is ever cranked hard. Not needed now.
