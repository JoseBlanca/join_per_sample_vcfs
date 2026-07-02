# Per-sample paralog summaries in the `.psp` — implementation plan

**Status:** draft, 2026-06-29, branch `tomato2-paralog-filter`. Turns the
settled architecture
([hidden_paralog_psp_integration.md](../architecture/hidden_paralog_psp_integration.md),
all four premises SETTLED) into build order. Design intent lives in the
spec [hidden_paralog_filter.md](../specs/hidden_paralog_filter.md) (§4,
§5, §8); the empirical prototype lives in
[`benchmarks/tomato2/`](../../../benchmarks/tomato2/).

## Domain intent

The hidden-paralog filter needs two per-sample summaries that are not
derivable from the VCF: a **coverage-by-GC** model (the introgression-
safe core, §4) and the **observed heterozygosity** prior (§5). Both are
single-sample quantities computable in the Stage-1 pileup walk, so we
produce them there and carry them in the `.psp` — var-calling then
consumes them with **no pre-pass**. This plan delivers the *production*
and *storage* of the two summaries plus the reader API to consume them.

The walk already visits every covered position single-threaded on the
Stage-1 consumer thread, with the `ref_fetcher` in hand, so both
accumulators ride that stream at O(1)/record with no new synchronization.

## Scope / non-goals

In scope: the `.psp` metadata-section format (container core), the SNP
summary payload (coverage histogram + het counts), the two accumulators,
the CLI knob, the Stage-1 wiring, the reader-side parse, and empirical
parity against the Python prototype.

**Out of scope** (separate follow-on plan): the model that *consumes*
the summaries — the depth∼GC curve / single-copy-scale fit, the H1-vs-H2
marginal LR (§5), the empirical-Bayes posterior and FDR cut (§6). This
plan stops at "the data is produced, stored, and readable."

## Principles (how the order was chosen)

- **Types first, then implementation**, within every step (project rule).
- **Container core before kind payload before wiring.** The schema-
  agnostic metadata-section mechanism (Milestone A) is independent of
  what SNP puts in it (B); both precede the Stage-1 plumbing (C).
- **Additive at the tail.** The metadata section sits between the block
  index and the trailer; block/index/trailer bytes are unperturbed and
  the 32-byte trailer stays byte-identical. A `.psp` with no section
  (zero-length gap) stays valid — the SSR path writes none.
- **Determinism + byte-identity.** Existing record/body bytes and the
  cross-thread/region byte-identity invariants must not change; only the
  new tail section appears. Pause between milestones (incremental rule).
- **Validate against the prototype.** The Rust accumulators must
  reproduce the tomato2 Python prototype's histogram + het before we
  build anything on top.

---

## Milestone A — `.psp` metadata section (schema-agnostic container core)

**A1. Metadata-section wire format + framing.** ☐
Types first: a `MetadataSection` value (opaque compressed bytes + its
decoded length) and the offset arithmetic. Then writer + reader framing.
- *Writer:* a generic `PspWriter::attach_metadata(bytes)` (or a
  `finish_with_metadata(bytes)` variant) that, at `finish()`, zstd-frames
  the bytes and writes them **between** the block index and the trailer.
  The 32-byte trailer encoding is unchanged.
- *Reader:* locate the section as the bytes between index-end
  (`index_offset + index_byte_length`) and trailer-start; a zero-length
  gap ⇒ `None`. Decompress; verify the zstd frame checksum. Drop the
  reader's current "index-end == trailer-start" assertion (relax to
  "≤ trailer-start").
- *Format version:* minor bump; document that old readers reject via the
  existing format-version gate.
- **Deliverable:** `src/psp/` framing + round-trip unit tests: write a
  section → read it back byte-for-byte; **empty case** (no section) reads
  `None` and the file is byte-identical to today's; truncated/corrupt
  section frame errors cleanly; SSR-kind file (no section) still reads.
- *Depends:* none. *Source:* arch Premise 4; `trailer.rs`, `index.rs`,
  `writer.rs::finish`.

---

## Milestone B — SNP summary payload + accumulators

**B1. SNP metadata payload types + TOML serde.** ☐
Types first: a `SnpSampleSummary` document — the coverage histogram
(flat row-major `u32` matrix + GC/depth bin schemes + window bp +
support counts) and the het block (`n_het_sites`, `n_variant_sites` +
the rough-genotype thresholds). Serialize/deserialize as TOML (one
serialization stack with the head + contamination artefact).
- **Deliverable:** the payload structs + `to_toml`/`from_toml` + a
  round-trip test (build a populated summary, serialize, parse, assert
  equality). No accumulation logic yet.
- *Depends:* A1 (it is the bytes A1 carries). *Source:* arch Premise 2.

**B2. Coverage accumulator.** ☐
Types first: a `CoverageByGcAccumulator` with its config (window bp, GC
bin scheme, depth bin scheme) and the per-tile running state
(`covered_count`, `gc_count`, `depth_sum`). Then the update/finalize.
- Per covered position: `gc_count += REF∈{G,C}`, `depth_sum += Σ
  allele-obs-count`, `covered_count += 1` (skip `N` REF). At each tile
  boundary: bump `histogram[gc_bin(gc_count/covered_count)][depth_bin(
  depth_sum/covered_count)]`, reset. Finalize → the `SnpSampleSummary`
  coverage fields.
- **Deliverable:** the accumulator type + unit tests over synthetic
  records (uniform-depth tiles land in the expected cell; `N` excluded;
  tile boundaries respected; multi-region feed accumulates into one
  histogram).
- *Depends:* B1 (output type). *Source:* arch Premise 1b/2/3.

**B3. Het accumulator.** ☐
Types first: a `HetAccumulator` with thresholds (min depth, minor-allele
VAF band defining "het", min minor-allele support defining "variant
site") and the two counters. Then the per-record rough genotype.
- Per record: compute minor-allele fraction from the allele obs counts;
  if it clears the variant-site threshold, `n_variant += 1` and
  `n_het += (VAF ∈ het band)`. Finalize → the het fields.
- **Deliverable:** the accumulator type + unit tests (clear het / clear
  hom-alt / sub-threshold sites counted correctly; multi-region feed).
- *Depends:* B1. *Source:* arch Premise 1b/2.

---

## Milestone C — wire into Stage 1

**C1. CLI flag + config plumbing.** ☐
Add `--gc-window-bp` (default 500) to the `pileup` subcommand; thread it
into the accumulator config; record it (and the het thresholds) in the
writer's `[writer.parameters]` for provenance.
- **Deliverable:** the flag parsed/validated, defaulted, surfaced in the
  `.psp` header parameters; CLI help text; a test that a non-default
  value reaches the accumulator and is recorded.
- *Depends:* B2/B3. *Source:* arch Premise 3; `cli.rs`,
  [pop_var_caller_pileup_cli.md](pop_var_caller_pileup_cli.md).

**C2. Hook the accumulators into the pileup→psp seam.** ☐
The accumulators live next to the shared `PspWriter` in the
`drive_*_to_psp` seam (not per-region), are fed each record, persist
across regions, and at end-of-run serialize the `SnpSampleSummary`
(B1) and attach it via A1's API before `finish()`.
- **Deliverable:** end-to-end — run the `pileup` subcommand on a small
  fixture, read the `.psp` back, assert the section is present and its
  histogram/het match an independent in-test accumulation; assert the
  **body + index + trailer bytes are unchanged** vs a build with the
  section disabled (additive-tail invariant).
- *Depends:* A1, B2, B3, C1. *Source:* arch Premise 1b/4;
  [`pileup_to_psp.rs`](../../../src/pileup/per_sample/pileup_to_psp.rs),
  [`stage1_pipeline.rs`](../../../src/pop_var_caller/stage1_pipeline.rs).

---

## Milestone D — reader API + validation

**D1. Reader-side summary accessor.** ☐
Expose the parsed `SnpSampleSummary` off the opened `.psp` (parse A1's
section bytes via B1) so a downstream consumer gets the histogram + het
counts without re-deriving them. Read-only; no model fit here.
- **Deliverable:** the accessor + a test reading a produced `.psp`;
  `None` when the file carries no section.
- *Depends:* A1, B1, C2. *Source:* arch Premise 4.

**D2. Empirical parity with the Python prototype.** ☐
Run the Rust producer on the tomato2 cohort and confirm the in-Rust
coverage histogram and het reproduce the validated prototype
([`build_gc_normalization.py`](../../../benchmarks/tomato2/src/build_gc_normalization.py)
and the het derivation) within tolerance.
- **Deliverable:** a parity check (script or test) + a short note in
  `ia/reports/implementations/`; any discrepancy resolved or documented.
- *Depends:* C2, D1. *Source:* spec §4, §10; `benchmarks/tomato2/`.

**D3. Cost check.** ☐
Measure the per-sample pileup wall/RSS delta from the two accumulators
on a representative run; confirm it is in the noise (arch expectation),
record the number.
- **Deliverable:** before/after numbers in the D2 report; a fast-path
  fix only if the delta is non-trivial.
- *Depends:* C2. *Source:* arch Premise 1b ("to be measured, not
  assumed"); [[project_pileup_thread_scaling]].

---

## Checkpoints

- After **A1**: the container can carry an opaque section, round-trips,
  empty-case byte-identical. Pause.
- After **C2**: a real `pileup` run emits both summaries in the `.psp`,
  additive-tail invariant proven. Pause — first end-to-end artefact.
- After **D2**: the produced numbers match the validated prototype — the
  data is trustworthy for the follow-on filter-model plan.

## Open items carried from the architecture doc

- Depth-axis bounding / bin-resolution defaults for the coverage
  histogram (tuning, recorded in the section — B1).
- The exact rough-genotype thresholds, and whether the het tally is
  coverage-gated to neutralize paralog inflation (B3; arch Premise 1b).
- TOML-for-dense-matrix wrinkle: revisit a binary sub-blob only if bin
  resolution is cranked hard (arch Premise 4).
