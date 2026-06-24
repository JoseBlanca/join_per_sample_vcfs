# SSR Stage 1 — `ssr-pileup` Mark-2 implementation report

**Date:** 2026-06-17. **Branch:** `ssr-pileup-mark2`. **Skill:**
rust-feature-implementation.

Built the **Mark-2 empirical-candidate** Stage-1 pipeline, replacing the Mark-1
reference-anchored rung model. Specs/plan:
[architecture/ssr_ladder_model.md](../../../doc/devel/architecture/ssr_ladder_model.md),
[architecture/ssr_pileup_mark2.md](../../../doc/devel/architecture/ssr_pileup_mark2.md),
[implementation_plans/ssr_pileup_mark2.md](../../../doc/devel/implementation_plans/ssr_pileup_mark2.md).

## Plan

Mark-2 makes candidate alleles **observed sequences** (not reference rungs): per
read, align it against the locus reference frame to **delimit** the repeat region,
**quality-gate** it, and tally **observed sequence → count** per locus. The
reference is only a coordinate frame; there is no on/off-ladder and no per-read
likelihood in Stage 1 (deferred to Stage 2). Built in 7 steps, each compiling +
tested, reusing the Mark-1 read-I/O path almost verbatim and swapping the per-read
compute + the on-disk schema.

## Assumptions / decisions made

- **`footprint.rs` kept as its own module** (not folded into `fetch_reads.rs` as
  the doc first said) — folding collided the two files' test helpers; the
  descriptive name addresses the original objection to the opaque `reach.rs`.
- **`passes_quality_gate` takes a scratch buffer** (`ViterbiScratch.qual_buf`) for
  the quantile select, rather than allocating per read.
- **No write-side obs-seq length cap.** obs-seq is a spanning-read tract (≪ the
  10 000-byte `MAX_ALLELE_SEQ_LEN`), so it is structurally unreachable; the decoder
  still caps it defensively for foreign files.
- **`BorderOffEnd` is detected positionally** (`tract_start == 0` or
  `tract_end == m`), relying on Stage 0's clean unique flanks + the reach gate
  rather than a flank-match-quality threshold (a v1 simplification).
- **Atomic cutover (steps 5–7).** `registry_ssr`'s only consumer was `ssr_mark1`,
  so the schema rewrite, the new driver, the CLI cutover, and the `ssr_mark1`
  deletion had to land in one commit (no compiling intermediate). User approved.

## Changes made

New `src/ssr/` (Mark-2), copied-then-adapted from `ssr_mark1` where noted:

- `types.rs` — `Locus`/`Motif` (copied); `Allele`/`NormalizedSeq` dropped.
- `catalog/` — Stage 0, copied verbatim (model-agnostic).
- `pileup/fetch_reads.rs` — reservoir depth cap + per-locus fetch (copied).
- `pileup/footprint.rs` — read-vs-locus geometry + the `reaches_locus` admission
  gate (lifted from Mark-1 `triage.rs`; rung machinery dropped).
- `pileup/alignment.rs` — **new:** the per-Q Viterbi (max-path) pair-HMM with a
  full backpointer matrix + traceback (`delimit_read` → `Region | BorderOffEnd`),
  tie-break `Match > Deletion > Insertion`, junction indel → preceding (5′) block;
  the first-quartile quality gate (`MIN_REGION_Q1 = 15`).
- `pileup/locus_tally.rs` — **new:** `ReadObs` / `SsrLocusObs` / `tally` (observed
  sequences → counts, sorted by bytes; lean QC scalar set).
- `pileup/driver.rs` — the Mark-2 driver (fetch → delimit + gate → tally → write),
  reusing the Mark-1 batched-`par_chunks` skeleton; `to_container_record`,
  `qc_counts`, header build.

Container + cutover:

- `psp/registry_ssr.rs` — **rewritten** to the Mark-2 schema (per-locus QC scalars
  + per-observation `obs-count`/`obs-seq-len`/`obs-seq` keyed by sequence). Rides
  the existing `Bytes`/varint codecs; no `schema_version` bump (pre-alpha).
- `psp/writer.rs` — SSR `validate_locus` drops the NaN-loglik sweep.
- `pop_var_caller/ssr_catalog.rs` + `ssr_pileup.rs` — repointed to `crate::ssr`;
  the `--window` knob removed.
- `lib.rs` — `pub mod ssr;` added, `ssr_mark1` removed.
- **Deleted:** `src/ssr_mark1/`, `benches/ssr_pileup_perf.rs`,
  `examples/profile_ssr_pileup.rs` + its `Cargo.toml` bench entry.

## Tests added/updated

- `alignment` (9): clean / longer / shorter / impure tract extracted verbatim
  between the flanks; missing left/right flank → `BorderOffEnd`; determinism; the
  gate keys on the first quartile not the min.
- `locus_tally` (4): sequence tally + byte-sort + dropout counts; order-independence.
- `registry_ssr` (4): multi-block round-trip (incl. a zero-observation locus, a
  many-sequence locus, a chromosome boundary); contig-end round-trip; wrong-kind
  poison; column well-formedness.
- `driver` (5): `qc_counts`; the name→id +1-shift adapter; end-to-end
  catalog→reference→BAM→`.ssr.psp` (clean read → one observed `CACACA`, count 1);
  **thread-count byte-identity**.
- `fetch_reads`/`footprint`/`types`/`catalog`: copied tests pass under the new path.

Anti-tautology: delimiter/tally tests build reads independently of the code under
test.

## Validation results

Run in the dev container:
- `cargo fmt --check` — clean.
- `cargo clippy --all-targets -- -D warnings` — clean.
- `cargo test` — **1116 lib tests + integration + doctests, 0 failed** (2 ignored).

## Tradeoffs and follow-ups

- **Calibration placeholders:** `MIN_REGION_Q1 = 15` (P1) and `MAX_READS_PER_LOCUS`
  need tuning on real data; measure the long-allele dropout the quality gate causes.
- **Unbanded DP (P6):** bounded per-locus by the reservoir cap; band later only if
  profiling shows it binds (the content pre-probe lever is gone with `ssr_mark1` —
  re-derive if needed).
- **No Mark-2 bench yet** — the Mark-1 `ssr_pileup_perf` bench was deleted; a Mark-2
  bench is a follow-up.
- **Stage 2 (`ssr-call`)** — candidate assembly (S1), stutter reachability between
  sequences (S2), the flat-error likelihood HMM (S3), EM, VCF — designed
  just-in-time after this.
- **Spec amendment owed:** spec §4.2/§4.3/§5.1 still describe the Mark-1 model.
