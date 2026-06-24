# SSR Stage-1 — `locus_record` aggregation (all-CSR)

**Date:** 2026-06-15 · **Skill:** rust-feature-implementation · **Branch:** `ssr-architecture`

Fold a locus's per-read `ReadOutcome`s into its evidence record (arch §8.2/§11).
Built after a design discussion that **revised the storage model** (recorded below).

## 1. Design decisions settled (with the user)

- **All-CSR, no histogram, no weight column.** Every spanning read is stored as
  one pruned `Qᵣ` profile; the spec §4.3 confident→histogram split is dropped.
  Rationale: per-locus reads are depth-capped (no in-memory problem), zstd
  compresses the redundant columns on disk, and the cohort step can summarize
  later — a measure-first call. The `hist_weight` base-qual aggregate is
  unnecessary now because, under realign-everything, **every read's forward
  log-lik already encodes its base quality** (and quality shows in the profile's
  peakedness). Diverges from spec §4.3 (histogram + weight) — amend like the
  realign-everything pivot.
- **Renormalized per-read log-probs**, not raw log-liks. Provably lossless for
  Stage 2: scaling a read's `Qᵣ` by `Z_r = Σ_L Qᵣ(L)` is a constant independent
  of allele `a` and stutter `θ` (`Qᵣ` is stutter-free), so it cancels in every
  genotype ratio/argmax. Quality is preserved in the profile shape.
- **Pruning kept** (`AMB_LL_DROP`): drop a read's candidates more than
  `AMB_LL_DROP` nats below its best, then renormalize the survivors — turns the
  dense `2W+1` profile into a sparse 1–3-length one. Placeholder value (≈4 nats),
  a calibration knob (arch §14).
- **QC scalars:** derive `n_spanning`/`n_flanking`/`n_frr` from the outcomes;
  take `depth`/`n_filtered`/`mapped_reads` from the fetcher (a `QcCounts` input).
  **Drop `n_flank_indel`** — a vestigial fast-path-gate-era count realign-everything
  no longer produces.
- **In-memory record**, profiles as `Vec`s; CSR flattening is the (deferred)
  container's job. Off-ladder evidence deferred (empty).

## 2. Changes made

[src/ssr/pileup/locus_record.rs](../../../src/ssr/pileup/locus_record.rs):
- `SsrLocusRecord` — QC scalars + `spanning: Vec<Vec<(u16, f32)>>` (one
  renormalized profile per spanning read). No `n_flank_indel`, no histogram, no
  weight.
- `QcCounts { depth, n_filtered, mapped_reads }` — the fetcher-supplied counts.
- `aggregate(locus, &[ReadOutcome], QcCounts) -> SsrLocusRecord`.
- `prune_and_renormalize` — `AMB_LL_DROP` pruning + log-sum-exp renormalization;
  `debug_assert`s off-ladder candidates don't appear yet.
- `AMB_LL_DROP = 4.0` (calibration placeholder).

## 3. Tests added (6)

The three classification counts derived from a mixed outcome list; locus coords +
`QcCounts` copied through; a sharp read prunes to a single length at log-prob 0;
a bimodal read keeps both lengths, normalized, ordered; pruning drops the far
tail; per-read profiles are normalized (mass sums to 1).

## 4. Validation results

Dev container (`./scripts/dev.sh`):
- `cargo fmt -- --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ssr::pileup::locus_record` — **6 passed, 0 failed**.

## 5. Tradeoffs and follow-ups

- Storage model diverges from spec §4.3 (all-CSR vs histogram+weight); amend the
  spec/arch (a doc follow-up, like the realign-everything pivot).
- `SsrLocusRecord` is in-memory; the container schema (deferred) will flatten the
  profiles to CSR columns. `SsrLocusRecord` placement may move to a top-level peer
  when the container consumer lands.
- Off-ladder evidence + the `count_repeats` fast-path shortcut remain deferred.
- **Next:** `fetch_reads` (the I/O fetcher: index walk, reservoir, bundles +
  `QcCounts`), then the driver — or the container refactor (the writer).
