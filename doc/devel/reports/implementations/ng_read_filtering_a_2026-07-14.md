# Implementation report: ng read filtering — Milestone A (types + scaffold)

**Date:** 2026-07-14
**Feature:** ng step 1 — read filtering
**Plan:** [read_filtering.md](../../ng/impl_plan/read_filtering.md) (Milestone A, steps A1–A3)
**Spec / arch:** [spec](../../ng/spec/read_filtering.md), [arch](../../ng/arch/read_filtering.md)

## 1. Plan

Land Milestone A of the read-filtering plan: the type foundations, no filtering
logic. Steps A1 (scaffold `read/`), A2 (seed `types.rs` scalar newtypes), A3
(step-1-local config/verdict/counts types). Because all three are pure-types /
no-logic under one checkpoint and tightly coupled, they were run as a single
plan-driven loop iteration (implement → review → apply → commit), named
explicitly in the commit — the sanctioned bundling for near-empty/type-only
steps.

## 2. Assumptions

- **Step bundling** A1+A2+A3 into one commit (above) — recorded here and in the
  commit message rather than committed silently.
- **`get()` on unconstrained newtypes** kept despite the duplication with the
  `pub` field, because `ng_step_interfaces.md` §1 mandates every newtype expose
  `.get()` for a uniform accessor surface.
- **Flat `ReadFilterConfig`** (`max_read_mismatch_fraction` + `mismatch_bq_floor`
  as sibling fields, not a nested `Option<MismatchFilter{…}>`) — deliberate
  mirror of the production `AlignmentMergedReaderConfig` per spec §4 / arch §2.2.
- **`DomainError` is `#[non_exhaustive]`** and holds one variant now
  (`MismatchFraction`); later constrained types add their own — matches the
  `RefSeqError` convention already in `ref_seq.rs`.

## 3. Changes made

- **`src/ng/read/mod.rs`** (new) — declares `pub mod filtering`; documents the
  deliberate `read/`-folder deviation (fixed prelude + future `ReadPrep`
  siblings).
- **`src/ng/read/filtering.rs`** (new) — step-1-local types only: `ReadFilterConfig`
  (+ `Default` = the production policy, thresholds from the reused `DEFAULT_*`
  constants), `FilterVerdict`, `DropReason` (9 variants, 1:1 with the counts),
  `ReadFilterCounts` (running tally). No logic.
- **`src/ng/types.rs`** — seeded the scalar newtypes read filtering touches:
  `MapQual`, `BaseQual`, `Bp` (unconstrained, `pub` field + `.get()`),
  `MismatchFraction` (constrained `[0,1]`, private field, checked `try_new`), and
  the ng-wide `DomainError` (first variant).
- **`src/ng/mod.rs`** — wired `pub mod read;`, re-exported the new vocabulary,
  refreshed the header doc.

## 4. Tests added

- `mismatch_fraction_accepts_boundary_values` — `0.0`, `1.0`, `0.10` accepted.
- `mismatch_fraction_rejects_out_of_range` — `-0.01`, `1.01`, `±INFINITY` rejected
  with `DomainError::MismatchFraction`.
- `mismatch_fraction_rejects_nan` — NaN rejected.
- `unconstrained_newtypes_expose_their_value` — `MapQual`/`BaseQual`/`Bp` `.get()`.
- `default_config_reproduces_the_production_filter_policy` — every default field
  cross-checked against the reused `DEFAULT_*` constants (the port anchor for A).
- `counts_default_is_all_zero` — explicit all-zero literal pins every counter.

## 5. Validation

Run in the dev container (`./scripts/dev.sh`):
- `cargo fmt -- --check` → in-scope ng files clean.
- `cargo clippy --lib` → no warnings on the in-scope files.
- `cargo test --lib -- ng::read::filtering ng::types` → 6 in-scope tests pass.

## 6. Tradeoffs and follow-ups

- **Deferred to Milestone B:** compile-time enforcement of the `DropReason` ↔
  `ReadFilterCounts` 1:1 mapping (the tally/`match` site lands with the cascade;
  A3 is deliberately no-logic). Acceptance criterion for B.
- The cascade (B), record-source seam (C), and `ReadFilter` iterator (D) follow.
