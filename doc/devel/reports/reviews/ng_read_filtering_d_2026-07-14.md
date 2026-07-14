# Code Review + Fixes: ng read filtering — Milestone D (the ReadFilter iterator)
**Date:** 2026-07-14
**Reviewer:** rust-code-review skill (orchestrator, 8 categories)
**Scope:** `ReadFilter` iterator + `ReadFilterError` + `record_drop` + `header()` + D tests, in `src/ng/read/filtering.rs`; re-exports in `src/ng/read/mod.rs`
**Status:** Approve-with-changes → fixes applied

---

## Verdict
Approve-with-changes. **1 Blocker, 1 Major, several Minor/Nits.** Reviewed at HEAD `535f715`. 8 categories dispatched (reliability, errors, naming, idiomatic, refactor_safety, module_structure, smells, extras); the state machine was traced correct (no correctness defect in the pipeline itself).

## Findings and resolution

**Blocker (errors) — `finish` is bypassable → silent fatal-error loss.** `ReadFilter: Iterator`, so `for read in filter {}` moves+drops the filter without calling `finish`, and a latched fatal error (truncated file / corrupt record / out-of-bounds fetch) then looks like a clean EOF. This touches the spec §5 `Item = MappedRead` decision (which explicitly accepted "nothing forces full consumption"), so rather than change the surface it is **strongly mitigated**: `finish` now `take`s the error (marking it observed) + `mem::take`s the counts; a `Drop` guard `debug_assert!`s that a latched error was observed (skipped while `panicking()`), catching the misuse in every dev/test run; `#[must_use]` on `finish`; and a prominent doc obligation on `ReadFilter`. **Residual** (release-build silent loss if `finish` is never called) is inherent to the spec's design and surfaced to the owner — the pipeline consumer must call `finish`.

**Major (errors) — `expect` without `// PANIC-FREE:`** on `u32::try_from(i)` in `new`. **Applied** — added the marker + rationale (`ref_id` is a 32-bit field; the index fits by construction).

**Minor (refactor_safety) — pre-decode used `if let Drop(..)`, not an exhaustive `match`,** so a future third `FilterVerdict` variant would silently fall through to decode/keep. **Applied** — made it an exhaustive `match` mirroring the post-decode arm.

**Minor (extras) — the fixture didn't exercise a read failing BOTH #9 and #8** (the only case where ng attribution diverges from production, and the stated reason the hand-count is "against the ng order"). **Applied** — the fixture's bad-CIGAR read is now also high-mismatch, so the exact counts (`bad_cigar==1`, `high_mismatch==1`, not `0`/`2`) now *discriminate* the #9-before-#8 order end-to-end.

**Minor (naming) — `fatal` → `fatal_error`, `buf` → `record_buf`** (bare adjective / unqualified, inconsistent with `ref_buf`). **Applied.**

**Minor (errors) — `Reference(#[from])`** was asymmetric with `Source`/`Decode`'s `#[source]` and its generated `From` was unused. **Applied** — changed to `#[source]`.

**Minor (module_structure) — the `record_source.rs` split is now warranted** (~264 lines of noodles adapters + traits in a 1770-line file). **Deferred** — a pure-move refactor, tracked as a follow-up commit.

**Minor (module_structure) — `FilterVerdict`/`DropReason` over-exposed** in the re-exports (cascade-internal, no public consumer). **Applied** — dropped from `mod.rs`'s `pub use`.

**Minor (smells) — stale module header** ("C and D follow"); **`one_contig_200_header` duplicated `one_contig_header`**; **latch-and-return repeated 3×** in `next()`. **Applied** — refreshed the header; extracted `contig_header(length)`; extracted a `fail()` helper.

**Reliability missing-tests — Applied:** empty-source (zero counts) + running-tally-before-exhaustion. (The remaining suggested classes — 0-`@SQ` header, decode-stays-stopped-with-trailing-record — deferred as low-value.)

**Deferred / not changed:** `SetupError` vs `RefSeqError` two-shapes (the arch §6 open item — revisit when a second setup failure mode appears); the `record_source.rs` split (follow-up); a couple of Nits (Source variant name, `R: RawRefSeq` bound on `counts`/`finish`).

## Validation
`cargo fmt -- --check` (ng) clean; `cargo clippy --lib` clean; `cargo test --lib -- ng::read::filtering` → **35 pass**.
