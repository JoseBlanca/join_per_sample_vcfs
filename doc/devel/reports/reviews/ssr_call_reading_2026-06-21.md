# Code Review: ssr-call reading & merge (Phases 0–3)

**Date:** 2026-06-21
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** SSR Stage-2 `ssr-call` cohort reader/merger — the new `src/ssr/cohort/` module + the enabling `OwnedRecordsIter` in `psp/reader.rs` + the `run_ssr_call` wiring (branch `ssr-cohort`, vs `main`).
**Status:** Request-changes

---

## 1. Scope

- **What:** a diff/PR — Phases 0–3 of the reading & merge layer.
- **Reviewed against:** branch `ssr-cohort` @ `2799005`.
- **In-scope files:**
  - [src/ssr/cohort/types.rs](../../../../src/ssr/cohort/types.rs)
  - [src/ssr/cohort/reader.rs](../../../../src/ssr/cohort/reader.rs)
  - [src/ssr/cohort/merge.rs](../../../../src/ssr/cohort/merge.rs)
  - [src/ssr/cohort/driver.rs](../../../../src/ssr/cohort/driver.rs)
  - [src/ssr/cohort/test_support.rs](../../../../src/ssr/cohort/test_support.rs)
  - [src/ssr/cohort/mod.rs](../../../../src/ssr/cohort/mod.rs)
  - [src/pop_var_caller/ssr_call.rs](../../../../src/pop_var_caller/ssr_call.rs)
  - `OwnedRecordsIter` + `into_records_of` in [src/psp/reader.rs](../../../../src/psp/reader.rs) + the export in [src/psp/mod.rs](../../../../src/psp/mod.rs)
- **Deliberately out of scope:** the genotyping EM/VCF (not built); the deferred threaded topology + prefetch pool (Phases 4–5); the rest of `psp/reader.rs` (`RecordsIter` and prior code); trivial CLI registration lines in `cli.rs`/`main.rs`/`mod.rs`.
- **Categories dispatched:** reliability, errors, naming, defaults, idiomatic, refactor_safety, module_structure, smells, extras. **Skipped:** `unsafe_concurrency` (no `unsafe`/threads/`Arc`/channels — single-threaded by design), `tooling` (no `Cargo.toml`/CI changes).

## 2. Verdict

**Request-changes.** 1 Blocker + 5 Major. The merge/cursor logic is correct for well-formed inputs and cleanly structured, but the cursor turns **data conditions into a `panic!` (B1) and a release-silent wrong result (M1)**, and decodes coordinates from disk without overflow protection (M2). These are all on the path that ingests arbitrary user `.ssr.psp`/catalog files. None is hard to fix.

## 3. Execution status

- `cargo fmt --check` → **exit 0**, clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → **exit 0**, `Finished` clean.
- `cargo doc --no-deps` → **exit 0**, generated clean.
- `cargo test --lib` → **`test result: ok. 1159 passed; 0 failed; 2 ignored`**.
- **Not run:** `cargo test --all-targets --all-features` (PROJECT_STATUS records a pre-existing, unrelated `psp_writer_perf` bench panic under `--all-targets`); `cargo audit` (no dependency changes). No findings labeled "Needs verification".

## 4. Open questions and assumptions

1. **Is "sample loci ⊆ catalog loci, and the catalog is coordinate-sorted" a guaranteed upstream precondition, or must `ssr-call` enforce it?** `ssr-pileup` produces a sample's loci by genotyping *against the catalog*, so in normal operation the subset + sort hold. But the same-catalog check only matches `reference_md5` (the reference, not the catalog *content*), and the catalog file is user-supplied (could be edited/unsorted). The answer sets the severity of **B1** and **M1**: if the caller must guarantee it, both drop a level and become "harden the boundary"; if `ssr-call` must tolerate arbitrary files, they stand.
2. **Are `.ssr.psp` inputs trusted (own-produced) or potentially corrupt/adversarial?** Determines whether **M2** (checked subtraction + typed error) is required or whether a `debug_assert` is acceptable. The crate is pre-alpha and own-produces these files today, but it is a cohort tool that opens many files by path.

## 5. Top 3 priorities

1. **B1** — `evidence_at` aborts the process (`panic!`) on a sample/catalog locus divergence that passes the md5 check; make it a typed `SsrCohortReadError`. [Findings → B1]
2. **M2** — `record.start - 1` / `record.end - 1` underflow to `u32::MAX` in release on a malformed `start`/`end == 0`; use `checked_sub` + a typed error. [Findings → M2]
3. **M4** — the duplicated `OwnedRecordsIter` state machine is only equivalence-tested on `SnpKind` + the happy path, never the `SsrKind` the cohort actually uses nor the error/poison path; extend the test (or factor a shared core) before the two copies drift. [Findings → M4]

## 6. Findings

### Blocker

**B1: [src/ssr/cohort/reader.rs:109](../../../../src/ssr/cohort/reader.rs#L109) — `evidence_at` panics on a data-reachable condition**
**Categories:** errors (primary); convergent: reliability, defaults, refactor_safety, smells, extras.
**Confidence:** High (that it panics). Severity contingent on Open Q1.

The `Ordering::Greater` arm `panic!`s when the merger asks for a locus *past* the cursor's held front — i.e. the cursor holds a stored locus the merger walked past without asking. The merger validates `catalog_reference_md5` and chromosome names, but **never checks that each sample record's `[start,end)` is a catalog locus**. A sample built against a catalog with the same `reference_md5` but different content (or any drift) makes this reachable from on-disk input, so a single divergent file aborts the whole cohort run with a panic instead of a typed error. The doc comment itself names the trigger "catalog/sample drift" — a data condition.

**Fix:** return a typed error instead of panicking. Add `SsrCohortReadError::UnexpectedLocus { held, asked }` and:
```rust
Some(Ordering::Greater) => Err(SsrCohortReadError::UnexpectedLocus {
    held: self.held.as_ref().map(|(f, _)| *f),
    asked: q,
}),
```
Keep a `debug_assert!`/comment documenting that the merger *should* prevent this, but do not trust on-disk data to uphold it. (See also M1 — the rewind sibling.) If Open Q1 resolves to "caller guarantees subset", downgrade to Major and harden as a `debug_assert` + one boundary test.

### Major

**M1: [src/ssr/cohort/reader.rs:93-99](../../../../src/ssr/cohort/reader.rs#L93-L99) — rewind guard is `debug_assert!`-only → silent wrong result in release**
**Categories:** reliability (primary); same root cause as B1.
**Confidence:** High.

The monotonic precondition is enforced only by `debug_assert!`. In a release build, a non-ascending `q` (e.g. an **unsorted user catalog**, or `LocusId` order diverging from catalog walk order) falls through to the `q < held` arm and returns `Ok(None)` — silently dropping the sample from a locus it actually covers, with no error. The only test (`rewind_panics`) pins debug behavior; release behavior is untested.

**Fix:** make the order violation a real error in all builds (a typed `SsrCohortReadError::NonAscendingQuery { prev, asked }` returned from `evidence_at`), or have the merger validate catalog monotonicity up front and document that as the guarantee. Add a test that exercises the release-surviving path. Resolve together with B1 — both are "the cursor trusts an ordering invariant that on-disk data can violate."

**M2: [src/ssr/cohort/reader.rs:140-144](../../../../src/ssr/cohort/reader.rs#L140-L144) — coordinate subtraction underflows on malformed input**
**Categories:** errors, extras, reliability, refactor_safety (convergent).
**Confidence:** High.

`adapt` computes `record.start - 1` and `record.end - 1` on values decoded from disk. The only guard is `debug_assert!(record.start >= 1)` — compiled out in release, and `end` has *no* guard at all. A malformed `start == 0` (or `end == 0`) wraps to `u32::MAX`, producing a bogus `LocusId` that then feeds the merge's ascending comparison → a silent wrong result or a downstream panic (B1), never a clean error. The `>= 1` invariant lives in the Stage-1 *writer*, not here.

**Fix:**
```rust
let start = record.start.checked_sub(1).ok_or(SsrCohortReadError::BadCoordinates { .. })?;
let end   = record.end.checked_sub(1).ok_or(SsrCohortReadError::BadCoordinates { .. })?;
```
plus a test feeding a `start == 0` / `end == 0` record. (Severity contingent on Open Q2 for the "adversarial" framing; the *asymmetry* — `start` guarded, `end` not — is a defect regardless.)

**M3: [src/ssr/cohort/reader.rs:143-144](../../../../src/ssr/cohort/reader.rs#L143-L144) ↔ [src/ssr/pileup/driver.rs:210-211](../../../../src/ssr/pileup/driver.rs#L210-L211) — the `+1`/`-1` coupling is prose-only**
**Categories:** refactor_safety.
**Confidence:** High.

The Stage-1 writer adds `+1` to both bounds; the cohort reader subtracts `-1`. They are linked only by parallel comments — no shared constant, no shared helper pair, no write→read round-trip test. Worse, the cohort tests *hand-write both sides* of the conversion (fixtures emit `start+1`, assertions expect `start`), so they cannot catch a drift on either side. A future change to one shift silently produces off-by-one loci.

**Fix:** add a round-trip test that writes a `SsrLocusObs` through the real Stage-1 `to_container_record`, reads it back through the cursor, and asserts the original 0-based coordinates — so the two shifts are pinned to each other by a single test rather than by comment. (A shared `psp`-level coordinate-shift helper used by both stages would be even safer.)

**M4: [src/psp/reader.rs](../../../../src/psp/reader.rs) `OwnedRecordsIter` — duplicated state machine, under-tested for drift**
**Categories:** refactor_safety (primary), smells, idiomatic (convergent).
**Confidence:** High.

`OwnedRecordsIter::{new, load_next_block, next}` is a near-verbatim copy of `RecordsIter`'s (the doc comments admit it). The deliberate choice to leave the production SNP path untouched is reasonable, but the only equivalence test (`owned_records_iter_matches_borrowing_across_blocks`) exercises **`SnpKind` on the clean happy path only** — never the `SsrKind` decoder the cohort actually drives it with, and never the shared error/poison branch. The two copies must now be hand-synchronized with no test that would catch divergence on the path that matters.

**Fix (pick one):** (a) extend the equivalence test to `SsrKind` *and* a mid-block decode-error case; or (b) factor the shared framing (`load_next_block` + the `next` loop) into one free function `fn step<R, S>(reader: &mut PspReader<R>, state: &mut StepState<S>, clamp: &RangeClamp)` that both iterators call, so the logic lives once. At minimum add `// KEEP IN SYNC WITH RecordsIter` markers at the two edit sites, not just in the type doc.

**M5: [src/pop_var_caller/ssr_call.rs:36-42](../../../../src/pop_var_caller/ssr_call.rs#L36-L42) / [src/ssr/cohort/driver.rs:46-53](../../../../src/ssr/cohort/driver.rs#L46-L53) — `--threads`/`--queue-depth` are silent no-ops with misleading help**
**Categories:** defaults (primary); convergent: reliability, smells.
**Confidence:** High.

`--threads` (default 4) and `--queue-depth` (default 4) parse, flow into `SsrCallConfig`, and are then never read — `driver::run` is single-threaded. A user passing `--threads 32` gets identical behavior with zero runtime feedback. The `--threads` help text ("Worker threads for the EM pool … output is identical … frozen-parameter genotyping") describes a phase that does not exist yet. `--help` is a primary default-visibility surface and currently misdescribes the running behavior.

**Fix:** emit a one-time `tracing::warn!` (or stderr note) when a non-default value is supplied for a reserved flag, and rewrite the help to say the flag is reserved/not-yet-effective. Optionally hide the flags until the EM lands.

### Minor

- **Mi1: [src/ssr/cohort/types.rs:85](../../../../src/ssr/cohort/types.rs#L85) — `CohortLocus.ref_tract` is misnamed.** Its own doc and the fixtures (`GGGGGGCACACATTTTTT` = 6 flank + 6 tract + 6 flank) confirm it holds tract **plus flanks**, while "tract" elsewhere on the struct means the repeat region only. Rename `ref_frame` / `ref_with_flanks`. (naming)
- **Mi2: [src/ssr/cohort/merge.rs:82](../../../../src/ssr/cohort/merge.rs#L82) — `SsrMergeError::Read(#[from] SsrCohortReadError)` collapses two origins** (cursor construction at line 163 vs. streaming at line 199) into one pathless variant labelled "streaming sample evidence" — unlike the sibling `Open`/`OpenPsp` which carry `path`. In an N-file cohort this loses *which* sample failed. Split the variants, or carry the sample label. (errors)
- **Mi3: [src/pop_var_caller/ssr_call.rs:5-8,31,35](../../../../src/pop_var_caller/ssr_call.rs#L5-L8) — docs contradict behavior.** The module banner still calls `run_ssr_call` a "stub" and `--output` help says "Output multi-sample VCF path", but the function now runs the merge and writes a **TSV**. The deferral is documented well in `driver.rs`, but the user-facing `--help` + banner disagree. Update them. (extras)
- **Mi4: [src/ssr/mod.rs:4](../../../../src/ssr/mod.rs#L4) — module-wide `#![allow(dead_code)]` masks the whole `ssr` tree** (catalog, pileup, cohort), not just the not-yet-wired cohort surface. Several `SsrQc` fields are written by the cursor and read by no non-test code, which this hides. Scope the allow to the specific forward-declared items, or narrow it to `cohort/`. (module_structure)
- **Mi5: `catalog_reference_md5` key is triplicated** — [merge.rs:31](../../../../src/ssr/cohort/merge.rs#L31) (production const), [test_support.rs:19](../../../../src/ssr/cohort/test_support.rs#L19) (duplicate const), [pileup/driver.rs:302](../../../../src/ssr/pileup/driver.rs#L302) (writer literal). The `test_support` copy is a silent-drift hazard: the merge tests pin to their own copy and could stop validating the real contract. Lift to one authority both stages import (e.g. `psp::registry_ssr`). (module_structure)
- **Mi6: [src/ssr/cohort/driver.rs:127-130](../../../../src/ssr/cohort/driver.rs#L127-L130) — `..SsrQc::default()` in a test literal** hides a future QC field from the compiler. Exhaustively construct it so a new field forces the test author to address it. (refactor_safety)
- **Mi7: missing `// PANIC-FREE:` invariant comments on `expect`** at [reader.rs:105](../../../../src/ssr/cohort/reader.rs#L105) and the `OwnedRecordsIter::new` decompressor `expect` in `psp/reader.rs` — both are genuinely sound, but the project convention wants the invariant spelled out. (errors)
- **Mi8: [src/ssr/cohort/merge.rs:138-148](../../../../src/ssr/cohort/merge.rs#L138-L148) — the catalog-md5 `_ =>` catch-all over the external `ParameterValue` enum** folds `None` and a non-`String` present value into "missing", with no `// REVIEW ON UPGRADE:` marker; a new `ParameterValue` variant would be absorbed silently. Add the marker (or match `None` explicitly and let a wrong-type value be its own error). (idiomatic)

### Nits

- [driver.rs `format_locus`](../../../../src/ssr/cohort/driver.rs#L91-L97) labels the per-sample field `alleles=`, but it reports the **distinct-observation count**, which the `SampleEvidence` doc explicitly says are "not alleles yet". Relabel `distinct=`/`nseq=`. (naming/extras)
- [driver.rs:77](../../../../src/ssr/cohort/driver.rs#L77) — `from_utf8_lossy` on the motif silently masks non-ASCII; motifs are ACGT, so a lossy substitution would hide corruption. (errors)
- [driver.rs:109](../../../../src/ssr/cohort/driver.rs#L109) — the `(u32, u32, usize)` test-helper tuple (`idx`, `depth`, `n_alleles`) is a transposition trap; a 3-field struct documents the positions. (smells)
- [merge.rs](../../../../src/ssr/cohort/merge.rs) — `path.display().to_string()` repeated within `open`; hoist once per file. (smells)
- [merge.rs:215-221](../../../../src/ssr/cohort/merge.rs#L215-L221) — `chrom_names()` relies on `name_to_global`'s values being a dense `0..n` permutation (true today via `enumerate`); state that invariant in a comment so a future change to the id assignment is on notice. (idiomatic — verified safe today)
- [reader.rs:143](../../../../src/ssr/cohort/reader.rs#L143) — the inline comment says "inclusive" for the `end` conversion, which is exclusive; fix the wording (the arithmetic is correct). (extras)

## 7. Out of scope observations

- None blocking. The Stage-1 writer's `to_container_record` (`pileup/driver.rs`) is the other half of M3's coupling and would be the natural home for a shared coordinate-shift helper if that fix is taken.

## 8. Missing tests to add now

Grouped by function under test.

- **`SampleEvidenceCursor::evidence_at`**
  - `evidence_at_returns_error_not_panic_when_asked_past_held` — a cursor holding locus L, asked for L' > L → typed error (after B1). Catches the process-abort.
  - `evidence_at_rejects_non_ascending_query_in_release` — two calls with descending `q` → typed error, not silent `None` (after M1). Must not be `#[cfg(debug_assertions)]`-only.
- **`SampleEvidenceCursor::adapt` / `new`**
  - `adapt_errors_on_zero_start_coordinate` and `adapt_errors_on_zero_end_coordinate` — malformed record → typed error, no `u32::MAX` wrap (M2).
  - `new_propagates_chrom_id_out_of_range` — a record whose `chrom_id` exceeds the remap → `ChromIdOutOfRange` (the variant is currently unconstructed in tests).
- **`CohortMerger::next_locus`**
  - `next_locus_propagates_a_corrupt_catalog_error` and `next_locus_propagates_a_cursor_decode_error` — error paths at merge.rs:178-199 are untested.
- **`OwnedRecordsIter`** (M4)
  - `owned_records_iter_matches_borrowing_for_ssr_kind` — equivalence on the `SsrKind` decoder, not just `SnpKind`.
  - `owned_records_iter_poisons_after_a_decode_error` and `owned_records_iter_empty_file_yields_none` — the error/poison and empty-file branches.
- **Coordinate round-trip** (M3)
  - `stage1_to_cohort_coordinate_round_trip` — write via the real `to_container_record`, read via the cursor, assert original 0-based coords.

## 9. What's good

- The `held` + `last_query` design (reader.rs) cleanly separates a *rewind* from a *sparse-Absent* with one extra field — the right model for the contract (the gap is enforcing it in release, M1, not the design).
- Module dependency direction is clean: `src/ssr/cohort/` imports only peer modules (`psp`, `ssr::types`, `ssr::catalog`) with no back-reference into any pipeline-stage module (verified by grep in module_structure).
- Adding `OwnedRecordsIter` instead of refactoring the production `RecordsIter` kept the byte-identity-critical SNP path untouched (psp/reader.rs) — good blast-radius discipline (it just needs better equivalence tests, M4).
- `write_dump` is generic over `W: Write` and `format_locus` is pure, so the output layer is unit-testable without touching the filesystem (driver.rs).
- `from_parts` (in-memory validation) is split from `open` (file I/O) and from the merge loop, which is exactly why the validation has 7 focused unit tests (merge.rs).

## 10. Commands to re-verify

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo doc --no-deps`
- `cargo test --lib`
- New (introduced by this review's recommendations): the §8 tests, then re-run `cargo test --lib`.

### Author response convention
Address each finding by id (B1, M1…Mi8) with `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the §4 open questions first — B1, M1, and M2 severities depend on them.
