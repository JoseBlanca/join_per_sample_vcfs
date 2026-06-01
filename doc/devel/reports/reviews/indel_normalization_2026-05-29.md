# Code Review: indel_normalization
**Date:** 2026-05-29
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Stage 1 indel left-alignment (normalization) feature on branch `indel-normalization`
**Status:** Approve-with-changes — **SUPERSEDED by a restructure (same day)**

> **Note (post-review).** This review covered the first implementation,
> which added indel normalization in the BAQ prep stage (`BaqEngine::process`).
> The review (and its sub-agents, scoped to the diff) **missed** that the
> BAM/CRAM input cascade already left-aligned every read via an always-on
> "F3" pass (`src/bam/alignment_input.rs`, pre-existing on `main`). The
> feature was subsequently restructured to **replace F3 with the GATK port**
> and revert the prep-stage detour. That restructure **mooted M2/M3/M4/M5**
> (artifacts of the reverted detour), **resolved M7** (canonical-form tests
> added), and **applied M1 + Mi1**. The findings below are retained as a
> record of the first implementation; a fresh review of the F3-replacement
> diff is the right next step.

---

## 1. Scope

- **What:** PR diff — the indel left-alignment feature (commits `6a06ab2` → `932aa02`) on branch `indel-normalization`.
- **Reviewed against:** `main` (merge-base `8d7dd17`); HEAD `932aa02`.
- **In-scope files:**
  - [src/pileup/walker/indel_norm.rs](../../../../src/pileup/walker/indel_norm.rs) (new — the port)
  - [src/pileup/walker/indel_norm/tests.rs](../../../../src/pileup/walker/indel_norm/tests.rs) (new)
  - [src/pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs) (module decl)
  - [src/pileup/per_sample/baq_engine.rs](../../../../src/pileup/per_sample/baq_engine.rs) (`process` + `apply_baq` + `raw_ref` + normalization)
  - [src/pileup/per_sample/baq_stream.rs](../../../../src/pileup/per_sample/baq_stream.rs) (`apply_baq` field/threading)
  - [src/pileup/per_sample/baq_tests.rs](../../../../src/pileup/per_sample/baq_tests.rs) (new tests + call-site updates)
  - [src/pop_var_caller/stage1_pipeline.rs](../../../../src/pop_var_caller/stage1_pipeline.rs) (unified BaqStream branch)
  - [src/pop_var_caller/var_calling_from_bam.rs](../../../../src/pop_var_caller/var_calling_from_bam.rs) (unified BaqStream branch)
- **Deliberately out of scope:**
  - [src/var_calling/posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs) — formatting-only reflow swept into commit `db4acc9` by `git add -A` + `cargo fmt`; no behavioral change (see Out of scope observations).
  - benches/baq_perf.rs, examples/dhat_baq.rs — trivial one-arg call-site updates for the new `apply_baq`/`BaqStream::new` signature.
  - doc/devel/implementation_plans/indel_normalization.md — the plan (spec, not code).
- **Categories dispatched:** reliability, errors, naming, defaults, idiomatic, refactor_safety, module_structure, unsafe_concurrency, smells, extras. (Tooling skipped — no `Cargo.toml`/build-config/feature changes in the diff.)

## 2. Verdict

**Approve-with-changes.** The algorithm is a faithful, well-documented port (GATK `leftAlignIndels`/`normalizeAlleles`, cross-checked vs freebayes), it is correct on well-formed input, the concurrency story is clean, the diff matches the stated intent with no scope creep, and the full suite is green (996 tests). The changes to make before merge cluster around **release-build robustness on untrusted CRAM/BAM input** (M1), an **observability regression introduced by the `--no-baq` unification** (M3), and **missing tests for the canonical-form logic the feature's value depends on** (M7). The remaining Majors are a refactor-safety hazard (M4), naming drift (M5), an error-classification gap (M2), and a module placement issue (M6).

## 3. Execution status

- `cargo fmt --check` → **exit 0**, clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → **clean**, no warnings.
- `cargo test --all-targets --all-features` → **exit 0**, **996 passed / 0 failed** across 8 suites.
- `cargo doc --no-deps` → **fails** with 13 `unresolved-intra-doc-link` errors. All 13 are **pre-existing project-wide** broken links (the crate denies `broken_intra_doc_links`); **none reference the in-scope new files** (verified by grepping the doc output for `indel_norm`/`baq_engine`/`baq_stream`). Noted, not attributed to this PR.
- `cargo audit` → **not run**: `cargo-audit` is not installed in this environment.
- Findings labeled "Needs verification": **2** (M1 reachability of malformed CIGAR-vs-seq input through the noodles decoder; M5 the exact `FilterCounts.baq_rejected` field name).

## 4. Open questions and assumptions

1. **Does the noodles CRAM/BAM decoder guarantee that a record's CIGAR read-consumption equals `seq.len()`?** If yes, M1's corrupt-CIGAR path is unreachable from production input and M1 drops to a defensive-hardening Minor. If no (the decoders in `src/bam/` add no such cross-check), M1 stands. Affects **M1**.
2. **Is the malformed-CIGAR corrupt output actually caught by the walker's `read.length()` check?** Left-alignment preserves read-consumption, so a malformed input (read-consume ≠ `seq.len`) yields an output that is also ≠ `seq.len`, which `WalkerState::admit_read` rejects as `MalformedRead` (a hard error, not silent corruption). The residual risk is a malformed input whose *rewritten* read-consumption coincidentally equals `seq.len` and is semantically wrong — not demonstrated. This is why M1 is filed Major-with-verification rather than Blocker. Affects **M1**.
3. **Is `--no-baq` expected to ever drop reads?** Pre-PR it never did (passthrough kept all reads). Post-PR, mandatory normalization means an indel read with an unfetchable window is skipped under `--no-baq`. Confirm this is acceptable (it is the price of mandatory normalization). Affects **M2, M3**.
4. **What is the canonical name for the rolled-up skip counter (`FilterCounts.baq_rejected`)?** Determines the blast radius of the M5 rename. Affects **M5**.

## 5. Top 3 priorities

1. **M1** — Validate the read-side footprint in `left_align_cigar` and make the read-accounting invariant release-safe; the only guards are `debug_assert`s compiled out of the shipping binary, on a path fed by untrusted CRAM/BAM. [indel_norm.rs:385,482]
2. **M3** — Stop reporting `baq_skip_counts = None` under `--no-baq`: mandatory prep now drops reads even with BAQ off, and those drops are silently discarded. [stage1_pipeline.rs:165]
3. **M7** — Add tests for `build_cigar`'s deletion-before-insertion reorder and indel-collision merge — the canonical-form logic that makes identical indels consolidate is entirely untested. [indel_norm.rs:262]

## 6. Findings

### Major

#### M1: src/pileup/walker/indel_norm.rs:385 — Read-side footprint never validated; release builds can emit a corrupt CIGAR on malformed input
**Categories:** extras, errors, reliability
**Confidence:** Medium (Needs verification — see Open questions 1–2)

`left_align_cigar`'s early-return guard validates only the reference side (`read_start + ref_length > ref_bases.len()`). It then seeds `read_indel_range` from `read_seq.len()` and shifts it by the CIGAR's *claimed* read lengths. When a CIGAR's read-consumption disagrees with `seq.len()` (no cross-check exists in `process` — it only checks `qual.len() == seq.len()` at baq_engine.rs:132 — nor in the decoders), `read_indel_range.start` ends non-zero. The `debug_assert_eq!(read_indel_range.start, 0, ...)` at [indel_norm.rs:482](../../../../src/pileup/walker/indel_norm.rs#L482) catches this in debug (the extras sub-agent reproduced `left: -1`), but in release the assert is gone and the routine emits a CIGAR whose op lengths no longer match the read (and `Range::size()`'s `(end-start) as usize as u32` truncates silently if a range went negative). This is untrusted-input territory.

**Backstop (why Major, not Blocker):** the walker's `admit_read` runs `read.length()` (cigar-read-consume vs `seq.len`) and hard-errors on mismatch, so the demonstrated case aborts the run rather than corrupting calls. The silent-corruption window requires a malformed input whose rewritten length coincidentally matches `seq.len` — unverified.

**Fix:** add the read-side clause to the existing fail-safe so normalization is a no-op on malformed reads (the walker then rejects them exactly as it does today):
```rust
let read_consumed: usize = cigar.iter().copied().map(length_on_read).sum();
if read_start + ref_length > ref_bases.len()
    || read_consumed != read_seq.len()
    || read_seq.is_empty()
{
    return LeftAlignResult { cigar: cigar.to_vec(), leading_deletion_bases_removed: 0 };
}
```
Add a unit test (`left_align_cigar_returns_input_on_malformed_read`). Consider promoting the `:482` invariant to a release-safe fail-safe (return input unchanged) keeping `debug_assert!(false, …)` so tests still fail loudly. **Verify Open questions 1–2 first** — if the decoder guarantees the invariant, downgrade to Minor (defensive hardening).

#### M2: src/pileup/per_sample/baq_engine.rs:194-203 — Genuine reference-fetch I/O failure is masked as a past-chrom-end clamp and the read is silently dropped
**Categories:** errors (+ reliability cross-cat)
**Confidence:** High. **Pre-existing** behavior (the arm existed before this PR), but the `--no-baq` unification **extends its blast radius**: indel reads under `--no-baq` now fetch a window where pre-PR they did not, so this drop path is newly reachable in no-BAQ mode.

`ref_fetcher.fetch(...)` returns a typed `ChromRefFetchError` with four variants (`OutOfBounds`, `InvalidStart`, `OutOfPattern`, `Io { source }`). The `Err(e)` arm collapses **all** into `BaqSkipReason::RefWindowPastChromEnd`, prints a `warning: … skipping read`, and drops the read. A deleted/truncated FASTA, unreadable `.fai`, or disk error is misattributed to a normal end-of-contig clamp and counted in `ref_window_past_chrom_end`. On a failing disk this sheds unbounded reads (missing coverage → missing calls) with no distinguishable signal.

**Fix:** match the variant — only `OutOfBounds` maps to the clamp-skip; `Io`/`OutOfPattern`/`InvalidStart` should propagate as a hard error (preferred) or at minimum route to a new `BaqSkipReason::RefFetchFailed` with its own counter and an `error:`-level message. (Diff in `tmp/review_2026-05-29_indel-norm/errors.md`.)

#### M3: src/pop_var_caller/stage1_pipeline.rs:165 — `--no-baq` reports `baq_skip_counts = None`, hiding prep-stage read drops
**Categories:** defaults
**Confidence:** High. **Introduced by this PR.**

Pre-PR, `--no-baq` used a passthrough that never dropped reads, so `None` ("baq: disabled") was accurate. Post-PR, `BaqEngine::process` is the mandatory prep stage and several `BaqSkipReason`s fire under `--no-baq` — `Unmapped`/`EmptyQuery`/`QualAbsent`/`PosOutOfRange`/`ReadTooLong`/`ChromIdOutOfRange` unconditionally, plus `NoMatchInCigar`/`ContainsRefSkip`/`RefWindowPastChromEnd` for any indel-bearing read (gated by `needs_norm`). Those drops are counted in `BaqStream::skip_counts` but [stage1_pipeline.rs:165](../../../../src/pop_var_caller/stage1_pipeline.rs#L165) throws the snapshot away as `None`. A user debugging "why is this indel missing under `--no-baq`" gets no signal.

**Fix:** always surface the counts (`let skip = Some(*baq_stream.skip_counts());`) and let the summary printer render "baq: disabled" from the `no_baq` flag it already holds. The counter now measures prep-stage skips, not BAQ skips, so the `Option` no longer carries useful semantics. (This couples with M5's naming.)

#### M4: src/pileup/per_sample/baq_stream.rs:132 — `apply_baq` is an inverted-polarity positional bool re-derived (`!no_baq`) at each call site with no single source of truth
**Categories:** refactor_safety (+ naming cross-cat)
**Confidence:** High.

The mode is `no_baq: bool` at the driver layer but flips to `apply_baq: bool` at the engine layer; the `!no_baq` inversion is performed independently at [stage1_pipeline.rs:149](../../../../src/pop_var_caller/stage1_pipeline.rs#L149) and [var_calling_from_bam.rs:756](../../../../src/pop_var_caller/var_calling_from_bam.rs#L756), each passing it as the **6th positional `bool`** to `BaqStream::new` (next to `chunk_size: usize`). Nothing forces the two sites to agree; dropping a `!` at one site still compiles and silently runs the wrong BAQ mode for an entire run — which changes called genotypes — with no test failure.

**Fix:** a two-variant `BaqMode` enum with a single `from_no_baq(no_baq)` conversion, threaded through `BaqStream::new` and `BaqEngine::process`; the `if apply_baq` blocks become `if mode.applies_hmm()`. Makes a transposed/un-negated argument a type error. (Full snippet in `tmp/review_2026-05-29_indel-norm/refactor_safety.md`.)

#### M5: src/pileup/per_sample/baq_engine.rs:23,76,120 — `Baq*` type names no longer describe the now-mandatory prep stage
**Categories:** naming, defaults
**Confidence:** Medium.

The module docs now call this "the per-read prep stage" and state normalization is "mandatory and independent of BAQ", yet `BaqEngine`/`BaqStream`/`BaqOutcome::Capped`/`BaqSkipReason`/`BaqSkipCounts` all claim BAQ semantics. Under `--no-baq` a read is normalized but never capped, yet it returns `BaqOutcome::Capped`, and prep-failure drops (`Unmapped`, …) roll up into a `baq_rejected`-named counter. The mismatch will keep accreting as more prep steps land here.

**Fix:** rename the prep cluster to its domain role — `BaqEngine→ReadPrepEngine`, `BaqStream→ReadPrepStream`, `BaqOutcome→ReadPrepOutcome`/`Capped→Prepared`, `BaqSkip*→ReadPrepSkip*`, files `baq_{engine,stream}.rs→read_prep_{engine,stream}.rs` — keeping genuinely BAQ-only names (`apply_baq`, `apply_baq_cap_into`, `BaqConfig`, `bq_baq`). Wide mechanical rename; if deferred, record a justification at the type definitions and **verify the `FilterCounts.baq_rejected` field name** (Open question 4).

#### M6: src/pileup/walker/mod.rs:21 — `indel_norm` is placed under `walker/` but its only consumer is the sibling `per_sample/` subtree
**Categories:** module_structure
**Confidence:** High.

A consumer sweep shows the sole non-test consumer of `indel_norm` is `per_sample/baq_engine` (one symbol, `left_align_prepared`); `walker` itself never calls it. The walker's module doc frames it as not knowing about BAQ/prep, yet hosts a module only the prep stage uses — a cross-subtree back-reference, and `crate::pileup::walker::indel_norm` advertises a walker relationship that doesn't exist.

**Fix:** `git mv` `indel_norm.rs` (+ tests) to `src/pileup/per_sample/`, declare it there, and change its imports to `use crate::pileup::walker::{CigarOp, PreparedRead}` (the types correctly stay in `walker`). `baq_engine.rs:13` becomes `use super::indel_norm;`. (A neutral top-level peer is *not* recommended — only one consumer.) Land the test file on the `per_sample` `baq_tests.rs`-style sibling convention while moving (folds in Mi11).

#### M7: src/pileup/walker/indel_norm.rs:262 — The canonical-form logic in `build_cigar` (D-before-I reorder, adjacent-op merge, indel collision-merge) is untested
**Categories:** reliability
**Confidence:** High.

`build_cigar` step (2) reorders an indel run to put the deletion before the insertion and step (3) merges adjacent identical ops; left-alignment can also make two indels collide and merge ("collisions merge", per intent). Every existing test asserts a single isolated indel — no input ever produces an I-before-D run, two merging indels, or adjacent same-ops needing a merge. A reorder/merge bug yields a non-canonical CIGAR that fails to bucket identical events — *the exact bug this feature exists to remove* — with no panic. The `unreachable!("run is all indels")` and the `splice`/`i += rep_len` bookkeeping at [indel_norm.rs:286-298](../../../../src/pileup/walker/indel_norm.rs#L286) are unexercised against a mixed run.

**Fix:** add `build_cigar_orders_deletion_before_insertion`, `build_cigar_merges_adjacent_identical_ops`, and `left_align_cigar_merges_colliding_indels` (see §8).

### Minor

- **Mi1: indel_norm.rs:143 — `Range::size()` wraps silently on an inverted range.** `(self.end - self.start) as usize` returns ~`2^64` if `end < start`; this feeds `CigarOp::Deletion(… as u32)`. Add `debug_assert!(self.end >= self.start)` in `size()` and/or emit op lengths via `u32::try_from` failing safe. *(extras, idiomatic — convergent.)*
- **Mi2: indel_norm.rs:367-489 — Per-indel-read hot-path allocations with no regression bench.** `left_align_cigar`/`build_cigar` allocate 4–5 fresh `Vec`s per indel-bearing read (`result_right_to_left`, `nonzero`, `merged`, `replacement`) plus drop the moved-out `cigar`; the codebase otherwise prefers reusable scratch (`BaqEngine` already owns `raw_ref`/`encoded_ref`). Indel-free reads are correctly on the fast path. No criterion bench covers the normalization path. Recommend: thread `&mut Vec` scratch from `BaqEngine` (mirroring `apply_baq_cap_into`) **and** add a `benches/` criterion fixture with a regression threshold before optimizing. *(extras.)*
- **Mi3: baq_engine.rs:437 vs indel_norm.rs:75 — `is_indel_op` duplicates `is_indel` byte-for-byte.** Two names for one predicate across the feature's two modules. Collapse to one (call `indel_norm::is_indel`). *(naming, idiomatic, reliability — convergent.)*
- **Mi4: indel_norm.rs:272-303 — `build_cigar` step (2) is a manual `splice`/index loop with an `unreachable!`.** Rebuild via a grouping fold instead of in-place splice (the output is consumed wholesale by step (3)); removes the `i += rep_len` arithmetic and the `unreachable!`. *(idiomatic.)*
- **Mi5: indel_norm.rs:543-545 — `count_mismatches` silently truncates (early `return`) when blocks run past the buffers**, masking the very discrepancy the guard exists to detect. If the guard is kept, make the out-of-bounds case an explicit mismatch/flag, not a quiet return. *(errors, reliability.)*
- **Mi6: indel_norm.rs:385-390 — The "window too short" fail-safe returns input unchanged with no log/telemetry.** A real upstream window-sizing regression would silently degrade indel recall. Add a one-line `eprintln!`/debug counter and a unit test driving the branch. *(errors.)*
- **Mi7: baq_engine.rs:442 — `read.pos as u32` in `mapped_to_prepared` is an unchecked re-narrowing** of the field already range-checked via `i32::try_from(read.pos)` at line 139 — the one line that breaks the file's stated "plain `as` casts wrap silently" discipline. Thread the validated `pos_0` through instead, or `debug_assert!` the invariant. *(idiomatic.)*
- **Mi8: indel_norm.rs:367 — `remove_deletions_at_ends: bool` is an opaque trailing positional bool** on a coordinate-affecting routine. One production caller (always `false`). Enum-ify (`EndDeletions::{Strip,Keep}`) or at least bind it to a named local at the call site. *(refactor_safety.)*
- **Mi9: indel_norm.rs:136 — `struct Range` shadows `std::ops::Range`** and is a shape-name for a signed mutable index window. Rename to `IndexRange` (matches the GATK source and the line-129 comment); consider `len` over `size`. *(naming.)*
- **Mi10: indel_norm.rs:172-193 — The `*_is_same` predicate trio triplicates the `zip().all()` skeleton** with only the index expression varying (six `as usize` casts). Factor into one helper taking an index selector. *(idiomatic.)*
- **Mi11: indel_norm/ holds only `tests.rs`** while the impl is the flat sibling `indel_norm.rs` — the only place in the tree using a directory purely for a test file. Match the `per_sample` `baq_engine.rs`/`baq_tests.rs` convention (folds into M6's move). *(module_structure.)*
- **Mi12: baq_engine.rs:104 — "normalization always on" is documented but not runtime-inspectable.** No log/counter records that left-alignment ran or how many reads it rewrote. Emit a one-line summary in the existing run-summary path. *(defaults.)*
- **Mi13: indel_norm.rs:367,29 — `left_align_cigar` and `LeftAlignResult` are `pub(crate)` but have no consumers outside the module** (only `left_align_prepared` is called externally; the other two are test-only). Tighten to private with `#[cfg(test)]` access. *(module_structure cross-cat.)*

### Nits
- indel_norm.rs:286 — `unreachable!("run is all indels")` lacks the project-required `// UNREACHABLE: …` justification comment.
- indel_norm.rs:510,529 — `count_mismatches` and its sole call site carry two coupled `#[cfg(debug_assertions)]` gates; fold into one debug-only helper so they can't drift.
- indel_norm.rs:536-538 — `sp`/`rp` cursor abbreviations; prefer `ref_pos`/`read_pos`.
- indel_norm/tests.rs:8-9 — `use super::*` plus `use crate::pileup::walker::CigarOp::{self,*}` double-imports `CigarOp`.
- indel_norm.rs:172,180,188 — `*_is_same` predicate names read awkwardly ("same as what?"); `all_seqs_agree_at_*` states the predicate more directly.

## 7. Out of scope observations

- **src/var_calling/posterior_engine.rs** — A pure `cargo fmt` reflow (40 lines, no behavior change) was swept into commit `db4acc9` by `git add -A`. Harmless, but it pollutes the feature commit. Follow-up: if commit hygiene matters, split it out; otherwise ignore.
- **var_calling_from_bam.rs:517-522** — the per-chrom path discards all Stage-1 `FilterCounts`/`BaqSkipCounts` (no cross-worker reducer); the code comment already flags this as a follow-up. Pre-existing observability gap, widened in relevance by M3.
- **baq_stream.rs:245-249 / var_calling_from_bam.rs:701** — non-test `.expect(...)` panics (rayon-worker fetcher construction; `pub`-reachable contig-bound). Pre-existing; carry justification comments. Worth converting to typed errors in a separate pass.

## 8. Missing tests to add now

Grouped by function. Names in `function_returns_expected_on_condition` form. Several reliability "Blocker" gaps are **already covered indirectly** by the engine-level tests (which call `left_align_cigar` with `read_start > 0` via the real BAQ window, e.g. `engine_left_aligns_homopolymer_deletion`); those are listed as direct-unit-test gaps, not behavior gaps.

**`build_cigar` (M7 — genuinely uncovered):**
- `build_cigar_orders_deletion_before_insertion` — input `[Match(2), Insertion(1), Deletion(1), Match(2)]` → `[Match(2), Deletion(1), Insertion(1), Match(2)]`. Catches a reorder bug that splits identical events.
- `build_cigar_merges_adjacent_identical_ops` — `[Match(1), Match(4)]` → `[Match(5)]`; `[Deletion(1), Deletion(2)]` → `[Deletion(3)]`.
- `left_align_cigar_merges_colliding_indels` — ref/read where an insertion and deletion roll into the same repeat block; assert a single merged canonical run.

**`left_align_cigar` / `normalize_alleles`:**
- `left_align_cigar_returns_input_on_malformed_read` — (M1) read-consume ≠ `seq.len`; assert input returned unchanged.
- `left_align_cigar_handles_nonzero_read_start` — direct unit test with left context (`read_start = 3`); covered indirectly by engine tests, lock it directly.
- `normalize_alleles_shifts_insertion_left_in_repeat` — empty (insertion) ref range + non-trivial `max_shift`; assert returned shifts/bounds.
- `left_align_cigar_is_idempotent` — `proptest` over random homopolymer/STR contexts: `left_align(out) == out` and mismatch count preserved.

**`count_mismatches` (debug guard):**
- `count_mismatches_counts_known_mismatches` — `count_mismatches(&[Match(4)], b"ACGT", b"AGGT", 0) == 1`.
- `count_mismatches_is_invariant_under_left_alignment`.
- `count_mismatches_early_returns_on_short_window` — lock/repair the Mi5 truncation behavior.

**`BaqEngine::process` / `BaqStream`:**
- `engine_applies_baq_and_left_aligns_same_indel_read` — `apply_baq=true` indel read in a repeat; assert both the rewritten CIGAR and `bq_baq.len() == seq.len()` with caps applied.
- `no_baq_indel_read_past_chrom_end_is_skipped` — documents the new `--no-baq` drop behavior (M2/M3).
- `no_baq_indel_free_read_does_not_require_reference` — proves the fast path skips the fetch.
- `stream_with_no_baq_left_aligns_and_keeps_raw_qual` — `BaqStream` with `apply_baq=false`, indel read → normalized CIGAR + raw `bq_baq`.
- `stream_bails_chunk_when_ref_id_out_of_range` — chunk with `ref_id` absent from `ContigList`; assert all skipped and `chrom_id_out_of_range` count.
- `stream_splits_chunk_on_chrom_transition` — interleave two `ref_id`s; assert order preserved and each read scored against its own contig.
- `engine_skip_chrom_id_out_of_range` — `ref_id > u32::MAX` (it is `usize`); assert `Skipped(ChromIdOutOfRange)`.

## 9. What's good

- **Concurrency is clean by construction** — `raw_ref` is per-worker `BaqEngine` scratch, `apply_baq` is a `Copy` captured immutably, and the rayon `par_drain`+`collect_into_vec` order-preservation is untouched; no new shared state, `unsafe`, or locks. ([baq_stream.rs:241-261](../../../../src/pileup/per_sample/baq_stream.rs#L241))
- **Reference window is genuinely reused** — one `fetch` fills `encoded_ref` (HMM) and/or `raw_ref` (norm) from the same bytes; no second fetch, and the indel-free fast path skips both. ([baq_engine.rs:177-204](../../../../src/pileup/per_sample/baq_engine.rs#L177))
- **`BaqConfig` is spelled field-by-field and every `CigarOp` match is exhaustive (no `_`)** — a new config field or CIGAR variant is a compile error at each site; the parity-drift hazard the PR worried about is actually closed. ([baq_engine.rs:214-218](../../../../src/pileup/per_sample/baq_engine.rs#L214))
- **Output determinism holds** — no `HashMap` iteration, pure per-read transform, order-preserving stream; the `.psp`/VCF stays reproducible.
- **The port is faithful and well-commented** — `left_align_cigar` mirrors GATK's right-to-left structure with per-step comments and cross-references to the spec, easing future audit against the reference implementations.

## 10. Commands to re-verify

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`  (currently 996 passed / 0 failed)
- New: the §8 tests once added; a `benches/` criterion fixture for the normalization hot path (Mi2).

### Author response convention
Address each finding by id (M1, Mi3, …) with `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer §4 open questions 1–4 first — M1's severity and M5's scope depend on them.
