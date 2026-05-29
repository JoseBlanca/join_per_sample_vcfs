# Code Review: indel_normalization — F3 replacement

**Date:** 2026-05-29
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** the diff that replaces the pre-existing "F3" indel left-aligner with the GATK port
**Status:** Approve-with-changes (most introduced defects fixed during review — see §6)

---

> Companion to the superseded [first review](indel_normalization_2026-05-29.md).
> That one reviewed an earlier BaqEngine-prep implementation and **missed**
> that the BAM/CRAM input cascade already left-aligned every read (the "F3"
> pass). This review covers the restructure (`aaf8685`): F3's single-pass
> shifter is deleted and its call site now invokes the GATK port
> `indel_norm::left_align_indels`.

## 1. Scope

- **What:** the F3-replacement diff vs `main` (the indel-normalization feature as it now stands).
- **In-scope files:**
  - [src/pileup/walker/indel_norm.rs](../../../../src/pileup/walker/indel_norm.rs) — the GATK port, now the **sole** indel-left-alignment implementation.
  - [src/pileup/walker/indel_norm/tests.rs](../../../../src/pileup/walker/indel_norm/tests.rs)
  - [src/pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs) — `pub(crate) mod indel_norm;`
  - [src/bam/alignment_input.rs](../../../../src/bam/alignment_input.rs) — the F3 deletion (~228 fn + ~358 test lines), the import (line 26), the call site (~1070-1082), and `fetch_ref_for_read` (~1136-1176). Rest of this file is pre-existing/out of scope.
- **Out of scope:** docs; [posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs) (incidental pre-existing fmt fix).
- **Categories dispatched:** reliability, errors, refactor_safety, extras, module_structure. (Style categories — naming/idiomatic/smells — were exhaustively run on `indel_norm.rs` in the first review and are unchanged here; their surviving Minors carry over. defaults/unsafe_concurrency/tooling N/A: no config/concurrency/build changes, and the reader merge loop is single-threaded.)

## 2. Verdict

**Approve-with-changes.** The port correctly replaces F3: the call-site integration is sound (`ref_seq` starts at the read's first aligned base, so `read_start = 0` is right; verified in `fetch_ref_for_read`), the reference fetch is reused from the F1 filter, F3 is fully deleted with no dead helpers, and the suite is green. The review surfaced two **defects this PR introduced** (a broken doc link, a now-false F3 comment) and one **resolved test gap** — all **fixed during the review** (§6). The substantive open items are a reader-level integration test (M2), a hot-path allocation regression + missing bench (M4), the module placement question (M5), and a behavioral question about port-produced adjacent indels vs the G2 filter (M6).

## 3. Execution status

- `cargo fmt --check`: exit 0, clean.
- `cargo clippy --lib --all-features -- -D warnings`: clean.
- `cargo test --lib`: 935+ pass (21 in `indel_norm` after the new wrapper tests).
- `cargo test --test '*'`: 46 pass.
- `cargo doc --no-deps`: the **new** in-scope broken link the review found (`[left_align_prepared]`) is **fixed**; the remaining errors are the 13 pre-existing project-wide ones.
- `cargo audit`: not installed.
- "Needs verification": 1 (M6 reachability).

## 4. Open questions

1. **Can the port produce an adjacent `D`/`I` pair that reaches the walker, and does the walker handle it?** The port merges colliding indels (F3 did not); G2 `cigar_is_bad` runs *before* the port, so a port-created adjacent pair is not re-checked. The reliability sub-agent could not construct a reachable non-cancelling collision (attempts trimmed to pure `M`). Affects **M6**.
2. **Where should `indel_norm` live** now its only consumer is `src/bam/`? Affects **M5**.

## 5. Top 3 priorities

1. **M2** — add a reader-level integration test (`AlignmentMergedReader` → `fetch_ref_for_read` → `left_align_indels`); a slice-offset/`read_start` mis-wiring would pass every current unit test.
2. **M4** — the port allocates 4–5 `Vec`s per indel read on the single-threaded reader loop where F3 was zero-alloc in-place; add scratch reuse + a criterion bench.
3. **M6** — confirm/handle the port-produced-adjacent-indel-vs-G2 question with a test.

## 6. Findings

### Fixed during review (introduced by this PR or trivially closable)
- **F1: src/bam/alignment_input.rs:1568 — stale F3 comment asserted "adjacent indels are invariant under F3"**, which is now false (the port merges colliding indels). It was the load-bearing rationale for ordering `cigar_is_bad` before alignment. *Categories: refactor_safety (Major), reliability (Minor).* **Fixed** — comment rewritten to explain rule 1 targets the aligner's original CIGAR and the port may legitimately produce a canonical adjacent pair.
- **F2: src/pileup/walker/indel_norm.rs:546 — broken intra-doc link `[left_align_prepared]`** (renamed to `left_align_indels`); a real `cargo doc` regression in an in-scope file. *Categories: all five.* **Fixed.**
- **F3: indel_norm.rs:289 — `unreachable!` lacked the required `// UNREACHABLE:` comment.** *Categories: errors, reliability.* **Fixed.**
- **M3 (resolved): the `left_align_indels` wrapper and its debug mismatch-count invariant were untested** — every shift test went through `left_align_cigar`, bypassing the wrapper. *Category: reliability (Major).* **Fixed** — added `wrapper_*` tests driving the in-place wrapper (also pins the `seq`/`ref_seq` arg order: a swap trips the invariant).

### Major (open)
- **M2: no end-to-end test through the reader.** *Category: reliability.* The ~15 deleted `f3_*` tests have no integration equivalent; the new tests cover the pure routine + wrapper but nothing exercises the `AlignmentMergedReader` wiring (the `pos-1` offset in `fetch_ref_for_read`, the `None`-on-overshoot skip). A wiring regression would be invisible. **Fix:** add a reader-level test with a small FASTA + a repeat-region indel read, asserting the emitted record's canonical anchor; assert an indel-free read and a window-past-contig-end read are handled.
- **M4: hot-path allocation regression + no perf guard.** *Category: extras.* F3 shifted in place with zero allocation; the port allocates `result_right_to_left`, `nonzero`, `merged`, `replacement` per indel-bearing read on the single-threaded `AlignmentMergedReader::next` loop, and there is no criterion bench for the normalization path. **Fix:** reuse scratch owned by the reader (or reuse `mapped.cigar`'s allocation via `mem::swap`); add a `benches/` fixture with a regression threshold. Carries over the first review's Mi2.
- **M5: `indel_norm` placement.** *Category: module_structure.* Its sole consumer is `src/bam/alignment_input.rs` and its only dependency is the walker's `CigarOp`; the walker never calls it. The sub-agent recommends **moving to `src/bam/indel_norm.rs`** (beside its consumer; no new cross-stage edge, since `bam` already imports `CigarOp` from `walker`), removing one `bam → pileup` edge. A neutral top-level peer is the runner-up if a second consumer appears.
- **M6 (Needs verification): port-produced adjacent indels bypass G2.** *Category: reliability.* The port can merge two G2-approved separated indels into an adjacent `D`/`I` pair that `cigar_is_bad` (run earlier) never re-checks. The cursor has a stacked-indel path and F1 is panic-safe, but this F3-never-produced shape has no end-to-end test. **Fix:** construct a colliding-indel read, confirm the emitted record is correct (or, if such a pair must be rejected, re-run the adjacency check post-alignment).

### Minor (open)
- **Mi1: indel_norm.rs:400 — the read-consume guard returns the CIGAR unchanged with no log/telemetry.** A corrupt-CRAM read silently skips normalization (then the walker length check drops it). Pure module has no logger; consider having the wrapper signal "skipped" so the reader can bump `filter_counts`. *errors.*
- **Mi2: visibility precision** — `left_align_cigar` and `LeftAlignResult` are `pub(crate)` but only `left_align_indels` is consumed externally; demote the other two to private (tests reach them via `super`). *module_structure.*
- **Mi3: test-layout** — `indel_norm/tests.rs` is a directory holding one test file beside the flat `indel_norm.rs`; the tree's convention is inline `mod tests` or a sibling `_tests.rs` (cf. `baq_engine.rs`/`baq_tests.rs`). *module_structure.*
- Carryover from the first review (still apply to the unchanged port): `Range`→`IndexRange` rename; factor the three `*_is_same` predicate helpers; the `build_cigar` splice loop could be a grouping fold.

## 7. Out of scope observations
- **alignment_input.rs:1162 — `expect("just set")` in `fetch_ref_for_read` lacks a `// PANIC-FREE:` comment.** Pre-existing (the F1 fetch helper; this PR only calls it). The invariant holds but is non-local. Suggest the comment or restructuring to return the value directly.
- **alignment_input.rs (~1090) — a read whose window overshoots the contig end is skipped for normalization** (`fetch_ref_for_read` → `None`), leaving the in-bounds portion un-canonicalised. Pre-existing (F3 had the same), and correct for the F1 filter that shares the fetch; note as a minor recall edge.

## 8. Missing tests to add now
- `reader_left_aligns_repeat_indel_end_to_end` — M2; FASTA + repeat-region indel read through `AlignmentMergedReader`; assert canonical anchor.
- `left_align_cigar_handles_nonzero_read_start` — a direct unit test (production uses `read_start = 0`, but the parameter is untested at non-zero).
- `left_align_cigar_multi_indel_each_aligned` — two separated indels in one CIGAR.
- `left_align_cigar_is_idempotent` — proptest over random repeat contexts (idempotence + mismatch-count invariance).
- `build_cigar_three_indel_run` — a 3+-indel run through steps (2)/(3).
- (added this review) `wrapper_*` tests for `left_align_indels`.

## 9. What's good
- **The call-site integration is verifiably correct** — `fetch_ref_for_read` returns `raw[pos-1 .. pos-1+ref_span]`, so `ref_seq[0]` is the read's first aligned base and `read_start = 0` is right; the reference fetch is shared with the F1 filter (no extra I/O). ([alignment_input.rs:1169-1175](../../../../src/bam/alignment_input.rs#L1169))
- **F3 deletion is complete and the port is a faithful drop-in** — no surviving `try_apply_indel_shift`/`m_op_len`/`max_left_shift_*`/`f3_*` references; `left_align_indels` has byte-identical arg order/types to the deleted F3.
- **The opaque-bool hazard is insulated** — `remove_deletions_at_ends` lives only on internal `left_align_cigar`; the wrapper hardcodes `false`, so the reader never sees it.
- **Determinism holds** — no `HashMap` iteration; pure right-to-left walk; reproducible output.
- **Structurally mandatory normalization** — running at the reader (upstream of BAQ and the `--no-baq` split) makes "always normalize" a property of placement, not of a flag.

## 10. Commands to re-verify
- `cargo fmt --check`; `cargo clippy --all-targets --all-features -- -D warnings`; `cargo test --all-targets --all-features`.
- New: the §8 tests; a normalization criterion bench (M4).

### Author response convention
Address each id (F1–F3 fixed; M2/M4/M5/M6; Mi1–Mi3) with `fixed in <commit>` / `disputed` / `deferred` / `won't fix`. Answer §4 first.
