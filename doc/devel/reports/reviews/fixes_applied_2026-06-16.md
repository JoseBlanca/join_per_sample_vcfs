# Fix Application Report: segment_reader_2026-06-16.md

**Date:** 2026-06-16
**Source review:** `doc/devel/reports/reviews/segment_reader_2026-06-16.md`
**Source state reviewed against:** branch `segment-read-fetcher` @ `e48b452`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

All six Major findings are resolved. M4/M5/M6 were applied earlier in the same
interactive session (during the open-question discussion); M1/M2/M3 plus the
two cheap Minors (Mi1/Mi2) were applied in this fix run. The remaining Minors
and Nits are deferred with reasons (naming calls, design decisions tied to the
not-yet-wired consumer, or fixture limitations).

### Review totals
- Blockers: 0
- Majors: 6 (M1–M6)
- Minors: 9 (Mi1–Mi9)
- Nits: grouped (unnumbered)

### Outcome totals
- Applied: 7 (M1, M2, M3, M4, M5, M6, Mi2)
- Applied with adaptation: 1 (Mi1 — field-doc rationale instead of an accessor)
- Already fixed: 0
- Deferred: 7 (Mi3, Mi4, Mi5, Mi6, Mi7, Mi8, Mi9) + Nits
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-targets --all-features` → lib (1059 pass, 1 ignored) + all
  integration suites pass; the only failure is the **pre-existing**
  `benches/psp_writer_perf.rs:386` panic (separate module, documented in
  PROJECT_STATUS, not touched by any fix here)
- `cargo doc --no-deps` → 0, clean
- `cargo audit` → not run (cargo-audit not installed; this run adds no
  dependencies, so the advisory surface is unchanged)
- Performance check → skipped (no `benches/` harness reaches
  `src/bam/segment_reader.rs`; the primitive has no production consumer yet)

### Unresolved high-priority findings
- None (all 6 Majors resolved).

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | Mutex-poison asymmetry | Apply | Applied | No | segment_reader.rs | Pass | No |
| M2 | Major | `from_input` catch-all absorbs new variant | Apply | Applied | No | segment_reader.rs | Pass | No |
| M3 | Major | Duplicated BAM/CRAM filter cascade | Apply | Applied | No | segment_reader.rs | Pass | No |
| M4 | Major | CRAM no early-stop | Apply | Applied | Q answered | segment_reader.rs | Pass | No |
| M5 | Major | Silent partial-config | Apply | Applied | Q answered | segment_reader.rs | Pass | No |
| M6 | Major | Broken impl-report link / report location | Apply | Applied | Q answered | ia/skills/* (main `8c32b00`), report `git mv`, PROJECT_STATUS.md | Pass | Other features' `ia/reports/` migration left as-is |
| Mi1 | Minor | `expect("handle held")` lacks rationale | Apply | Applied with adaptation | No | segment_reader.rs | Pass | Guard-wrapper encapsulation deferred |
| Mi2 | Minor | Module-wide `#![allow(dead_code)]` too broad | Apply | Applied | No | segment_reader.rs | Pass | Remove entirely when #4 consumer lands |
| Mi3 | Minor | `get_reads_from_segment` `get_` prefix / diverges from `query` | Defer | Deferred | No | None | N/A | Naming call; revisit at #5 convergence |
| Mi4 | Minor | Pool never shrinks — doc note | Defer | Deferred | No | None | N/A | Trivial doc; fold into #4 wiring |
| Mi5 | Minor | Inline `Io`/`MalformedRecord` bypass helper | Defer | Deferred | No | None | N/A | `MalformedRecord` dedup subsumed by M3; residual `Io` sites are reader-scoped |
| Mi6 | Minor | Cast-spelling + boundary test | Defer | Deferred | No | None | N/A | Cosmetic; no behavior change |
| Mi7 | Minor | CRAM container-skip `span==0` untested | Defer | Deferred | No | None | N/A | Fixture writer can't force multi-container |
| Mi8 | Minor | `from_input` trusts path↔header pairing | Defer | Deferred | No | None | N/A | Validate when #4 consumer lands |
| Mi9 | Minor | Test-only back-ref into `pileup::per_sample::cram_files` | Defer | Deferred | No | None | N/A | Crate-wide fixture-relocation follow-up |
| Nits | Nit | naming/doc nits | Defer | Deferred | No | None | N/A | Optional |

## 3. Questions asked and answers

1. **M5 (Open question 1)** — Should the primitive apply the reference-dependent
   mismatch filters, or treat them as the consumer's job?
   - **Answer:** Apply only the cheap filters; make the reference-dependent
     fields unrepresentable (option A). Implemented as `SegmentReadFilter`.
2. **M4 (Open question 2)** — Is the CRAM full-contig crai walk an accepted
   deferral, or is the per-record early-stop worth doing now?
   - **Answer:** Do the early-stop now ("stop as soon as we pass the useful
     containers").
3. **M6 (Open question 3)** — Where should reports live?
   - **Answer:** `doc/devel/reports/`. Fixed the four report-writing skills on
     `main` and moved this feature's report.

## 4. Per-finding log

### M1 — Mutex-poison asymmetry
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `borrow_handle` panicked on a poisoned pool (`.expect`) while
  `return_handle` silently dropped the reader (`if let Ok`). Unified both
  (BAM + CRAM) on a `lock_pool()` helper that recovers the guard via
  `PoisonError::into_inner`, so a poison originating elsewhere neither cascades
  into a panic-storm nor silently shrinks the pool.
- **Implementation summary:** Added `BamFile::lock_pool` / `CramFile::lock_pool`
  (recovering `unwrap_or_else(|p| p.into_inner())`); routed
  `borrow_handle` / `return_handle` / `pool_len` through it for both formats.
- **Review suggestion used verbatim?:** No (skill: a finding is not a patch).
- **Adaptation:** Centralized in a `lock_pool` helper rather than inlining
  `into_inner` at six sites.
- **Verification performed:** Test-first — `bam_pool_survives_a_poisoned_lock`
  poisons the pool (panic while holding the guard) then asserts a query still
  returns its reads and the pool returns to resting size. Confirmed the test
  panics against the pre-fix `.expect`, passes after.
- **Files changed:** `src/bam/segment_reader.rs`
- **Tests added or modified:** `bam_pool_survives_a_poisoned_lock` (new)
- **Validation:** `cargo test --lib segment_reader` → 0, pass; full suite below.
- **Follow-up:** None. **Residual risk:** None (lock guards only `Vec` ops).

### M2 — `from_input` catch-all absorbs a new `#[non_exhaustive]` variant
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Replaced the `(kind, index) =>` wildcard with the two genuine
  mismatch pairs enumerated explicitly (`(Cram, BamCsi|BamBai)` and
  `(Bam, Crai)`), mirroring `AlignmentMergedReader::query`. A new
  `AlignmentIndex` variant now forces a new arm instead of being silently
  served as "format mismatch".
- **Implementation summary:** Rewrote the trailing match arm; relies on
  `AlignmentFileKind: Copy` (same as the sibling) so `kind.display_name()`
  stays valid in the arm body.
- **Review suggestion used verbatim?:** Yes (structure matches the sibling).
- **Adaptation:** None.
- **Verification performed:** Existing `from_input_rejects_format_index_mismatch`
  (Bam+Crai) plus new `from_input_rejects_cram_with_bam_index` (Cram+BamCsi)
  cover both arms; compile confirms exhaustiveness.
- **Files changed:** `src/bam/segment_reader.rs`
- **Tests added or modified:** `from_input_rejects_cram_with_bam_index` (new)
- **Validation:** `cargo test --lib from_input_rejects` → 0, 2 pass.
- **Follow-up:** None. **Residual risk:** None.

### M3 — Duplicated per-record filter cascade (BAM vs CRAM)
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Extracted the shared cascade (contig check → overlap →
  `classify_pre_decode` → min-length → convert, with one `MalformedRecord`
  construction) into a free `classify_segment_record(...) -> Result<Option<MappedRead>>`
  that both iterators call. A future filter change now has one site, so BAM and
  CRAM cannot silently diverge. The BAM-only sorted early-stop stays in the BAM
  loop (it ends the scan; CRAM has no per-record sort guarantee).
- **Implementation summary:** New `classify_segment_record`; BAM `next` and
  CRAM `refill` rewired to call it; added `use std::path::Path`. Both
  `next`/`refill` shrank well under the length/nesting guidelines as a result.
- **Review suggestion used verbatim?:** No (adapted the signature to the crate).
- **Adaptation:** Helper returns `Result<Option<_>>` (keep/drop/err) and takes
  `&Path` for the error path; the early-stop is kept out of the helper.
- **Verification performed:** Behavior-preserving refactor guarded by the
  existing overlap / boundary / target-contig / MAPQ / min-length / spanning
  tests on **both** formats, plus the M4 early-stop test — all pass unchanged.
- **Files changed:** `src/bam/segment_reader.rs`
- **Tests added or modified:** None new (existing both-format suite is the net).
- **Validation:** `cargo test --lib segment_reader` → 21 pass.
- **Follow-up:** None. **Residual risk:** None (output verified identical).

### M4 — CRAM has no early-stop
- **Severity:** Major
- **Initial decision:** Apply (user-confirmed)
- **Final status:** Applied
- **Reasoning:** In `refill`, a container whose start is past `segment.end` now
  ends the walk (`return Ok(true)`) instead of `continue`-ing to the end of the
  contig's `.crai`. The `.crai` is coordinate-ordered, so nothing later
  overlaps. Turns the per-call cost from O(index-tail) to O(near-locus),
  matching the BAM path. Container *cache* (plan §8) stays deferred.
- **Implementation summary:** Split the container-disjoint check; added the
  `container_start > segment.end → return Ok(true)` stop before the span check.
- **Verification performed:** `cram_early_stops_when_container_starts_past_segment`
  asserts a before-the-container segment yields nothing and a reaching segment
  still returns the reads (stop not too eager).
- **Files changed:** `src/bam/segment_reader.rs`
- **Tests added or modified:** `cram_early_stops_when_container_starts_past_segment` (new)
- **Validation:** `cargo test --lib segment_reader` → pass.
- **Follow-up:** Binary-search the crai *head* (a separate, measure-first
  optimization; discussed and intentionally not done now). **Residual risk:** None.

### M5 — Silent partial-config
- **Severity:** Major
- **Initial decision:** Apply (user-confirmed, option A)
- **Final status:** Applied
- **Reasoning:** Replaced the `AlignmentMergedReaderConfig` parameter with a new
  `SegmentReadFilter` exposing only the cheap fields (`min_mapq`,
  `min_read_length`, `drop_qc_fail`, `drop_duplicate`). The reference-dependent
  fields are now unrepresentable, so a caller cannot set one and have it
  silently ignored. `classify_pre_decode` is still reused via
  `SegmentReadFilter::pre_decode_config()` (Tier-2 forced inert), and that
  literal breaks if a new `AlignmentMergedReaderConfig` field is added — forcing
  a decision. Added `# Filtering` docs to `from_input` / `get_reads_from_segment`.
- **Implementation summary:** New `SegmentReadFilter` (+ `Default`,
  `pre_decode_config`); `from_input` and the file structs take/store `filter`;
  both iterators read `self.file.filter`. `alignment_input`/merged reader
  untouched (no byte-identity risk).
- **Verification performed:** Existing filter tests adapted to `SegmentReadFilter`
  (incl. the MAPQ-isolation test); all pass.
- **Files changed:** `src/bam/segment_reader.rs`
- **Tests added or modified:** `permissive_config()` + helpers retyped; no new test.
- **Validation:** `cargo test --lib segment_reader` → pass.
- **Follow-up:** None. **Residual risk:** None.

### M6 — Broken impl-report link / report location
- **Severity:** Major
- **Initial decision:** Apply (user-confirmed)
- **Final status:** Applied
- **Reasoning:** The report was at `ia/reports/` but linked from
  `doc/devel/reports/` (the project convention). Fixed the convention at the
  source (all four report-writing skills) on `main`, then moved this feature's
  report and repointed PROJECT_STATUS.
- **Implementation summary:** (1) On `main`, commit `8c32b00` fixes
  `feature_implementation_skill.md`, `code_review_skill.md`,
  `code_review_fixes.md`, `performance_review_skill.md` to save under
  `doc/devel/reports/` (+ relative-link-depth and `ls` chronology updates;
  historical `ia/reviews/perf_baq` citation left intact). (2) On this branch,
  `git mv` the impl report into `doc/devel/reports/implementations/`, deepened
  its internal links by one `../`, and repointed both PROJECT_STATUS links.
- **Verification performed:** Grep confirmed no bare `reviews/`/`ia/reports`
  save-paths remain in the four skills and no double-substitution; the report's
  links are now 4-deep; the two PROJECT_STATUS segment_reader links resolve.
- **Files changed:** `ia/skills/*` (on `main`), report move + PROJECT_STATUS.md
  (this branch).
- **Tests added or modified:** None (docs/tooling).
- **Validation:** N/A (no code). **Follow-up:** Other features' reports still in
  `ia/reports/` are left where their PROJECT_STATUS links point — migrating them
  would touch other features' blocks (out of scope). **Residual risk:** None.

### Mi1 — `expect("handle held")` lacks rationale
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** The review suggested a private accessor or a guard wrapper. A
  `fn reader(&mut self)` accessor borrows *all* of `self`, which breaks the
  disjoint-field borrows `next` relies on (it reads `self.file` while mutating
  `self.handle`). So the safe minimal fix is to document the invariant where it
  lives rather than restructure ownership.
- **Implementation summary:** Expanded the `handle: Option<…>` field docs on
  both iterators to state the `Some`-until-`Drop` invariant and why the
  `expect` cannot fire.
- **Adaptation:** Field-doc rationale instead of an accessor (accessor breaks
  borrows); the guard-wrapper encapsulation is left as a follow-up.
- **Files changed:** `src/bam/segment_reader.rs`
- **Validation:** compile + full suite. **Follow-up:** optional guard wrapper.
  **Residual risk:** None.

### Mi2 — Module-wide `#![allow(dead_code)]` too broad
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Swapped `#![allow(dead_code)]` for
  `#![cfg_attr(not(test), allow(dead_code))]` (the sibling convention in
  `psp/reader.rs` / `vcf/sink.rs`): the non-test build is still suppressed
  (everything is dead without a consumer), but the test build now flags any
  genuinely-dead helper added later.
- **Implementation summary:** One-line attribute change + comment.
- **Verification performed:** `cargo clippy --all-targets -D warnings` still
  clean (the non-test lib compile remains suppressed).
- **Files changed:** `src/bam/segment_reader.rs`
- **Follow-up:** Remove entirely when #4 lands. **Residual risk:** None.

### Mi3 — `get_reads_from_segment` naming
- **Severity:** Minor — **Initial decision:** Defer — **Final status:** Deferred
- **Reasoning:** Naming judgment with no single right answer (the review itself
  offered `query_segment` vs `reads_in_segment`); renaming the primary method
  is best decided when the primitive converges with `AlignmentMergedReader::query`
  at #5. No external callers, so cost of waiting is nil.
- **Files changed:** None. **Follow-up:** decide at #5.

### Mi4 — Pool never shrinks (doc note)
- **Severity:** Minor — **Initial decision:** Defer — **Final status:** Deferred
- **Reasoning:** Trivial doc note; folds naturally into the #4 wiring where the
  resting-size behavior becomes observable to a real caller. **Files changed:** None.

### Mi5 — Inline `Io`/`MalformedRecord` bypass the helper
- **Severity:** Minor — **Initial decision:** Defer — **Final status:** Deferred
- **Reasoning:** M3 already collapsed the duplicated `MalformedRecord`
  construction into `classify_segment_record`. The residual inline `Io` sites
  are reader-scoped (seek / `read_container`) and constructing them through a
  `&self` helper would re-borrow `self` while the reader borrow is live — the
  same borrow conflict as Mi1. Not worth a restructure now. **Files changed:** None.

### Mi6 — Cast-spelling + boundary test
- **Severity:** Minor — **Initial decision:** Defer — **Final status:** Deferred
- **Reasoning:** Cosmetic (`usize::from(p) as u64` spelling) with no behavior
  change; the boundary test is nice-to-have. **Files changed:** None.

### Mi7 — CRAM container-skip `span==0` untested
- **Severity:** Minor — **Initial decision:** Defer — **Final status:** Deferred
- **Reasoning:** The skip is an optimization (correctness preserved by the
  per-record filter); a true multi-container test needs the fixture writer to
  emit small containers, which noodles' writer Builder does not expose.
  **Files changed:** None. **Follow-up:** add a fixture knob.

### Mi8 — `from_input` trusts path↔header pairing
- **Severity:** Minor — **Initial decision:** Defer — **Final status:** Deferred
- **Reasoning:** A precondition for the not-yet-written #4/#5 callers; validating
  (re-reading the header or asserting `@SQ` agreement) is best added when the
  consumer exists. **Files changed:** None.

### Mi9 — Test-only back-ref into `pileup::per_sample::cram_files`
- **Severity:** Minor — **Initial decision:** Defer — **Final status:** Deferred
- **Reasoning:** `#[cfg(test)]`-only and the house convention (`cram_input.rs`,
  `alignment_input.rs` do the same); relocating fixtures is a crate-wide
  follow-up, not this finding's scope. **Files changed:** None.

### Nits
- **Final status:** Deferred. `cfg`-vs-`config`, match-arm binding inconsistency,
  raw `reader` field name, `io_error` constructor naming, stale `query`
  doc-attribution on `ContigNotInList`, `refill -> Result<bool>` meaning. All
  optional; grouped for a future cosmetic pass. **Files changed:** None.

## 5. Deferred findings to carry forward
- Mi3, Mi4, Mi5, Mi6, Mi7, Mi8, Mi9 + Nits — see per-finding reasons above.

## 6. Disputed findings to return to reviewer
- None.

## 7. Failed-validation findings
- None.

## 8. Blocked-by-context-mismatch findings
- None.

## 9. Performance check
- **Triggered:** no.
- **Baseline saved:** no — not needed.
- **Benches run:** none.
- **Outcome:** Skipped — no `Apply` touched perf-sensitive code reachable from a
  `benches/` harness (`segment_reader` has no bench and no production consumer).
- **Notes:** The `--all-targets` run does surface the pre-existing
  `benches/psp_writer_perf.rs:386` panic; it is unrelated (separate module,
  untouched) and documented in PROJECT_STATUS.

## 10. Commands run
- `cargo test --lib segment_reader`
- `cargo test --lib bam_pool_survives_a_poisoned_lock` (test-first, M1)
- `cargo test --lib from_input_rejects`
- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo doc --no-deps`

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-targets --all-features` → 0 for lib (1059 pass / 1 ignored)
  + all integration suites; pre-existing `psp_writer_perf` bench panic only
- `cargo doc --no-deps` → 0, clean
- `cargo audit` → not run (not installed; no dependency changes)

## 12. Notes
- M4/M5/M6 were applied during the interactive open-question discussion that
  preceded this formal fix run; recorded here as Applied with their work.
- M1/M3/Mi1/Mi5 all bump against the same constraint — the iterators rely on
  disjoint-field borrows (`self.file` read while `self.handle` is mutated), so
  any "extract a `&self`/`&mut self` helper" refactor that re-borrows the whole
  receiver is off the table; helpers were shaped (free function for M3, field
  doc for Mi1) to respect that.
