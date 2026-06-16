# Code Review: segment_reader
**Date:** 2026-06-16
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** PR diff on branch `segment-read-fetcher` (commit `e48b452`) — the new `src/bam/segment_reader.rs` segment-read primitive plus its supporting visibility/error edits
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** PR diff (the single commit `e48b452` vs `main`).
- **Reviewed against:** branch `segment-read-fetcher` @ `e48b452`.
- **In-scope files:**
  - [src/bam/segment_reader.rs](../../../../src/bam/segment_reader.rs) — the new module (primary), including its `#[cfg(test)]` suite.
  - [src/bam/errors.rs](../../../../src/bam/errors.rs) — the two new variants `MissingCramReference`, `InvalidSegment`.
  - [src/bam/mod.rs](../../../../src/bam/mod.rs) — module registration.
  - [src/bam/alignment_input.rs](../../../../src/bam/alignment_input.rs), [src/bam/bam_input.rs](../../../../src/bam/bam_input.rs), [src/bam/cram_input.rs](../../../../src/bam/cram_input.rs) — `pub(super)` visibility lifts only.
  - [PROJECT_STATUS.md](../../../../PROJECT_STATUS.md) — status block + current-focus edit.
- **Deliberately out of scope:** the unchanged function bodies of the sibling decoders (`OwnedIndexed{Bam,Cram}Records` etc. — pre-existing, reviewed); the deferred consumers (#4 SSR fetcher, #5 SNP `--regions` retrofit); `Cargo.toml` (unchanged).
- **Categories dispatched:** reliability (always), errors (always), naming (always), defaults (config-acting values), idiomatic (always), refactor_safety (always), module_structure (multi-file), unsafe_concurrency (`Mutex`/`Arc`/`Sync`/rayon), smells (always), extras (PR + hot path + file input), tooling (crate).

## 2. Verdict

**Approve-with-changes.** No Blockers. The primitive's concurrency design is sound (verified: `Sync` for the right reasons, no lock held across I/O, no handle leak on the error paths) and the happy path is well tested on both formats. The Major findings are robustness/maintainability/perf/docs issues — chiefly two that could later cause *silent* wrong results (M3 BAM/CRAM filter divergence, M5 silently-ignored config fields), one perf gap that undercuts the primitive's stated purpose (M4), one extensibility-safety regression vs a documented sibling convention (M2), and a convergent poison-handling asymmetry (M1).

## 3. Execution status

- `cargo fmt --check` — exit 0, clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0, clean.
- `cargo test --lib --tests` — 1056 lib tests pass (1 ignored) + all integration suites pass (3+11+6+5+10+7+4+0).
- `cargo doc --no-deps` — exit 0.
- `cargo audit` — **not run** (cargo-audit not installed on this machine). Diff adds no dependencies, so the advisory surface is unchanged; listed under §10.
- Findings labeled "Needs verification": 0 (all cited locations were read; one perf claim, M4, is reasoned from code structure not a benchmark).

## 4. Open questions and assumptions

1. **Is this primitive *meant* to apply the mismatch-fraction filters (`max_read_mismatch_fraction` / `mismatch_bq_floor`)?** It accepts the full `AlignmentMergedReaderConfig` but applies only the `classify_pre_decode` gates + `min_read_length`, matching the merged reader (which also defers F1/F3 to the downstream read-processing stage). If "no, the consumer applies them" is intended, M5 is satisfied by documenting it at the call site; if "yes", it is a correctness gap. Affects **M5**, **M3** (cross-cut).
2. **Is the CRAM full-contig crai walk an accepted deferral?** The plan §8 defers the CRAM *container cache*, but the per-record *early-stop* (M4) is a separate, ~6-line change. Affects **M4**.
3. **Where should impl/review reports live — `ia/reports/` or `doc/devel/reports/`?** The feature-implementation skill says `ia/reports/`, but every PROJECT_STATUS link and the majority of reports use `doc/devel/reports/`. The impl report was saved to `ia/` and linked to `doc/`, so the link is broken. Affects **M6**.

## 5. Top 3 priorities

1. **M3** — extract the duplicated BAM/CRAM per-record filter cascade into one chokepoint, so a future filter change can't silently diverge the two formats' read sets. [Findings → M3](#m3)
2. **M1** — reconcile the Mutex-poison policy (`borrow_handle` panics, `return_handle` silently drops the reader). Convergent across ~10 categories. [Findings → M1](#m1)
3. **M4** — add the CRAM per-record early-stop so tiny-locus queries don't scan the contig's crai tail on every one of ~10⁶ calls. [Findings → M4](#m4)

## 6. Findings

### Major

<a id="m1"></a>
- [src/bam/segment_reader.rs:330](../../../../src/bam/segment_reader.rs#L330) (+ :343, :533, :546) — **[Major]** Mutex-poison handling is asymmetric and untested
- **Confidence:** High
- **Categories:** unsafe_concurrency, errors, reliability, defaults, idiomatic, refactor_safety, smells, naming, extras, tooling (convergent — surfaced by 10 of 11 sub-agents)
- **Problem:** `borrow_handle` (BAM [:330](../../../../src/bam/segment_reader.rs#L330), CRAM [:533](../../../../src/bam/segment_reader.rs#L533)) does `.lock().expect("reader pool not poisoned")` — it *panics* on a poisoned pool. `return_handle` (BAM [:343](../../../../src/bam/segment_reader.rs#L343), CRAM [:546](../../../../src/bam/segment_reader.rs#L546)) does `if let Ok(mut pool) = self.readers_pool.lock()` — it *silently drops* the reader on poison. Today the lock only guards `Vec::pop`/`push`/`len` (no unwinding work under it), so poisoning is effectively unreachable — but the two halves disagree, and a future edit doing fallible work under the lock would arm it: one poison would then make every new query panic *and* every dropping iterator silently shrink the pool. No test pins either behavior.
- **Why it matters:** This is a shared, `Sync`, many-threaded primitive; a one-panic-into-cascade hazard plus a silent reader leak is exactly the latent concurrency bug worth closing while the surface is small. The asymmetry also contradicts the module's own "typed-error-not-panic" stance.
- **Suggested fix:** Adopt one policy on both sides — recover the guard rather than panic or discard:
  ```rust
  let mut pool = self.readers_pool.lock().unwrap_or_else(|e| e.into_inner());
  ```
  on both borrow and return (poison becomes non-fatal; handles are never silently lost). Add a `// PANIC-FREE:`-style comment naming why a poison is recoverable here, and a test that pins the chosen behavior.

<a id="m2"></a>
- [src/bam/segment_reader.rs:205](../../../../src/bam/segment_reader.rs#L205) — **[Major]** `from_input`'s `(kind, index) =>` catch-all silently absorbs a new `#[non_exhaustive]` variant
- **Confidence:** High
- **Categories:** refactor_safety
- **Problem:** `AlignmentIndex` is `#[non_exhaustive]` ([index_preflight.rs:94](../../../../src/bam/index_preflight.rs#L94)). `from_input` matches the three valid `(kind, index)` pairs, then collapses *everything else* into one `(kind, index) =>` arm returning `AlignmentIndexFormatMismatch`. The sibling `AlignmentMergedReader::query` ([alignment_input.rs:964-981](../../../../src/bam/alignment_input.rs#L964-L981)) handles the identical enum the opposite way — it enumerates each genuine mismatch pair explicitly, with a comment stating it does so *precisely* to force a new compile arm when a variant is added. The new code reintroduces the anti-pattern the existing code documents itself as avoiding: a fourth index variant would compile clean and be reported as a "format mismatch" at runtime instead of being served.
- **Why it matters:** This primitive is the planned future home of both the SSR fetcher and the SNP `--regions` path; a new index format failing with a misleading error and no compiler help is a real extensibility trap, and two functions over the same enum now disagree on safety.
- **Suggested fix:** Mirror `query`'s arm structure — enumerate the genuine mismatch pairs (`(Cram, BamCsi|BamBai) | (Bam, Crai)`) explicitly with `index @ ...` bindings so the wildcard disappears and a new variant forces a new arm. Full snippet in `tmp/review_2026-06-16_segment-reader/refactor_safety.md`.

<a id="m3"></a>
- [src/bam/segment_reader.rs:434-472](../../../../src/bam/segment_reader.rs#L434-L472) and [:663-687](../../../../src/bam/segment_reader.rs#L663-L687) — **[Major]** The per-record filter cascade is copy-pasted between the BAM and CRAM iterators with no shared chokepoint
- **Confidence:** High
- **Categories:** refactor_safety, smells, reliability, errors, idiomatic (cross)
- **Problem:** Both iterators run the same cascade per record — contig check → `segment.overlaps_record` → `classify_pre_decode` → `min_read_length` → `record_buf_to_mapped_read` with the same `MalformedRecord` construction — in two textually-independent copies (BAM in `next`, CRAM in `refill`). They already differ subtly (BAM has the sorted early-stop; CRAM doesn't), which is exactly where copy-paste drift hides. A future change — wiring in the currently-ignored mismatch filters (see M5), or reordering a gate — applied to one arm and not the other would compile clean and make BAM and CRAM return *different reads for the same data*. Nothing cross-checks BAM-result == CRAM-result.
- **Why it matters:** Silent BAM/CRAM divergence in a read source feeding a variant/SSR caller is a wrong-results class of bug a refactor would introduce with no compile error. It also pushes both functions over the length/nesting guidelines (`next` ~93 lines, `refill` ~100 lines / 4 levels) — both shrink under threshold once the cascade is extracted.
- **Suggested fix:** Extract one pure `classify_segment_record(cfg, target, segment, path, source_file_index, &record) -> Result<Option<MappedRead>, AlignmentInputError>` that both arms call (`Ok(Some)` keep / `Ok(None)` drop / `Err` malformed). Keep the BAM-only sorted early-stop in the BAM arm; share only the filter tail. Worth doing now, before the #5 convergence with `OwnedIndexed*Records`.

<a id="m4"></a>
- [src/bam/segment_reader.rs:589-692](../../../../src/bam/segment_reader.rs#L589-L692) — **[Major]** CRAM path has no early-stop; every call scans all target-contig containers to end-of-contig
- **Confidence:** High
- **Categories:** reliability, extras
- **Problem:** `BamSegmentReads::next` latches `done` once `alignment_start > segment.end` ([:441](../../../../src/bam/segment_reader.rs#L441)) — the sort-order stop. `CramSegmentReads::refill` has no equivalent: it walks the `.crai` cursor to the end of the index, and a container that *starts* after `segment.end` is `continue`-skipped ([:614](../../../../src/bam/segment_reader.rs#L614)) rather than stopping the walk. For a tiny SSR locus on a large contig, each of ~10⁶ calls iterates every later target-contig crai record. The work is cheap per record (a `crai::Record` clone + integer compare, no decode) but it is O(containers-after-locus) per call vs the BAM path's O(1)-after-stop — directly opposed to the primitive's stated "10⁶ tiny loci without re-scanning" purpose.
- **Why it matters:** On a human-genome crai (thousands of containers/contig) this is a measurable per-call tax the BAM path does not pay, in the exact workload the primitive exists to serve.
- **Suggested fix:** In `refill`, once an index record's container start exceeds `segment.end`, `return Ok(true)` (the crai is coordinate-ordered on a contig, so nothing later overlaps) instead of `continue`. This subsumes the existing `container_start > segment.end` skip arm. Snippet in `tmp/review_2026-06-16_segment-reader/extras.md`. (The plan §8 defers the CRAM container *cache*; this early-stop is a distinct, cheap change — confirm against Open question 2 before deferring.)

<a id="m5"></a>
- [src/bam/segment_reader.rs:160-167](../../../../src/bam/segment_reader.rs#L160-L167), [:213-245](../../../../src/bam/segment_reader.rs#L213-L245), [:449-459](../../../../src/bam/segment_reader.rs#L449-L459), [:669-677](../../../../src/bam/segment_reader.rs#L669-L677) — **[Major]** Silent partial-config: `max_read_mismatch_fraction` / `mismatch_bq_floor` are accepted but ignored, undocumented at the call site
- **Confidence:** High
- **Categories:** defaults
- **Problem:** `from_input` takes the whole `AlignmentMergedReaderConfig` by value, but the per-record pipeline honors only `min_mapq`/`drop_qc_fail`/`drop_duplicate` (via `classify_pre_decode`) + `min_read_length`. `max_read_mismatch_fraction` and `mismatch_bq_floor` are silently discarded. The omission is explained in the *module header* ([:40-48](../../../../src/bam/segment_reader.rs#L40-L48)) but not on `from_input` or `get_reads_from_segment` — the call-site surface. A caller reusing the same config type they pass to `AlignmentMergedReader` (the natural thing) gets a different read set with no signal.
- **Why it matters:** An honored-subset-of-a-richer-config is a hidden default (the effective `max_read_mismatch_fraction` is forced "off" regardless of input) — a correctness-relevant surprise for the SSR consumer, invisible at every call site. See Open question 1.
- **Suggested fix:** Add a `# Filtering` doc section to both `from_input` and `get_reads_from_segment` naming the honored and ignored fields explicitly. Stronger: a `debug_assert!`/typed rejection when an ignored field is set to a non-inert value. (If the intent is that this type is the wrong surface, consider a narrower segment-reader config exposing only the honored fields.)

<a id="m6"></a>
- [PROJECT_STATUS.md](../../../../PROJECT_STATUS.md) — **[Major]** The impl-report link is broken; the report is stored off the project's report-directory convention
- **Confidence:** High
- **Categories:** extras (surfaced as "impl report does not exist")
- **Problem:** PROJECT_STATUS links the report at `doc/devel/reports/implementations/segment_read_fetcher_2026-06-16.md`, but the file is actually at [ia/reports/implementations/segment_read_fetcher_2026-06-16.md](../../../../ia/reports/implementations/segment_read_fetcher_2026-06-16.md) (verified: absent at the linked path, present at the `ia/` path). Every other PROJECT_STATUS report link and the review-skill's own path references use `doc/devel/reports/`. The `extras` sub-agent, following the link, concluded the report did not exist at all and could not verify the documented assumptions.
- **Why it matters:** PROJECT_STATUS is the project's navigation index; a broken link makes the report (which documents the deferral decisions behind M4 and the silent choices behind M5) unreachable, and the off-convention location means future tooling/readers won't find it.
- **Suggested fix:** `git mv ia/reports/implementations/segment_read_fetcher_2026-06-16.md doc/devel/reports/implementations/` so the file matches the PROJECT_STATUS link and the dominant convention (then the link resolves). Reconcile the feature-implementation skill's `ia/reports/` instruction with the project's `doc/devel/reports/` reality separately.

### Minor

<a id="mi1"></a>
- [src/bam/segment_reader.rs:401](../../../../src/bam/segment_reader.rs#L401) (+ :413, :423, :623) — **[Minor]** `expect("handle held")` repeated 5× with no `// PANIC-FREE:` rationale
- **Confidence:** High — **Categories:** errors, smells, reliability, idiomatic
- **Problem:** The invariant (handle is `Some` for the iterator's life, `take`n only in `Drop`, and `next`/`refill` never run during `Drop`) is genuinely sound, but the call sites duplicate the bare `expect` with no comment, and `idiomatic` notes the `Option<Handle>` exists solely to satisfy `Drop`.
- **Suggested fix:** A small private `fn reader(&mut self) -> &mut ...` accessor that centralizes the `expect` + carries one `// PANIC-FREE:` comment; or a `PooledReader` guard wrapper that owns the `Option`/`Drop` and lets the iterators hold the reader directly (removes the per-iterator `Drop` impl too).

<a id="mi2"></a>
- [src/bam/segment_reader.rs:58](../../../../src/bam/segment_reader.rs#L58) — **[Minor]** Module-wide `#![allow(dead_code)]` is broader than needed
- **Confidence:** High — **Categories:** tooling, smells, defaults, unsafe_concurrency (cross)
- **Problem:** The rationale comment is exemplary, but the blanket scope will also silence *genuinely* dead helpers added to the file later. The codebase already uses a narrower form that still lets tests police dead code: `#[cfg_attr(not(test), allow(dead_code))]` ([src/psp/reader.rs:339](../../../../src/psp/reader.rs#L339), [src/vcf/sink.rs:41](../../../../src/vcf/sink.rs#L41)).
- **Suggested fix:** Swap to `#[cfg_attr(not(test), allow(dead_code))]` (file- or item-level), and remove it entirely when the #4 consumer lands.

<a id="mi3"></a>
- [src/bam/segment_reader.rs:231](../../../../src/bam/segment_reader.rs#L231) — **[Minor]** `get_reads_from_segment` uses the non-idiomatic `get_` prefix and diverges from the sibling verb `query`
- **Confidence:** Medium — **Categories:** naming
- **Problem:** Rust API guidelines (C-GETTER) discourage `get_`; more importantly, the same module's `AlignmentMergedReader::query` is the established verb for "index-seek and return reads for a region", and this primitive is slated to back that very path at #5. No current callers, so renaming is free now.
- **Suggested fix:** Rename to `query_segment` (or `reads_in_segment`); update the two private per-format twins and the tests.

<a id="mi4"></a>
- [src/bam/segment_reader.rs:312-355](../../../../src/bam/segment_reader.rs#L312-L355) — **[Minor]** The reader pool grows to peak concurrency and never shrinks
- **Confidence:** High — **Categories:** unsafe_concurrency
- **Problem:** Each pooled handle holds an open FD + decode buffers for the `AlignmentFile`'s lifetime. Fine for the documented "one segment per thread" SSR driver, but the primitive is `pub(crate)` and could be misused.
- **Suggested fix:** A one-line doc note on `AlignmentFile`/the pool field stating the resting size equals peak concurrent borrows and persists for the file's life.

<a id="mi5"></a>
- [src/bam/segment_reader.rs:310](../../../../src/bam/segment_reader.rs#L310) (+ :624, :632 vs the `io_error` helper at :377/:580) — **[Minor]** `Io`/`MalformedRecord` constructed inline in several places, bypassing the `io_error` helper
- **Confidence:** High — **Categories:** refactor_safety, smells
- **Problem:** Two spellings of the same `Io { path, source }` error coexist (the helper, and inline forms in `get_reads_from_segment` and `refill`). Not a field-shape hole, but invites a future edit touching only some copies.
- **Suggested fix:** Route all `Io` construction through one helper (a free `fn io_err(path, source)` or a `path`-closure usable before the iterator exists).

<a id="mi6"></a>
- [src/bam/segment_reader.rs:439](../../../../src/bam/segment_reader.rs#L439) (+ :607, :674) — **[Minor]** `Position → u64` cast spellings vary, and the `u32`-end boundary is untested
- **Confidence:** Medium — **Categories:** idiomatic, reliability
- **Problem:** Three spellings of the lossless `usize::from(pos) as u64` cast coexist across the module (and `.get() as u64` / `usize::from(first) as u64` in the sibling). No test exercises a segment end at `u32::MAX` or a read footprint at the contig end.
- **Suggested fix:** Unify on one spelling/helper; add a boundary test (segment end = contig length; read ending exactly at the segment end).

<a id="mi7"></a>
- [src/bam/segment_reader.rs:606-617](../../../../src/bam/segment_reader.rs#L606-L617) — **[Minor]** CRAM container-skip predicate untested; `span == 0` silently disables the skip
- **Confidence:** Medium — **Categories:** reliability
- **Problem:** The `if span > 0` guard means a zero-span crai record falls through to a full seek + decode (correctness preserved by the per-record filter, but the optimization is silently off). The only test that touches the skip passes regardless of its correctness (overlapping containers). 
- **Suggested fix:** Add `cram_stops_or_completes_with_downstream_containers` (reads clustered at `[1,20]` and `[150,170]` to force separate containers; query `chr1 5 15`; assert only the first cluster returns). See §8.

<a id="mi8"></a>
- [src/bam/segment_reader.rs:139-211](../../../../src/bam/segment_reader.rs#L139-L211) — **[Minor]** `from_input` trusts the `path`↔`header` pairing without validation
- **Confidence:** Medium — **Categories:** extras
- **Problem:** The doc calls the header "already-validated", but nothing ties it to `path`. A future #4/#5 caller that mispairs them would get wrong reads (wrong ref-id resolution) silently.
- **Suggested fix:** Document the precondition prominently now; when the consumer lands, consider re-reading the header from `path` (cheap, once) or asserting `@SQ` agreement.

<a id="mi9"></a>
- [src/bam/segment_reader.rs:781](../../../../src/bam/segment_reader.rs#L781) — **[Minor]** Test-only back-reference into the pipeline-stage module `crate::pileup::per_sample::cram_files`
- **Confidence:** High — **Categories:** module_structure
- **Problem:** The CRAM fixtures import from a `pileup`-stage module — a peer→stage coupling. It is `#[cfg(test)]`-only and the house convention ([cram_input.rs:425](../../../../src/bam/cram_input.rs#L425), [alignment_input.rs](../../../../src/bam/alignment_input.rs)), so not a blocker.
- **Suggested fix:** Crate-wide follow-up: relocate the CRAM/FASTA fixture builders to a `crate::bam` test-support module. Not for this PR.

### Nits

Grouped (no individual action required beyond a small pass):
- `cfg` field/param abbreviates where the owning type spells it `config` everywhere in the sibling file — align on `config`.
- Match-arm binding inconsistency: `Self::Bam(file)` in `get_reads_from_segment` dispatch vs `Self::Bam(reads)` in the `Iterator` dispatch — pick one role-name.
- The raw noodles `reader` field could carry the library-name prefix per the project convention for raw dependency-type bindings at layer transitions.
- `io_error` is a noun for an error-*constructing* method (constructor-naming nit).
- [src/bam/errors.rs:135-142](../../../../src/bam/errors.rs#L135-L142) — `ContigNotInList`'s doc attributes the error to `AlignmentMergedReader::query`, but `resolve_segment` now raises it too; refresh the attribution.
- `CramSegmentReads::refill -> Result<bool>` — the `bool`'s meaning (`true` = exhausted) lives only in the doc comment; a 2-variant enum or a comment at the call site would be clearer.

## 7. Out of scope observations

- The pre-existing `AlignmentInputError` enum-wide pattern of concatenating `{source}` into `Display` (e.g. `Io`, `OpenFailed`, `MalformedRecord`) duplicates the `#[source]` chain in the message — pre-existing, not introduced here; the two *new* variants correctly avoid it. Follow-up in a separate errors-polish pass if desired.
- The post-open cursor-position contract on the lifted openers (`open_{bam,cram}_reader_with_header` must leave the reader at a record/container boundary) is now load-bearing across two modules but enforced only by docs — note for the #5 convergence when these helpers may be refactored.

## 8. Missing tests to add now

From the `reliability` challenge-tests pass (specs/bodies in `tmp/review_2026-06-16_segment-reader/reliability.md`):

- `get_reads_from_segment_does_not_leak_handle_on_index_query_error` — BAM index `.query()` returns `Err` (the early-return path *before* `borrow_handle`); pins the "borrow only after the fallible query" ordering so a refactor can't leak a pooled reader. If `query` can't be made to fail deterministically, document the invariant at the call site instead.
- `bam_next_returns_none_after_error_is_fused` — drive `BamSegmentReads::next` to a `Some(Err(_))` (truncated bgzf after header, or a `record_buf_to_mapped_read`-rejecting record) and assert the next poll is `None`; pins the `done = true` fuse the happy-path tests never execute.
- `bam_min_read_length_filtering_everything_yields_empty_and_returns_handle` — `min_read_length` exceeds every read; assert empty result and `pool_len == 1` (distinct from the empty-*segment* test).
- `cram_min_read_length_filtering_everything_yields_empty` — same for the CRAM decode path (its own copy of the filter — guards M3 drift).
- `cram_unknown_contig_and_invalid_segment_are_rejected` — the CRAM dispatch has no error-path tests; mirror the BAM ones and assert `pool_len == 0` (nothing borrowed on the rejected path).
- `cram_stops_or_completes_with_downstream_containers` — reads at `[1,20]` and `[150,170]` in separate containers; query `chr1 5 15`; assert only the first cluster returns (exercises the keep-walking path and motivates M4).
- `resolve_segment_accepts_single_base_segment` / boundary unit test — `start == end`, `start == 1`, and a segment end at the contig length.

## 9. What's good

- The `Sync` design is verified correct *for the right reasons* — the only interior mutability is the pool `Mutex`, no lock is held across I/O, and a `assert_sync::<AlignmentFile>()` compile-time test guards it ([segment_reader.rs](../../../../src/bam/segment_reader.rs)).
- No handle leak on any error path: the BAM index query runs *before* the handle is borrowed, and `Drop` returns the handle on every iterator-construction success path including the error-fused path — a genuinely careful ownership design.
- Reuse over forking: the per-record decode delegates to the existing `classify_pre_decode` / `record_buf_to_mapped_read` / `query_interval` helpers rather than reimplementing noodles, with visibility lifts kept to `pub(super)`.
- All call-site failure modes are typed (`MissingCramReference`, `InvalidSegment`, `ContigNotInList`, `AlignmentIndexFormatMismatch`) with context, and the iterators fuse after the first error on every path.
- `permissive_config()` is an exhaustive struct literal (no `..Default`), so a new config field forces a test update — the right refactor-safety pattern (in contrast to several spots where M2/M3 want the same compiler-forcing applied to production matches).

## 10. Commands to re-verify

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib --tests` (or the targeted `cargo test --lib segment_reader`)
- `cargo doc --no-deps`
- **Author must run:** `cargo install cargo-audit && cargo audit` (could not run here; no deps changed, but it is a CI gate).

### Author response convention
Address each finding by identifier (M1–M6, Mi1–Mi9) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the §4 open questions first (they gate M4/M5/M6).
