# Fix Application Report: ssr_pileup_2026-06-17.md

**Date:** 2026-06-17
**Source review:** `doc/devel/reports/reviews/ssr_pileup_2026-06-17.md`
**Source state reviewed against:** branch `ssr-pileup-review` @ `9c47399` (working tree at `354802b`; the parallel session added two `driver.rs` parallelism commits — `pair_hmm.rs` / `segment_reader.rs` / `fetch_reads.rs` were unchanged since review, so all targeted findings context-matched cleanly).
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

Severities below are the **review-discussion re-grades** (recorded with the user this session): on learning this is the first working version with **no released `.ssr.psp` baseline and no committed byte-snapshot test** (verified), M1→Minor (reword), M2→Minor (latent; folded into M5), M4→Nit (doc), Mi5→Nit. The must-fix set was **B1, M3, M5** (with M2 folded into M5).

### Review totals
- Blockers: 1
- Majors: 5 (M1–M5; re-graded: M1→Minor, M2→Minor, M4→Nit)
- Minors: 12 (Mi1–Mi12)
- Nits: grouped

### Outcome totals
- Applied: 13 (B1, M2, M5, Mi1, Mi3, Mi4, Mi5, Mi6, Mi7, Mi8, Mi9, Mi10, + doc nits)
- Applied with adaptation: 3 (M1, M3, M4)
- Already fixed: 0
- Deferred: 3 (Mi2, Mi11, Mi12) + the cosmetic name nits
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → **0**, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → **0**, clean (one round of `redundant_closure` fixed: `|e| to_io(e)` → `&to_io`)
- `cargo test --lib` → **0**, `1171 passed; 0 failed; 3 ignored` (+4 new tests vs the 1167 baseline)
- `cargo test --tests` → **0**, all suites pass (lib 1171, main 3, integration 11+6+5+11+7+4+1)
- `cargo doc --no-deps` → **0**, `Generated …/doc/pop_var_caller/index.html` (**B1 gate restored**)
- `cargo audit` → not run (cargo-audit unavailable in-container per CLAUDE.md/PROJECT_STATUS; no new deps this run)
- Performance check → **skipped** — no `Apply` changed hot-path executable code covered by a bench (the only benched file touched, `pair_hmm.rs`, got a test + a comment; all executable edits are in `segment_reader.rs`, which is not reachable from `benches/`).

Note on `--all-targets`: not run to completion (it includes `benches/psp_writer_perf.rs`, whose runtime panic at `:386` is the pre-existing out-of-scope failure documented in PROJECT_STATUS). All targets **compile** (`clippy --all-targets` green, including the edited `benches/ssr_pileup_perf.rs` and `examples/profile_ssr_pileup.rs`); the CI gate `--lib --tests` is green.

### Unresolved high-priority findings
- None. All Blocker/Major findings are Applied or Applied-with-adaptation. One follow-up carried (M3's multi-container *fetch* test — needs a fixture capability; the eviction/dedup gap it shared is closed).

## 2. Findings table

| ID | Severity (re-graded) | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | Broken `cargo doc` — unresolved `AlignmentFile` link | Apply | **Applied** | fetch_reads.rs | doc 0 |
| M1 | Minor (was Major) | "byte-identical to main" math claim overclaims | Apply | **Applied w/ adaptation** | PROJECT_STATUS.md | n/a (doc) |
| M2 | Minor (was Major) | `decoded==0` walk divergence from `refill` | Apply (via M5) | **Applied** | segment_reader.rs | test/clippy 0 |
| M3 | Major | `ContainerCache` eviction + multi-container untested | Apply | **Applied w/ adaptation** | segment_reader.rs | test 0 |
| M4 | Nit (was Major) | cache cap not runtime-inspectable | Apply (doc) | **Applied w/ adaptation** | segment_reader.rs | doc 0 |
| M5 | Major | `.crai` walk + decode duplicated with `refill` | Apply | **Applied** | segment_reader.rs | full suite 0 |
| Mi1 | Minor | stale `window=10` in bench/example/harness | Apply | **Applied** | bench_harness.rs, benches/…, examples/… | clippy 0 |
| Mi2 | Minor | `bench_harness` pub seam not cfg-gated | Defer | **Deferred** | None | n/a |
| Mi3 | Minor | `.expect("opened by ensure_open")` lacks PANIC-FREE note | Apply (via M5) | **Applied** | segment_reader.rs | clippy 0 |
| Mi4 | Minor | single-candidate / `n==lcp` bit-identity test gap | Apply | **Applied** | pair_hmm.rs | test 0 |
| Mi5 | Nit (was Minor) | manual env-var tools mislabeled | Apply | **Applied** | driver.rs, catalog/io.rs | doc 0 |
| Mi6 | Minor | inline `Io` vs `io_error` + misleading comment | Apply (via M5) | **Applied** | segment_reader.rs | clippy 0 |
| Mi7 | Minor | `#[allow(too_many_arguments)]` unjustified | Apply | **Applied** | pair_hmm.rs | clippy 0 |
| Mi8 | Minor | `Arc<Vec<RecordBuf>>` → `Arc<[RecordBuf]>` | Apply | **Applied** | segment_reader.rs | clippy/test 0 |
| Mi9 | Minor | no compile-time `Send` guard on `WorkerReader` | Apply | **Applied** | segment_reader.rs | test 0 |
| Mi10 | Minor | `ContainerCache.cap` vs `Reservoir.capacity` | Apply | **Applied** | segment_reader.rs | clippy 0 |
| Mi11 | Minor | new bench has no regression/baseline guard | Defer | **Deferred** | None | n/a |
| Mi12 | Minor | driver batch-fill error lacks positional context | Defer | **Deferred** | None | n/a |
| Nits | Nit | ensure_open doc / Arc-shared wording / linear-scan note | Apply | **Applied** | segment_reader.rs | doc 0 |
| Nits | Nit | `cands`/`out`/`v` names; mimalloc `psp_writer_perf` doc | Defer | **Deferred** | None | n/a |

## 3. Questions asked and answers

None this run — all policy decisions were settled in the prior review discussion (recorded in §4 of the source review and the session transcript). Key resolutions carried in:
1. **M1** — first working version, no `.ssr.psp` baseline/golden → reword the claim, keep the faster math (Option C). The ~0.5%-class drift is acceptable to the user.
2. **M4** — cache cap is not output-affecting → keep hardcoded, do **not** add to the header (the review's header recommendation was withdrawn); document as intentionally fixed.

## 4. Per-finding log

### B1 — Broken `cargo doc` (unresolved `AlignmentFile` link)
- **Severity:** Blocker · **Initial:** Apply · **Final:** Applied
- **Reasoning:** verified the gate still failed on current code; the diff orphaned the module-doc link when the import moved `AlignmentFile`→`WorkerReader`.
- **Implementation:** re-pointed the `[`AlignmentFile`]` link to `[`WorkerReader`]` and reworded the bullet to describe the decode-caching CRAM / pooled BAM source.
- **Review suggestion used verbatim?:** Yes (intent).
- **Files changed:** `src/ssr/pileup/fetch_reads.rs`
- **Validation:** `cargo doc --no-deps` → 0, `Generated …/index.html` (was a hard error before).
- **Follow-up / Residual risk:** None.

### M1 — "byte-identical to main" math claim overclaims
- **Severity:** Minor (re-graded from Major) · **Initial:** Apply (reword) · **Final:** Applied with adaptation
- **Reasoning:** Verification showed the overclaim lives **only in `PROJECT_STATUS.md` and commit messages** — the code's sole bit-identity claim (`pair_hmm.rs:592`, "bit-identical to `forward`") is the *correct* P1 claim and was kept. With no `.ssr.psp` baseline or golden test (verified), the ULP-difference is inert. So the fix is a doc reword, not a code change.
- **Adaptation:** reworded the PROJECT_STATUS perf bullet rather than touching code (no code overclaim existed); commit messages are immutable history and left as-is.
- **Files changed:** `PROJECT_STATUS.md` (perf-applied bullet now distinguishes P1 = provably bit-identical from H1/H2 = algebraically equivalent / within-ULP, absorbed by f32 storage; notes no baseline exists).
- **Validation:** n/a (doc).
- **Residual risk:** if a committed golden `.ssr.psp` is ever added, regenerate it after H1/H2 (none exists today).

### M2 — `decoded==0` walk divergence from `refill`
- **Severity:** Minor (re-graded; latent — unreachable for well-formed CRAM, as `read_container` returns 0 only at the EOF marker which the `.crai` never indexes) · **Initial:** Apply via M5 · **Final:** Applied
- **Implementation:** the caching walk now `break`s on an empty decode, matching `refill`'s `Ok(true)` stop. Folded into M5: the shared helper returns `Ok(None)` at EOF, which `refill` maps to its historical stop and the caching path maps to `break`.
- **Files changed:** `src/bam/segment_reader.rs` (`fetch_mapped_reads`).
- **Validation:** `caching_cram_reader_matches_per_call_path` + the 21 CRAM in-module tests pass (clippy/test 0).
- **Residual risk:** None (a decoded-but-record-less container is not a real CRAM shape; documented inline).

### M3 — `ContainerCache` eviction + multi-container untested
- **Severity:** Major · **Initial:** Apply · **Final:** Applied with adaptation
- **Reasoning:** the highest-value gap (FIFO eviction order + insert dedup never executed) is closed with direct unit tests; the multi-container *fetch* test is blocked on a fixture capability.
- **Implementation:** added `container_cache_evicts_oldest_over_capacity` (insert 4 distinct offsets into a cap-3 cache; assert the oldest is evicted, the rest by record-count) and `container_cache_insert_is_idempotent_on_repeat_offset` (repeat-offset insert is dropped; first wins).
- **Adaptation:** the review also asked for a multi-container `fetch_mapped_reads` vs per-call diff. The existing CRAM test fixture (`cram_alignment_file`) builds a **single** container and exposes no knob to force multiple containers/slices, so that test is **deferred** (same multi-container-fixture limitation already tracked as the segment-reader review's Mi7). The unit tests fully cover the eviction/dedup logic that had zero coverage.
- **Files changed:** `src/bam/segment_reader.rs` (tests).
- **Validation:** `cargo test --lib` → both new tests pass (1171 total).
- **Follow-up:** multi-container/eviction-under-load `fetch_mapped_reads` byte-identity test, when a multi-container CRAM fixture knob lands.

### M4 — cache cap not runtime-inspectable
- **Severity:** Nit (re-graded from Major) · **Initial:** Apply (doc) · **Final:** Applied with adaptation
- **Reasoning:** per the discussion, the cap is **not output-affecting** (decode count varies, decoded records do not), so the review's "record it in the `.ssr.psp` header" recommendation was **withdrawn** (recording a non-result knob in result-provenance would mislead).
- **Adaptation:** documented `DEFAULT_MAX_CACHED_CONTAINERS` as an intentionally fixed, non-tunable, non-recorded speed/memory constant — instead of exposing/recording it.
- **Files changed:** `src/bam/segment_reader.rs` (const doc).
- **Validation:** `cargo doc` → 0.
- **Residual risk:** None.

### M5 — `.crai` walk + decode duplicated with `refill`
- **Severity:** Major · **Initial:** Apply · **Final:** Applied
- **Implementation:** extracted a module-level `decode_cram_container(reader, repository, header, path, offset) -> Result<Option<Vec<RecordBuf>>>` that performs the seek + `read_container` + per-slice decode once, returning `Ok(None)` at EOF. Both `CramSegmentReads::refill` (per-call) and `CachingCramReader::decode_container` (decode-once) now call it, so the seek/decode path and its EOF semantics cannot drift. This is what made M2 + Mi3 + Mi6 fall out for free (single `to_io` error closure; PANIC-FREE-commented `expect`). The two now-dead `io_error` methods were removed.
- **Adaptation note:** the helper returns `Option` (not a bare `Vec`) specifically to preserve `refill`'s exact EOF-vs-empty distinction while letting the caching path collapse to a slice — so the refactor is behaviour-preserving on the out-of-scope `refill` path, validated by its 21 existing tests.
- **Review suggestion used verbatim?:** No — the review sketched separate `overlapping_container_offsets` + `decode_container_records` helpers; I unified only the decode (the walk control flow differs structurally between the lazy Iterator and the eager loop, and unifying it would have meant a riskier rewrite of out-of-scope code for no behavioural gain).
- **Files changed:** `src/bam/segment_reader.rs`.
- **Validation:** full suite green — `cargo test --lib` 1171, `--tests` all pass, clippy `--all-targets` 0. The unchanged 21 CRAM tests passing is the byte-identity evidence for the `refill` rewrite.
- **Residual risk:** the `.crai` *walk* (offset selection) is still written twice (Iterator state machine vs `for` loop) but is now behaviourally identical (M2 closed the one divergence); unifying that too is a larger, lower-value refactor left as a potential follow-up.

### Mi1 — stale `window=10`
- **Final:** Applied. Updated `bench_harness.rs` doc (default is `DEFAULT_WINDOW` = 6), the bench prose + the fixed/swept window values (`10`→`6`, sweep `[3,10,15]`→`[3,6,15]`), and the example's default + doc (`10`→`6`). Files: `bench_harness.rs`, `benches/ssr_pileup_perf.rs`, `examples/profile_ssr_pileup.rs`. Validation: clippy `--all-targets` 0.

### Mi3 — PANIC-FREE comment on `expect`
- **Final:** Applied (via M5). `decode_container`'s `self.reader.as_mut().expect("opened by ensure_open")` now carries `// PANIC-FREE: ensure_open populated self.reader immediately above.` Validation: clippy 0.

### Mi4 — single-candidate bit-identity test gap
- **Final:** Applied. Added `score_candidates_single_candidate_matches_forward_bit_for_bit` (window 0 → one rung → `n == lcp`, empty tail loop), asserting `to_bits()` equality vs `forward`. Catches a seam regression in the least-exercised corner. File: `pair_hmm.rs`. Validation: test passes.

### Mi5 — manual tools mislabeled
- **Final:** Applied. Both `ssr_psp_concordance` and `filter_catalog_to_regions` doc comments now lead with "Manual … — **not a CI regression test**". Files: `driver.rs`, `catalog/io.rs`. Validation: doc 0.

### Mi6 — inline `Io` vs helper + misleading comment
- **Final:** Applied (via M5). The shared helper uses a single `to_io` error closure; the inline `AlignmentInputError::Io {…}` duplicates and the misleading "disjoint fields, so the inline error construction is allowed" comment are gone. File: `segment_reader.rs`.

### Mi7 — `#[allow(too_many_arguments)]` unjustified
- **Final:** Applied. Added a hot-path-row-kernel justification comment above the `allow` on `fill_row`. File: `pair_hmm.rs`. Validation: clippy 0.

### Mi8 — `Arc<Vec<RecordBuf>>` → `Arc<[RecordBuf]>`
- **Final:** Applied. `ContainerCache` and `get_or_decode` now use `Arc<[RecordBuf]>` (one fewer indirection on the served slice; `Vec` → `Arc<[_]>` via `.into()`). File: `segment_reader.rs`. Validation: clippy/test 0.

### Mi9 — no compile-time `Send` guard on `WorkerReader`
- **Final:** Applied. Added `assert_send` helper + `worker_reader_is_send` test (`assert_send::<WorkerReader<'static>>()`), mirroring `alignment_file_is_sync`. File: `segment_reader.rs`. Validation: test passes.

### Mi10 — `ContainerCache.cap` → `capacity`
- **Final:** Applied. Renamed the field + `new` param to `capacity` (matches `Reservoir.capacity`). File: `segment_reader.rs`. Validation: clippy 0.

### Doc nits
- **Final:** Applied. Fixed the `ensure_open` doc (was describing `decode_container`), corrected the `CachingCramReader` "`Arc`-shared" wording (repository is clone-shared, filter is `Copy`), and documented `ContainerCache`'s intentional linear scan. File: `segment_reader.rs`.

## 5. Deferred findings to carry forward
- **Mi2** — `bench_harness` `pub mod` not `cfg`-gated. Deferred: it's a deliberate, doc-flagged seam on an internal CLI; `#[cfg(any(test, feature = "bench-internals"))]` gating is a build-config decision, not a defect. Low priority.
- **Mi11** — new bench has no regression/baseline guard. Deferred: criterion benches are not CI-gated in this project; adding a threshold is a design choice.
- **Mi12** — driver batch-fill error lacks positional context. Deferred: `driver.rs` was restructured to `par_chunks` by the parallel session (context shifted); low value, re-assess against the new loop.
- **M3 (partial)** — multi-container `fetch_mapped_reads` byte-identity test. Deferred: needs a multi-container CRAM fixture knob the test harness lacks (shared with segment-reader Mi7). The eviction/dedup half is Applied.
- Cosmetic name nits (`cands`/`out`/`v`) and the harmless `psp_writer_perf` name in the `alloc-mimalloc` feature doc — Deferred (low value).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** No.
- **Baseline saved:** No.
- **Reason:** no `Apply` changed hot-path executable code reachable from `benches/`. The benched `pair_hmm.rs` got only a new test + an `#[allow]` justification comment; all executable edits are in `segment_reader.rs` (the CRAM fetch path, not reachable from any bench). Per the skill, skipped.

## 10. Commands run
- `git status` / `git diff main..HEAD` / `git diff --stat 9c47399`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` (×2)
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo test --tests`
- `./scripts/dev.sh cargo doc --no-deps`

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0 (after fixing 7 `redundant_closure` lints: `|e| to_io(e)` → `&to_io`)
- `cargo test --lib` → 0, `1171 passed; 0 failed; 3 ignored`
- `cargo test --tests` → 0, all suites pass (lib 1171 + main 3 + integration 45)
- `cargo doc --no-deps` → 0, generated index.html (B1 gate restored)

## 12. Notes
- The M5 extraction touched out-of-scope pre-existing code (`CramSegmentReads::refill`); it is byte-identity-gated by `refill`'s 21 existing CRAM tests, all of which pass — the refactor is behaviour-preserving on that path by construction (`Option`-returning helper preserves the EOF-vs-empty distinction).
- M2's "break on empty" closes the only *behavioural* divergence between the two CRAM walks; the remaining duplication is the walk's offset-selection control flow, which is now identical in behaviour and left as an optional follow-up.
- New test count: 1167 → 1171 (`worker_reader_is_send`, `container_cache_evicts_oldest_over_capacity`, `container_cache_insert_is_idempotent_on_repeat_offset`, `score_candidates_single_candidate_matches_forward_bit_for_bit`).
