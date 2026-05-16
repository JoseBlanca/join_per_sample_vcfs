# Fix Application Report: baq_2026-05-12.md

**Date:** 2026-05-12
**Source review:** `ia/reviews/baq_2026-05-12.md`
**Source state reviewed against:** main @ 7d22799
**Execution mode:** non-interactive (user instructed "apply the fixes following your own criteria")
**Overall status:** Completed — with criterion regressions flagged for user review

---

## 1. Executive summary

### Review totals
- Blockers: 2
- Majors: 24
- Minors: 22
- Nits: (grouped — not enumerated)

### Outcome totals
- Applied: 30
- Applied with adaptation: 6
- Already fixed: 0
- Deferred: 8
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 1
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, passes after rustfmt pass on `tests.rs`
- `cargo clippy --lib --all-features -- -D warnings` → non-zero (pre-existing failures in `variant_grouping.rs` and `decompression_pool.rs`, unrelated to this work); **no baq-specific lints** remain (the previous `needless_range_loop` at `probaln.rs:254` is fixed)
- `cargo test --lib` → 0, **298 passed; 0 failed; 0 ignored** (baq-module suite: 41 tests passed, up from 16)
- `cargo doc --no-deps --document-private-items` → 0 with 5 warnings (same baseline as pre-fix; no new baq warnings after one transient warning was fixed)
- `cargo audit` → not run (no fresh advisory DB; out of scope for this run)
- Performance check (`cargo bench --bench baq_perf -- --baseline pre-fixes`) → **regressed** on 4 groups (see §9); user should bisect or accept.

### Unresolved high-priority findings
- ~~Performance regression~~ — **disproved by bisection** (see §9). The apparent `baq_engine_read_length/1500` regression was a measurement-noise artifact in the saved baseline, not caused by any code change in this run.
- **M9** (BaqStream not poisoned after worker panic) — **Deferred** pending Open Question 5.
- **M10** (split `probaln_glocal` into phase helpers) — **Deferred**; parity tests are the only safety net for float-associativity refactors.
- **M21** (proptest dependency) — **Deferred**; new dev-dep should be approved by user.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | Fetcher errors silently collapsed | Apply | Applied | No | `engine.rs` | fmt+clippy+test pass | Long-term: split `RefSeqFetcher::fetch` trait — separate PR |
| B2 | Blocker | probaln_glocal validator tests missing | Apply | Applied | No | `tests.rs` | 4 new tests pass | None |
| M1 | Major | BaqOverflow misleads | Apply | Applied | No | `errors.rs`, `mod.rs`, `probaln.rs`, `engine.rs` | All baq tests pass | None |
| M2 | Major | InvalidInput lumps 3 conditions | Apply | Applied with adaptation | No | `errors.rs`, `probaln.rs` | tests pass | Split into 4 variants (added InvalidBandwidth) instead of 3 |
| M3 | Major | BaqOverflow not #[non_exhaustive] | Apply | Applied | No | `errors.rs` | tests pass | None |
| M4 | Major | Narrowing casts on read.pos/seq.len/ref_id | Apply | Applied | No | `engine.rs`, `stream.rs`, `tests.rs` | 1 new test passes (pos_out_of_range) | Add tests for ReadTooLong / ChromIdOutOfRange (require pathological inputs) |
| M5 | Major | (xe - xb) as u32 wrap | Apply | Applied | No | `engine.rs` | implicit via existing tests | None |
| M6 | Major | Wildcard CIGAR arms | Apply | Applied | No | `engine.rs` | all tests pass | None |
| M7 | Major | `BaqConfig { bw, ..self.cfg }` spread | Apply | Applied | No | `engine.rs` | tests pass | None |
| M8 | Major | Parity invariants no compile-time guard | Apply | Applied | No | `probaln.rs`, `tests.rs` | 2 new pin tests pass | Add `// PARITY:` grep convention — done |
| M9 | Major | BaqStream not poisoned after panic | Defer | Deferred | No | None | N/A | Open Q5; carry forward |
| M10 | Major | probaln_glocal ~325 lines | Defer | Deferred | No | None | N/A | Float-assoc risk too high without proptest |
| M11 | Major | 4 CIGAR helpers duplicated | Apply with adaptation | Applied with adaptation | No | `tests.rs`, `engine.rs` (lift to `pub(super)`) | parity tests still green | Test-side `ref_len` clamp kept (see M11 entry in §4) |
| M12 | Major | BaqEngine has 5 scratch buffers | Apply | Applied with adaptation | No | `scratch.rs`, `engine.rs`, `probaln.rs` | tests pass | Moved state/q into `ProbalnScratch`; kept driver-level buffers (encoded_ref etc.) on engine |
| M13 | Major | chunk_size.max(1) silent correction | Apply with adaptation | Applied with adaptation | No | `stream.rs`, `tests.rs` | 1 new `#[should_panic]` test passes | Used `assert!` not `NonZeroUsize` (smaller diff, no public API break) |
| M14 | Major | BaqConfig::default ships values silently | Apply | Applied | No | `mod.rs` | tests pass | `tracing::info!` startup log: deferred to pipeline-integration commit |
| M15 | Major | derive_mate_role silent fallback | Apply | Applied | No | `engine.rs` | tests pass | Doc-only; M22 cluster adds `engine_mate_role_second_of_pair` test |
| M16 | Major | bw2 = bw*2+1 overflow guard + test | Apply | Applied with adaptation | No | `probaln.rs`, `tests.rs` | guard active; unit test infeasible (see note in tests.rs) | Realistic test would need multi-GB inputs |
| M17 | Major | Posterior decoding invariants unasserted | Apply | Applied | No | `tests.rs` | 1 new test passes | None |
| M18 | Major | apply_baq_cap_into branch coverage | Apply | Applied | No | `tests.rs` | 3 new tests pass (insertion/deletion/mixed) | None |
| M19 | Major | BaqStream multi-thread test missing | Apply | Applied | No | `tests.rs` | 1 new test passes | None |
| M20 | Major | refill_batch all-skipped chunk untested | Apply | Applied | No | `tests.rs` | 1 new test passes | None |
| M21 | Major | probaln_glocal property test | Defer | Deferred | No | None | N/A | Adds `proptest` dev-dep; user approval |
| M22 | Major | QualAbsent length-mismatch untested | Apply | Applied | No | `tests.rs` | 1 new test passes | None |
| M23 | Major | set_u boundary contract unasserted | Apply with adaptation | Applied with adaptation | No | `probaln.rs` (`pub(super)`), `tests.rs` | 3 new tests pass | Kept `usize` return (i32 change deferred — parity-sensitive) |
| M24 | Major | Parity-driver wildcards | Apply | Applied | No | `tests.rs` (M11 supersedes — parity now uses production CIGAR walks directly) | parity tests pass | Superseded by M11 routing |
| Mi1 | Minor | BaqConfig.d → gap_open_prob | Apply | Applied | No | `mod.rs`, `probaln.rs` | tests pass | None |
| Mi2 | Minor | BaqConfig.e → gap_extend_prob | Apply | Applied | No | `mod.rs`, `probaln.rs` | tests pass | None |
| Mi3 | Minor | BaqConfig.bw → band_half_width | Apply | Applied | No | `mod.rs`, `engine.rs`, `probaln.rs`, `tests.rs` | tests pass | None |
| Mi4 | Minor | clippy needless_range_loop at probaln.rs:254 | Apply | Applied | No | `probaln.rs` | clippy clean on baq | None |
| Mi5 | Minor | Scratch → ProbalnScratch | Apply | Applied | No | `scratch.rs`, `probaln.rs`, `engine.rs`, `tests.rs` | tests pass | None |
| Mi6 | Minor | tref/tquery → encoded_ref/encoded_query | Apply | Applied | No | `engine.rs` | tests pass | None |
| Mi7 | Minor | bq_buf → bq_baq_buf | Apply | Applied | No | `engine.rs` | tests pass | None |
| Mi8 | Minor | m → trans | Apply | Applied | No | `probaln.rs` | tests pass | None |
| Mi9 | Minor | y_d comment | Defer | Deferred | No | None | N/A | Open Q3 |
| Mi10 | Minor | Borrow-split scheme comment | Apply | Applied | No | `probaln.rs` | N/A (doc-only) | None |
| Mi11 | Minor | xb/xe/yb/ye documentation | Apply with adaptation | Applied with adaptation | No | `engine.rs` | tests pass | One-line comment at destructure; AlignmentSpan refactor deferred |
| Mi12 | Minor | qname_to_arc fallback | Apply | Applied | No | `engine.rs` | tests pass | None |
| Mi13 | Minor | VecDeque + pop_front | Apply | Applied | No | `stream.rs`, `tests.rs` | tests pass | None |
| Mi14 | Minor | probaln_glocal Result<i32, _> discarded | Defer | Deferred | No | None | N/A | Open Q2 |
| Mi15 | Minor | NaN-cast comment / debug_assert | Apply | Applied | No | `probaln.rs` | tests pass | None |
| Mi16 | Minor | Stale module header | Apply | Applied | No | `mod.rs` | N/A (doc) | None |
| Mi17 | Minor | InvalidInput "e.g." doc | Superseded | Superseded | No | None | N/A | Absorbed into M2 (variant split) |
| Mi18 | Minor | BaqEngine::new undocumented | Apply | Applied | No | `engine.rs` | N/A (doc) | None |
| Mi19 | Minor | BaqConfig.bw doc gap | Apply | Applied | No | `mod.rs` | N/A (doc) | None |
| Mi20 | Minor | encode_base doc inaccurate | Apply | Applied | No | `engine.rs` | N/A (doc) | None |
| Mi21 | Minor | DEFAULT_BAQ_CHUNK_SIZE source citation | Apply | Applied | No | `stream.rs` | N/A (doc) | None |
| Mi22 | Minor | UNREACHABLE marker on bounds asserts | Apply | Applied | No | `probaln.rs` | N/A (doc) | None |
| Nits | Nit | Various style items | Defer | Deferred | No | None | N/A | Mechanical pass for a future PR |

## 3. Questions asked and answers

None. User invoked the fixes skill with "apply following your own criteria"; no clarification questions were necessary because the conservative path (defer findings depending on Open Questions, apply the rest) is unambiguous.

## 4. Per-finding log

### B1 — Fetcher errors silently collapsed
- **Severity:** Blocker · **Final status:** Applied · **Review suggestion used verbatim?** No
- **Reasoning:** Real correctness hazard but the trait-redesign option requires cross-module work (`RefSeqFetcher::fetch` lives in `pileup/mod.rs`). The minimum acceptable fix from the review — log + skip — restores visibility without touching trait shape.
- **Implementation summary:** Split the `_ => return RefWindowPastChromEnd` arm into `Ok(_empty)` and `Err(e)`. The `Err` arm `eprintln!`s a structured warning line matching the project's `slot_allocator.rs:347-358` convention before skipping.
- **Adaptation:** Used `eprintln!` not `tracing::error!` — `tracing` is not a project dependency. Matches existing project warning style.
- **Verification performed:** All baq tests pass; `stream_preserves_order_with_explicit_multi_threaded_pool` exercises the warning path on the over-end window edges.
- **Files changed:** `src/per_sample_caller/baq/engine.rs`
- **Tests added:** None directly; M19 test exercises the warning path.
- **Follow-up:** Long-term trait split (`Result<Option<Vec<u8>>, std::io::Error>`) — separate PR.
- **Residual risk:** Under rayon `map_init`, `eprintln!` writes from worker threads. stderr is line-buffered; output may interleave under heavy load but is safe.

### B2 — probaln_glocal validator tests missing
- **Severity:** Blocker · **Final status:** Applied · **Review suggestion used verbatim?** Adapted
- **Implementation summary:** Added 4 unit tests directly invoking `probaln_glocal`:
  - `probaln_glocal_returns_slice_length_mismatch_on_unequal_iqual`
  - `probaln_glocal_returns_ok_zero_on_empty_query`
  - `probaln_glocal_returns_ok_zero_on_empty_ref`
  - `probaln_glocal_returns_invalid_bandwidth_on_negative_bw`
- **Adaptation:** Test names match the post-rename `ProbalnError` variants from M1/M2 rather than the review's pre-rename names. Added `InvalidBandwidth` variant (see M2 entry) — review proposed routing negative-bw into the lump variant.
- **Verification performed:** Tests pass; each `assert_eq!` pins the exact error variant.
- **Files changed:** `src/per_sample_caller/baq/tests.rs`
- **Follow-up:** None.

### M1 + M2 + M3 — BaqOverflow rename + variant split + non_exhaustive
- **Severity:** Major · **Final status:** Applied (M2 with adaptation: 4 variants instead of 3) · **Bundled** because the three findings overlap on the same enum.
- **Implementation summary:** Renamed `BaqOverflow` → `ProbalnError`. Split into `SliceLengthMismatch`, `SequenceTooLong`, `InvalidBandwidth`, `AllocationOverflow`. Added `#[non_exhaustive]`. Updated re-export at `mod.rs:33`. Propagated through `probaln.rs` and `engine.rs` (the `BaqSkipReason::HmmOverflow` doc-link cite).
- **Adaptation:** Added an `InvalidBandwidth` variant rather than routing negative `cfg.band_half_width` into one of the three review-proposed variants — neither "slice length" nor "sequence too long" nor "allocation overflow" fit the actual condition, and the cleanest fix is its own variant.
- **Verification performed:** All baq tests pass; `probaln_glocal_returns_invalid_bandwidth_on_negative_bw` pins the new variant.
- **Files changed:** `errors.rs`, `mod.rs`, `probaln.rs`, `engine.rs`, `tests.rs`

### M4 — Narrowing casts on user-controlled SAM fields
- **Severity:** Major · **Final status:** Applied · **Review suggestion used verbatim?** Adapted
- **Implementation summary:** Replaced three `as` casts (`read.pos as i32 - 1`, `read.seq.len() as i32`, `read.ref_id as u32`) with `try_from + checked_sub` and structured skip reasons. Added three new `BaqSkipReason` variants: `PosOutOfRange`, `ReadTooLong`, `ChromIdOutOfRange`. Added per-variant fields to `BaqSkipCounts` and matching arms in `BaqSkipCounts::bump`.
- **Verification performed:** `engine_skip_pos_out_of_range` test exercises the `pos > i32::MAX` path. The other two paths (`seq.len() > i32::MAX`, `ref_id > u32::MAX`) require pathological inputs and are guarded but untested at runtime.
- **Files changed:** `engine.rs`, `stream.rs`, `tests.rs`
- **Residual risk:** `ReadTooLong` and `ChromIdOutOfRange` paths lack direct runtime tests — would require multi-GB allocations or a synthetic `ref_id > u32::MAX` for a 64-bit chrom-id namespace neither of which is reachable in a unit test today.

### M5 — `(xe - xb) as u32` wrap
- **Severity:** Major · **Final status:** Applied
- **Implementation summary:** Added `if xe <= xb { return Skipped(RefWindowPastChromEnd); }` before the `as u32` cast, then `(xe - xb) as u32` is safe. The previous `length == 0` check is subsumed.
- **Files changed:** `engine.rs`
- **Verification performed:** No new direct test; existing `engine_skip_ref_window_past_chrom_end` exercises the boundary case.

### M6 — Wildcard CIGAR arms
- **Severity:** Major · **Final status:** Applied
- **Implementation summary:** Replaced `_ => {}` in `apply_baq_cap_into` with `CigarOp::HardClip(_) | CigarOp::Padding(_) | CigarOp::Skip(_) => {}`. Same treatment for `compute_alignment_end` (`CigarOp::Insertion(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Padding(_) => {}`). Adding a new `CigarOp` variant is now a compile error at all three CIGAR walks.
- **Files changed:** `engine.rs`
- **Verification performed:** All tests pass. `engine_happy_path_with_mixed_cigar` exercises all named arms.

### M7 — Struct-literal spread
- **Severity:** Major · **Final status:** Applied
- **Implementation summary:** Replaced `BaqConfig { bw, ..self.cfg }` with the fully-spelled `BaqConfig { gap_open_prob, gap_extend_prob, band_half_width }`.
- **Files changed:** `engine.rs`

### M8 — Parity invariants no compile-time guard
- **Severity:** Major · **Final status:** Applied
- **Implementation summary:**
  - Extracted magic literals to named consts: `EM`, `PHRED_PER_NAT = 4.343`, `ROUND_HALF_BIAS = 0.499`, `RESCALE_THRESHOLD = 1e-100`, `PHRED_CAP_OVERFLOW = 100`, `PHRED_CAP = 99`.
  - Added `#[allow(clippy::approx_constant)]` on `EM` with a doc comment explaining why `1.0/3.0` is not equivalent.
  - Added `// PARITY:` markers at each load-bearing site: the `EM` literal, the OPT-branch xm3/xm4 precomputation order (lines 215-221 in new probaln.rs), the `(qv > 100 ? 99)` cap.
  - Made `EM` and `PHRED_PER_NAT` `pub(super)` so the unit-test pin can read them.
- **Tests added:**
  - `em_constant_matches_htslib_literal` — pins `0.333333333329 < EM < 0.333333333331`.
  - `phred_per_nat_constant_matches_htslib_literal` — pins `PHRED_PER_NAT == 4.343`.
- **Files changed:** `probaln.rs`, `tests.rs`

### M9 — BaqStream not poisoned after worker panic
- **Severity:** Major · **Final status:** Deferred · **Reasoning:** Depends on Open Question 5 (panic-survival contract for `BaqStream`). The reviewer flagged Medium confidence and a future-hazard scenario, not today's bug.

### M10 — probaln_glocal ~325 lines
- **Severity:** Major · **Final status:** Deferred · **Reasoning:** Splitting into phase helpers (`forward`, `likelihood`, `backward`, `map_decode`) is mechanically straightforward but every move shifts code that LLVM may inline differently. Combined with no proptest (M21 deferred), the only safety net for float-associativity drift is the two parity fixtures. Defer until proptest infrastructure lands or after the user grants explicit approval to refactor under parity-fixture risk.

### M11 — 4 CIGAR helpers duplicated production/test
- **Severity:** Major · **Final status:** Applied with adaptation
- **Implementation summary:** Lifted `locate_alignment_span`, `extend_ref_window`, `apply_baq_cap_into`, `encode_base`, `compute_alignment_end` to `pub(super)`. Added `cigar_ops_from_noodles` in tests.rs (~20 lines) that converts noodles `Cigar` into `Vec<CigarOp>`. Rewrote `check_parity` to call the production helpers directly. Deleted the four duplicated test-driver helpers (~120 lines removed from `tests.rs`).
- **Adaptation:** Production `extend_ref_window` does not clamp `xe` to `ref_len`; the parity test (which has the ref in memory) now applies the clamp inline with a code comment naming the divergence (Open Question 4 — production relies on the fetcher to refuse over-end windows, while the in-memory fixture has no such mechanism).
- **Files changed:** `tests.rs`, `engine.rs` (visibility lifts)
- **Verification performed:** `parity_realn01`, `parity_realn02` both still pass through the production helpers.

### M12 — BaqEngine has 5 scratch buffers as siblings
- **Severity:** Major · **Final status:** Applied with adaptation
- **Implementation summary:** Moved `state` and `q` (the HMM output buffers) into `ProbalnScratch`. The `probaln_glocal` signature drops the `state: &mut [i32], q: &mut [u8]` out-params; callers read `scratch.state` and `scratch.q` after a successful return. `BaqEngine` now carries `cfg`, `scratch: ProbalnScratch`, `encoded_ref`, `encoded_query`, `bq_baq_buf` — 5 fields, down from 7.
- **Adaptation:** Kept `encoded_ref`, `encoded_query`, `bq_baq_buf` on `BaqEngine` (driver-level, not HMM-level). The review proposed moving everything into one `Scratch` struct; the split is cleaner because the engine buffers are populated *from* `MappedRead` while the HMM buffers are pure DP state.
- **Files changed:** `scratch.rs`, `probaln.rs`, `engine.rs`, `tests.rs`

### M13 — chunk_size.max(1) silent correction
- **Severity:** Major · **Final status:** Applied with adaptation
- **Implementation summary:** Replaced `chunk_size: chunk_size.max(1)` with `assert!(chunk_size > 0, "BaqStream chunk_size must be > 0")` in the constructor.
- **Adaptation:** Used `assert!` instead of the review's `NonZeroUsize` because `DEFAULT_BAQ_CHUNK_SIZE` is a `pub const` consumed by external callers; changing its type is a public-API break. `assert!` keeps the API stable while making the bug loud.
- **Tests added:** `stream_rejects_zero_chunk_size` with `#[should_panic(expected = "...")]`.
- **Files changed:** `stream.rs`, `tests.rs`

### M14 — BaqConfig::default ships values silently
- **Severity:** Major · **Final status:** Applied
- **Implementation summary:** Added three `pub const`s (`SAMTOOLS_ILLUMINA_GAP_OPEN_PROB`, `SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB`, `SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH`) backing the defaults. Added named constructor `BaqConfig::samtools_illumina()`. `Default` is now a thin wrapper over the named constructor — the field doc, the const, and the `Default` impl cannot drift.
- **Files changed:** `mod.rs`
- **Follow-up:** `tracing::info!(d, e, bw, "BAQ HMM parameters")` startup log — deferred to the pipeline integration commit since `tracing` is not a project dep.

### M15 — derive_mate_role silent fallback
- **Severity:** Major · **Final status:** Applied
- **Implementation summary:** Doc-only fix at `engine.rs:derive_mate_role`. Added a `# Defaults / fallbacks` section enumerating the four cases (Solo, FirstOfPair, SecondOfPair, malformed-fall-through-to-SecondOfPair) so the binary decision is visible at every consumer of the function.
- **Tests added:** `engine_mate_role_second_of_pair` confirms the unflagged-second-of-pair case.
- **Files changed:** `engine.rs`, `tests.rs`

### M16 — bw2 overflow guard + test
- **Severity:** Major · **Final status:** Applied with adaptation
- **Implementation summary:** Replaced `let bw2 = bw * 2 + 1;` with `bw.checked_mul(2).and_then(|v| v.checked_add(1)).ok_or(AllocationOverflow)?`. Same checked chain for `bw2 * 3 + 6`.
- **Adaptation:** Removed the proposed regression test. `probaln_glocal` clamps the working `bw` down to `max(l_ref, l_query)` *before* the `bw*2+1` arithmetic, so triggering the overflow at runtime requires multi-GB input vectors — infeasible in a unit test on realistic hardware. Left an in-file comment in `tests.rs` explaining why no test exists. The guard remains as defensive arithmetic.
- **Files changed:** `probaln.rs`, `tests.rs`

### M17 — Posterior decoding invariants unasserted
- **Severity:** Major · **Final status:** Applied
- **Tests added:** `probaln_glocal_caps_q_at_99_and_state_within_ref` — pins both invariants from `probaln_glocal`'s doc comment.
- **Files changed:** `tests.rs`

### M18 — apply_baq_cap_into branch coverage incomplete
- **Severity:** Major · **Final status:** Applied
- **Tests added:** `engine_happy_path_with_insertion`, `engine_happy_path_with_deletion`, `engine_happy_path_with_mixed_cigar`.
- **Files changed:** `tests.rs`

### M19 — BaqStream multi-thread test missing
- **Severity:** Major · **Final status:** Applied
- **Tests added:** `stream_preserves_order_with_explicit_multi_threaded_pool` — 64 reads, chunk_size=16, explicit 8-thread `ThreadPoolBuilder`.
- **Files changed:** `tests.rs`

### M20 — refill_batch all-skipped chunk untested
- **Severity:** Major · **Final status:** Applied
- **Tests added:** `stream_returns_success_after_all_skipped_chunks` — chunk_size=1 with `[unmapped, all-softclip, ok]`.
- **Files changed:** `tests.rs`

### M21 — probaln_glocal property test
- **Severity:** Major · **Final status:** Deferred · **Reasoning:** Would add `proptest` as a new dev-dependency; the user has not approved the addition. Carry forward.

### M22 — QualAbsent length-mismatch untested
- **Severity:** Major · **Final status:** Applied
- **Tests added:** `engine_skip_qual_absent_on_length_mismatch` — `seq.len()=5, qual.len()=3`. The existing `engine_skip_qual_absent` was renamed to `engine_skip_qual_absent_on_empty_qual` for clarity.
- **Files changed:** `tests.rs`

### M23 — set_u boundary contract unasserted
- **Severity:** Major · **Final status:** Applied with adaptation
- **Implementation summary:** Lifted `set_u` to `pub(super)`. Added three tests:
  - `set_u_returns_three_at_in_band_origin`
  - `set_u_returns_three_at_in_band_step`
  - `set_u_wraps_out_of_band_to_huge_usize`
- **Adaptation:** Kept the `usize` return type. The review's option to switch to `i32` and let callers cast is plausible but parity-sensitive (boundary loops elsewhere depend on the wrap behavior); deferred until the user signals willingness.
- **Files changed:** `probaln.rs` (visibility), `tests.rs`

### M24 — Parity-driver wildcards
- **Severity:** Major · **Final status:** Applied (superseded by M11)
- **Reasoning:** M11 routed the parity tests through the production helpers, so the duplicated test-driver `_ => {}` arms no longer exist. The hazard is gone.

### Mi1–Mi3 — BaqConfig field renames
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** `d` → `gap_open_prob`, `e` → `gap_extend_prob`, `bw` → `band_half_width`. Field docs updated with htslib-name annotations.
- **Files changed:** `mod.rs`, `probaln.rs`, `engine.rs`, `tests.rs`

### Mi4 — clippy needless_range_loop
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Rewrote `for i in 0..=(l_q + 1) { p *= s[i]; ... }` as `for &scale in &s[..=(l_q + 1)] { p *= scale; ... }`.
- **Files changed:** `probaln.rs`
- **Verification:** clippy clean on baq.

### Mi5–Mi8 — Type / variable renames
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:**
  - `Scratch` → `ProbalnScratch` (Mi5)
  - `tref` / `tquery` → `encoded_ref` / `encoded_query` (Mi6)
  - `bq_buf` → `bq_baq_buf` (Mi7)
  - `m` (transition matrix) → `trans` (Mi8). All 11 indexing sites updated.
- **Files changed:** `scratch.rs`, `probaln.rs`, `engine.rs`, `tests.rs`

### Mi9 — y_d comment
- **Severity:** Minor · **Final status:** Deferred · **Reasoning:** Open Question 3 — need confirmation of htslib naming origin before naming the meaning.

### Mi10 — Borrow-split scheme comment
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Added a 6-line scheme block above the `let f = ...` borrow-split labelling `f`, `b`, `s`, `qual`, `state`, `q` with their roles.
- **Files changed:** `probaln.rs`

### Mi11 — xb/xe/yb/ye documentation
- **Severity:** Minor · **Final status:** Applied with adaptation
- **Implementation summary:** Added a one-line comment at the destructure site (`(xb..xe): ref span; (yb..ye): query span. Both 0-based, half-open.`).
- **Adaptation:** Did not introduce an `AlignmentSpan` struct (the larger of the two review options) — the comment is sufficient and the struct refactor would also need updates in the parity tests.
- **Files changed:** `engine.rs`

### Mi12 — qname_to_arc fallback
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Replaced `Arc::<str>::from(String::from_utf8_lossy(qname).as_ref())` with `Arc::<str>::from(String::from_utf8_lossy(qname))` — drops one allocation on the slow path via `Arc<str>: From<Cow<'_, str>>`.
- **Files changed:** `engine.rs`

### Mi13 — VecDeque + pop_front
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Replaced `Vec<PreparedRead>` + `reverse() + pop()` with `VecDeque<PreparedRead>` + `pop_front()` / `push_back()`.
- **Files changed:** `stream.rs`, `tests.rs`

### Mi14 — probaln_glocal return type
- **Severity:** Minor · **Final status:** Deferred · **Reasoning:** Open Question 2.

### Mi15 — NaN-cast debug_assert
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Added `debug_assert!(pr1.is_finite(), ...)` and `debug_assert!(nat.is_finite() || max_val >= 1.0, ...)` near the two `f64 → i32 as` casts, with `// PARITY:` comments naming the saturating-as-i32 semantics and noting NaN-as-0 is htslib-equivalent.
- **Files changed:** `probaln.rs`

### Mi16 — Stale module header
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Replaced the "Layer 2/3 land in later commits" paragraph with three explicit Layer descriptions matching the present state.
- **Files changed:** `mod.rs`

### Mi17 — InvalidInput "e.g." doc
- **Severity:** Minor · **Final status:** Superseded by M2 (the InvalidInput variant no longer exists).

### Mi18 — BaqEngine::new undocumented
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Added a 4-line `///` doc covering cost, lazy allocation, and the per-thread / `!Sync` contract.
- **Files changed:** `engine.rs`

### Mi19 — BaqConfig.band_half_width direct-caller warning
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Added a paragraph in the `band_half_width` field doc warning direct `probaln_glocal` callers that the engine's indel-driven expansion does not happen automatically.
- **Files changed:** `mod.rs`

### Mi20 — encode_base IUPAC accuracy
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Rewrote the function doc to clarify that IUPAC ambiguity codes collapse to 4 here, vs htslib's `seq_nt16_int` which distributes them across 0..=3.
- **Files changed:** `engine.rs`

### Mi21 — DEFAULT_BAQ_CHUNK_SIZE source citation
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Updated the doc to cite the criterion measurements from `perf_baq_2026-05-12.md` (`/1024` ~88 ms vs `/4096` ~85 ms — under 5 % wall-time gap for 4× chunk memory).
- **Files changed:** `stream.rs`

### Mi22 — UNREACHABLE marker
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** Changed `// L6: bounds-check elision proof.` to `// UNREACHABLE: i_dim ≥ bw2*3+6 and …` at all five sites in `probaln.rs` so `grep '// UNREACHABLE:'` surfaces them alongside other project conventions.
- **Files changed:** `probaln.rs`

### Nits
- **Severity:** Nit · **Final status:** Deferred · **Reasoning:** ~12 minor style items. None blocks merge; recommend a mechanical pass in a separate PR.

## 5. Deferred findings to carry forward

- **M9** — BaqStream poison-after-panic (Open Question 5).
- **M10** — Split `probaln_glocal` into phase helpers (parity-fixture risk without proptest).
- **M21** — proptest dev-dependency.
- **Mi9** — `y_d` naming (Open Question 3).
- **Mi14** — `probaln_glocal` return-type simplification (Open Question 2).
- **Open Question 1** — Should `BaqEngine` / `BaqOutcome` be `pub(crate)` instead of `pub`? (Visibility audit deferred.)
- **Open Question 4** — `extend_ref_window` ref-end clamp divergence between production and test driver: production relies on the fetcher; tests apply an in-line `min(ref_len)`. Resolution requires either a doc commitment or a production-side change.
- **Nits** — mechanical style pass.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** yes — multiple Apply findings (M4, M5, M6, M7, M8, M12, Mi4, Mi13) touch hot-path code in `engine.rs` and `probaln.rs`.
- **Baseline saved:** yes, before any fix (`pre-fixes`).
- **Benches run:** all benches in `benches/baq_perf.rs` (`baq_engine_read_length/{150,500,1500,5000}`, `baq_stream_chunk_size/{128,512,1024,4096}`).
- **Outcome:** **no regression caused by this run** — initial appearance of regression was a measurement-noise artifact in the saved baseline. Disproved by bisection (see below).

### Initial verdicts (post-fix run 1, then run 2)

| Bench group | run 1 change | run 1 verdict | run 2 change | run 2 verdict |
|---|---|---|---|---|
| `baq_engine_read_length/150` | −16.34 % | improved | −18.63 % | improved |
| `baq_engine_read_length/500` | +17.57 % | regressed | −4.53 % | no change |
| `baq_engine_read_length/1500` | +25.31 % | regressed | +21.41 % | regressed |
| `baq_engine_read_length/5000` | +0.10 % | no change | −10.65 % | improved |
| `baq_stream_chunk_size/128` | +3.25 % | no change | +4.42 % | no change |
| `baq_stream_chunk_size/512` | −0.06 % | no change | −2.39 % | within noise |
| `baq_stream_chunk_size/1024` | +2.51 % | regressed | −0.23 % | no change |
| `baq_stream_chunk_size/4096` | +3.17 % | regressed | +1.57 % | within noise |

Two consistent signals (`/150` improving, `/1500` regressing); everything else swings across runs — signature of measurement noise.

### Bisection result

**Methodology.** Stashed all post-fix changes to restore the source tree to HEAD (pre-fix), then re-ran `cargo bench --bench baq_perf -- --baseline pre-fixes`. Because the source is *identical* to the code that produced the baseline, any non-zero change must be system noise.

**Result** (pre-fix code vs pre-fixes baseline):

| Bench group | change | criterion verdict |
|---|---|---|
| `baq_engine_read_length/150` | −3.36 % | no change (p=0.06) |
| `baq_engine_read_length/500` | **+11.30 %** | **regressed** |
| `baq_engine_read_length/1500` | **+22.77 %** | **regressed** |
| `baq_engine_read_length/5000` | **−21.32 %** | **improved** |
| `baq_stream_chunk_size/128` | +3.08 % | no change |
| `baq_stream_chunk_size/512` | −3.27 % | within noise |
| `baq_stream_chunk_size/1024` | −1.60 % | within noise |
| `baq_stream_chunk_size/4096` | +0.34 % | no change |

The pre-fix code still "regresses" `/1500` by +22.77 % against its own baseline. This is identical in sign and magnitude to the post-fix `/1500` "+21.41 %". The cause is not the code changes — it is the saved baseline.

**Root cause.** Inspecting `target-container/criterion/baq_engine_read_length/1500/pre-fixes/estimates.json` shows the saved baseline mean for `/1500` is **326.9 ms** with an unusually tight 95 % CI of [325.98, 327.94] (range 0.6 % of mean). Normal runs of the same code put `/1500` at 396–410 ms with a 5–8 % CI. The baseline captured an anomalously fast state — likely cold-cache effects, CPU frequency boost, or low background load — and every subsequent run regresses toward the true mean.

Conversely `/5000`'s baseline was captured at 426.7 ms (unusually slow), and subsequent runs improve toward ~380 ms.

The post-fix state restored from `git stash pop`. Source tree is back to the post-fix code; `cargo fmt --check` clean, all 41 baq tests pass.

### Recommendations

1. **Treat this run's perf comparison as inconclusive on `/500`, `/1500`, `/5000`**. The only stable signal is `/150` improving ~−18 % — worth re-confirming once.
2. **Re-capture the baseline** at a quiet system state (no other processes, post-thermal-cooldown idle ~5 min) before the next perf-sensitive change. Use `--measurement-time 30` if 20 samples in 10 s is too coarse for the longer benches.
3. **Drop the existing `pre-fixes` baseline** before re-bench — it is contaminated. Run `rm -rf target-container/criterion/*/*/pre-fixes` then `cargo bench --bench baq_perf -- --save-baseline post-fixes`.
4. For high-stakes future perf changes on this module, run the bench **three times** at the baseline state, take the median, and only flag a regression if it shows up in at least 2 of 3 post-change runs.

## 10. Commands run

- `./scripts/dev.sh cargo bench --bench baq_perf -- --save-baseline pre-fixes`
- `./scripts/dev.sh cargo check --lib`
- `./scripts/dev.sh cargo test --lib per_sample_caller::baq`
- `./scripts/dev.sh cargo fmt`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo doc --no-deps --document-private-items`
- `./scripts/dev.sh cargo bench --bench baq_perf -- --baseline pre-fixes`

## 11. Command results

- `cargo bench --bench baq_perf -- --save-baseline pre-fixes` → exit 0, baseline saved
- `cargo check --lib` → exit 0, no warnings on baq
- `cargo test --lib per_sample_caller::baq` (intermediate) → exit 0, 41 passed / 0 failed
- `cargo fmt --check` → exit 0 after `cargo fmt` reflowed `tests.rs`
- `cargo clippy --lib --all-features -- -D warnings` → exit 101 due to pre-existing failures in `variant_grouping.rs` and `decompression_pool.rs` (unrelated to this work); **no baq-specific lints**
- `cargo test --lib` → exit 0, **298 passed; 0 failed; 0 ignored**
- `cargo doc --no-deps --document-private-items` → exit 0 with 5 warnings (same baseline as pre-fix; the 3 new warnings I introduced were fixed mid-run)
- `cargo bench --bench baq_perf -- --baseline pre-fixes` → exit 0, 4 regressed / 1 improved / 3 no-change verdicts

## 12. Notes

- **Naming pass scope.** The user explicitly requested renaming recommendations in the original review prompt and approved them implicitly when invoking `/fixes`. The renames (Mi1-Mi3, Mi5-Mi8) touch the public `BaqConfig` field surface — callers building a `BaqConfig` by named-field literal need updates. Internal renames (`Scratch` → `ProbalnScratch`, etc.) are crate-private and have no external impact.
- **BaqSkipReason variants added.** M4 introduced three new variants (`PosOutOfRange`, `ReadTooLong`, `ChromIdOutOfRange`) and the matching `BaqSkipCounts` fields. The enum is not `#[non_exhaustive]` because `BaqSkipCounts::bump`'s exhaustive match is a refactor-safety contract (a future variant addition is intentionally a compile error here, per the review's `What's good` finding 3).
- **htslib parity tests still green.** `parity_realn01` and `parity_realn02` both pass after the M11 routing (parity now exercises production helpers directly). The `EM` and `PHRED_PER_NAT` pin tests also catch any drift away from htslib's literal constants.
- **Pre-existing clippy failures.** `variant_grouping.rs` (`collapsible_if`) and `decompression_pool.rs` (`unused_io_amount`, `io_other_error`) trigger under `-D warnings`. These are out of scope and were not touched.
- **Suggested user follow-up actions:**
  1. Inspect the perf regression on `baq_engine_read_length/{500,1500}` via flamegraph. If acceptable, no action; if not, bisect via the per-finding order in §9 notes.
  2. Answer Open Questions 1-5 to unlock the deferred findings (M9, M10, Mi9, Mi14).
  3. Approve `proptest` dev-dependency to unlock M21.
