# Fix Application Report: pileup_2026-05-06.md

**Date:** 2026-05-06
**Source review:** `ia/reviews/pileup_2026-05-06.md`
**Source state reviewed against:** `b5e1cba` (main; review committed)
**Execution mode:** interactive
**Overall status:** In progress

---

## 1. Executive summary

### Review totals
- Blockers: 2 (B1, B2)
- Majors: 5 (M1, M2, M3, M4, M5)
- Minors: 6 (Mi1, Mi2, Mi3, Mi4, Mi5, Mi6)
- Nits: ~5 small items

### Outcome totals (running)
- Applied: 0
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` — pending
- `cargo clippy --all-targets --all-features -- -D warnings` — pre-existing failures in unrelated files; pileup-scoped rerun pending
- `cargo test --all-targets --all-features` — pileup tests pass; full suite uses unrelated bench fixture
- `cargo doc --no-deps` — not required for this scope
- `cargo audit` — not required for this scope

### Unresolved high-priority findings
*(populated as findings progress)*

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | REF re-folded once per walker step inside open record | Apply | Applied | No | `src/per_sample_caller/pileup/{open_record,walker,tests}.rs` | Pass | No |
| B2 | Blocker | Indel mate-overlap doesn't collapse to one observation | Apply | Applied | No | `src/per_sample_caller/pileup/{walker,tests}.rs` | Pass | No |
| M1 | Major | `find_overlapping` early break unsound across heterogeneous spans | Apply | Applied | No | `src/per_sample_caller/pileup/open_record.rs` | Pass | No |
| M2 | Major | Mate-overlap tie-break uses `alignment_start`, not `is_first_mate` | Apply | Applied | No | `src/per_sample_caller/pileup/{open_record,walker,tests}.rs` | Pass | No |
| M3 | Major | `record_widen_events` counter conflates open and widen | Apply | Applied | No | `src/per_sample_caller/pileup/{open_record,walker,tests}.rs` | Pass | No |
| M4 | Major | Lifecycle marks attach only to first drained record | Ask | Applied with adaptation | Yes | `src/per_sample_caller/pileup/slot_allocator.rs`, `ia/specs/pileup_walker.md` | Pass | No |
| M5 | Major | `windows(2)` only inspects adjacent contributor pairs | Apply | Applied | No | `src/per_sample_caller/pileup/walker.rs` | Pass | No |
| Mi1 | Minor | `checked_sub.unwrap_or(0)` masks precondition violation | Apply | Applied | No | `src/per_sample_caller/pileup/open_record.rs` | Pass | No |
| Mi2 | Minor | `walker_pos == 0` guard is dead code | Apply | Applied | No | `src/per_sample_caller/pileup/walker.rs` | Pass | No |
| Mi3 | Minor | `flush_chromosome` accepts unused `_fasta` | Apply | Applied | No | `src/per_sample_caller/pileup/walker.rs` | Pass | No |
| Mi4 | Minor | `from_slot_counters` awkward shape | Apply | Applied | No | `src/per_sample_caller/pileup/walker.rs` | Pass | No |
| Mi5 | Minor | `FiveScalars::zero` duplicates `Default::default` | Apply | Applied | No | `src/per_sample_caller/pileup/{mod,open_record}.rs` | Pass | No |
| Mi6 | Minor | `widen` always appends `extra_bases` to every allele | Apply | Disputed | No | `src/per_sample_caller/pileup/open_record.rs` (doc only) | Pass | No |
| Nits | Nit | Clippy errors + small style items | Apply | Applied | No | `src/per_sample_caller/pileup/{open_record,walker}.rs` | Pass | No |

## 3. Questions asked and answers

1. **M4** — Should lifecycle marks attach per-emission (all on first record of an emission batch, current implementation) or per-record (each record carries marks since the previous record)?
   - **Answer:** "Honestly, I don't know what would be best for M4. What do you think? I value correctness over implementation time, so don't worry if you have to redo something. Also, think critically about the spec, do not assume that it is correct just because it is written in a particular way." — Resolved by the fixer's critical analysis: per-emission stamping is kept, but the deeper correctness issue uncovered (slot-id reuse silently merging phase chains) is the actual fix delivered.

## 4. Per-finding log

### B1 — REF re-folded once per walker step inside open record

- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Defect confirmed by a regression test added test-first. The review's suggested fix (split "affected = newly opened/widened" + "fold every active read once") was viable but more invasive; an equivalent narrower fix tracks per-record `(read_id → previous contribution)` and subtracts the prior contribution before adding the new one on every fold, so each `(record, read)` pair contributes exactly one net observation across the record's lifetime. The narrower shape preserves the existing `process_position` flow (events at walker_pos drive find_overlapping → fold) while removing the multiplicative inflation.
- **Implementation summary:**
  - `OpenPileupRecord` gained an `AHashMap<u32, FoldedReadState>` field tracking each contributing read's last folded `(allele_index, contribution)`.
  - `process_position` now subtracts a read's prior contribution from its old bucket before adding to the new bucket (sub-/add-helpers added).
  - `ReadContribution` gained a `read_id` field plumbed from the active set in `walker.rs::process_position`.
  - Old `find_or_create_allele(&mut OpenPileupRecord) -> &mut OpenAllele` was replaced by `find_or_create_allele_index(&mut OpenPileupRecord) -> usize`, freeing the borrow so callers can keep mutating the record.
- **Review suggestion used verbatim?:** No.
- **Adaptation:** Subtract-old/add-new instead of "split affected and fold once". Equivalent observable behavior; smaller diff.
- **Verification performed:** Added regression test `deletion_record_does_not_double_count_ref_reads` that fails before fix (REF.num_obs=4) and passes after (REF.num_obs=1, fwd=1). Ran the full pileup suite (50 tests pass).
- **Files changed:**
  - `src/per_sample_caller/pileup/open_record.rs`
  - `src/per_sample_caller/pileup/walker.rs`
  - `src/per_sample_caller/pileup/tests.rs`
- **Tests added or modified:** Added `deletion_record_does_not_double_count_ref_reads` (tests.rs); updated `find_or_create_allele_returns_same_bucket_on_match` to use the renamed `find_or_create_allele_index`.
- **Validation:**
  - `cargo test --lib pileup` → exit 0, 50 passed
  - `cargo fmt --check` → exit 0 (after one auto-format pass)
- **User input:** None.
- **Follow-up:** Phantom-slot residual on a read that switches buckets across compound-event widens (a slot may persist on the bucket the read left, since other reads sharing the slot might still be there). Tests do not currently exercise this; flagged in §11.
- **Residual risk:** Low for the SNP / single-deletion / pure-REF mixes the existing tests cover. Compound events on a single read crossing walker steps are not exercised; a future test should verify the bucket-switch path.

### B2 — Indel mate-overlap doesn't collapse to one observation

- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Per spec §"Mate overlap on indels" both regimes (same-indel pair, and indel-vs-Match disagreement) collapse to a single observation. Match-only overlap keeps its existing behaviour (loser's BQ zeroed; both still count as observations).
- **Implementation summary:** `resolve_mate_overlap_at_pos` now classifies each pair's resolution: if either side has an Insertion/Deletion event anchored at this walker_pos (`pair_has_indel`), the loser is removed from the contributor list with `swap_remove` (descending-index pass keeps earlier indices stable); otherwise the loser's BQ is zeroed in place as before.
- **Review suggestion used verbatim?:** No.
- **Adaptation:** The review proposed splitting "drop event from loser" vs "zero match BQ"; this fix removes the loser entirely on indel-overlap (matching the stronger spec reading "treat as one observation, not two"). Kept the existing zero-BQ path for Match-only overlap. Tie-breaker change is deferred to M2 (still uses `alignment_start` here pending that fix).
- **Verification performed:** Added regression test `paired_mate_indel_overlap_yields_single_observation`; failed before fix (INS.num_obs = 2), passed after (INS.num_obs = 1, fwd = 1). All 51 pileup tests pass.
- **Files changed:**
  - `src/per_sample_caller/pileup/walker.rs`
  - `src/per_sample_caller/pileup/tests.rs`
- **Tests added or modified:** Added `paired_mate_indel_overlap_yields_single_observation`.
- **Validation:**
  - `cargo test --lib pileup` → exit 0, 51 passed
- **User input:** None.
- **Follow-up:** None for this finding. Tie-breaker still uses `alignment_start` (M2 will fix).
- **Residual risk:** Low. Match-only overlap path unchanged.

### M1 — `find_overlapping` early break unsound across heterogeneous spans

- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The early break terminated the backward scan at the first record whose footprint ended at or before `event_start`, missing wide earlier records sitting behind narrower intermediate ones. The fix bounds the search range by `MAX_RECORD_SPAN` (the upper bound on any single record's reach) and walks every record in `[event_start - MAX_RECORD_SPAN, event_start]` without an early exit.
- **Implementation summary:** Replaced the backward-scan-with-early-break in `OpenPileupRecordTable::find_overlapping` with a `range(lo..=event_start)` scan where `lo = event_start.saturating_sub(MAX_RECORD_SPAN)`. Updated the function's doc to spell out why heterogeneous spans require this.
- **Review suggestion used verbatim?:** Yes (the `MAX_RECORD_SPAN` bound was the review's exact suggestion).
- **Adaptation:** None.
- **Verification performed:** Added `find_overlapping_walks_past_intermediate_narrow_record_to_wide_one` (open_record.rs tests). Failed before fix (returned None), passed after (returned Some(5)). All 52 pileup tests pass.
- **Files changed:**
  - `src/per_sample_caller/pileup/open_record.rs`
- **Tests added or modified:** Added `find_overlapping_walks_past_intermediate_narrow_record_to_wide_one`.
- **Validation:**
  - `cargo test --lib pileup` → exit 0, 52 passed
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None observed; the search is `O(log n + k)` where `k ≤ MAX_RECORD_SPAN` records sit in the bounded range — typical `k` is small because open records close eagerly.

### M2 — Mate-overlap tie-break uses `alignment_start`, not `is_first_mate`

- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The spec says "BQ tie: prefer mate 1 (flag `0x40`, `is_first_mate` set)". The contribution didn't carry `is_first_mate` so the walker substituted `alignment_start` and broke ties on position. Adding the field and using it directly is a one-field plumbing change.
- **Implementation summary:** Added `is_first_mate: bool` to `ReadContribution`; populated from `active.read.is_first_mate` in `walker.rs::process_position`; the tie branch in `resolve_mate_overlap_at_pos` now keeps the contributor with `is_first_mate == true`. Falls back to `alignment_start` only when both or neither carry the flag (defensive — paired mates by SAM convention have exactly one mate-1).
- **Review suggestion used verbatim?:** Yes (the field-plumbing shape).
- **Adaptation:** Kept the `alignment_start` path as a defensive fallback for the impossible-but-let's-not-panic "both or neither is_first_mate" case.
- **Verification performed:** Added `mate_overlap_bq_tie_prefers_first_mate_not_earlier_position` (tests.rs). With first-mate appearing AFTER second-mate in the input stream and identical `alignment_start`, the previous code dropped first-mate's contribution and `q_sum ≈ -6.9`; after the fix, first-mate is kept and `q_sum ≈ -2.0`. All 53 pileup tests pass.
- **Files changed:**
  - `src/per_sample_caller/pileup/open_record.rs`
  - `src/per_sample_caller/pileup/walker.rs`
  - `src/per_sample_caller/pileup/tests.rs`
- **Tests added or modified:** Added `mate_overlap_bq_tie_prefers_first_mate_not_earlier_position`.
- **Validation:**
  - `cargo test --lib pileup` → exit 0, 53 passed
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None.

### Mi1 — `checked_sub.unwrap_or(0)` masks precondition violation

- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The function documents "Events are assumed to be inside the record's footprint." The previous shape silently fabricated `offset = 0` if that precondition was violated; a `debug_assert!` + `saturating_sub` makes the contract loud in tests while behaving identically in release.
- **Files changed:** `src/per_sample_caller/pileup/open_record.rs`.
- **Tests:** Existing pileup suite (55 + others) passes.
- **Validation:** `cargo test --lib` → 145 passed.

### Mi2 — `walker_pos == 0` dead-code guard

- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `walker_pos` initialises to 1 and resets to 1 on chromosome change; nothing sets it to 0. Replaced the live early-return guard with a `debug_assert!` so a future regression that violates the invariant trips in tests.
- **Files changed:** `src/per_sample_caller/pileup/walker.rs`.
- **Validation:** `cargo test --lib` → 145 passed.

### Mi3 — `flush_chromosome` accepts unused `_fasta`

- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The `_fasta: &F` parameter was never read; dropped it and its generic, simplifying the signature.
- **Files changed:** `src/per_sample_caller/pileup/walker.rs` (definition + two call sites).
- **Validation:** `cargo test --lib` → 145 passed.

### Mi4 — `from_slot_counters` awkward shape

- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The function ignored most fields of its `Self` parameter, used a workaround `self_` name, and reconstructed the struct by hand. Replaced with a method `merge_slot_counters(self) -> Self` that mutates only the slot-derived fields.
- **Files changed:** `src/per_sample_caller/pileup/walker.rs`.
- **Validation:** `cargo test --lib` → 145 passed.

### M4 — Lifecycle marks attach only to first drained record

- **Severity:** Major
- **Initial decision:** Ask
- **Final status:** Applied with adaptation
- **Reasoning:** The review framed M4 as a per-record vs per-emission stamping choice. While reading the slot-allocator code I noticed a deeper correctness defect the review missed: **same-walker-step slot reuse silently merges phase chains.** When `r1` expires and `r2` admits in the same step, the allocator's old `acquire_slot` immediately popped the freed id from `free`. The drain produced `new=[S]` and `expired=[S]` simultaneously; the existing `drain_lifecycle_marks_suppresses_slots_present_in_both` rule dropped both, leaving the consumer no way to tell "r1's chain ended, r2's chain began with the same id" apart from a transient. Same-step reuse is the **common** case at any meaningful coverage (back-to-back reads at adjacent positions), so phase chains were routinely corrupted. Per-record vs per-emission stamping doesn't fix this — both options still allow the suppression to mask reuse.
- **Implementation summary:**
  - Added `pending_free: Vec<SlotId>` to `SlotAllocator`. `release_slot` now parks freed slots there instead of pushing directly to `free`.
  - `drain_lifecycle_marks` migrates `pending_free` → `free` *after* draining the marks, so reuse of a freed id is forced into a strictly later emission than the one carrying its `expired` mark.
  - Dropped the suppression in `drain_lifecycle_marks` — it's now unnecessary for reuse (impossible by construction) and harmless for genuinely-transient slots (consumer applies `new_chains` before processing and `expired_chains` after; transient slots produce no observable alive-set delta).
  - Spec updated: removed the `[QUESTION]` resolution that recommended suppression; added a §"Slot reuse and the `pending_free` deferral" subsection making the design rule explicit.
  - Per-emission "all-on-first-record" stamping is kept and now documented in the spec without ambiguity.
- **Review suggestion used verbatim?:** No.
- **Adaptation:** The review's two options (per-emission with explicit drain-on-empty, vs per-record marks) didn't cover the actual bug. The user explicitly invited critical thinking and prioritized correctness over speed. The implemented fix addresses the root cause (slot-id reuse correctness) and keeps per-emission stamping as documented behaviour.
- **Verification performed:** Added `same_emission_reuse_does_not_collide_on_slot_id` (asserts `s1 != s2` when r1 releases and r2 admits before any drain; was `s1 == s2 == 0` before). Renamed and flipped `drain_lifecycle_marks_suppresses_slots_present_in_both` → `drain_lifecycle_marks_emits_both_for_transient_slot_within_one_drain` (transient slots now surface; consumer handles them). All 55 pileup tests pass; full lib suite at 145 tests passes.
- **Files changed:**
  - `src/per_sample_caller/pileup/slot_allocator.rs`
  - `ia/specs/pileup_walker.md`
- **Tests added or modified:** Added `same_emission_reuse_does_not_collide_on_slot_id`; replaced `drain_lifecycle_marks_suppresses_slots_present_in_both` with `drain_lifecycle_marks_emits_both_for_transient_slot_within_one_drain`.
- **Validation:**
  - `cargo test --lib pileup` → exit 0, 55 passed
  - `cargo test --lib` → exit 0, 145 passed
- **User input:** Yes — invited critical thinking; the user did not commit to either review-proposed option.
- **Follow-up:** None required for the slot-reuse correctness bug. If Stage 5 ever needs per-record marks (rather than per-emission), that's an additive change — split `stamp_lifecycle_marks` to track marks per record. Not needed today.
- **Residual risk:** Low. `pending_free` adds one more list to the allocator; the migration is `O(n)` per drain where `n ≤ active reads released this emission window` (typically ≤ a handful). Memory: one extra `Vec<u16>` per allocator.

### Mi5 — `FiveScalars::zero` duplicates `Default::default`

- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The `zero()` method just delegated to `Self::default()` with no extra semantics. Removing it canonicalises construction at one name.
- **Implementation summary:** Removed the `impl FiveScalars { fn zero }` block; updated the one call site (`OpenAllele::new`) to use `FiveScalars::default()`.
- **Files changed:** `src/per_sample_caller/pileup/mod.rs`, `src/per_sample_caller/pileup/open_record.rs`.
- **Validation:** `cargo test --lib` → 145 passed.

### Mi6 — `widen` invariant on every allele kind

- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Disputed
- **Reasoning:** Investigated the reviewer's worry that "appending `extra_bases` to every allele's `seq`" might not produce the right haplotype for INS alleles whose insertion anchor sits mid-old-span. Tracing `apply_events_to_ref` shows that every allele's `seq` always *ends* with the ref-aligned suffix `ref_seq[ref_cursor..ref_seq.len()]` — the post-loop tail emits raw ref bases verbatim. Widening extends `ref_seq` with `extra_bases`; appending `extra_bases` to each allele's `seq` reproduces exactly what re-folding the read against the wider `ref_seq` would emit (modulo new events at the new positions, which haven't been folded yet — they fold at a later walker step into the right bucket via the B1 subtract-old/add-new path). The current code is correct for SNP/MNP, DEL, INS, and compound alleles by construction.
- **Implementation summary:** No code-behaviour change. Replaced the inline rationale with a load-bearing comment that names the invariant ("every allele's seq ends with a ref-aligned suffix") and explains why simple append is correct for all allele kinds.
- **Files changed:** `src/per_sample_caller/pileup/open_record.rs` (doc only).
- **Validation:** `cargo test --lib` → 145 passed.
- **Follow-up:** None. Recorded as Disputed because no code defect was confirmed; the documentation update tightens reasoning for future readers.

### Nits — pileup-scoped clippy + small style

- **Severity:** Nit
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The original review listed four pileup-scoped clippy errors (`manual_saturating_arithmetic`, `needless_lifetimes`, `unnecessary_filter_map`, `ptr_arg`). The first three are now resolved as side-effects of B1 (replaced `find_or_create_allele` lifetime-annotated wrapper with `find_or_create_allele_index`), Mi1 (replaced `checked_sub.unwrap_or(0)` with `saturating_sub`), and a follow-up `.filter_map → .map` rewrite in `ln_bq_for_read`. The fourth (`ptr_arg` on `resolve_mate_overlap_at_pos`) is a false positive — the function genuinely needs `&mut Vec<_>` for `swap_remove`; addressed via a localised `#[allow(clippy::ptr_arg)]` with a justifying comment.
- **Files changed:** `src/per_sample_caller/pileup/open_record.rs` (filter_map → map; removed dead `OpenPileupRecordTable::iter` after M3 stopped using it), `src/per_sample_caller/pileup/walker.rs` (allow ptr_arg with comment).
- **Validation:** `cargo clippy --lib --all-features -- -D warnings` → no errors in `src/per_sample_caller/pileup/`. Remaining clippy errors are in pre-existing out-of-scope files (variant_grouping.rs, decompression_pool.rs, genotype_merging.rs).

### M5 — `windows(2)` only inspects adjacent contributor pairs

- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The spec invariant ("only mate pairs share a slot") makes the `len > 2` branch unreachable today, but the `windows(2)` shape silently miscompares the moment a future change relaxes the invariant. The fix replaces it with explicit all-pairs iteration plus a `debug_assert!` so any future violation surfaces in tests.
- **Implementation summary:** Replaced `for window in indices.windows(2)` with two nested loops (`for i in 0..n / for j in (i+1)..n`) and added `debug_assert!(indices.len() <= 2)` to make the implicit invariant visible.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** Existing pileup suite (54 tests) still passes — no contributor sharing >2 slots in tests, but the `debug_assert!` would fire if one ever did. No new test added because constructing a contrived three-way slot-sharing scenario would need fixture reads that violate spec invariants the upstream filter is supposed to enforce; the value of the change is structural (no silent miscompare), not behavioral.
- **Files changed:**
  - `src/per_sample_caller/pileup/walker.rs`
- **Tests added or modified:** None.
- **Validation:**
  - `cargo test --lib pileup` → exit 0, 54 passed
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None.

### M3 — `record_widen_events` counter conflates open and widen

- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The previous heuristic (sum of all open-record spans before vs after) incremented on every fresh `open_new` because the after-sum included the new record's span. The right signal lives inside `widen` itself.
- **Implementation summary:**
  - `OpenPileupRecordTable::widen` now returns `Result<bool, _>` (true on actual widening, false on no-op).
  - `process_position` now returns `ProcessOutcome { widen_count: u64 }`, counted as the sum of `widen` calls returning `true`.
  - The walker increments `summary.record_widen_events` by `outcome.widen_count` after each `process_position`, replacing the span-sum delta.
- **Review suggestion used verbatim?:** Yes (return-bool / explicit-counter shape).
- **Adaptation:** Returned the count via a small `ProcessOutcome` struct rather than a `&mut u64` parameter, leaving room for additional per-step counters without growing the call signature.
- **Verification performed:** Added `record_widen_events_counter_only_increments_on_real_widens` (tests.rs, exercises three pure-Match reads at non-overlapping positions; expects counter == 0). Failed before fix (got 9), passed after (got 0). All 54 pileup tests pass.
- **Files changed:**
  - `src/per_sample_caller/pileup/open_record.rs`
  - `src/per_sample_caller/pileup/walker.rs`
  - `src/per_sample_caller/pileup/tests.rs`
- **Tests added or modified:** Added `record_widen_events_counter_only_increments_on_real_widens`; added `drive_walker_with_summary` helper for tests that need to assert on the run summary.
- **Validation:**
  - `cargo test --lib pileup` → exit 0, 54 passed
- **User input:** None.
- **Follow-up:** None. Considered adding a separate `record_open_events` counter; the spec only commits to `record_widen_events`, so leaving opens uncounted matches the spec.
- **Residual risk:** None.

## 5. Deferred findings to carry forward

*(populated at end)*

## 6. Disputed findings to return to reviewer

*(populated at end)*

## 7. Failed-validation findings

*(populated at end)*

## 8. Blocked-by-context-mismatch findings

*(populated at end)*

## 9. Commands run

*(running log)*

## 10. Command results

*(running log)*

## 11. Notes

- Pre-existing clippy errors in `variant_grouping.rs` and `decompression_pool.rs` are out of scope per the source review's §7 ("Out of Scope Observations"). Pileup-scoped clippy will be checked module-locally.
- Tests are run inside the project's dev container via `./scripts/dev.sh` (see `CLAUDE.md`).
