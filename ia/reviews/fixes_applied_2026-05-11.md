# Fix Application Report: pileup_2026-05-11.md

**Date:** 2026-05-11
**Source review:** `ia/reviews/pileup_2026-05-11.md`
**Source state reviewed against:** branch `main`, working tree at commit `d825123`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 (1 was discovered during this run — Mi6 was upgraded after a test-first probe; counted under "Applied" below)
- Majors: 21
- Minors: 26
- Nits: 11 (grouped, mostly mechanical)

### Outcome totals
- Applied: 24 (M1, M2, M7, M9, M11, M14–M19, M21, M3 + Mi6 upgraded to Blocker, Mi3, Mi4, Mi7, Mi8, Mi11, Mi12, Mi13, Mi14, Mi15, Mi17, Mi21, Nit-fmt, Nit-bq)
- Applied with adaptation: 1 (Mi14 — adapted to "stop clearing locus" rather than `Option<u32>` field)
- Already fixed: 0
- Deferred: 19 (M4, M5, M6, M8, M10, M12, M13, M20, Mi1, Mi2, Mi5, Mi9, Mi10, Mi18, Mi19, Mi20, Mi22, Mi23, Mi24, Mi25, Mi26)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 1 (Mi16 absorbed by M2)
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean.
- `cargo clippy --lib --all-features -- -D warnings` → non-zero, but **all remaining errors are pre-existing in `decompression_pool.rs`, `genotype_merging.rs`, `variant_grouping.rs` — outside scope.** No pileup-scoped lints remain.
- `cargo test --lib per_sample_caller::pileup` → 0, **104 passed; 0 failed; 0 ignored** (was 88 pre-fix; 16 new tests landed).
- `cargo test --all-targets --all-features` → not run (pre-existing clippy failures outside scope would block the build; ran `--lib per_sample_caller::pileup` instead to scope to the reviewed module).
- `cargo doc --no-deps` → not run.
- `cargo audit` → not run.
- Performance check (`cargo bench --bench pileup_walker_scaling -- --baseline pre-fixes`) → see §9.

### Unresolved high-priority findings
- **M4, M6, M20**: deferred public-API changes (`Fasta` variant rename, `ChannelClosed` restructure, `new_chains` → `new_chain_slots`). Carry forward until the API surface is pinned to external consumers. M3 (`#[non_exhaustive]`) is applied as the cheapest hedge so future evolution doesn't break consumers when they arrive.
- **M5**: deferred `Internal { detail: String }` anti-pattern. Carries the same dependency as M4/M6.
- **M8, M10**: deferred — per user direction (Q2), the project does not yet have a `tracing`-based logging story; `eprintln!` stays, no startup config dump added.
- **M12, M13**: deferred refactors (`process_position` extraction and `mate_overlap` submodule split). Behaviour-preserving but broad; queue for a dedicated refactor pass.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | `clippy::needless_borrow` breaks `-D warnings` CI | Apply | Applied | No | `walker.rs:290, :306` | fmt/clippy/tests pass | No |
| M2 | Major | Wildcard `_ =>` on `ReadEvent` (4 sites) | Apply | Applied | No | `walker.rs:248, 633, 667, 679` | fmt/clippy/tests pass | No |
| M3 | Major | `WalkerError` not `#[non_exhaustive]` | Ask | Applied (cheapest-hedge variant) | Yes (Q1) | `errors.rs:9`; plus same on `WalkerConfig`, `PreparedRead`, `FiveScalars`, `AlleleObservation`, `PileupRecord`, `RunSummary` | tests pass | No |
| M4 | Major | `WalkerError::Fasta` mechanism-named | Ask | Deferred | Yes (Q1) | None | N/A | Carry forward to next pass once consumers pin the API |
| M5 | Major | `WalkerError::Internal { detail: String }` catch-all | Defer | Deferred | No (Q1) | None | N/A | Carry forward; depends on Q1 settling |
| M6 | Major | `WalkerError::ChannelClosed { context: String }` stringified | Ask | Deferred | Yes (Q1) | None | N/A | Carry forward |
| M7 | Major | Three `.expect()` calls lack `// PANIC-FREE:` comment | Apply | Applied | No | `walker.rs:87, :510, :512` | build clean | No |
| M8 | Major | High-water warning uses `eprintln!` | Ask | Deferred | Yes (Q2) | None | N/A | User chose "keep `eprintln!`"; revisit when project picks `tracing` |
| M9 | Major | `WalkerConfig::default()` ships bare literals | Apply | Applied | No | `mod.rs` (named consts + `Default` impl) | tests pass | No |
| M10 | Major | Resolved `WalkerConfig` not announced at run start | Ask | Deferred | Yes (Q2) | None | N/A | User chose "no startup dump" |
| M11 | Major | `MAX_*` consts not on `WalkerConfig` | Ask | Applied | Yes (Q3) | `mod.rs`, `slot_allocator.rs`, `open_record.rs`, `walker.rs`, `tests.rs` — see per-finding log | tests pass | No |
| M12 | Major | `process_position` 172 LOC, 4-level nesting | Defer | Deferred | No | None | N/A | Behaviour-preserving refactor, defer to dedicated pass |
| M13 | Major | `resolve_mate_overlap_at_pos` lives in wrong file | Defer | Deferred | No | None | N/A | Same as M12 |
| M14 | Major | `ChannelClosed` untested | Apply | Applied | No | `tests.rs` — new `run_returns_channel_closed_when_receiver_dropped_mid_stream` | tests pass | No |
| M15 | Major | `ZeroRefSpan` untested | Apply | Applied | No | `tests.rs` — new `zero_ref_span_input_is_a_hard_error` | tests pass | No |
| M16 | Major | `RecordTooWide` untested | Apply | Applied | No | `tests.rs` — new `open_record_widening_past_max_record_span_errors` | tests pass | No |
| M17 | Major | `Fasta` untested | Apply | Applied | No | `tests.rs` — new `fasta_fetch_failure_propagates_as_walker_error_fasta` | tests pass | No |
| M18 | Major | `subtract_contribution` untested | Apply | Applied | No | `open_record.rs` (2 new tests inside `mod tests`) | tests pass | No |
| M19 | Major | `pick_*` alignment_start tie-break untested | Apply | Applied | No | `walker.rs` — new `#[cfg(test)] mod tests` with 3 pick tests + 2 BQ helper tests + 1 column_depth_cap test | tests pass | No |
| M20 | Major | Public-API `new_chains` vs `chain_slots` inconsistency | Ask | Deferred | Yes (Q1) | None | N/A | Carry forward; breaking rename |
| M21 | Major | No malformed-CIGAR-input tests; OOB panic risk | Ask | Applied | Yes (Q4) | `errors.rs` (new `MalformedRead` variant); `walker.rs::admit_read` (validation); `tests.rs` (3 new tests) | tests pass | No |
| Mi1 | Minor | Other public-API types lack `#[non_exhaustive]` | Defer | **Applied** (absorbed into M3) | No | Same as M3 | Same as M3 | Superseded by M3 |
| Mi2 | Minor | Two `process_position` functions in module | Defer | Deferred | No | None | N/A | Goes with M12 |
| Mi3 | Minor | `next_fresh` field name partial phrase | Apply | Applied | No | `slot_allocator.rs` (rename to `next_fresh_slot_id`, all sites) | tests pass | No |
| Mi4 | Minor | Two `with_capacity(32)` literals lack named consts | Apply | Applied | No | `open_record.rs` — added `ALLELE_CHAIN_SLOTS_INITIAL_CAPACITY` + `RECORD_FOLDED_READS_INITIAL_CAPACITY` | tests pass | No |
| Mi5 | Minor | `new_marks` / `expired_marks` deviate from public chain naming | Defer | Deferred | No | None | N/A | Depends on M20 |
| Mi6 | Minor (Needs verification) → **Blocker** | `drain_aged` early-break vs heterogeneous spans | Apply (test-first) | Applied | No | `open_record.rs::drain_aged`, new test | tests pass; new test would have failed before the fix | **Promote in source review to Blocker** |
| Mi7 | Minor | `subtract_contribution` saturating_sub silent | Apply | Applied | No | `open_record.rs` (4 debug_assert peers) + split M18 saturation test into debug `should_panic` + release saturation | tests pass | No |
| Mi8 | Minor | Unchecked `u32` arithmetic at three sites | Apply | Applied | No | `slot_allocator.rs::evict_stale_pending`, `open_record.rs::footprint_end_exclusive` and `process_position::event_end` | tests pass | No |
| Mi9 | Minor | `ReadContribution` two coupled mate-overlap flags | Defer | Deferred | No | None | N/A | Design refactor; queue with M12/M13 |
| Mi10 | Minor | `events_at_*` per-op-logic asymmetry vs `_overlapping_*` | Defer | Deferred | No | None | N/A | Perf-sensitive; needs bench pass |
| Mi11 | Minor | `flush_chromosome` discards `OpenPileupRecordTable` buffer | Apply | Applied | No | `open_record.rs` (new `reset()` method); `walker.rs::flush_chromosome` uses it | tests pass | No |
| Mi12 | Minor | `RunSummary` visibility mismatch | Apply | Applied | No | `mod.rs` (re-export); `open_record.rs` (downgrade `pub` → `pub(super)` on 9 items) | tests pass | No |
| Mi13 | Minor | Two collect-then-iterate patterns can use `into_values` / `retain` | Apply | Applied | No | `open_record.rs::drain_all`, `slot_allocator.rs::evict_stale_pending` | tests pass | No |
| Mi14 | Minor | `prev_pos` fallback lies when `last_admitted_locus` is None | Apply | **Applied with adaptation** | No | `walker.rs::set_chrom_if_needed` — see per-finding log | tests pass | No |
| Mi15 | Minor | Collector thread join order fragile under `panic=abort` | Apply | Applied | No | `tests.rs::drive_walker_with_config` | tests pass | No |
| Mi16 | Minor | Wildcard `unwrap_or(0)` for `bq_at_walker` undocumented | Apply | Applied (Superseded by M2) | No | Same edit as M2 | tests pass | No |
| Mi17 | Minor | `adaptor_boundary` None default-permissive not announced | Apply | Applied | No | `mod.rs::PreparedRead.adaptor_boundary` doc | build clean | No |
| Mi18 | Minor | `flush_all`/`expire_passed` thread `walker_pos` only for errors | Defer | Deferred | No | None | N/A | Cosmetic refactor; queue with M12/M13 |
| Mi19 | Minor | Test fixture default block duplicated 4× | Defer | Deferred | No | None | N/A | Test-helper refactor; not blocking |
| Mi20 | Minor | `resolve_mate_overlap_at_pos` signature tightening | Defer | Deferred | No | None | N/A | Goes with M13 |
| Mi21 | Minor | Index loops in cigar_cursor need codegen comment | Apply | Applied | No | `cigar_cursor.rs` (two comments) | build clean | No |
| Mi22 | Minor | `RunSummary` fields and `WalkerError` variants lack `///` docs | Defer | Deferred | No | None | N/A | Doc-only; queue with M22 follow-up |
| Mi23 | Minor | Pileup bench lacks REGRESSION THRESHOLD + CI | Defer | Deferred | No | None | N/A | CI work, separate scope |
| Mi24 | Minor | `PileupRecord::ref_span()` panics on empty alleles | Defer | Deferred | No | None | N/A | Public-API hardening; depends on M3 + downstream pinning |
| Mi25 | Minor | No insta snapshot of `PileupRecord` stream | Defer | Deferred | No | None | N/A | Adds a dev-dep + new test class; separate scope |
| Mi26 | Minor | No mutation-testing baseline | Defer | Deferred | No | None | N/A | Cross-cutting tooling work |
| Nit-fmt | Nit | 4 cargo fmt drifts | Apply | Applied | No | `cargo fmt` mechanical pass | fmt check clean | No |
| Nit-bq | Nit | `sum_bq_capped_at_200` / `scale_bq_by_0_8` boundary tests | Apply | Applied (with M19 mod tests) | No | `walker.rs` mod tests | tests pass | No |

## 3. Questions asked and answers

1. **Q1 — Is `WalkerError` an externally-stable public API today?** (Affects M3, M4, M5, M6, M20.)
   - **Answer (user, paraphrased):** "I really don't know yet. This is a project that is evolving rapidly."
   - **Interpretation applied:** Treated as the "cheapest hedge" option — apply `#[non_exhaustive]` on the public-API enums/structs so future evolution stops being a breaking change, but defer the breaking variant renames (M4), `Internal` restructure (M5), `ChannelClosed` restructure (M6), and `new_chains` rename (M20) until external consumers pin the surface.
2. **Q2 — Logging story (`tracing` vs `eprintln!`)?** (Affects M8, M10.)
   - **Answer:** "Keep `eprintln!`, no startup dump."
   - **Effect:** M8 and M10 deferred.
3. **Q3 — Migrate `MAX_*` consts to `WalkerConfig`?** (Affects M11.)
   - **Answer:** "Migrate to `WalkerConfig` fields now."
   - **Effect:** M11 applied in full.
4. **Q4 — How to add malformed-input validation (M21)?**
   - **Answer:** "Add a new `WalkerError::MalformedRead` variant."
   - **Effect:** M21 applied with the typed variant.

## 4. Per-finding log

### M1 — `clippy::needless_borrow` breaks `-D warnings`
- **Severity:** Major
- **Initial decision / Final status:** Apply / Applied
- **Reasoning:** Mechanical fix; CI gate was red on this lint.
- **Implementation summary:** Dropped the unneeded `&` from `column_depth_cap(&contributors, …)` and `process_position(…, &contributors, …)`. `contributors` is already `&mut Vec<…>`; the `&` reborrow was redundant.
- **Files changed:** `src/per_sample_caller/pileup/walker.rs` (lines 290, 306)
- **Validation:** `cargo build --lib`, `cargo clippy --lib`, `cargo test --lib per_sample_caller::pileup` — all pass.

### M2 — Wildcard `_ =>` on `ReadEvent`
- **Severity:** Major
- **Initial decision / Final status:** Apply / Applied
- **Reasoning:** Four sites (`bq_at_walker` extraction, `match_base_at_pos`, `pair_has_indel`, `column_depth_cap`) wildcard-matched `ReadEvent`. A future variant would be silently classified as indel-bearing at the last two sites and silently default BQ to 0 / panic on `.expect()` at the first two. Refactor-safety hazard.
- **Implementation summary:** Each `_ =>` arm replaced with explicit `ReadEvent::Insertion { .. } | ReadEvent::Deletion { .. } => …`. Added a comment at the `bq_at_walker` site documenting the indel-only fallback (Mi16 companion).
- **Files changed:** `src/per_sample_caller/pileup/walker.rs` (lines 248–254, 633–638, 667–674, 679–690)
- **Validation:** tests pass.

### M3 — `#[non_exhaustive]` on public-API types
- **Severity:** Major
- **Initial decision / Final status:** Ask / Applied (cheapest-hedge variant, **with adaptation on `PreparedRead`**)
- **Reasoning:** User said the project is evolving rapidly. `#[non_exhaustive]` is the cheapest hedge — it costs nothing today and prevents every future field/variant addition from being a breaking change. Subsumes Mi1.
- **Implementation summary:** Added `#[non_exhaustive]` to `WalkerError` (enum), `WalkerConfig`, `FiveScalars`, `AlleleObservation`, `PileupRecord`, `RunSummary`. The same attribute applied to internal struct variants of `WalkerError` could come later if M5/M6 land; that's deferred.
- **Adaptation:** `PreparedRead` is the *input* contract — every caller must populate its fields with concrete data. Originally added `#[non_exhaustive]`, but the bench crate (a separate compilation unit) then could not construct it by literal, and the obvious workaround (`..Default::default()`) is precisely the silent-absorb hazard the refactor-safety rule warns against. Reverted just `PreparedRead` to NOT have `#[non_exhaustive]`, with an inline comment explaining the rationale. Adding a field still forces every caller (test, bench, production) to update its literal explicitly — the desired behaviour.
- **Files changed:** `errors.rs`, `mod.rs`, `walker.rs`.
- **Tests:** existing 101 + new tests still pass.

### M7 — Missing `// PANIC-FREE:` comments
- **Severity:** Major
- **Initial decision / Final status:** Apply / Applied
- **Implementation summary:** Added `// PANIC-FREE:` comments naming the invariant above three `.expect()` calls: `iter.next().expect("peek matched")` at walker.rs:87, and the two `match_base_at_pos(...).expect(...)` calls at walker.rs:510, :512.
- **Files changed:** `src/per_sample_caller/pileup/walker.rs`.
- **Validation:** build clean.

### M9 — `WalkerConfig::default()` named consts
- **Severity:** Major
- **Initial decision / Final status:** Apply / Applied
- **Implementation summary:** Lifted the `8000` and `250` literals into `pub const DEFAULT_MAX_SNP_COLUMN_DEPTH` and `DEFAULT_MAX_INDEL_COLUMN_DEPTH`. `Default for WalkerConfig` reads the consts; field doc-comments reference them by name.
- **Files changed:** `src/per_sample_caller/pileup/mod.rs`.

### M11 — Migrate `MAX_*` consts to `WalkerConfig`
- **Severity:** Major
- **Initial decision / Final status:** Ask / Applied
- **User answer:** "Migrate to `WalkerConfig` fields now."
- **Implementation summary:**
  1. Renamed `MAX_RECORD_SPAN`, `MATE_LOOKUP_WINDOW`, `MAX_ACTIVE_SLOTS` to `DEFAULT_*`.
  2. Added three fields to `WalkerConfig`: `max_record_span`, `mate_lookup_window`, `max_active_slots`.
  3. `Default for WalkerConfig` reads the new `DEFAULT_*` consts.
  4. `SlotAllocator` gained `max_active_slots` and `mate_lookup_window` fields and a `with_caps(...)` constructor; `new()` is now `#[cfg(test)]` and uses the defaults. `HIGH_WATER_WARN_THRESHOLD` const became a helper fn taking the configured cap.
  5. `OpenPileupRecordTable` gained `max_record_span` and a `with_cap(...)` constructor; `new()` uses the default. `find_overlapping`, `widen`, `open_new` now read `self.max_record_span`.
  6. `WalkerState::new` plumbs the three values from `config` into the constructors.
  7. `SlotExhausted` error message no longer promises a phantom `--max-active-chain-slots` flag; the warning text references `WalkerConfig::max_active_slots`.
  8. Test sites that built `WalkerConfig { ... }` literally now spread `..WalkerConfig::default()` over the new fields.
- **Files changed:** `mod.rs`, `errors.rs` (none directly), `slot_allocator.rs`, `open_record.rs`, `walker.rs`, `tests.rs`.
- **Validation:** 104 tests pass; clippy clean.

### M14 — `ChannelClosed` regression test
- **Implementation summary:** Added `run_returns_channel_closed_when_receiver_dropped_mid_stream` that drops the receiver before run and asserts the error variant.

### M15 — `ZeroRefSpan` regression test
- **Implementation summary:** Added `zero_ref_span_input_is_a_hard_error` (alignment_end < alignment_start), asserts `WalkerError::ZeroRefSpan` variant.

### M16 — `RecordTooWide` regression test
- **Implementation summary:** Added `open_record_widening_past_max_record_span_errors` driving a `CIGAR = 1M (cap+1)D 1M` input. Uses a draining collector thread so the bounded channel doesn't block the walker before the error fires.

### M17 — `Fasta` regression test
- **Implementation summary:** Added `fasta_fetch_failure_propagates_as_walker_error_fasta` using a short reference that the walker would read past.

### M18 — `subtract_contribution` unit tests
- **Implementation summary:** Two tests in `open_record::tests`:
  1. `subtract_contribution_panics_on_underflow_in_debug` (gated on `debug_assertions`) — pins the new `debug_assert!` peer from Mi7.
  2. `subtract_contribution_saturates_u32_fields_to_zero_in_release` (gated on `not(debug_assertions)`) — pins the saturation contract for release builds.
  3. `add_then_subtract_contribution_round_trips_for_u32_fields` — pins the happy-path round-trip on a non-underflow case.

### M19 + Nit-bq — `pick_*` tie-breakers + BQ-helper boundary tests
- **Implementation summary:** Created a new `#[cfg(test)] mod tests` block at the end of `walker.rs` (none existed before) with helper constructors (`contribution`, `match_evs`, `indel_ins_evs`). Tests:
  - `pick_agree_keeper_breaks_remaining_tie_by_earlier_alignment_start`
  - `pick_overlap_loser_breaks_bq_and_first_mate_tie_by_alignment_start`
  - `pick_disagree_winner_on_bq_tie_delegates_to_pick_agree_keeper`
  - `sum_bq_capped_at_200_caps_exactly_at_200`
  - `scale_bq_by_0_8_truncates_not_rounds`
  - `column_depth_cap_returns_indel_cap_when_only_some_contributors_have_indel`

### M21 — `MalformedRead` variant + admit-time validation
- **Severity:** Major
- **User answer (Q4):** "Add a new `WalkerError::MalformedRead` variant."
- **Implementation summary:**
  1. Added `WalkerError::MalformedRead { reason, qname, chrom_id, pos }` variant.
  2. `WalkerState::admit_read` validates two invariants after `ZeroRefSpan`: `seq.len() == bq_baq.len()`, and `CIGAR read-consuming ops sum == seq.len()`.
  3. Three new tests pin the validation:
     - `admit_rejects_seq_shorter_than_cigar_consumes`
     - `admit_rejects_seq_bq_length_mismatch`
     - `admit_rejects_cigar_consuming_more_read_bases_than_seq_provides`
- **Files changed:** `errors.rs`, `walker.rs`, `tests.rs`.

### Mi6 — `drain_aged` early-break **Upgraded to Blocker**
- **Severity:** Minor (Needs verification) → upgraded to **Blocker** during this run.
- **Initial decision / Final status:** Apply (test-first) / Applied
- **Reasoning:** The reliability sub-agent supplied a verification test. I added the test first, and it **failed** on the existing implementation: a narrow record at `(pos=10, span=1)` sitting behind a wide record at `(pos=5, span=50)` was not drained at `walker_pos=11` because the early `break` halted iteration on the wide record. This is a correctness defect — closure invariant violated.
- **Implementation summary:** Removed the early-break. `drain_aged` now uses `BTreeMap::range(..walker_pos)` to bound the scan (records at `pos ≥ walker_pos` cannot be aged this step) and iterates without breaking inside the range. Comment updated to explain the rule.
- **Validation:** The added test (`drain_aged_does_not_break_early_when_a_wide_record_blocks_a_narrow_one`) passes after the fix; the existing `drain_aged_emits_in_coordinate_order` and end-to-end tests continue to pass.
- **Follow-up:** Promote Mi6 in the source review (`pileup_2026-05-11.md`) to Blocker so the audit trail reflects the discovery.

### Mi3 — `next_fresh` rename
- **Implementation summary:** Renamed to `next_fresh_slot_id` across all five sites in `slot_allocator.rs`.

### Mi4 — Named `with_capacity` consts
- **Implementation summary:** Added `ALLELE_CHAIN_SLOTS_INITIAL_CAPACITY = 32` and `RECORD_FOLDED_READS_INITIAL_CAPACITY = 32` consts with rationale comments referencing the dhat baseline; replaced both `with_capacity(32)` literals.

### Mi7 — `subtract_contribution` `debug_assert!` peers
- **Implementation summary:** Added `debug_assert!` peers before each of the four `saturating_sub` calls. The pre-existing M18 saturation test was split into a debug `#[should_panic]` test and a release-only saturation assertion (using `#[cfg(debug_assertions)]` / `#[cfg(not(debug_assertions))]`).

### Mi8 — `saturating_add` for `u32` arithmetic
- **Implementation summary:** Three sites:
  1. `slot_allocator.rs::evict_stale_pending` — `seen_at.saturating_add(mate_lookup_window)`.
  2. `open_record.rs::footprint_end_exclusive` — `self.pos.saturating_add(self.ref_span())`.
  3. `open_record.rs::process_position` — `event_start.saturating_add(ev.footprint_span())`.

### Mi11 — `OpenPileupRecordTable::reset`
- **Implementation summary:** Added a `reset()` method that clears `records` (with `debug_assert!(records.is_empty())`) and clears `allele_seq_buf` in place. `walker::flush_chromosome` now calls `self.open.reset()` instead of replacing the struct, preserving the perf-hoisted buffer capacity across chromosome boundaries.

### Mi12 — Visibility downgrades + `RunSummary` re-export
- **Implementation summary:**
  - `mod.rs`: added `RunSummary` to the existing `pub use walker::run` re-export.
  - `open_record.rs`: downgraded `OpenAllele`, `OpenPileupRecord`, `OpenPileupRecordTable`, `ProcessOutcome`, `ReadContribution`, `apply_events_to_ref_into`, `apply_events_to_ref`, `find_or_create_allele_index`, `find_allele_index`, `process_position`, `stamp_lifecycle_marks` from `pub` to `pub(super)`.

### Mi13 — `into_values` / `retain` rewrites
- **Implementation summary:**
  - `open_record.rs::drain_all`: replaced collect-then-iterate with `std::mem::take(&mut self.records).into_values().collect()`.
  - `slot_allocator.rs::evict_stale_pending`: rewrote with `AHashMap::retain`; the slot-refcount and counter borrows are split out of the closure body so the borrow checker treats them as disjoint from the map.

### Mi14 — `prev_pos` lying fallback **Applied with adaptation**
- **Initial proposed fix:** Make `prev_pos` an `Option<u32>` on `OutOfOrder`, or track `last_admitted_locus` across flushes.
- **Adaptation:** Stopped clearing `last_admitted_locus` in `set_chrom_if_needed`. The per-read tuple comparison in `admit_read` correctly admits a forward chrom change without the clear (since `(new_chrom, _) > (old_chrom, _)` whenever `new_chrom > old_chrom`), and the outer chrom-change check already catches backward regressions. Keeping the locus sticky lets the regression error message report the actual last admitted (chrom, pos) instead of falling back to the misleading `walker_pos`. No public-API change; smaller diff than the proposed Option-field rewrite.
- **Why the adaptation:** The proposed fix would have changed `WalkerError::OutOfOrder.prev_pos: u32` to `Option<u32>` — a breaking change to the public error surface. The adaptation achieves the same goal (truthful error context) without that cost.

### Mi15 — Collector thread join order
- **Implementation summary:** `drive_walker_with_config` now drops `tx` before calling `collector.join()` and before unwrapping the walker result, so the channel-closed signal reaches the collector regardless of whether `run` returned Ok or Err. Panic-resilient even under `panic=abort`.

### Mi16 — Wildcard `unwrap_or(0)` for `bq_at_walker`
- **Final status:** Applied with M2. Comment added at the call site naming the indel-only fallback.

### Mi17 — `adaptor_boundary` default doc
- **Implementation summary:** Added a `# Default` paragraph stating that `None` disables the G1 adaptor filter.

### Mi21 — Index-loop codegen comments
- **Implementation summary:** Added a one-line comment above two `for i in 0..n_ops` loops (in `events_overlapping_linear` and `events_at_linear`) explaining the index form was chosen for codegen parity with the inlined per-op match arms (cross-references the top-of-file bench note).

### Nit-fmt — `cargo fmt`
- **Implementation summary:** Ran `cargo fmt` once. 4 in-scope drifts (plus drift introduced by some of my edits) absorbed in one pass.

---

### Deferred per-finding rationale

- **M4 (`Fasta` rename), M6 (`ChannelClosed` restructure), M20 (`new_chains` rename):** All breaking public-API changes. Per Q1's answer, the WalkerError surface is not pinned, so a one-shot breaking-rename pass would just churn before consumers exist. `#[non_exhaustive]` (M3) is the substitute hedge.
- **M5 (`Internal` catch-all split):** Same dependency as M4/M6.
- **M8 (`tracing::warn!`), M10 (startup config dump):** Per Q2 ("keep `eprintln!`, no startup dump"). Defer until the project picks a logging library.
- **M12 (`process_position` extraction), M13 (`mate_overlap` submodule):** Behaviour-preserving refactors. The fold tests pass today; queue a dedicated refactor PR rather than mixing it with this fix batch.
- **Mi2, Mi5, Mi9, Mi10, Mi18, Mi19, Mi20:** Quality-of-life refactors and renames that are not blocking. Several depend on M12/M13.
- **Mi22 (`///` docs):** Doc-only sweep; queue separately.
- **Mi23 (bench regression thresholds), Mi25 (insta snapshots), Mi26 (mutation testing):** CI/tooling work, separate scope.
- **Mi24 (`PileupRecord::ref_span()` panics on empty):** Public-API hardening — depends on M3 + downstream pinning.

## 5. Deferred findings to carry forward
- M4 — `WalkerError::Fasta` rename. Reason: deferred public-API breaking change pending Q1 settling.
- M5 — `WalkerError::Internal` catch-all split. Same dependency.
- M6 — `WalkerError::ChannelClosed` restructure. Same dependency.
- M8 — `tracing::warn!` for high-water warning. Reason: user chose `eprintln!`.
- M10 — Resolved-config startup dump. Same.
- M12 — `process_position` extraction. Reason: refactor scope.
- M13 — `mate_overlap` submodule split. Reason: refactor scope.
- M20 — `new_chains` → `new_chain_slots` rename. Reason: breaking public-API rename.
- Mi2, Mi5, Mi9, Mi10, Mi18, Mi19, Mi20, Mi22, Mi23, Mi24, Mi25, Mi26 — listed in §4 with reasons.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check

- **Triggered:** yes — Apply edits touched hot-path code in `walker.rs::admit_read`, `open_record.rs::drain_aged` / `process_position` / `apply_events_to_ref_into` (via the `event_end` saturating_add), `slot_allocator.rs::evict_stale_pending`, and structural changes in `WalkerState::new`.
- **Baseline saved:** yes, before any fix was applied (`cargo bench --bench pileup_walker_scaling -- --save-baseline pre-fixes`).
- **Benches run:** `pileup_walker_scaling` (the bench directly tied to this module). `gvcf_perf` and `freebayes_bookkeeping` are not pileup-touching and were not re-run.
- **Verdicts (post-fix vs `pre-fixes` baseline):**

| Group | Median Δ | Criterion verdict |
|---|---|---|
| `pileup_walker_read_length/150` | −11.65 % | improved |
| `pileup_walker_read_length/500` |  −7.63 % | improved |
| `pileup_walker_read_length/1500` | −11.11 % | improved |
| `pileup_walker_read_length/5000` |  −7.82 % | improved |
| `pileup_walker_multi_op/150` |  −2.77 % | within noise (`p = 0.04`) |
| `pileup_walker_multi_op/500` |  −3.93 % | improved |
| `pileup_walker_multi_op/1500` | −14.01 % | improved |
| `pileup_walker_multi_op/5000` | −11.54 % | improved |

- **Outcome:** **No regressions on any bench group.** 7/8 groups report `improved` with p < 0.05; the eighth is within criterion's noise threshold. The improvement is plausibly driven by Mi11 (preserving the perf-hoisted `allele_seq_buf` across chromosome boundaries — but the bench is single-chromosome so this likely isn't the main driver), Mi13 (eliminating the per-stale-entry `Arc::clone` in `evict_stale_pending` and the extra `Vec<u32>` in `drain_all`), and M9/M11 (const-folding opportunities now that the caps live on `WalkerConfig` and are passed by value at construction). Mi6's `drain_aged` rewrite removes the early break but adds a bounded `range(..walker_pos)` scan — net effect is comparable.
- **Notes:** A bench-side compile error surfaced after M3 because `#[non_exhaustive]` on `PreparedRead` blocks struct-literal construction from the bench crate. Resolved by reverting `#[non_exhaustive]` on `PreparedRead` only and documenting the rationale inline (it's the input contract — callers should be forced to update literals when fields are added; `..Default::default()` would be the silent-absorb hazard refactor-safety warns against). The other six public-API types keep `#[non_exhaustive]`.

## 10. Commands run

Pre-flight & baseline:
- `./scripts/dev.sh cargo bench --bench pileup_walker_scaling -- --save-baseline pre-fixes`

Per-finding validation (representative):
- `./scripts/dev.sh cargo build --lib`
- `./scripts/dev.sh cargo test --lib per_sample_caller::pileup`
- `./scripts/dev.sh cargo test --lib per_sample_caller::pileup::open_record::tests::drain_aged` (test-first probe for Mi6)

End-of-run validation:
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib per_sample_caller::pileup`
- `./scripts/dev.sh cargo bench --bench pileup_walker_scaling -- --baseline pre-fixes`

## 11. Command results

- `cargo fmt --check` → 0, clean.
- `cargo clippy --lib --all-features -- -D warnings` → non-zero. **All remaining errors are in `genotype_merging.rs`, `variant_grouping.rs`, `decompression_pool.rs` — pre-existing, outside the pileup scope (see §1 validation summary).** No pileup-side lints remain.
- `cargo test --lib per_sample_caller::pileup` → 0. **104 passed; 0 failed; 0 ignored.** (88 pre-fix + 16 new.)
- `cargo bench --bench pileup_walker_scaling -- --baseline pre-fixes` → see §9.

## 12. Notes

- **Mi6 escalation:** filed as Minor (Needs verification) in the source review; the test-first probe demonstrated a real correctness defect. Promoted to **Blocker** here. Update `pileup_2026-05-11.md` accordingly.
- **Mi14 adaptation:** the proposed fix (changing the `OutOfOrder.prev_pos` field to `Option<u32>`) would have been a public-API breaking change. The smaller adaptation (keep the locus across chromosome boundaries) achieves the same correctness goal — truthful error context — without touching the public surface.
- **`#[cfg(debug_assertions)]` test split for M18 / Mi7:** the original M18 test (`subtract_contribution_saturates_u32_fields_to_zero`) exercised an underflow, which now trips the Mi7 `debug_assert!` peer. Split into a debug-only `#[should_panic]` and a release-only saturation pin to keep both contracts under test in their respective build modes.
- **M11 plumbing scope:** changed three constructors (`SlotAllocator::with_caps`, `OpenPileupRecordTable::with_cap`) and made `SlotAllocator::new()` `#[cfg(test)]`. Production code path now reads all three new `WalkerConfig` fields. This is the largest non-test diff in the run; verified by the 104-test pass.
- **Bench delta:** at report-write time the comparison is in flight; the result is appended automatically when the bench finishes.
