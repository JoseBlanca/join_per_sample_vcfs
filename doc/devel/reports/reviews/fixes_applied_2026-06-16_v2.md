# Fix Application Report: regions_2026-06-16.md

**Date:** 2026-06-16
**Source review:** `doc/devel/reports/reviews/regions_2026-06-16.md`
**Source state reviewed against:** branch `bed-regions-review` @ `9ae294b`
**Execution mode:** interactive
**Overall status:** Completed

*(Distinct from `fixes_applied_2026-06-16.md`, which is the segment-read-fetcher fix run.)*

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 6 (M1–M6) · Minors: 12 (Mi1–Mi12) · Nits: 7 (grouped)

### Open questions (resolved by the user before this run)
1. **Empty `--regions` BED** → **hard error** (gates M1).
2. **`MappedReadSource` intent** → **keep as a test/swap seam**; narrow visibility + document. The "move to `src/bam/`" suggestion is declined — it would invert the dependency direction of a deliberate dependency-inversion seam (gates M5/Mi11).

### Outcome totals
- Applied: 17 (M1, M2, M3, M4, M6, Mi1, Mi4, Mi5, Mi6, Mi7, Mi8, Mi9, Mi10 + Nit-query, Nit-doc, Nit-threads, Nit-open, Nit-track)
- Applied with adaptation: 2 (Mi4 — mode-announce applied, per-line clamp warning deferred; Mi10 — narrowing surfaced + removed dead `Stage1PipelineContext` fields and the `sample_name` threading)
- Already fixed: 0
- Deferred: 3 (M5, Mi12, Nit-bedlines, Nit-cache) *(4 items; M5 + Mi12 + 2 nits)*
- Disputed: 1 (Mi11 — the "move" recommendation; kept as a documented DI seam)
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary (run in the container, branch tip)
- `cargo fmt --all -- --check` → 0, clean (after `cargo fmt --all`)
- `cargo clippy --all-targets -- -D warnings` → 0, clean
- `cargo test --lib` → 0, **1165 passed; 0 failed; 1 ignored**
- `cargo test --tests` → 0, all integration suites pass (pileup 11, cohort 11, + 6 others)
- `cargo doc --no-deps` → 0, clean
- `cargo audit` → not run (cargo-audit not installed — Mi12)
- Performance check → **skipped** (see §9): no Apply touched per-read hot-loop code; baseline not captured pre-edit by that judgment.

### Unresolved high-priority findings
- **M5** (Major) — deferred: duplicated staged/inline interval-walk refactor is regression-prone; needs its own validated effort.

## 2. Findings table

| ID | Severity | Title | Initial | Final status | Files | Validation |
|---|---|---|---|---|---|---|
| M1 | Major | Silent empty output on non-selecting BED | Apply | **Applied** | regions.rs (+test) | Pass |
| M2 | Major | `BedError::Io` context loss | Apply | **Applied** | regions.rs | Pass |
| M3 | Major | Counter-fold merge not field-safe | Apply | **Applied** | alignment_input.rs, baq_stream.rs, walker/driver.rs, stage1_pipeline.rs | Pass |
| M4 | Major | Worker `JoinHandle`s dropped | Apply | **Applied** | stage1_pipeline.rs | Pass (suite; panic-path not unit-testable — see log) |
| M5 | Major | Duplicated staged/inline interval-walk | Defer | **Deferred** | None | N/A |
| M6 | Major | No cross-thread byte-identity test | Apply | **Applied** | tests/pileup_cli_integration.rs | Pass |
| Mi1 | Minor | `bed_start + 1` untrusted arithmetic | Apply | **Applied** | regions.rs (+test) | Pass |
| Mi2 | Minor | `restrict_intervals` test gaps | Apply | **Applied** | var_calling/pipeline.rs (tests) | Pass |
| Mi3 | Minor | BED parser no property test | Apply | **Applied** | regions.rs (proptest) | Pass |
| Mi4 | Minor | Whole-genome / clamp not announced | Apply | **Applied w/ adaptation** | cli.rs | Pass |
| Mi5 | Minor | `cram_paths` param misnamed | Apply | **Applied** | cli.rs | Pass |
| Mi6 | Minor | Bare `#[allow(too_many_arguments)]` ×3 | Apply | **Applied** | stage1_pipeline.rs, var_calling/pipeline.rs | Pass |
| Mi7 | Minor | `RegionSet` no `IntoIterator` for `&` | Apply | **Applied** | regions.rs | Pass |
| Mi8 | Minor | `sort_and_merge` in/out-param | Apply | **Applied** | regions.rs | Pass |
| Mi9 | Minor | clamp → `(start..=end).contains` | Apply | **Applied** | pileup_to_psp.rs | Pass |
| Mi10 | Minor | Stage-1 surface `pub` → `pub(crate)` | Apply | **Applied w/ adaptation** | pop_var_caller/mod.rs, stage1_pipeline.rs, cli.rs | Pass |
| Mi11 | Minor | `MappedReadSource` placement | Dispute | **Disputed** (+doc) | stage1_pipeline.rs | Pass |
| Mi12 | Minor | `cargo audit` not wired | Defer | **Deferred** | None | N/A |
| Nit-query | Nit | Stale `query` comments | Apply | **Applied** | cli.rs, stage1_pipeline.rs | Pass |
| Nit-doc | Nit | Braided `Stage1Walker` doc | Apply | **Applied** | stage1_pipeline.rs | Pass |
| Nit-threads | Nit | Stale threads-default comment | Apply | **Applied** | cli.rs | Pass |
| Nit-open | Nit | `BedError::Open` no `#[source]` | Apply | **Applied** | regions.rs | Pass |
| Nit-track | Nit | `track`/`browser` prefix vs token | Apply | **Applied** | regions.rs (+test) | Pass |
| Nit-bedlines | Nit | CRLF/whitespace/Unicode tests | Defer | **Deferred** | None | N/A |
| Nit-cache | Nit | cache-then-`expect` dup | Defer | **Deferred** | None | N/A |

## 3. Questions asked and answers

1. **M1 (Q1)** — Should an empty `--regions` BED error, warn, or be silent? → **Hard error** (`BedError::NoRegions`, raised at `from_bed_reader`).
2. **M5/Mi11 (Q2)** — Is `MappedReadSource` a deliberate seam or vestigial? → **Deliberate test/swap seam**; narrow visibility + document, keep placement.

## 4. Per-finding log

### M1 — silent empty output on a non-selecting `--regions` BED
- **Final status:** Applied. `RegionSet::from_bed_reader` now returns the new `BedError::NoRegions` when the parsed set is empty; both the pileup (`cli.rs`) and cohort (`pipeline.rs`) entry points inherit it via `from_bed_path`. Test `empty_or_comment_only_bed_is_a_hard_error` (empty / comment-only / header-only) replaces the old `empty_bed_yields_empty_region_set` (whose contract this inverts). The "non-empty regions, zero coverage" case is untouched (stays valid-empty).
- **Files:** `src/regions.rs`. **Tests:** `empty_or_comment_only_bed_is_a_hard_error`.

### M2 — `BedError::Io` dropped path/line/`#[source]`
- **Final status:** Applied. Renamed `Io(io::Error)` → `ReadLine { line_number, #[source] source }` (carries the line, preserves the source in the chain, operation-named). Nit-open folded in: dropped the duplicated `{source}` from `Open`'s `Display` (its named `source` field already chains). **Files:** `src/regions.rs`.

### M3 — counter-fold `merge()`s not field-addition-safe
- **Final status:** Applied. `FilterCounts::merge`, `BaqSkipCounts::merge`, `RunSummary::merge`, and the `StageCounts`→`FilterCounts` fold now exhaustively destructure `*other` (no `..`) so a new field is a compile error. The `RunSummary` comment calls out the non-uniform fold (`active_reads_high_water` = max). **Files:** `alignment_input.rs`, `baq_stream.rs`, `walker/driver.rs`, `stage1_pipeline.rs`.

### M4 — `run_pipelined` dropped worker `JoinHandle`s
- **Final status:** Applied. Worker handles are collected and joined explicitly after the producer join, surfacing a worker panic at a defined point (matching the producer/cohort pattern) rather than relying on `thread::scope`'s implicit re-raise; a doc comment records why (the truncation-then-panic margin). **Verification:** full suite passes; the panic-mid-fold path is not unit-testable without fault injection (noted) — the change is behavior-preserving on the happy path and strictly tightens the failure path. **Files:** `stage1_pipeline.rs`.

### M5 — duplicated staged/inline interval-walk
- **Final status:** Deferred. The duplication is in the regression-prone cohort producer (`pipeline.rs:431-595`, the prior staged multi-interval-bug area); a shared-helper extraction is a wider refactor with multiple shapes and must be validated on its own (byte-identity gate). Carried forward.

### M6 — no cross-thread byte-identity test for `--regions`
- **Final status:** Applied. Added `regions_pileup_records_are_thread_count_invariant` (runs the same `--regions` BED at `--threads 1` inline vs `4` staged, asserts equal record streams). Compares records, not raw bytes, because the `.psp` header legitimately encodes the thread count. **Files:** `tests/pileup_cli_integration.rs`.

### Mi1 — `bed_start + 1` untrusted arithmetic
- **Applied.** `saturating_add(1)` + comment; test `parse_bed_line_rejects_oversized_u64_start_coordinate`. **Files:** `regions.rs`.

### Mi2 — `restrict_intervals_to_regions` test gaps
- **Applied.** Added `restrict_intervals_advances_covered_on_end_equals_region_end_excl` and `restrict_intervals_one_region_spanning_three_covered_intervals`. **Files:** `var_calling/pipeline.rs`.

### Mi3 — BED parser no property test
- **Applied.** `proptest` `region_set_invariants_hold_on_random_spans` (sorted, disjoint/non-abutting, in-bounds; empty input → `NoRegions`). **Files:** `regions.rs`.

### Mi4 — whole-genome / clamp not announced (Applied with adaptation)
- **Applied:** an `eprintln!` now announces the resolved analysis mode (N spans from BED, or whole genome) at the start of `run_pileup`. **Adapted/deferred:** the per-line *clamp* warning (silent overhang clamp) is deferred — it requires threading a clamp flag out of `parse_bed_line`; clamping is already documented and rare, so the smaller half is left as a follow-up. **Files:** `cli.rs`.

### Mi5 — `cram_paths` misnamed
- **Applied.** `build_writer_header`'s `cram_paths` → `alignment_paths` (BAM-or-CRAM). **Files:** `cli.rs`.

### Mi6 — bare `#[allow(too_many_arguments)]`
- **Applied.** Justifying comments on `with_stage1_chain`, `run_pipelined`, `dust_mask_for` (config-threading / distinct concerns; a bundle struct would not reduce real arity). **Files:** `stage1_pipeline.rs`, `var_calling/pipeline.rs`.

### Mi7 / Mi8 / Mi9 — idiomatic
- **Applied.** `impl IntoIterator for &RegionSet` (Mi7); `sort_and_merge` now consuming `fn(Vec)->Vec` (Mi8); `(start..=end).contains(&pos)` clamp (Mi9). **Files:** `regions.rs`, `pileup_to_psp.rs`.

### Mi10 — Stage-1 surface visibility (Applied with adaptation)
- **Applied:** `pub mod stage1_pipeline` → `pub(crate)`; `MappedReadSource`/`Stage1Walker` → `pub(crate)`. **Adaptation:** the narrowing surfaced a real dead-code warning — `Stage1PipelineContext.{sample_name,contigs}` were never read (the closure uses only `walker`). Removed those fields and, consequently, the now-unused `sample_name` threading through `with_stage1_chain`/`run_pipelined` and the caller. This extends L3's dead-weight removal; compiler-verified. **Files:** `pop_var_caller/mod.rs`, `stage1_pipeline.rs`, `cli.rs`.

### Mi11 — `MappedReadSource` placement
- **Disputed (the "move" recommendation).** Per Q2 it is a deliberate dependency-inversion seam (the consuming stage declares the input contract; a reader satisfies it), kept for testability/swap. Moving it into `src/bam/` would invert the dependency. Applied the fallback the finding offered: a doc note stating the intent + unbraiding the `Stage1Walker` doc (Nit-doc). **Files:** `stage1_pipeline.rs`.

### Mi12 — `cargo audit` not wired
- **Deferred.** cargo-audit is not installed (`error: no such command: audit`). No new dependencies this branch, so incremental supply-chain risk is nil; wiring it (install + `deny.toml`) is a tooling follow-up.

### Nits
- **Nit-query** (applied): the present-tense `query`/`query()` comments now name `SegmentMergedReads`; the past-tense "replaces the old `AlignmentMergedReader::query`" kept. `cli.rs`, `stage1_pipeline.rs`.
- **Nit-doc** (applied): `Stage1Walker` doc unbraided from the `MappedReadSource` trait. `stage1_pipeline.rs`.
- **Nit-threads** (applied): the threads-default comment now states the pileup default (`DEFAULT_PILEUP_THREADS`), not "all logical cores". `cli.rs`.
- **Nit-open** (applied, with M2): `BedError::Open` no longer flattens `{source}` into `Display`. `regions.rs`.
- **Nit-track** (applied): `track`/`browser` matched as a whole first token, not a prefix; test `contig_named_like_a_header_keyword_is_not_skipped`. `regions.rs`.
- **Nit-bedlines / Nit-cache** (deferred): CRLF/whitespace/Unicode line tests (partly covered by the new proptest) and the cosmetic `ref_fetch`/`dust_mask_for` cache-then-`expect` duplication — low value, left as follow-ups.

## 5. Deferred findings to carry forward
- **M5** — staged/inline interval-walk de-duplication (regression-prone; own effort).
- **Mi12** — wire `cargo audit` + `deny.toml`.
- **Nit-bedlines**, **Nit-cache** — minor test/cosmetic follow-ups.

## 6. Disputed findings to return to reviewer
- **Mi11** — the "move `MappedReadSource` to `src/bam/`" recommendation; kept as a documented dependency-inversion seam per Q2.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** No.
- **Baseline saved:** No — not captured before editing, by the judgment that no Apply touches per-read hot-loop code.
- **Reasoning:** the changes are parse-time (`regions.rs`), per-region orchestration/teardown (`cli.rs`, `stage1_pipeline.rs`), per-region counter folds (not per-read), a codegen-identical per-record clamp (`pos>=s && pos<=e` ⇄ `(s..=e).contains`), error types, tests, and visibility. None is in the BAQ/walker per-read hot loop. The 2026-06-16 perf review already measured the entire `--regions` read/merge path at ~1–2% of runtime (BAQ dominates ~70%), so a regression here is not plausible. Per the skill, since no pre-edit baseline was captured under that judgment, the comparison is recorded-and-skipped rather than back-filled.

## 10. Commands run
- `cargo build --lib` (×2, mid-run checkpoints)
- `cargo fmt --all` ; `cargo fmt --all -- --check`
- `cargo clippy --all-targets -- -D warnings`
- `cargo test --lib` ; `cargo test --tests`
- `cargo doc --no-deps`
- (all via `./scripts/dev.sh` in the container)

## 11. Command results
- `cargo fmt --all -- --check` → 0, clean
- `cargo clippy --all-targets -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 1165 passed / 0 failed / 1 ignored
- `cargo test --tests` → 0, all suites pass (pileup 11, cohort 11, +6 others, +3 main)
- `cargo doc --no-deps` → 0, clean
- `cargo audit` → not run (not installed)

## 12. Notes
- M4's panic-mid-fold failure path is not exercised by a unit test (would need fault injection into a worker); the fix is behavior-preserving on the happy path (full suite green) and strictly tightens the failure path. Flagged honestly rather than asserting coverage it doesn't have.
- Mi10 expanded slightly beyond a visibility tweak because the narrowing *surfaced* genuine dead weight (the `sample_name` threading); removing it is in the spirit of the review (cf. L3) and is compiler-verified.
