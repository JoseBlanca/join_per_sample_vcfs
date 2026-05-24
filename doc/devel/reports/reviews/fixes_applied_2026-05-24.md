# Fix Application Report: bam_input_support_2026-05-24

**Date:** 2026-05-24
**Source review:** `doc/devel/reports/reviews/bam_input_support_2026-05-24.md`
**Source state reviewed against:** commit `1c59d61` (post-review HEAD)
**Execution mode:** interactive (two AskUserQuestion decisions: CSI depth, scope)
**Overall status:** _in progress — preflight complete; commits 1–10 pending_

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 19 (M1–M19)
- Minors: 22 (Mi1–Mi22)
- Nits: grouped (doc-sweep + small idiomatic items)

### Outcome totals (to be filled in as work proceeds)
- Applied: 0 → target ~31
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 0 → target 7 (Mi5, Mi7, Mi12, Mi13, Mi15, Mi20, Mi21; plus M17's redesign half)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0 (two questions answered in preflight)

### Validation summary (to be filled in at the end)
- `cargo fmt --check` → _pending_
- `cargo clippy --all-targets --all-features -- -D warnings` → _pending_
- `cargo test --all-targets --all-features` → _pending_
- `cargo doc --no-deps` → _expected to remain at 14 pre-existing warnings; will assert no new ones_
- `cargo audit` → not run (not in project verification list)
- Performance check (`cargo bench -- --baseline pre-fixes`) → **skipped** — no `Apply` finding touches the hot per-record decode path; the merge / per-chrom worker hot path receives only startup-time changes (helper extraction, dispatch enumeration) and one error-path enum addition. Baseline was not captured (skill rule: do not stash/revert to back-fill).

### Unresolved high-priority findings (to be filled in)
- _pending_

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | Error vocabulary still says "CRAM" on BAM inputs | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M2 | Major | `VarCallingFromBamCliError::Io(#[from] io::Error)` collapses multiple origins | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M3 | Major | `load_alignment_index` returns bare `io::Error` | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M4 | Major | `AlignmentInputError` missing `#[non_exhaustive]` | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M5 | Major | `load_per_input_headers` duplicates opener + skips CRAM-version gate | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M6 | Major | Catch-all `(file_kind, index)` arm absorbs future variants | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M7 | Major | `BamIndex` enum missing `#[non_exhaustive]` | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M8 | Major | `PileupCliError` missing `#[non_exhaustive]` | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M9 | Major | `AlignmentIndexFormatMismatch` has no test | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M10 | Major | `AlignmentInputError::UnsupportedExtension` test gap on pileup path | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M11 | Major | `OwnedIndexedBamRecords::next` chunk-walking + IO-error untested | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M12 | Major | `load_alignment_index` no direct unit tests | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M13 | Major | `.csi`/`.bai` policy not in a named constant | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M14 | Major | `build_csi_for_bam` uses `Indexer::default()`; depth hidden | Apply | _pending_ | **Yes — depth=6** | _pending_ | _pending_ | No |
| M15 | Major | `--build-map-file-index` help text doesn't say `.csi`-only build | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M16 | Major | `MissingMapFileIndex` Display lists only one path for BAM | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M17 | Major | `MixedAlignmentFileFormats` dual-surfaced; only one path tested | Apply (test only) | _pending_ | No | _pending_ | _pending_ | Yes — defer the shared-enum redesign |
| M18 | Major | `.unwrap()` defensive panic after classify pre-pass | Apply | _pending_ | No | _pending_ | _pending_ | No |
| M19 | Major | Impl-report commit table + deferred-list inaccurate | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi1 | Minor | `from_open_crams` method + param `open_crams` stale | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi2 | Minor | Stale local var names in dispatch loops | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi3 | Minor | `cram_cfg` + `cram_config_from_args` stale | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi4 | Minor | `CramHeader` struct name | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi5 | Minor | `process_one_chromosome_from_bam` name covers CRAM too | **Defer** | Deferred | No | None | N/A | Yes — held deliberately per impl-report; subcommand-name coupling |
| Mi6 | Minor | `validate_fasta_agreement` param + detail strings | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi7 | Minor | `input_crams` field/local mixed-format | **Defer** | Deferred | No | None | N/A | Yes — PSP schema field reaches into out-of-scope `src/psp/header.rs` |
| Mi8 | Minor | Doc-comment + module-doc sweep | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi9 | Minor | `.csi`/`.bai` scan-and-index loop triplicated | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi10 | Minor | `index_display_name` free function → impl method | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi11 | Minor | `load_alignment_index` re-encodes csi/bai policy | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi12 | Minor | `MixedAlignmentFileFormats` / `UnsupportedExtension` duplicated | **Defer** | Deferred | No | None | N/A | Yes — folded into M17's redesign half (deferred together) |
| Mi13 | Minor | `From<AlignmentIndexError>` 4-arm passthrough | **Defer** | Deferred | No | None | N/A | Yes — design call; the explicit rename is load-bearing for the CLI vocabulary |
| Mi14 | Minor | `pub mod cram_input` / `pub mod bam_input` → `pub(crate) mod` | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi15 | Minor | `crate::bam` module name rename | **Defer** | Deferred | No | None | N/A | Yes — large blast radius; decide before next API consumer outside `crate::bam` |
| Mi16 | Minor | `BamIndex::Csi` arm has no isolated test | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi17 | Minor | `OwnedBamRecords` `Err(e)` arm has no test | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi18 | Minor | `Stage1` variant doc says "CRAM-input validation" | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi19 | Minor | `BinnedIndex` vs `LinearIndex` undocumented | Apply | _pending_ | No | _pending_ | _pending_ | No |
| Mi20 | Minor | `OwnedIndexedBamRecords::next` long phase-numbered loop | **Defer** | Deferred | No | None | N/A | Yes — cosmetic; revisit when a second use case appears |
| Mi21 | Minor | Per-record `RecordBuf::default()` allocation | **Defer** | Deferred | No | None | N/A | Yes — perf-tuning workstream candidate; same pattern as CRAM analogue |
| Mi22 | Minor | `Cargo.toml` missing justification comment | Apply | _pending_ | No | _pending_ | _pending_ | No |

## 3. Questions asked and answers

1. **M14** — what CSI depth should `build_csi_for_bam` use?
   - **Answer:** Bump to depth=6 (~2^32 ≈ 4.3 Gbp addressable). Safe for wheat / large plant genomes; future-proofs against plant references the project doesn't target today but might tomorrow.

2. **Scope** — how much of the review should the fix run cover?
   - **Answer:** All Majors + coupled Minors + the 6 missing tests. Defer Mi5, Mi7, Mi12, Mi13, Mi15, Mi20, Mi21, and M17's redesign half.

3. **M5 / review §4 Q3** — was `load_per_input_headers` skipping the CRAM-version check intentional?
   - **Answer:** No — the skip was accidental. Fix it. Prioritised M5 to the first commit so the bug is closed early.

## 4. Per-finding log

_To be filled in as fixes are applied. One section per finding._

## 5. Deferred findings to carry forward

- **Mi5** — `process_one_chromosome_from_bam` function-name rename. Held deliberately per the original impl-report follow-ups; the subcommand-name coupling means renaming either both or neither.
- **Mi7** — `input_crams` field-name rename. Reaches into out-of-scope `src/psp/header.rs`; PSP schema field rename has its own design pass.
- **Mi12** — `MixedAlignmentFileFormats` + `UnsupportedExtension` duplication across the two error enums. Folded into M17's deferred redesign half.
- **Mi13** — `From<AlignmentIndexError>` 4-arm passthrough. Design call; the explicit rename is load-bearing for the CLI's `--flag` vocabulary.
- **Mi15** — `crate::bam` module rename to `crate::alignment`. Large blast radius; decide before the next API consumer outside `crate::bam` is added.
- **Mi20** — `OwnedIndexedBamRecords::next` phase-extraction. Cosmetic; revisit when a second use case appears.
- **Mi21** — Per-record `RecordBuf::default()` allocation. Folds into the standing parallelisation-tuning workstream.
- **M17 (redesign half only)** — Lift `MixedAlignmentFileFormats` into a shared `AlignmentFormatError` sub-enum. Apply only the test-add piece in this run.

## 6. Disputed findings to return to reviewer

_To be filled in._ Expected: none.

## 7. Failed-validation findings

_To be filled in._ Expected: none.

## 8. Blocked-by-context-mismatch findings

_To be filled in._ Expected: none.

## 9. Performance check

- **Triggered:** no
- **Baseline saved:** no — see Notes.
- **Benches run:** none.
- **Verdicts:** N/A.
- **Outcome:** skipped.
- **Notes:** All planned `Apply` fixes touch startup-time / error-path / dispatch-time code. The per-record hot path in `OwnedBamRecords::next` / `OwnedIndexedBamRecords::next` / merge loop is unchanged by any fix. Baseline not captured per the skill's "do not stash/revert to back-fill" rule.

## 10. Commands run

_Pending._ Each commit's validation block (cargo fmt / clippy / test) will be listed here.

## 11. Command results

_Pending._

## 12. Notes

- The original review's §4 open questions Q1 (rename strategy) and Q3 (CRAM-version-gating skip) are resolved by the applied fixes: M1's rename pass continues the format-neutral rename already in motion; M5's helper extraction restores the CRAM-version gate. Q2 (`#[non_exhaustive]` omission) is treated as an oversight per the plan's explicit call (closed by M4/M7/M8). Q4 is answered above.
- The two test files `tests/pileup_cli_integration.rs` and `tests/cohort_cli_integration.rs` assert on the format-name strings (`"CRAM"`, `"BAM"`) — these stay stable across the M1 rename pass.
- The `cargo doc` failure is pre-existing per `git blame`; PROJECT_STATUS.md already tracks it. This fix run will not address the pre-existing intra-doc-link warnings (out of scope).
