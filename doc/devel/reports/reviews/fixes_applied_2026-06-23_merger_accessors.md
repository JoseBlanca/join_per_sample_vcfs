# Fix Application Report: ssr_call_merger_accessors_2026-06-23.md

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_merger_accessors_2026-06-23.md`
**Source state reviewed against:** commit `f01c14e`, branch `ssr-cohort`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 1 (Mi1) · Nits: 3

### Outcome totals
- Applied: 1 (md5-assertion Nit) · Deferred: 1 (Mi1 → H3/H4) · No action: 2 Nits

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib ssr::cohort::merge` → 0, 11 passed
- `cargo doc --no-deps` / `cargo audit` → not run (no doc-link or dependency change)
- Performance check → not applicable (accessors + tests, not benched)

### Unresolved high-priority findings
- None. (Mi1 deferred to H3/H4 — the header-build site is the right place to enforce uniqueness; recorded in the plan.)

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|
| Mi1 | Minor | `sample_names()` basenames can collide → duplicate VCF columns | Defer | Deferred | None (plan note) | N/A | H3/H4 (validate uniqueness) |
| Nit-a | Nit | Assert contig `md5` in the accessor test | Apply | Applied | `merge.rs` | Pass | No |
| Nit-b | Nit | `chromosomes()` owned clone vs slice | — | No action | None | N/A | Revisit if H4 awkward |
| Nit-c | Nit | `sample_name_from_label` doesn't strip doubled extensions | — | No action | None | N/A | Correct for well-formed inputs |

## 3. Questions asked and answers
None — Mi1 is a placement decision (enforce at the header-build site), not a user question.

## 4. Per-finding log

### Mi1 — sample-name collision
- **Final status:** Deferred. `sample_names()` is a pure derivation; uniqueness must be enforced where all names are assembled into the header (H3 `write_vcf_header` or H4 driver), returning a typed error before writing a malformed VCF. Added a carry-over note to the H3 plan step.
- **Files changed:** `doc/devel/implementation_plans/ssr_call_driver.md` (plan note).
- **Residual risk:** until H3/H4, two inputs sharing a file name across directories would yield duplicate `#CHROM` columns — but no VCF is emitted yet (H2 only adds the accessor), so the risk is not live.

### Nit-a — assert contig md5
- **Final status:** Applied. Added `assert_eq!(chroms[0].md5, "0".repeat(32))` to `chromosomes_and_sample_names_expose_the_header_table`.
- **Files changed:** `src/ssr/cohort/merge.rs`.
- **Validation:** `cargo test --lib ssr::cohort::merge` → 0, 11 passed.

### Nit-b — owned clone
- **Final status:** No action. Called once per run; the clone is negligible and matches `chrom_names()`'s owned-return convention. Revisit only if H4 finds the owned return awkward.

### Nit-c — doubled extension
- **Final status:** No action. Correct for well-formed `.ssr.psp` inputs; stripping arbitrary doubled extensions is out of contract.

## 5. Deferred findings to carry forward
- Mi1 — enforce sample-name uniqueness when building the VCF header (H3/H4).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no Apply touched perf-sensitive (benched) code.

## 10. Commands run
- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo test --lib ssr::cohort::merge`

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib ssr::cohort::merge` → 0, 11 passed

## 12. Notes
- Reviewed and fixed inline (sub-agent fan-out overloaded on earlier steps); small contained diff.
