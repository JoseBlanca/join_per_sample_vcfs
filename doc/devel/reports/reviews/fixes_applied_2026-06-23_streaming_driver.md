# Fix Application Report: ssr_call_streaming_driver_2026-06-23.md

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_streaming_driver_2026-06-23.md`
**Source state reviewed against:** commit `d17f524`, branch `ssr-cohort`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 3 (Mi1, Mi2, Mi3) · Nits: 2

### Outcome totals
- Applied: 2 (Mi2 cross-thread determinism test; Mi1 doc + carry-over) · Applied with adaptation: 0
- Already fixed: 0 · Deferred: 1 (Mi2 filtered-record e2e) · Disputed: 0
- Nits: 2 applied (sweep-order comment; the 7-arg one = no action by the review itself)
- Mi3: Applied (documented in the `UnresolvedSamples` error)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, **1271 passed, 0 failed, 2 ignored** (+1 determinism test)
- `cargo doc --no-deps` / `cargo audit` → not run (doc + test changes; no new dep)
- Performance check → not applicable (driver pipeline is not covered by a `benches/` harness)

### Unresolved high-priority findings
- None. (Mi1's real fix — representative subset selection — is the deferred calibration item; documented in code + plan.)

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | Positional subset biases burn-in + couples to decision E | Apply (doc) | Applied | `driver.rs` (doc), plan note | Pass |
| Mi2 | Minor | E2E gaps: cross-thread determinism + filtered record | Apply / Defer | Applied (determinism) + Deferred (filtered e2e) | `driver.rs` (test) | Pass |
| Mi3 | Minor | Zero-covered-loci surfaces as `UnresolvedSamples` | Apply (doc) | Applied | `driver.rs` (doc) | Pass |
| Nit-a | Nit | Note serial sweep preserves order (`_seq`) | Apply | Applied | `driver.rs` (comment) | Pass |
| Nit-b | Nit | `genotype_locus` 7 args | — | No action | None | N/A |

## 3. Questions asked and answers
- **OQ1 (decision-E on the subset):** acknowledged; the fix is representative selection (calibration). Documented in the `BURN_IN_MAX_LOCI` caveat and carried as a plan/PROJECT_STATUS open item.
- **OQ2 (H1-Mi1 backfill universe):** resolved as accept-the-documented-boundary (benign while `fallback_p == EM default`). No code change.

## 4. Per-finding log

### Mi1 — positional subset / decision-E coupling
- **Final status:** Applied (documentation + carry-over). Extended the `BURN_IN_MAX_LOCI` doc with the decision-E coupling caveat; the substantive fix (reservoir/stratified selection) stays the deferred calibration item. No behaviour change (correct at the cap for real genomes; the integration cohort is far under the cap).
- **Files changed:** `src/ssr/cohort/driver.rs` (doc); PROJECT_STATUS open item.

### Mi2 — e2e test gaps
- **Final status:** Applied (cross-thread determinism) + Deferred (filtered-record e2e). Added `run_is_byte_identical_across_thread_counts` (run at `threads = 1` and `4`, assert identical VCF bytes) — proves the headline property end to end. A filtered-record e2e needs a deterministically non-PASS-admitting locus; deferred to calibration (the FILTER formatting is covered by `vcf_out::formats_a_no_call_filtered_record`).
- **Files changed:** `src/ssr/cohort/driver.rs` (test). **Tests added:** `run_is_byte_identical_across_thread_counts`.
- **Validation:** `cargo test --lib ssr::cohort::driver` → 0, 11 passed.

### Mi3 — zero-covered-loci clarity
- **Final status:** Applied. Documented in the `UnresolvedSamples` error doc that a cohort with no covered loci surfaces here (every sample unresolved). A distinct error was judged not worth the surface for an edge case.
- **Files changed:** `src/ssr/cohort/driver.rs` (doc).

### Nit-a — sweep order comment
- **Final status:** Applied. Added a comment that the serial sweep preserves catalog order so the `_seq` reorder is Milestone J's concern.

### Nit-b — `genotype_locus` arg count
- **Final status:** No action. 7 args is under clippy's threshold; bundling the cfgs is premature.

## 5. Deferred findings to carry forward
- Mi2 (filtered-record e2e) — defer to calibration; FILTER formatting already unit-tested.
- Mi1 (representative subset selection) — the deferred calibration item (reading Q-R5), now also motivated by the decision-E coupling.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — the driver pipeline is not reachable from any `benches/` harness.

## 10. Commands run
- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo test --lib`

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 1271 passed / 0 failed / 2 ignored

## 12. Notes
- Reviewed and fixed inline (sub-agent fan-out repeatedly overloaded). The integration test is the genuine DoD check; the determinism test now also pins the headline byte-identity end to end.
