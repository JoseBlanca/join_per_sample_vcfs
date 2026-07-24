# Fix Application Report: ng_normalizer_d1_2026-07-24.md

**Date:** 2026-07-24
**Source review:** `doc/devel/reports/reviews/ng_normalizer_d1_2026-07-24.md`
**Source state reviewed against:** branch `ng-locus-evidence`, uncommitted D1 diff
**Execution mode:** non-interactive (plan-driven)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1 · Majors: 2 · Minors: 3 · Nits: ~3

### Outcome totals
- Applied: 6 · Disputed/recorded (Nits): 3 · Failed validation: 0

### Validation summary
- `cargo fmt --check` → clean · `cargo clippy --all-targets --all-features -- -D warnings` → clean
- `cargo test --example ng_normalizer_screen` → 5/5 · `cargo test --example ng_synthesize_stress_reads` → 3/3 · `cargo test --lib` → 2331 passed
- Performance check → skipped (example binaries; not on a bench path). The screen ran the whole GIAB HG002 10x BAM in 8.3 s.

### Unresolved high-priority findings
None. The reviewer confirmed the measurement logic and generator are correct and reproduced the reported numbers by hand.

## 2. Findings table

| ID | Severity | Title | Decision | Final status |
|---|---|---|---|---|
| B1 | Blocker | `build_loci` untested | Apply | Applied (3 tests) |
| M1 | Major | flank-non-`A` invariant only in a comment | Apply | Applied |
| M2 | Major | D==20 cap boundary untested | Apply | Applied |
| Mi1 | Minor | `Tally` derives `Copy` | Apply | Applied |
| Mi2 | Minor | bare `{error}` drops source chain | Apply | Applied |
| Mi3 | Minor | `matched` bare participle | Apply | Applied (refactor removed it) |
| N* | Nit | `Err(_)` detail, `main` length, canonical-vs-raw | Record | Recorded (no change) |

## 3. Questions asked and answers
None (the adoption/prune direction was raised with the owner separately at Checkpoint D, not blocking D1).

## 4. Per-finding log (condensed)

- **B1 Applied:** `every_read_is_well_formed_and_placed_rightmost` (read/ref-consumption balance + all-`A` middle), `the_shift_distances_straddle_the_cap_as_reported` (recomputes 36 cap-hits / 28 disagreements from the length list), `the_flanks_do_not_extend_the_homopolymer`.
- **M1/M2 Applied:** the flank invariant test above; `an_indel_exactly_at_the_cap_hits_it_but_still_agrees` (shift == 20 → cap hit but leftmost, so agree — the reason cap-hits exceed disagreements).
- **Mi1 Applied:** `Tally` is now `Clone` (not `Copy`) — a mutable accumulator.
- **Mi2 Applied:** an `error_chain` helper walks `source()` so a failure names its cause.
- **Mi3 Applied:** `build_loci` refactored into `append_locus` + `read_sequence`, de-duplicating the deletion/insertion arms and removing the `matched` binding.
- **Nits recorded:** the `Err(_)` recovery keeps counts only (0 observed on GIAB); `main` stays (the tested core `screen_one` is factored out); canonical window fetch is correct against the uppercase read and identical across the three normalizers.

## 5–8.
Deferred: none. Disputed: none. Failed validation: none. Blocked: none.

## 9. Performance check
Skipped — example binaries. (Owner-raised fast-path note: the recommended production normalizer 1a already early-returns on a no-indel read; the screen pre-filters them; documented on 1a in the follow-up default-decision commit.)

## 10–12. Commands / notes
- `cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --example …`, `cargo test --lib`; both examples built release and run on the host.
- Note: the screen was run on GIAB HG002 10x (0 disagreements / 63,757 indel reads) and on the synthetic stress set (28/64 disagreements, all from 1b's cap) — the recorded Checkpoint D result.
