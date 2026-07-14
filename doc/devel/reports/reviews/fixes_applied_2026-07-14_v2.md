# Fix Application Report: ng read filtering — Milestone B

**Date:** 2026-07-14
**Source review:** [ng_read_filtering_b_2026-07-14.md](ng_read_filtering_b_2026-07-14.md)
**Source state reviewed against:** working tree on `main`, HEAD `e3a26e6`
**Execution mode:** interactive (one owner decision on M1)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 1 · Minors: 6 · Nits: 4

### Outcome totals
- Applied: 5 · Applied with adaptation: 1 (M1 — kept fatal per owner + added test/doc) · Already fixed: 0 · Deferred: 1 · Disputed: 2 · Failed validation: 0

### Validation summary
- `cargo fmt -- --check` (ng) → clean
- `cargo clippy --lib` → clean
- `cargo test --lib -- ng::read::filtering` → 18 tests pass
- Performance check → Skipped (no `benches/` path touched; #8 already allocation-free per read)

### Unresolved high-priority findings
None. M1 resolved by owner (keep fatal).

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed |
|---|---|---|---|---|---|
| M1 | Major | #8 OutOfBounds fetch fatal vs skip-and-keep | Ask | Applied with adaptation | `read/filtering.rs` (test + doc) |
| Mi1 | Minor | unchecked `read.pos`/`read.ref_id` u32 casts | Apply | Applied | `read/filtering.rs` |
| Mi2 | Minor | DropReason↔counts unenforced | Defer | Deferred (→ D) | None |
| Mi3 | Minor | verdict_* noun-named | Dispute | Disputed | None |
| Mi4 | Minor | stale module doc header | Apply | Applied | `read/filtering.rs` |
| Mi5 | Minor | pre-decode order partially pinned | Apply | Applied | `read/filtering.rs` |
| Mi6 | Minor | mismatch_bq_floor co-dependent field | Dispute | Disputed | None |
| Nit-a | Nit | redundant `None` override | Apply | Applied | `read/filtering.rs` |
| Nit-b | Nit | poly-A ref literal ×7 | Apply | Applied | `read/filtering.rs` |
| Nit-c/d | Nit | terse helpers / `..default()` in tests | Defer | Deferred (kept) | None |

## 3. Questions asked and answers

1. **M1** — On a #8 `OutOfBounds` fetch (read window past contig end), abort the run (fatal) / skip-#8-and-keep / drop the read?
   - **Answer:** Abort the run (keep the fatal model).

## 4. Per-finding log

### M1 — #8 OutOfBounds fetch is fatal
- **Final status:** Applied with adaptation. No behavior change (owner kept fatal). Added `high_mismatch_fetch_past_contig_end_is_fatal` (pins the propagated `OutOfBounds`) and expanded the `verdict_post_decode` doc to state that an out-of-bounds window signals a malformed record and is fatal by design (spec §7).

### Mi1 — unchecked u32 casts
- **Final status:** Applied. `ContigId(u32::try_from(read.ref_id).expect("ref_id fits u32"))` and `u32::try_from(read.pos).expect("read position fits u32")`, with a `PANIC-FREE` note: real contigs are < 4.29 Gbp, so the conversions never truncate on valid input, and a value that didn't fit is a corrupt record to fail loudly on (consistent with M1's fatal model) rather than silently truncate.

### Mi2 — DropReason↔counts unenforced
- **Final status:** Deferred to Milestone D. The exhaustive `match` belongs at the tally site (the iterator), which lands in D. Recorded as a D acceptance criterion.

### Mi3 — verdict_* noun-named
- **Final status:** Disputed. The arch doc §3 mandates `verdict_pre_decode` / `verdict_post_decode`; the code matches its stated interface. Renaming would edit the design doc — out of scope.

### Mi4 — stale module doc header
- **Final status:** Applied. Header now reads "Milestones A … and B … have landed; C and D follow."

### Mi5 — pre-decode order partially pinned
- **Final status:** Applied. Added `pre_decode_charges_filters_in_full_cascade_order`, peeling every filter off in sequence (duplicate → low-MAPQ → supplementary → secondary → unmapped → QC-fail).

### Mi6 — mismatch_bq_floor co-dependent
- **Final status:** Disputed. Deliberate flat mirror of production `AlignmentMergedReaderConfig` (spec §4); same resolution as Milestone A's Mi3.

### Nit-a / Nit-b
- **Final status:** Applied. Removed the redundant `max_read_mismatch_fraction: None`; added `poly_a_ref(n)` and used it across the seven post-decode tests.

### Nit-c / Nit-d
- **Final status:** Deferred (kept as-is — terse helpers are section-commented; `..default()` acceptable in tests).

## 5. Deferred findings to carry forward
- Mi2 — DropReason↔counts exhaustive enforcement → Milestone D.
- Out-of-scope: reused-buffer no-stale assertion → D (plan D3); differential-vs-`classify_pre_decode` → D fixture.

## 6. Disputed findings to return to reviewer
- Mi3 — verdict_* names are arch-mandated.
- Mi6 — flat config is the spec-mandated production mirror.

## 7–8. Failed-validation / blocked
None.

## 9. Performance check
Skipped — no `Apply` touched a `benches/` path; #8's reference read uses a reused scratch buffer (no per-read allocation).

## 10. Commands run
- `./scripts/dev.sh cargo fmt -- src/ng/read/filtering.rs` / `-- --check`
- `./scripts/dev.sh cargo clippy --lib`
- `./scripts/dev.sh cargo test --lib -- ng::read::filtering`

## 11. Command results
- fmt → clean; clippy `--lib` → clean; test → 18 ng::read::filtering pass.

## 12. Notes
- M1 is the one owner-facing decision this milestone; recorded above. The #8/#9 ordering and left-alignment-deferral questions were settled with the owner before implementation (see the commit + the ng step-1 `PROJECT_STATUS.md` block).
