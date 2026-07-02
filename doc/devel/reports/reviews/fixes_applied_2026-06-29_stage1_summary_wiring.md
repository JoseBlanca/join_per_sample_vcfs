# Fix Application Report: stage1_summary_wiring_2026-06-29.md

**Date:** 2026-06-29
**Source review:** `doc/devel/reports/reviews/stage1_summary_wiring_2026-06-29.md`
**Source state:** branch `tomato2-paralog-filter`, C1+C2 working diff
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers 1 (B1); Majors 2 (M1, M2); Minors 2 (Mi1, Mi2); Nits 1.

### Outcome totals
- Applied: 6 (B1, M1, M2, Mi1, Mi2 + the now-moot Nit dropped by B1's fix).
- Deferred / Disputed: 0.

### Validation summary
- `cargo fmt --check` → exit 0
- `cargo clippy --lib --tests --bins -- -D warnings` → clean on scope
  (pre-existing `vcf/writer.rs` doc lints + pre-existing broken
  examples/benches excepted — neither in this diff)
- `cargo test --lib` → 1393 passed, 0 failed
- `cargo test --test pileup_cli_integration --test cohort_cli_integration`
  → 22 passed
- `cargo doc --no-deps` → pre-existing `ClassicStutterModel` failure only

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files |
|---|---|---|---|---|---|
| B1 | Blocker | coverage double-count / order-assert on multi-base records | Apply | Applied | pileup_to_psp.rs |
| M1 | Major | `--gc-window-bp 0` release panic | Apply | Applied | cli.rs |
| M2 | Major | upstream-error masked by finalise failure | Apply | Applied | cli.rs |
| Mi1 | Minor | `gc_window_bp` absent from writer params | Apply | Applied | cli.rs |
| Mi2 | Minor | `observe_record` empty-allele index panic | Apply | Applied | pileup_to_psp.rs |
| Nit | Nit | `offset as u32` cast | Apply | Superseded | pileup_to_psp.rs |

## 3. Questions asked and answers
None.

## 4. Per-finding log (key items)

### B1 — coverage double-count / coordinate-order violation
Applied. `observe_record` now attributes **one** coverage observation per
record, at its anchor position (`record.pos`, REF base = `alleles[0].seq[0]`),
instead of iterating the whole REF span. The walker already emits a separate
per-position record for every covered base, so anchor-only gives exactly one
observation per covered position, in coordinate order — no double-count, no
backwards position. Regression test
`observe_record_handles_multibase_record_in_coordinate_order` feeds a span-4
record then a span-1 record and asserts the fold stays in one tile (the pre-fix
span loop would have tripped the order assert).

### M1 — `--gc-window-bp 0` panic
Applied. `value_parser = clap::value_parser!(u32).range(1..)` makes `0` a clap
diagnostic rather than reaching the accumulator's release `assert!`.

### M2 — error masking on the upstream-break path
Applied. The stashed upstream read error is now surfaced (and the `.tmp`
removed) **before** `summary_acc.finish()` / `attach_metadata` / `finish` /
`finalise_output`, so a secondary failure on the doomed file cannot shadow it.

### Mi1 — provenance
Applied. `gc_window_bp` threaded through `build_writer_header` →
`effective_parameters` and inserted into `[writer.parameters]`; the
`effective_parameters_records_every_knob` test asserts it. The remaining
summary params live in the metadata-section document (noted in the
`effective_parameters` body comment).

### Mi2 — empty-allele guard
Applied. `record.alleles.first()` early-return + `seq.first().copied()`
fallback; no index panic on a malformed record.

### Nit — superseded
The `offset as u32` cast vanished with B1's removal of the span loop.

## 5–8. Deferred / Disputed / Failed / Blocked
None.

## 9. Performance check
Skipped here — the per-record fold is O(1) and on the serial consumer thread;
the end-to-end cost delta is D3's job (a measured cost check on a real run).

## 10. Commands run
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo test --test pileup_cli_integration --test cohort_cli_integration`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --tests --bins --all-features -- -D warnings -A clippy::doc_lazy_continuation`

## 11. Command results
- `cargo test --lib` → 1393 passed, 0 failed
- integration → 22 passed
- `cargo clippy` → clean on scope

## 12. Notes
- B1 is the kind of bug the per-step adversarial review exists to catch: it
  passed 1392 lib tests because no test drove a multi-base record through the
  new seam. The fix is simpler than the original code.
- C1 and C2 were merged into one cycle: a `--gc-window-bp` flag wired to a
  config that nothing consumed (C1 alone) is a dead-code intermediate, so the
  flag + the consumer landed together.
