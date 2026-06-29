# Fix Application Report: psp_metadata_section_2026-06-29.md

**Date:** 2026-06-29
**Source review:** `doc/devel/reports/reviews/psp_metadata_section_2026-06-29.md`
**Source state reviewed against:** branch `tomato2-paralog-filter`, A1 working diff
**Execution mode:** non-interactive (autonomous loop)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 3 (M1, M2, M3)
- Minors: 5 (Mi1–Mi5)
- Nits: 2

### Outcome totals
- Applied: 8 (M1, M2, M3, Mi1, Mi2, Mi3, Mi4, + 2 Nits as one)
- Already fixed: 0
- Deferred: 2 (Mi5; the out-of-scope `write_io_err` follow-up)
- Disputed: 0
- Won't fix: 2 (errors advisory; compress→`Io` Display cosmetic)

### Validation summary
- `cargo fmt --check` → exit 0, clean
- `cargo clippy --lib -- -D warnings` → 5 errors, all pre-existing in
  `src/vcf/writer.rs` (`doc_lazy_continuation`); clean with that lint
  allowed (`-A clippy::doc_lazy_continuation` → exit 0). None in scope.
- `cargo test --lib` → 1355 passed, 0 failed (+7 vs pre-fix 1348)
- `cargo doc --no-deps` → fails only on a pre-existing broken link
  (`ClassicStutterModel`, in `src/ssr/`); no psp/metadata doc issues
- `cargo audit` → not run (no dependency change)

### Unresolved high-priority findings
- None. Mi5 deferred by design (cold path, documented cap).

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| M1 | Major | Trailing bytes after frame unchecked | Apply | Applied | metadata.rs, errors.rs | Pass |
| M2 | Major | `pub` in `pub(crate)` mod, doc-link | Apply | Applied | metadata.rs, errors.rs | Pass |
| M3 | Major | No reader-level corruption tests | Apply | Applied | writer.rs (tests) | Pass |
| Mi1 | Minor | Malformed-input test matrix incomplete | Apply | Applied | metadata.rs (tests) | Pass |
| Mi2 | Minor | Frame determinism untested | Apply | Applied | metadata.rs (tests) | Pass |
| Mi3 | Minor | `cap as u64 + 1` overflow footgun | Apply | Applied | metadata.rs | Pass |
| Mi4 | Minor | attach-empty vs no-section untested | Apply | Applied | writer.rs (tests) | Pass |
| Mi5 | Minor | Cap not inspectable at call site | Defer | Deferred | None | N/A |
| Nits | Nit | use-in-fn; `// Mi5` above `//!` | Apply | Applied | metadata.rs | Pass |
| OOS | — | crate-wide `write_io_err` helper | Defer | Deferred | None | N/A |
| adv | — | Decoder::new vs read both `Zstd` | Dispute | Won't fix | None | N/A |
| adv | — | compress→`Io` Display "writing" | Dispute | Won't fix | None | N/A |

## 3. Questions asked and answers

None (non-interactive).

## 4. Per-finding log (applied findings)

### M1 — Trailing bytes after the zstd frame unchecked
- **Final status:** Applied. **Test-first.**
- **Implementation:** `decompress_metadata_capped` now validates the
  region is exactly one frame via
  `zstd::zstd_safe::find_frame_compressed_size(frame) == frame.len()`
  before decoding; a shortfall returns the new
  `PspReadError::MetadataTrailingBytes { trailing }`. This is more robust
  than the reviewer's `Cursor::position` suggestion (which is unreliable
  because the streaming decoder may read-ahead past the frame) — the
  frame API gives the exact compressed length. Truncated frames are also
  rejected here (`find_frame_compressed_size` needs a complete frame).
- **Tests:** `rejects_trailing_bytes_after_frame`,
  `truncated_frame_errors_cleanly`, `non_frame_bytes_error_cleanly`.

### M2 — `pub` codec/const in a `pub(crate)` module
- **Final status:** Applied. Downgraded `MAX_METADATA_DECOMPRESSED_BYTES`,
  `compress_metadata`, `decompress_metadata` to `pub(crate)` (their real
  reach — only in-crate callers). Changed the `errors.rs` doc reference
  from an intra-doc link to a plain code span (no private-intra-doc-link).

### M3 — reader-level corruption tests
- **Final status:** Applied. `reader_rejects_index_overrunning_trailer`
  (past-trailer + `u64`-overflow arms, via in-place trailer rewrite) and
  `reader_rejects_corrupt_metadata_frame` (byte-flip in the frame region,
  caught at `PspReader::new`).

### Mi1 / Mi2 / Mi3 / Mi4 / Nits
- **Mi1:** truncated-frame + non-frame-bytes tests (see M1).
- **Mi2:** `compress_metadata_is_deterministic` byte-equality test.
- **Mi3:** `(cap as u64).saturating_add(1)` in the cap streamer.
- **Mi4:** `attach_empty_metadata_is_some_empty_distinct_from_none`.
- **Nits:** moved `use std::io::Read` to the module import group; moved
  the `// Mi5` shorthand out from above the `//!` module doc.

## 5. Deferred findings to carry forward
- **Mi5** — surface the 64 MiB cap at the call site (accessor / `Debug` /
  `tracing`). Cold open path; the cap is a named, documented `const` and
  the `MetadataSectionTooLarge` error names it. Low value.
- **OOS** — crate-wide `write_io_err(context)` helper to collapse the
  ≥3 repeated `PspWriteError::Io { .. }` closures. Pre-existing-code
  refactor, separate follow-up.

## 6. Disputed findings to return to reviewer
- The two advisory items (both error variants of a frame failure map to
  `Zstd`; a compression failure maps to `Io`) are type-faithful and the
  conflation is harmless / effectively unreachable. **Won't fix.**

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no `Apply` touched a hot path covered by `benches/` (the
section codec runs once per file open, on kilobytes).

## 10. Commands run
- `./scripts/dev.sh cargo test --lib psp::`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings [-A clippy::doc_lazy_continuation]`
- `./scripts/dev.sh cargo doc --no-deps --lib`

## 11. Command results
- `cargo test --lib` → 1355 passed, 0 failed
- `cargo fmt --check` → exit 0
- `cargo clippy` → 5 pre-existing errors in vcf/writer.rs only
- `cargo doc` → pre-existing `ClassicStutterModel` link failure only

## 12. Notes
- M1's fix changes the metadata section's read contract: the region must
  be exactly one zstd frame. This only tightens validation; no valid file
  the writer produces is affected (it emits exactly one frame).
