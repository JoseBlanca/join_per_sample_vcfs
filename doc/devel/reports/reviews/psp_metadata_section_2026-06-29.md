# Code Review: psp-metadata-section (A1)

**Date:** 2026-06-29
**Reviewer:** rust-code-review skill (orchestrator + 7 category sub-agents)
**Scope:** feature A1 — optional metadata section in the `.psp` format
**Status:** Approve-with-changes

---

## 1. Scope

- Reviewed: working diff for A1 on branch `tomato2-paralog-filter`.
- In-scope files: `src/psp/metadata.rs` (new), `src/psp/writer.rs`
  (`attach_metadata` + `finish()` section write + tests), `src/psp/reader.rs`
  (index-end relaxation + section read + `metadata()`), `src/psp/errors.rs`
  (3 new variants), `src/psp/mod.rs` (module decl).
- Out of scope: rest of crate. Pre-existing clippy `doc_lazy_continuation`
  errors in `src/vcf/writer.rs` are unrelated.
- Categories dispatched: reliability, errors, refactor_safety, defaults,
  extras, idiomatic+naming+smells, module_structure.

## 2. Verdict

**Approve-with-changes.** Implementation logic is correct (gap arithmetic
overflow-safe, cursor positioning valid, trailer pointers intact, zip-bomb
cap correct). One real correctness gap (trailing bytes after the frame),
one visibility-honesty issue, and trust-boundary test gaps.

## 3. Execution status

- `cargo test --lib` → 1348 passed, 0 failed.
- `cargo fmt --check` → exit 0.
- `cargo clippy --lib -- -D warnings` → 5 errors, all pre-existing in
  `src/vcf/writer.rs` (`doc_lazy_continuation`); none in scope (verified
  clean with that lint allowed).
- Needs-verification findings: 0.

## 4. Open questions and assumptions

None blocking. The relaxed reader check (M1/refactor) is the intended
behavior change; the residual risk it documents is closed by M1's fix.

## 5. Top 3 priorities

1. **M1** — reject trailing bytes after the metadata frame (unchecked
   region on the hostile-input boundary).
2. **M3** — add reader-level corruption tests (`IndexOverrunsTrailer`,
   corrupt metadata frame through `PspReader::new`).
3. **M2** — fix the `pub`-in-`pub(crate)` visibility/doc-link mismatch.

## 6. Findings

### Major

**M1: src/psp/metadata.rs — trailing bytes after the zstd frame are silently accepted (unchecked gap region).**
**Categories:** extras, refactor_safety (convergent). Confidence High.
The reader hands the whole gap `[index_end, trailer_start)` to a streaming
decoder, which stops at the first complete frame and ignores trailing
bytes. The trailer's `index_checksum` covers only the index; the zstd
content checksum covers only the frame payload — so `<frame><garbage>`
decodes "successfully" and the post-frame bytes are unchecked. Fix:
validate the gap is exactly one frame via
`zstd::zstd_safe::find_frame_compressed_size(frame) == frame.len()`, else
a new `PspReadError::MetadataTrailingBytes { trailing }` (mirrors
`IndexTrailingBytes`).

**M2: src/psp/metadata.rs, src/psp/mod.rs — `pub` codec/const inside a `pub(crate)` module, not re-exported.**
**Category:** defaults. Confidence High.
`MAX_METADATA_DECOMPRESSED_BYTES` and the codec fns are written `pub` but
the module is `pub(crate) mod metadata` with no re-export, so they are
effectively `pub(crate)` and the `errors.rs` intra-doc link to the const
only resolves in-crate. Fix: downgrade the items to `pub(crate)` to match
their real reach; make the `errors.rs` reference a plain code span (not an
intra-doc link) to avoid a private-intra-doc-link warning.

**M3: src/psp/reader.rs — trust-boundary paths lack reader-level tests.**
**Category:** reliability. Confidence High.
No test drives `IndexOverrunsTrailer` (overflow + past-trailer arms) or a
corrupt metadata frame through `PspReader::new`. A predicate inversion or
a regression to `.unwrap()` on hostile input would pass the current suite.
Fix: add reader-level corruption tests.

### Minor

**Mi1: src/psp/metadata.rs — malformed-input test matrix incomplete.**
(extras) Add truncated-frame and non-frame-bytes cases (each asserting
`Zstd`), and the trailing-bytes case (asserting M1's new variant).

**Mi2: src/psp/metadata.rs — written-frame determinism asserted only indirectly.**
(extras) Add `compress_metadata(p) == compress_metadata(p)` byte-equality.

**Mi3: src/psp/metadata.rs — `cap as u64 + 1` overflows at `usize::MAX`.**
(reliability) Latent footgun in the private test-facing helper (returns
empty silently). Use `(cap as u64).saturating_add(1)`.

**Mi4: src/psp/writer.rs/reader.rs — attached-empty vs no-section untested.**
(reliability) Add a writer→reader test pinning attach(empty) ⇒ `Some([])`
distinct from no-attach ⇒ `None`.

**Mi5: src/psp/metadata.rs, reader.rs — 64 MiB cap not inspectable at the call site.**
(defaults) No accessor / `Debug` / `tracing` for the effective cap.
**Deferred** — cold open path; the cap is a named, documented `const` and
the `MetadataSectionTooLarge` error names it. Not worth a tracing/accessor
addition.

### Nits

- `src/psp/metadata.rs` — `use std::io::Read;` inside the function body;
  move to the module import group.
- `src/psp/metadata.rs:1-3` — the `// Mi5:` shorthand comment sits above
  the `//!` module doc; move below or fold in.

## 7. Out of scope observations

- Duplicated `PspWriteError::Io { context, block_index: None, column_tag:
  None, source }` closure recurs ≥3× crate-wide (writer.rs 248/649/713/
  752/759/854/905/914 + the 2 new A1 sites). A crate-local
  `write_io_err(context)` helper mirroring `reader.rs::io_err` would
  collapse them — a pre-existing-code refactor, separate follow-up.

## 8. Missing tests to add now

- `decompress_metadata_rejects_trailing_bytes_after_frame` (M1).
- `decompress_metadata_rejects_truncated_frame` / `_non_frame_bytes` (Mi1).
- `compress_metadata_is_deterministic` (Mi2).
- `reader_rejects_index_overrunning_trailer` (M3).
- `reader_rejects_corrupt_metadata_frame` (M3).
- `attach_empty_metadata_is_some_empty_distinct_from_none` (Mi4).

## 9. What's good

- Gap arithmetic is overflow-safe by construction
  (`checked_add(...).filter(end <= trailer_start)` → no underflow in the
  later subtraction).
- Zip-bomb guard uses a streaming decoder + `take(cap+1)` instead of
  `decode_all`, bounding memory before materialisation.
- No-section byte-identity is preserved and directly asserted
  (`index_offset + index_byte_length == trailer_start`).
- The two `IndexTrailingBytes` uses (index-internal vs the now-replaced
  positional one) are cleanly separated; `index.rs` keeps its meaning.
- Visibility is minimised and the new module sits correctly beside the
  other `pub(crate)` wire-codec submodules.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --lib psp::`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings -A clippy::doc_lazy_continuation`
