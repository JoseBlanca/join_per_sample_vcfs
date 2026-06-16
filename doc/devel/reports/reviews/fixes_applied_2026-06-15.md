# Fix Application Report: psp_container_generalization_2026-06-15.md

**Date:** 2026-06-15
**Source review:** `doc/devel/reports/reviews/psp_container_generalization_2026-06-15.md`
**Source state reviewed against:** branch `ssr-architecture`, `git diff aa6a105 f0dfffc -- src/psp/`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 5 (M1–M5)
- Minors: 11 (Mi1–Mi11)
- Nits: grouped (record_interval tuple, KNOWN_KINDS table, projected_bytes consts, spanning type, UNREACHABLE markers, dead-code expect, private_bounds, unknown-column skip alloc, stale `key_ssr` doc, frr/amb glossary)

### Outcome totals
- Applied: 14 (M1, M2, M3, M4, M5, Mi2, Mi3, Mi5, Mi6, Mi7, Mi8, Mi9, Mi10, Mi11)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 2 (Mi1, Mi4) + 4 deferred Nits
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

Three Majors/Minors that hinged on design/policy choices were resolved by the user up front (see §3): **M4** → poison-the-iterator, **M5** → `column_key!` macro, **Mi3** → require the `kind` tag.

### Validation summary
- `cargo fmt --check` → 0, clean (after one `cargo fmt` to wrap a long `matches!` in a new test).
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean (`Finished … in 3.49s`, 0 warnings).
- `cargo test --lib` → 0, `test result: ok. 1142 passed; 0 failed; 1 ignored` (was 1138; +4 net new psp tests).
- `cargo doc --no-deps` → 0, `Finished` (no warnings).
- `cargo audit` → not run — `cargo-audit` is not installed in the dev container (`no such command: audit`). No dependencies were added, so the advisory surface is unchanged.
- Performance check (`cargo bench -- --baseline pre-fixes`) → **skipped** — no pre-fix baseline was captured (preflight step 6 missed); the skill forbids reverting to back-fill one. See §9.

### Unresolved high-priority findings
- None. All five Majors are Applied and validated. Two Minors deferred (Mi1, Mi4) are non-blocking maintainability follow-ups.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | SSR write-then-can't-read (exclusive `last_pos`) | Apply | Applied | No | `registry_ssr.rs` | Pass (+test) | No |
| M2 | Major | Missing `sum(n_spanning)==n_total_alleles` | Apply | Applied | No | `registry_ssr.rs`, `errors.rs` | Pass | test=follow-up |
| M3 | Major | SSR structural errors mislabeled as `Io` | Apply | Applied | No | `registry_ssr.rs`, `errors.rs` | Pass | test=follow-up |
| M4 | Major | `records_of` no kind guard | Ask→Apply | Applied | Yes (poison) | `reader.rs`, `errors.rs` | Pass (+test) | No |
| M5 | Major | `from_tag` weakens M4 compile-time guarantee | Ask→Apply | Applied | Yes (macro) | `kind.rs`, `registry.rs`, `registry_ssr.rs` | Pass | No |
| Mi1 | Minor | `n_total_alleles` SNP-name on generic field | Defer | Deferred | No | None | N/A | Yes |
| Mi2 | Minor | Lost eager block-buffer drop | Apply | Applied | No | `kind.rs`, `reader.rs`, `registry_ssr.rs` | Pass | No |
| Mi3 | Minor | Missing `kind` silently defaults to snp | Ask→Apply | Applied | Yes (require) | `header.rs` | Pass (+test) | No |
| Mi4 | Minor | `columns_for_kind` module placement | Defer | Deferred | No | None | N/A | Yes |
| Mi5 | Minor | `pub mod registry_ssr` over-exposed | Apply | Applied | No | `mod.rs` | Pass | No |
| Mi6 | Minor | Dead `AllelesLessThanRecords` variant | Apply | Applied | No | `errors.rs` | Pass | No |
| Mi7 | Minor | `PosOutOfRange` wrong range for SSR `end` | Apply | Applied | No | `errors.rs`, `writer.rs` | Pass (+test) | No |
| Mi8 | Minor | NaN logliks round-trip silently | Apply | Applied | No | `errors.rs`, `writer.rs` | Pass (+test) | No |
| Mi9 | Minor | `span[i] as u32` unchecked truncation | Apply | Applied | No | `registry_ssr.rs` | Pass | test=follow-up |
| Mi10 | Minor | Stale module docs (kind.rs, mod.rs) | Apply | Applied | No | `kind.rs`, `mod.rs` | Pass (doc) | No |
| Mi11 | Minor | `INITIAL_*_HINT` magic capacities | Apply | Applied | No | `registry_ssr.rs` | Pass (doc) | No |
| Nits | Nit | markers / glossary / dead doc | Partial | Applied (4) / Deferred (4) | No | `writer.rs`, `reader.rs`, `registry_ssr.rs` | Pass | Yes |

## 3. Questions asked and answers

1. **M4** — How should `records_of::<S>()` guard against a schema/kind mismatch?
   - **Answer:** *Poison iterator on mismatch* — keep the infallible signature; the iterator yields `KindMismatch` on its first `next()`, then ends.
2. **M5** — How should the `from_tag` exhaustiveness gap be closed?
   - **Answer:** *`column_key!` macro (single source)* — generate the enum + `tag()` + `from_tag()` from one variant↔tag list, restoring a compile-time bijection.
3. **Mi3** — What should happen for a header missing the `kind` tag?
   - **Answer:** *Require kind (reject missing)* — drop the serde default; a tag-less header is a hard parse error.

## 4. Per-finding log

### M1 — SSR writer emits a file its own reader rejects (exclusive `last_pos`)
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Correctness bug, High confidence, one clearly-correct path (the review-recommended inclusive store, which also matches the SNP point semantics and the index.rs doc). No policy invention.
- **Implementation summary:** `SsrBlock::append` now stores `self.last_pos = self.last_pos.max(record.end - 1)` (inclusive last-touched base) instead of the exclusive `record.end`. `end > start >= 1` is guaranteed by `validate_locus`, so `end - 1` never underflows.
- **Review suggestion used verbatim?:** Yes (the `max(end - 1)` form).
- **Verification performed:** Added `ssr_locus_at_contig_end_round_trips` (writes a locus with `end == chrom.length + 1`, reopens, asserts equality). Confirmed it exercises the defect: before the fix, `PspReader::new` rejected the file with `BlockIndexPosOutOfRange`.
- **Files changed:** `src/psp/registry_ssr.rs`
- **Tests added or modified:** `ssr_locus_at_contig_end_round_trips`
- **Validation:** `cargo test --lib psp` → 0, 204 passed. Full `cargo test --lib` → 0, 1142 passed.
- **User input:** None · **Follow-up:** None · **Residual risk:** None.

### M2 — Missing `sum(n_spanning) == n_total_alleles` check
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Robustness/correctness on untrusted input; the SNP path has the symmetric B3 check. One correct path. Uses the M3 typed variant.
- **Implementation summary:** `SsrDecoder::decode_block` now reconciles `self.n_spanning.iter().sum::<u64>()` against `header.n_total_alleles`, returning the new `PspReadError::SsrProfileCountMismatch { expected, got }` on disagreement (catches under-count = silent profile loss, and over-count). The `next_record` over-run guard is retained as defense-in-depth.
- **Review suggestion used verbatim?:** No — adapted: a dedicated `SsrProfileCountMismatch` variant (per M3) instead of the review's placeholder reuse of `AllelesLessThanRecords`.
- **Verification performed:** Logic verified by reading the decode path; the existing `ssr_round_trips_through_the_container` (which includes a zero-spanning locus) still passes, confirming the check accepts well-formed blocks where `sum == n_total_alleles`.
- **Files changed:** `src/psp/registry_ssr.rs`, `src/psp/errors.rs`
- **Tests added or modified:** None new (see Follow-up).
- **Validation:** `cargo test --lib` → 0, 1142 passed; `cargo clippy … -D warnings` → 0.
- **User input:** None
- **Follow-up:** A regression test asserting `SsrProfileCountMismatch` needs a hand-built malformed block (the writer cannot produce an inconsistent block via the public API). Tracked as a follow-up test fixture.
- **Residual risk:** The guard is unit-tested only indirectly (well-formed acceptance); the rejection path awaits a raw-block fixture.

### M3 — SSR structural errors mislabeled as `PspReadError::Io`
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Error-modeling defect; `PspReadError` is `#[non_exhaustive]`, so adding variants is non-breaking. One correct path.
- **Implementation summary:** Added `PspReadError::BlockStructureInvalid { context: &'static str }`. Routed both SSR sites through it — the CSR-offset-disagreement check in `decode_block` and the `profile_end >= lengths_offsets.len()` over-run guard in `next_record` — replacing the synthesized `Io { source: io::Error::other(...) }`.
- **Review suggestion used verbatim?:** No — adapted: one `BlockStructureInvalid` variant for both structural sites (plus `SsrProfileCountMismatch` for the count axis), rather than the review's two finer-grained names.
- **Verification performed:** Build + clippy clean; the offset-disagreement message is preserved verbatim in `context`.
- **Files changed:** `src/psp/registry_ssr.rs`, `src/psp/errors.rs`
- **Tests added or modified:** None new (see Follow-up).
- **Validation:** `cargo test --lib` → 0, 1142 passed.
- **User input:** None
- **Follow-up:** A `decode_block_returns_error_on_csr_offset_disagreement` test needs a hand-built block (parallel CSR columns with mismatched row boundaries). Tracked as a follow-up fixture.
- **Residual risk:** None for SNP; SSR rejection paths await the raw-block fixture.

### M4 — `records_of::<S>()` no kind guard (silent cross-schema misdecode)
- **Severity:** Major · **Initial decision:** Ask → Apply · **Final status:** Applied
- **Reasoning:** User chose poison-the-iterator over a fallible signature, keeping `records_of`/`records`/`region_records` infallible.
- **Implementation summary:** Added `PspReadError::KindMismatch { expected: &'static str, found: String }`. `RecordsIter::new` sets a `pending_error: Option<PspReadError>` when `S::KIND != reader.header().kind`; `Iterator::next` surfaces it once (then poisons). Works for all entry points — `records()`/`region_records()` (default `SnpKind`) never flag an snp file, but an snp decoder over an ssr file (or vice-versa) now fails loudly.
- **Review suggestion used verbatim?:** No — the user-chosen poison mechanism instead of the review's `Result` return.
- **Verification performed:** Added `records_of_wrong_kind_yields_kind_mismatch_then_poisons` (reads an ssr file via the default SNP `records()`, asserts `KindMismatch { expected: "snp", found: "ssr" }` then `None`).
- **Files changed:** `src/psp/reader.rs`, `src/psp/errors.rs`
- **Tests added or modified:** `records_of_wrong_kind_yields_kind_mismatch_then_poisons`
- **Validation:** `cargo test --lib psp` → 0, 204 passed; full suite 1142 passed.
- **User input:** Poison-on-mismatch (see §3).
- **Follow-up:** None. (Generalizing `region_records` over `S` remains the separate out-of-scope item; `region_records` on a wrong-kind file is now also poisoned.)
- **Residual risk:** None.

### M5 — `from_tag` array weakens the M4 compile-time guarantee
- **Severity:** Major · **Initial decision:** Ask → Apply · **Final status:** Applied
- **Reasoning:** User chose the `column_key!` macro, which restores a true compile-time single source for both `ColumnKey` and `SsrColumnKey`.
- **Implementation summary:** Added a `column_key!` `macro_rules!` in `kind.rs` (`pub(crate) use`) that generates the enum + `pub const fn tag()` + `pub fn from_tag()` as two exhaustive `match`es from one `Variant = 0xNN` list. Replaced both hand-written enums + `tag()`/`from_tag()` impls with macro invocations. A forgotten variant now cannot compile (it must appear in the list to be emitted at all). The old runtime drift-guard tests (`column_keys_are_unique_and_cover_v1_0`, `ssr_columns_are_well_formed`) are kept as a check that `V1_0_COLUMNS`/`SSR_COLUMNS` tags match the keys.
- **Review suggestion used verbatim?:** Adapted — the review's illustrative macro generated `from_tag` as a `match`; the implemented macro additionally attaches enum/variant doc-attrs (`$(#[$meta])*`) and constrains tags to `:literal` so they double as patterns.
- **Verification performed:** Build + clippy + full suite clean; the macro expansion is exercised by every existing SNP/SSR encode/decode test.
- **Files changed:** `src/psp/kind.rs`, `src/psp/registry.rs`, `src/psp/registry_ssr.rs`
- **Tests added or modified:** None new (existing key tests cover it).
- **Validation:** `cargo test --lib` → 0, 1142 passed; `cargo clippy … -D warnings` → 0.
- **User input:** `column_key!` macro (see §3).
- **Follow-up:** None · **Residual risk:** None.

### Mi1 — `BlockHeader.n_total_alleles` SNP-name on generic field
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** A pure but broad cross-file rename (`block.rs` wire struct + `writer.rs` flush + `reader.rs` SNP decode + `registry_ssr.rs`); the review itself offers deferral. Out of minimal-diff discipline for a correctness-focused pass, deferred to a dedicated rename commit.
- **Implementation summary:** None.
- **Files changed:** None · **Tests:** None · **Validation:** N/A
- **Follow-up:** Rename `n_total_alleles` → `n_entries` (varint-positional on disk, so behaviour-preserving). · **Residual risk:** Cosmetic only.

### Mi2 — Lost eager block-buffer drop (memory-residency regression)
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Restores the documented pre-refactor contract; ties to the project's RAM-for-scaling thesis. Small, clear.
- **Implementation summary:** Added `BlockDecoder::unload(&mut self)` to the trait. `SnpDecoder::unload` sets `cur_block = None` (frees the per-record/per-allele fixed-width columns; the H1 CSR reuse slabs are kept). `SsrDecoder::unload` clears the per-locus column vecs and sets `loaded = false` (CSR reuse slabs kept). `Iterator::next` calls `self.decoder.unload()` when a block is exhausted, before advancing. Also corrected the stale `materialise_next_record` doc (it claimed `mem::take` ownership-move + wholesale drop; the H1 rework reads by index/copy, and the drop is now via `unload`).
- **Review suggestion used verbatim?:** Yes (the `unload()` trait-method shape).
- **Verification performed:** Full suite (1142) passes — `next_record` is only reachable while `cur_block_loaded`/`loaded`, so freeing on exhaust is behaviour-neutral; the existing sequential/region read tests cover the exhaust+advance path.
- **Files changed:** `src/psp/kind.rs`, `src/psp/reader.rs`, `src/psp/registry_ssr.rs`
- **Tests added or modified:** None new (covered by existing multi-block read tests).
- **Validation:** `cargo test --lib` → 0, 1142 passed.
- **User input:** None · **Follow-up:** None · **Residual risk:** The memory delta is not asserted by a test (would need an RSS/live-heap probe); behaviour is correctness-neutral.

### Mi3 — Missing `kind` silently defaults to snp
- **Severity:** Minor · **Initial decision:** Ask → Apply · **Final status:** Applied
- **Reasoning:** User chose to make the tag mandatory (no on-disk back-compat requirement).
- **Implementation summary:** Removed `#[serde(default = "default_kind")]` from `WireHeader::kind` and deleted `fn default_kind`. A tag-less header now fails TOML deserialization (`PspReadError::HeaderToml`). Inverted the test `missing_kind_defaults_to_snp` → `missing_kind_is_rejected`.
- **Review suggestion used verbatim?:** Yes (drop the default; invert the test).
- **Verification performed:** `missing_kind_is_rejected` asserts `PspReadError::HeaderToml` for a header with the `kind` line stripped.
- **Files changed:** `src/psp/header.rs`
- **Tests added or modified:** `missing_kind_is_rejected` (replaces `missing_kind_defaults_to_snp`)
- **Validation:** `cargo test --lib psp` → 0, 204 passed.
- **User input:** Require kind (see §3).
- **Follow-up:** None · **Residual risk:** Pre-tag step-1a/1b `.psp` files (if any survive) are now unreadable — acceptable under the stated no-back-compat constraint.

### Mi4 — `columns_for_kind` module placement
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** A structural refactor (relocate `SNP_KIND`/`KNOWN_KINDS`/`columns_for_kind` into `kind.rs` and dispatch via the trait). Not required by any correctness fix; deferred to keep the diff focused. Mi5 (visibility) was applied independently.
- **Implementation summary:** None.
- **Files changed:** None · **Validation:** N/A
- **Follow-up:** Move the kind taxonomy to a neutral home and dispatch through `SnpKind::columns()`/`SsrKind::columns()` to drop the `registry ⇄ registry_ssr` edge. · **Residual risk:** None.

### Mi5 — `pub mod registry_ssr` over-exposed
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** One-line visibility tightening to match the sibling `registry` (which is `pub(crate)`); grep confirmed no consumer outside `src/psp/`.
- **Implementation summary:** `pub mod registry_ssr;` → `pub(crate) mod registry_ssr;`.
- **Review suggestion used verbatim?:** Yes.
- **Verification performed:** Build clean — all `registry_ssr` users (`writer.rs`, `registry.rs`, its own test) are in-crate.
- **Files changed:** `src/psp/mod.rs`
- **Validation:** `cargo build --lib` → 0; full suite 1142 passed.
- **User input:** None · **Follow-up:** None · **Residual risk:** None.

### Mi6 — Dead `AllelesLessThanRecords` variant
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Co-locate the retention rationale at the declaration (the review preferred either removal or a doc note; a doc note is the minimal, format-stable choice since the enum is part of the public error surface).
- **Implementation summary:** Added a doc comment on `BlockHeaderInvariantKind::AllelesLessThanRecords` explaining it is intentionally no longer raised by the generic validator and where SNP integrity now lives.
- **Review suggestion used verbatim?:** Yes (the doc-note option).
- **Files changed:** `src/psp/errors.rs` · **Validation:** `cargo doc --no-deps` → 0; full suite 1142.
- **User input:** None · **Follow-up:** None · **Residual risk:** None.

### Mi7 — `PosOutOfRange` wrong range / conflated failures for SSR `end`
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Diagnostic correctness; `#[non_exhaustive]` enum, non-breaking.
- **Implementation summary:** Added `PspWriteError::LocusEndOutOfRange { record_index, chrom_id, start, end, chrom_length }` with an honest `(start, length+1]` message; `validate_locus` now raises it (instead of `PosOutOfRange`) for a bad `end`.
- **Review suggestion used verbatim?:** Adapted — kept the single `end <= start || end > length+1` condition with the honest-range message (the review's "fold it into the message" option) rather than splitting into two arms.
- **Verification performed:** Added `write_locus_rejects_end_le_start_with_locus_end_variant` (asserts `LocusEndOutOfRange { start: 50, end: 50, .. }` for `end == start`).
- **Files changed:** `src/psp/errors.rs`, `src/psp/writer.rs`
- **Tests added or modified:** `write_locus_rejects_end_le_start_with_locus_end_variant`
- **Validation:** `cargo test --lib psp` → 0, 204 passed.
- **User input:** None · **Follow-up:** None · **Residual risk:** None.

### Mi8 — NaN logliks round-trip silently
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** NaN is never a valid probability; rejecting it at write is the conservative, clearly-correct guard (the review's "reject NaN, permit -inf" option). No design ambiguity.
- **Implementation summary:** Added `PspWriteError::NonFiniteLoglik { record_index, profile_index }`. `validate_locus` rejects any `loglik.is_nan()` while permitting `±inf` (legitimate `log(0)`).
- **Review suggestion used verbatim?:** Yes (writer-side NaN-only reject).
- **Verification performed:** Added `write_locus_rejects_nan_loglik_but_permits_neg_inf` (NaN → `NonFiniteLoglik { profile_index: 0 }`; `-inf` accepted).
- **Files changed:** `src/psp/errors.rs`, `src/psp/writer.rs`
- **Tests added or modified:** `write_locus_rejects_nan_loglik_but_permits_neg_inf`
- **Validation:** `cargo test --lib psp` → 0, 204 passed.
- **User input:** None · **Follow-up:** None · **Residual risk:** None.

### Mi9 — `span[i] as u32` unchecked truncation
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Make the narrowing symmetric with the adjacent `delta-start` `u32::try_from` guard.
- **Implementation summary:** `SsrDecoder::next_record` now does `u32::try_from(self.span[i])` → `ColumnElementDecode { column: "span", source: VarintOverflow }` on overflow, instead of `as u32`.
- **Review suggestion used verbatim?:** Yes.
- **Verification performed:** Build + full suite clean; the well-formed round-trip tests still pass (spans are small).
- **Files changed:** `src/psp/registry_ssr.rs`
- **Tests added or modified:** None new (see Follow-up).
- **Validation:** `cargo test --lib` → 0, 1142 passed.
- **User input:** None
- **Follow-up:** The overflow path needs a hand-built block (the writer caps span ≤ contig length); tracked with the M2/M3 fixture work.
- **Residual risk:** Defensive guard, unit-tested only via well-formed acceptance.

### Mi10 — Stale module docs (kind.rs, mod.rs)
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Doc-only accuracy.
- **Implementation summary:** `kind.rs` module doc no longer says "writer-side surface only" (it now mentions the read-side `BlockDecoder` mirror). `mod.rs` reframed as a "generic `.psp` columnar container" hosting `snp`+`ssr`, qualified `registry` as the SNP registry, and added `kind`/`registry_ssr` to the layout.
- **Files changed:** `src/psp/kind.rs`, `src/psp/mod.rs` · **Validation:** `cargo doc --no-deps` → 0.
- **User input:** None · **Follow-up:** None · **Residual risk:** None.

### Mi11 — `INITIAL_*_HINT` magic capacities
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Doc-only — anchor the numbers.
- **Implementation summary:** Documented the basis (≈one block-window of loci × ~16 profiles × ~2 candidates) and why they're fixed rather than block-size-derived, with inline annotations on each const.
- **Files changed:** `src/psp/registry_ssr.rs` · **Validation:** build clean.
- **User input:** None · **Follow-up:** None · **Residual risk:** None.

### Nits
- **Applied (4):** `// UNREACHABLE:` markers on the four `from_tag(...).expect()` sites (writer.rs, reader.rs, registry_ssr.rs ×2); fixed the stale `key_ssr` doc in `ssr_columns_are_well_formed`; added an `amb`/`frr` glossary to the `registry_ssr` module doc. (The `from_tag` array-scan / dead-code-expect nits are subsumed by M5's macro and Mi6's note.)
- **Deferred (4):** named `RecordInterval` struct for `record_interval`; `KNOWN_KINDS` derived-from-a-table; `LadderEntry`/`Profile` named type for `spanning`; named consts for `SsrBlock::append`'s `projected_bytes` arithmetic. All optional readability changes with no behavioural impact.
- **Files changed:** `src/psp/writer.rs`, `src/psp/reader.rs`, `src/psp/registry_ssr.rs` · **Validation:** clippy/doc/test clean.

## 5. Deferred findings to carry forward
- **Mi1** — rename `BlockHeader.n_total_alleles` → `n_entries` (broad mechanical rename).
- **Mi4** — relocate the kind taxonomy (`columns_for_kind`) out of the SNP registry into a neutral home.
- **Deferred Nits** — `RecordInterval` struct, `KNOWN_KINDS` table, `LadderEntry`/`Profile` type, `projected_bytes` named consts.
- **Follow-up tests** (need hand-built malformed-block fixtures, not producible via the writer): `decode_block` rejection paths for M2 (`SsrProfileCountMismatch`), M3 (`BlockStructureInvalid` on offset disagreement), and Mi9 (`span` overflow). The fixes themselves are in and the well-formed acceptance paths are tested; only the rejection assertions await a raw-block harness.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** Yes — `src/psp/reader.rs` (changed by Mi2/M4/M5) is reachable from `benches/psp_reader_perf.rs`.
- **Baseline saved:** No — editing began before a criterion baseline was captured (preflight step 6 was missed). The skill forbids reverting/un-editing to back-fill a baseline, so the comparison is skipped rather than fabricated.
- **Benches run:** None.
- **Verdicts:** N/A.
- **Outcome:** Skipped — no pre-fix baseline.
- **Notes:** The hot-path changes are micro and neutral-to-faster: (a) one `Option::take()` per `next()` (None after the first call), (b) one `unload()` (a field-to-`None` write) per exhausted block, (c) `from_tag` rewritten as a `match` instead of a 12-element `.into_iter().find()` array scan (per column per block — if anything cheaper). No algorithmic change to the per-record decode. A baseline-anchored comparison is a recommended follow-up if the SSR/cohort read path is profiled.

## 10. Commands run
- `./scripts/dev.sh cargo build --lib`
- `./scripts/dev.sh cargo test --lib psp`
- `./scripts/dev.sh cargo fmt`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo doc --no-deps`

## 11. Command results
- `cargo build --lib` → 0, `Finished` (clean).
- `cargo test --lib psp` → 0, `204 passed; 0 failed`.
- `cargo fmt --check` → 0, clean (after one `cargo fmt`).
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, `Finished … in 3.49s` (0 warnings).
- `cargo test --lib` → 0, `1142 passed; 0 failed; 1 ignored`.
- `cargo doc --no-deps` → 0, `Finished` (no warnings).
- `cargo audit` → not run (not installed in the container).

## 12. Notes
- All five Majors are Applied and validated; the SNP byte-identity / behaviour claims confirmed in the review are unaffected (no SNP semantics changed — M5's macro reproduces the same tag↔key bijection; Mi2's `unload` is behaviour-neutral; the new SNP-path code is the `pending_error` kind check, inert for snp-on-snp).
- The three malformed-input rejection tests (M2/M3/Mi9) are the only deliberate test gap; they require a raw-block fixture the writer can't emit, so the rejection assertions are deferred while the fixes and well-formed acceptance are in place.
