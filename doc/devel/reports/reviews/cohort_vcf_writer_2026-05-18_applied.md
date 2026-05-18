# Fix Application Report: cohort_vcf_writer_2026-05-18

**Date:** 2026-05-18
**Source review:** [cohort_vcf_writer_2026-05-18.md](cohort_vcf_writer_2026-05-18.md)
**Source state reviewed against:** `review/vcf-writer` @ `9cf0010`
**Execution mode:** interactive (4 open questions resolved by user before editing began)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 15
- Minors: 19
- Nits: grouped (~6 items)

### Outcome totals
- Applied: 1 Blocker + 14 Major + 16 Minor = **31** (plus the grouped Nits cluster, applied alongside the Minors)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 4 (M15 + Mi11 + Mi12 + Mi14)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --all -- --check` → exit non-zero, **but** every in-scope file (`src/var_calling/vcf_writer/*`, `tests/cohort_vcf_writer_integration.rs`) passes; the 19 remaining hunks are all in parallel SIMD WIP files (`src/var_calling/posterior_engine/*`, `benches/var_calling_perf.rs`, `examples/profile_posterior_engine.rs`) on a different branch and out of scope.
- `cargo clippy --lib --tests --all-features -- -D warnings` → exit 0, clean.
- `cargo test --all-targets` → exit 0. Lib: **778 passed / 0 failed** (was 756 before the run — +22 new tests). Integration crates: 3+4+5+25+26+17+2+7+4+8+12 = 113 passed / 0 failed.
- `cargo doc --no-deps` → **not run.** The parallel `posterior_engine` WIP has broken intra-doc links unrelated to vcf_writer that would dominate the output. Re-run once the parallel branch lands. The fixes in scope did not add or change intra-doc links; the rustdoc burden lands on the next implementer who touches `posterior_engine`.
- `cargo audit` → **not run** (not part of this repo's standard pre-merge gate).
- Performance check (`cargo bench -- --baseline pre-fixes`) → **not applicable** — no `vcf_writer` bench exists (this is itself M15, deferred).

### Unresolved high-priority findings
- **M15** — deferred. Adding the criterion bench is its own slice; the writer ships now without a baseline. Listed in §5.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | `tally_called_alleles` silently drops out-of-bounds ALT allele indices | Apply | Applied | No | `errors.rs`, `record_encode.rs` | unit test + full suite | No |
| M1 | Major | `Io(#[from] io::Error)` collapses 8 distinct call sites | Apply | Applied | No | `errors.rs`, `sink.rs`, `writer.rs` | full suite | No |
| M2 | Major | `Encode(String)` is a catch-all stringly-typed variant | Apply | Applied | Q2 → box | `errors.rs`, `header.rs`, `record_encode.rs` | full suite | No |
| M3 | Major | public `VcfWriteError` not `#[non_exhaustive]` | Apply | Applied | Q4 → yes | `errors.rs` | full suite | No |
| M4 | Major | `Display` messages prepend mechanism + flatten source | Apply | Applied | No | `errors.rs` (folded into M1/M2 rewrite) | full suite | No |
| M5 | Major | `SampleCountMismatch` message inverts field semantics | Apply | Applied | No | `errors.rs`, `record_encode.rs` (test) | unit test + full suite | No |
| M6 | Major | Multiple unchecked slice indexes panic on malformed `PosteriorRecord` | Apply | Applied | No | `errors.rs`, `record_encode.rs` | unit tests + full suite | No |
| M7 | Major | DP silent saturation `unwrap_or(i32::MAX)` | Apply | Applied | No | `errors.rs`, `record_encode.rs` | unit test + full suite | No |
| M8 | Major | `num_obs as i32` wrapping cast | Apply | Applied | No | `record_encode.rs` | unit test + full suite | No |
| M9 | Major | GQ cap `99.0` magic + drift vs engine cap | Apply | Applied | No | `record_encode.rs` | drift test + full suite | No |
| M10 | Major | `finish` doesn't fsync the parent directory | Apply | Applied | Q3 → yes | `errors.rs`, `sink.rs` | smoke test + full suite | No |
| M11 | Major | `Default for WriterConfig` empty path + behavioural defaults | Apply | Applied | Q1 → drop | `mod.rs` (+ test updates in 5 files) | full suite | No |
| M12 | Major | `CohortVcfWriter` config not inspectable | Apply | Applied | No | `writer.rs` | accessor test + full suite | No |
| M13 | Major | `<output>.tmp` leak on header-write failure during `new()` | Apply | Applied | No | `writer.rs` | covered by error type + path; explicit test deferred (needs sink fault-injection seam) | follow-up: add fault-injection test once a sink seam exists |
| M14 | Major | `CohortVcfWriter` missing `#[must_use]` | Apply | Applied | No | `writer.rs` | clippy + full suite | No |
| M15 | Major | No criterion bench for `vcf_writer` despite hot-path designation | Defer | Deferred | No | None | N/A | separate slice |
| Mi1 | Minor | Missing test that `last_locus` doesn't advance on error | Apply | Applied | No | `writer.rs` (test) | unit test | No |
| Mi2 | Minor | `path_is_bgzf` case-sensitive, undocumented | Apply | Applied | No | `sink.rs` | unit test | No |
| Mi3 | Minor | Empty `contigs` silently accepted, neither documented nor tested | Apply | Applied | No | `header.rs` (doc + test) | unit test | No |
| Mi4 | Minor | Duplicated overflow checks; `ContigLengthOverflow` no test | Apply | Applied | No | `header.rs` | unit test | No |
| Mi5 | Minor | No test for empty `tool_string` / empty `command_line` skip paths | Apply | Applied | No | `header.rs` (tests) | 2 unit tests | No |
| Mi6 | Minor | Duplicate `##source` / `##commandline` insert blocks | Apply | Applied | No | `header.rs` (extracted `insert_unstructured`) | full suite | No |
| Mi7 | Minor | Duplicated genotype-table lookup with same 5-field error | Apply | Applied | No | `record_encode.rs` (extracted `lookup_genotype`) | full suite | No |
| Mi8 | Minor | `gt_buf` clone+reuse defeats amortisation | Apply | Applied | No | `record_encode.rs` | full suite | No |
| Mi9 | Minor | `metadata.contigs.clone()` after consume | Apply | Applied | No | `writer.rs` | full suite | No |
| Mi10 | Minor | `..Default::default()` in one test breaks the compile-fence | Apply | Applied | No | `record_encode.rs` (folded into M11 — `Default` is gone so this is moot) | full suite | No |
| Mi11 | Minor | `WriterConfig` should be `CohortVcfWriterConfig` | Defer | Deferred | No | None | N/A | public-API rename — pair with Mi12 in a coordinated naming pass |
| Mi12 | Minor | `tool_string` field name encodes type, not domain | Defer | Deferred | No | None | N/A | public-API rename — pair with Mi11 |
| Mi13 | Minor | `ref_allele` used for ALT slots in integration test | Apply | Applied | No | `tests/cohort_vcf_writer_integration.rs` | full integration suite | No |
| Mi14 | Minor | Three-way fixture duplication across test modules | Defer | Deferred | No | None | N/A | needs shared `test_fixtures` module — separate consolidation slice |
| Mi15 | Minor | Missing `# Errors` rustdoc sections on `pub fn`s | Apply | Applied | No | `writer.rs` | (doc-only) |  No |
| Mi16 | Minor | BGZF EOF magic duplicated | Apply | Applied (adaptation) | No | `sink.rs` (BGZF_EOF promoted to `pub(crate)` const; in-crate test deduplicated; integration test keeps its copy because a separate test crate can't reach `pub(crate)`) | sink unit test | follow-up: add a `pub` re-export or test-support feature flag if the duplication becomes a problem |
| Mi17 | Minor | `genotype_order` rebuilt per record | Apply | Applied | No | `writer.rs` (HashMap cache keyed by `(ploidy, n_alleles)`) | cache-shape test + full suite | No |
| Mi18 | Minor | Missing malformed-input tests for `Encode` (UTF-8) + `GenotypeIndexOutOfBounds` | Apply | Applied | No | `record_encode.rs` (2 new tests) | unit tests | No |
| Mi19 | Minor | Plan rustdoc link is `https://example.invalid` | Apply | Applied | No | `mod.rs` | (doc-only) | No |
| Nits | Nit | Grouped: `let _ = write!`, `Vec<String>` for `Keys`, `alt_allele` aliases, closure type annotation, `matches!` bare `..`, variant doc-comment coverage | Apply | Applied (bundled with Minor cleanups) | No | `record_encode.rs`, `header.rs`, `sink.rs`, `tests/...` | full suite | No |

## 3. Questions asked and answers

Four open questions from §4 of the review were posed up-front and answered before any code was edited:

1. **Q1 — Should `Default for WriterConfig` exist at all?**
   - **Answer:** Drop `impl Default`. Force explicit `output` via a constructor. Affects M11, M12, Mi11.

2. **Q2 — Should `VcfWriteError::Encode` expose noodles types in its `#[source]` chain, or stay opaque via `Box<dyn Error>`?**
   - **Answer:** Box behind `Box<dyn Error + Send + Sync>`. The crate is destined for crates.io publication; insulate the public API from noodles' fast semver cadence. Affects M2, M4.

3. **Q3 — Is a parent-directory fsync required for the documented durability contract?**
   - **Answer:** Yes — include the parent-dir fsync in `finish`. Workflow-manager use cases depend on "after `finish()` Ok, the file exists across crash". Affects M10.

4. **Q4 — Add `#[non_exhaustive]` to `VcfWriteError`?**
   - **Answer:** Yes. Zero runtime cost, future variants are non-breaking minor bumps once the crate ships. Affects M3.

## 4. Per-finding log

### B1 — `tally_called_alleles` silently drops out-of-bounds ALT allele indices
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Convergent finding across reliability + 5 cross-cat notes; the silent shape produces VCFs that parse cleanly but violate the `sum(AC) + REF == AN` invariant.
- **Implementation summary:** Added typed `VcfWriteError::AlleleIndexOutOfBounds { chrom_id, pos, sample_idx, allele_idx, n_alleles }` variant. Replaced the `if alt_pos < ac_per_alt.len()` silent guard in `tally_called_alleles` with an explicit out-of-range check that returns the new variant. Reordered the bounds check to happen *before* `an_total += 1` so the AN tally is consistent with what was tallied (i.e. we error before incrementing, so the partial tally is discarded with the error).
- **Review suggestion used verbatim?:** No — adapted: the suggestion's check was after `an_total += 1`; we moved it earlier for tally consistency.
- **Verification performed:** Unit test `tally_called_alleles_errors_on_out_of_bounds_allele_index` drives a synthetic record whose `genotype_order(2, 3)[3] == [0, 2]` references allele 2 against a 2-allele record; asserts the typed error.
- **Files changed:** `src/var_calling/vcf_writer/errors.rs`, `src/var_calling/vcf_writer/record_encode.rs`.
- **Tests added or modified:** `tally_called_alleles_errors_on_out_of_bounds_allele_index`.
- **Validation:** `cargo test --lib var_calling::vcf_writer` → 41 passed.

### M1 — `Io(#[from] io::Error)` collapses 8 distinct call sites
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Operators tracing a failure can't tell `create-tmp` from `rename` from `fsync` etc.
- **Implementation summary:** Removed `#[from] io::Error`. Added per-operation variants — `CreateTmp`, `WriteHeader`, `WriteRecord`, `FinishBgzf`, `FsyncFile`, `FsyncDir`, `Rename` — each carrying the relevant path(s) and `#[source] io::Error`. Updated every `?` over an `io::Error` in `sink.rs` and `writer.rs` to attach the matching variant via `.map_err`.
- **Review suggestion used verbatim?:** No — adapted: added `FsyncDir` (a new operation introduced by M10), used `WriteHeader` for the `BufWriter::into_inner` flush-on-finish path (closest semantic match; documented inline).
- **Files changed:** `errors.rs`, `sink.rs`, `writer.rs`.
- **Validation:** `cargo test --all-targets` → all green.

### M2 — `Encode(String)` is a catch-all stringly-typed variant
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Per user's Q2 answer, box behind `Box<dyn Error + Send + Sync>` to insulate the public API from noodles semver.
- **Implementation summary:** New shape: `Encode { operation: &'static str, source: Box<dyn std::error::Error + Send + Sync> }`. The `operation` tag is a project-side `&'static str` (e.g. `"##source key parse"`, `"allele bytes UTF-8"`) so we keep operation context. All six call sites updated to `.map_err(|e| VcfWriteError::Encode { operation: "…", source: Box::new(e) })`.
- **Review suggestion used verbatim?:** No — adapted: the suggested form was `Encode(#[source] Box<…>)` (single positional); we used a struct variant with an explicit `operation` tag because losing the operation context would have made the error harder to triage.
- **Files changed:** `errors.rs`, `header.rs`, `record_encode.rs`.

### M3 — `VcfWriteError` not `#[non_exhaustive]`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Per user's Q4 answer, future variants should be non-breaking.
- **Implementation summary:** Added `#[non_exhaustive]` on the enum. No call-site updates needed because all in-crate `match`es include either `..` patterns or are exhaustive against variants the writer guarantees to use.
- **Files changed:** `errors.rs`.

### M4 — `Display` messages prepend mechanism + flatten source
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Tools that walk the error chain print the cause twice with mechanism-prefix `Display`.
- **Implementation summary:** Folded into the M1/M2 redesign: every variant's `#[error("…")]` now describes the operation that failed, without interpolating the source. `Display` no longer contains `{0}` or `{source}` for variants that carry a `#[source]`; the cause is reached via `source()`.
- **Files changed:** `errors.rs` (folded into the same rewrite as M1 + M2).

### M5 — `SampleCountMismatch` message inverts the meaning of its two count fields
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Wrong-direction diagnostic actively misleads operators.
- **Implementation summary:** Swapped the message wording to "cohort metadata names X samples but the posterior arrays carry Y", matching the call-site population. Added a wording-lock test.
- **Files changed:** `errors.rs`, `record_encode.rs` (test).
- **Tests added:** `sample_count_mismatch_message_names_cohort_first`.

### M6 — Multiple unchecked slice indexes panic on malformed `PosteriorRecord`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Six unchecked indexes panic on a malformed upstream; same pattern as the existing `SampleCountMismatch`/`UnknownChromId` guards but missing.
- **Implementation summary:** New `VcfWriteError::InconsistentRecord { chrom_id, pos, field, expected, actual }` variant. New `validate_record_shape(record, expected_samples, table)` helper at the top of `encode` that checks every dependent vector length (`best_genotype`, `gq_phred`, `allele_frequencies`, `compound_frequencies`, `scalars`, `posteriors`, `chain_anchor_flags`) plus the `alleles` non-empty guard and the table-length consistency. Single early return on any mismatch.
- **Review suggestion used verbatim?:** No — adapted: added `compound_frequencies` and the table-length check (omitted from the suggested shape).
- **Files changed:** `errors.rs`, `record_encode.rs`.
- **Tests added:** `encode_errors_on_empty_alleles`, `encode_errors_on_undersized_best_genotype`.

### M7 — DP silent saturation `unwrap_or(i32::MAX)`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Silent saturation is data corruption.
- **Implementation summary:** New `VcfWriteError::DepthOverflow { chrom_id, pos, sample_idx: Option<usize>, depth: u64 }` variant. Replaced `unwrap_or(i32::MAX)` at the cohort-DP-total site with a typed `map_err`. (Per-sample DP fixed in M8's pass.)
- **Files changed:** `errors.rs`, `record_encode.rs`.
- **Tests added:** `encode_errors_when_per_sample_depth_overflows_i32`.

### M8 — `num_obs as i32` wrapping cast
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** A single allele's `num_obs > i32::MAX` would wrap to a negative AD silently.
- **Implementation summary:** Replaced the wrapping `as i32` with `i32::try_from(...).map_err(|_| DepthOverflow { sample_idx: Some(s), depth: u64::from(num_obs), .. })?` inside the AD-cell loop. Threads the same `DepthOverflow` variant introduced by M7. Also the per-sample DP-total cast is now `try_from`-with-error along the same shape.
- **Files changed:** `record_encode.rs`.
- **Tests added:** `encode_errors_when_single_allele_depth_overflows_i32`.

### M9 — GQ cap `99.0` magic + drift vs engine cap
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Two independent ceilings on the same field with no shared symbol.
- **Implementation summary:** Added `pub(super) const GQ_MAX: f32 = 99.0;` next to `QUAL_MAX` with the same doc style. Clamp updated to `clamp(0.0, GQ_MAX as f64)`. Drift test `gq_writer_cap_is_at_least_engine_cap` asserts `GQ_MAX >= PosteriorEngineConfig::default().max_gq_phred`.
- **Files changed:** `record_encode.rs`.
- **Tests added:** `gq_writer_cap_is_at_least_engine_cap`.

### M10 — `finish` doesn't fsync the parent directory after the atomic rename
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Per user's Q3 answer, the durability contract requires it for workflow-manager safety.
- **Implementation summary:** Added `sync_parent_dir(final_path)` helper that opens the parent (defaulting to `.` when `final_path.parent()` is `None`/empty) and calls `sync_all`. Called from `SinkKind::finish` immediately after `fs::rename`. Errors surface as the typed `VcfWriteError::FsyncDir { final_path, source }` (new variant introduced in M1).
- **Files changed:** `errors.rs`, `sink.rs`.
- **Tests added:** `finish_renames_and_fsyncs_parent_dir` (happy-path smoke that confirms the new code is exercised; the crash-recovery property is covered by the existence of the call).

### M11 — `Default for WriterConfig` empty path + behavioural defaults
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Per user's Q1 answer, drop the `Default` impl entirely.
- **Implementation summary:** Removed `impl Default for WriterConfig`. Added `pub const DEFAULT_FILTER_PASS: bool = true;` and `pub const DEFAULT_EMIT_GP: bool = false;` at module scope. Added `WriterConfig::new(output: PathBuf) -> Self` constructor that fills the two booleans from the consts. Updated every `WriterConfig::default()` call site (9 in `header.rs`, `record_encode.rs`, `writer.rs`) to either `WriterConfig::new(...)` or an explicit struct literal.
- **Files changed:** `mod.rs`, `header.rs` (tests), `record_encode.rs` (tests), `writer.rs` (tests), `tests/cohort_vcf_writer_integration.rs`.

### M12 — `CohortVcfWriter` config not inspectable
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `pub fn config(&self) -> &WriterConfig` accessor. Also added a manual `impl Debug for CohortVcfWriter` with exhaustive destructure (so a new field surfaces at compile time) and `finish_non_exhaustive` (so the noodles writer / header internals stay hidden).
- **Files changed:** `writer.rs`.
- **Tests added:** `config_accessor_exposes_frozen_settings`.

### M13 — `<output>.tmp` leak on header-write failure during `new()`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Constructor failures need to clean up the tmp file so they're distinguishable from "user forgot to call finish".
- **Implementation summary:** In `CohortVcfWriter::new`, wrap the `inner.write_header(&header)` call in an `if let Err(...)` that calls `std::fs::remove_file(<tmp>)` (best-effort) before returning the typed error. Surfaces as `VcfWriteError::WriteHeader { tmp_path, source }`.
- **Review suggestion used verbatim?:** No — adapted: did not add the fault-injection regression test because the writer has no sink-injection seam today; would require either a `#[cfg(test)]` private constructor or a `FaultyWriter` mock. **Follow-up:** add this test when the sink-injection seam exists (likely as part of the CLI slice, where the writer is wired up against a `Box<dyn Write>` sink).
- **Files changed:** `writer.rs`.

### M14 — `CohortVcfWriter` missing `#[must_use]`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `#[must_use = "CohortVcfWriter must be finalised by calling \`.finish()\`; …"]` on the struct. Also added `#[must_use]` on `config()` (read-only accessor; the rule applies).
- **Files changed:** `writer.rs`.

### M15 — No criterion bench for `vcf_writer` despite hot-path designation
- **Severity:** Major
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** Adding a `criterion` bench is its own slice — fixture shape, sample-size sweep, allocator choice, regression-threshold setup are decisions that fit poorly inside a fixes-applied run. The writer's hot-path concern is real (per-record `Vec<Vec<Option<SampleValue>>>` allocations, format-key reuse) but a single bench addition is more thoughtful than a quick add.
- **Follow-up:** A separate slice should add `benches/vcf_writer_perf.rs` with the three measurement points the review suggests (write_record throughput at 1000 samples × biallelic SNPs, same with `emit_gp = true`, encode-only).

### Mi1 — Missing test that `last_locus` doesn't advance on error
- Applied. Test `out_of_order_does_not_advance_last_locus` in `writer.rs`.

### Mi2 — `path_is_bgzf` case-sensitive
- Applied. `path_is_bgzf` now lower-cases the file-name before suffix matching. Doc-comment on the function updated. Test `uppercase_suffixes_select_bgzf_sink`.

### Mi3 — Empty `contigs` silently accepted
- Applied (documented direction). `CohortMetadata::contigs` doc-comment now says "empty vector is accepted; operationally meaningless but structurally valid". Test `empty_contigs_accepted_produces_header_with_no_contig_lines`.

### Mi4 — Duplicated overflow checks; `ContigLengthOverflow` no test
- Applied. Dropped the dead `usize::try_from(u32)` arm; only the `> i32::MAX as u32` check remains. Test `contig_length_above_i32_max_rejected`.

### Mi5 — No test for empty `tool_string` / empty `command_line` skip paths
- Applied. `header_omits_source_when_tool_string_empty` and `header_omits_commandline_when_command_line_empty`.

### Mi6 — Duplicate `##source` / `##commandline` insert blocks
- Applied. Extracted `insert_unstructured(builder, key, value)` helper. Tagged via a tiny `header_op_for_key` mapping so the `Encode { operation }` tag stays a `&'static str`.

### Mi7 — Duplicated genotype-table lookup with same 5-field error
- Applied. Extracted `lookup_genotype(table, record, sample_idx, gt_idx) -> Result<&[u8], _>` helper. Both `tally_called_alleles` and `build_samples` now call it.

### Mi8 — `gt_buf` clone+reuse defeats amortisation
- Applied. Dropped the `gt_buf.clear()` + `.clone()` pattern. Now allocates a fresh `String::with_capacity(2 * gt.len())` per iteration and moves it into the `SampleValue::String` — same allocation count, clearer code.

### Mi9 — `metadata.contigs.clone()` after consume
- Applied. Moved `metadata.contigs` instead of cloning.

### Mi10 — `..Default::default()` in `format_keys_track_emit_gp` test
- Applied. Folded into M11: `Default` is gone, so the test uses an explicit struct literal via the `cfg_emit_gp_on()` helper.

### Mi11 — `WriterConfig` → `CohortVcfWriterConfig`
- Deferred. Public-API rename without strong forcing function; pair with Mi12 in a coordinated naming pass.

### Mi12 — `tool_string` → `source_label`/`tool_name`
- Deferred. Same as Mi11.

### Mi13 — `ref_allele` used for ALT slots
- Applied. Added `alt_allele` alias to `tests/cohort_vcf_writer_integration.rs`; updated four ALT-position call sites.

### Mi14 — Three-way fixture duplication
- Deferred. Needs a shared `pub(crate) test_fixtures` module (or `tests/common/mod.rs`) — bigger than the slot.

### Mi15 — Missing `# Errors` rustdoc sections
- Applied. `CohortVcfWriter::new`, `write_record`, and `finish` now each carry a `# Errors` block enumerating the variants they may return.

### Mi16 — BGZF EOF magic duplicated
- Applied (adaptation). Promoted `BGZF_EOF` to `pub(crate) const` in `sink.rs` so the in-crate test reuses it. Added `#[cfg_attr(not(test), allow(dead_code))]` so the constant doesn't trigger the dead-code warning in production builds. The integration test in `tests/` keeps its own copy because a separate test crate can't reach `pub(crate)`; doc-comment on the const explains the two-place truth. **Follow-up:** the full deduplication requires either a `pub` re-export (exposes wire-format magic in the public API — undesirable for a published crate) or a `test-support` feature flag; not done in this run.

### Mi17 — `genotype_order` rebuilt per record
- Applied. Added `genotype_tables: HashMap<(u8, usize), Vec<Vec<u8>>>` field on `CohortVcfWriter`. `write_record` now does `entry((ploidy, n_alleles)).or_insert_with(|| genotype_order(...))` then clones the table out to pass to `encode`. The `encode` signature gained a `table: &[Vec<u8>]` parameter (was internally calling `genotype_order` itself). Test `genotype_table_cache_is_keyed_by_ploidy_and_n_alleles` exercises the cache shape.
- Minor caveat: each record still clones the table out of the HashMap (lifetime of the entry vs. `encode` borrow). A future refactor could replace the clone with a `Cow<'_, [Vec<u8>]>` or compute the table inline once and reuse — but the per-record cost is one `Vec<Vec<u8>>` clone of a few hundred bytes, dwarfed by the noodles encode allocations. Worth revisiting if M15's bench shows it matters.

### Mi18 — Missing malformed-input tests
- Applied. `invalid_utf8_allele_surfaces_encode_error` + `out_of_range_best_genotype_surfaces_typed_error`. (`ContigLengthOverflow` covered by Mi4.)

### Mi19 — Placeholder URL
- Applied. Module doc now reads `Plan: doc/devel/implementation_plans/cohort_vcf_writer.md.` (no broken link, no `example.invalid`).

### Nits — grouped
- `let _ = write!`: replaced with `write!(out, "{a}").expect("write to String is infallible")` — same shape, names the discarded `fmt::Error` honestly.
- `Vec<String>` for `Keys`: left as-is — the intermediate allocation is one per writer instance, not a hot path; the proposed `iter::once`/`chain` form is readability-neutral.
- `alt_allele` aliases (`writer.rs:140`, `record_encode.rs:348`): retained — both functions encode intent at the call site (REF vs ALT); deleting them was the alternative the Nit offered.
- Closure type annotation (`header.rs:61-65,72-76`): replaced with turbofish `parse::<…::Other>()`.
- `matches!` bare `..`: left as-is — the test patterns survive field renames better with `..` than with `field: _` (the latter would force a test edit on every variant field addition).
- Variant doc-comment coverage: every variant on `VcfWriteError` now has a doc comment naming the failure mode.

## 5. Deferred findings to carry forward

- **M15** — no criterion bench for `vcf_writer`. Add `benches/vcf_writer_perf.rs` as a focused follow-up slice with the three measurement points the review suggests.
- **Mi11** — `WriterConfig` → `CohortVcfWriterConfig` rename. Pair with Mi12.
- **Mi12** — `tool_string` → `source_label`/`tool_name` rename. Pair with Mi11.
- **Mi14** — Three-way test-fixture deduplication. Needs a `pub(crate) test_fixtures` module or `tests/common/mod.rs` plus a small `test-support` feature gate to reach the integration test crate.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** No — no `vcf_writer` bench exists (this is M15, deferred).
- **Baseline saved:** N/A.
- **Benches run:** N/A.
- **Verdicts:** N/A.
- **Outcome:** Skipped — no Apply touched perf-sensitive code covered by an existing bench.
- **Notes:** Adding a `vcf_writer_perf` bench is itself M15.

## 10. Commands run

- `cargo build` (multiple incremental — final clean)
- `cargo test --lib var_calling::vcf_writer`
- `cargo test --test cohort_vcf_writer_integration`
- `cargo fmt --all -- --check`
- `cargo fmt` (in-scope files only; out-of-scope SIMD WIP files reverted via `git checkout HEAD --`)
- `cargo clippy --lib --tests --all-features -- -D warnings`
- `cargo test --all-targets`

## 11. Command results

- `cargo fmt --all -- --check` (final) → 19 hunks, **all in out-of-scope parallel-SIMD-WIP files**; in-scope vcf_writer files pass.
- `cargo clippy --lib --tests --all-features -- -D warnings` → exit 0, clean.
- `cargo test --lib var_calling::vcf_writer` → 41 passed / 0 failed.
- `cargo test --test cohort_vcf_writer_integration` → 4 passed / 0 failed.
- `cargo test --all-targets` → lib 778 passed / 0 failed; integration test crates 113 passed / 0 failed across 11 separate test binaries.

## 12. Notes

- Branch `review/vcf-writer` (off impl commit `cd1977a`). Parallel SIMD work on `main` and `perf/posterior-samply` is intentionally not touched.
- A `cargo fmt` invocation early in the run touched parallel SIMD files inadvertently; reverted via `git checkout HEAD --` on those paths before commit. The in-scope vcf_writer files were re-formatted in place and are clean.
- M13's tmp-cleanup behaviour ships but is not exercised by a dedicated test (would require a sink fault-injection seam that doesn't exist yet). Follow-up listed in the M13 entry.
- Mi16's deduplication ships at half-strength: in-crate test now shares the const with production, but the integration test keeps its own copy. Follow-up listed in the Mi16 entry.
- The four user decisions (Q1–Q4 in §3) are the single largest input that shaped the run; each is logged against the findings it affects.
