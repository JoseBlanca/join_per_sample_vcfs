# Fix Application Report: ssr_pileup_mark2_2026-06-17.md

**Date:** 2026-06-17
**Source review:** `doc/devel/reports/reviews/ssr_pileup_mark2_2026-06-17.md`
**Source state reviewed against:** branch `ssr-pileup-mark2` @ `5cce21f`
**Execution mode:** interactive
**Overall status:** Partial (all 3 Blockers + 9 of 12 Majors applied; remainder deferred with reasons)

---

## 1. Executive summary

### Review totals
- Blockers: 3
- Majors: 12
- Minors: 14 (Mi1–Mi14 + the grouped "Mi-docs")
- Nits: 1 group

### Outcome totals
- Applied: 18
- Applied with adaptation: 2 (B3, M9)
- Already fixed: 0
- Deferred: 9 (M4, M12, Mi2, Mi8, Mi9, Mi10, Mi11, Mi12, Nits)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean (after one `cargo fmt` pass over the new code).
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, no warnings.
- `cargo test --lib --tests --all-features` → 0, **1129 lib + all integration tests pass, 0 failed** (2 ignored = live trf-mod + one pre-existing).
- `cargo test --all-targets --all-features` → non-zero, **pre-existing out-of-scope** `psp_writer_perf` bench panic only (`benches/psp_writer_perf.rs:386`, index OOB — documented in PROJECT_STATUS as the known `--all-targets` failure; not touched by this work). All lib/integration tests passed before the bench ran.
- `cargo doc --no-deps` → 0, **now green** (was the B1 failure).
- `cargo audit` → not run (cargo-audit not installed in the container; no dependency changes this run).
- Performance check → skipped — no Apply touched a hot path covered by `benches/` (the Mark-2 `ssr-pileup` bench was deleted at the cutover; no bench reaches the changed code).

### Unresolved high-priority findings
- None at Blocker/Major confidence that block merge. Deferred Majors: **M4** (HMM-model provenance — schema design call) and **M12** (DP-state enum — determinism-critical refactor, now guarded by the B3 tests). Both carry written justification below.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | Broken intra-doc link fails `cargo doc` | Apply | Applied | No | mod.rs, fetch_reads.rs | doc green | No |
| B2 | Blocker | Length-inconsistent record panics the run | Ask→Apply | Applied | Yes (drop+count) | fetch_reads.rs, footprint.rs, driver.rs | Pass | No |
| B3 | Blocker | Byte-identity untested on cap/multi-chunk paths | Apply | Applied with adaptation | No | driver.rs, fetch_reads.rs | Pass | Yes (allele-diverse cap test) |
| M1 | Major | `Io(#[from])` collapses 5 I/O sites | Apply | Applied | No | driver.rs | Pass | No |
| M2 | Major | Zero-length-flank delimitation untested/reachable | Apply | Applied | Yes (gate at Stage-0) | postprocess.rs | Pass | No |
| M3 | Major | `MIN_FLANK_BP` absent from `.ssr.psp` header | Apply | Applied | No | driver.rs | Pass | No |
| M4 | Major | HMM model constants not recorded | Defer | Deferred | No | None | N/A | Yes |
| M5 | Major | Cross-platform float determinism unscoped | Apply | Applied | No | driver.rs | Pass | No |
| M6 | Major | `..FilterCounts::default()` masks new buckets | Apply | Applied | No | driver.rs | Pass | No |
| M7 | Major | `TimestampFormat` discards cause | Apply | Applied | No | driver.rs | Pass | No |
| M8 | Major | `ref_to_read` indel branches untested | Apply | Applied | No | footprint.rs | Pass | No |
| M9 | Major | `process_locus` outcome routing untested | Apply | Applied with adaptation | No | driver.rs | Pass | Yes (BorderOffEnd e2e) |
| M10 | Major | `build_ssr_writer_header` errors untested | Apply | Applied | No | driver.rs | Pass | No |
| M11 | Major | `passes_quality_gate` boundary untested | Apply | Applied | No | alignment.rs | Pass | No |
| M12 | Major | DP state primitive, not enum | Defer | Deferred | No | None | N/A | Yes |
| Mi1 | Minor | `pick` empty-slice guard | Apply | Applied | No | alignment.rs | Pass | No |
| Mi2 | Minor | `MAX_READS_PER_LOCUS` test-only default | Defer | Deferred | No | None | N/A | Yes |
| Mi3 | Minor | Q1 index undocumented | Apply | Applied | No | alignment.rs | Pass | No |
| Mi4 | Minor | `qc_counts` unchecked casts | Apply | Applied | No | driver.rs | Pass | No |
| Mi5 | Minor | `cfg.cap as i64` unchecked | Apply | Applied | No | driver.rs | Pass | No |
| Mi6 | Minor | Reservoir "unbiased" doc overstates | Apply | Applied | No | fetch_reads.rs | Pass | No |
| Mi7 | Minor | `input_crams` misleading name | Apply | Applied (local) | No | driver.rs | Pass | Yes (psp field) |
| Mi8 | Minor | `delimit_read` too long | Defer | Deferred | No | None | N/A | Yes |
| Mi9 | Minor | `build_ssr_writer_header` too long | Defer | Deferred | No | None | N/A | Yes |
| Mi10 | Minor | Duplicated `Locus` test fixtures | Defer | Deferred | No | None (partial: BAM helpers refactored) | N/A | Yes |
| Mi11 | Minor | `LazyLock` → `const` scalars | Defer | Deferred | No | None | N/A | Yes |
| Mi12 | Minor | `passes_quality_gate` takes whole scratch | Defer | Deferred | No | None | N/A | Yes |
| Mi-docs | Minor | Stale "still to land" module docs | Apply | Applied | No | mod.rs, fetch_reads.rs | doc green | No |
| Nits | Nit | M/I/D casts, `hap`/`nn`/`prev`/`cur`, region comment | Defer | Deferred | No | None | N/A | Yes |

## 3. Questions asked and answers

1. **B2** — When a fetched read is length-inconsistent (empty/short `QUAL`, or a CIGAR that consumes more read bases than `seq` holds), how should ssr-pileup handle it?
   - **Answer:** Drop + count locally (guard in the SSR fetch/process path, tally as a filtered drop, no change to the shared segment reader).

(M2's "gate at Stage-0" decision was taken in the prior conversation turn before this skill run; recorded here for completeness.)

## 4. Per-finding log

### B1 — Broken intra-doc link fails `cargo doc`
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Verified `cargo doc` failure reproduced; the rebuild renamed `locus_record`→`locus_tally`/`aggregate`→`tally`. The broken link lived inside a stale "still to land" paragraph, so fixed together with Mi-docs.
- **Implementation summary:** Rewrote the `mod.rs` and `fetch_reads.rs` module docs to the landed shape; the dead `super::locus_record::aggregate` link is replaced by a live `super::locus_tally::tally` reference.
- **Review suggestion used verbatim?:** No (rewrote the surrounding prose, not just the link).
- **Verification performed:** `cargo doc --no-deps` now exits 0 and generates docs.
- **Files changed:** src/ssr/pileup/mod.rs, src/ssr/pileup/fetch_reads.rs
- **Tests added or modified:** None (doc-only).
- **Validation:** `cargo doc --no-deps` → 0, Generated index.html.
- **Follow-up:** None. **Residual risk:** None.

### B2 — Length-inconsistent record panics the run
- **Severity:** Blocker
- **Initial decision:** Ask → Apply (drop + count locally)
- **Final status:** Applied
- **Reasoning:** Confirmed at the code level (no seq/qual/CIGAR-length reconciliation on the SSR fetch path) and **empirically reachable**: the new test shows noodles round-trips a `*`-QUAL BAM read as an empty `qual` buffer, which the old `&read.qual[region]` slice would panic on. User chose drop+count.
- **Implementation summary:** In `fetch_locus_reads`, drop a read where `qual.len() != seq.len()` or `cigar_read_len(cigar) != seq.len()`, before reservoir admission; count it in a new `LocusReads.malformed` field, folded into `n_filtered` by `qc_counts`. Added `footprint::cigar_read_len`. Added a defense-in-depth `r_start.min(read_len)` clamp in `extract_region` so an out-of-range index can never reach the slice.
- **Review suggestion used verbatim?:** No — adapted (the review sketched either a `process_locus` guard or a `record_buf_to_mapped_read` change; chose the fetch-boundary drop+count, which keeps the change in-scope and counts cleanly).
- **Verification performed:** `fetch_drops_length_inconsistent_reads_and_counts_them` (empty-QUAL read dropped, `malformed == 1`, clean read still admitted); `qc_counts_folds_malformed_reads_into_n_filtered`; `cigar_read_len_sums_read_consuming_ops`. Test-first: the empty-QUAL fixture demonstrates the exact input class the review flagged.
- **Files changed:** src/ssr/pileup/fetch_reads.rs, src/ssr/pileup/footprint.rs, src/ssr/pileup/driver.rs
- **Tests added or modified:** `fetch_drops_length_inconsistent_reads_and_counts_them`, `qc_counts_folds_malformed_reads_into_n_filtered`, `cigar_read_len_sums_read_consuming_ops` (+ `aln_record_no_qual` helper).
- **Validation:** `cargo test --lib ssr::pileup` → 0, 51 passed.
- **Follow-up:** None required. **Residual risk:** A CIGAR that over-consumes the read cannot be written by noodles (it validates on write), so that trigger is only reachable from a corrupt file; the drop covers it and the clamp is the backstop.

### B3 — Byte-identity untested on cap/multi-chunk paths
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** The determinism *mechanism* was verified sound in review; the gap was test coverage on the reservoir-eviction and multi-chunk paths. Added all three suggested tests.
- **Implementation summary:** (1) `reservoir_keeps_a_deterministic_subset_when_offered_far_past_capacity` — offers 10 000 into cap 8, asserts subset + cap-sized + cross-run identity. (2) `run_output_is_thread_invariant_across_multiple_chunks` — 80 loci (> `MIN_FETCH_CHUNK`), T=1 vs T=8, asserts identical records **and** strictly ascending order (chunk-order preserved). (3) `run_output_is_thread_invariant_when_the_reservoir_cap_bites` — 50 reads at one locus, cap 4, T=1 vs T=4, asserts identical records, `depth == 50`, `observed` capped at 4. Refactored the BAM fixtures (`write_bam` + `spanning_read`, contig length + reps parameterized) and added `run_config_with_cap`.
- **Review suggestion used verbatim?:** No.
- **Adaptation:** The cap-bites e2e test uses identical-allele reads, so it exercises the eviction branch end-to-end and pins depth/observed-count thread-invariance, but it does not vary *which allele* survives (all reads show `CACACA`). The allele-subset ordering invariant is covered instead by the reservoir unit test (`is_deterministic_for_a_fixed_seed_and_order` over a distinguishable 1..100 stream + the new subset test). An allele-diverse cap-bites fixture is deferred (fixture cost; see Follow-up).
- **Verification performed:** The three tests pass; multi-chunk confirmed to split (80 > 64).
- **Files changed:** src/ssr/pileup/driver.rs, src/ssr/pileup/fetch_reads.rs
- **Tests added or modified:** the three named above + fixture refactor.
- **Validation:** `cargo test --lib ssr::pileup` → 0, 51 passed.
- **Follow-up:** An allele-diverse cap-bites e2e (reads showing mixed `CACACA`/`CACACACA`) would tighten the subset-ordering assertion end-to-end. **Residual risk:** Low — the per-locus-atomic + seeded-reservoir + multi-chunk-order properties are all now tested; only the e2e allele-subset selection is covered indirectly.

### M1 — `Io(#[from])` collapses 5 I/O sites
- **Severity:** Major — **Initial:** Apply — **Final:** Applied
- **Reasoning/Implementation:** Replaced the blanket `Io(#[from] io::Error)` with five operation-named variants (`OpenCatalog`, `CreateOutput`, `FlushOutput`, `SyncOutput`, `RenameOutput`), each carrying its path(s) and `#[source]`; mapped each `?` site in `run`.
- **Review suggestion used verbatim?:** Yes (variant set as suggested).
- **Files changed:** src/ssr/pileup/driver.rs — **Tests:** None (error-shape; exercised by the happy-path e2e). **Validation:** clippy + tests pass. **Follow-up/Residual:** None.

### M2 — Zero-length-flank delimitation (reachable)
- **Severity:** Major (→ Blocker once reachability confirmed) — **Initial:** Apply — **Final:** Applied
- **Reasoning:** Open Question 1 resolved "yes" — `finish_locus` clamps flanks at contig ends with no minimum, so a tract at contig position 0 / last base yields an empty flank, which `Locus::new` permits and the delimiter mishandles. User chose to gate at Stage-0.
- **Implementation summary:** `finish_locus` now returns `None` when `ref_start == new_start` (empty left) or `ref_end == new_end` (empty right); updated the module-doc step list (added step 6). (Applied in the prior turn.)
- **Files changed:** src/ssr/catalog/postprocess.rs — **Tests:** `build_loci_drops_locus_with_empty_left_flank`, `..._right_flank`; fixed the collateral `build_loci_drops_period_one_homopolymer_and_spares_the_neighbour_ssr` fixture (added a right flank). **Validation:** `cargo test --lib ssr::catalog` → 0, 25 passed. **Follow-up/Residual:** None.

### M3 — `MIN_FLANK_BP` absent from header
- **Severity:** Major — **Initial:** Apply — **Final:** Applied
- **Implementation:** Added `reach_min_flank_bp` to the `.ssr.psp` provenance parameters (imported `footprint::MIN_FLANK_BP`). **Files:** driver.rs. **Tests:** asserted present in `build_ssr_writer_header_errors_on_overflow_and_missing_md5` (M10). **Validation:** pass. **Follow-up/Residual:** None.

### M4 — HMM model constants not recorded
- **Severity:** Major — **Initial:** Defer — **Final:** Deferred
- **Reasoning:** Recording the gap-open/extend/emission model in the header is a provenance-schema design decision: `ParameterValue` is Integer/String/Boolean, so a float gap-open (`2.9e-5`) needs a scaling or string convention, and a model-version tag is a new invention. The model is frozen for the `format_version`. Deferring rather than inventing a schema on the fly (no clearly-correct single path).
- **Files changed:** None. **Follow-up:** Decide a `delimiter_model` provenance convention (string tag + format_version binding, or a scaled-int gap-open) and record it. **Residual risk:** Two `.ssr.psp` from differently-tuned builds are byte-different with no in-artifact discriminator — acceptable while the model is frozen and pre-alpha.

### M5 — Cross-platform float determinism unscoped
- **Severity:** Major — **Initial:** Apply — **Final:** Applied
- **Implementation:** Added a paragraph to the driver module doc scoping byte-identity to a fixed target+toolchain and stating cross-platform identity is not guaranteed (transcendental `f64` ties). **Files:** driver.rs. **Tests:** None (doc). **Validation:** doc green. **Follow-up/Residual:** If cross-platform identity is later required, a golden `.ssr.psp` + fp-contract pinning is the path (noted in the doc).

### M6 — `..FilterCounts::default()` masks new buckets
- **Severity:** Major — **Initial:** Apply — **Final:** Applied
- **Implementation:** Named all ten `FilterCounts` fields in the test literal (dropped `..default()`), matching the `FilterCounts::merge` exhaustive-destructure convention. **Files:** driver.rs (test). **Tests:** `qc_counts_excludes_dups_from_coverage_but_keeps_them_filtered` rewritten. **Validation:** pass. **Follow-up/Residual:** None.

### M7 — `TimestampFormat` discards cause
- **Severity:** Major — **Initial:** Apply — **Final:** Applied
- **Implementation:** `TimestampFormat { value: String, #[source] source: DatetimeParseError }`, populated from the parsed string. **Files:** driver.rs. **Validation:** clippy + tests pass. **Follow-up/Residual:** None.

### M8 — `ref_to_read` indel branches untested
- **Severity:** Major — **Initial:** Apply — **Final:** Applied
- **Implementation:** `extract_region_maps_window_across_an_internal_deletion` (M8 D2 M10 → 0..16) and `..._insertion` (M8 I2 M10 → 0..20). **Files:** footprint.rs (tests). **Validation:** pass. **Follow-up/Residual:** None.

### M9 — `process_locus` outcome routing untested
- **Severity:** Major — **Initial:** Apply — **Final:** Applied with adaptation
- **Reasoning/Adaptation:** Added `process_locus_routes_sequence_and_low_quality_outcomes` (one locus, a clean read + a read with the tract dimmed below Phred 15 → `depth == 2`, `n_low_quality == 1`, `observed == [(CACACA, 1)]`). This covers the Sequence vs LowQuality gate routing end-to-end. The `BorderOffEnd` arm is **not** in this test — a synthetic read that reaches the delimiter and yields `BorderOffEnd` requires a soft-clipped fixture whose outcome is fiddly to pin without becoming tautological; it is covered at the unit level by `alignment::missing_left_flank_is_border_off_end` / `missing_right_flank_is_border_off_end`, and the `process_locus` arm is a trivial passthrough.
- **Files changed:** driver.rs (test). **Validation:** pass. **Follow-up:** A soft-clipped off-end e2e read to cover the `BorderOffEnd` routing arm. **Residual risk:** Low.

### M10 — `build_ssr_writer_header` errors untested
- **Severity:** Major — **Initial:** Apply — **Final:** Applied
- **Implementation:** `build_ssr_writer_header_errors_on_overflow_and_missing_md5` — hand-built `ContigList` asserting `ContigLengthOverflow` (length > u32::MAX), `MissingMd5` (md5 None), and that a valid contig records `quality_q1_threshold`/`reservoir_cap`/`reach_min_flank_bp`. **Files:** driver.rs (test). **Validation:** pass. **Follow-up/Residual:** None.

### M11 — `passes_quality_gate` boundary untested
- **Severity:** Major — **Initial:** Apply — **Final:** Applied
- **Implementation:** `quality_gate_is_inclusive_at_the_threshold_and_handles_singletons` — `[15]`→pass, `[14]`→fail, `[15,15]`→pass, plus a 4-element region whose lower-quartile element is exactly 15. **Files:** alignment.rs (test). **Validation:** pass. **Follow-up/Residual:** None.

### M12 — DP state primitive, not enum
- **Severity:** Major (Medium confidence, idiomatic) — **Initial:** Defer — **Final:** Deferred
- **Reasoning:** Converting `M`/`I`/`D` to a `#[repr(u8)] enum` touches ~14 sites in the hottest, most determinism-critical function and its backpointer matrix. It is behavior-preserving but a sizeable refactor better done on its own, now guarded by the new B3 determinism tests. Not implemented to keep this run's diff minimal and verifiable.
- **Files changed:** None. **Follow-up:** Focused `DpState` refactor gated by `run_output_is_thread_invariant_*`. **Residual risk:** None functional; the `_` traceback arm and casts remain a maintainability hazard.

### Minor findings
- **Mi1 (Applied):** `pick` — added `# Panics` doc + `debug_assert!(!cands.is_empty())`. (alignment.rs)
- **Mi2 (Deferred):** `MAX_READS_PER_LOCUS` doc vs the off-file CLI default depends on Q3 (CLI default), which spans an out-of-scope file. Doc-only Minor.
- **Mi3 (Applied):** Documented the Q1 nearest-rank lower-quartile index `⌊(n-1)/4⌋`. (alignment.rs)
- **Mi4 (Applied):** `qc_counts` now uses `saturating_add` + `u32::try_from(..).unwrap_or(u32::MAX)`. (driver.rs)
- **Mi5 (Applied):** `reservoir_cap` parameter uses `i64::try_from(cfg.cap).unwrap_or(i64::MAX)`. (driver.rs)
- **Mi6 (Applied):** Reservoir doc reworded to "effectively-uniform (modulo bias ≤ seen/2^64), deterministic" with the rationale. (fetch_reads.rs)
- **Mi7 (Applied, local):** Renamed the local `input_crams` → `input_alignment_files`; the `WriterProvenance.input_crams` field rename spans `src/psp/` (out of scope) and is left as a follow-up.
- **Mi8 / Mi9 (Deferred):** Splitting `delimit_read` / `build_ssr_writer_header` are refactors not required by any higher fix.
- **Mi10 (Deferred, partial):** The BAM fixtures were refactored (`write_bam` + `spanning_read`) as a side effect of B3; the cross-module `Locus` fixture dedup touches `src/ssr/types.rs` (out of scope) and is left.
- **Mi11 (Deferred):** `LazyLock` `GAP_EXTEND_PROB`/`INS_EMIT_LN` → `const` literals carries a byte-identity risk (a literal may differ from `.exp()`/`.ln()` by 1 ULP and shift output); deliberately not changed without a round-trip-verified literal.
- **Mi12 (Deferred):** `passes_quality_gate` taking the whole scratch is a Minor API refactor not required here.
- **Mi-docs (Applied):** Stale "still to land" prose in `mod.rs`/`fetch_reads.rs` rewritten (bundled with B1).

### Nits (Deferred)
The M/I/D cast count (subsumed by M12), `hap`/`nn`/`prev`/`cur` naming, and the `region`/`r` relative-range comment were not addressed — cosmetic, no behavioral impact.

## 5. Deferred findings to carry forward
- **M4** — record an HMM-model provenance tag (schema decision).
- **M12** — `DpState` enum refactor (guarded by the new B3 tests).
- **Mi2** — `MAX_READS_PER_LOCUS` doc vs CLI default (Q3).
- **Mi8 / Mi9** — `delimit_read` / `build_ssr_writer_header` function splits.
- **Mi10** — cross-module `Locus` test-fixture dedup (in `src/ssr/types.rs`).
- **Mi11** — `LazyLock` → `const` (needs ULP round-trip check).
- **Mi12** — `passes_quality_gate` scratch parameter.
- **Mi7 (residual)** — `WriterProvenance.input_crams` crate-wide rename (`src/psp/`).
- **B3 / M9 (follow-up)** — allele-diverse cap-bites e2e; `BorderOffEnd` routing e2e.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** No.
- **Baseline saved:** No.
- **Benches run:** None.
- **Outcome:** Skipped — no Apply touched a hot path covered by `benches/`. The Mark-2 `ssr-pileup` bench was deleted at the cutover, and no remaining bench reaches `src/ssr/pileup/` or `src/ssr/catalog/postprocess.rs`. The B2 per-read length check is a cheap comparison on a non-benched path; the `extract_region` clamp is one `.min()`.
- **Notes:** The pre-existing `psp_writer_perf` bench panic is unrelated (a bench-harness off-by-one, documented in PROJECT_STATUS).

## 10. Commands run
- `cargo fmt` / `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib ssr::pileup`
- `cargo test --lib ssr::catalog`
- `cargo test --lib --tests --all-features`
- `cargo test --all-targets --all-features`
- `cargo doc --no-deps`
(all via `./scripts/dev.sh` in the dev container.)

## 11. Command results
- `cargo fmt --check` → 0, clean (after one `cargo fmt`).
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, no warnings.
- `cargo test --lib ssr::pileup` → 0, 51 passed.
- `cargo test --lib --tests --all-features` → 0, 1129 lib + integration passed, 0 failed (2 ignored).
- `cargo test --all-targets --all-features` → non-zero, only the pre-existing `psp_writer_perf` bench panic.
- `cargo doc --no-deps` → 0, generated.
- `cargo audit` → not run (not installed; no dep changes).

## 12. Notes
- B2's empty-`QUAL` test empirically confirmed the review's reachability assumption: noodles decodes a `*`-quality BAM record to an empty `qual` buffer, the exact input that panicked the pre-fix slice.
- M2 was applied in the conversation turn immediately before this skill run (Stage-0 flank gate); recorded here so the source review's findings are fully accounted for.
- Nothing was committed; all changes are in the working tree on branch `ssr-pileup-mark2`.
