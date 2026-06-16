# Fix Application Report: ssr_types_2026-06-12.md

**Date:** 2026-06-12
**Source review:** `doc/devel/reports/reviews/ssr_types_2026-06-12.md`
**Source state reviewed against:** commit `74b5a2d` on branch `ssr-architecture`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

All review findings are resolved. The one Blocker (B1) and one Major (M1) are
Applied with new regression tests; all five Minors are Applied; both Nits are
accounted for (one kept as defensible, one verified already-consistent). The
change is confined to `src/ssr/mod.rs` + `src/ssr/types.rs`.

### Review totals
- Blockers: 1 (B1)
- Majors: 1 (M1)
- Minors: 5 (Mi1–Mi5)
- Nits: 2

### Outcome totals
- Applied: 7 (B1, M1, Mi1, Mi2, Mi3, Mi4, Mi5)
- Applied with adaptation: 0
- Already fixed: 1 (Nit2 — `purity_fraction` term verified consistent)
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0
- Kept as-is (Nit): 1 (Nit1 — flank accessors, premature-but-defensible)

### Validation summary
- `cargo fmt --check` → **ssr files clean.** The only diffs are the two
  pre-existing out-of-scope `src/vcf/*` drifts noted in review §7 (Rust 1.95.0
  toolchain); not touched by this change.
- `cargo clippy --lib --all-features` → **0 warnings in `ssr`.** Under
  `--all-targets -- -D warnings` the build fails on 2 pre-existing out-of-scope
  errors in `src/vcf/qual_refine.rs` (`needless_range_loop` L93,
  `excessive_precision` L147) — Rust 1.95.0 drift, review §7, not this change.
- `cargo test --lib ssr::` → **16 passed, 0 failed** (6 original + 10 new).
- `cargo test --lib` (whole suite, `--skip` the two pre-existing heavy
  `vcf::record_encode` i32-overflow tests) → **1039 passed, 0 failed** (exit 0).
- `cargo doc --no-deps` → **passes** (B1 was the only failure; now clean).
- `cargo audit` → not run (no manifest/dependency change in this slice; the
  dependency tree is unchanged from the reviewed state).
- Performance check → **Not applicable.** The new `src/ssr/` module is not
  reachable from any `benches/` harness.

### Unresolved high-priority findings
- None. (B1 and M1 both Applied and validated.)

### Note on the full-suite run
Two pre-existing tests in `src/vcf/record_encode.rs`
(`encode_errors_when_per_sample_depth_overflows_i32`,
`encode_errors_when_single_allele_depth_overflows_i32`) feed `support(u32::MAX)`
into the encoder and are pathologically slow (>60 s each) under the 1.95.0
container — they stalled a naive `cargo test --lib` for >11 min. They live in
`src/vcf/`, are untouched by this `src/ssr/`-only change (which has no non-test
consumers), and were `--skip`ped for the clean full-suite confirmation. They are
an out-of-scope environment/perf issue, not a regression from these fixes.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | Broken `cargo doc` gate — unresolved intra-doc links | Apply | Applied | No | `src/ssr/mod.rs` | `cargo doc --no-deps` passes | No |
| M1 | Major | `Locus` invariants unenforced → release panic | Apply | Applied | Yes (Q1: private fields + `Locus::new`) | `src/ssr/types.rs` | 10 new tests pass; clippy/doc clean | No |
| Mi1 | Minor | `MotifError` not `#[non_exhaustive]`, positional field | Apply | Applied | No | `src/ssr/types.rs` | tests pass | No |
| Mi2 | Minor | `Motif` Eq/Hash zeroed-tail coupling needs guard | Apply | Applied | No | `src/ssr/types.rs` | `motif_is_hashable_and_dedups_in_a_set` | No |
| Mi3 | Minor | Duplicated tract-range arithmetic | Apply (bundled w/ M1) | Applied | No | `src/ssr/types.rs` | tract/flank tests pass | No |
| Mi4 | Minor | `MAX_MOTIF_LEN` doc lacks source | Apply | Applied | No | `src/ssr/types.rs` | `cargo doc` clean | No |
| Mi5 | Minor | `pub` vs `pub(crate)` surface | Apply | Applied | Yes (Q2: `pub(crate)`) | `src/ssr/types.rs`, `src/ssr/mod.rs` | clippy clean (+ `#![allow(dead_code)]`) | Remove the allow when Stage 0 wires the types up |
| Nit1 | Nit | Flank/tract accessors have no caller yet | Keep (defensible) | Kept | No | None | covered by new tests | No |
| Nit2 | Nit | `purity_fraction` term consistency | Already fixed | Already fixed | No | None | spec verified | No |

## 3. Questions asked and answers

1. **M1 (Open Q1)** — How should the `Locus` invariant fix be structured: full
   encapsulation (private fields + `Locus::new`) or keep `pub` fields + a
   `validate()` method?
   - **Answer:** Private fields + `Locus::new`.
2. **Mi5 (Open Q2)** — Should the SSR shared types be `pub` or `pub(crate)`?
   - **Answer:** `pub(crate)`.
3. **Nit2 / Mi4-naming (Open Q3)** — Is `purity_fraction` the exact field name
   the catalog spec uses end-to-end?
   - **Answer (resolved by inspection):** Yes —
     `doc/devel/specs/ssr_genotyping.md:61` names the field `purity_fraction`
     explicitly; the catalog format header row (§"catalog format") and the
     shared-types doc use the same term. No rename needed.

## 4. Per-finding log

### B1 — `cargo doc` build gate broken by unresolved intra-doc links
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The design-doc paths in `src/ssr/mod.rs:13-14` were written as
  `` [`...`] `` which rustdoc parses as intra-doc links; the crate denies broken
  links, aborting `cargo doc` crate-wide. Pure markup fix.
- **Implementation summary:** Dropped the `[ ]` brackets on the two doc-path
  lines so they render as inline code spans (matching the existing `types.rs:7`
  style).
- **Review suggestion used verbatim?:** Yes (the review's suggested replacement).
- **Verification performed:** Reproduced the failure before the fix
  (`unresolved link to doc/devel/specs/ssr_genotyping.md` at `mod.rs:13`);
  `cargo doc --no-deps` passes after.
- **Files changed:** `src/ssr/mod.rs`
- **Tests added or modified:** None (doc-gate fix; verified by `cargo doc`).
- **Validation:** `cargo doc --no-deps` → exit 0, "Generated …/doc/…".
- **Follow-up:** None.

### M1 — `Locus` documented invariants unenforced; release builds underflow / slice OOB
- **Severity:** Major
- **Initial decision:** Apply (after Q1)
- **Final status:** Applied
- **Reasoning:** `Locus` had all-`pub` fields and no constructor, so the
  documented coordinate/purity invariants could be violated by a struct-literal
  caller (or the future Stage-0 builder fed bad TRF/reference input), producing a
  release-mode `u32` underflow in `tract_offset` or an out-of-bounds slice in
  `ref_tract`/`right_flank`. Inputs are external, so "the caller upholds it" is
  unsafe. User chose full encapsulation.
- **Implementation summary:**
  - Made all `Locus` fields private.
  - Added `Locus::new(chrom, start, end, motif, purity_fraction, ref_bytes,
    ref_bytes_start) -> Result<Self, LocusError>` validating
    `ref_bytes_start <= start <= end`, `end <= ref_bytes_start + ref_bytes.len()`
    (computed in `u64` to avoid its own overflow), and `purity_fraction` finite
    ∈ `[0.0, 1.0]`.
  - Added a typed `LocusError` (`#[non_exhaustive]`, `thiserror`) with
    `BadCoordinates`, `TractBeyondRefBytes`, `BadPurity` variants.
  - Added infallible field accessors (`chrom`, `start`, `end`, `motif`,
    `purity_fraction`, `ref_bytes`, `ref_bytes_start`) since the fields are now
    private and downstream stages read them.
  - Kept the `debug_assert!` (now in `tract_range`) as documentation of the
    lower-bound precondition the accessors rely on.
- **Review suggestion used verbatim?:** No — adapted.
- **Adaptation:** The review's illustrative `new` signature returns
  `Result<Self, LocusError>`; implemented as such but with named-field error
  variants (consistent with Mi1) and computed `ref_bytes_end` in `u64`. Added
  field accessors the review implied but didn't enumerate.
- **Verification performed:** Test-first — added boundary/error tests that
  exercise the exact release-mode failure modes (`ref_bytes_start > start`,
  `start > end`, `end > ref_bytes_end`, non-finite/out-of-range purity) and the
  clamped-flank edges; all pass. Confirmed accessors compile away the prior
  fallibility.
- **Files changed:** `src/ssr/types.rs`
- **Tests added or modified:** `locus_new_rejects_ref_bytes_start_after_start`,
  `locus_new_rejects_start_after_end`, `locus_new_rejects_end_beyond_ref_bytes`,
  `locus_new_rejects_non_finite_or_out_of_range_purity`,
  `locus_new_accepts_boundary_purity`,
  `locus_ref_tract_returns_tract_when_left_flank_empty`,
  `locus_right_flank_empty_when_clamped_at_contig_end`,
  `locus_is_perfect_false_on_nan_purity`; updated `sample_locus` helper and
  `locus_is_perfect_thresholds_on_purity` to build via `Locus::new`.
- **Validation:** `cargo test --lib ssr::` → 16 passed; `cargo clippy --lib` →
  0 ssr warnings; `cargo doc --no-deps` → passes.
- **Follow-up:** None. (The Stage-0 builder, when written, must route through
  `Locus::new` — now structurally enforced.)
- **Residual risk:** None identified.

### Mi1 — `MotifError` not `#[non_exhaustive]`; positional field
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Cheap hardening on brand-new internal code: `#[non_exhaustive]`
  future-proofs adding variants; a named field is more self-documenting.
- **Implementation summary:** Added `#[non_exhaustive]` to `MotifError`; changed
  `BadLength(usize)` → `BadLength { len: usize }` and updated the `#[error]`
  format string and the two existing match sites in tests.
- **Review suggestion used verbatim?:** Yes.
- **Files changed:** `src/ssr/types.rs`
- **Tests added or modified:** `motif_rejects_out_of_range_lengths` updated to
  the named-field variant.
- **Validation:** `cargo test --lib ssr::` → pass.
- **Follow-up:** None.

### Mi2 — `Motif` Eq/Hash zeroed-tail coupling
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The derived `Eq`/`Hash` correctness depends on the zeroed tail;
  a future constructor or `pub buf` could break it silently. Chose the
  guard-comment option (keeps the `Copy`, allocation-free form, per the review's
  lean) plus a regression test.
- **Implementation summary:** Added an `INVARIANT (...)` paragraph to the `Motif`
  doc stating every constructor must zero the unused tail; added a hashing test.
- **Review suggestion used verbatim?:** Yes (comment option).
- **Files changed:** `src/ssr/types.rs`
- **Tests added or modified:** `motif_is_hashable_and_dedups_in_a_set`.
- **Validation:** `cargo test --lib ssr::` → pass.
- **Follow-up:** None.

### Mi3 — Duplicated tract-range arithmetic
- **Severity:** Minor
- **Initial decision:** Apply (bundled with M1)
- **Final status:** Applied
- **Reasoning:** `ref_tract`/`left_flank`/`right_flank` each recomputed the
  offset and tract length. Bundled with M1 because both restructure the same
  accessor block.
- **Implementation summary:** Replaced `tract_offset()` with
  `tract_range(&self) -> Range<usize>`; the three accessors are now one-liners
  slicing with it. Added `use std::ops::Range;`.
- **Review suggestion used verbatim?:** Yes.
- **Files changed:** `src/ssr/types.rs`
- **Tests added or modified:** Covered by existing + new flank/tract tests.
- **Validation:** `cargo test --lib ssr::` → pass; clippy clean.
- **Follow-up:** None.

### Mi4 — `MAX_MOTIF_LEN` doc lacks source
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** No-magic-numbers convention — document why 6.
- **Implementation summary:** Extended the doc to name the microsatellite
  period ≤ 6 scope and cite `doc/devel/specs/ssr_genotyping.md` (glossary).
- **Review suggestion used verbatim?:** Adapted (cited the glossary rather than
  inventing a §1.3 section number that doesn't exist).
- **Files changed:** `src/ssr/types.rs`
- **Validation:** `cargo doc --no-deps` → passes.
- **Follow-up:** None.

### Mi5 — `pub` vs `pub(crate)` surface
- **Severity:** Minor
- **Initial decision:** Apply (after Q2)
- **Final status:** Applied
- **Reasoning:** User chose the more minimal `pub(crate)`, matching the recent
  `var_calling` demotion. No cross-crate consumer exists.
- **Implementation summary:** `MAX_MOTIF_LEN`, `Motif`, `MotifError`, `Locus`,
  `LocusError` are now `pub(crate)`. Because the SSR stage consumers aren't built
  yet, the `pub(crate)` surface has no in-crate caller and would trip
  `dead_code` under `-D warnings`; added a single module-level
  `#![allow(dead_code)]` in `src/ssr/mod.rs` with a comment to remove it once
  Stage 0 wires the types up — preferred over scattering per-item `#[allow]`s.
- **Review suggestion used verbatim?:** N/A (resolved via Q2).
- **Files changed:** `src/ssr/types.rs`, `src/ssr/mod.rs`
- **Validation:** `cargo clippy --lib` → 0 ssr warnings.
- **Follow-up:** Remove `#![allow(dead_code)]` when `ssr-catalog` consumes these
  types.
- **Residual risk:** Low — the allow is module-scoped and explicitly temporary;
  it could mask an unrelated future dead item in `ssr` until removed.

### Nit1 — Flank/tract accessors have no caller in this phase
- **Severity:** Nit
- **Final status:** Kept (premature-but-defensible, per the review)
- **Implementation summary:** None. Their contract is now documented by the new
  §8 boundary tests.

### Nit2 — `purity_fraction` term consistency
- **Severity:** Nit
- **Final status:** Already fixed (verified consistent)
- **Verification performed:** `doc/devel/specs/ssr_genotyping.md:61` names the
  field `purity_fraction`; catalog format and shared-types docs agree. No change.

## 5. Deferred findings to carry forward
None.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** No — the new `src/ssr/` module is not reachable from any
  `benches/` harness (verified: no SSR bench exists yet).
- **Baseline saved:** No (not applicable).
- **Outcome:** Skipped — no Apply touched perf-sensitive code.

## 10. Commands run
- `./scripts/dev.sh cargo doc --no-deps`
- `./scripts/dev.sh cargo test --lib ssr::`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo clippy --lib --all-features`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh rustfmt --edition 2024 src/ssr/types.rs src/ssr/mod.rs`
- `./scripts/dev.sh cargo test --lib -- --skip encode_errors_when_per_sample_depth_overflows_i32 --skip encode_errors_when_single_allele_depth_overflows_i32`

## 11. Command results
- `cargo doc --no-deps` → exit 0; generated docs (B1 fixed).
- `cargo test --lib ssr::` → exit 0; **16 passed, 0 failed**.
- `cargo clippy --lib --all-features` → 0 ssr warnings (2 pre-existing
  out-of-scope `src/vcf/qual_refine.rs` warnings remain).
- `cargo clippy --all-targets -- -D warnings` → fails on 2 pre-existing
  out-of-scope errors in `src/vcf/qual_refine.rs` (review §7); ssr clean.
- `cargo fmt --check` → ssr clean; 2 pre-existing out-of-scope `src/vcf/*` drifts.
- `cargo test --lib --skip <2 heavy vcf overflow tests>` → exit 0;
  **1039 passed, 0 failed, 1 ignored, 2 filtered out** (the 2 skipped), 32.98s.

## 12. Notes
- The full-suite run was stalled by two pre-existing, out-of-scope
  `vcf::record_encode` i32-overflow tests that materialize `u32::MAX`
  observations and take >60 s each under the 1.95.0 container. They are
  unrelated to this `src/ssr/`-only change and were `--skip`ped for a clean
  confirmation.
- `cargo audit` not run: this slice changes no manifest or dependency.
- The `#![allow(dead_code)]` in `src/ssr/mod.rs` is the one piece of debt this
  run introduces, and it is intentional and temporary (tracked in Mi5).
</content>
