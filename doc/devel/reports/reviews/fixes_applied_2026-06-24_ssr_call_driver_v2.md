# Fix Application Report: ssr_call_driver_2026-06-24.md (v2 — deferred-Minor follow-up)

**Date:** 2026-06-24
**Source review:** `doc/devel/reports/reviews/ssr_call_driver_2026-06-24.md`
**Source state reviewed against:** branch `ssr-cohort`, HEAD `f79ca8e` (after the v1 fix run)
**Execution mode:** interactive
**Overall status:** Completed (6 of the 8 carried-over Minors Applied; 2 Deferred by user decision)

---

## 0. Context

This is the **second** fix run against the `ssr_call_driver_2026-06-24` review. The
[v1 run](fixes_applied_2026-06-24_ssr_call_driver.md) applied the Blocker (B1), all six
Majors (M1–M6), and 8 Minors (Mi6–Mi12, Mi15), deferring 8 Minors as follow-ups: **Mi1,
Mi2, Mi3, Mi4, Mi5, Mi13, Mi14, Mi16**. This run picks up that deferred set. Three of the
eight needed a product/policy or measure-first decision; those were put to the user before
any code changed (§3). The user chose: Mi14 → numeric QUAL; Mi13 → reject non-ACGTN; Mi3 /
Mi4 → keep deferred until an SSR bench exists.

## 1. Executive summary

### Carried-over findings (the v1-deferred set)
- Minors carried over: 8 (Mi1, Mi2, Mi3, Mi4, Mi5, Mi13, Mi14, Mi16)

### Outcome totals (this run)
- Applied: 6 (Mi1, Mi2, Mi5, Mi13, Mi14, Mi16)
- Applied with adaptation: 0 (Mi1 adapted the helper's return shape — noted per finding)
- Already fixed: 0
- Deferred: 2 (Mi3, Mi4 — user decision: bench-gated, no SSR bench yet)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary (final full gate, host toolchain)
- `cargo fmt --check` → 0
- `cargo clippy --all-targets --all-features -- -D warnings` → 0
- `cargo test --lib --all-features` → 0, `1292 passed; 0 failed; 2 ignored` (+5 vs the 1287 v1 baseline)
- `cargo test --all-targets --all-features` → the **only** failure is the pre-existing
  `benches/psp_writer_perf.rs:386` panic (baseline; not this work). All lib + integration tests pass.
- `cargo doc --no-deps` → 0
- `cargo audit` → not run (cargo-audit not installed; diff adds no dependencies — `Cargo.toml`/`Cargo.lock` untouched)
- Performance check: **skipped** — no bench under `benches/` references `ssr`/`cohort::`, so no Apply touched a bench-covered hot path.

### Unresolved high-priority findings
- None. The carried-over set is entirely Minor. The two Deferred (Mi3/Mi4) are bench-gated
  allocation levers, held by explicit user decision until an SSR cohort bench exists.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation |
|---|---|---|---|---|---|---|---|
| Mi1 | Minor | duplicated read-to-nearest-allele attribution | Apply | Applied | No | attribution.rs (new), em.rs, prepass.rs, vcf_out.rs, mod.rs | Pass |
| Mi2 | Minor | `compute_data_ll` 11-arg signature | Apply | Applied | No | em.rs | Pass |
| Mi3 | Minor | per-round alloc churn in `compute_data_ll` | Ask → Defer | Deferred | Yes | None | N/A |
| Mi4 | Minor | per-locus `f_present` alloc in the sweep | Ask → Defer | Deferred | Yes | None | N/A |
| Mi5 | Minor | tuple bins primitive obsession | Apply | Applied | No | param_estimation.rs, prepass.rs | Pass |
| Mi13 | Minor | `from_utf8_lossy` alleles | Ask → Apply | Applied | Yes | merge.rs, vcf_out.rs | Pass |
| Mi14 | Minor | `QUAL=.` for a variable-but-zero locus | Ask → Apply | Applied | Yes | vcf_out.rs | Pass |
| Mi16 | Minor | once-per-run `level_per_group.clone()` | Apply | Applied | No | driver.rs | Pass |

## 3. Questions asked and answers

1. **Mi14** — A variable locus whose site-QUAL proxy clamps to exactly `0.0` currently
   prints `QUAL=.`; keep `.` (documented) or always print numeric?
   - **Answer:** Always print numeric.
2. **Mi13** — REF/ALT allele bytes go through `from_utf8_lossy` (a bad byte → silent U+FFFD).
   Apply M5-style reject-with-typed-error to alleles, or defer?
   - **Answer:** Reject non-ACGTN.
3. **Mi3/Mi4** — Apply the bench-gated scratch-reuse now (no SSR bench to measure the win),
   or defer until a bench exists?
   - **Answer:** Defer until a bench exists.

## 4. Per-finding log

### Mi1 — duplicated read-to-nearest-allele attribution
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The "find the nearest parent allele by `|read_units − allele_units|`, sign
  the slip" kernel was re-derived at three in-scope sites (`em::attribute_locus`,
  `prepass::accumulate_locus`, `vcf_out::allele_balance`) with the tie-break implicit in each
  (`min_by_key` first-min twice; `<=`-to-first in `allele_balance`). All three resolve a tie
  toward the first/lowest-index parent today, so a single shared primitive is behaviour-
  preserving and removes the drift risk the review flagged.
- **Implementation summary:** new `src/ssr/cohort/attribution.rs` with
  `nearest_parent(read_units, parent_units) -> Option<(index, delta)>`, the tie-break
  (`min_by_key` keeps the first minimum → lowest index) documented as the one shared rule.
  Wired into all three sites: `attribute_locus` (parents = `genotype_units`),
  `accumulate_locus` (parents = the peaks' `repeat_len`, collected once per sample),
  `allele_balance` (parents = the two `genotype_units`). Each site's downstream logic is
  unchanged: `delta != 0` ⇔ the old `read_units != parent`; `parent_idx == 0` ⇔ the old
  `<=`-to-`a`.
- **Review suggestion used verbatim?:** No (adapted).
- **Adaptation:** the review proposed a helper returning `(delta, parent, count)`; the shared
  decision is only the tie-broken nearest parent, so the primitive returns `(index, delta)`
  and each caller derives `parent`/`count` locally (the count never leaves the caller). This
  keeps the helper minimal and each site's accounting visible.
- **Verification performed:** the EM and driver **byte-identity-across-thread-count** tests
  (`cohort_em_is_byte_identical_across_thread_counts`, `run_is_byte_identical_across_thread_counts`)
  and every genotype-correctness/allele-balance test pass unchanged — the tie-break is
  preserved at all three sites. Added 3 direct unit tests on the primitive (nearest pick +
  signed delta; the first-parent tie-break, including the order-sensitivity; empty input).
- **Files changed:** `src/ssr/cohort/attribution.rs` (new), `src/ssr/cohort/em.rs`,
  `src/ssr/cohort/prepass.rs`, `src/ssr/cohort/vcf_out.rs`, `src/ssr/cohort/mod.rs`.
- **Tests added or modified:** `attribution::tests::{picks_the_nearest_parent_and_signs_the_delta,
  breaks_a_distance_tie_toward_the_first_parent, returns_none_for_no_parents}` (new).
- **Validation:** `cargo test --lib ssr::cohort` → 0, `154 passed`; fmt/clippy/doc → 0.
- **Follow-up:** the soft per-read responsibility split (the deferred refinement the modules
  document) would replace this hard assignment in every caller at once — now a single edit.
- **Residual risk:** None observed — byte-identity preserved.

### Mi2 — `compute_data_ll` 11-arg signature
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `compute_data_ll` took 11 positional args (two adjacent `f64`s and two
  adjacent slices invited a silent transposition) and was called identically at two sites.
  The iteration-invariant locus shaping is the natural bundle.
- **Implementation summary:** new `struct LocusModel<'a> { locus, candidates, cand_units,
  genotypes, distinct }`, built once in `run_locus_em_with` and passed by reference. This drops
  `compute_data_ll` to 7 args `(model, params, level_per_group, theta, level_multiplier,
  lambda, scratch)` — below clippy's `too_many_arguments` threshold, so its
  `#[allow(clippy::too_many_arguments)]` was removed. Both call sites now pass `&model`.
- **Review suggestion used verbatim?:** Yes (the `LocusModel` bundle; additionally folded
  `locus` into it so the count drops under the lint and the allow is deleted).
- **Adaptation:** included `locus` in the bundle (the review listed only the four shaping
  fields) so the arg count falls to 7 and the `#[allow]` is removed outright.
- **Verification performed:** all 18 `ssr::cohort::em` tests pass unchanged (the EM output is
  identical — pure signature refactor); byte-identity tests green.
- **Files changed:** `src/ssr/cohort/em.rs`.
- **Tests added or modified:** None (behaviour-identical; existing EM tests cover it).
- **Validation:** `cargo test --lib ssr::cohort::em` → 0, `18 passed`; fmt/clippy → 0.
- **Follow-up:** None.
- **Residual risk:** None.

### Mi3 — per-round alloc churn in `compute_data_ll`
- **Severity:** Minor
- **Initial decision:** Ask (perf, bench-gated)
- **Final status:** Deferred
- **Reasoning:** The review says bench-gate the scratch-reuse, and no bench under `benches/`
  exercises the SSR cohort path, so the perf win cannot be measured (only correctness could
  be validated). Put to the user; the user chose to defer until an SSR bench exists.
- **Implementation summary:** None.
- **Follow-up:** revisit with Mi4 once an SSR cohort bench lands (then measure the wall cost).
- **Residual risk:** None (no change).

### Mi4 — per-locus `f_present` alloc in the parallel sweep
- **Severity:** Minor
- **Initial decision:** Ask (perf, bench-gated)
- **Final status:** Deferred
- **Reasoning:** Same bench-gated rationale as Mi3; deferred by the same user decision.
- **Implementation summary:** None.
- **Follow-up:** `map_init` scratch in `write_genotyped_chunk`'s `par_iter`, measured against
  a future SSR bench.
- **Residual risk:** None (no change).

### Mi5 — tuple bins primitive obsession
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `SampleStutterStats.by_length: Vec<(u16, u64, u64)>` (repeated in `add_bin`,
  `merge_sample_stats`, `fit_level`) and the per-locus `Vec<(u16, u64)>` allele tally
  (`add_allele_copies`) carried meaning only in doc comments (`bin.1`, `bin.2`).
- **Implementation summary:** added `struct LengthBin { length, faithful, slipped }`
  (param_estimation.rs) and `struct AlleleCopies { length, copies }` (prepass.rs); converted
  the field, `add_bin`, `merge_sample_stats`, `fit_level`, `accumulate_locus`,
  `add_allele_copies`, and the one test reader to named-field access. The `Vec` + linear scan
  is kept (a one-line note records that the distinct-length cardinality is tiny).
- **Review suggestion used verbatim?:** Yes for `LengthBin`; additionally named the
  `(u16, u64)` allele tally (`AlleleCopies`) the review listed as the same smell.
- **Adaptation:** None beyond also doing the second tuple the review named.
- **Verification performed:** the prepass + inbreeding determinism/level-fit tests pass
  unchanged (pure field-naming refactor; same arithmetic).
- **Files changed:** `src/ssr/cohort/param_estimation.rs`, `src/ssr/cohort/prepass.rs`.
- **Tests added or modified:** the `het … bins at both alleles` test updated to `bin.length`.
- **Validation:** `cargo test --lib ssr::cohort` → 0, `151 passed`; fmt/clippy → 0.
- **Follow-up:** None.
- **Residual risk:** None.

### Mi13 — `from_utf8_lossy` alleles
- **Severity:** Minor
- **Initial decision:** Ask → Apply (user chose reject)
- **Final status:** Applied
- **Reasoning:** REF/ALT bytes come from untrusted catalog reference tracts and `.ssr.psp`
  observed sequences; a non-nucleotide byte became a silent U+FFFD in the VCF. The user chose
  to reject (consistent with M5's name validation), so malformed input fails loud.
- **Implementation summary:** new `SsrMergeError::InvalidAlleleByte { chrom, start, byte }`
  and a `first_non_acgtn` helper (case-insensitive A/C/G/T/N — SSR loci may be soft-masked
  lowercase; control chars / non-UTF8 / IUPAC ambiguity rejected). Validated at the merge
  boundary in `next_locus`: the catalog reference tract, then each present sample's observed
  sequences — the single chokepoint where untrusted bytes enter, and it propagates through
  both passes via the existing `?` on merger items. `vcf_out::format_vcf_record`'s two
  `from_utf8_lossy` sites become `std::str::from_utf8(...).expect(ACGTN_GUARANTEE)` with a
  `// PANIC-FREE:` note citing the now-enforced invariant (a broken guarantee surfaces loudly
  rather than as a silent replacement char).
- **Review suggestion used verbatim?:** No (validate-or-escape → reject-with-typed-error, the
  crate's loud-failure style and the M5 precedent).
- **Adaptation:** validation at the merge boundary (covers ref tract + observed seqs once);
  case-insensitive ACGTN so soft-masked tracts are not rejected.
- **Verification performed:** test-first — `merge_rejects_a_non_nucleotide_observed_sequence`
  asserts the exact `InvalidAlleleByte { chrom: "chr1", start: 16, byte: b'X' }`; a unit test
  pins `first_non_acgtn` (accepts upper/lowercase ACGTN; flags tab, IUPAC `R`, `0xFF`).
- **Files changed:** `src/ssr/cohort/merge.rs`, `src/ssr/cohort/vcf_out.rs`.
- **Tests added or modified:** `merge::tests::{first_non_acgtn_accepts_nucleotides_and_flags_others,
  merge_rejects_a_non_nucleotide_observed_sequence}` (new).
- **Validation:** `cargo test --lib ssr::cohort::merge` → 0, `13 passed`; fmt/clippy/doc → 0.
- **Follow-up:** None — Mi13 was the remaining untrusted-output gap after M5.
- **Residual risk:** Low — accepts the A/C/G/T/N alphabet (incl. soft-masked); a stricter VCF
  4.4 allele grammar (e.g. rejecting `N` in ALT) could be layered later if needed.

### Mi14 — `QUAL=.` for a variable-but-zero locus
- **Severity:** Minor
- **Initial decision:** Ask → Apply (user chose numeric)
- **Final status:** Applied
- **Reasoning:** A locus emitted *because* variable but whose `site_qual` proxy clamped to
  `0.0` printed `QUAL=.`, which downstream QUAL filters read as "missing" rather than
  "lowest." `site_qual` always returns a finite value in `[0, qual_cap]`, so there is no
  genuinely-unscored case on this path. The user chose always-numeric.
- **Implementation summary:** `format_vcf_record`'s `qual_field` conditional (`if qual > 0.0
  { numeric } else { "." }`) is now an unconditional `{qual:.1}`, with a comment that `.` is
  reserved for a genuinely unscored site this path never produces.
- **Review suggestion used verbatim?:** Yes (the "always print numeric" branch of the review's
  either/or).
- **Adaptation:** None.
- **Verification performed:** updated `formats_a_no_call_filtered_record` to assert
  `cols[5] == "0.0"` (was `"."`) — the pin that a clamped-zero QUAL prints numerically.
- **Files changed:** `src/ssr/cohort/vcf_out.rs`.
- **Tests added or modified:** `formats_a_no_call_filtered_record` assertion updated.
- **Validation:** `cargo test --lib ssr::cohort::vcf_out` → 0, `10 passed`; fmt/clippy → 0.
- **Follow-up:** None.
- **Residual risk:** None.

### Mi16 — once-per-run `level_per_group.clone()`
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Inside the burn-in `pool.install` closure, `build_param_set`'s borrow of
  `grouped` ends before `run_cohort_em` is called, and `grouped` is unused afterward, so the
  `grouped.level_per_group.clone()` argument can move instead.
- **Implementation summary:** `grouped.level_per_group.clone()` → `grouped.level_per_group`
  (a move), with a one-line note that the borrow has ended and `grouped` is dead after the
  call. (The *other* `.clone()` at `build_param_set`'s `level_seed` stays — that one borrows
  `&grouped`.)
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** compiles (the borrow checker proves `grouped` is unused after);
  driver tests pass unchanged.
- **Files changed:** `src/ssr/cohort/driver.rs`.
- **Tests added or modified:** None.
- **Validation:** `cargo test --lib ssr::cohort` → 0; fmt/clippy → 0.
- **Follow-up:** None.
- **Residual risk:** None.

## 5. Deferred findings to carry forward

- **Mi3** — hoist `compute_data_ll`'s per-round `Vec<Vec<f64>>`/`obs_qr` allocations into reused
  scratch. Bench-gated (no SSR bench yet); deferred by user decision.
- **Mi4** — per-locus `f_present` allocation in the parallel sweep → `map_init` scratch. Same
  bench-gated rationale; deferred by the same decision.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** No — no bench under `benches/` references `ssr`/`cohort::`, so no Apply
  changed a bench-covered hot path.
- **Baseline saved:** N/A.
- **Outcome:** Skipped — no Apply touched perf-sensitive code. (Mi3/Mi4, the two findings that
  *would* be perf-sensitive, were deferred precisely because there is no bench to measure them.)

## 10. Commands run

- `cargo build --lib` (per finding)
- `cargo test --lib ssr::cohort` / `…ssr::cohort::{em,merge,vcf_out}` (per finding)
- `cargo fmt` / `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings` (final)
- `cargo test --lib --all-features` (final)
- `cargo test --all-targets --all-features` (final)
- `cargo doc --no-deps` (final)

## 11. Command results

- `cargo fmt --check` → 0
- `cargo clippy --all-targets --all-features -- -D warnings` → 0
- `cargo test --lib --all-features` → 0, `1292 passed; 0 failed; 2 ignored`
- `cargo test --all-targets --all-features` → 101 (pre-existing `psp_writer_perf` bench panic only)
- `cargo doc --no-deps` → 0

## 12. Notes

- All six Applied findings are behaviour-preserving except the two intended changes: Mi13 adds
  a new typed error for malformed allele bytes (loud failure on untrusted input), and Mi14
  changes a `QUAL=.` to `QUAL=0.0` for a clamped-zero emitted locus. Mi1/Mi2/Mi5/Mi16 are
  refactors with identical arithmetic/output — guarded by the byte-identity-across-thread-count
  tests.
- With this run, the only `ssr_call_driver_2026-06-24` findings still open are Mi3 and Mi4,
  both bench-gated by explicit user decision.
