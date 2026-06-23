# Code Review: ssr-call CohortMerger accessors (Step H2)

**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator; inline synthesis)
**Scope:** commit `f01c14e` on branch `ssr-cohort` — `CohortMerger::chromosomes()` / `sample_names()` + second-pass re-open
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** PR diff, commit `f01c14e` (Step H2).
- **In-scope files:** [merge.rs](../../../../src/ssr/cohort/merge.rs) — `chromosomes` field + capture in `from_parts`, `chromosomes()` / `sample_names()` accessors, `sample_name_from_label`, 3 new tests.
- **Out of scope:** the merge core (unchanged), the consumed `ParsedChromosome` type.
- **Categories considered:** reliability, errors, naming, idiomatic, refactor_safety, extras. `unsafe_concurrency` n/a.

## 2. Verdict

**Approve-with-changes.** The accessors are correct, the chromosome table is captured before `inputs` is consumed, and the two-pass re-open is proven. No Blocker/Major. One Minor (sample-name uniqueness is not enforced — a VCF-validity concern to close where the header is built) plus two Nits.

## 3. Execution status

- `cargo fmt --check` — clean. `cargo clippy --lib -D warnings` — clean. `cargo test --lib` — **1268 passed, 0 failed, 2 ignored** (+3 H2 tests).
- Sub-agent fan-out not used (small contained diff; prior fan-outs overloaded). "Needs verification": 0.

## 4. Open questions and assumptions

1. **Sample-name uniqueness.** `sample_names()` derives the basename per input; two inputs in different directories with the same file name (`a/s.ssr.psp`, `b/s.ssr.psp`) collapse to the same name, which would emit duplicate `#CHROM` sample columns — an invalid VCF. The accessor is a pure derivation; uniqueness should be enforced where all names are known (the H3 header writer or H4 driver). See Mi1.

## 5. Top 3 priorities

1. **Mi1** — enforce sample-name uniqueness (error on collision) when building the VCF header (H3/H4).
2. **Nit** — assert `md5` in the `chromosomes()` test.
3. **Nit** — `chromosomes()` returns an owned clone; a slice would avoid the copy (called once, so negligible).

## 6. Findings

### Minor

**Mi1: [merge.rs](../../../../src/ssr/cohort/merge.rs) `sample_names` — basename derivation can collide, producing duplicate VCF sample columns.**
Two inputs with the same file name in different directories both reduce to one name; the VCF `#CHROM … <samples>` line requires unique sample IDs, so a collision yields a malformed VCF that downstream tools reject. Confidence: High (mechanical consequence of the basename rule).
*Fix:* not in `sample_names()` itself (a pure derivation). When the header is assembled (H3 `write_vcf_header` or H4 driver), check `sample_names()` for duplicates and return a typed error (e.g. a new `SsrCallError::DuplicateSampleName { name }`) before writing. Recorded as a carry-over to H3/H4.

### Nits

- [merge.rs](../../../../src/ssr/cohort/merge.rs) — the `chromosomes_and_sample_names_expose_the_header_table` test asserts `name` and `length` but not `md5`; add one `md5` assertion so a future contig-table mismapping is caught (the fixture sets `"0".repeat(32)`).
- [merge.rs](../../../../src/ssr/cohort/merge.rs) — `chromosomes()` returns an owned `Vec<ParsedChromosome>` clone; since it is called once per run a `&[ParsedChromosome]` borrow would avoid the copy, but the clone is harmless and matches `chrom_names()`'s owned-return convention. Leave unless H4 finds it awkward.
- `sample_name_from_label` won't strip a doubled/altered extension (`foo.ssr.psp.bak`); correct for well-formed inputs, no action.

## 7. Out of scope observations
None. (Duplicate *contig* names would also desync `chrom_names()`/`chromosomes()` length, but that is pre-existing malformed-input territory shared with `chrom_names()`, not introduced here.)

## 8. Missing tests to add now
- `chromosomes_carries_contig_md5` — extend the existing accessor test with `assert_eq!(chroms[0].md5, "0".repeat(32))`.
- (At H3/H4, once Mi1 is implemented) `*_errors_on_duplicate_sample_names`.

## 9. What's good
- The chromosome table is cloned from the first input **before** the `for (path, reader) in inputs` loop consumes `inputs`, reusing the already-required `first` borrow — no second pass, no lifetime gymnastics. [merge.rs from_parts](../../../../src/ssr/cohort/merge.rs)
- The re-open test asserts the *stream*, not just that `open` succeeds twice — it pins the actual two-pass contract the driver relies on.

## 10. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --lib --all-features -- -D warnings` · `cargo test --lib ssr::cohort::merge`

### Author response convention
Address Mi1 + Nits by id; answer open question 1 first.
