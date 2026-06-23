# Code Review: ssr-call VCF header writer (Step H3)

**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator; inline synthesis)
**Scope:** commit `2ff5b1b` on branch `ssr-cohort` — `write_vcf_header` + `FILTER_DESCRIPTIONS`
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** PR diff, commit `2ff5b1b` (Step H3).
- **In-scope files:** [vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs) — `write_vcf_header`, `FILTER_DESCRIPTIONS`, the `std::io::Write` + `ParsedChromosome` imports, 3 new tests.
- **Out of scope:** `format_vcf_record` (the data-line writer, unchanged); the sample-name uniqueness check (H4, where `SsrCallError` lives — avoids a `vcf_out → driver` back-reference).
- **Categories considered:** reliability, errors, naming, idiomatic, defaults, extras (stable-output). `unsafe_concurrency`/`module_structure`/`tooling` n/a.

## 2. Verdict

**Approve-with-changes.** A correct, dedicated SSR header writer. The FILTER IDs are sourced from `filter_text` (single source of truth with the data lines), `REPCN` is `Number=.` `Integer` to match the comma-joined units `format_vcf_record` emits, and the column line is built correctly. No Blocker/Major/Minor — only Nits about documenting the input contract.

## 3. Execution status

- `cargo fmt --check` clean · `cargo clippy --lib -D warnings` clean · `cargo test --lib` = **1271 passed, 0 failed, 2 ignored** (+3 H3 tests).
- Sub-agent fan-out not used (small formatting diff). "Needs verification": 0.

## 4. Open questions and assumptions

1. **Header-string inputs are assumed VCF-clean.** `write_vcf_header` interpolates `warnings`, contig names, and sample names verbatim. A warning containing a newline, or a sample/contig name containing a tab or `"`, would produce a malformed header. In practice warnings come from `f_is_warning` (single line), contig names from the `.ssr.psp` header, and sample names from path basenames — all VCF-clean — and sample-name *validity/uniqueness* is H4's check. So this is a documented contract, not a live bug. See Nits.

## 5. Top 3 priorities

Nits only — see below. The material follow-up (sample-name uniqueness) is already tracked for H4.

## 6. Findings

### Nits

- [vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs#L81-L90) — the writer assumes each `warning` is a single line (a newline would split it into a bogus meta line). True for `f_is_warning`, but worth a one-line doc note on `write_vcf_header` stating the single-line contract (and that VCF-validity of sample/contig names is the caller's responsibility, enforced in H4).
- [vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs#L91-L93) — `##contig` emits `ID` + `length` only; the `md5` carried by `ParsedChromosome` is dropped. This matches arch §5 (`ID, length`) and is a fine default, but the md5 is useful reference provenance and is in hand — note it as a deliberate omission so a future reader knows it was a choice, not an oversight.
- Description strings and the `Number=.`/`Number=1` choices are correct; no escaping needed for the static Descriptions. No action.

## 7. Out of scope observations
None.

## 8. Missing tests to add now
None required — the three tests cover all sections, the FILTER/`filter_text` consistency (and PASS omission), and the no-warning path. (A multi-line-warning guard test would only be warranted if a guard is added; not necessary given the documented contract.)

## 9. What's good
- `FILTER_DESCRIPTIONS` pairs each `Admission` with its description and pulls the ID from `filter_text`, so the header declarations and the per-record FILTER values cannot drift — and a test asserts exactly that. [vcf_out.rs:55-69](../../../../src/ssr/cohort/vcf_out.rs#L55-L69)
- `REPCN` declared `Number=.` `Type=Integer` correctly matches the comma-joined units the data line emits, rather than a single opaque String.
- Generic `W: Write` keeps the writer testable in memory (the tests write to a `Vec<u8>`), no file needed.

## 10. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --lib --all-features -- -D warnings` · `cargo test --lib ssr::cohort::vcf_out`

### Author response convention
Address the Nits by id; answer open question 1 first (it drives the doc-note Nit).
