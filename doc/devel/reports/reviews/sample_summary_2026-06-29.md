# Code Review: sample-summary data model (B1)

**Date:** 2026-06-29
**Reviewer:** rust-code-review skill (orchestrator + 3 consolidated category sub-agents)
**Scope:** feature B1 — `src/sample_summary/mod.rs` (per-sample summary TOML data model) + `pub mod` line in `src/lib.rs`
**Status:** Approve-with-changes

---

## 1. Scope
- Reviewed: B1 working diff on branch `tomato2-paralog-filter`.
- In-scope: `src/sample_summary/mod.rs` (new), `src/lib.rs` (one line).
- Out of scope: rest of crate; pre-existing `vcf/writer.rs` clippy lints.
- Categories (consolidated for a low-risk serde data model): reliability,
  errors, defaults, idiomatic, naming, smells, module_structure,
  refactor_safety.

## 2. Verdict
**Approve-with-changes.** Clean, well-documented, fully-validated data
model. Findings: a version under-enforcement, a version-policy doc/impl
mismatch, an error-API hardening gap, and serde test-coverage gaps.

## 3. Execution status
- `cargo test --lib sample_summary` → 7 passed (pre-fix).
- `cargo fmt --check` → exit 0. `cargo clippy --lib` → clean on scope.

## 4. Findings

### Major
- **M1 (reliability): `validate()` accepts `version == 0`.** Only `> 1`
  was rejected; `0` is the serde/zeroed-blob default and never written.
  Fix: refuse `version == 0` too. *(Applied.)*
- **M2 (defaults): version-policy doc says "newer major refused" but the
  field is a flat `u16` with a `>` gate.** Reword to the honest strict
  policy (every bump breaking; accept `1..=VERSION`); drop "major" from
  the `UnsupportedVersion` text. *(Applied.)*
- **M3 (errors): `SampleSummaryError` not `#[non_exhaustive]`.** Adding a
  variant would be a break. Fix: `#[non_exhaustive]` on the enum (and, for
  convention parity with `ContaminationArtefactError`, the structs).
  *(Applied.)* The companion suggestion to wrap the `toml` source types
  is **disputed**: `ContaminationArtefactError` exposes them directly and
  the crate is unpublished, so there is no external-API coupling.
- **M4 (reliability): no property/round-trip test for the serde pair.**
  The checklist requires proptest for serializers. *(Applied — proptest
  over random valid documents incl. boundary values.)*

### Minor
- **Mi1 (errors): UTF-8 failure dropped its `Utf8Error` source and was
  mislabelled `InvalidField`.** Fix: dedicated `NotUtf8 { source }`.
  *(Applied.)*
- **Mi2 (reliability): untested boundary/error branches** (zero
  gc_bins/window_bp/depth_bins; invalid-UTF-8; malformed-TOML;
  parse-path value-invariant). *(Applied — tests added.)*
- **Mi3 (naming/errors nit): `Serialize`/`Deserialize` variants are
  mechanism-named.** *(Applied — `SerializeToml`/`ParseToml`.)*
- **Mi4 (smells): the `bad` closure is duplicated across two `validate`
  bodies (2×, below the ≥3 extract threshold).** *(Deferred — revisit if
  B2/B3 add a third.)*

### Cross-category / deferred
- `to_string_pretty` over a large `counts` `Vec<u32>` and the absence of
  an upper bound on `gc_bins`/`depth_bins` (efficiency). **Deferred —**
  already bounded upstream by A1's 64 MiB decompressed-section cap.
- `sample_summary/` is a directory with only `mod.rs` today. Justified by
  the imminent B2/B3 accumulator submodules; no action.

## 9. What's good
- Invariants enforced on both serialise and parse paths via one
  `validate()`; the doc-comment states the contract.
- Mirrors the established `contamination_artefact.rs` serde-artefact
  pattern (derive + `to_string_pretty`/`from_str` + thiserror + validate).
- Kebab-case wire keys pinned by a test (incl. snake_case leak guard).

## 10. Commands to re-verify
- `./scripts/dev.sh cargo test --lib sample_summary`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings -A clippy::doc_lazy_continuation`
