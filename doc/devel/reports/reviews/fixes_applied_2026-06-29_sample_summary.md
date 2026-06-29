# Fix Application Report: sample_summary_2026-06-29.md

**Date:** 2026-06-29
**Source review:** `doc/devel/reports/reviews/sample_summary_2026-06-29.md`
**Source state:** branch `tomato2-paralog-filter`, B1 working diff
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers 0; Majors 4 (M1–M4); Minors 4 (Mi1–Mi4); Nits 1 (folded into Mi3).

### Outcome totals
- Applied: 7 (M1, M2, M3, M4, Mi1, Mi2, Mi3)
- Disputed: 1 (M3's toml-wrapping half)
- Deferred: 2 (Mi4; efficiency/upper-bound cross-cat)

### Validation summary
- `cargo fmt --check` → exit 0
- `cargo clippy --lib -- -D warnings` → clean on scope (5 pre-existing
  `vcf/writer.rs` lints only)
- `cargo test --lib` → 1368 passed, 0 failed (+13 for B1, +6 vs pre-fix)
- `cargo doc --no-deps` → pre-existing `ClassicStutterModel` failure only
- `cargo audit` → not run (no dependency change)

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files |
|---|---|---|---|---|---|
| M1 | Major | `version == 0` accepted | Apply | Applied | sample_summary/mod.rs |
| M2 | Major | version-policy doc/impl mismatch | Apply | Applied | sample_summary/mod.rs |
| M3 | Major | error enum not `#[non_exhaustive]` | Apply | Applied | sample_summary/mod.rs |
| M3b | Major | wrap `toml` source types | Dispute | Disputed | None |
| M4 | Major | no proptest round-trip | Apply | Applied | sample_summary/mod.rs |
| Mi1 | Minor | UTF-8 error dropped source | Apply | Applied | sample_summary/mod.rs |
| Mi2 | Minor | untested boundary/error branches | Apply | Applied | sample_summary/mod.rs |
| Mi3 | Minor | mechanism-named variants | Apply | Applied | sample_summary/mod.rs |
| Mi4 | Minor | duplicated `bad` closure (2×) | Defer | Deferred | None |
| CC | — | large-counts / bin upper bound | Defer | Deferred | None |

## 3. Questions asked and answers
None (non-interactive).

## 4. Per-finding log (key items)

### M1 — `version == 0` accepted
Applied. `validate()` now refuses `version == 0 || version > VERSION`.
Test `validate_rejects_version_zero`.

### M2 — version-policy doc/impl mismatch
Applied. Reworded the const doc, the struct field doc, and the
`UnsupportedVersion` message + enum doc to the honest flat-counter policy
(`accept 1..=VERSION`; every bump breaking). The `>` operator (refuse
future) was correct and kept.

### M3 — `#[non_exhaustive]`
Applied to `SampleSummaryError` and the three on-disk structs (convention
parity with `ContaminationArtefactError`). In-crate construction (B2/B3,
tests) is unaffected; only external struct literals are blocked, which is
the intended hardening.

### M3b — wrap `toml` source types — **Disputed**
The crate's own convention (`ContaminationArtefactError`, lines 397–407)
exposes `toml::de::Error`/`toml::ser::Error` directly, and the crate is
not a published library, so there is no downstream version-coupling to
protect. Kept the direct `#[source]` exposure; the value (preserved cause
chain) is retained without the newtype boilerplate.

### M4 — proptest round-trip
Applied. `round_trips_arbitrary_valid_summary` generates random valid
documents (boundary `gc_bins`/`depth_bins` at 1, `u32::MAX` counts,
`lo == hi`, `min_depth = 0`) and asserts `from_toml_bytes(to_toml_bytes)`
round-trips.

### Mi1 / Mi2 / Mi3
- Mi1: new `NotUtf8 { source: Utf8Error }` variant; `from_toml_bytes`
  maps to it (cause preserved). Test `from_toml_bytes_rejects_invalid_utf8`.
- Mi2: `validate_rejects_zero_bin_dimensions`,
  `from_toml_bytes_rejects_malformed_toml`,
  `from_toml_bytes_rejects_value_invariant_violation`.
- Mi3: `Serialize`/`Deserialize` → `SerializeToml`/`ParseToml`.

## 5. Deferred findings to carry forward
- **Mi4** — extract the `bad` closure to a shared `invalid_field`
  constructor once a third `validate` body (B2/B3) makes it ≥3 uses.
- **CC** — explicit upper bound on `gc_bins * (depth_bins + 1)`; deferred
  because A1's 64 MiB decompressed-section cap already bounds the parsed
  matrix size.

## 6. Disputed findings to return to reviewer
- **M3b** — wrapping `toml` errors. Rationale above.

## 7. Failed-validation findings
None.

## 9. Performance check
Skipped — no `Apply` touched a `benches/`-covered hot path (TOML
serde of a per-sample summary, once per file).

## 10. Commands run
- `./scripts/dev.sh cargo test --lib sample_summary`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings [-A clippy::doc_lazy_continuation]`

## 11. Command results
- `cargo test --lib` → 1368 passed, 0 failed
- `cargo fmt --check` → exit 0
- `cargo clippy` → 5 pre-existing vcf/writer.rs errors only

## 12. Notes
- M3b dispute is the only non-applied Major; it is a deliberate
  convention-match, not a deferral of work.
