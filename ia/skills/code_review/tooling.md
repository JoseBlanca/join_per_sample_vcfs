# Dependencies, tooling, observability

**Purpose.** Project-level discipline that spans the crate: which dependencies to take, which lints to enforce, what to log, what CI guarantees, what publishable metadata is required.

**Triggers.** Read `Cargo.toml`, the `[lints]` table, `lib.rs` / `main.rs` crate-level attributes, every `println!` / `eprintln!` / `tracing` call, every CI configuration file in the repo. Check for `forbid(unsafe_code)`, `deny(missing_docs)`, the `[lints.clippy]` table.

**Skip when.** The scope is a snippet without a `Cargo.toml` or known crate context. A few rules still apply (e.g. `println!` discipline, structured logging) — apply only those if the orchestrator dispatches anyway.

## Rules

- **Dependency justification.** Each new dependency justified by (a) what it provides, (b) why writing it ourselves is not preferable (LOC, maintenance cost, correctness risk), (c) its health: last release, open issues, maintainer count, transitive dep count, `cargo audit` status. Prefer `std` first, then small focused crates, then frameworks. A dependency added for a single function is almost always a finding — vendor the function instead.
- **Version constraints.** Libraries use the loosest range that compiles (`serde = "1"`); binaries commit `Cargo.lock` and accept the loose range. Pinning a library to an exact patch version is a finding unless justified by a known incompatibility.
- **Features.** Listed in crate-level `///` docs with a one-line description each. Additive only — mutually-exclusive features are a finding. Tested in CI under `--no-default-features`, default, and `--all-features`. Used to gate optional dependencies (`dep:foo`) where possible.
- **No `println!` / `eprintln!` for diagnostics.** Use `tracing` (preferred) or `log`. `println!` is acceptable only as a CLI binary's primary output. Library code never uses it.
- **Structured logging.** `tracing` events use structured fields (`info!(user_id = %id, "request received")`), not interpolated strings. Spans wrap operations whose duration matters. Levels: `error` for actionable failures, `warn` for recovered/degraded, `info` for lifecycle, `debug`/`trace` for development.
- **Observability completeness.** Long-running services expose structured logs, metrics (counter per error variant, histogram per operation latency, gauge per pool/queue), and health/readiness endpoints. Every public-API error variant is countable in metrics — what can't be measured can't be alerted on.
- **SemVer.** Public API changes classified per [SemVer for Rust](https://doc.rust-lang.org/cargo/reference/semver.html). Breaking changes require a major bump and `CHANGELOG.md` entry. CI runs `cargo semver-checks` (or `cargo public-api`) to detect unintended breakage. Added enum variants require `#[non_exhaustive]` to remain non-breaking.
- **Crate-level policy declarations.** Crates that aim to be `unsafe`-free declare `#![forbid(unsafe_code)]`. Panic-free crates document this and enforce via `#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic)]`. Once declared, CI enforces.
- **CI enforcement.** On every PR and on `main`: `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --all-targets --all-features`, `cargo doc --no-deps` with `RUSTDOCFLAGS="-D warnings"`, `cargo audit`, `cargo deny check`. MSRV declared in `Cargo.toml` (`rust-version`) and tested in a separate CI job. Failures block merge.
- **Curated clippy lint set declared in the crate.** `-D warnings` only catches lints clippy already raises by default. Opt into a correctness-focused set in `[lints.clippy]` in `Cargo.toml` (Rust ≥ 1.74) or at the crate root. Recommended baseline:
  - `indexing_slicing = "deny"` — pairs with the *Refactor safety* slice-matching rule.
  - `unwrap_used = "deny"`, `expect_used = "deny"` — pair with *Error handling*.
  - `panic = "deny"` on panic-free crates.
  - `fallible_impl_from = "deny"` — pairs with the `TryFrom` rule.
  - `wildcard_enum_match_arm = "deny"` — pairs with exhaustive matching.
  - `fn_params_excessive_bools = "deny"` — pairs with the boolean-parameters smell.
  - `must_use_candidate = "warn"` — surfaces missing `#[must_use]` annotations.

  Per-call-site suppression with `#[allow(...)]` requires the justification comment already mandated under code smells.
- **Publishable metadata.** Crates intended for publication declare `license`, `repository`, `description`, `categories`, `keywords`. `cargo deny check licenses` enforces the allowed-license list across transitive deps.
