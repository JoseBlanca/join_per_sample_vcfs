# Additional recommendations

**Purpose.** Higher-bar techniques that aren't always justified but pay off on critical-path code, public crates, or releases. Apply per item — not all-or-nothing.

**Triggers.** Look for: parser/validator/security-boundary code (mutation testing); code that accepts untrusted input (malformed-input tests); code producing stable output (golden vs snapshot); hot paths (perf guards); public items missing docs; PR description mismatch with diff; release-related changes.

**Skip when.** None of the above apply to the in-scope code. The orchestrator may dispatch this category just for "Diff matches stated intent" if the scope is a PR.

## Rules

- **Mutation testing for critical logic.** Run `cargo mutants` on parsers, validators, billing, security boundaries, or other correctness-critical modules. Every surviving mutant is either killed by a new test or documented as semantically equivalent — surviving mutants without justification are a Major finding.
- **Malformed-input tests as first-class cases.** Wherever the code accepts data it didn't produce itself, test: invalid UTF-8, truncated input, oversized input, duplicated fields, out-of-order fields, missing required fields, present-but-null fields, integer overflow in length prefixes, deeply nested structures. Each malformed input has a test asserting the *specific error variant* returned, not just "an error".
- **Golden vs. snapshot tests.** Golden tests (committed `tests/golden/` files) for any output whose stability is part of an external contract — wire formats, public JSON APIs, log formats consumed downstream. Snapshot tests (`insta`) for intentionally versioned, easy-to-rebless output — error messages, debug output, generated code. Never snapshot external-contract output.
- **Performance guards on hot paths.** `criterion` benchmarks under `benches/`, each with a `// REGRESSION THRESHOLD: N%` comment. CI compares against the `main` baseline and fails on regressions beyond threshold. Benchmark inputs are realistic and documented.
- **Documentation completeness.** Every `pub` item has `///` docs with a one-line summary, `# Errors` if it returns `Result`, `# Panics` if it can panic, `# Safety` if `unsafe`, and at least one `# Examples` block that runs as a doc test. Libraries enforce via `#![deny(missing_docs)]`.
- **Diff matches stated intent.** The PR description and the diff agree. Unrelated refactors mixed into a "fix bug X" PR are a Minor finding (request a split) and a real risk that one of them is the actual bug.
- **Release checklist for public crates.** `RELEASING.md` covers: `cargo semver-checks` clean against the previous version, `CHANGELOG.md` updated (`Added`/`Changed`/`Deprecated`/`Removed`/`Fixed`/`Security`), deprecation warnings added at least one minor version before removal, MSRV changes announced, migration guide for breaking changes.
- **Supply-chain checks beyond `cargo audit`.** `cargo deny` for licenses, banned crates, and duplicate versions. `cargo outdated` reviewed quarterly. New dependencies vetted for typo-squatting and maintenance status (last release > 18 months requires justification).
- **Reproducible builds (where applicable).** Pinned toolchain via `rust-toolchain.toml`, `Cargo.lock` committed for binaries, no network access in `build.rs` beyond `cargo`'s own.
