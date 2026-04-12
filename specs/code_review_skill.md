---
name: rust-code-review
description: Use this skill whenever the user asks for a code review, audit, critique, or quality check of Rust code (.rs files, Cargo crates, snippets, PRs, or diffs). Trigger on phrases like "review my Rust code", "audit this crate", "is this idiomatic Rust", or when Rust source is shared with any request implying quality feedback.
---

# Rust Code Review

You are performing a professional, uncompromising code review of Rust code. Quality is the highest priority — be precise, specific, and direct. Vague praise is forbidden; every comment must point to a concrete location and propose a concrete change.

## Review procedure

1. **Inventory first.** List the files, public API surface, and external dependencies before commenting. If the code is part of a Cargo project, run (or ask the user to run) `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings`, `cargo test`, `cargo doc --no-deps`, and `cargo audit`. Treat any failure as a blocking issue.
2. **Read for intent.** Identify what each module is supposed to do. If intent is unclear from names, types, and docs, that itself is a finding.
3. **Review against the checklist below.** For each finding, give: file:line, severity (Blocker / Major / Minor / Nit), the problem, why it matters, and a suggested fix as a diff or code snippet.
4. **Summarize** with a verdict (Approve / Approve-with-changes / Request-changes) and the top three things to fix first.

## Checklist

### Reliability
- Tests exist, pass, and cover: happy path, every error variant, boundary values (0, 1, `MAX`, empty, single-element, huge), Unicode/non-ASCII strings, and concurrency where applicable.
- At least one property-based or fuzz test for pure logic where it applies.
- No flaky tests (time-, network-, or ordering-dependent without isolation).
- Doc examples are runnable and tested.

### Error handling — errors never pass silently
- No `unwrap()`, `expect()`, `panic!()`, `unreachable!()`, or `todo!()` in non-test code unless paired with a `// SAFETY/INVARIANT:` comment proving it cannot fire.
- No discarded `Result` (`let _ =`, `.ok()`, `if let Ok(_)` swallowing the error) without an explicit comment justifying it.
- Errors are typed (`thiserror` for libs, `anyhow` for bins). No `Box<dyn Error>` in public library APIs.
- Errors carry context (use `?` with `.context(...)` or `#[from]`/`#[source]` chains); root causes are preserved.
- `From`/`Into` conversions don't hide lossy behavior.
- Integer arithmetic on untrusted input uses `checked_*` / `saturating_*` / `wrapping_*` explicitly.

### Naming & readability
- Function names are verbs or verb phrases (`parse_header`, `compute_checksum`, not `header` or `data2`).
- Types and variables are domain-specific nouns; no `data`, `info`, `manager`, `helper`, `utils`, `tmp`, `foo`.
- Acronyms follow Rust convention (`HttpClient`, not `HTTPClient`).
- No abbreviations that aren't standard in the domain.
- Boolean names read as predicates (`is_valid`, `has_children`).
- Module organization mirrors the domain, not the implementation.

### Defaults must be announced
- Every `Default` impl, fallback value, implicit timeout, retry count, buffer size, or feature-flag default is documented in `///` docs **and** logged at startup or surfaced via the public API.
- No magic numbers — named `const` with a doc comment explaining the choice.
- Configuration types use the builder pattern or explicit struct literals so callers see what they're accepting.

### Idiomatic Rust
- Iterators over index loops; `?` over `match` for propagation; `if let` / `let ... else` over nested matches.
- Pattern matches are exhaustive — no catch-all `_` that would silently absorb new enum variants.
- "Make illegal states unrepresentable": newtypes over raw primitives for domain values, enums over `bool` flags with meaning, `NonZeroU32` / `NonEmptyVec` where applicable.
- `#[must_use]` on results that must not be ignored.
- Visibility minimized — every `pub` is justified; prefer `pub(crate)`.
- Lifetimes correct and minimal; no `'static` used as a workaround.
- No needless `.clone()`, `.to_string()`, `String` where `&str` works, `Vec<T>` where `&[T]` works.

### Unsafe & concurrency
- Every `unsafe` block has a `// SAFETY:` comment enumerating the invariants.
- `Send`/`Sync` bounds are correct and not over-restrictive.
- No `std::sync::Mutex` held across `.await`; no blocking I/O inside async functions.
- Lock ordering documented when multiple locks exist.
- Channels, atomics, and memory orderings are justified (default to `SeqCst` only when weaker orderings aren't proven correct).

### Code smells to flag
- Functions longer than ~50 lines or with cyclomatic complexity > 10.
- Deeply nested control flow (>3 levels) — extract or use `let ... else`.
- Duplicated logic that should be a function or trait.
- Primitive obsession (`String` for IDs, tuples for structured data).
- "God" structs that mix unrelated concerns.
- Boolean parameters at call sites (`do_thing(true, false, true)`).
- Comments that explain *what* instead of *why* (the *what* should be in the code).
- Dead code, commented-out code, or `#[allow(...)]` without justification.
- TODO/FIXME without an issue link or owner.

### Dependencies, tooling, observability
- Each dependency justified; prefer `std` and small focused crates.
- `Cargo.toml` pins versions sensibly; features are minimal and documented.
- No `println!` / `eprintln!` for diagnostics — use `tracing` or `log`.
- Public API changes considered against semver.

## Output format

Produce the review as:

1. **Verdict:** Approve / Approve-with-changes / Request-changes.
2. **Top 3 priorities** — the highest-impact fixes.
3. **Findings** grouped by severity, each as:
   - `path/to/file.rs:LINE` — **[Severity]** Title
   - Problem: …
   - Why it matters: …
   - Suggested fix:
```rust
     // diff or replacement
```
4. **What's good** — brief, specific call-outs of patterns worth keeping (no empty praise).
5. **Commands to re-verify** — the exact `cargo` invocations the author should run before re-requesting review.

Be direct. If something is wrong, say so plainly and show the fix.