---
name: rust-code-review
description: Use this skill whenever the user asks for a code review, audit, critique, or quality check of Rust code (.rs files, Cargo crates, snippets, PRs, or diffs). Trigger on phrases like "review my Rust code", "audit this crate", "is this idiomatic Rust", or when Rust source is shared with any request implying quality feedback.
---

# Rust Code Review

You are performing a professional, uncompromising code review of Rust code. Quality is the highest priority — be precise, specific, and direct. Vague praise is forbidden; every comment must point to a concrete location and propose a concrete change.

## Review principles (must always hold)

- **Correctness over style.** Prioritize behavioral correctness, and failure transparency over formatting or personal preferences.
- **No silent assumptions.** When the code leaves something unspecified — invariants on inputs, call-site guarantees, threading context, whether a collection is sorted, whether a value can be zero or empty — do not silently pick an answer and review as if it were fact. State the assumption explicitly in the finding and lower its severity until it can be verified.
- **Evidence-first findings.** Do not invent file paths, line numbers, logs, or command results. If you cannot verify, label it as "Needs verification".
- **Actionability.** Every non-trivial finding must include a concrete fix (diff/snippet, test, or refactor step), or — if the right fix depends on intent the reviewer cannot infer — a specific question whose answer would determine the fix.
- **User-visible defaults.** Any implicit default (timeout, retry count, buffer size, fallback value, feature-flag state) must be discoverable, preferably in the public API, but if that is too cumbersome, at least in the documentation or at runtime (logged at startup, returned from a config accessor, or surfaced in --help).
- **Scope discipline**. Review what was asked. For a diff or PR, focus on changed lines and their direct callers/callees; flag pre-existing issues in untouched code separately under "Out of scope observations" rather than mixing them with diff findings. Exception: pre-existing Blocker-severity issues (security, data loss, undefined behavior) should be raised prominently regardless of scope.

## Review procedure

1. **Establish scope**. Determine whether the review covers a full crate, a diff/PR, or a snippet. State the scope at the top of the review. For diffs, identify the changed files and the direct callers/callees of changed items; these define the in-scope surface.
2. **Inventory first.** List files, public API surface, error types, concurrency primitives, and external dependencies before commenting.
3. **Run verification commands** If a real execution environment is available, run: cargo fmt --check, cargo clippy --all-targets --all-features -- -D warnings, cargo test --all-targets --all-features, cargo doc --no-deps, and cargo audit. Quote actual output verbatim. **Never simulate, paraphrase, or guess at command output** — if the commands cannot be run, list them under "Commands the author must run" and proceed with static review only. Any real failure is at least **Major**; correctness-impacting failures are **Blocker**.
4. **Determine intent before judging correctness**. Before reviewing whether the code works, establish what it is meant to do — its purpose, its contract with callers, the inputs it is built to handle, and the invariants it maintains. Correctness is meaningless without intent. If intent is unclear from names, types, tests, and docs, file it as a finding against documentation or naming — whichever would have prevented the ambiguity — at Minor severity, or Major if the ambiguity could plausibly cause misuse.
5. **Review against the checklist below.** For each finding provide: file:line, severity, confidence (High/Medium/Low), assumptions (if any), problem, impact, and concrete fix.
6. **Challenge tests.** For each non-trivial function in the changed code, identify at least one input class not covered by existing tests, name the specific bug that input class would expose, and provide the test as code.
7. **Summarize** with verdict (Approve / Approve-with-changes / Request-changes), top three priorities, and a re-verification checklist.

## Severity rubric

Each finding is filed at one of four severity levels. The italicized clause is the decision rule; the rest are characteristic examples.

- **Blocker**: *Must be fixed before merge.* Likely wrong results, data loss or corruption, unsoundness or undefined-behavior risk, security vulnerabilities (path traversal, deserialization of untrusted input, timing-leaky secret comparison, injection), swallowed errors that hide failures, or absence of tests for any code path that, if broken, would produce wrong results without panicking.
- **Major**: *Should be fixed before merge; may be deferred only with explicit written justification from the author.* High risk of regressions, poor or missing error context, ambiguous or undocumented defaults, concurrency hazards (lock held across `.await`, incorrect `Send`/`Sync` bounds, unjustified memory ordering), or API design that invites misuse.
- **Minor**: *Should be fixed soon but need not block merge.* Readability or maintainability issues with moderate long-term cost: unclear naming, oversized functions, duplication, primitive obsession, missing non-critical docs.
- **Nit**: *Optional; author's discretion.* Small style or consistency issues that do not materially affect behavior.

### Severity interacts with confidence

A finding may be filed at **Blocker** only at **High** confidence. Lower-confidence findings about Blocker-class issues are filed at **Major** and paired with the specific verification step (test, call-site audit, author question) that would confirm or refute them. Once verified, the severity is upgraded in a follow-up.

### Volume guidance for Nits

Nit findings are grouped into a single "Nits" section, not enumerated as separate top-level findings. If there are more than ~5 nits of the same kind (formatting, naming convention, import ordering), do not list them individually — recommend a single mechanical fix instead (`cargo fmt`, a specific `clippy` lint, or a rename pass).

## Checklist

### Reliability

- Tests exist for every public item and every non-trivial private function.
- All tests pass under `cargo test --all-targets --all-features`.
- Test coverage includes: happy path, every error variant, boundary values (0, 1, `MAX`, empty, single-element, very large), Unicode and non-ASCII strings, malformed input, and concurrent access where shared state exists.
- Property-based test (`proptest`/`quickcheck`) or fuzz target required for: parsers, serializers/deserializers, any pure function over a structured input domain, and any function whose correctness depends on an algebraic law (associativity, idempotence, round-tripping).
- Concurrency tests required when the code uses `unsafe`, `Arc`, `Mutex`/`RwLock`, atomics, channels, or shared `async` state.
- No flaky tests. Time-dependent tests use an injected clock; network-dependent tests use a mocked transport or are gated behind a feature flag excluded from default CI; randomness uses a seeded RNG with the seed printed on failure; ordering-dependent tests use explicit synchronization, not `sleep`.
- Doc examples compile and run as doc tests. Use of `no_run` or `ignore` requires a comment explaining why.
- Every doc-comment invariant, every previously-fixed bug, and every `unsafe` safety condition has a named regression test that would fail if the invariant were violated.
- Test names describe the behavior under test and the expected outcome (`parse_returns_error_on_empty_input`, not `test_parse_2`). The name alone should communicate the bug on a CI failure.

### Error handling — errors never pass silently

- No `unwrap()` or `expect()` in non-test code unless paired with a `// PANIC-FREE:` comment naming the invariant that guarantees the call cannot fire. Prefer restructuring to eliminate the call entirely.
- No `panic!()` or `unreachable!()` in non-test code unless the branch is genuinely unreachable given the type system, with a `// UNREACHABLE:` comment explaining why. Prefer making the state unrepresentable.
- No `todo!()` or `unimplemented!()` in code submitted for review.
- No discarded `Result` or `#[must_use]` value (`let _ = …`, `.ok()` without recovery, `if let Ok(_) = …` ignoring the error) without a comment naming the specific failure mode being discarded and why it is safe to ignore.
- Errors are typed (`thiserror` for libs, `anyhow` for bins). No `Box<dyn Error>` in public library APIs. Error enums are scoped to the module or operation that produces them, not crate-wide; catch-all `Other(String)` variants require justification.
- Errors carry context. Each `?` across a meaningful boundary adds context identifying *what was being attempted* and *which input caused it*. Root causes are preserved via `#[source]` or `#[from]`, never flattened into a string.
- `From`/`Into` conversions are total and lossless. Lossy or fallible conversions use `TryFrom`/`TryInto`.
- Integer arithmetic on values derived from untrusted input uses `checked_*`, `saturating_*`, or `wrapping_*` explicitly. Plain operators are only exceptionally acceptable on values whose range is provably bounded by the type system; in release builds they wrap silently, which is data corruption, not a panic.
- Every recovery point — retry, fallback, default branch, caught-and-ignored error — emits a `tracing` event at `warn` or higher with structured fields (error chain, causing input, recovery action). Silent recovery is silent failure.
- `Drop` implementations do not panic. Failures during cleanup are logged via `tracing` and swallowed.

### Naming & readability

- Function names are verbs or verb phrases (`parse_header`, `compute_checksum`). Pure accessors returning a field-like value are an exception and use the noun (`len`, `name`, `as_bytes`) per Rust convention. Avoid placeholders (`header`, `data2`, `do_stuff`).
- Types, fields, and variables name what they *are* in the domain, not their generic shape. Forbidden when used alone: `data`, `info`, `manager`, `handler`, `helper`, `util`/`utils`, `tmp`, `foo`/`bar`, `thing`, `object`, `item`, `value`, `result`. Acceptable as suffixes when the prefix carries the meaning (`parse_result`).
- Acronyms follow Rust convention: `HttpClient`, not `HTTPClient`; `parse_url`, not `parse_URL`.
- No abbreviations except (a) Rust standard conventions (`len`, `ctx`, `buf`, `cfg`), (b) acronyms universally known in the crate's domain (`HTTP`, `AST`), or (c) variables whose lifetime is bounded to ~3 lines (`i` in a loop). Domain-specific abbreviations require a glossary entry in the crate-level docs.
- Boolean variables, fields, and `-> bool` functions read as predicates: `is_valid`, `has_children`, `contains_key`, `should_retry`. No noun-form booleans (`validity`, `children`, `retry`).
- The same concept uses the same name everywhere in the crate. `user_id` in one signature is `user_id` in all signatures — not `uid`, `user`, or `id`.
- Module names reflect the concept they own (`parser`, `connection`, `tariff`), not the layer or pattern (`models`, `services`, `helpers`, `core`, `common`). A reader scanning `src/` should be able to predict each module's contents from its name alone.
- Names are as long as they need to be to be unambiguous in their scope, and no longer.

### Defaults must be announced

The principle: a default a caller cannot see at the call site is not announced, no matter how well-documented it is elsewhere. The API itself is the primary surface for default visibility; docs and logs are reinforcements, not substitutes.

- **Visible at the call site.** Functions whose behavior depends on implicit values take an explicit configuration parameter. Never bury a behaviorally-significant default behind a no-argument constructor: `Client::new(url)` hides the timeout; `Client::builder(url).build()` or `Client::new(url, ClientConfig::default())` exposes it.
- **`Default` is discouraged.** A type may implement `Default` only if (a) the defaults are behaviorally inert (zero, empty, identity) or (b) an equivalent explicit constructor exists whose docs enumerate every default chosen. `#[derive(Default)]` on a config struct whose fields drive runtime behavior is a finding.
- **No hidden `Option` fallbacks.** A parameter typed `Option<T>` that silently falls back to a non-`None` default is a hidden default. Either take `T` directly, take a documented config struct, or name the parameter to reveal the fallback and document the default in a `# Defaults` section of the function's `///` docs.
- **Inspectable at runtime.** Every effective default is recoverable from a running instance via a public accessor or structured dump (`Debug`, `Serialize`, startup `tracing` event listing all resolved values). "Read the source" is not an acceptable answer to "what timeout is this using?".
- **Logged when applied, not only at startup.** When a default is *used* — including lazy, per-request, and fallback defaults — emit a `tracing` event at `debug` or higher identifying the field, the value, and why no explicit value was provided.
- **No magic numbers.** Each default lives in a named `pub const` with a doc comment explaining the choice, the units, and the source (spec, measurement, vendor recommendation). Public API docs reference the constant by name (`/// Defaults to [`DEFAULT_TIMEOUT`]`) so the doc and the value cannot drift.
- **Builder visibility.** Configuration types use the builder pattern or explicit struct literals. `.build()` methods document every field filled from defaults; the builder's `Debug` impl distinguishes explicitly-set fields from defaulted ones.
- **Environment-dependent defaults.** Document precedence (CLI > env > config file > built-in), test it (one test per layer plus one per pairwise override), and emit a startup `tracing` event recording each effective value and the source it came from. An operator must be able to answer "where did this value come from?" from logs alone.

### Idiomatic Rust

- Iterators over index loops; `?` over `match` for propagation; `if let` / `let ... else` over nested matches.
- Pattern matches on internal enums are exhaustive — no `_` arm that would silently absorb a new variant. For `#[non_exhaustive]` external enums or primitives where exhaustive matching is impossible, the catch-all includes a `// REVIEW ON UPGRADE:` comment.
- **Make illegal states unrepresentable.** Triggers a finding when: (a) a function takes a primitive with domain meaning (use a newtype: `Cents`, `UserId`); (b) two `bool` parameters at the same call site (use an enum); (c) a struct's fields are only valid in certain combinations (split into an enum); (d) a return type allows an "empty" value where empty is invalid (use `NonEmpty<T>` or a constructor-validated newtype); (e) a numeric field has a documented range (use a newtype with a fallible constructor).
- `#[must_use]` is required on builders, status-check functions, RAII guards, transaction handles, and any type whose construction implies an obligation.
- Visibility is minimized: private by default; `pub(crate)` for cross-module use; `pub` only for the external API. A `pub` item with no caller outside its module, or with no `///` doc comment, is a finding.
- Lifetimes are correct and minimal: no `'static` added to silence the borrow checker (use `Arc`, owned data, or a different ownership boundary); no elidable annotations written explicitly; no lifetime parameter appearing only once in a signature.
- Function parameters use the most general bound the implementation requires: `impl AsRef<Path>`, `impl IntoIterator<Item = T>`, `impl Into<String>`. Over-specific signatures force needless conversions on callers.
- No needless allocation: `.clone()` where a borrow suffices, `.to_string()` / `.to_owned()` where `&str` works, `String`/`Vec<T>` parameters where `&str`/`&[T]` work, `format!` where `Display` or `write!` would do, `.collect::<Vec<_>>()` followed by another iteration, `Box<T>` where `T` fits by value.
- `Copy` is implemented only on small (≤ ~16 bytes), value-like types where bitwise copy is semantically meaningful. Each `.clone()` at a call site is justified by ownership or replaced with a borrow.
- Imports are grouped (`std`, external, internal) and individual `use` statements are preferred over globs (except preludes and test modules). Re-exports are intentional.

### Unsafe & concurrency

- `unsafe` is a last resort. Each `unsafe` block requires a justification explaining why the safe alternative is insufficient (benchmarked perf gap, FFI, or an `std` API that is itself `unsafe`). "Faster" without a measurement is a finding.
- Every `unsafe` block has a `// SAFETY:` comment that names (a) every precondition the called function requires, (b) the specific reason each precondition holds here, referencing the code that establishes it, and (c) for `unsafe fn` and `unsafe impl`, the caller obligations — documented in `///` docs as `# Safety`. "Pointer is valid" is not a safety comment.
- `Send`/`Sync` are correct and minimal. Manual `unsafe impl Send`/`Sync` requires a safety comment proving the type contains no trait-violating shared mutability. Auto-derived `Send`/`Sync` on types containing `*mut T`, `Cell`, `RefCell`, or FFI handles is a finding — explicit, justified impls (positive or negative) are required.
- No `static mut`. Mutable globals use `Mutex`, `RwLock`, `OnceLock`, or atomics. Never roll your own initialization with `unsafe`.
- Inside `async` code: no non-`Send` guard (`std::sync::MutexGuard`, `RwLockReadGuard`/`WriteGuard`, `RefMut`, `Ref`) held across `.await`. Use `tokio::sync` primitives or scope the `std` lock to a sync block that drops the guard before awaiting. No blocking I/O (`std::fs`, `std::net`, `std::thread::sleep`, blocking drivers) — use the runtime's primitives or `spawn_blocking`.
- When multiple locks are acquired anywhere in the crate, the crate-level docs declare a total ordering. Every acquisition site is reviewed against it — inconsistent ordering across call sites is a deadlock and a Blocker. Prefer designs that hold one lock at a time.
- Atomics: every memory ordering is justified in a comment naming the synchronization relationship (which load pairs with which store, what data is protected). Default `SeqCst`; weaker orderings require explicit argument. `Relaxed` used for synchronization (rather than independent counters) is a Blocker until proven correct.
- Channel choice is justified: bounded vs. unbounded (unbounded is a memory hazard under load and requires justification), `mpsc` vs. `broadcast` vs. `oneshot`, sync vs. async. Sender and receiver drop semantics are handled on both ends.
- Thread and async task panics are not silently ignored. `JoinHandle::join` (threads) and `JoinHandle` / spawned-task results (async) are awaited and surfaced as typed errors. Detached spawns require justification.
- Functions used in `tokio::select!` arms or that may be cancelled are documented as cancellation-safe or not. State that must survive cancellation is held to allow resumption or rollback.

### Code smells to flag

Code smells are heuristics, not absolute rules — each one is a prompt to look, not a guaranteed finding. When a smell is flagged, the reviewer proposes a concrete refactor or accepts the existing code with a one-line justification.

- **Long or branchy functions.** Longer than ~75 lines, or with more than ~4 distinct decision points (`if`, `match` arm, `?`, loop condition). Use `cargo clippy -- -W clippy::cognitive_complexity` to surface the worst. A long flat state machine may be clearer than its decomposed alternative — justify and move on.
- **Deeply nested control flow** (>3 levels). Extract, invert with early `return`, or use `let ... else`.
- **Duplicated logic** in three or more places. Extract to a function, trait, or macro. Two occurrences may be left with a follow-up note.
- **"God" structs.** Triggers when (a) different methods operate on disjoint subsets of fields, (b) more than ~7 fields with no clear grouping, or (c) the struct requires `Arc<Mutex<Self>>` because too many things mutate too many fields. Split along the field-disjointness boundary.
- **Primitive obsession.** `String` for IDs/emails/paths/URLs, `u64` for money, tuples for structured data. Each occurrence is a finding pointing to a newtype.
- **Boolean parameters at call sites.** `do_thing(true, false)` is unreadable. Replace with an enum or named struct fields.
- **Stringly-typed APIs.** `String`/`&str` parameters or returns representing values from a finite set. Use an enum.
- **Long `Result`/`Option` `match` chains** where `?`, `.map`, `.and_then`, `.ok_or`, or `let ... else` would read better. Conversely, deeply nested `.and_then` chains where sequential `?` would be clearer.
- **Macros where a generic function would work.** `macro_rules!` is justified for syntax functions cannot express, not for avoiding type parameters.
- **Comments explain *why*, not *what*.** The *what* belongs in code — better names, smaller functions, self-documenting types. Inverse smell: missing `///` docs on `pub` items, or invariants that live only in the author's head.
- **Dead code, commented-out code, or `#[allow(...)]` without justification.** `#[allow(...)]` requires a comment naming (a) the specific lint, (b) why suppression is correct here, (c) the condition for removal. `#[allow(dead_code)]` on genuinely unused items is a finding — delete them.
- **`TODO`/`FIXME`/`XXX`** without an issue reference, owner, and removal condition.

### Dependencies, tooling, observability

- **Dependency justification.** Each new dependency justified by (a) what it provides, (b) why writing it ourselves is not preferable (LOC, maintenance cost, correctness risk), (c) its health: last release, open issues, maintainer count, transitive dep count, `cargo audit` status. Prefer `std` first, then small focused crates, then frameworks. A dependency added for a single function is almost always a finding — vendor the function instead.
- **Version constraints.** Libraries use the loosest range that compiles (`serde = "1"`); binaries commit `Cargo.lock` and accept the loose range. Pinning a library to an exact patch version is a finding unless justified by a known incompatibility.
- **Features.** Listed in crate-level `///` docs with a one-line description each. Additive only — mutually-exclusive features are a finding. Tested in CI under `--no-default-features`, default, and `--all-features`. Used to gate optional dependencies (`dep:foo`) where possible.
- **No `println!` / `eprintln!` for diagnostics.** Use `tracing` (preferred) or `log`. `println!` is acceptable only as a CLI binary's primary output. Library code never uses it.
- **Structured logging.** `tracing` events use structured fields (`info!(user_id = %id, "request received")`), not interpolated strings. Spans wrap operations whose duration matters. Levels: `error` for actionable failures, `warn` for recovered/degraded, `info` for lifecycle, `debug`/`trace` for development.
- **Observability completeness.** Long-running services expose structured logs, metrics (counter per error variant, histogram per operation latency, gauge per pool/queue), and health/readiness endpoints. Every public-API error variant is countable in metrics — what can't be measured can't be alerted on.
- **SemVer.** Public API changes classified per [SemVer for Rust](https://doc.rust-lang.org/cargo/reference/semver.html). Breaking changes require a major bump and `CHANGELOG.md` entry. CI runs `cargo semver-checks` (or `cargo public-api`) to detect unintended breakage. Added enum variants require `#[non_exhaustive]` to remain non-breaking.
- **Crate-level policy declarations.** Crates that aim to be `unsafe`-free declare `#![forbid(unsafe_code)]`. Panic-free crates document this and enforce via `#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic)]`. Once declared, CI enforces.
- **CI enforcement.** On every PR and on `main`: `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --all-targets --all-features`, `cargo doc --no-deps` with `RUSTDOCFLAGS="-D warnings"`, `cargo audit`, `cargo deny check`. MSRV declared in `Cargo.toml` (`rust-version`) and tested in a separate CI job. Failures block merge.
- **Publishable metadata.** Crates intended for publication declare `license`, `repository`, `description`, `categories`, `keywords`. `cargo deny check licenses` enforces the allowed-license list across transitive deps.

## Additional recommendations

- **Mutation testing for critical logic.** Run `cargo mutants` on parsers, validators, billing, security boundaries, or other correctness-critical modules. Every surviving mutant is either killed by a new test or documented as semantically equivalent — surviving mutants without justification are a Major finding.
- **Malformed-input tests as first-class cases.** Wherever the code accepts data it didn't produce itself, test: invalid UTF-8, truncated input, oversized input, duplicated fields, out-of-order fields, missing required fields, present-but-null fields, integer overflow in length prefixes, deeply nested structures. Each malformed input has a test asserting the specific error variant returned, not just "an error".
- **Golden vs. snapshot tests.** Golden tests (committed `tests/golden/` files) for any output whose stability is part of an external contract — wire formats, public JSON APIs, log formats consumed downstream. Snapshot tests (`insta`) for intentionally versioned, easy-to-rebless output — error messages, debug output, generated code. Never snapshot external-contract output.
- **Performance guards on hot paths.** `criterion` benchmarks under `benches/`, each with a `// REGRESSION THRESHOLD: N%` comment. CI compares against the `main` baseline and fails on regressions beyond threshold. Benchmark inputs are realistic and documented.
- **Documentation completeness.** Every `pub` item has `///` docs with a one-line summary, `# Errors` if it returns `Result`, `# Panics` if it can panic, `# Safety` if `unsafe`, and at least one `# Examples` block that runs as a doc test. Libraries enforce via `#![deny(missing_docs)]`.
- **Diff matches stated intent.** The PR description and the diff agree. Unrelated refactors mixed into a "fix bug X" PR are a Minor finding (request a split) and a real risk that one of them is the actual bug.
- **Release checklist for public crates.** `RELEASING.md` covers: `cargo semver-checks` clean against the previous version, `CHANGELOG.md` updated (`Added`/`Changed`/`Deprecated`/`Removed`/`Fixed`/`Security`), deprecation warnings added at least one minor version before removal, MSRV changes announced, migration guide for breaking changes.
- **Supply-chain checks beyond `cargo audit`.** `cargo deny` for licenses, banned crates, and duplicate versions. `cargo outdated` reviewed quarterly. New dependencies vetted for typo-squatting and maintenance status (last release > 18 months requires justification).
- **Reproducible builds (where applicable).** Pinned toolchain via `rust-toolchain.toml`, `Cargo.lock` committed for binaries, no network access in `build.rs` beyond `cargo`'s own.

## Output format

Produce the review in the following order. Use the section headings verbatim so the format is machine-readable.

### 1. Scope
- What was reviewed: full crate / PR diff / snippet.
- Reviewed against: commit hash, branch, or "as-provided".
- In-scope files (list).
- Deliberately out of scope (list, with reason).

### 2. Verdict
Approve / Approve-with-changes / Request-changes.

### 3. Execution status
- Commands run, with exit code and one-line result.
- Commands not run, and why.
- Count of findings labeled "Needs verification".

### 4. Open questions and assumptions
Numbered. Each entry references the findings it affects. The author resolves these before responding to individual findings.

### 5. Top 3 priorities
The highest-impact fixes, with one-line rationale and pointer to the full finding.

### 6. Findings
Grouped by severity (**Blocker** → **Major** → **Minor** → **Nits**). Within each severity, ordered by confidence (High → Low), then by file. Nits collected into a single sub-section, not enumerated.

Each finding:
- `path/to/file.rs:LINE` — **[Severity]** Title
- **Confidence:** High / Medium / Low
- **Assumptions (if any):** …
- **Problem:** …
- **Why it matters:** …
- **Suggested fix:** unified diff, complete replacement snippet, or numbered refactor procedure. Self-contained — applicable without reading the rest of the review.
```rust
  // diff or replacement
```

### 7. Out of scope observations
Pre-existing issues in untouched code, surfaced but not blocking. Each: file, brief description, suggested follow-up (separate PR or issue). Pre-existing Blocker-severity issues (security, data loss, UB) appear under Findings instead, marked "pre-existing".

### 8. Missing tests to add now
Each: proposed test name in `function_returns_expected_on_condition` form, input class covered, specific bug it would catch, and the test as code or specification. Grouped by function under test.

### 9. What's good
Up to 5 specific, transferable patterns worth keeping, each one sentence with a file reference. No general praise. Skip the section entirely if nothing specific qualifies.

### 10. Commands to re-verify
- Commands the reviewer ran (re-run to confirm they still pass).
- New commands or test invocations the review introduced.

### Author response convention
Address each finding by its identifier (e.g., "B2", "M5") with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer open questions from section 4 first.

---

Be direct. If something is wrong, say so plainly and show the fix. Vague praise and vague criticism are equally useless.

## Reusable prompt template

Use this to invoke the skill consistently. Fill in the Context block; the rest defers to the skill body.

> Perform a Rust code review per the **rust-code-review** skill. Follow its principles, procedure, severity rubric, and output format in full — do not abbreviate or skip sections.
>
> **Context**
> - **Scope:** <files / module / PR diff / branch comparison>
> - **Domain intent:** <one-paragraph description of what this code is meant to do>
> - **Audience:** <internal service / public library / CLI / embedded>
> - **Constraints:** <performance budgets, MSRV, `no_std`, target platforms, deadline pressure>
> - **Out of scope:** <anything the reviewer should not comment on — legacy modules being deleted, generated code, vendored deps>
> - **Prior review history:** <previously reviewed? known tracked issues?>
>
> **Anti-hallucination contract.** Quote tool output verbatim. If a command was not run, list it under "Commands not run". If a file or line cannot be located, say so. Never invent file paths, line numbers, error messages, clippy warnings, test results, or behavior. Findings without verifiable evidence are labeled "Needs verification" per the skill.
>
> **Reminders of the most-violated rules** (not a substitute for the full skill):
> 1. Reliability first — verify behavior, not style.
> 2. Errors must never pass silently.
> 3. Tests cover edge cases and every error path.
> 4. Names are precise, domain-relevant, verb-based for functions.
> 5. Defaults are visible at the call site, in docs, and at runtime.
> 6. Code smells get concrete refactors, not vague complaints.

---

## Saving Code Review Reports

### Directory Structure

All code review reports are saved to the `reviews/` directory at the crate root:
```
project-root/
├── src/
├── tests/
├── reviews/
│   ├── gvcf_parser_2026-04-13.md
│   ├── genotype_merging_2026-04-13.md
│   └── ...
├── Cargo.toml
└── ...
```

### File Naming Convention

Save each review report as:
```
reviews/<module_or_filename>_<YYYY-MM-DD>.md
```

Examples:
- `reviews/gvcf_parser_2026-04-13.md` — review of `src/gvcf_parser.rs` and its test file
- `reviews/genotype_merging_2026-04-13.md` — review of `src/genotype_merging.rs`
- `reviews/lib_changes_2026-04-13.md` — review of multiple files or a PR diff

If multiple reviews are created on the same date, append a version suffix:
- `reviews/gvcf_parser_2026-04-13_v1.md`
- `reviews/gvcf_parser_2026-04-13_v2.md`

### Report Format

Each review report is a Markdown file structured exactly as follows (use section headings verbatim):

#### 1. **Document Header**
```markdown
# Code Review: <module-or-file-name>
**Date:** <YYYY-MM-DD>  
**Reviewer:** <name or "GitHub Copilot">  
**Module:** <module-name>  
**Status:** <Approve / Approve-with-changes / Request-changes>

---
```

#### 2. **Sections (in order, use headings verbatim)**
1. **Scope** — describe what was reviewed
2. **Verdict** — one-line summary of approval status
3. **Execution Status** — table of commands run and results
4. **Open Questions & Assumptions** — numbered list (or "none")
5. **Top 3 Priorities** — list of highest-impact fixes
6. **Findings** — grouped by severity (Blocker → Major → Minor → Nits)
7. **Out of Scope Observations** — pre-existing issues not in scope
8. **Missing Tests to Add Now** — test code or specs
9. **What's Good** — transferable patterns worth keeping
10. **Commands to Re-verify** — instructions for follow-up verification
11. **Author Response Template** (optional) — response format for the author

#### 3. **Findings Template (for each finding)**

Each finding uses this format:

```markdown
#### <SEVERITY_CODE>: [file.rs](../path/file.rs#L123–L456) — Short title

- **Confidence:** High / Medium / Low
- **Assumptions:** None. / Listed here.
- **Problem:**  
  [2–4 sentences describing the defect, with code snippets if needed.]

- **Why it matters:**  
  [1–2 sentences explaining impact: data loss, performance, correctness, maintainability.]

- **Suggested fix:**  
  [Concrete code diff, full replacement function, or numbered refactor steps.]
  
  ```rust
  // example code
  ```
```

#### 4. **Severity Codes**
Use prefixes in findings:
- `B` for Blocker (e.g., B1, B2, B3)
- `M` for Major (e.g., M1, M2)
- `Mi` for Minor (e.g., Mi1, Mi2)
- Group Nits under a single "Nits" section without enumeration

#### 5. **File Links**
All file references use relative Markdown links with line numbers:
- Link to line: `[file.rs](../path/file.rs#L123)`
- Link to range: `[file.rs](../path/file.rs#L123–L456)`
- Display text should be the path or a brief description, but not backticks

### Commands Before Saving

Before writing the review report, the reviewer must run (where applicable):

1. `cargo fmt --check` — format check
2. `cargo clippy --all-targets --all-features -- -D warnings` — lint check
3. `cargo test --all-targets --all-features` — test execution
4. `cargo doc --no-deps` — documentation build (optional for libs)
5. `cargo audit` — security audit (optional)

Quote output **verbatim**. If a command cannot be run, list it under **Execution Status** as "not run" with reason.

### Review Checklist Before Saving

- [ ] All findings reference specific line numbers with `file.rs#Lxxx` links
- [ ] Each Blocker finding has High confidence
- [ ] Each Major/Minor finding has a concrete fix (code, test, or specific question)
- [ ] Test recommendations include test names matching `test_behavior_on_condition` pattern
- [ ] "What's Good" section has 3–5 specific, transferable patterns (not vague praise)
- [ ] File paths are relative to project root (e.g., `src/file.rs`, not `/home/...`)
- [ ] Severity codes (B1, M1, etc.) are used consistently
- [ ] Open questions are numbered and linked from findings
- [ ] Author response template provided (if review blocks merge)