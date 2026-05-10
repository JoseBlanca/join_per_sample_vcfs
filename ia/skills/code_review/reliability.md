# Reliability checklist

**Purpose.** Behavioural correctness verified by tests — coverage, robustness, regression protection, test naming.

**Triggers.** Look at every public item, every non-trivial private function, every previously-fixed bug, every `unsafe` safety condition, every doc-comment invariant. Read the existing test files for the scope and identify gaps.

**Skip when.** Never skipped — even snippets without tests get findings (the missing tests *are* the findings).

## Rules

- **Tests exist** for every public item and every non-trivial private function. Each missing test is a finding with a proposed test name in `function_returns_expected_on_condition` form.
- **All tests pass** under `cargo test --all-targets --all-features`. If the orchestrator already ran this and quoted output, do not re-run; use that output. A failing test is at least Major; a failure on a correctness-critical path is Blocker.
- **Coverage classes.** Every function under review is covered against: happy path, every error variant it can return, boundary values (`0`, `1`, `MAX`, empty, single-element, very large), Unicode and non-ASCII strings, malformed input, and concurrent access where shared state exists. Each missing class is a finding.
- **Property-based or fuzz tests** are required for: parsers, serializers/deserializers, any pure function over a structured input domain, and any function whose correctness depends on an algebraic law (associativity, idempotence, round-tripping). Use `proptest` / `quickcheck` or a `cargo-fuzz` target.
- **Concurrency tests** are required when the code uses `unsafe`, `Arc`, `Mutex`/`RwLock`, atomics, channels, or shared `async` state.
- **No flaky tests.** Time-dependent tests use an injected clock; network-dependent tests use a mocked transport or a feature flag excluded from default CI; randomness uses a seeded RNG with the seed printed on failure; ordering-dependent tests use explicit synchronization, not `sleep`. Each violation is at least Major.
- **Doc examples** compile and run as doc tests. Use of `no_run` or `ignore` requires a comment explaining why.
- **Regression tests** exist for every doc-comment invariant, every previously-fixed bug, and every `unsafe` safety condition. The test would fail if the invariant were violated.
- **Test names** describe the behavior under test and the expected outcome (`parse_returns_error_on_empty_input`, not `test_parse_2`). The name alone communicates the bug on a CI failure.

## Challenge tests (additional pass)

For each non-trivial function in the changed code, identify at least one input class **not** covered by existing tests, name the specific bug that input class would expose, and provide the test as code. Add these under a `## Missing tests` heading at the bottom of your output file (in addition to per-rule findings). Each entry:

- Proposed test name in `function_returns_expected_on_condition` form
- Input class covered (e.g. "empty input", "header longer than 4 KiB")
- Specific bug it would catch
- Test body as code (or an unambiguous specification)
