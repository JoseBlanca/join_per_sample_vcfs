---
name: rust-feature-implementation
description: Use this skill whenever the user asks to implement a new Rust feature, add new behavior, create a new module/API, or write production code beyond small fixes. Trigger on requests like "implement feature", "add support for", "build", "create", "develop", "write Rust code", or "add tests" in .rs crates.
---

# Rust Feature Implementation

Implement new Rust features with production-grade quality.

This skill is for writing new code, not just patching existing defects. It emphasizes clear design, readable domain language, correctness, and strong tests, with explicit attention to high-throughput data processing workloads.

## Priorities

- **Correctness first.** Code must satisfy functional requirements and domain invariants.
- **Tests are mandatory.** Every feature ships with meaningful tests.
- **Readable code.** Names communicate domain intent; structure is easy to follow.
- **Idiomatic Rust.** Use standard patterns and ownership/lifetime design that fit Rust best practices.
- **Performance-aware design.** Optimize where it matters, especially in hot paths that process large datasets.

## Default workflow

1. Clarify the feature contract.
2. Design the shape of the solution (types, module boundaries, data flow, errors).
3. Share a short implementation plan with the user and ask for comments when interactive mode is desired.
4. Add or update tests to represent expected behavior (prefer red/green when practical).
5. Implement the feature with minimal, focused changes.
6. Validate with Rust tooling.
7. Report what changed, why, and how it was validated.

## Interactive planning gate

Before coding, default to this behavior unless the user explicitly asks to skip it:

- Present a concise plan (scope, affected files/modules, test strategy, and risks).
- Ask for confirmation or comments.
- If feedback is provided, incorporate it and restate the adjusted plan briefly.
- If no feedback is requested by the user context, proceed after presenting the plan.

This keeps the process collaborative without blocking execution.

## Feature contract checklist

Before coding, ensure the following is known or explicitly assumed:

- Inputs and outputs (types, format, ordering, cardinality).
- Edge cases and failure modes.
- Performance constraints (latency/throughput/memory).
- API compatibility requirements (internal crate vs public API).
- Concurrency expectations (single-threaded, rayon, async, lock behavior).

If critical information is missing and changes the implementation direction, ask focused questions before coding.

## Rust design rules

- Make illegal states unrepresentable when reasonable (newtypes, enums, validated constructors).
- Keep public APIs explicit and unsurprising.
- Use typed errors (`thiserror` for libs, `anyhow` for app boundaries).
- Keep visibility minimal (`pub(crate)`/private by default).
- Prefer borrowing over cloning when ownership does not require a copy.
- Avoid hidden defaults for behaviorally significant values.
- Keep functions cohesive; split long functions by intent.

## Data-processing style guidelines

The codebase commonly processes very large datasets composed of smaller independent items.

- Prefer iterator pipelines and functional transformations (`map`, `filter`, `fold`, `try_fold`, `flat_map`) where they improve clarity.
- Use chunked/batched processing when full materialization would increase memory pressure.
- Avoid unnecessary intermediate allocations and copies.
- Mutate state only when it measurably improves performance or clarity versus a pure transform.
- For parallel processing, ensure deterministic behavior where required and avoid hidden shared mutable state.

Do not force strict functional programming. Using structs/enums plus `impl` blocks is encouraged when it improves domain modeling and maintainability.

## Naming and readability standards

- Functions should be verbs or verb phrases (`parse_header`, `merge_genotypes`).
- Types should be domain nouns (`VariantGroup`, `GenotypePosterior`).
- Avoid vague names: `item`, `data`, `value`, `thing`, `obj`, `tmp`, `helper`.
- Keep names explicit enough to be unambiguous in context.
- Prefer straightforward control flow with early returns over deep nesting.
- Add comments only for non-obvious invariants, algorithms, or performance tradeoffs.

## Testing policy (mandatory)

Every implemented feature must include tests. The exact order is flexible:

- Prefer red/green TDD when it accelerates understanding.
- If not using strict TDD, add tests immediately after implementation before considering the task done.

Minimum test coverage expectations:

- Happy path behavior.
- Relevant edge cases and boundary values.
- Error paths and invalid/malformed input handling.
- Regression test for any bug discovered during implementation.

When appropriate, add:

- Property-based tests for parser/transformer logic.
- Integration tests for end-to-end feature behavior.
- Benchmarks for suspected hot paths (`benches/` with `criterion`) when performance is a requirement.

## Performance and scalability checks

For data-intensive features:

- Estimate algorithmic complexity in practical terms.
- Check memory behavior (streaming vs buffering).
- Validate that I/O and parsing are not accidentally serialized if parallelism is intended.
- Avoid micro-optimizations without evidence.

Only add low-level optimization complexity when supported by profiling or benchmark results.

## Validation commands

Run and report real command results when environment allows:

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo test <targeted_test_name>` for newly added scenarios

If a command cannot be run, state that explicitly and do not fabricate output.

## Definition of done

A feature is complete only when all of the following are true:

- Plan was shared up front, and user feedback was incorporated when provided.
- Requirements are implemented.
- Code is readable and domain-named.
- New tests exist and pass.
- Existing tests pass (or failures are explained and unrelated).
- Formatting/linting expectations are met or explicitly documented.
- Any tradeoffs (performance, API decisions, deferred work) are clearly documented.

## Response format for implementation tasks

When using this skill, structure the response in this order:

1. **Plan**: short implementation plan.
2. **Changes made**: files and behavior added/modified.
3. **Tests added/updated**: what each new test validates.
4. **Validation results**: commands run and outcomes.
5. **Tradeoffs and follow-ups**: explicit non-goals, deferred improvements.

## Reusable prompt template

Use this template to invoke the skill consistently:

> Implement this feature using the **rust-feature-implementation** skill.
>
> **Feature request**
> - Goal: <what should be added>
> - Context: <domain/module>
> - Inputs/Outputs: <contracts>
> - Constraints: <performance, memory, API, MSRV>
> - Testing expectations: <unit/integration/property/bench>
> - Out of scope: <what should not be changed>
>
> Requirements:
> 0. Present a short plan first and ask for comments before implementation (unless explicitly skipped).
> 1. Write readable, idiomatic Rust.
> 2. Use domain-specific naming.
> 3. Favor iterator/functional approaches for item-wise processing where clear.
> 4. Include tests that prove behavior and edge cases.
> 5. Report commands actually run and real outcomes.