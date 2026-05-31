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

0. **Read `PROJECT_STATUS.md` at the project root.** Use the "About this project" paragraph for design context, "Current focus" to confirm the task matches the project's direction, and the in-scope feature's block (if any) to locate its spec, plan, and prior reports. See *Project status protocol* below.
1. Clarify the feature contract.
2. Design the shape of the solution (types, module boundaries, data flow, errors).
3. Share a short implementation plan with the user and ask for comments when interactive mode is desired.
4. Add or update tests to represent expected behavior (prefer red/green when practical).
5. Implement the feature with minimal, focused changes.
6. Validate with Rust tooling.
7. Report what changed, why, and how it was validated.
8. **Update `PROJECT_STATUS.md`** — append the new implementation report to the feature's block (create the block if first time touching this feature), update `Status:`, refresh `Last completed task`. See *Project status protocol* below.

## Interactive planning gate

Before coding, default to this behavior unless the user explicitly asks to skip it:

- Present a concise plan (scope, affected files/modules, test strategy, and risks).
- Ask for confirmation or comments.
- If feedback is provided, incorporate it and restate the adjusted plan briefly.
- If no feedback is requested by the user context, proceed after presenting the plan.

This keeps the process collaborative without blocking execution.

## Feature contract checklist

Before coding, ensure the following is known or explicitly assumed:

- **Domain intent.** A one-paragraph statement of what the feature is *meant* to do — its purpose, the contract it offers callers, and the invariants it maintains. Correctness claims are meaningless without intent; write intent down before judging design choices against it.
- Inputs and outputs (types, format, ordering, cardinality).
- Edge cases and failure modes.
- Performance constraints (latency/throughput/memory).
- API compatibility requirements (internal crate vs public API).
- Concurrency expectations (single-threaded, rayon, async, lock behavior).

If critical information is missing and changes the implementation direction, ask focused questions before coding.

**No silent assumptions.** When the spec leaves something unspecified and you proceed anyway — picking a default, choosing an ordering, deciding what an empty input means, assuming a collection is sorted — write the assumption down in the plan and in the final report. Silent choices are where features end up technically working but contractually wrong; surfacing them gives the user a chance to redirect before the choice is buried in code.

## Rust design rules

- Make illegal states unrepresentable when reasonable (newtypes, enums, validated constructors).
- Keep public APIs explicit and unsurprising.
- Use typed errors (`thiserror` for libs, `anyhow` for app boundaries).
- Keep visibility minimal (`pub(crate)`/private by default).
- Prefer borrowing over cloning when ownership does not require a copy.
- Avoid hidden defaults for behaviorally significant values.
- Keep functions cohesive; split long functions by intent.

## Module structure and refactoring process

These rules govern *where* code lives and *how* to move it. The corresponding review-time checklist is `code_review/module_structure.md`; the rules here are the writer-side counterpart — apply them when creating a new module, lifting a type, or restructuring the tree.

- **Crate-absolute paths over `super::super::*`.** When referencing items more than one module level up, use `crate::<module>::<item>` rather than chained `super::*`. The crate-absolute form survives file moves without breaking, reads self-evidently at the use-site, and shows up in editor go-to-definition. Single-level `super::*` is fine.

- **New top-level modules: flat `.rs` unless submodules will follow.** A leaf module belongs at `src/<name>.rs`. Create `src/<name>/mod.rs` only when the module already has (or imminently will have) siblings inside the directory. The import path `crate::<name>::*` is identical either way; the directory form just costs a directory entry for nothing.

- **Where to put a new data type.** If a type will be consumed by 3+ modules across pipeline-stage boundaries — even prospectively — define it as its own top-level peer (`src/<type>.rs`), not inside the producer's module. Otherwise the producer's module silently becomes the home for the pipeline's interchange data, and refactor pressure builds. Pre-checking at write time saves the lift later.

- **Back-reference sanity check.** Before declaring a new top-level peer module (data model, algorithm, file format, utility) done, grep its `use crate::*` lines. Any import from a pipeline-stage module is a back-reference and probably means: lift the depended-on type too, move the new module inside the stage, or decouple the consumer via a trait. The trap is *moving a directory without addressing what it imports* — the new location misrepresents the dependencies.

- **`pub` consumers of pipeline-stage types: consider a trait at the boundary.** When a `pub` module takes a concrete record type from a pipeline-stage module (e.g., a writer that consumes a specific `Record`), the natural cost is a back-reference into that stage. If the consumer's data needs are stable and read-only, a trait with read-only accessors monomorphises to zero runtime cost and removes the back-reference from the public API. Worth proposing only when (a) the consumer is `pub` and (b) the input type is stage-internal. Otherwise the trait is over-abstraction with one implementor.

- **`clippy::module_inception` when promoting a flat file.** Promoting `foo.rs` into `foo/` with submodules leaves the original content at `foo/foo.rs`, which the lint fires on. Rename the inner file to something specific (`driver.rs`, `core.rs`, `impl.rs`) as part of the same change — not as a follow-up. The new name also tends to clarify the file's role.

- **`git mv` for renames.** Use `git mv old new`, not `mv` + `git add new` + `git rm old`. Git's rename detection is heuristic and reliable only when both files are present and the content is close enough; explicit `git mv` works at any similarity level and keeps `git log --follow` cheap.

- **Bulk path substitutions: `sed` is acceptable for mechanical renames.** When a refactor touches 20+ files with a purely-mechanical path change (e.g. `crate::old_path::` → `crate::new_path::`), a single `sed -i` pass across the affected files is more reviewable than 150 individual `Edit` calls — the diff at the end is the unit of review. Reserve `Edit` for changes that need per-call judgment. Restrict `sed` to substitutions whose correctness is obvious from the pattern alone.

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
- Types name *shapes*; bindings (variables, fields, arguments) name *roles*. Prefer `type Locus = (usize, u64)` with `per_file_prev_locus: Vec<Option<Locus>>` over a role-named alias like `type PerFileOrder = Option<(usize, u64)>` — the binding carries the role, the type stays reusable across other roles. Newtype wrappers (`struct UserId(u64)`) are exempt; their job is precisely to mint type-level role distinctions.
- Avoid vague names: `item`, `data`, `value`, `thing`, `obj`, `tmp`, `helper`.
- Prefer self-explaining names over short cryptic ones, especially for private items where length costs nothing: `ReadFingerprintWithSourceFile` over `WindowEntry`, `per_file_prev_locus` over `prev_per_file`. The test is "would a reader hitting this identifier cold know what it is, without looking elsewhere?", not "is this name long?".
- At internal layer transitions, prefix raw dependency-type bindings with the library name (`noodles_cram_reader` rather than `cram_reader`) when the value is about to be wrapped by project code. Once wrapped, the project-named wrapper carries the layer. Private-code convention only — public APIs still don't expose dependency types.
- Prefer straightforward control flow with early returns over deep nesting.
- Add comments only for non-obvious invariants, algorithms, or performance tradeoffs.
- Type and field doc comments lead with what the value *is* (shape, parts, invariants), then add *why* (role, lifecycle, rationale) only when it would surprise a reader who understood the shape. Inline comments still follow the *why* not *what* rule — the difference is that doc comments on type-level declarations have no surrounding code to read the *what* off.
- In design-explaining doc comments, prefer plain English to Rust jargon ("borrowing would force every consumer to copy the bytes itself" rather than "we avoid hidden clones"). The jargon gestures; the plain version names. Mechanical code annotations can still use jargon where it's load-bearing.

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

If a command cannot be run, state that explicitly and do not fabricate output. The same rule applies to file paths, line numbers, type signatures, and behavior claims: cite only what was read or verified. When something cannot be checked, label it "needs verification" rather than asserting it as fact.

## Definition of done

A feature is complete only when all of the following are true:

- Plan was shared up front, and user feedback was incorporated when provided.
- Requirements are implemented.
- Code is readable and domain-named.
- New tests exist and pass.
- Existing tests pass (or failures are explained and unrelated).
- Formatting/linting expectations are met or explicitly documented.
- Any tradeoffs (performance, API decisions, deferred work) are clearly documented.
- The feature's block in `PROJECT_STATUS.md` reflects the new state (per *Project status protocol*).

## Response format for implementation tasks

When using this skill, structure the response in this order:

1. **Plan**: short implementation plan.
2. **Assumptions**: silent choices made when the spec was incomplete — what was unspecified, what was chosen, and why. Omit only if the spec was fully concrete.
3. **Changes made**: files and behavior added/modified.
4. **Tests added/updated**: what each new test validates.
5. **Validation results**: commands run and outcomes.
6. **Tradeoffs and follow-ups**: explicit non-goals, deferred improvements.

## Saving the implementation report

For non-trivial features, save a report mirroring the *Response format* sections to:

```
ia/reports/implementations/<feature-slug>_<YYYY-MM-DD>.md
```

Append `_v<N>` if a report for the same slug and date already exists. File references inside the report use relative Markdown links from the report's directory (e.g. `[file.rs](../../../src/file.rs#L42)`). Skip the saved report only for changes small enough that the commit message carries the full context.

## Project status protocol

The project tracks the lifecycle of every feature in `PROJECT_STATUS.md` at the project root. It is a navigation aid, not a source of truth — use it to find the relevant spec, plan, and prior reports for the in-scope feature, then verify against current code as usual.

**At task start.** Read `PROJECT_STATUS.md`. The immutable "About this project" paragraph (delimited by `ABOUT-PARAGRAPH-START` / `ABOUT-PARAGRAPH-END` HTML comments) gives the design context and points at the authoritative spec; "Current focus" confirms the project's direction and last-completed work; the per-feature blocks point at the plan, prior impl reports, and reviews for the in-scope feature.

**At task end** (after the implementation report has been saved): update only the in-scope feature's block in `PROJECT_STATUS.md`.

- If the block already exists: append a link to the new implementation report under `Impl report:` (or `Impl reports:` if the feature has several), refresh `Status:`, close `Open:` items the run resolved, add any new `Open:` items the run surfaced.
- If the block does not exist yet: create one in the matching pipeline-stage section, using the format of existing blocks (do not improvise a new format).
- Refresh the **Current focus** paragraph: rewrite `Last completed task` to name this implementation and link the new report. Touch `Next task` only if the human PM has not already set one — otherwise leave it alone, optionally appending `(suggested follow-up: …)` after the existing text.

**Status vocabulary** (one per block; use the closest match):

| Status            | Meaning                                                   |
|---|---|
| `planned`         | Plan written, no code yet                                 |
| `in-flight`       | Code being written; no implementation report yet          |
| `implemented`     | Implementation report exists; no review yet               |
| `reviewed`        | Review exists; fixes not yet applied                      |
| `fixes-applied`   | Fix-application report exists for the latest review       |
| `shipped`         | Reviewed, fixes applied, no open follow-ups               |
| `superseded`      | Replaced by another feature (link the replacement)        |

After a feature_implementation run, the typical new status is `implemented` (no review yet).

**Hard rules.**

- Do not edit the **About this project** paragraph or anything between the `ABOUT-PARAGRAPH-START` / `ABOUT-PARAGRAPH-END` comments.
- Do not modify another feature's block.
- Do not summarize report content inside the block — the block is a list of pointers.
- Prefer in-place updates of an existing bullet over accumulating a long history (`git log` and `ls reviews/` carry chronology).
- If `PROJECT_STATUS.md` and the current code disagree, trust the code; the status file is stale and should be updated, not relied on.

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