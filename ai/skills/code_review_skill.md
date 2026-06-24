---
name: rust-code-review
description: Use this skill whenever the user asks for a code review, audit, critique, or quality check of Rust code (.rs files, Cargo crates, snippets, PRs, or diffs). Trigger on phrases like "review my Rust code", "audit this crate", "is this idiomatic Rust", or when Rust source is shared with any request implying quality feedback.
---

# Rust Code Review

You are performing a professional, uncompromising code review of Rust code. Quality is the highest priority — be precise, specific, and direct. Vague praise is forbidden; every comment must point to a concrete location and propose a concrete change.

The review is split across focused per-category checklists in `ia/skills/code_review/`. **You are the orchestrator**: you triage which categories apply to the scope, dispatch one sub-agent per category in parallel, then synthesize their findings into a single report. Per-category rules are not duplicated in this file — read each category file when dispatching the corresponding sub-agent.

## Review principles (must always hold)

- **Correctness over style.** Prioritize behavioral correctness and failure transparency over formatting or personal preferences.
- **No silent assumptions.** When the code leaves something unspecified — invariants on inputs, call-site guarantees, threading context, whether a collection is sorted, whether a value can be zero or empty — do not silently pick an answer and review as if it were fact. State the assumption explicitly in the finding and lower its severity until it can be verified.
- **Evidence-first findings.** Do not invent file paths, line numbers, logs, or command results. If you cannot verify, label as "Needs verification".
- **Actionability.** Every non-trivial finding must include a concrete fix (diff/snippet, test, or refactor step), or — if the right fix depends on intent the reviewer cannot infer — a specific question whose answer would determine the fix.
- **Scope discipline.** Review what was asked. For a diff or PR, focus on changed lines and their direct callers/callees; flag pre-existing issues in untouched code separately under "Out of scope observations". Exception: pre-existing Blocker-severity issues (security, data loss, undefined behavior) are raised under Findings regardless of scope.

The severity rubric and per-finding format are defined in `ia/skills/code_review/_finding_format.md`. Read it once at the start of every review — both you (for synthesis) and every sub-agent you dispatch will follow it.

## Review procedure

You run each numbered step once per review. Only step 6 fans out into parallel sub-agents.

### 0. Read `PROJECT_STATUS.md`

Read `PROJECT_STATUS.md` at the project root to orient on what the project is, the current focus, and the in-scope feature's prior artefacts (plan, implementation report, prior reviews, prior fix-applied reports). See *Project status protocol* below for what to read and what to update at the end.

### 1. Establish scope

Determine whether the review covers a full crate, a diff/PR, or a snippet. State the scope. For diffs, identify changed files and the direct callers/callees of changed items; these define the in-scope surface.

### 2. Inventory

List files, public API surface, error types, concurrency primitives, and external dependencies. Detect category triggers: does this code use `unsafe` / `async` / `Arc` / `Mutex` / atomics / channels? Does it have a public API? Is there a `Cargo.toml`? Is it a parser / validator / security boundary? Is it on a hot path?

### 3. Run verification commands

If a real execution environment is available, run, in the project's container per `CLAUDE.md`:

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo doc --no-deps`
- `cargo audit`

Quote actual output verbatim. **Never simulate, paraphrase, or guess at command output.** If commands cannot be run, list them under "Commands the author must run" and proceed with static review only. Any real failure is at least **Major**; correctness-impacting failures are **Blocker**. The verbatim output is passed into each sub-agent's prompt at step 6 so they do not re-run.

### 4. Determine intent

Before judging whether the code works, establish what it is meant to do — its purpose, its contract with callers, the inputs it is built to handle, the invariants it maintains. Correctness is meaningless without intent. If intent is unclear from names, types, tests, and docs, file it as a finding against documentation or naming at Minor severity (Major if the ambiguity could plausibly cause misuse).

Write a one-paragraph "domain intent" summary; it is passed into each sub-agent's prompt.

### 5. Triage categories

Decide which per-category checklists apply. Each lives at `ia/skills/code_review/<category>.md`.

| Category | Apply when |
|---|---|
| `reliability` | Always. Snippets without test files: still flag missing tests. |
| `errors` | Always. |
| `naming` | Always. |
| `defaults` | Scope contains public API, configuration, or any default-acting value. Skip pure-internal snippets with no parameters. |
| `idiomatic` | Always. |
| `refactor_safety` | Always. |
| `module_structure` | Scope spans multiple files (crate / PR / module-level review). Skip pure single-file snippets. |
| `unsafe_concurrency` | Code uses `unsafe`, `Arc`, `Mutex`/`RwLock`, atomics, channels, `async`, or thread spawning. Skip otherwise. |
| `smells` | Always. |
| `tooling` | Scope is a crate (has `Cargo.toml`). Skip pure snippets. |
| `extras` | Scope contains parsers/validators/security boundaries; accepts untrusted input; produces stable output; is on a hot path; is a public crate; or is a PR (for "Diff matches stated intent"). Apply per item. |

When in doubt, dispatch — a sub-agent that finds nothing applicable writes `No findings.` and is cheap.

### 6. Dispatch sub-agents in parallel

Create the scratch directory: `tmp/review_<YYYY-MM-DD>_<scope-slug>/` (append `_v<N>` if it already exists). The slug is a short kebab-case identifier of the reviewed module or PR (e.g. `gvcf_parser`, `pr-142`).

Project rule: scratch space is project-local `tmp/`, never `/tmp`. Add `tmp/` to `.gitignore` if it is not already covered by the existing target ignores.

For each selected category, dispatch a `general-purpose` sub-agent **in parallel** — issue a single message with multiple Agent tool calls. Each sub-agent prompt:

> Run the **<category>** checklist on the following Rust code review scope.
>
> **Scope:** <full crate / PR diff / snippet>
> **Domain intent:** <one paragraph from step 4>
> **In-scope files (full paths):** <list>
> **Out of scope:** <list, with reasons>
> **Verification command output:** <verbatim quotes from step 3, or "not run, because …">
>
> **Instructions:**
> 1. Read `ia/skills/code_review/<category>.md` for the rules to apply.
> 2. Read `ia/skills/code_review/_finding_format.md` for the severity rubric and finding format.
> 3. Read each in-scope file.
> 4. Apply each rule and produce findings in the specified format.
> 5. Write findings to `tmp/review_<date>_<slug>/<category>.md`. If no findings apply, write only the line `No findings.`
> 6. Do not invent file paths, line numbers, command output, or behavior. Cite only locations you have read.
> 7. Stay within the category. Issues that belong elsewhere go under a `## Cross-category observations` heading at the bottom of your file.

Substitute `<category>` and the scope fields for each dispatch. Do **not** assign severity codes (B1, M1, Mi1, …) inside sub-agents — that happens at synthesis.

### 7. Collect findings and route cross-category observations

For each `tmp/review_<date>_<slug>/<category>.md`, in turn:

1. **Read the file.** If it contains only `No findings.`, mark the category as clean and continue.
2. **Tally findings** by severity. Preserve each sub-agent's per-finding text — you will paste it during synthesis (step 9), with a severity code prepended.
3. **Read the `## Cross-category observations` section at the bottom.** Sub-agents are instructed to note issues that belong to other categories there rather than filing them out-of-category. For each note, decide:
   - **Promote to a finding** if the issue clearly matches a rule in the destination category and was not already raised there. Add it to that category's tally; during synthesis, file it under the destination category's severity and cite that multiple categories surfaced it ("convergent finding").
   - **Merge** if the same issue is raised in multiple places (its own category's findings *and* one or more cross-category notes from other agents). File once during synthesis, citing every category that surfaced it. Convergent evidence raises confidence; do not duplicate the entry.
   - **Defer** if the note is genuinely out of scope (different module, different concern) — record it for "Out of scope observations" in the synthesized report.
4. **Sanity-check sub-agent output.** If a sub-agent appears to have skipped its scope, cited locations it did not read (a hallucination risk), or produced findings without concrete fixes, redispatch *that one category* with an explicit instruction to fix the gap. Do not redo the others.

The point of this step is to do the routing once, deterministically, before synthesis — not to leave it as an exercise during the writeup.

### 8. Verify the test challenge

The `reliability` sub-agent runs the "challenge tests" pass for every non-trivial function. Spot-check that it did so for the changed-code surface; if a non-trivial function was missed, supplement with the missing entries before synthesis.

### 9. Synthesize the unified report

Compose the report using the *Output format* below. Verdict, top 3, and "What's good" need the full picture and are produced by you. Assign severity codes during synthesis: `B1, B2, …` for Blocker; `M1, M2, …` for Major; `Mi1, Mi2, …` for Minor; Nits stay grouped without numbering.

Each finding is filed once. Issues that were merged in step 7 (raised by multiple categories) carry a `**Categories:** <a>, <b>` line citing every sub-agent that surfaced them — this is convergent evidence and should be visible to the author, not hidden behind deduplication.

Save to `doc/devel/reports/reviews/<module-slug>_<YYYY-MM-DD>.md` per the saving conventions below. Leave the per-category files in `tmp/` as an audit trail.

### 10. Update `PROJECT_STATUS.md`

After the review report is saved, update the in-scope feature's block in `PROJECT_STATUS.md` to point at it. See *Project status protocol* below for the rules.

## Output format

Produce the synthesized report in the following order. Use the section headings verbatim so the format is machine-readable.

### 1. Scope
- What was reviewed: full crate / PR diff / snippet.
- Reviewed against: commit hash, branch, or "as-provided".
- In-scope files (list).
- Deliberately out of scope (list, with reason).
- Categories dispatched (list, each with a one-line reason).

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
Grouped by severity (**Blocker** → **Major** → **Minor** → **Nits**). Within each severity, ordered by confidence (High → Low), then by file. Nits collected into a single sub-section, not enumerated. Each finding follows the format defined in `ia/skills/code_review/_finding_format.md`, with the severity code (B1, M1, …) prepended to the title — e.g. `B1: src/parser.rs:42 — Title`.

### 7. Out of scope observations
Pre-existing issues in untouched code, surfaced but not blocking. Each: file, brief description, suggested follow-up (separate PR or issue). Pre-existing Blocker-severity issues (security, data loss, UB) appear under Findings instead, marked "pre-existing".

### 8. Missing tests to add now
Each: proposed test name in `function_returns_expected_on_condition` form, input class covered, specific bug it would catch, and the test as code or specification. Grouped by function under test. The `reliability` sub-agent's challenge-tests output feeds this section directly.

### 9. What's good
Up to 5 specific, transferable patterns worth keeping, each one sentence with a file reference. No general praise. Skip the section entirely if nothing specific qualifies.

### 10. Commands to re-verify
- Commands the reviewer ran (re-run to confirm they still pass).
- New commands or test invocations the review introduced.

### Author response convention
Address each finding by its identifier (e.g., "B2", "M5") with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer open questions from section 4 first.

---

Be direct. If something is wrong, say so plainly and show the fix. Vague praise and vague criticism are equally useless.

## Saving the report

### Directory and filename

Save to the project's `doc/devel/reports/reviews/` directory:

```
doc/devel/reports/reviews/<module-slug>_<YYYY-MM-DD>.md
```

Examples:

- `doc/devel/reports/reviews/gvcf_parser_2026-04-13.md`
- `doc/devel/reports/reviews/genotype_merging_2026-04-13.md`
- `doc/devel/reports/reviews/pr-142_2026-04-13.md`

If a review for the same scope and date already exists, append `_v<N>`:

- `doc/devel/reports/reviews/gvcf_parser_2026-04-13_v2.md`

### Document header

```markdown
# Code Review: <module-slug>
**Date:** <YYYY-MM-DD>
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** <one-line description>
**Status:** <Approve / Approve-with-changes / Request-changes>

---
```

The body is sections 1–10 of *Output format* above, in order, with verbatim headings.

### File links inside findings

References to source files use relative Markdown links from the `doc/devel/reports/reviews/` directory:

- Single line: `[file.rs](../../../../path/file.rs#L123)`
- Range: `[file.rs](../../../../path/file.rs#L123-L456)`

Display text is the path (no backticks).

### Pre-save checklist

- [ ] Every Blocker finding has High confidence; lower-confidence Blocker-class issues are filed at Major with a verification step.
- [ ] Every Major and Minor finding has a concrete fix (code, test, or specific question).
- [ ] Every cited file:line was actually read (no invented locations).
- [ ] Test recommendations use `function_returns_expected_on_condition` naming.
- [ ] "What's good" has 3–5 specific patterns or is omitted.
- [ ] File paths are repo-relative (`src/foo.rs`, not `/home/...`).
- [ ] Severity codes (B1, M1, Mi1) are consistent and dense (no gaps).
- [ ] Open questions are numbered and referenced from the findings they affect.
- [ ] Per-category files in `tmp/review_<date>_<slug>/` are left in place as an audit trail.
- [ ] `PROJECT_STATUS.md` updated (per *Project status protocol*).

## Project status protocol

The project tracks the lifecycle of every feature in `PROJECT_STATUS.md` at the project root. It is a navigation aid, not a source of truth — use it to find the relevant spec, plan, and prior reports for the in-scope feature, then verify against current code as usual.

**At task start.** Read `PROJECT_STATUS.md`. The immutable "About this project" paragraph (delimited by `ABOUT-PARAGRAPH-START` / `ABOUT-PARAGRAPH-END` HTML comments) gives the design context and points at the authoritative spec; "Current focus" confirms the project's direction and last-completed work; the per-feature blocks point at the plan, prior impl reports, and reviews for the in-scope feature.

**At task end** (after the review report has been saved): update only the in-scope feature's block in `PROJECT_STATUS.md`.

- Append a link to the new review under `Latest review:` (or `Latest reviews:` if multiple). Prefer to replace the previous link rather than accumulate a long list — `git log` and `ls doc/devel/reports/reviews/` carry chronology.
- Update `Status:` to `reviewed` (or keep `fixes-applied` / `shipped` if the review found nothing material; explain in one trailing word: `shipped (re-reviewed 2026-MM-DD)`).
- Add any new `Open:` items the review surfaced; do not close existing ones (a review surfaces, it does not resolve — that is the apply-fixes skill's job).
- Refresh **Current focus** — rewrite `Last completed task` to name this review and link the report. Touch `Next task` only if the human PM has not already set one; otherwise leave it alone, optionally appending `(suggested follow-up: apply fixes)`.

**Status vocabulary:** `planned` / `in-flight` / `implemented` / `reviewed` / `fixes-applied` / `shipped` / `superseded`. After a code-review run, the typical new status is `reviewed`.

**Hard rules.**

- Do not edit the **About this project** paragraph or anything between the `ABOUT-PARAGRAPH-START` / `ABOUT-PARAGRAPH-END` comments.
- Do not modify another feature's block.
- Do not summarize the review's findings inside the block — the block is a list of pointers; findings live in the saved report.
- If the in-scope feature has no block yet, create one using the format of existing blocks.
- If `PROJECT_STATUS.md` and the current code disagree, trust the code; the status file is stale and should be updated, not relied on.

## Reusable prompt template

Use this to invoke the skill consistently. Fill in the Context block; the rest defers to the skill body.

> Perform a Rust code review per the **rust-code-review** skill. Follow its principles, procedure, severity rubric, and output format in full — do not abbreviate or skip sections.
>
> **Context**
> - **Scope:** <files / module / PR diff / branch comparison>
> - **Domain intent:** <one-paragraph description of what this code is meant to do>
> - **Audience:** <internal service / public library / CLI / embedded>
> - **Constraints:** <performance budgets, MSRV, `no_std`, target platforms, deadline pressure>
> - **Out of scope:** <legacy modules being deleted, generated code, vendored deps>
> - **Prior review history:** <previously reviewed? known tracked issues?>
>
> **Anti-hallucination contract.** Quote tool output verbatim. If a command was not run, list it under "Commands not run". If a file or line cannot be located, say so. Never invent file paths, line numbers, error messages, clippy warnings, test results, or behavior. Findings without verifiable evidence are labeled "Needs verification" per the skill.
>
> **Reminders of the most-violated rules** (not a substitute for the per-category checklists):
> 1. Reliability first — verify behavior, not style.
> 2. Errors must never pass silently.
> 3. Tests cover edge cases and every error path.
> 4. Names are precise, domain-relevant, verb-based for functions.
> 5. Defaults are visible at the call site, in docs, and at runtime.
> 6. Code smells get concrete refactors, not vague complaints.
> 7. Make the compiler flag refactors — no `..Default::default()`, exhaustive destructures.
