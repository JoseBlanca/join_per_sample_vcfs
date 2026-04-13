---
name: apply-code-review-fixes
description: Implement approved fixes from a completed Rust code review report conservatively and verifiably. Input: a review report plus access to the source tree. Output: narrowly-scoped source changes, targeted regression tests where needed, and a mandatory fix-application report that accounts for every finding.
---

# Apply Code Review Fixes

Use this skill to implement fixes from a completed Rust code review report.

This skill is downstream of the code-review workflow. Its job is not to perform a new review, and it is not a mechanical patch applicator. Its job is to take reviewed findings, decide which ones are safe and implementation-ready, apply only those changes that can be justified and verified, and leave a precise written record of what was done and what was not done.

A review report is evidence and guidance, not an executable patch.

## What this skill does

- reads a completed review report
- extracts and tracks every finding in that report
- decides, for each finding, whether to apply it, ask the user for clarification, defer it, dispute it, or mark it already fixed
- implements only those fixes that are clearly correct and sufficiently verified
- adds or strengthens targeted tests when needed
- validates the resulting code with the appropriate Rust tooling
- writes a mandatory fix-application report that accounts for every finding

## Core principles

- **Correctness over speed.** Prefer a smaller, slower, more verifiable change over a broad or speculative one.

- **A review finding is not a patch.** Suggested code in a review may be illustrative, incomplete, outdated, over-broad, or wrong in detail. Treat it as input to implementation, not as authoritative text to splice into the source tree.

- **Never silently resolve ambiguity.** If a finding depends on an assumption, unresolved invariant, or policy choice, do not guess. Either ask the user a focused question or defer the finding.

- **No silent omissions.** Every finding in the review must be accounted for explicitly. Nothing disappears just because it was inconvenient, ambiguous, or already fixed.

- **Minimal diff discipline.** Change only what is required to resolve the finding and its direct fallout. Do not mix in unrelated cleanup, renames, refactors, or linter-driven rewrites.

- **Reproduce or verify before fixing when practical.** For correctness defects, prefer to demonstrate the problem with a targeted test or direct verification before editing production code.

- **Validation must be real.** Do not claim a fix is complete unless it has been validated with actual commands and their real outputs. Never invent successful runs, line matches, or tool results.

- **Tests are part of the fix when needed.** If a finding reveals missing coverage, weak assertions, or an unprotected invariant, the fix is incomplete without the appropriate test improvement.

- **Public behavior is not changed casually.** If a fix would alter public API, semver expectations, or externally visible behavior, that impact must be recognized explicitly rather than smuggled in as an implementation detail.

- **Interactive clarification is allowed, but controlled.** Ask only when a short, concrete user decision would unlock a correct implementation of an important finding.

- **No work without accounting.** Every run must end with a fix-application report describing what was applied, adapted, deferred, disputed, blocked, or already fixed.

## Inputs and outputs

### Required inputs

This skill expects the following inputs:

- a completed Rust code review report in the established review format
- access to the relevant source files
- access to the relevant test files
- access to the current working tree or branch state

If available, the skill should also use:

- author responses to the original review
- prior fix-application reports for the same review
- any user decisions already made about disputed or deferred findings

### Minimum input assumptions

The skill may proceed only if:

- the review report is readable
- individual findings can be identified from the report
- the referenced source files are available
- the current code can be inspected against the report

If those conditions are not met, do not improvise. Report the limitation clearly and do not claim that fixes were applied.

### Required outputs

Every run of this skill must produce the following outputs:

1. **Modified source files**, if one or more findings are applied.
2. **Modified or added tests**, when needed to reproduce, verify, or guard the fix.
3. **A fix-application report**, written in the designated `reviews/` directory. This report, and in particular its Executive Summary section, is the authoritative output to the user. At the end of the run, tell the user where the report was written and list any unresolved high-priority findings that require attention.

### Mandatory reporting rule

The fix-application report is mandatory even if:

- no code changes were made
- all findings were deferred
- the report was disputed
- the issues were already fixed
- the run stopped pending user clarification

A run without a report is incomplete.

### Output quality requirements

All outputs produced by this skill must satisfy these requirements:

- changes are narrowly scoped to the findings being implemented
- tests added or edited are relevant to the actual defect
- validation claims are backed by real command execution
- every finding from the source review is accounted for in the report
- unresolved work is explicitly carried forward rather than implied or forgotten

## Preflight and finding extraction

Before making any code change, perform a preflight pass over the review and the current source tree.

The purpose of preflight is to determine whether the review can be acted on safely, what findings exist, whether the current code still matches the report, and which findings are actually implementation-ready.

### Preflight steps

1. Read the review header and scope.
2. Read the `Open Questions & Assumptions` section.
3. Read the `Top 3 Priorities` section.
4. Read all findings in full, not just the summary lines.
5. Identify every finding ID and extract its metadata.
6. Locate the referenced files and inspect the cited code or its current equivalent.
7. Check whether the current code still matches the report closely enough to proceed.
8. Complete the initial findings table in the report before editing any code.

Do not modify source files before this preflight is complete.

### Metadata to extract for each finding

For each finding, extract and record at least:

- finding ID
- title
- severity
- confidence
- assumptions
- affected file or files
- cited location, if present
- summary of the problem
- summary of the suggested fix
- whether the finding appears to affect public API or externally visible behavior
- whether the finding appears testable directly
- whether the current code still appears to exhibit the issue

### Verifying the current code against the review

A review may be partially stale. Before acting on a finding, verify whether the current code still matches the reported problem.

Use semantic context, not raw line numbers alone. Line numbers may drift.

A finding is considered **context-matched** only if:

- the relevant code still exists
- the reported behavior or missing guard still appears present
- the review’s reasoning still fits the current implementation structure

If the code has changed enough that the finding cannot be safely mapped to the current implementation, do not guess. Mark it later as `Blocked by context mismatch` unless a careful re-verification resolves it.

### When preflight should stop early

Stop and report rather than continuing if any of the following hold:

- the review file cannot be interpreted reliably
- findings cannot be extracted clearly
- the relevant source files are unavailable
- the codebase has drifted so far that most findings no longer map safely
- the report’s scope is too ambiguous to determine what is in or out of scope

In that case, still produce the fix-application report and record why the run could not proceed.

### Preflight output

The output of preflight is the initial state of the report's findings table, with every finding from the review recorded and each assigned:

- extracted metadata
- current context status
- an initial decision state:
  - `Apply`
  - `Ask`
  - `Defer`
  - `Dispute`
  - `Already fixed`

This initial state must be written to the report before any code edits begin.

## Decision model: Apply / Ask / Defer / Dispute / Already fixed

Every finding must be classified before implementation begins.

Classification is not based on severity alone. A finding is not ready to implement just because it is important. It must also be sufficiently clear, sufficiently supported by the current code, and sufficiently testable.

### Apply

Choose `Apply` only when all of the following hold:

- the finding is important enough to act on in this run, normally **Blocker** or **Major**
- confidence is **High**
- the current code still exhibits the issue
- assumptions are either `None` or already resolved explicitly
- no open question materially affects the implementation
- there is one clearly correct implementation path
- the fix does not require inventing new policy
- the change does not introduce unapproved public API or semver-significant behavior changes
- the fix can be validated directly

If any of these conditions fail, do not mark the finding `Apply`.

### Ask

Choose `Ask` when:

- the defect is credible and important
- the current code supports the finding
- the implementation depends on one or a small number of explicit user decisions
- the likely answers would lead to materially different implementations
- a short, focused question could resolve the ambiguity

Use `Ask` for cases such as:

- reject invalid input vs normalize it
- keep current API vs change it
- document a narrow contract vs implement broader support
- choose between two concrete behavioral policies

`Ask` is a temporary working state, not a terminal status.

### Defer

Choose `Defer` when:

- the finding is real or plausible, but not implementation-ready
- the ambiguity is too broad to resolve with a short user question
- the fix would require wider redesign
- the fix would widen scope materially beyond the reviewed surface
- the change affects public API or compatibility and no explicit approval exists
- multiple reasonable implementations exist and choosing among them would amount to designing the feature

Deferred means the finding remains open and must be carried forward explicitly.

### Dispute

Choose `Dispute` when verification against the current code or behavior indicates that the review finding is wrong, obsolete, or based on a misunderstanding.

A disputed finding requires evidence. Do not dispute casually.

At minimum, a disputed finding must record:

- what was checked
- what was observed
- why that observation conflicts with the review

### Already fixed

Choose `Already fixed` when the current code already satisfies the finding.

This still requires verification. Do not mark a finding already fixed just because the code looks different. Confirm that the underlying issue is no longer present.

### General decision rules

- Severity affects priority, not certainty.
- Confidence affects whether implementation may proceed without further clarification.
- A high-severity finding may still be deferred if it is not implementation-ready.
- A lower-severity finding may be applied if it is necessary to complete a higher-severity fix correctly.
- Do not leave findings unclassified.

### Terminal statuses

Once a finding is implemented or otherwise resolved, assign one of these terminal statuses:

#### Applied
The finding was implemented materially as intended and validated successfully.

#### Applied with adaptation
The finding was implemented successfully, but the final implementation differed materially from the review's suggested fix.

Use this when:
- the review proposed a broader change than necessary
- the suggested patch did not fit the current code structure
- an existing crate pattern or type was used instead of the exact proposed mechanism
- the test strategy or normalization point was changed for good reason

This status must explain the adaptation.

#### Already fixed
The underlying defect is no longer present in the current code.

#### Deferred
The finding remains open but was not implemented in this run.

#### Disputed
Verification contradicts the review finding.

#### Failed validation
Implementation was attempted, then reverted because validation failed.

#### Blocked by context mismatch
The current code cannot be safely mapped to the report finding.

#### Superseded
Another implemented finding fully absorbed this one and no separate handling remains necessary. Must reference the absorbing finding.

## Interactive clarification policy

Interactive clarification is permitted, but only in a narrow and controlled way.

The purpose of interaction is to unblock important findings whose implementation depends on one or two concrete user decisions. It is not a substitute for analysis, and it must not turn the skill into an open-ended design interview.

### When to ask

Ask the user only when all of the following hold:

- the finding is important, normally **Blocker** or **Major**, or is required to complete one of those correctly
- the code supports the finding well enough that the ambiguity is about implementation policy, not about whether the issue exists
- a short question can resolve the ambiguity
- the answer is likely to produce a concrete implementation path

Examples of acceptable questions:

- strict rejection or permissive handling of malformed input
- fixed behavior or configurable behavior
- preserve current public API or allow a breaking change
- support only the currently documented domain or broaden support now

### When not to ask

Do **not** ask when:

- the answer should be inferable from the code, tests, docs, or review
- the ambiguity is too broad and would still leave design work after the answer
- the finding is low priority
- the question is vague or open-ended
- asking would merely offload implementation thinking to the user

Bad question:
- `How should I fix this?`

Good question:
- `B2: For malformed genotype strings, should the parser reject them with a typed error, or keep permissive parsing and surface a warning?`

### Question budget

Keep the number of questions tightly bounded.

- Prefer at most **two questions per finding**.
- Prefer batching related questions into a short numbered list.
- If safe progress would require more than **2–3 user decisions** before any meaningful implementation can start, stop and report that the fix set is not implementation-ready.

### Effect of user answers

When the user answers:

- record the answer in the fix-application report
- update the finding's status in the report
- move the finding from `Ask` to one of:
  - `Apply`
  - `Defer`
  - `Dispute`

Do not continue as if an answer had been given when it has not.

## Implementation workflow

Once preflight is complete and the initial findings table is written to the report, implement fixes in a controlled order.

The workflow is designed to prevent speculative edits, hidden scope creep, and untraceable partial work.

### Order of work

Process findings in this general order:

1. `Apply` findings that are **Blocker**
2. `Apply` findings that are **Major**
3. lower-severity findings only if they are required to complete a higher-severity fix correctly
4. findings in `Ask` state only after the user has answered
5. do not implement deferred or disputed findings in the same run unless their status is explicitly changed first

Severity determines priority, but not permission. Only findings already classified as implementation-ready may be edited.

### One finding at a time

Implement findings one at a time unless two findings are tightly coupled and cannot be separated safely.

Default rule:
- one finding
- one scoped patch
- one validation pass
- one report update

Bundling multiple findings into one edit is allowed only when:
- they affect the same narrow code path
- validating them separately would be artificial
- the report clearly explains that they are interdependent

If findings are bundled, the report must say so explicitly.

### Verify before editing

Immediately before editing a file for a given finding:

- locate the relevant code
- confirm the current code still exhibits the issue or still lacks the required guard/test
- confirm the planned fix still fits the current implementation structure

If that verification fails, do not edit speculatively. Update the finding status to `Blocked by context mismatch`, `Already fixed`, or `Disputed`, whichever is justified.

### Prefer test-first for correctness bugs

For findings about correctness, malformed input, invariants, boundaries, or regression risk, prefer this order:

1. add or strengthen a targeted test
2. confirm that it fails or otherwise demonstrates the current problem
3. change production code
4. rerun the targeted test
5. rerun broader validation

Use this especially for:
- parser acceptance/rejection defects
- overflow/underflow defects
- invalid state handling
- bugs that current tests fail to catch
- defects where the review explicitly called out missing coverage

### Implement narrowly

When changing production code:

- change the smallest surface that resolves the finding
- preserve surrounding behavior unless the finding requires changing it
- avoid opportunistic cleanup
- avoid unrelated renames
- avoid broad rewrites unless the structure itself is the problem
- preserve existing crate conventions unless they are part of the defect

Do not replace an entire function just because the review included a full replacement example. Use the smallest correct implementation.

### Adapt review suggestions deliberately

If the review’s suggested fix is not the best implementation in the current codebase, adapt it deliberately.

Acceptable reasons to adapt include:
- the suggestion is broader than necessary
- the suggestion does not fit the crate’s actual error model
- the suggestion would introduce avoidable API churn
- the current code allows a smaller, clearer, safer fix
- the current test structure calls for a different regression strategy

Whenever this happens, the final status should usually be `Applied with adaptation`, and the report must describe the difference.

### Update the report immediately

After each finding is handled:

- update its entry in the report's findings table
- record files changed
- record tests added or modified
- record validation commands run
- record the resulting status

Do not wait until the end of the run to reconstruct what happened.

### Revert failed attempts completely

If a finding is implemented and then fails validation:

- revert that finding’s code changes fully
- do not leave partial edits in place
- mark the finding `Failed validation`
- record the failing command and the relevant output
- move on only if the remaining findings can still be handled cleanly

### End-of-run reconciliation

Before closing the run:

- ensure every finding from the review appears exactly once in the report
- ensure every finding has either a terminal status or an active `Ask` state awaiting user input
- ensure unresolved high-severity findings are clearly carried forward


## Validation policy

A fix is not complete just because the code compiles or because one test passes. Validation must be strong enough to support the claim that the finding was actually resolved.

### Validation principles

- Validate the specific defect first, then the surrounding surface, then the wider project.
- Use real command execution only. Never invent or paraphrase successful results.
- Validation must be at least as strict as the original review standard, unless the project state makes that impossible.
- If full validation cannot be run, say so explicitly and record exactly what was run instead.

### Layered validation order

After implementing a finding, validate in layers:

1. the targeted test or direct check for that finding
2. any affected module or integration tests
3. formatting check
4. clippy on affected targets or the project
5. full project test, lint, doc, and security validation where available

This layered order catches local regressions early while still requiring broad confirmation before claiming completion.

### Standard validation commands

Use the following commands where the environment and project allow:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
cargo doc --no-deps
cargo audit
```

If the original review used narrower commands because the broader project state was already broken or out of scope, record that explicitly. Do not claim full-project success when only partial validation was run.

### Validation requirements by outcome

#### Applied

For a finding with final status Applied, the report must include:

- the targeted validation for that finding
- the broader validation run after the change
- the commands executed
- exit codes
- concise results

#### Applied with adaptation

For a finding with final status Applied with adaptation, the report must include:

- the targeted validation for that finding
- the broader validation run after the change
- the commands executed
- exit codes
- concise results
- a short note explaining why the implementation differed materially from the review suggestion

#### Failed validation

For a finding with final status Failed validation, the report must include:

- the command that failed
- the exit code
- the relevant output
- confirmation that the attempted code change was reverted

#### Already fixed

For a finding with final status Already fixed, the report must include:

- the verification step performed
- what was observed
- why the original finding no longer applies
- confirmation that no code change was made

#### Disputed

For a finding with final status Disputed, the report must include:

- the verification step performed
- what was observed
- why that observation contradicts the review finding
- confirmation that no code change was made

#### Deferred

For a finding with final status Deferred, the report must include:

- why the finding was not implementation-ready
- whether user clarification was attempted
- what follow-up is required before implementation can proceed


#### Blocked by context mismatch

For a finding with final status Blocked by context mismatch, the report must include:

- what code was inspected
- why the report could not be mapped safely to the current implementation
- why proceeding would have required speculation

### No weakening of standards

Do not weaken validation just to get a finding over the line.

Examples of unacceptable behavior:

- running only a targeted test when the change affects shared code
- skipping clippy because a new warning is inconvenient
- counting compilation as proof of correctness
- claiming a fix is complete when the relevant regression test was never added or run

### Validation under pre-existing project failures

If unrelated pre-existing failures prevent clean full-project validation:

- record that fact explicitly
- still run the strongest relevant validation possible
- distinguish clearly between:
  - failures introduced by the attempted fix
  - failures already present before the fix

Do not let pre-existing project noise become an excuse for weak local verification.

### Validation and the report

A finding must not be given final status Applied or Applied with adaptation unless the validation required for that finding has been recorded in the fix-application report.

At minimum, each applied finding must record:

- the validation commands run
- the exit codes
- a concise result summary
- whether the validation was targeted, broader, or full-project

If this information is missing, the finding is not fully accounted for.

## Test policy

Tests are not optional support work. When a finding exposes weak coverage, a missing regression guard, or unverified behavior, the fix is incomplete without the necessary test changes.

### When tests must be added or strengthened

Add or strengthen tests when a finding involves any of the following:

- incorrect behavior
- malformed input handling
- boundary values
- overflow or underflow risk
- invariant enforcement
- public contract expectations
- a previously untested error path
- a bug that current tests failed to catch
- a review finding explicitly identifying missing or weak tests

### Test-first preference

Prefer a test-first workflow for correctness findings whenever practical:

1. add or refine a targeted test
2. show that it fails or otherwise demonstrates the defect
3. implement the fix
4. rerun the test to show the defect is now covered

If a failing test cannot be added first, the fix-application report must explain why.

### What a good test must do

Tests added or edited by this skill must:

- verify the actual bug, not a loose proxy
- assert the expected behavior precisely
- assert the exact error variant when the contract is typed and stable enough for that assertion
- use inputs that directly exercise the reported defect
- remain narrow enough that a failure points clearly to the regression
- avoid depending accidentally on unrelated behavior

Weak assertions such as `is_err()` are often insufficient when the contract expects a specific typed error or specific failure mode.

### Naming expectations

Test names should communicate the behavior under test and the condition being exercised.

Preferred style:

- `parse_genotype_returns_error_on_invalid_character`
- `get_span_returns_error_on_position_overflow`

Do not force a naming scheme that conflicts with an already consistent and clear crate convention. Prefer clarity and local consistency over mechanical renaming.

### Using review-proposed tests

If the review includes sample test code:

- do not copy it blindly
- adapt it to the crate’s actual API and conventions
- verify that it genuinely demonstrates the defect
- tighten weak assertions before accepting it
- ensure the test still makes sense in the current codebase

A review-proposed test is input to implementation, not automatically the final version.

### Scope discipline for tests

Test edits must stay aligned with the finding being implemented.

Do not:

- add unrelated test cleanup in the same patch
- rewrite large parts of the test suite unnecessarily
- change the crate’s overall test style while fixing a narrow defect
- widen the patch into general test refactoring unless that refactoring is directly required for the finding

The goal is to protect the specific defect and its direct invariants, not to redesign the test suite.

### Regression protection

If a finding concerns a bug that could plausibly recur, add a regression test that would fail if the same defect returned.

The fix-application report should note explicitly when a test was added as a named regression guard.

### Relationship between code changes and test changes

When a finding reveals both a production defect and insufficient tests:

- the production fix and the test improvement belong to the same implementation unit
- do not treat the code fix as complete if the defect remains unguarded
- do not mark the finding `Applied` or `Applied with adaptation` unless the necessary test coverage has been added or the report explains why no new test was needed

### When no new test is needed

A finding may be fixed without adding a new test only if at least one of the following is true:

- an existing test already covers the defect precisely once corrected
- the finding concerns a non-behavioral issue that does not require regression protection
- the relevant behavior is already covered by stronger existing tests and the fix simply aligns implementation with that coverage

In such cases, the report must state why no new or modified test was required.

### Test validation requirements

For every test added or modified as part of a fix, the report must record:

- which test was added or changed
- what defect or invariant it covers
- whether it failed before the fix, or how the defect was otherwise demonstrated
- the result after the fix

If the test could not be shown failing first, the report must explain why that was not practical.

## API and public contract

If a fix would change a public function signature, alter the semantics of a public return value, add or remove public error variants, change which inputs are accepted or rejected in a way callers can observe, or otherwise affect the externally visible behavior of the crate: recognize it explicitly rather than treating it as an implementation detail.

Default rule: do not apply such changes automatically. Use `Ask` if one focused user decision unlocks a correct implementation. Use `Defer` if the change requires broader compatibility or release discussion. Do not smuggle semver-significant changes into a narrowly-scoped correctness patch.

If a public-impacting change is applied, the implementation is incomplete unless the corresponding doc comments, error documentation, or changelog entry are also updated.

When in doubt about whether a fix changes the public contract, treat it as if it does.

## Mandatory fix-application report

Every run of this skill must produce a fix-application report in the `reviews/` directory.

This report is mandatory even if:

- no code changes were made
- all findings were deferred
- all relevant issues were already fixed
- the run stopped awaiting user clarification
- the review was disputed
- the current code could not be mapped safely to the report

A run without a fix-application report is incomplete.

### Purpose of the report

The fix-application report exists to provide:

- **tracking** — every finding from the source review is accounted for
- **auditability** — later readers can reconstruct what was changed, what was not changed, and why
- **verification** — validation commands and outcomes are recorded
- **handoff** — unresolved findings are carried forward explicitly
- **control** — the report prevents silent omissions

### Accounting rule

Before editing code, complete the findings table in the report with every finding ID from the source review.

At the end of the run, every finding must appear in the report with either:

- a terminal status:
  - `Applied`
  - `Applied with adaptation`
  - `Already fixed`
  - `Deferred`
  - `Disputed`
  - `Failed validation`
  - `Blocked by context mismatch`
  - `Superseded`

or:

- an active temporary state:
  - `Ask`, only if the run is awaiting user clarification

No finding may disappear from the report.

### Completion rule

The run is not complete until all of the following are true:

- every finding from the source review appears in the report
- no finding has more than one terminal status
- every applied finding records files changed and validation performed
- every non-applied finding records a reason
- unresolved high-severity findings are listed explicitly
- the report has been written to disk

### Report file naming

Write the report to:

- `reviews/fixes_applied_<YYYY-MM-DD>.md`

If multiple runs occur on the same date, append a version suffix:

- `reviews/fixes_applied_<YYYY-MM-DD>_v1.md`
- `reviews/fixes_applied_<YYYY-MM-DD>_v2.md`

### Reconciliation check

Before finishing, verify all of the following:

- every finding ID from the review is present in the report
- no finding has more than one final status
- every `Applied` or `Applied with adaptation` finding lists changed files and validation
- every `Deferred`, `Disputed`, `Failed validation`, and `Blocked by context mismatch` finding includes a reason
- every `Ask` finding includes the question asked and either the answer received or a note that the run is awaiting answer
- every `Superseded` finding names the finding that absorbed it

If any of these checks fail, the report is incomplete.

### The report is the control artifact

The fix-application report is not a cosmetic summary written after the real work. It is one of the main control artifacts of the workflow. Fill it in incrementally as you work — do not reconstruct it from memory at the end of the run.

### Report required even for blocked runs

If the run cannot proceed because of ambiguity, missing files, context drift, validation failure, or pending user answers, still write the report.

In those cases, the report must state clearly:

- what prevented progress
- which findings remain unresolved
- what user input or follow-up work is required
- whether any code was changed before the run stopped

### No silent carry-forward

A finding may be carried forward to a later run, but only if the report states that explicitly.

For every carried-forward finding, the report must record:

- the finding ID
- its final status in this run
- why it was not completed
- what would be needed to complete it in a later run

Do not rely on memory or implicit handoff.

## Report template

Use the following structure exactly for the fix-application report.

# Fix Application Report: <original-review-report.md>

**Date:** YYYY-MM-DD  
**Source review:** `reviews/<review-file>.md`  
**Source state reviewed against:** <commit / branch / as-provided>  
**Execution mode:** <interactive / non-interactive>  
**Overall status:** <Completed / Partial / Blocked>

---

## 1. Executive summary

### Review totals
- Blockers: N
- Majors: N
- Minors: N
- Nits: N

### Outcome totals
- Applied: N
- Applied with adaptation: N
- Already fixed: N
- Deferred: N
- Disputed: N
- Failed validation: N
- Blocked by context mismatch: N
- Superseded: N
- Awaiting user answer: N

### Validation summary
- `cargo fmt --check` → <exit code / not run>, <result>
- `cargo clippy --all-targets --all-features -- -D warnings` → <exit code / not run>, <result>
- `cargo test --all-targets --all-features` → <exit code / not run>, <result>
- `cargo doc --no-deps` → <exit code / not run>, <result>
- `cargo audit` → <exit code / not run>, <result>

### Unresolved high-priority findings
- <finding ID> — <short reason>
- <finding ID> — <short reason>

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | <title> | Apply | Applied | No | `src/...`, `tests/...` | Pass | No |
| B2 | Blocker | <title> | Ask | Deferred | Yes | None | N/A | Yes |
| M1 | Major | <title> | Apply | Applied with adaptation | No | `src/...` | Pass | No |

## 3. Questions asked and answers

1. **B2** — <question text>  
   - **Answer:** <user answer / awaiting answer>

2. **M4** — <question text>  
   - **Answer:** <user answer / awaiting answer>

If no questions were asked, write:

- None.

## 4. Per-finding log

### <Finding ID> — <title>
- **Severity:** <Blocker / Major / Minor / Nit>
- **Initial decision:** <Apply / Ask / Defer / Dispute / Already fixed>
- **Final status:** <Applied / Applied with adaptation / Already fixed / Deferred / Disputed / Failed validation / Blocked by context mismatch / Superseded / Ask>
- **Reasoning:** <why this finding received its initial decision and final status>
- **Implementation summary:** <what changed, or state "None">
- **Review suggestion used verbatim?:** <Yes / No / Not applicable>
- **Adaptation:** <if applicable; otherwise "None">
- **Verification performed:** <what was checked before deciding or changing code>
- **Files changed:** <file list or "None">
- **Tests added or modified:** <test names or "None">
- **Validation:**
  - `<command>` → <exit code>, <result>
  - `<command>` → <exit code>, <result>
- **User input:** <question and answer, or "None">
- **Follow-up:** <what remains to be done, or "None">
- **Residual risk:** <if any; otherwise "None">

Repeat this subsection for every finding in the source review.

## 5. Deferred findings to carry forward
- <finding ID> — <short reason>
- <finding ID> — <short reason>

If none, write:

- None.

## 6. Disputed findings to return to reviewer
- <finding ID> — <short reason>
- <finding ID> — <short reason>

If none, write:

- None.

## 7. Failed-validation findings
- <finding ID> — <failing command and short reason>
- <finding ID> — <failing command and short reason>

If none, write:

- None.

## 8. Blocked-by-context-mismatch findings
- <finding ID> — <short reason>
- <finding ID> — <short reason>

If none, write:

- None.

## 9. Commands run

List the exact commands actually executed, in order.

Example:

- `cargo test --test gvcf_parser_test`
- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`

## 10. Command results
- `<command>` → <exit code>, <one-line result>
- `<command>` → <exit code>, <one-line result>
- `<command>` → <exit code>, <one-line result>

## 11. Notes
- <important assumption resolved>
- <important adaptation relative to the review suggestion>
- <follow-up work required>

If none, write:

- None.

### Template rules

- Every finding from the source review must appear exactly once in the findings table and exactly once in the per-finding log.
- Use explicit status names only:
  - `Applied`
  - `Applied with adaptation`
  - `Already fixed`
  - `Deferred`
  - `Disputed`
  - `Failed validation`
  - `Blocked by context mismatch`
  - `Superseded`
  - `Ask` only when awaiting a user answer
- Do not use vague labels such as `skipped`, `handled`, or `partially done`.
- If no code changes were made, the report is still required and must explain why.
- If the run is blocked awaiting user input, keep the affected finding in `Ask` state and record the exact question.

