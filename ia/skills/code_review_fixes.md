---
name: apply-code-review-fixes
description: Implement approved fixes from a completed Rust code review report conservatively and verifiably. Input: a review report plus access to the source tree. Output: narrowly-scoped source changes, targeted regression tests where needed, and a mandatory fix-application report that accounts for every finding.
---

# Apply Code Review Fixes

Implement fixes from a completed Rust code review report. This skill is downstream of the code-review workflow — its job is to take reviewed findings, decide which are safe and implementation-ready, apply only changes that can be justified and verified, and leave a precise written record.

A review report is evidence and guidance, not an executable patch. Suggested code may be illustrative, incomplete, or wrong in detail. Treat it as input to implementation, not as authoritative text to splice in.

## Core principles

- **Correctness over speed.** Prefer smaller, more verifiable changes over broad or speculative ones.
- **A review finding is not a patch.** Suggested code in a review may be illustrative, incomplete, outdated, over-broad, or wrong in detail.
- **Never silently resolve ambiguity.** If a finding depends on an assumption or policy choice, ask the user or defer. Do not guess.
- **No silent omissions.** Every finding must be accounted for explicitly in the report.
- **Minimal diff discipline.** Change only what the finding requires. No unrelated cleanup, renames, or refactors.
- **Validation must be real.** Never invent successful runs or tool results. Only claim what was actually executed.
- **Tests are part of the fix.** If a finding reveals missing coverage or an unprotected invariant, the fix is incomplete without the test.
- **Public behavior is not changed casually.** Fixes that alter public API or externally visible behavior must be recognized explicitly — use `Ask` or `Defer`, not silent inclusion.

## Preflight

Before making any code change, perform a preflight pass:

1. Read the review header, scope, open questions, and top priorities.
2. Read all findings in full — not just summaries.
3. For each finding, extract: ID, title, severity, confidence, assumptions, affected files, problem summary, suggested fix summary, whether it affects public API, whether it's testable, whether the current code still exhibits the issue.
4. Verify findings against current code using semantic context (line numbers may have drifted). A finding is **context-matched** only if the relevant code still exists and the reported behavior is still present.
5. Write the initial findings table to the report — every finding recorded with metadata and an initial decision — before editing any code.

**Stop early** if: the review can't be interpreted reliably, source files are unavailable, or the codebase has drifted so far that most findings no longer map safely. Still produce the report explaining why.

## Decision model

Every finding must be classified before implementation. Classification depends on clarity and implementability, not severity alone.

### Apply

All of the following must hold:

- severity is normally **Blocker** or **Major**; confidence is **High**
- current code still exhibits the issue
- assumptions are `None` or already resolved
- one clearly correct implementation path exists
- no new policy invention required
- no unapproved public API or semver-significant changes
- the fix can be validated directly

### Ask

Use when: the defect is credible and important, the code supports the finding, but implementation depends on a concrete user decision (e.g., reject vs. normalize invalid input, keep vs. change current API). Ask only focused questions that would produce materially different implementations.

**Question discipline:** at most 2 questions per finding; batch related questions; if more than 2–3 user decisions are needed before any work can start, report that the fix set is not implementation-ready. Bad: "How should I fix this?" Good: "B2: Should the parser reject malformed genotype strings with a typed error, or keep permissive parsing with a warning?"

`Ask` is temporary — after the user answers, move the finding to `Apply`, `Defer`, or `Dispute`. Record answers in the report.

### Defer

Use when: the finding is real but not implementation-ready — ambiguity too broad for a short question, wider redesign needed, public API impact without explicit approval, or multiple reasonable implementations requiring design decisions.

### Dispute

Use when: verification shows the finding is wrong, obsolete, or based on a misunderstanding. Requires evidence: what was checked, what was observed, why it contradicts the review.

### Already fixed

Use when: the current code already satisfies the finding. Requires verification — don't mark it just because the code looks different.

### Terminal statuses

After implementation or resolution, assign one of:

| Status | Meaning |
|---|---|
| **Applied** | Implemented as intended, validated |
| **Applied with adaptation** | Implemented successfully but differed materially from suggestion — must explain why |
| **Already fixed** | Defect no longer present |
| **Deferred** | Not implemented, remains open |
| **Disputed** | Verification contradicts the finding |
| **Failed validation** | Attempted, failed, reverted |
| **Blocked by context mismatch** | Current code can't be safely mapped to the finding |
| **Superseded** | Another finding fully absorbed this one — must reference which |

## Implementation workflow

Process findings in order: Blocker `Apply` → Major `Apply` → lower-severity only if required by a higher fix → `Ask` findings only after user answers. Never implement `Deferred` or `Disputed` findings.

### One finding at a time

Default: one finding, one scoped patch, one validation pass, one report update. Bundle only when findings are tightly coupled and can't be validated separately — the report must explain the bundling.

### Verify before editing

Immediately before editing, confirm the code still exhibits the issue and the planned fix still fits. If not, update the status (`Blocked by context mismatch`, `Already fixed`, or `Disputed`) instead of editing speculatively.

### Test-first for correctness bugs

For findings about correctness, malformed input, invariants, boundaries, or regression risk:

1. Add or strengthen a targeted test
2. Confirm it fails or demonstrates the defect
3. Change production code
4. Rerun the targeted test
5. Rerun broader validation

If test-first isn't practical, the report must explain why.

### Test quality requirements

Tests must: verify the actual bug (not a loose proxy), assert expected behavior precisely (including exact error variants when the contract is typed and stable), use inputs that exercise the reported defect, and remain narrow enough that failure points to the regression. `is_err()` alone is usually insufficient when a specific error variant is expected.

Don't copy review-proposed tests blindly — adapt to the crate's actual API, verify they demonstrate the defect, and tighten weak assertions.

### Implement narrowly

Change the smallest surface that resolves the finding. Do not replace entire functions just because the review included a full replacement example. Preserve surrounding behavior and crate conventions. If you adapt the review's suggestion, use status `Applied with adaptation` and explain the difference.

### Revert failed attempts

If validation fails: revert completely, mark `Failed validation`, record the failing command and output. Don't leave partial edits.

### Update the report incrementally

After each finding is handled, update the report immediately — don't reconstruct from memory at the end.

## Validation

Validate in layers after each fix:

1. Targeted test for the specific finding
2. Affected module/integration tests
3. `cargo fmt --check`
4. `cargo clippy --all-targets --all-features -- -D warnings`
5. `cargo test --all-targets --all-features`
6. `cargo doc --no-deps`
7. `cargo audit`

**Rules:** Use real command execution only. Do not weaken validation to get a finding over the line (e.g., skipping clippy, counting compilation as correctness). If pre-existing failures prevent clean full-project validation, record that fact and distinguish pre-existing failures from newly introduced ones.

A finding cannot receive status `Applied` or `Applied with adaptation` unless the report records: commands run, exit codes, and concise results.

## Fix-application report

**Every run must produce a report** in `reviews/`, even if no code changed, all findings were deferred, or the run was blocked. A run without a report is incomplete.

**File naming:** `reviews/fixes_applied_<YYYY-MM-DD>.md`, with `_v1`, `_v2` suffixes for multiple runs on the same date.

**Accounting rule:** every finding from the source review must appear exactly once in the findings table with a terminal status or an active `Ask` state. No finding may disappear.

**Completion rule:** the run is not complete until: every finding is in the report, every applied finding records files changed and validation, every non-applied finding records a reason, and unresolved high-severity findings are listed explicitly.

At the end of the run, tell the user where the report was written and list any unresolved high-priority findings.

### Report template

```markdown
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

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | <title> | Apply | Applied | No | `src/...` | Pass | No |
| B2 | Blocker | <title> | Ask | Deferred | Yes | None | N/A | Yes |

## 3. Questions asked and answers

1. **B2** — <question text>
   - **Answer:** <user answer / awaiting answer>

If none: "None."

## 4. Per-finding log

### <Finding ID> — <title>
- **Severity:** <Blocker / Major / Minor / Nit>
- **Initial decision:** <Apply / Ask / Defer / Dispute / Already fixed>
- **Final status:** <terminal status or Ask>
- **Reasoning:** <why this decision and status>
- **Implementation summary:** <what changed, or "None">
- **Review suggestion used verbatim?:** <Yes / No / N/A>
- **Adaptation:** <explanation if adapted; otherwise "None">
- **Verification performed:** <what was checked>
- **Files changed:** <file list or "None">
- **Tests added or modified:** <test names or "None">
- **Validation:**
  - `<command>` → <exit code>, <result>
- **User input:** <question and answer, or "None">
- **Follow-up:** <what remains, or "None">
- **Residual risk:** <if any; otherwise "None">

*(Repeat for every finding)*

## 5. Deferred findings to carry forward
- <finding ID> — <short reason>

If none: "None."

## 6. Disputed findings to return to reviewer
- <finding ID> — <short reason>

If none: "None."

## 7. Failed-validation findings
- <finding ID> — <failing command and short reason>

If none: "None."

## 8. Blocked-by-context-mismatch findings
- <finding ID> — <short reason>

If none: "None."

## 9. Commands run
- `<command>`
- `<command>`

## 10. Command results
- `<command>` → <exit code>, <one-line result>

## 11. Notes
- <important assumption, adaptation, or follow-up>

If none: "None."
```

### Template rules

- Every finding appears exactly once in the findings table and per-finding log.
- Use only the explicit status names defined in this skill. No vague labels like `skipped` or `partially done`.
- If no code changes were made, the report must explain why.
- If blocked awaiting user input, keep affected findings in `Ask` state with the exact question recorded.
