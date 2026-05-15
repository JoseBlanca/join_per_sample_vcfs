# Finding format and severity rubric (shared)

This file is read by every category sub-agent dispatched by the rust-code-review orchestrator. It defines how findings are written so the orchestrator can synthesize them.

## Severity rubric

Each finding is filed at one of four severity levels.

- **Blocker** — *Must be fixed before merge.* Wrong results, data loss or corruption, unsoundness or undefined-behavior risk, security vulnerabilities (path traversal, deserialization of untrusted input, timing-leaky secret comparison, injection), swallowed errors that hide failures, or absence of tests for any code path that, if broken, would produce wrong results without panicking.
- **Major** — *Should be fixed before merge; deferral requires written justification.* High regression risk, missing or poor error context, ambiguous or undocumented defaults, concurrency hazards (lock held across `.await`, incorrect `Send`/`Sync` bounds, unjustified memory ordering), or API design that invites misuse.
- **Minor** — *Should be fixed soon; need not block merge.* Readability or maintainability issues with moderate long-term cost: unclear naming, oversized functions, duplication, primitive obsession, missing non-critical docs.
- **Nit** — *Optional.* Small style or consistency issues that do not materially affect behavior.

### Severity interacts with confidence

A finding may be filed at **Blocker** only at **High** confidence. Lower-confidence findings about Blocker-class issues are filed at **Major** and paired with the specific verification step (test, call-site audit, author question) that would confirm or refute them.

### Volume guidance for Nits

Nits are grouped under a single `## Nits` heading, not enumerated. If there are more than ~5 nits of the same kind (formatting, naming convention, import ordering), do not list them individually — recommend a single mechanical fix instead (`cargo fmt`, a specific clippy lint, or a rename pass).

## Output format per finding

Each finding uses this exact format:

```
- `path/to/file.rs:LINE` — **[Severity]** Title
- **Confidence:** High / Medium / Low
- **Assumptions (if any):** ...
- **Problem:** 2–4 sentences with code snippet if needed.
- **Why it matters:** 1–2 sentences naming the impact.
- **Suggested fix:** unified diff, complete replacement snippet, or numbered refactor steps. Self-contained.
  ```rust
  // diff or replacement
  ```
```

Group findings by severity within your output file: **Blocker** first, then **Major**, then **Minor**, then a single `## Nits` section. Do **not** assign severity codes (B1, M1, …) — the orchestrator does that during synthesis.

If no findings apply, write a single line and stop:

```
No findings.
```

## Cross-category observations

If you notice an issue that clearly belongs in a different category, do not file it as a finding in your own output. Add a `## Cross-category observations` heading at the bottom of your file with a one-line note (`<file:line> — looks like an X issue, see Y category`). The orchestrator routes these.

## Anti-hallucination contract

- Cite only file paths and line numbers you have read. Use the Read tool first.
- Quote tool output verbatim. If you ran `cargo test`, paste the actual output. If you did not run a command, do not claim a result.
- If you cannot verify something the rule asks about, downgrade severity, set `**Confidence:** Low`, and attach an explicit verification step under `**Suggested fix:**`.
- Stay within the category. The rule of thumb: if you found yourself reading a file outside the in-scope list to make a finding, you are probably out of scope — note it and stop.
