# Finding format and severity rubric (shared)

This file is read by every category sub-agent dispatched by the rust-performance-review orchestrator. It defines how findings are written so the orchestrator can synthesize them.

## Severity rubric

Performance findings are not correctness defects; severity tracks **how confidently we expect the change to produce a measurable, worth-the-complexity gain**.

- **Hot-path** — *Strong evidence the site is on a hot path (named in a profile, top of a flamegraph, or in a function backed by a reproducible benchmark) and the proposed change has a clear mechanism for a measurable win.* The candidate must account for a meaningful share of total time / allocations / lock-wait / syscalls in the quoted measurement.
- **Likely** — *Pattern-matched improvement with plausible call frequency but no profile yet.* Worth running the proposed measurement before merging. Most performance findings start here.
- **Speculative** — *The pattern exists but it is unclear the site is hot or the gain is real.* Filed so the team is aware; do not act without an experiment that contradicts the null hypothesis "this won't matter".
- **Note** — *Educational pointer.* Not actionable as written; reserved for patterns the team should recognize for future code, or for cold-path occurrences flagged for completeness.

### Severity interacts with confidence

A finding may be filed at **Hot-path** only at **High** confidence and only with quoted profile or benchmark output that names the site. Pattern-only findings cap at **Likely**. Findings about cold or once-per-run code (CLI parsing, startup, error reporting) cap at **Note** unless the orchestrator explicitly requested a sweep of cold code.

### Volume guidance

If the same pattern recurs more than ~5 times in similar form (every `format!` in a logging module, every `Vec::new` in a builder), do not enumerate. Report once with a recommendation for a mechanical sweep.

## Output format per finding

Each finding uses this exact format:

```
- `path/to/file.rs:LINE` — **[Severity]** Title
- **Confidence:** High / Medium / Low
- **Hot-path evidence:** quoted profile / benchmark / "pattern-match only"
- **Pattern matched:** which rule in this category's checklist the candidate matches
- **Mechanism:** 1–3 sentences naming the specific cost (allocation, cache miss, lock wait, syscall, branch mispredict, vtable, etc.) and why the proposed fix removes it.
- **Measurement plan:** the benchmark or profile that would confirm or refute the gain. Concrete: command, metric, and the threshold that makes the change worth merging.
- **Complexity cost:** what the fix adds (new type, lifetime, dependency, `unsafe`, build flag, extra invariant). Be honest.
- **Suggested experiment / fix:** unified diff, replacement snippet, or numbered steps. Self-contained.
  ```rust
  // diff or replacement
  ```
```

Group findings by severity within your output file: **Hot-path** first, then **Likely**, then **Speculative**, then a single **Note** section. Do **not** assign severity codes (H1, L1, …) — the orchestrator does that during synthesis.

If no findings apply, write a single line and stop:

```
No findings.
```

## Cross-category observations

If you notice an issue that clearly belongs in a different category, do not file it in your own output. Add a `## Cross-category observations` heading at the bottom of your file with a one-line note (`<file:line> — looks like an X issue, see Y category`). The orchestrator routes these.

Common cross-category routings:

- A `Mutex` lock taken inside a tight loop is both an allocation/clone smell and a contention smell — file under `concurrency`, mention from `allocations` if the call shape suggests cloning.
- False sharing on atomics is `data_layout`, not `concurrency` — but the trigger is concurrent state, so `concurrency` may notice it first.
- A `format!` allocation in a hot serializer is `allocations`; the choice between `format!` and `write!` is `hot_loops`.

## Anti-hallucination contract

- Cite only file paths and line numbers you have read. Use the Read tool first.
- Quote tool output verbatim. If you ran a benchmark or read a profile, paste the actual output. If you did not, say "pattern-match only" — never fabricate measurements.
- Do not estimate speedups in percentage terms or as multipliers unless you have run the experiment yourself. Use qualitative language ("removes a per-iteration allocation", "moves the lock out of the inner loop") instead.
- Cold code, error-handling paths, and once-per-run setup are not hot paths. Findings filed against them are downgraded to **Note** unless the orchestrator marks them in scope.
- Stay within the category. The rule of thumb: if you found yourself reading a file outside the in-scope list to make a finding, you are probably out of scope — note it and stop.
