---
name: rust-performance-review
description: Use this skill whenever the user asks for a performance review, hot-path audit, profiling guidance, or to "find places to make this faster" in Rust code. Trigger on phrases like "look for performance improvements", "is this on the hot path", "where is this allocating", "review for cache locality", "is there lock contention here", or when the user shares a flamegraph / criterion result and asks for next steps.
---

# Rust Performance Review

You are scanning Rust code for performance-improvement candidates. You are **not** rewriting the code, and you are **not** asserting wins; you are surfacing places where a measurable gain is plausible and proposing the experiment that would confirm or refute it.

This work has a hostile prior: most "obvious" optimizations either do not matter (cold code), do not help (the compiler already handles it), or trade real complexity for an imagined win. Be skeptical. Each candidate must answer: where on the call graph does this run? what would we measure? how much complexity does the fix add?

The review is split across focused per-category checklists in `ia/skills/performance_review/`. **You are the orchestrator**: you triage which categories apply to the scope, dispatch one sub-agent per category in parallel, then synthesize their findings into a single report. Per-category rules are not duplicated in this file — read each category file when dispatching the corresponding sub-agent.

## Review principles (must always hold)

- **Profile first; pattern-match second.** A finding without a profile, benchmark, or strong call-graph evidence that the code is on a hot path is filed at low priority and paired with the measurement that would promote it. Do not propose a rewrite of code you cannot show is hot.
- **Every candidate has a measurement plan.** State the experiment: which benchmark, which profiler, which metric (wall time / allocations / cache misses / lock-wait time / syscalls), and what threshold makes the fix worth merging. "It will be faster" is not a finding.
- **Complexity is a cost.** Every fix names the complexity it introduces (extra type, lifetime gymnastics, `unsafe`, build-config knob, dependency) and weighs it against the expected gain. A fix that doubles maintenance for a 2% wall-time win is a bad trade in critical lab code.
- **Hot-path discipline.** Optimizations in code that runs once at startup, in a CLI flag handler, or in error-handling paths are noise. Be explicit about the call-frequency assumption and downgrade severity when call frequency is unverified.
- **Correctness first.** Never recommend an optimization that weakens invariants, introduces `unsafe`, or relaxes atomic ordering without a separate, justified safety review. Performance findings that touch correctness boundaries are flagged as such and held to the evidence bar of a correctness review.
- **One change per measurement.** Bundling allocator switch + LTO + a code refactor in the same PR produces an unreadable result. Each PR names the single hypothesis being tested.

The severity rubric and per-finding format are defined in `ia/skills/performance_review/_finding_format.md`. Read it once at the start of every review — both you (for synthesis) and every sub-agent you dispatch will follow it.

## Review procedure

You run each numbered step once per review. Only step 5 fans out into parallel sub-agents.

### 0. Read `PROJECT_STATUS.md`

Read `PROJECT_STATUS.md` at the project root to orient on what the project is, the current focus, and the in-scope feature's prior artefacts (plan, implementation report, prior reviews, prior perf reviews, prior benches). See *Project status protocol* below for what to read and what to update at the end.

### 1. Establish scope and call frequency

Determine whether the review covers a full crate, a module, a PR diff, or a single function under benchmark. State the scope. For each in-scope file, identify which functions are believed to be on the hot path and on what evidence: a profile, a benchmark, a call-graph argument, or "no evidence yet". Code with no hot-path evidence is reviewed at lower priority — flag this as a constraint passed into every sub-agent prompt.

### 2. Inventory existing measurement

List every existing benchmark (`benches/`), profile artifact (`flamegraph.svg`, `*.profraw`, `dhat-*.json`, `samply.json.gz`), and any quoted measurement the user has provided. Note what is missing — which functions have no benchmark, no profile, or only synthetic timings. **If no measurement exists for the in-scope code, the first deliverable is a measurement plan, not code rewrites** — say so explicitly in the verdict.

If real benchmark or profile output is available, quote it verbatim. The verbatim output is passed into each sub-agent's prompt at step 5 so they do not re-run.

**If a sampling profile is impossible to collect in the user's environment, raise the alarm before moving on.** This is not a "we'll work around it with other tools" situation — sampling profile access is a load-bearing input to a useful review, and proceeding without it significantly degrades the review's actionable output.

Common blockers (diagnose precisely — the fix depends on the OS):

**Linux** — `perf_event_open` is the kernel ABI all the major sampling tools use (`perf record`, `cargo flamegraph`, `samply`).

- `kernel.perf_event_paranoid >= 3` (Debian's default) blocks `perf_event_open` for unprivileged users. Fix: `sudo sysctl kernel.perf_event_paranoid=2` (persist via `/etc/sysctl.d/99-perf.conf`).
- Rootless container without `CAP_PERFMON`. Adding `--cap-add=PERFMON` on the container invocation only helps in *rootful* podman/docker — in rootless mode the kernel's `capable(CAP_PERFMON)` check runs against the host UID's caps, which the user namespace cannot grant. Lower the host paranoid instead, or run perf against the release binary from the host directly (`target/.../bench_binary --bench ... --profile-time 30`).
- Restricted CI / sandboxed runners. May require running benches on a privileged worker, or accepting that any code-level findings will be marked Speculative.

**macOS (incl. Apple Silicon M-series)** — There is no `perf_event_paranoid` knob; `perf` is Linux-specific. Profiling your own Rust process generally *just works* without elevated permissions. The native toolchain:

- `samply record cargo bench --bench <name> -- <filter> --profile-time 30` — same crate as Linux, cross-platform; outputs to a self-contained viewer (https://profiler.firefox.com). Easiest first try.
- `cargo flamegraph --bench <name> -- <filter> --profile-time 30` — on macOS this auto-switches to a `dtrace` backend.
- **Instruments.app** (ships with Xcode Command Line Tools, free) — the canonical GUI profiler; CPU Profiler template uses M-series hardware counters natively. CLI variant: `xcrun xctrace record --template 'Time Profiler' --launch -- <binary>`.
- `sample <pid>` — the simple, always-available CLI sampler that ships with the OS.

macOS blockers that *do* exist (less common, mention only if hit):

- System Integrity Protection (SIP) restricts some kernel-level DTrace probes; user-space sampling against your own process is unaffected. If `dtrace` complains about probe registration, disabling specific SIP features requires a Recovery-mode `csrutil` change — flag to user, don't recommend casually.
- Profiling a signed third-party binary (cross-process attach) needs the binary's entitlements to allow it. Not relevant for our own `cargo bench` output.

#### Tools verified available on this project's machine (macOS / Apple Silicon dev host)

Confirmed present as of 2026-06-01 (re-check with `command -v` if a run is far in the future — installs drift):

- **`samply`** (`~/.cargo/bin/samply`) — sampling profiler. `samply record --save-only -o out.json.gz -r 2000 -- <binary> …`. **Gotcha:** `--save-only` writes *unsymbolicated* module-relative addresses; symbolication happens in the Firefox-profiler UI (or via a symbol server samply runs). For a *headless* self-time ranking, prefer `sample` (below), which symbolicates inline.
- **`sample` (macOS built-in)** — attach to a running PID: `sample <pid> <secs> -file out.txt`. Produces a **symbolicated** call tree directly. This is the path that worked headlessly for the 2026-06-01 `var_calling` review. Note rayon pool threads show up as idle `Sleep::sleep` frames with the full sample count — filter them out.
- **`xctrace` (Instruments)** — `xcrun xctrace record --template 'Time Profiler' …`. GUI-oriented; `.trace` output is awkward to parse headlessly. Use `sample`/`samply` first.
- **host `cargo` + `rustc` (native arm64-apple-darwin)** — the dev container builds *Linux* binaries (`target-container/`), which cannot be profiled natively on the macOS host and whose in-VM `perf` targets the VM kernel. So for any host profiling you must **build natively on the host** (`cargo bench --no-run` / `cargo build --release --example …`), which lands in `target/` (gitignored) — this is the sanctioned exception to "cargo runs in the container". First native build is ~1 min.
- **DHAT** — no install needed; wired via the `dhat` crate behind the `dhat-heap` feature. `cargo run --release --example dhat_var_calling --features dhat-heap -- …` writes `dhat-heap.json` (open at <https://nnethercote.github.io/dh_view/dh_view.html>, or parse offline — the stacks are deep, so attribute to the first `src/<crate-module>/…rs:line` frame, skipping the `examples/dhat_*` alloc-hook frame and the alloc/core/BTreeMap internals). Sibling examples: `dhat_pileup`, `dhat_psp_reader`, `dhat_psp_writer`, `dhat_baq`. Run on the **host** when the fixture needs `$HOME/genomes/…` (the macOS container mounts only the project tree).
- **`cargo-show-asm`** (`cargo asm …`, installed 2026-06-01) — codegen inspection for the `hot_loops` category: confirm autovectorization / bounds-check elision on a named hot function (`cargo asm --lib --simplify "<crate>::path::to::fn"`). Functions that the profiler sees as distinct symbols are not fully inlined and `cargo asm` will find them; fully-inlined helpers must be inspected at their caller.
- **criterion benches** — `./scripts/dev.sh cargo bench --bench <name> -- <filter>` for reproducible wall numbers (run in the container for the committed baseline; run the native host build only when you also want to `samply`/`sample` it).
- **Not installed on the host** (it's in the container instead — see below): `cargo-flamegraph`, `perf`, `valgrind`, `coz`, `hyperfine`.

#### Tools in the dev container (`./scripts/dev.sh` → Debian 12 / aarch64 — the production target)

Verified 2026-06-06. The container is **Linux aarch64** (the production build target), runs as **root**, `kernel.perf_event_paranoid=2`. **It is ephemeral** — every `./scripts/dev.sh …` is a fresh container; persist any built artifact or output under the bind-mounted project tree (`tmp/`), never `/tmp` (which is wiped). Already baked into the image (`Containerfile`): `perf` (`linux-perf`), `valgrind`, `hyperfine`, `cargo-flamegraph`, `samply`, `cargo-show-asm`. `coz` is in the `Containerfile` but **not yet in the live image** (its build layer lands on the next image rebuild — see note below); until then use the pre-built copy at `tmp/coz/install/bin/coz`.

> **PATH gotcha — use `bash -c`, not `bash -lc`.** A *login* shell sources `/etc/profile`, which resets `PATH` and **drops `/usr/local/cargo/bin`** — so `samply` / `cargo-flamegraph` / `cargo-asm` (all cargo-installed) silently vanish from a `bash -lc` command and appear "missing". `./scripts/dev.sh cargo …` and `./scripts/dev.sh bash -c '…'` both keep the right PATH. Reference genomes/PSPs are reachable via the same absolute host paths (HOME is propagated).
>
> **Image rebuilds:** `./scripts/rebuild-image.sh` (now runtime-aware: podman or Apple `container`). A *warm-cache* rebuild only re-runs changed layers (≈ minutes). A **cold** buildkit cache (fresh builder VM) forces a full from-scratch rebuild of this heavy image (samtools/bcftools/freebayes/GATK/nodejs + the Rust profilers + coz/LIEF) ≈ **40 min** — so don't trigger one casually; the coz/LIEF layer also bundles mbedTLS, which is the bulk of its build time.

- **`perf` CPU self-time sampling — WORKS, and on the production arch.** `perf record -e cpu-clock -F 499 -g -- <binary>` then `perf report`/`perf script`. Use the software `cpu-clock` event (not the default `cycles`).
- **`cargo flamegraph`** (already in the image) — wraps `perf`, so it needs the same `cpu-clock` override (its default `cycles` event returns `<not supported>` here). Working invocation (note `bash -c`, not `-lc`): `./scripts/dev.sh bash -c 'cargo flamegraph -c "record -e cpu-clock -F 997 --call-graph dwarf -g" --bin pop_var_caller -- var-calling --reference … --output /dev/null --threads 1 …'` → writes `flamegraph.svg`. The container is the right home for it (runs as root, no `sudo`; on the macOS host `cargo flamegraph` uses `dtrace`, which needs `sudo` and is awkward headless — prefer host `sample`/`samply` there instead). `samply` likewise wraps `perf` and takes the same `-e cpu-clock`.
- **`coz`** — **causal profiler**; the right tool for the producer→workers→writer pipeline ("speeding up which line actually moves end-to-end throughput?"). Add `#include <coz.h>` + a `COZ_PROGRESS` throughput point (e.g. per chunk written), build with `-g`, run `coz run --- <binary>`; reads `profile.coz`. **It samples on `PERF_TYPE_SOFTWARE`/`SW_TASK_CLOCK`, so it works here despite the missing PMU** (verified end-to-end). Until the next image rebuild bakes it in, run it from `tmp/coz/install/bin/coz` (set `COZ_LIBCOZ_PATH=tmp/coz/install/lib/aarch64-linux-gnu/libcoz.so` if invoked outside its prefix).
- **`valgrind`** — callgrind gives deterministic instruction-count ranking for A/B (slow, ~10–50×); cachegrind is unreliable on Apple Silicon — avoid.
- **`hyperfine`** — statistical CLI wall-clock benchmarking (reference-tool comparisons).

#### What NONE of these environments can do (the PMU is not virtualized in the Apple-`container` VM)

`perf stat -e cache-misses,instructions,…` returns `<not supported>` in the container; consequently **`perf c2c` (false-sharing), `perf sched` / off-CPU blocking analysis, and hardware counters (IPC, cache/branch-miss) are unavailable** on this machine (host or container). This bites the concurrency category most (the "why doesn't wall scale past T≥2" question). Cover it by either:
1. **Instruments on the host** — the CPU Profiler / System Trace templates *do* read M-series hardware counters and thread-state (off-CPU), but are GUI-oriented and awkward to parse headlessly.
2. **Manual channel instrumentation** — time spent blocked in crossbeam `send` (workers can't keep up) vs `recv` (producer can't keep up); a few lines, always works, names the bottleneck stage directly.
3. **A bare-metal Linux box** with full PMU + sched tracepoints if `perf c2c` / `perf sched` are truly needed.

The reference FASTA for tomato fixtures lives at `~/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa`; real per-sample `.psp` cohorts are under `benchmarks/tomato1/results/ours/cohort/psp/` (and `/Users/jose/devel/pop_var_caller/tmp/aligned_psp/` for the 50-sample set).

**Other OSes / WSL** — Treat as the closest match (WSL2 ≈ Linux; FreeBSD has `dtrace`). If you cannot match it, file the gap and ask the user.

**Other profilers are not substitutes.** Each of them answers one narrow question:

- DHAT — allocations only. A function can be 30 % of CPU and 0 % of allocations.
- `cargo asm` — codegen inspection. Tells you whether bounds checks elided or autovec happened, but not whether the affected line matters.
- valgrind / callgrind — deterministic instruction counts. Useful for ranking changes against each other but ~10–50× slower than native; the instruction count is *not* the same as the CPU-cycle cost (cache misses, branch mispredicts, FMA throughput all invisible).
- Criterion wall-clock — cross-run variance is typically 5–10 %, larger than the effect size of individual changes. Within-run confidence intervals look impressive but are tighter than reality.

None of these tells you which line of code holds 30 % of the CPU. Only a sampling profile does. A review without one will produce a long list of pattern-matched candidates of which most will show **no measurable gain when applied**, because the underlying site was not in the top 20 % of self-time. This failure mode is documented end-to-end in `ia/reviews/perf_baq_2026-05-12.md` (round one applied seven Likely findings against pattern-match alone; the criterion bench could not detect a net change; round two used a real profile and got measurable wins on the same code base).

When sampling is blocked, **tell the user explicitly**, in this order:

1. Name what is blocked and why (kernel paranoid setting, missing capabilities, sandboxed runner — diagnose precisely).
2. Name the concrete fix to apply (the `sysctl` line, the container flag, "run on bare metal", etc.).
3. State the cost of *not* fixing it: the review's actionable output is significantly degraded; the candidate list will not be rankable by hot-path evidence and many candidates that look promising will fall through.
4. Ask whether they want to fix the access issue first, or proceed at lower-quality output.

Do not silently route around this with "we'll use DHAT instead" or "we'll use cargo asm and infer." If the user accepts proceeding without sampling, downgrade every code-level finding's severity to **Likely** (no **Hot-path** is reachable without profile evidence per the rubric) and call out the constraint at the top of the report.

### 3. Determine intent and targets

Before judging whether the code is fast enough, establish what it is meant to do — its purpose, its expected input sizes, its latency or throughput target, its target hardware. Performance review is meaningless without these numbers. If the targets are not stated, ask before continuing or file the gap as a Note finding.

Write a one-paragraph "performance intent" summary; it is passed into each sub-agent's prompt.

### 4. Triage categories

Decide which per-category checklists apply. Each lives at `ia/skills/performance_review/<category>.md`.

| Category | Apply when |
|---|---|
| `methodology` | Always. Establishes that profiling, benchmarking, and `Cargo.toml` are sound before any code-level finding is acted on. |
| `allocations` | Code on the hot path constructs `Vec` / `String` / `Box` / `Arc` / map types, clones owned data, calls `format!` / `to_owned` / `to_string`, or has unbounded buffer growth. |
| `data_layout` | Hot path iterates over collections of structs, holds shared atomic state across threads, defines a `pub struct` with several fields, or there is suspicion of cache-miss pressure or false sharing. |
| `concurrency` | Code uses `Arc`, `Mutex`, `RwLock`, atomics, channels, `rayon`, `tokio`, `async fn`, or `spawn`. Skip otherwise. |
| `hot_loops` | Hot path contains tight numeric or byte-processing loops, slice indexing, iterator chains, generic dispatch, `format!`, or anywhere autovectorization or branch layout could plausibly matter. |
| `io_and_syscalls` | Code performs file or socket I/O, reads/writes large data, or makes per-record syscalls. |

When in doubt, dispatch — a sub-agent that finds nothing applicable writes `No findings.` and is cheap.

### 5. Dispatch sub-agents in parallel

Create the scratch directory: `tmp/perf_review_<YYYY-MM-DD>_<scope-slug>/` (append `_v<N>` if it already exists). The slug is a short kebab-case identifier of the reviewed module or PR (e.g. `gvcf_parser`, `pr-142`).

Project rule: scratch space is project-local `tmp/`, never `/tmp`. Add `tmp/` to `.gitignore` if it is not already covered by the existing target ignores.

For each selected category, dispatch a `general-purpose` sub-agent **in parallel** — issue a single message with multiple Agent tool calls. Each sub-agent prompt:

> Run the **<category>** checklist on the following Rust performance-review scope.
>
> **Scope:** <full crate / PR diff / module / single function>
> **Performance intent and targets:** <one paragraph from step 3, including throughput/latency targets and target hardware>
> **In-scope files (full paths):** <list>
> **Hot-path evidence:** <verbatim profile output / benchmark numbers / "none — pattern-match only">
> **Out of scope:** <list, with reasons>
>
> **Instructions:**
> 1. Read `ia/skills/performance_review/<category>.md` for the rules to apply.
> 2. Read `ia/skills/performance_review/_finding_format.md` for the severity rubric and finding format.
> 3. Read each in-scope file.
> 4. Apply each rule and produce findings in the specified format. For every candidate, propose the **measurement plan** (benchmark or profile that confirms the gain) and the **complexity cost** of the fix. Findings without a measurement plan are downgraded.
> 5. Write findings to `tmp/perf_review_<date>_<slug>/<category>.md`. If no findings apply, write only the line `No findings.`
> 6. Do not invent profile output, benchmark numbers, or call frequencies. Cite only what was provided in the prompt or what you read in the source.
> 7. Stay within the category. Issues that belong elsewhere go under a `## Cross-category observations` heading at the bottom of your file.

Substitute `<category>` and the scope fields for each dispatch. Do **not** assign severity codes (H1, L1, …) inside sub-agents — that happens at synthesis.

### 6. Collect findings

When all sub-agents complete, read each `tmp/perf_review_<date>_<slug>/<category>.md`. Tally findings, note any cross-category observations, and decide whether each cross-category note becomes its own finding or merely informs synthesis.

Promote findings whose measurement plans converge (multiple categories agree the same site is hot). Demote findings that turn out to be cold-path on closer inspection. If a sub-agent appears to have skipped its scope or produced unverifiable findings, redispatch that one category with an explicit instruction to fix the gap.

### 7. Synthesize the unified report

Compose the report using the *Output format* below. Verdict, measurement plan, and "What's already good" need the full picture and are produced by you. Assign severity codes during synthesis: `H1, H2, …` for Hot-path; `L1, L2, …` for Likely; `S1, S2, …` for Speculative; Notes stay grouped without numbering. Save to `doc/devel/reports/reviews/perf_<module-slug>_<YYYY-MM-DD>.md` per the saving conventions below. Leave the per-category files in `tmp/` as an audit trail.

### 8. Update `PROJECT_STATUS.md`

After the perf review report is saved, update the in-scope feature's block in `PROJECT_STATUS.md` to point at it. See *Project status protocol* below for the rules.

## Output format

Produce the synthesized report in the following order. Use the section headings verbatim so the format is machine-readable.

### 1. Scope and constraints
- What was reviewed: full crate / module / PR diff / single function.
- Reviewed against: commit hash, branch, or "as-provided".
- Throughput / latency targets, expected input sizes, target hardware.
- Hot-path evidence available (profile, benchmark, or "none — pattern-match only").
- In-scope files (list).
- Deliberately out of scope (list, with reason).
- Categories dispatched (list, each with a one-line reason).

### 2. Verdict
One of:
- **Profile first** — there is not enough hot-path evidence to recommend code changes. Section 3 lists the measurements to take.
- **Run experiments** — candidates exist; their priority and order are listed in section 5.
- **Apply the listed wins** — at least one candidate is well-evidenced (matching profile, plausible mechanism, contained complexity); apply with the proposed measurements as gates.

### 3. Measurement plan
The benchmarks and profiles to add or run, in the order they unblock other findings. Each entry: command, expected output, what threshold answers what question. If the verdict is *Profile first*, this is the primary deliverable.

### 4. Build / toolchain configuration
LTO, codegen-units, panic, opt-level, debug, allocator, target-cpu — anything in `Cargo.toml` or `.cargo/config.toml` that should change before code-level work. Driven by the `methodology` sub-agent's findings.

### 5. Code-level findings
Grouped by severity (**Hot-path** → **Likely** → **Speculative** → **Note**). Within each severity, ordered by confidence (High → Low), then by file. Each finding follows the format defined in `ia/skills/performance_review/_finding_format.md`, with the severity code (H1, L1, …) prepended to the title — e.g. `H1: src/parser.rs:42 — Title`.

### 6. Out-of-scope observations
Performance smells in untouched code, surfaced but not blocking. Each: file, brief description, suggested follow-up (separate PR or issue).

### 7. What's already good
Up to 3 specific, transferable patterns the code is already getting right (e.g., "uses `with_capacity` everywhere it knows the size", "shards by sample to avoid `Mutex` contention"). One sentence each, with file references. No general praise. Skip the section entirely if nothing specific qualifies.

### Author response convention
Address each finding by its identifier (e.g., "H1", "L2") with one of: `applied in <commit>` / `experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. The "experiment shows no gain" path is expected and welcome — that is what the measurement plan is for.

---

Be direct. If something is plausibly hot, say so plainly and propose the experiment. Vague "could be faster" and vague "looks fine" are equally useless.

## Saving the report

### Directory and filename

Save to the project's `doc/devel/reports/reviews/` directory:

```
doc/devel/reports/reviews/perf_<module-slug>_<YYYY-MM-DD>.md
```

Examples:

- `doc/devel/reports/reviews/perf_gvcf_parser_2026-05-10.md`
- `doc/devel/reports/reviews/perf_pipeline_2026-05-10.md`
- `doc/devel/reports/reviews/perf_pr-142_2026-05-10.md`

If a review for the same scope and date already exists, append `_v<N>`:

- `doc/devel/reports/reviews/perf_gvcf_parser_2026-05-10_v2.md`

### Document header

```markdown
# Performance Review: <module-slug>
**Date:** <YYYY-MM-DD>
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** <one-line description>
**Verdict:** <Profile first / Run experiments / Apply the listed wins>
**Hot-path evidence:** <profile / benchmark / pattern-match only>

---
```

The body is sections 1–7 of *Output format* above, in order, with verbatim headings.

### File links inside findings

References to source files use relative Markdown links from the `doc/devel/reports/reviews/` directory:

- Single line: `[file.rs](../../../../path/file.rs#L123)`
- Range: `[file.rs](../../../../path/file.rs#L123-L456)`

Display text is the path (no backticks).

### Pre-save checklist

- [ ] Every Hot-path finding has High confidence and quotes the profile / benchmark output that names the site.
- [ ] Every Likely finding has a measurement plan that would confirm or refute the gain.
- [ ] Every finding has a complexity cost named honestly.
- [ ] Every cited file:line was actually read (no invented locations).
- [ ] No fabricated percentages or speedup multipliers anywhere in the report.
- [ ] Cold-path findings are downgraded to Note unless the orchestrator marked them in scope.
- [ ] Severity codes (H1, L1, S1) are consistent and dense (no gaps).
- [ ] Per-category files in `tmp/perf_review_<date>_<slug>/` are left in place as an audit trail.
- [ ] Build configuration findings (section 4) are separated from code-level findings (section 5).
- [ ] `PROJECT_STATUS.md` updated (per *Project status protocol*).

## Project status protocol

The project tracks the lifecycle of every feature in `PROJECT_STATUS.md` at the project root. It is a navigation aid, not a source of truth — use it to find the relevant spec, plan, and prior reports for the in-scope feature, then verify against current code as usual.

**At task start.** Read `PROJECT_STATUS.md`. The immutable "About this project" paragraph (delimited by `ABOUT-PARAGRAPH-START` / `ABOUT-PARAGRAPH-END` HTML comments) gives the design context and points at the authoritative spec; "Current focus" confirms the project's direction and last-completed work; the per-feature blocks point at the plan, prior impl reports, and prior reviews (correctness and performance) for the in-scope feature.

**At task end** (after the perf review report has been saved): update only the in-scope feature's block in `PROJECT_STATUS.md`.

- Append a link to the new perf review under `Latest reviews:` (perf reviews live alongside correctness reviews in the same bullet group; the filename prefix `perf_` is enough to tell them apart). Prefer to replace the previous perf review link rather than accumulate a long list.
- Leave `Status:` unchanged — perf reviews do not advance the correctness lifecycle. Exception: if the verdict is `Profile first` and that blocks further code-level work, add an `Open:` item naming the measurement to take.
- Add `Open:` items for any Hot-path or Likely findings the run surfaced; do not close existing items.
- Refresh **Current focus** — rewrite `Last completed task` to name this perf review and link the report. Touch `Next task` only if the human PM has not already set one; otherwise leave it alone, optionally appending `(suggested follow-up: …)` after the existing text.

**Status vocabulary:** `planned` / `in-flight` / `implemented` / `reviewed` / `fixes-applied` / `shipped` / `superseded`. (Perf reviews do not have their own status — they annotate whichever status the feature is already in.)

**Hard rules.**

- Do not edit the **About this project** paragraph or anything between the `ABOUT-PARAGRAPH-START` / `ABOUT-PARAGRAPH-END` comments.
- Do not modify another feature's block.
- Do not summarize the perf findings inside the block — the block is a list of pointers; numbers and measurement plans live in the saved report.
- If the in-scope feature has no block yet, create one using the format of existing blocks.
- If `PROJECT_STATUS.md` and the current code disagree, trust the code; the status file is stale and should be updated, not relied on.

## Reusable prompt template

Use this to invoke the skill consistently. Fill in the Context block; the rest defers to the skill body.

> Perform a Rust performance review per the **rust-performance-review** skill. Follow its principles, procedure, severity rubric, and output format in full — do not abbreviate or skip sections.
>
> **Context**
> - **Scope:** <files / module / PR diff / branch>
> - **Performance intent and targets:** <what the code does, on what input sizes, with what latency / throughput target, on what hardware>
> - **Hot-path evidence:** <flamegraph path / criterion bench output / DHAT report / "none — pattern-match only">
> - **Constraints:** <MSRV, `no_std`, target platforms, deadline pressure>
> - **Out of scope:** <legacy modules being deleted, vendored deps>
> - **Prior review history:** <previously reviewed? known tracked issues?>
>
> **Anti-hallucination contract.** Quote tool output verbatim. If a benchmark or profile was not run / read, say so under "Hot-path evidence" and downgrade findings accordingly. Never invent measurements, percentages, file paths, or line numbers.
>
> **Reminders of the most-violated rules** (not a substitute for the per-category checklists):
> 1. Profile first; pattern-match second.
> 2. Every finding has a measurement plan and a complexity cost.
> 3. Cold code is not the hot path; downgrade findings filed against it.
> 4. Correctness wins over speed. Route correctness-adjacent findings through a separate review.
> 5. One change per measurement — do not bundle allocator + LTO + refactor.
> 6. **If the environment blocks sampling profilers (`perf record` / `cargo flamegraph` / `samply` / Instruments / `dtrace`), halt and tell the user before continuing.** This is not a "we'll use DHAT instead" situation — DHAT, `cargo asm`, valgrind, and wall-clock criterion each answer one narrow question and **none** of them tells you which line owns the CPU. Without a sampling profile, most pattern-matched Likely findings will show no measurable gain when applied. The fix depends on OS (Linux: `sudo sysctl kernel.perf_event_paranoid=2`; macOS: usually no fix needed — `samply` / Instruments work without elevated perms); see step 2 for the full per-OS handling.

---

## Sources

The two articles that motivated this skill:

- *Inside Rust's std and parking_lot Mutexes: Who Wins?* — https://blog.cuongle.dev/p/inside-rusts-std-and-parking-lot-mutexes-who-win
- *About memory pressure, lock contention and data-oriented design* — https://mnt.io/articles/about-memory-pressure-lock-contention-and-data-oriented-design/

Background reading the per-category rules draw on:

- *The Rust Performance Book* (Nicholas Nethercote) — https://nnethercote.github.io/perf-book/
- *How to avoid bounds checks in Rust (without unsafe!)* — https://shnatsel.medium.com/how-to-avoid-bounds-checks-in-rust-without-unsafe-f65e618b4c1e
- *Why my Rust benchmarks were wrong* (Guillaume Endignoux on `black_box`) — https://gendignoux.com/blog/2022/01/31/rust-benchmarks.html
- *Optimization adventures: making a parallel Rust workload faster with (or without) Rayon* — https://gendignoux.com/blog/2024/11/18/rust-rayon-optimized.html
- `crossbeam_utils::CachePadded` docs — https://docs.rs/crossbeam-utils/latest/crossbeam_utils/struct.CachePadded.html
- `bumpalo` — https://github.com/fitzgen/bumpalo
- `dashmap` — https://github.com/xacrimon/dashmap
- `rustc-hash` (FxHash) — https://github.com/rust-lang/rustc-hash
- `mimalloc` — https://github.com/microsoft/mimalloc
- `flamegraph-rs` — https://github.com/flamegraph-rs/flamegraph
- `criterion.rs` — https://github.com/bheisler/criterion.rs
