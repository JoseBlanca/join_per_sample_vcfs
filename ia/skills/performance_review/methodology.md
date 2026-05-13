# Methodology & build configuration checklist

**Purpose.** Before any code-level finding is acted on, the team has a benchmark that exercises the function, a profile that confirms it is hot, and a build configuration that is not silently dropping speed on the floor. Findings in this category often outweigh anything in the code-level categories — they are the foundation, run them first.

**Triggers.** Always.

**Skip when.** Never skipped.

## Rules

- **No optimization without a benchmark.** A function that is supposedly hot but has no `criterion` benchmark in `benches/` (or equivalent harness) is a finding. The benchmark is the first deliverable; code changes follow it. A consistent naming convention (`<module>_<workload>_bench`) makes baselines easy to track.
- **Use `std::hint::black_box`.** Benchmarks must wrap inputs in `black_box(...)` and consume outputs the same way; constant-folded benchmarks are the most common cause of "wins" that disappear in production. Any `b.iter(|| f(constant))` without `black_box` is a Likely finding.
- **Profiles must come from release builds with debug info.** Add to `Cargo.toml`:
  ```toml
  [profile.release]
  debug = true
  ```
  or a dedicated `[profile.bench-with-debug]`. Profiling a stripped `--release` binary gives unreadable flamegraphs.
- **Use sampling profilers, not just timers.** `perf record` + `flamegraph` (Linux) or `samply` (cross-platform) is the default for "where does time go". `cargo flamegraph --bench <name>` is the shortest path. Wall-time deltas alone do not localize the bottleneck.
- **Use heap profilers for allocation hypotheses.** `dhat-rs` or `heaptrack` answer "how many allocations, of what size, from where". Wall time alone cannot tell you whether a fix removed allocator pressure or moved it.
- **Benchmark before *and* after, on the same machine, in the same configuration.** Quote both numbers in the PR. Reviewers cannot validate "this is faster" without the before number; criterion's `--save-baseline` and `--baseline` are the supported flow.
- **A single tight within-run CI is not a trustworthy cross-commit number.** Run-to-run variance from system load can shift criterion's median by 20–25 % even with the within-run CI in the 1 % range — criterion's CI describes one run, not the population of runs. Any cross-commit comparison `≥ 5 %` that would trigger a code rollback needs *either* two clean measurements per side, *or* a **revert experiment**: actually revert the change, re-bench, and confirm both directions agree on sign and magnitude. The revert experiment is the cheapest insurance against rolling back a real win that happens to compare against an unlucky-fast prior measurement.
- **Build configuration is part of the program.** Audit `Cargo.toml` and `.cargo/config.toml` for `[profile.release]` settings. Whole-program inlining via `lto = "fat"` and `codegen-units = 1` can produce wins that dwarf any code-level tweak; do not debate inline hints when these are unset.
  ```toml
  [profile.release]
  lto = "fat"          # whole-program inlining; +compile time
  codegen-units = 1    # one unit, more inlining opportunity
  panic = "abort"      # smaller, slightly faster, but no unwinding
  opt-level = 3
  debug = true         # for profiling; otherwise debug = "line-tables-only"
  ```
  Each of these is a separate experiment. `lto = "thin"` is a cheaper middle ground.
- **Profile-guided optimization (PGO) is on the table for stable critical paths.** For a binary that runs the same workload repeatedly (lab pipeline), PGO is the highest-leverage build-time change after LTO. The recipe is documented in `rustc --profile-use` / `cargo-pgo`; cost is a more elaborate build flow that runs the binary on a representative workload between two compiles.
- **Allocator choice is a one-line change.** Trying `mimalloc` or `jemallocator` as the global allocator behind a feature flag is cheap. The mechanism behind the win on multithreaded workloads is per-thread caches that eliminate lock contention inside the allocator; the size of the win is workload-specific. A workload with heavy allocation traffic and no allocator A/B benchmark is a Likely finding.
  ```rust
  #[global_allocator]
  static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
  ```
- **`target-cpu` is set when the deployment target is known.** `RUSTFLAGS="-C target-cpu=native"` or a fixed microarchitecture (`x86-64-v3`, `znver3`, …) in `.cargo/config.toml`. Required for AVX2/AVX-512 autovectorization. **Do not** use `native` for a binary distributed to other machines.
- **Microbenchmarks lie about cache behavior.** A benchmark that calls one function repeatedly with the same input fits in L1 and looks artificially fast. Real workloads thrash. A single-input criterion bench standing in for a streaming workload is a Likely finding with the recommendation to add a streaming bench.
- **`iter_batched` setup is *sampled* under `--profile-time`.** Criterion times only the body, but a sampling profiler (`samply`, `perf record`, `cargo flamegraph`) records every on-CPU instruction in the whole bench routine — including per-iter `setup()`. A heavy fixture clone in setup will dominate the profile and starve the actually-hot code of samples. Prefer building fixtures once outside `b.iter` and borrowing them per iteration. If a per-iter setup is unavoidable, run a contamination sanity-check: the top-level inclusive-time hit on the bench thread should be the function under test (`≥ 50 %`), not the setup closure (`c_with_alloca`, `iter_batched::{{closure}}`, `Routine::profile`). Anything else means the profile is unusable for code-level findings.
- **Benches that align with the API's internal state need an in-bench verification assertion.** A prime-loop that "fills up to a flush boundary", a setup that "warms past the first reallocation", an `iter_batched` body that "triggers exactly one syscall" — each is making a load-bearing claim about internal state, and off-by-one boundary logic produces nominal numbers that look reasonable while measuring the wrong work. The failure is silent unless something dramatic moves the number (e.g. a 1000× speedup that catches your eye). Verify in the bench: assert on a return value, on a `#[doc(hidden)] pub` peek into internal state, or on an after-the-fact counter — whatever proves the body did what the bench's name claims. The assertion belongs in the bench file, not in a manual one-off check, because it's the only thing that survives the next refactor.
- **Compiler version drift is real.** Pin `rust-toolchain.toml` if absolute numbers matter; nightly autovectorization decisions change between versions.
- **One change per measurement.** Each PR names the single hypothesis being tested. Bundling allocator + LTO + refactor produces an unreadable result.

## Output

Findings in this category often go in section 4 of the synthesized report ("Build / toolchain configuration") rather than section 5 ("Code-level findings"). Flag which is which in your output: prepend `[build]` to the title for build-config findings, `[bench]` for missing-or-wrong-benchmark findings, `[profile]` for missing-profile findings.
