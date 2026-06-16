# Performance Review: ssr-pileup

**Date:** 2026-06-16
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** `src/ssr/pileup/` — the SSR Stage-1 pileup (per-sample evidence extraction; the pair-HMM realignment hot path)
**Verdict:** Apply the listed wins — **wins applied (see Author response): ≈−54% on the realignment hot path, byte-identical**
**Hot-path evidence:** criterion bench (`benches/ssr_pileup_perf.rs`) + symbolicated `sample` profile (host macOS M-series) — both new this review

> **Update 2026-06-16 (applied).** H1, H2, and a new finding **P1 (shared-prefix
> DP)** that emerged while implementing them are applied to
> [pair_hmm.rs](../../../../src/ssr/pileup/pair_hmm.rs); each is gated by a
> bit-identity test (`score_candidates_is_bit_identical_to_per_candidate_forward`)
> and the criterion bench. Combined effect on `ssr_realign/window/10`:
> 276.9 ms → ~125.9 ms (**−54%**); `period/tetra` 430 → ~200 ms. L2 and an
> `−∞`-skip variant were tried and reverted (no gain / slower). See the **Author
> response** section at the end for the per-finding outcome and measured deltas.

---

## 1. Scope and constraints

- **What was reviewed:** the module `src/ssr/pileup/` — driver, fetch, triage, candidate generation, the pair-HMM forward scorer, count_repeats, and the per-locus aggregator.
- **Reviewed against:** branch `ssr-pileup-review` off `main` (`b8ad8b4`).
- **Throughput / input sizes / hardware:** Stage 1 genotypes ONE sample's BAM/CRAM against an SSR catalog of ~1–2 M microsatellite loci genome-wide. Each locus is fetched (indexed query), reservoir depth-capped (default cap 1000), and **every spanning read** (typical WGS depth ~20–50×) is realigned against `2·window + 1` candidate haplotype rungs (default `window = 10` → 21 rungs) with a log-space pair-HMM forward. Dev host: macOS Apple Silicon (6 fast cores). Production target: 32-core aarch64 Linux. No formal throughput target is stated; the realignment is the dominant CPU cost (this review establishes that).
- **Hot-path evidence available:** Yes — this review created the measurement (there was none before). A criterion microbench over the realignment core (`bench_harness` seam) and a 30 s `sample` profile of a tight-loop driver (`examples/profile_ssr_pileup.rs`). Verbatim numbers in §3.
- **In-scope files:** `driver.rs`, `fetch_reads.rs`, `triage.rs`, `read_analysis.rs`, `candidate_generation.rs`, `pair_hmm.rs`, `count_repeats.rs`, `locus_record.rs`.
- **Deliberately out of scope:** `src/bam/segment_reader.rs` (the fetch backend — read enough to judge the fetch/concurrency findings, flagged not patched), the SSR catalog (Stage 0) and `.ssr.psp` container schema, and the deferred off-ladder candidate path (not yet wired).
- **Categories dispatched:** methodology (always), hot_loops (the DP — primary), allocations (per-read candidate Vecs), data_layout (the `[f64;3]` DP cells), concurrency (rayon batch model), io_and_syscalls (per-locus indexed fetch).

## 2. Verdict

**Apply the listed wins.** The profile is decisive and one-sided: ~74.6 % of self-time is in `exp` + `log` from `ln_sum_exp2` / `ln_sum_exp3` in [pair_hmm.rs](../../../../src/ssr/pileup/pair_hmm.rs). Two of those calls are provably redundant work (H1, H2) — local, low-complexity, and very close to byte-identical. They are well-evidenced enough to implement now with the bench as the gate. Everything else (banding, bounds-check elision, the allocation cleanups, the driver batch-barrier, the missing aarch64-Linux SIMD floor, the absent end-to-end/fetch bench) is a Likely/Speculative experiment with a named measurement.

The single largest *raw* lever — the `2·window + 1` candidate multiplier (cost is perfectly linear in rung count) — is a **design/calibration decision, not a perf fix**: it touches the realign-everything contract and Stage 2's stutter support points. It is recorded (§6) for the SSR design owners, not actioned here.

## 3. Measurement plan

The measurement now exists; this section is what to run to confirm each finding.

**Baseline already captured (host, native `--release`, M-series):**

```
# criterion: ssr_realign/* — one synthetic locus, depth 30, flank 50, through bench_harness
ssr_realign/period/homopolymer  214.44 ms   ssr_realign/window/3   93.49 ms (7 rungs)
ssr_realign/period/di           273.79 ms   ssr_realign/window/10 278.98 ms (21 rungs)
ssr_realign/period/tri          351.86 ms   ssr_realign/window/15 411.64 ms (31 rungs)
ssr_realign/period/tetra        433.71 ms   ssr_realign/units/8   204.47 ms
                                            ssr_realign/units/30  421.10 ms
# Perfectly linear in rung count: 93.49/7 = 13.36, 278.98/21 = 13.28, 411.64/31 = 13.28 ms/rung.

# sample profile (30 s, dinuc/window 10/depth 30), symbolicated self-time, 25627 in-loop samples:
exp  (libsystem_m.dylib)   9573   37.4%
log  (libsystem_m.dylib)   9545   37.2%   → exp+log ≈ 74.6% (ln_sum_exp2/ln_sum_exp3)
malloc/free/memmove           ~5   ~0.02% (negligible — allocation is cold on the clean path)
```

Commands:
- Wall: `./scripts/dev.sh cargo bench --bench ssr_pileup_perf -- <group>` for the committed (container) baseline; native host build (`cargo bench --bench ssr_pileup_perf --no-run` then run the binary) when also sampling.
- Profile: `cargo build --release --example profile_ssr_pileup`, run it, `sample $! 25 -file tmp/ssr_pileup_sample.txt`.
- Codegen (L2): `./scripts/dev.sh cargo asm` on `forward` to confirm bounds-check branches.
- Allocations (L3/L4): add `examples/dhat_ssr_pileup.rs` (model on `dhat_baq.rs`) driving `analyze_workload` at high `depth`.
- End-to-end / fetch (L5/L7/L9): the missing bench — build `driver::run` over a real catalog + indexed CRAM, report loci/s; add a 3-phase `Instant` accumulator (read/map/write) for L5 and a distinct-chunk-per-locus counter for L7.

**For every H/L finding that changes the DP math, the gate is the same:** re-run the bench (expect the named symbol's share to drop) **and** diff the `.ssr.psp` against the pre-change build. H1/H2 are expected byte-identical-or-last-ULP (confirm the f32 cast in `prune_and_renormalize` absorbs it; if not, it is a QUAL-only delta). Banding/approx-math/linear-space are *not* byte-identical and must gate on genotype concordance over a real catalog.

## 4. Build / toolchain configuration

`[profile.release]` is already correct: `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`, `debug = "line-tables-only"`; `[profile.bench]` and `[profile.profiling]` are sound. `alloc-mimalloc` is wired as an opt-in feature, but the profile shows malloc ~0.02 % — no allocator headroom on this workload; PGO is premature until the end-to-end bench exists.

**One build finding (B1, Likely, High):** [.cargo/config.toml](../../../../.cargo/config.toml) sets a `target-cpu` floor for `(x86_64, linux) → x86-64-v3` and `(aarch64, macos) → apple-m1`, but **no floor for `(aarch64, linux)` — the production target**. The prod binary therefore builds generic armv8-a and the pair-HMM forward misses Neoverse SIMD, while the macOS dev box that produced these numbers gets `apple-m1` — so the bench *flatters* production. The file's own comment says the floors exist precisely for "the HMM inner loops". Add a block matched to the fleet's oldest node (Graviton2/Ampere Altra = `neoverse-n1`, Graviton3 = `neoverse-v1`); confirm the floor is ≤ the oldest prod node before committing. Measure on the aarch64 dev box: build the bench with and without `RUSTFLAGS="-C target-cpu=neoverse-n1"`, compare `ssr_realign/window/10`.

## 5. Code-level findings

### Hot-path

**H1: [pair_hmm.rs:137-146](../../../../src/ssr/pileup/pair_hmm.rs#L137-L146) — `ln_sum_exp2` computes a redundant `exp(0)` on every call**
- **Confidence:** High.
- **Evidence:** `exp (libsystem_m.dylib) 9573 samples = 37.4%` of self-time. `ln_sum_exp2` is the *only* `exp` caller in the module (`ln_sum_exp3` is two nested `ln_sum_exp2`).
- **Mechanism:** after `let m = a.max(b)`, the line `(a - m).exp() + (b - m).exp()` always evaluates `exp(0.0) = 1.0` for whichever argument equals `m` — one of the two `exp` calls is wasted. Rewrite as `m + (other - m).exp().ln_1p()` where `other = a.min(b)`: one `exp` instead of two, and `ln_1p` is more accurate near zero. Roughly half of all `exp` calls are the redundant `exp(0)`.
- **Measurement plan:** re-run `ssr_realign/window/10`; expect a measurable per-rung wall drop. Re-profile; expect the `exp` share (37.4 %) to fall toward ~half. Verify `.ssr.psp` byte-identity (last-ULP `ln_1p` vs `ln(1+x)` should be absorbed by the downstream f32 cast — confirm).
- **Complexity cost:** none structural — a 3-line rewrite of one `#[inline] fn`; keep the two `NEG_INFINITY` short-circuits.
```rust
#[inline]
fn ln_sum_exp2(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY { return b; }
    if b == f64::NEG_INFINITY { return a; }
    let m = a.max(b);
    let other = a.min(b);
    m + (other - m).exp().ln_1p()   // exp(m-m)=1 folded in; one exp, not two
}
```

**H2: [pair_hmm.rs:149-152](../../../../src/ssr/pileup/pair_hmm.rs#L149-L152), [212-217](../../../../src/ssr/pileup/pair_hmm.rs#L212-L217) — `ln_sum_exp3` does 2 `log` (and up to 4 `exp`) per M-cell; a single-pass 3-way reduction cuts it to 1 `log`**
- **Confidence:** High.
- **Evidence:** `log (libsystem_m.dylib) 9545 samples = 37.2%`. The M-cell recurrence (the densest term in the `O(m·n)` DP) calls `ln_sum_exp3` once per interior cell, and `ln_sum_exp3 = ln_sum_exp2(ln_sum_exp2(a,b), c)` — two `.ln()` calls.
- **Mechanism:** take `max(a,b,c)` once, then `m + ln(exp(a−m)+exp(b−m)+exp(c−m))` with the max term folded to `+1.0` (per H1) — one `log` and two `exp` instead of two `log` and up to four `exp`. Halves the `log` calls on the dominant cell. Apply the same to the final `ln_sum_exp3(last[M], last[I], last[D])` (once per read — fold for consistency, not for the win). Mind the `-∞` bookkeeping for three operands.
- **Measurement plan:** after H1+H2, re-profile; expect the `log` share to drop materially, biggest on `period/tetra` and `units/30` (more M-cells). Same byte-identity check.
- **Complexity cost:** one small function rewrite carrying its own three-operand `-∞` handling; no new types/deps.

### Likely

**L1: [pair_hmm.rs:184-228](../../../../src/ssr/pileup/pair_hmm.rs#L184-L228) — band the DP to the diagonal (arch §5.5 `PAIR_HMM_BAND_BP`)**
- **Confidence:** Medium. **Evidence:** cost grows with haplotype length (the `period/*` and `units/*` benches); transcendentals scale with cell count, and the module doc (`pair_hmm.rs:19-25`) explicitly defers banding "until the slow path is shown to bind" — the profile now shows it binds.
- **Mechanism:** the read's extracted region and the candidate `left_flank + motif×L + right_flank` are near-equal length and co-aligned, so forward mass sits near the main diagonal. Restricting `j` to `[i−band, i+band]` skips far-off-diagonal cells and their `ln_sum_exp` transcendentals; saving `∝ (1 − band/n)`, largest on the longest haplotypes (tetra, units-30).
- **Measurement plan:** add a `band` parameter, re-run `period/*` and `units/*`. **Not byte-identical** (drops off-diagonal paths) — gate on genotype concordance over a real catalog, not byte diff.
- **Complexity cost:** a correctness-adjacent approximation + a band parameter to thread and calibrate + boundary handling (out-of-band cells must read as `-∞`). More than a local rewrite; the module doc anticipates it.

**L2: [pair_hmm.rs:206-228](../../../../src/ssr/pileup/pair_hmm.rs#L206-L228) — elide per-cell bounds checks in the inner `for j in 1..=n` loop**
- **Confidence:** Medium. **Evidence:** the ~25 % non-transcendental remainder is the inlined DP cell arithmetic; four indexed slice accesses per cell (`prev[j-1]`, `prev[j]`, `cur[j-1]`, `candidate_seq[j-1]`) carry bounds checks LLVM may not prove redundant (the `Vec`s are grown by `resize_for` but their `len()` isn't visible at the loop).
- **Mechanism:** a leading `assert!(scratch.prev.len() > n && scratch.cur.len() > n && candidate_seq.len() >= n)`, or per-row re-borrowed slices (`&scratch.prev[..=n]`), lets LLVM drop the per-cell checks and straight-line the arithmetic. The `[M]/[I]/[D]` static indices are already fine.
- **Measurement plan:** `cargo asm` on `forward` before/after to confirm `panic_bounds_check` branches disappear; then `ssr_realign/window/10` for a small-but-real delta. Byte-identical (pure codegen).
- **Complexity cost:** low — one `assert!` (a documented invariant) or per-row re-borrow; no `unsafe`. The `mem::swap` rolling rows means the re-borrow must happen inside the `i` loop after the swap.

**L3: [candidate_generation.rs:60](../../../../src/ssr/pileup/candidate_generation.rs#L60), [73-77](../../../../src/ssr/pileup/candidate_generation.rs#L73-L77) — `build_rungs` allocates ~42 Vecs per spanning read, freed next read**
- **Confidence:** Medium. **Evidence:** allocation is cold in the clean-read profile (~0.02 %), so this is **not** hot-path; the bar is a high-depth (near cap 1000) or many-distinct-loci run, named in the measurement plan.
- **Mechanism:** per read, `build_rungs` loops `2·window+1 = 21` rungs; each does `allele.to_sequence(locus)` (a throwaway `tract` Vec) then a `candidate_seq` Vec — ~42 alloc/free per read, ~21 000 avoidable per cap-1000 locus. The `tract` Vec never escapes (it is `extend_from_slice`'d into `candidate_seq` and dropped).
- **Measurement plan:** add `examples/dhat_ssr_pileup.rs` at high depth; confirm `build_rungs`/`to_sequence` as top alloc-blocks sites; act if they are a double-digit % of blocks.
- **Complexity cost / fix:** Step 1 (byte-identical, cheap): an `append_sequence(locus, &mut candidate_seq)` on `Allele` removes the `tract` temporary. Step 2: recycle the inner `candidate_seq` buffers across reads (truncate-don't-drop) since `build_rungs` only appends/clears — collapses per-locus allocs from `O(depth·rungs)` to `O(rungs)`; adds a "inner Vecs are reused" invariant. Gate on the existing determinism tests.

**L4: [locus_record.rs:73](../../../../src/ssr/pileup/locus_record.rs#L73), [105](../../../../src/ssr/pileup/locus_record.rs#L105) — `aggregate`/`prune_and_renormalize` per-read Vec churn**
- **Confidence:** Medium. **Evidence:** pattern-match (not in the profile — the harness exercises scoring, not the aggregator).
- **Mechanism:** `spanning` is `Vec::new()` (reallocs toward `n_spanning`); `prune_and_renormalize` collects a `survivors` Vec then reads it 3× and collects a second output Vec per read.
- **Measurement plan / fix:** pre-size `spanning` to `outcomes.len()` (one-liner, the bound is exact); fuse `survivors` away by accumulating `sum_exp` and the output in one pass (or a `LocusScratch` scratch Vec). Validate via the `assert_normalized` tests.
- **Complexity cost:** low; the fused form adds a scratch parameter.

**L5: [driver.rs:387-411](../../../../src/ssr/pileup/driver.rs#L387-L411) — read→map→write batch barrier idles the pool during serial phases**
- **Confidence:** Medium. **Evidence:** pattern-match (the single-locus bench cannot see a cross-batch effect; no threads sweep exists).
- **Mechanism:** each loop iteration runs serially: fill a ≤8192-locus batch (`read_locus`) → `pool.install(par_iter)` → write the batch in catalog order. While the main thread reads the next batch or writes the previous one, the whole worker pool is parked (and vice-versa) — a fixed serial fraction per batch that does not shrink with thread count (Amdahl), echoing the MEMORY note's "~40 % serial floor, bulk-synchronous" on the SNP path. Stronger at the 32-core prod target than on 6 dev cores.
- **Measurement plan:** end-to-end threads sweep (1/2/4/6) over a multi-batch catalog; wrap the three phases in `Instant` accumulators (no PMU needed). Act if `read+write` is >~15 % of wall at the prod thread count.
- **Complexity cost:** a pipelined producer (read *k+1* / write *k−1* overlapping map *k*) is a topology change that must preserve catalog-order writes (the byte-identity guarantee). Defer behind the sweep.

**L6: [driver.rs:401](../../../../src/ssr/pileup/driver.rs#L401) — `map_init(LocusScratch::new, …)` re-inits the grow-and-keep scratch per rayon job, not per worker**
- **Confidence:** Medium. **Evidence:** pattern-match; contradicts the module's own "each worker reuses its own `LocusScratch` … reuse survives parallelization" claim (`driver.rs:8`, `90-93`).
- **Mechanism:** rayon calls the `map_init` closure once per split/steal *job*, not once per thread, so `PairHmmScratch`'s grow-and-keep rows reset at every split boundary — the first loci of each job re-`resize`. Bounded churn; likely small (the profile shows allocation cold), so this is speculative-in-practice — measure before changing.
- **Measurement plan:** alloc-count over a multi-batch run; or bump an `AtomicU64` in `LocusScratch::new` and compare to thread count (if ≫ threads, confirmed).
- **Complexity cost:** moderate — a `thread_local!<RefCell<LocusScratch>>` or an explicit per-worker fan-out. Only act if the alloc profile disagrees with "cold".

**L7: [fetch_reads.rs:179-190](../../../../src/ssr/pileup/fetch_reads.rs#L179-L190) — one indexed query per locus; adjacent clustered loci re-seek/re-decode the same bgzf block / CRAM container**
- **Confidence:** Medium. **Evidence:** pattern-match (fetch path not profiled).
- **Mechanism:** `get_reads_from_segment` runs ~1–2 M tiny queries genome-wide, each a seek + per-chunk/container decode. Microsatellite loci cluster, so neighbouring loci frequently decode the *same* block twice. Coalescing a run of adjacent loci sharing an index chunk into one query (decode once, demux to per-locus reservoirs by overlap) removes the duplication. The catalog is read sorted so queries are monotonic — this is a duplicated-decode problem, not random-thrash (favourable).
- **Measurement plan:** cheap pre-check first — instrument `get_reads_from_segment` to count distinct `chunk.start()` vs locus count; a ratio ≪ 1 quantifies the headroom. Then `strace -c`/`dtruss` on a real run. Act if seek+inflate ≥ ~15 % of fetch wall and duplication is high.
- **Complexity cost:** a real new `fetch_locus_group` component crossing the catalog-walk/worker boundary; must preserve the per-locus deterministic reservoir seeding + `--threads`-invariant order. Measure before building.

(Build finding **B1** is in §4.) **L9 — no end-to-end / fetch bench** is in §6 as it is a measurement gap, not a code change.

### Speculative

- **S1: [pair_hmm.rs:137-152](../../../../src/ssr/pileup/pair_hmm.rs#L137-L152) — fast/approximate `exp`/`ln`.** The scores are pruned at `AMB_LL_DROP = 4.0` nats then cast to f32, so the precision floor is far above f64 epsilon — a low-degree minimax `exp`/`ln` could replace the libm calls (the 74.6 %). Highest-risk lever: breaks byte-identity, needs a genotype-concordance study, adds a hand-rolled kernel to audit. The `wide` crate is already a dependency (lane-parallel exp/ln) — worth evaluating, but the D-recurrence is serially dependent across `j`, so SIMD-across-cells is blocked; the gain would come from a faster scalar kernel, not vectorization. Separate numerical review.
- **S2: [pair_hmm.rs:165-235](../../../../src/ssr/pileup/pair_hmm.rs#L165-L235) — linear-space forward with per-column rescaling.** Eliminates per-cell `exp`/`log` entirely (multiply-add + ~1 division per column). The most thorough fix and the most invasive — a full recurrence rewrite + rescaling protocol + correctness re-derivation. The code deliberately chose log-space; revisiting it is a design decision. Gate on concordance.
- **S3: [pair_hmm.rs:114-235](../../../../src/ssr/pileup/pair_hmm.rs#L114-L235) — f32 DP cells.** The forward output is already cast to f32 downstream, so the row footprint could halve — **but** the log-sum-exp accumulation over a ~100+ bp haplotype needs the mantissa bits, and pruning compares against a 4-nat threshold, so a few ULP could flip survivors. A correctness/precision question routed to numerical review, not a layout win.
- **S4: [pair_hmm.rs:246](../../../../src/ssr/pileup/pair_hmm.rs#L246) — `score_candidates` result Vec + `Allele::clone()` per candidate.** Free today (on-ladder `Allele` is a `u16` copy), but becomes a `Box<[u8]>` heap clone once off-ladder lands. Do not pre-optimize; re-measure with an off-ladder workload.
- **S5: [driver.rs:398-407](../../../../src/ssr/pileup/driver.rs#L398-L407) — load imbalance from cap-1000 hypervariable loci** in a `par_iter` that assumes uniform per-item cost; a single heavy locus extends the batch's `install` (compounding L5). The cap bounds the worst case (hence Speculative). Mitigations (cost-aware batch ordering) are deferred behind an end-to-end thread-states capture.
- **S6: [pair_hmm.rs:114-116](../../../../src/ssr/pileup/pair_hmm.rs#L114-L116) — `[f64;3]` AoS DP rows.** Examined and concluded **AoS is correct**: the inner loop consumes the whole M/I/D triple of each neighbour, and the D-recurrence reads the cell just written (serial dependence), so SoA's only upside (autovec) is blocked. Filed as a *measured-negative* experiment to document the decision, not assume it.

## 6. Out-of-scope observations

- **Design lever (not a perf fix) — the `2·window+1` candidate multiplier.** [triage.rs:240](../../../../src/ssr/pileup/triage.rs#L240) / [candidate_generation.rs:60-83](../../../../src/ssr/pileup/candidate_generation.rs#L60-L83): every spanning read runs `2·window+1` full forwards regardless of how clean it is; cost is perfectly linear in rung count, so `window` (default 10) is the single biggest raw lever. A clean pure-tiling read has an exact count and `count_pure_tiling` ([count_repeats.rs:43-96](../../../../src/ssr/pileup/count_repeats.rs#L43-L96)) computes it without DP — but it is deliberately **unwired** under the realign-everything design, and narrowing the window touches Stage 2's stutter support-point contract. Route to the SSR design owners, not a codegen change.
- **L9 / measurement gap — no end-to-end bench.** The new bench measures only the realignment core (one locus held hot in L1); the fetch path, the driver batch loop, and cross-loci cache behaviour are unmeasured. Add a `driver::run`-over-a-real-fixture bench (loci/s). The driver tests already build a tiny BAM + hand-written catalog in-process (`driver.rs` `stage1_fixture`) — reuse that to get a dependency-free end-to-end harness without `trf-mod`; a pre-built real catalog is the realism upgrade.
- **`src/bam/segment_reader.rs:614`** — `RecordBuf::default()` allocated fresh inside the BAM per-record decode loop (per read, on the worker path); a reusable cleared `RecordBuf` would remove it. Decode-side, out of module; flagged not patched.
- **[fetch_reads.rs:175](../../../../src/ssr/pileup/fetch_reads.rs#L175), [192](../../../../src/ssr/pileup/fetch_reads.rs#L192)** — `Reservoir::new` allocates `Vec::with_capacity(1000)` per locus even though typical depth fills a few %; a per-worker reusable buffer drained per locus would remove ~1–2 M oversized allocations. Allocations follow-up, behind the high-depth DHAT run.
- **[locus_record.rs:86-97](../../../../src/ssr/pileup/locus_record.rs#L86-L97)** — `chrom: locus.chrom().into()` is one `Box<str>` per locus (~1–2 M small allocs of a handful of distinct contig names); intern or carry a `chrom_id` if a many-loci DHAT run shows it.

## 7. What's already good

- **Grow-and-keep scratch done right** ([pair_hmm.rs:107-132](../../../../src/ssr/pileup/pair_hmm.rs#L107-L132), `PairHmmScratch::resize_for` only grows; `LocusScratch` reused per worker) — `forward` allocates nothing beyond the reused rows, confirmed by the ~0.02 % malloc in the profile.
- **The reservoir is allocation-clean on the hot fetch loop** ([fetch_reads.rs:89-125](../../../../src/ssr/pileup/fetch_reads.rs#L89-L125)) — pre-sized to the cap, replaces in place, moves the `Vec` out; no per-read allocation.
- **Contention-free shared reader path** ([segment_reader.rs](../../../../src/bam/segment_reader.rs)) — the `Mutex<Vec<reader>>` pool is locked only for pop/push at iterator create/drop, never across the scan; index/header are immutable `Arc`s, each call borrows its own reader (no shared cursor). The right use of `std::sync::Mutex`.
- **Correct durable-write shape** ([driver.rs:367](../../../../src/ssr/pileup/driver.rs#L367), [407-413](../../../../src/ssr/pileup/driver.rs#L407-L413)) — 64 KiB `BufWriter`, batched writes, single end-of-stage `fsync` + atomic rename (no per-record syscall or fsync).

### Author response

Measured on the host (M-series, native `--release`), `bench_harness` workload.
Baselines: `base` = pre-change; `h1` = post-H1. All applied changes are gated by
the bit-identity test `score_candidates_is_bit_identical_to_per_candidate_forward`
plus the existing `pair_hmm` / `driver` suites.

- **H1 — applied.** Fold the redundant `exp(0)` in `ln_sum_exp2` (now `m + (other−m).exp().ln_1p()`). Measured vs `base`: window/10 **−8.2%**, di −5.0%, tetra −4.5%. Also more accurate (`ln_1p`). Kept.
- **H2 — applied.** Single-pass `ln_sum_exp3` (one `ln_1p`, max term folded). Measured vs `h1`: window/10 **−10.1%**, di −10.0%, tetra **−14.8%** (largest where haplotypes are longest, as predicted). Kept.
- **P1 (new) — applied. Shared-prefix DP in `score_candidates`.** All `2·window+1` rungs are `left_flank + motif×L + right_flank` and share a long common prefix (left flank + the shorter rungs' tract); the forward DP over a prefix is independent of the bytes after it, so the prefix is scored **once** (its seam column saved per read row) and every rung's tail continues from it. Bit-identical to per-candidate `forward` (proven by the new test). Measured vs `h1`: window/10 **−50.3%**, di −50.4%, tetra **−51%**. This is the single biggest win. Kept. *(Follow-up not yet done: incremental motif-by-motif seam advance would also share the per-rung tract tail — a further ~20%, more complex; deferred.)*
- **L2 — experiment shows no gain, reverted.** Slicing the DP rows to elide per-cell bounds checks measured identical to H2 alone (vs `h1`: −10.0%/−14.9%/−10.1%) — fat-LTO already elided them. An `−∞`-skip variant of `ln_sum_exp3` was *slower* (the branch costs more than the cheap `exp(−∞)` it saves) — reverted.
- **L1 (banding), S1 (approx `exp`/`ln`), S2 (linear-space), S3 (f32 cells) — not applied.** All are approximations that break byte-identity and need genotype-concordance validation over a real catalog; held as flagged experiments (correctness-first). With P1 landed, the remaining `exp` (~48% of self-time) is the natural target for these.
- **L3/L4/L6 (allocations), L5/L7 (driver/fetch), B1 (build floor), L9 (end-to-end bench) — open**, per §4/§5/§6. Allocation is cold on the profiled path (~0.02%); these await the end-to-end / high-depth measurements named in their plans.

Net: combined H1+H2+P1 ≈ **−54%** on the realignment core, all byte-identical and test-gated.
