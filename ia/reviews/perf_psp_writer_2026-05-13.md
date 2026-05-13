# Performance Review: psp_writer
**Date:** 2026-05-13
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the `.psp` per-sample-pileup writer module
**Verdict:** Apply the listed wins
**Hot-path evidence:** samply sampling profile + criterion bench (`benches/psp_writer_perf.rs`)

---

## 1. Scope and constraints

- **What was reviewed:** the on-disk `.psp` writer at
  [src/per_sample_caller/psp/](../../src/per_sample_caller/psp/), with
  the writer in [writer.rs](../../src/per_sample_caller/psp/writer.rs)
  as the focus and the supporting block / varint / index / trailer /
  registry / header build-side modules in scope wherever the writer
  hot path touches them.
- **Reviewed against:** commit `326b58c` on branch `main` (working
  tree: the new bench at
  [benches/psp_writer_perf.rs](../../benches/psp_writer_perf.rs) added
  for this review; three `pub fn new` constructors added to
  `PileupRecord` / `AlleleObservation` / `AlleleSupportStats` so the
  bench can build records from outside the crate).
- **Throughput / latency targets, input sizes, target hardware:**
  the writer is Stage-1's sink. Production input: ~3 × 10⁹
  per-position pileup records per human WGS sample at 30× coverage.
  Target throughput: ≥ 1 M records/s sustained (to keep up with the
  upstream pileup walker). Current bench measures 6.6–6.8 M records/s
  on a SNP-only single-chromosome workload — comfortable headroom
  *today*, but the bench shape is narrow. Target hardware:
  developer-class Linux/x86_64, AVX2 floor already configured via
  `.cargo/config.toml`'s `target-cpu=x86-64-v3`.
- **Hot-path evidence available:**
  - Criterion bench result for `psp_writer/snp_typical_3_3M`:
    ```
    time:   [482.65 ms 491.64 ms 498.62 ms]
    thrpt:  [6.6183 Melem/s 6.7123 Melem/s 6.8373 Melem/s]
    ```
  - samply sampling profile of the bench binary in
    `tmp/psp_writer_perf.samply.json.gz` (26 335 on-CPU samples on
    the bench thread, ~25 s wall). Full resolved hot-list in
    [tmp/perf_review_2026-05-13_psp_writer/_evidence.md](../../tmp/perf_review_2026-05-13_psp_writer/_evidence.md);
    the dominant share inclusive on the writer's call graph is
    35.28 % `RawVecInner::grow_amortized` / 35.39 % `__libc_realloc`
    / 37.70 % `Vec::push_mut` underneath `write_record`, 17.43 %
    `zstd::CCtx::compress_stream` (with 7.12 % `compressBegin_internal`
    + 5.68 % `cwksp_clean_tables` on the per-call setup path), and
    10.09 % self-time inside `PspWriter::validate_record`.
- **In-scope files:**
  - [src/per_sample_caller/psp/writer.rs](../../src/per_sample_caller/psp/writer.rs)
  - [src/per_sample_caller/psp/block.rs](../../src/per_sample_caller/psp/block.rs)
  - [src/per_sample_caller/psp/varint.rs](../../src/per_sample_caller/psp/varint.rs)
  - [src/per_sample_caller/psp/index.rs](../../src/per_sample_caller/psp/index.rs)
  - [src/per_sample_caller/psp/trailer.rs](../../src/per_sample_caller/psp/trailer.rs)
  - [src/per_sample_caller/psp/registry.rs](../../src/per_sample_caller/psp/registry.rs)
  - [src/per_sample_caller/psp/header.rs](../../src/per_sample_caller/psp/header.rs) (writer-side build path only)
  - [src/per_sample_caller/psp/mod.rs](../../src/per_sample_caller/psp/mod.rs)
  - [benches/psp_writer_perf.rs](../../benches/psp_writer_perf.rs)
- **Deliberately out of scope:**
  - The reader path — not implemented yet.
  - The pileup walker and BAQ engine — separate prior reviews.
  - The header TOML parser/serializer — one-shot per file, cold.
  - Test code under `#[cfg(test)]`.
  - CRAM input layer / FASTA fetcher — not on the writer call graph.
- **Categories dispatched:**
  - `methodology` — always; the bench was added this session and its
    coverage / build configuration / heap-profiler tooling were never
    audited before.
  - `allocations` — profile shows 35 %+ inclusive realloc traffic under
    `write_record`.
  - `data_layout` — `Vec<Vec<SlotId>>` ragged 2-D in `BlockAccumulator`
    plus the SoA-vs-AoS question across the 12 column buffers.
  - `hot_loops` — `validate_record` 10.09 % self, `encode_u64_leb128`
    1.22 % self, `encode_list_column` 1.14 % self.
  - `io_and_syscalls` — `__brk` / `__glibc_morecore` / `systrim` cluster
    visible in the profile is driven by the per-column zstd encoder
    lifecycle, and the writer has no production caller yet so the
    BufWriter contract needs to be stated before one is wired up.
  - `concurrency` — **skipped.** The writer is fully single-threaded;
    no `Arc`, `Mutex`, atomics, channel, `rayon`, or `async fn`. The
    pipeline as a whole is parallel at the per-sample level (one
    writer per worker), but each `PspWriter` instance is owned by one
    thread.

## 2. Verdict

**Apply the listed wins** — H1–H4 are all High-confidence, well-evidenced
candidates with contained complexity and clear measurement plans. They
cover three independent axes (Vec capacity, ragged-2D collapse, zstd
encoder lifecycle) so they can land as three separate PRs each gated
by the existing bench, plus a fourth PR for the `validate_record`
ACGTN lookup table. Before any of them lands, the bench harness needs
the two additions in §3 below (phase-chain-heavy workload + DHAT
example) so the inclusive-realloc and zstd-setup costs are visible in
both wall-clock and allocation-count terms.

## 3. Measurement plan

Run before any code-level work; the listed wins below cite the same
benches.

1. **Add `examples/dhat_psp_writer.rs`** mirroring
   [examples/dhat_pileup.rs](../../examples/dhat_pileup.rs) — the
   `dhat-heap` feature is already wired in
   [Cargo.toml](../../Cargo.toml#L34-L38). Run inside container:
   `./scripts/dev.sh cargo run --release --example dhat_psp_writer
   --features dhat-heap`. Save the baseline `dhat-heap.json`. This is
   the primary instrument for H1, H2, and the L-series scratch-buffer
   findings — wall-clock alone will not show whether a per-record
   `Vec::push` site was eliminated, only DHAT names it.

2. **Add phase-chain-heavy and multi-allele bench workloads** to
   [benches/psp_writer_perf.rs](../../benches/psp_writer_perf.rs):
   - `phase_chain_heavy_1M`: 1 M records, average 5–20 active chain
     markers per record across `new_chains` / `expired_chains` /
     per-allele `chain_slots`. Exercises the three `Vec<Vec<SlotId>>`
     columns and the per-record `clone()` calls — both invisible on
     the current SNP-only bench.
   - `multi_allele_500k`: 500 k records with 2–4 alleles each, indels
     mixed in. Exercises the variable-byte bytes column, the per-allele
     scalars, and longer-allele branches of `validate_record`.
   - `streaming_bufwriter_file_1M`: same content as `snp_typical_3_3M`
     against `BufWriter::with_capacity(64 KiB, tempfile)` and a
     sibling against raw `File::create(...)`. Exposes the per-block
     `write_all` syscall cost that `io::sink()` hides today.

3. **Add a `write_record`-only and a `flush_block`-only sub-bench** so
   the per-record append cost is separable from the per-flush
   encode+compress cost. The current single bench averages four
   distinct phases into one number, which dilutes any one-phase fix.
   Sketch: pre-fill a writer to one record short of an auto-flush,
   then `b.iter` either 100 k more `write_record`s (steady state) or
   one boundary-crossing `write_record` (flush triggered).

4. **Capture a fresh samply profile after each change** that claims a
   win. Compare verbatim hot-list against
   [_evidence.md](../../tmp/perf_review_2026-05-13_psp_writer/_evidence.md):
   - H1 should drop `RawVecInner::grow_amortized` inclusive from
     35 % to a small single-digit percentage and `__libc_realloc`
     inclusive correspondingly.
   - H3 should drop `__brk` self from 4.70 % to < 1 % and
     `ZSTD_cwksp_clean_tables` / `ZSTD_compressBegin_internal`
     inclusive from 5.68 % / 7.12 % to near zero.
   - H4 should drop `validate_record` self from 10.09 % to ~3–5 %
     on the SNP-only workload, and scale roughly with the *byte
     count* not the *allele count* on the multi-allele workload.

5. **Run an allocator A/B** once H1+H2 have landed (see Build /
   toolchain Section 4 below). mimalloc behind a feature flag, save
   criterion baselines pre/post, threshold to merge: ≥ 5 %
   wall-clock improvement on `snp_typical_3_3M`.

## 4. Build / toolchain configuration

- **`Cargo.toml` — Global allocator A/B (mimalloc / jemalloc).**
  35 %+ inclusive realloc + ~18 % across `_int_free` / `__brk` /
  `__glibc_morecore` / `systrim` says the workload is heavily
  allocator-bound. The code-level fixes in §5 reduce *traffic* but
  what is left still pays system-malloc's per-call cost. Add behind a
  feature, A/B-bench against the H1+H2+L1+L2 baseline. Complexity:
  one optional dep, one feature, one `#[global_allocator]` in the
  bench / example shim. Cost is sticky: changes per-build
  reproducibility.
- **`Cargo.toml:6-13` — try `[profile.release-abort]` with `panic =
  "abort"`.** With `unwind`, every `?` and every implicit unwind edge
  in `validate_record` (10.09 % self) and `BlockAccumulator::append_record`
  carries landing-pad metadata. Use a derived profile so test
  workflows don't change. Threshold: ≥ 2 % wall-clock improvement.
- **No `rust-toolchain.toml`.** Reproducibility risk for the
  baselines saved this session, not a perf win itself. Add a pin
  before saving any further criterion baselines that will be cited
  externally.
- **PGO** is on the table for the writer (single binary, single
  dominant workload), but defer until the §5 code-level wins are
  exhausted; PGO adds non-trivial build-pipeline complexity for
  ≤ 5–10 % typically.
- **Already correct, no action needed (called out so they are not
  re-litigated):**
  - `[profile.release]`: `lto = "fat"`, `codegen-units = 1`,
    `debug = "line-tables-only"`.
  - `[profile.bench] inherits = "release"; debug = true`.
  - [.cargo/config.toml](../../.cargo/config.toml): `target-cpu =
    x86-64-v3` on Linux/x86_64 (AVX2 + FMA + BMI2 floor) — already
    matches the AVX2 expectation the profile shows in
    `__memcpy_avx_unaligned_erms` and `__memset_avx2_unaligned_erms`.

## 5. Code-level findings

### Hot-path

- [src/per_sample_caller/psp/writer.rs:470-481](../../src/per_sample_caller/psp/writer.rs#L470-L481) — **H1: `BlockAccumulator::new` builds twelve `Vec::new()` columns that grow under per-record `push`**
  - **Confidence:** High
  - **Hot-path evidence:**
    ```
    37.70%   9929  alloc::vec::Vec<T,A>::push_mut [std]   (inclusive)
    35.39%   9319  [libc] __GI___libc_realloc            (inclusive)
    35.28%   9291  alloc::raw_vec::RawVecInner<A>::grow_amortized [std]
    34.69%   9135  PspWriter<W>::write_record [crate]    (inclusive)
    ```
    Most of the realloc / memcpy inclusive cost is underneath
    `write_record`, not `flush_block`. The accumulator's 12 columns
    start at capacity 0 and grow under the amortised-doubling rule;
    for ~330 000 records per 16 MiB block × 12 columns that is ≥ 17
    reallocations per column per block.
  - **Pattern matched:** allocations checklist — *Pre-size containers
    when the size is known or bounded*; *Allocations belong outside
    hot loops*.
  - **Mechanism:** each `Vec::push` on a `Vec::new()`-initialised
    column hits `RawVecInner::grow_amortized` once per power-of-two
    boundary, which is a `realloc(...)` and a `memcpy` of the prior
    contents. Initial capacity sized from `DEFAULT_TARGET_BLOCK_BYTES`
    and the per-record / per-allele projection turns ~17 reallocations
    per column per block into one (or zero, if shape (b) below is
    chosen and the accumulator is reused across flushes).
  - **Measurement plan:** `cargo bench --bench psp_writer_perf --
    --save-baseline pre-H1` before, re-run after with
    `--save-baseline post-H1`. Threshold to merge: ≥ 5 % wall-clock
    improvement on `snp_typical_3_3M`. Pair with a DHAT pass (from §3
    item 1) — number of `realloc` calls should drop by an order of
    magnitude even if wall-clock gains are modest on the small bench.
  - **Complexity cost:** *Shape (a)* — three named constants near
    `DEFAULT_TARGET_BLOCK_BYTES` (`INITIAL_RECORDS_HINT`,
    `INITIAL_ALLELES_HINT`, `INITIAL_ALLELE_SEQ_BYTES_HINT`) and 11
    `Vec::with_capacity(...)` swaps. *Shape (b)* — promote
    `BlockAccumulator` from `Option<...>` to an always-resident
    accumulator with `is_open: bool`, add a `reset(...)` method that
    calls `.clear()` on every column (capacity is retained). Shape (b)
    is structurally cleaner and amortises the first-block
    high-water-mark growth across the whole file.
  - **Suggested experiment / fix:**
    ```rust
    // writer.rs (Shape (a) — minimum diff)
    const INITIAL_RECORDS_HINT: usize = 350_000;
    const INITIAL_ALLELES_HINT: usize = 360_000;
    const INITIAL_ALLELE_SEQ_BYTES_HINT: usize = 1_000_000;

    fn new(chrom_id: u32, first_pos: u32, snapshot_active_slots: Vec<SlotId>) -> Self {
        Self {
            chrom_id, first_pos, last_pos: first_pos, snapshot_active_slots,
            delta_pos: Vec::with_capacity(INITIAL_RECORDS_HINT),
            n_alleles: Vec::with_capacity(INITIAL_RECORDS_HINT),
            new_chain_slots: Vec::with_capacity(INITIAL_RECORDS_HINT),
            expired_chain_slots: Vec::with_capacity(INITIAL_RECORDS_HINT),
            allele_seq_len: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_seq_bytes: Vec::with_capacity(INITIAL_ALLELE_SEQ_BYTES_HINT),
            allele_obs_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_q_sum_log: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_fwd_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_placed_left_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_placed_start_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_chain_slots: Vec::with_capacity(INITIAL_ALLELES_HINT),
            projected_bytes: 0,
        }
    }
    ```
    Land Shape (a) as PR #1; Shape (b) is a follow-up PR if the bench
    still shows allocator pressure under the next flushes.

- [src/per_sample_caller/psp/writer.rs:447-457](../../src/per_sample_caller/psp/writer.rs#L447-L457) + [writer.rs:495-508](../../src/per_sample_caller/psp/writer.rs#L495-L508) + [writer.rs:565-583](../../src/per_sample_caller/psp/writer.rs#L565-L583) — **H2: Collapse the three `Vec<Vec<SlotId>>` list columns to a flat CSR (`data: Vec<SlotId>`, `offsets: Vec<u32>`) representation**
  - **Confidence:** High
  - **Hot-path evidence:** same inclusive-realloc cluster as H1 — even
    the all-empty-list bench workload pushes one empty
    `record.new_chains.clone()` and one `record.expired_chains.clone()`
    per record into the outer `Vec<Vec<SlotId>>`, plus one
    `allele.chain_slots.clone()` per allele. On phase-chain-heavy
    workloads each non-empty inner Vec is a real `malloc`. Encode-time
    pointer-chase visible too: `encode_column` at writer.rs:565-583
    materialises a `Vec<&[u16]>` per list column per flush just to
    satisfy the `&[&[u16]]` signature of `encode_list_column`.
  - **Pattern matched:** data_layout checklist — *ragged `Vec<Vec<T>>`
    is the worst of both worlds (outer SoA spine + inner per-row heap
    allocation)*; allocations checklist — *Clones are evidence, not
    noise*.
  - **Mechanism:** standard CSR layout — replace each
    `Vec<Vec<SlotId>>` with `{ data: Vec<SlotId>, offsets: Vec<u32> }`
    where `offsets[i]..offsets[i+1]` is row `i`. Append becomes
    `data.extend_from_slice(&src); offsets.push(data.len() as u32);`
    — zero allocations per record on the steady-state path. Encode
    becomes one streaming pass over `data` with the varint counts
    derived from `offsets.windows(2)`, removing both the inner-Vec
    pointer chase and the `Vec<&[u16]>::collect()` materialisation in
    `encode_column`. On-disk format is unchanged (the wire is already
    flat per `block::encode_list_column`).
  - **Measurement plan:** primary test on the
    `phase_chain_heavy_1M` workload from §3 item 2. The empty-list
    `snp_typical_3_3M` should be unchanged or slightly improved (one
    fewer outer-Vec `push_mut` per record). DHAT confirms per-record
    `Vec` allocations for the three list columns drop to zero.
    Threshold to merge: phase-chain workload improves measurably **and**
    SNP workload is not slower at criterion's 95 % CI.
  - **Complexity cost:** one new `ListColumn` helper struct
    (~25 lines) in writer.rs; three list-column field replacements;
    a new `encode_list_column_csr(data: &[T], offsets: &[u32], out: &mut Vec<u8>)`
    in block.rs (or fold into a generic that accepts
    `IntoIterator<Item = &[T]>`). No `unsafe`, no new crate dep. The
    wire format is unchanged.
  - **Suggested experiment / fix:**
    ```rust
    // writer.rs — new helper
    struct ListColumn {
        data: Vec<SlotId>,
        offsets: Vec<u32>,   // len = n_rows + 1; offsets[0] = 0
    }
    impl ListColumn {
        fn with_capacity(n_rows: usize, n_slots: usize) -> Self {
            let mut offsets = Vec::with_capacity(n_rows + 1);
            offsets.push(0);
            Self { data: Vec::with_capacity(n_slots), offsets }
        }
        fn push_row(&mut self, row: &[SlotId]) {
            self.data.extend_from_slice(row);
            self.offsets.push(self.data.len() as u32);
        }
        fn rows(&self) -> impl Iterator<Item = &[SlotId]> + '_ {
            self.offsets.windows(2).map(move |w| &self.data[w[0] as usize..w[1] as usize])
        }
    }
    // BlockAccumulator: three list columns replace the Vec<Vec<SlotId>> fields.
    // append_record: replace `.push(record.new_chains.clone())` with
    //                `.push_row(&record.new_chains)` etc.
    // encode_column 0x20/0x21/0x22: encode_list_column(list_column.rows(), &mut out)
    ```
    If H2 is split (just the encoder-side change without the
    accumulator refactor), the minimum-effort version is to make
    `encode_list_column` generic over `IntoIterator<Item = &[T]>` and
    call `.iter().map(Vec::as_slice)` at the call site — that removes
    the `Vec<&[u16]>::collect()` alone. Land the full CSR version if
    the phase-chain bench shows the inner-Vec mallocs matter.

- [src/per_sample_caller/psp/block.rs:375-385](../../src/per_sample_caller/psp/block.rs#L375-L385) — **H3: Reuse one `zstd::bulk::Compressor` across all columns and all blocks**
  - **Confidence:** High
  - **Hot-path evidence:**
    ```
    17.43%   4589  zstd_safe::CCtx::compress_stream [dep]
     7.12%   1875  ZSTD_compressBegin_internal [dep]
     7.12%   1875  ZSTD_CCtx_init_compressStream2 [dep]
     5.68%   1495  ZSTD_cwksp_clean_tables [dep]
     4.70%   1238  [libc] __brk
     4.70%   1239  [libc] __glibc_morecore
     4.62%   1217  [libc] systrim
    ```
    A fresh `zstd::Encoder` per column per flush (~72 encoder
    lifecycles per bench iteration; thousands per real file)
    allocates and frees the CCtx workspace each time. At level 9 the
    workspace is large enough that glibc cannot satisfy it from
    fastbins, so it `sbrk`s up on construction and `systrim`s back on
    drop — the per-call setup is the entire `__brk` / `morecore` /
    `systrim` cluster, and the `compressBegin_internal` /
    `cwksp_clean_tables` samples are the re-initialisation work the
    persistent compressor would skip.
  - **Pattern matched:** io_and_syscalls checklist — *Compressed I/O:
    decompress once into a reusable buffer; per-record decompression
    context creation is a hidden allocation* (symmetrically for
    compression); allocations checklist — *Allocations belong outside
    hot loops*.
  - **Mechanism:** allocate one `zstd::bulk::Compressor` in
    `PspWriter::new`, set `include_checksum(true)` once, call
    `compress_to_buffer(input, &mut out)` per column. The CCtx
    workspace and tables are reused; the per-call cost drops to the
    actual compression work plus a frame-header reset. Pairs naturally
    with finding L1 (writer-owned scratch buffers) — the same
    refactor wires both.
  - **Measurement plan:** same harness; threshold to merge: a fresh
    samply capture shows `__brk` self < 1 % and the
    `compressBegin_internal` / `cwksp_clean_tables` lines fall below
    the 0.5 % SELF reporting bar.
  - **Complexity cost:** one extra `Compressor<'static>` field on
    `PspWriter`; one regression test that compresses a known payload
    twice (fresh-encoder vs reused-context) and asserts the resulting
    bytes are byte-identical — verifies the per-frame content
    checksum (currently set by `include_checksum(true)`) survives
    context reuse. No `unsafe`, no new crate dep.
  - **Suggested experiment / fix:**
    ```rust
    // block.rs
    pub fn zstd_compress_into(
        compressor: &mut zstd::bulk::Compressor<'static>,
        input: &[u8],
        out: &mut Vec<u8>,
    ) -> io::Result<()> {
        out.clear();
        let bound = zstd_safe::compress_bound(input.len());
        out.resize(bound, 0);
        let n = compressor.compress_to_buffer(input, out).map_err(io::Error::other)?;
        out.truncate(n);
        Ok(())
    }
    // PspWriter::new — construct once, store on Self.
    let mut compressor = zstd::bulk::Compressor::new(ZSTD_COMPRESSION_LEVEL)?;
    compressor.include_checksum(true)?;
    ```

- [src/per_sample_caller/psp/writer.rs:248-281](../../src/per_sample_caller/psp/writer.rs#L248-L281) — **H4: `validate_record` per-allele ACGTN inner loop — replace `matches!(b, b'A'|...)` with a `[bool; 256]` lookup table walked via `.all()`**
  - **Confidence:** Medium
  - **Hot-path evidence:**
    ```
    10.09%   2658  PspWriter<W>::validate_record [crate]  (self)
    ```
    On the current bench, every allele is one byte and the loop
    iterates once per allele — so the ACGTN check alone cannot
    *dominate* the 10.09 % (other validation costs share). The
    function is on the per-record path; the byte check is the only
    one inside it whose cost scales with allele length. On
    multi-allele / indel workloads (the `multi_allele_500k` workload
    from §3 item 2) the per-byte branch-chain dominates linearly
    with no protection.
  - **Pattern matched:** hot_loops checklist — *Match on bytes for
    ASCII parsers* (chained-compare → cmp/je chain);
    *Autovectorization needs the compiler's confidence* (data-
    dependent early return through `Result::Err` blocks autovec).
  - **Mechanism:** `matches!(b, b'A'|b'C'|b'G'|b'T'|b'N')` compiles to
    a chain of compare-and-branch instructions, and the `for (j, &b)
    in allele.seq.iter().enumerate() { ... return Err(...); }` shape
    forces a per-iteration early-exit branch. Replacing with `static
    ALLOWED: [bool; 256]` walked via
    `allele.seq.iter().all(|&b| ALLOWED[b as usize])` is one
    load + one branch per byte, the loop body has no early exit, and
    LLVM can unroll/vectorise it. The error path re-scans to recover
    the offending byte index for the message (cold).
  - **Measurement plan:** a microbench comparing the current loop
    against the lookup-table walk for allele lengths 1 / 50 / 500;
    inspect `cargo asm` to confirm the longer-allele path autovectorises
    on `target-cpu=x86-64-v3`. Then re-run `cargo bench --bench
    psp_writer_perf` on `snp_typical_3_3M` (expect a small win) and
    on `multi_allele_500k` (expect a large win, proportional to
    average allele length). Side benefit of the microbench: time
    `validate_record` with the ACGTN body replaced by a no-op to
    attribute how much of the 10.09 % is the byte check vs the rest
    of the function.
  - **Complexity cost:** ~10 lines: a `const ALLOWED: [bool; 256]`
    table, the `.all()` call, and a `#[cold]` re-scan helper for the
    error message.
  - **Suggested experiment / fix:**
    ```rust
    const ALLOWED: [bool; 256] = {
        let mut t = [false; 256];
        t[b'A' as usize] = true;
        t[b'C' as usize] = true;
        t[b'G' as usize] = true;
        t[b'T' as usize] = true;
        t[b'N' as usize] = true;
        t
    };
    // replaces writer.rs:258-265
    if !allele.seq.iter().all(|&b| ALLOWED[b as usize]) {
        let (j, &b) = allele.seq.iter().enumerate()
            .find(|(_, &b)| !ALLOWED[b as usize])
            .expect(".all() failed means a non-ACGTN byte exists");
        return Err(PspWriteError::InvalidRecord {
            record_index,
            reason: format!("allele {i} byte {j} = {b:#04x} (only A/C/G/T/N allowed)"),
        });
    }
    ```

### Likely

- [src/per_sample_caller/psp/writer.rs:368-419](../../src/per_sample_caller/psp/writer.rs#L368-L419) + [block.rs:375-385](../../src/per_sample_caller/psp/block.rs#L375-L385) — **L1: Promote per-flush scratch buffers to writer-scoped reusable `Vec`s**
  - **Confidence:** Medium
  - **Hot-path evidence:** pattern-match on top of H3's allocator
    pressure. `flush_block` 7.17 % inclusive. The
    `__brk` / `morecore` / `systrim` cluster is driven primarily by
    the encoder workspace (H3), but the per-column compressed-buffer
    `Vec<u8>`s, the `payloads`, the `manifest`, the `header_bytes`,
    and the per-column uncompressed `Vec<u8>` in `encode_column`
    (writer.rs:542) are 5–6 small allocations per column ×12 columns
    per flush ×N flushes per file — secondary contributors.
  - **Pattern matched:** allocations checklist — *Allocations belong
    outside hot loops*.
  - **Mechanism:** add five reusable `Vec<u8>` fields on `PspWriter`
    (uncompressed scratch, per-column compressed buffers as
    `Vec<Vec<u8>>` of length `V1_0_COLUMNS.len()`, payloads
    scratch, manifest scratch, header_bytes scratch). `Vec::clear`
    keeps capacity. Pair with H3 — `zstd_compress_into` writes into a
    borrowed `Vec<u8>` directly.
  - **Measurement plan:** primary signal is DHAT alloc-count under
    `flush_block`; threshold to merge: alloc count drops by ≥ a
    factor of 5 in flush-attributed allocations. Wall-clock change
    likely small (paired with H3 will compound visibly).
  - **Complexity cost:** ~5 new fields on `PspWriter`. The
    accompanying `encode_column_into(&def, &block, &mut out)`
    signature change (returning unit, taking `out` by `&mut`) is
    mechanical.
  - **Suggested experiment / fix:** see the sketch in
    [tmp/perf_review_2026-05-13_psp_writer/allocations.md](../../tmp/perf_review_2026-05-13_psp_writer/allocations.md#L264-L284)
    — applied together with H3 in one PR.

- [src/per_sample_caller/psp/writer.rs:378-385](../../src/per_sample_caller/psp/writer.rs#L378-L385) — **L2: Fold the manifest construction into the single column-pass loop in `flush_block`**
  - **Confidence:** Medium
  - **Hot-path evidence:** pattern-match. No isolated samply line —
    the cost is one `.iter().map().collect()` over 12 entries per
    flush, hard to separate. Filed because the fix is mechanical and
    composes with L1.
  - **Pattern matched:** allocations checklist — *`.collect::<Vec<_>>()`
    followed by another iteration is wasted work*. `payloads` is
    iterated twice (once to build `manifest`, once to write to the
    sink).
  - **Mechanism:** push one `ColumnManifestEntry` into a writer-scoped
    `manifest_scratch` Vec inside the same `for column_def in
    V1_0_COLUMNS` loop that builds the compressed payload. Removes
    the second pass.
  - **Measurement plan:** roll into the same DHAT pass as L1.
    Wall-clock impact small.
  - **Complexity cost:** body of `flush_block` becomes one loop
    instead of two passes; ~10-line diff.
  - **Suggested experiment / fix:** see sketch in
    [allocations.md](../../tmp/perf_review_2026-05-13_psp_writer/allocations.md#L373-L388).

- [src/per_sample_caller/psp/writer.rs:393](../../src/per_sample_caller/psp/writer.rs#L393) — **L3: Replace `block.snapshot_active_slots.clone()` with `std::mem::take`**
  - **Confidence:** High
  - **Hot-path evidence:** pattern-match only; the snapshot is small
    and empty in the bench workload. Listed because the `block` is
    consumed by `flush_block` immediately on this line — the clone is
    purely cosmetic.
  - **Pattern matched:** allocations checklist — *Clones are evidence,
    not noise*.
  - **Mechanism:** `let mut block = self.block.take()...` then
    `active_chain_slots: std::mem::take(&mut block.snapshot_active_slots)`.
    Zero behaviour change.
  - **Measurement plan:** bench-neutral. Fold into the H3+L1+L2 PR.
  - **Complexity cost:** none. One-character `let` → `let mut` and
    `.clone()` → `mem::take(...)`.
  - **Suggested experiment / fix:**
    ```rust
    let mut block = self.block.take().expect("flush_block called with no open block");
    // ...
    active_chain_slots: std::mem::take(&mut block.snapshot_active_slots),
    ```

- [src/per_sample_caller/psp/varint.rs:30-36](../../src/per_sample_caller/psp/varint.rs#L30-L36) — **L4: Split `encode_u64_leb128` into an inlined `< 0x80` fast path and a `#[cold] #[inline(never)]` multi-byte slow path**
  - **Confidence:** Medium
  - **Hot-path evidence:** `1.22 % self` in `encode_u64_leb128`. On
    this workload `delta-pos`, `n-alleles`, and `allele-seq-len` are
    all 1, so every varint payload is one byte; ~36 M calls total
    per bench iteration.
  - **Pattern matched:** hot_loops checklist — *Cold paths get
    `#[cold]`*; fast-path / slow-path split.
  - **Mechanism:** the current `while value >= 0x80 { ... }` always
    enters the loop test before the final push. Specialising
    `if value < 0x80 { out.push(value as u8); } else { cold(...); }`
    lets LLVM lay the cold body out-of-line and tightens the inlined
    fast-path into one push.
  - **Measurement plan:** a microbench against 1 M `u64` values drawn
    mostly from `0..128`, plus `cargo asm` on the rewritten function.
    Re-run `cargo bench --bench psp_writer_perf`; merge if ≥ 1 %
    improvement on `snp_typical_3_3M`.
  - **Complexity cost:** one extra fn (the `#[cold]` slow path); no
    new types, no `unsafe`.
  - **Suggested experiment / fix:**
    ```rust
    #[inline]
    pub fn encode_u64_leb128(value: u64, out: &mut Vec<u8>) {
        if value < 0x80 { out.push(value as u8); } else { encode_u64_leb128_cold(value, out); }
    }
    #[cold]
    #[inline(never)]
    fn encode_u64_leb128_cold(mut value: u64, out: &mut Vec<u8>) {
        while value >= 0x80 { out.push((value as u8) | 0x80); value >>= 7; }
        out.push(value as u8);
    }
    ```

- [src/per_sample_caller/psp/writer.rs:248-296](../../src/per_sample_caller/psp/writer.rs#L248-L296) — **L5: Mark `validate_record` error-construction sites `#[cold]`**
  - **Confidence:** Medium
  - **Hot-path evidence:** indirect — `validate_record` is 10.09 %
    self. The cold arms inflate the function's hot-path icache
    footprint.
  - **Pattern matched:** hot_loops checklist — *Cold paths get `#[cold]`*.
  - **Mechanism:** every error arm builds the error with `format!`
    + `String` and never fires in a valid-input run, but LLVM has no
    hint to lay it outside the hot fall-through. Extract each call
    into a `#[cold] #[inline(never)] fn err_invalid(...) ->
    PspWriteError` helper (or wrap in `#[cold]` closures).
  - **Measurement plan:** `cargo bench --bench psp_writer_perf`
    before/after, optional `cargo asm` to confirm the hot block is
    contiguous. Merge if ≥ 2 %. Composes with H4 (same function).
  - **Complexity cost:** a handful of `#[cold]` helpers; no
    semantic change.
  - **Suggested experiment / fix:**
    ```rust
    #[cold] #[inline(never)]
    fn err_invalid(record_index: u64, reason: String) -> PspWriteError {
        PspWriteError::InvalidRecord { record_index, reason }
    }
    // call sites: `return Err(err_invalid(record_index, format!("...")));`
    ```

- [src/per_sample_caller/psp/block.rs:248-255](../../src/per_sample_caller/psp/block.rs#L248-L255) — **L6: Slab-cast list-column emission on little-endian (`bytemuck::cast_slice` or hand-rolled `unsafe` raw-bytes)**
  - **Confidence:** Medium
  - **Hot-path evidence:** `encode_list_column` 1.14 % self —
    essentially all from the per-entry varint(0) on the empty-list
    bench. On non-empty `chain_slots` workloads, the per-element
    `extend_from_slice(&u16::to_le_bytes())` path becomes the inner
    loop.
  - **Pattern matched:** hot_loops checklist — *Autovec needs the
    compiler's confidence*. Per-element 2-byte `extend_from_slice` is
    per-element capacity bookkeeping.
  - **Mechanism:** for `u16` (and the future `u32` etc.) list
    columns, on `target_endian = "little"` the entire inner list is a
    single memcpy via `out.extend_from_slice(bytemuck::cast_slice::<u16,
    u8>(list))`. Hoists cap-test and length update out of the
    per-element loop. Gate on `cfg(target_endian = "little")`; keep
    the per-element loop for big-endian.
  - **Measurement plan:** depends on the `phase_chain_heavy_1M`
    workload from §3 item 2. Microbench `encode_list_column` over
    inner-list lengths 0 / 4 / 32. Merge if the non-empty case
    improves ≥ 10 % and the empty case is no slower.
  - **Complexity cost:** either a `bytemuck` dependency (cheap,
    widely-used) or a localised `unsafe { from_raw_parts(...) }` with
    a `// SAFETY:` comment asserting `T: Pod` and the LE invariant.
  - **Suggested experiment / fix:** see sketch in
    [hot_loops.md](../../tmp/perf_review_2026-05-13_psp_writer/hot_loops.md#L264-L282)
    — gate on `cfg(target_endian = "little")`. Defer until H2 lands
    so the slab call site is `&block.<col>.data[off_a..off_b]` (CSR),
    not `&[&[u16]]` (current).

- (project root) — **L7: Add `examples/dhat_psp_writer.rs` mirroring `examples/dhat_pileup.rs`**
  - **Confidence:** High
  - **Hot-path evidence:** the `dhat-heap` feature is already wired
    in [Cargo.toml:34-38](../../Cargo.toml#L34-L38); examples exist for
    `dhat_pileup` and `dhat_baq` but not for the writer. The writer's
    allocator pressure is 35 % inclusive — DHAT is the right
    instrument to gate H1, H2, L1, L2.
  - **Pattern matched:** methodology checklist — *Use heap profilers
    for allocation hypotheses; wall time alone cannot tell you
    whether a fix removed allocator pressure or moved it*.
  - **Mechanism:** N/A — this is infrastructure.
  - **Measurement plan:** N/A. The example *is* the measurement plan
    for H1/H2/L1/L2.
  - **Complexity cost:** ~30 lines mirroring the existing examples.
  - **Suggested experiment / fix:** copy
    [examples/dhat_pileup.rs](../../examples/dhat_pileup.rs) and adapt
    its harness to drive `PspWriter::new` + 1 M `write_record` calls
    + `finish` inside a `dhat::Profiler::new_heap()` block. Size the
    workload so ≥ 2 block flushes are triggered (i.e. one full block
    + one partial).

- [benches/psp_writer_perf.rs](../../benches/psp_writer_perf.rs) — **L8: Add `phase_chain_heavy_1M`, `multi_allele_500k`, and `streaming_bufwriter_file_1M` workloads**
  - **Confidence:** High
  - **Hot-path evidence:** the bench file currently has one workload
    (`snp_typical_3_3M`). H2's payoff lives on the phase-chain
    workload; H4's payoff lives on the multi-allele workload; the I/O
    findings (S2) only become visible with a real-file sink. Without
    these, a measurable win from H2 or H4 cannot be confirmed on
    `snp_typical_3_3M` alone.
  - **Pattern matched:** methodology checklist — *Microbenchmarks lie
    about cache behavior; a single-input criterion bench standing in
    for a streaming workload is a Likely finding*.
  - **Mechanism:** N/A — this is bench coverage.
  - **Measurement plan:** add workloads, save baselines as
    `phase_chain_heavy_1M`, `multi_allele_500k`,
    `streaming_bufwriter_file_1M`. Per-workload pass criterion: each
    must clear 1 M records/s (the upstream walker rate) at the lower
    95 % CI for the "writer keeps up" guarantee to hold across the
    three axes.
  - **Complexity cost:** ~100 lines of fixture builders; the bench
    runtime grows from 1 × 30 s to ~4 × 30 s (3 workloads + the
    file-sink one likely sample-size-tuned down).
  - **Suggested experiment / fix:** sketch in
    [methodology.md](../../tmp/perf_review_2026-05-13_psp_writer/methodology.md#L122-L144).

- [benches/psp_writer_perf.rs:138-143](../../benches/psp_writer_perf.rs#L138-L143) — **L9: Wrap bench inputs in `black_box`**
  - **Confidence:** Medium
  - **Hot-path evidence:** the bench wraps the *output* of `write_all`
    in `black_box` but not the inputs (`&records`, `header.clone()`).
    With `lto = "fat"` on `[profile.bench]`, the optimiser sees
    `records` as a loop-invariant slice across `b.iter` calls.
  - **Pattern matched:** methodology checklist — *Use
    `std::hint::black_box`. Constant-folded benchmarks are the most
    common cause of "wins" that disappear in production*.
  - **Mechanism:** harden the bench against loop-invariant
    optimisations that may hide future regressions.
  - **Measurement plan:** apply, rerun, expect no change today — this
    is a correctness-of-measurement fix, not an optimisation.
  - **Complexity cost:** zero.
  - **Suggested experiment / fix:**
    ```rust
    b.iter(|| black_box(write_all(black_box(&records), black_box(header.clone()))));
    ```

- [benches/psp_writer_perf.rs](../../benches/psp_writer_perf.rs) — **L10: Add `write_record`-only and `flush_block`-only sub-benches**
  - **Confidence:** Medium
  - **Hot-path evidence:** the existing bench averages four phases
    (`new + write_record* + flush_block* + finish`). H1 targets only
    the per-record path; H3 targets only the per-flush path. Without
    separable sub-benches, a single-axis fix is diluted in the
    headline number.
  - **Pattern matched:** methodology checklist — *One change per
    measurement*.
  - **Mechanism:** N/A — bench coverage.
  - **Measurement plan:** sketch in §3 item 3. Threshold: each
    sub-bench must produce a stable mean (criterion CV < 5 %).
  - **Complexity cost:** ~40 lines; the `flush_block`-only bench
    needs a `#[doc(hidden)] pub` helper on `PspWriter` to query
    block-fill progress.
  - **Suggested experiment / fix:** add to the existing benchmark
    group; pattern follows
    [benches/pileup_walker_scaling.rs](../../benches/pileup_walker_scaling.rs)'s
    multi-sub-bench layout.

- [src/per_sample_caller/psp/writer.rs:72](../../src/per_sample_caller/psp/writer.rs#L72) — **L11: Add a rustdoc note on `PspWriter::new` about wrapping `File` sinks in `BufWriter`**
  - **Confidence:** Medium
  - **Hot-path evidence:** `grep -rn PspWriter::new` shows the only
    non-test caller is the bench's `io::sink()`. When production
    wiring lands, the writer emits 1 header `write_all` + 12 column
    `write_all`s per flush + 1 index `write_all` + 1 trailer
    `write_all` (32 bytes) per file; against a raw `File` each is one
    `write(2)` syscall. The 32-byte trailer is a particularly bad
    small-write to leave unbuffered.
  - **Pattern matched:** io_and_syscalls checklist — *Wrap in
    `BufWriter::with_capacity(N, file)`*.
  - **Mechanism:** the writer itself should *not* internally wrap its
    sink (then `finish` cannot return a `BufWriter` for the caller to
    `flush()` and `into_inner()` to surface the trailing buffered
    bytes' error). The hook is at the call site; the writer docs
    should advertise the requirement.
  - **Measurement plan:** validated by the
    `streaming_bufwriter_file_1M` workload from L8 + `strace -c -e
    trace=write`. `write` count bounded by `file_bytes / buf_size +
    O(1)`, not by `12 × n_blocks`.
  - **Complexity cost:** one-line rustdoc.
  - **Suggested experiment / fix:**
    ```rust
    /// ...
    /// **Buffering.** This writer issues many small `write_all`
    /// calls per flush (one block header + one per column payload +
    /// one trailer). Wrap a real-file sink in
    /// `BufWriter::with_capacity(64 * 1024, file)` (or larger)
    /// before passing it here. `io::sink()` and `Cursor<Vec<u8>>`
    /// are in-memory and don't need wrapping.
    pub fn new(mut sink: W, header: WriterHeader) -> Result<Self, PspWriteError> { ... }
    ```

### Speculative

- [src/per_sample_caller/psp/writer.rs:59](../../src/per_sample_caller/psp/writer.rs#L59) — **S1: `BTreeSet<SlotId>` running active-slot set → sorted `Vec<u16>` or `SmallVec<[u16; 16]>`**
  - **Confidence:** Low (workload relevance today)
  - **Hot-path evidence:** pattern-match only — the bench workload
    keeps `active_slots` empty throughout, so this set does not
    appear in the profile. Filed at the orchestrator's request for
    phase-chain-heavy workloads.
  - **Pattern matched:** data_layout checklist — *Pointer chases
    fragment the cache*; allocations checklist — *Tiny, often-empty
    collections are SmallVec candidates — sometimes*.
  - **Mechanism:** `BTreeSet<u16>` stores keys in heap-allocated
    B-tree nodes; `contains` / `insert` / `remove` each pay a cache
    miss to reach the root. A sorted `Vec<u16>` of length `n ≤ 32`
    lives in one allocation and fits in one cache line; `contains`
    is `binary_search.is_ok()`. `SmallVec<[u16; 16]>` keeps it
    stack-resident in the common case.
  - **Measurement plan:** depends on the `phase_chain_heavy_1M`
    workload (L8) where ~8 active chains overlap at any locus.
    Threshold: ≥ 5 % improvement on the heavy workload **and** ≤ 1 %
    regression on `snp_typical_3_3M`.
  - **Complexity cost:** sorted `Vec<u16>` requires open-coding
    `insert` / `remove` via `partition_point` (~10 lines + tests);
    `SmallVec` is a new dep. Either preserves the ascending-iteration
    property used by the block snapshot.
  - **Suggested experiment / fix:** sketch in
    [data_layout.md](../../tmp/perf_review_2026-05-13_psp_writer/data_layout.md#L218-L234).

- [src/per_sample_caller/psp/writer.rs:138](../../src/per_sample_caller/psp/writer.rs#L138) — **S2: Pre-size the active-slot snapshot Vec**
  - **Confidence:** Low
  - **Hot-path evidence:** pattern-match only; cold relative to the
    per-record path (once per block, never on this workload).
  - **Pattern matched:** allocations checklist — *Pre-size containers
    when the size is known or bounded*.
  - **Mechanism:** replace `self.active_slots.iter().copied().collect()`
    with `let mut v =
    Vec::with_capacity(self.active_slots.len()); v.extend(self.active_slots.iter().copied()); v`.
    Marginal even on a phase-chain workload.
  - **Measurement plan:** not bench-visible on the current workload.
    Roll into S1 if it lands.
  - **Complexity cost:** none.
  - **Suggested experiment / fix:** as above.

- [src/per_sample_caller/psp/writer.rs:541-594](../../src/per_sample_caller/psp/writer.rs#L541-L594) — **S3: `encode_column` 12-arm `match def.tag` dispatch — confirm jump-table codegen, no action**
  - **Confidence:** Low
  - **Hot-path evidence:** `encode_column` 2.47 % inclusive, not in
    the > 0.5 % SELF table. Runs ~72 times per bench iteration —
    too cold to matter even if the dispatch is suboptimal.
  - **Pattern matched:** hot_loops checklist — *Static dispatch in
    hot loops* (with the caveat that this is not in a hot loop).
  - **Mechanism:** with sparse non-contiguous tag values
    (0x01..0x04, 0x10..0x14, 0x20..0x22) LLVM often emits a
    binary-search chain rather than a jump table. At 72 calls per
    bench × ≤ 1 ns per dispatch ≈ ≤ 100 ns / 480 ms, immeasurable.
  - **Measurement plan:** `cargo asm` on `encode_column` to record
    the codegen shape; revisit if v1.x adds columns and the match
    grows past ~24 arms.
  - **Complexity cost:** none proposed today.
  - **Suggested experiment / fix:** none — record the codegen shape
    for future reference.

- [src/per_sample_caller/psp/writer.rs:273-296](../../src/per_sample_caller/psp/writer.rs#L273-L296) — **S4: `windows(2)` walks in `validate_record` — confirm cost on a chain-heavy workload before rewriting**
  - **Confidence:** Low
  - **Hot-path evidence:** workload shape: empty `chain_slots` /
    `new_chains` / `expired_chains` short-circuits all three loops
    on this bench. Cannot conclude hot from current profile.
  - **Pattern matched:** hot_loops checklist — *Iterators are
    zero-cost most of the time*; don't rewrite without measurement.
  - **Mechanism:** `slots.windows(2)` over a non-empty slice produces
    `(slots[i], slots[i+1])` pairs; the comparison is branchless on
    a register. `<[T]>::is_sorted` would be a drop-in but same
    codegen; not the win to chase.
  - **Measurement plan:** depends on `phase_chain_heavy_1M` (L8). If
    these windows show > 0.5 % self on that workload, upgrade.
  - **Complexity cost:** N/A.
  - **Suggested experiment / fix:** hold until L8 lands.

- [src/per_sample_caller/psp/writer.rs:152-187](../../src/per_sample_caller/psp/writer.rs#L152-L187) — **S5: Document `BufWriter::into_inner` and `File::sync_all` discipline at end-of-stage**
  - **Confidence:** Low
  - **Hot-path evidence:** pattern-match only; correctness-and-future-API
    trap, not a speed win.
  - **Pattern matched:** io_and_syscalls checklist — *`BufWriter::drop`
    may swallow flush errors. For correctness-critical writes, call
    `flush()` explicitly and propagate the error*; *`fsync` is not
    free … at end-of-stage, not per-record*.
  - **Mechanism:** when the production caller wires `BufWriter<File>`
    into `PspWriter::new`, `PspWriter::finish` returns the
    `BufWriter`; if the caller drops it instead of
    `flush()`+`into_inner()`, any data still in the buffer flushes in
    `BufWriter::drop` and any error is swallowed. For billions-of-records
    files this can silently truncate the trailer.
  - **Measurement plan:** N/A.
  - **Complexity cost:** trivial: end-of-stage:
    ```rust
    let buf = writer.finish()?;          // PSP-level flush
    let file = buf.into_inner()?;        // surface BufWriter errors
    file.sync_all()?;                    // durability
    ```
  - **Suggested experiment / fix:** add the snippet as a rustdoc
    example on `PspWriter::finish`.

### Note

- [src/per_sample_caller/psp/writer.rs:198-244](../../src/per_sample_caller/psp/writer.rs#L198-L244) — the chrom-id bounds check at line 205 immediately precedes the index `&self.header.chromosomes[record.chrom_id as usize]` at line 212, so LLVM already elides the second check. Confirmed by reading; nothing to fix.
- [src/per_sample_caller/psp/block.rs:65-91](../../src/per_sample_caller/psp/block.rs#L65-L91) — `WireScalar::encode_le` per-element `extend_from_slice` of a 2/4/8-byte slice has the same slab opportunity as L6 (`encode_list_column`). Below the 0.5 % SELF reporting bar on this workload; file only as part of L6's diff to keep codecs symmetric.
- [src/per_sample_caller/psp/writer.rs:243,253,262,269,277,292](../../src/per_sample_caller/psp/writer.rs#L243) — error-path `format!()` allocations in `validate_record`. Cold by construction (only fires on producer bugs). Listed only because L5 marks them `#[cold]`.
- [src/per_sample_caller/psp/registry.rs](../../src/per_sample_caller/psp/registry.rs) — `ColumnDef`'s `name`, `description`, and `length_column` are all `&'static str`, so registry walks on `flush_block` do not allocate. Good baseline.
- [src/per_sample_caller/psp/writer.rs:437-461](../../src/per_sample_caller/psp/writer.rs#L437-L461), [src/per_sample_caller/psp/index.rs:30-36](../../src/per_sample_caller/psp/index.rs#L30-L36) — `BlockAccumulator` and `BlockIndexEntry` are both `repr(Rust)`; the compiler reorders for minimum padding. No layout fix.
- [src/per_sample_caller/psp/writer.rs:55](../../src/per_sample_caller/psp/writer.rs#L55), [src/per_sample_caller/psp/writer.rs:349](../../src/per_sample_caller/psp/writer.rs#L349) — `Option<BlockAccumulator>` is niche-optimised (compiler stuffs the `None` discriminant into the inner `Vec`'s nullable pointer). No size or indirection cost; the per-record branch is predictable.
- [src/per_sample_caller/psp/writer.rs:160-181](../../src/per_sample_caller/psp/writer.rs#L160-L181) — `encode_index` + `encode_trailer` single-shot writes are fine (one syscall each on a buffered sink, zero against `io::sink()`).
- `_evidence.md` workload shape — `XXH64_round` and `xxhash` self-samples (0.94 % combined) confirm the per-column zstd content checksum is on the measured path. Any change to L1's buffer reuse must preserve the `include_checksum(true)` semantics; H3's regression test verifies this.

## 6. Out-of-scope observations

- [benches/psp_writer_perf.rs](../../benches/psp_writer_perf.rs) — the bench was added in this session and is intentionally minimal. L7–L10 trace the gaps; landing them is the next-PR menu before any §5 H/L fix is merged.
- The reader is not yet implemented. Several findings here (`BlockAccumulator` capacity heuristics, scratch-buffer reuse pattern, `zstd::bulk::Compressor` shape) have natural mirrors on the reader side (`zstd::bulk::Decompressor` reuse, capacity-pre-sized column buffers). When the reader lands, this same review should be repeated against it.
- The eventual CLI / per-sample pipeline wiring — when written — must wrap `File` sinks in `BufWriter::with_capacity(64 KiB, file)` and end-of-stage with the `into_inner()? → sync_all()?` discipline (S5). Filed as a follow-up issue, not blocking this review's findings.
- `samply` could not run inside the rootless container (CAP_PERFMON unavailable) — the profile this review cites was captured by running the container-built bench binary on the host with host samply (`perf_event_paranoid=1`). Documented for the next reviewer.

## 7. What's already good

- **SoA scalar columns.** `BlockAccumulator` already keeps each
  fixed-width column as its own `Vec<T>` ([writer.rs:445-457](../../src/per_sample_caller/psp/writer.rs#L445-L457)), exactly the right shape for the slab-encode and per-column zstd pipeline — the remaining work (H1, H2) is just per-buffer capacity and the ragged-2D collapse, not a re-architecture.
- **`ColumnDef` is fully `&'static`.** [registry.rs](../../src/per_sample_caller/psp/registry.rs) holds `name`, `description`, and `length_column` as `&'static str`; the writer's per-flush walk of `V1_0_COLUMNS` is allocation-free.
- **Build configuration is tight.** [Cargo.toml:6-13](../../Cargo.toml#L6-L13) has `lto = "fat"`, `codegen-units = 1`, and `debug = "line-tables-only"` on `[profile.release]`, with `[profile.bench]` inheriting and adding `debug = true`. [.cargo/config.toml](../../.cargo/config.toml) sets `target-cpu = x86-64-v3` on Linux/x86_64 — AVX2 + FMA + BMI2 floor, which matches the profile's `__memcpy_avx_unaligned_erms` / `__memset_avx2_unaligned_erms` use.

### Author response convention

Address each finding by its identifier (e.g., "H1", "L2") with one
of: `applied in <commit>` / `experiment shows no gain — closing` /
`disputed because …` / `deferred to <issue>` / `won't fix because …`.
The "experiment shows no gain" path is expected and welcome — that
is what the measurement plan is for.
