# Performance Review: psp_reader
**Date:** 2026-05-13
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the `.psp` per-sample-pileup reader module
**Verdict:** Run experiments
**Hot-path evidence:** criterion bench only (`benches/psp_reader_perf.rs`, added this session); **sampling profile NOT available** in this environment — see §1 *Constraints*.

---

## 1. Scope and constraints

- **What was reviewed:** the on-disk `.psp` reader at
  [src/per_sample_caller/psp/](../../src/per_sample_caller/psp/),
  with [reader.rs](../../src/per_sample_caller/psp/reader.rs) as the
  focus and the supporting block / varint / index / trailer / registry
  decoder primitives in scope wherever the reader hot path touches
  them.
- **Reviewed against:** commit `d4dc55f` on branch `main` (working
  tree: the new bench at
  [benches/psp_reader_perf.rs](../../benches/psp_reader_perf.rs) added
  for this review; one new `[[bench]]` entry in
  [Cargo.toml](../../Cargo.toml)).
- **Throughput / latency targets, input sizes, target hardware:** the
  reader is Stage 2's source. Production input: ~3 × 10⁹ per-position
  pileup records per WGS sample at 30× coverage. Target throughput:
  ≥ 1 M records/s sustained on the consumer side (matches the writer
  target — symmetric pipeline). Today's bench shows 11.4 / 5.9 /
  5.0 Mrec/s on snp / phase / multi streaming workloads and
  1.27 Mrec/s on a 100 k-record region query — comfortably over the
  floor on streaming, but the region path is single-block-bound.
  Target hardware: developer-class Linux/x86_64, AVX2 floor configured
  via `.cargo/config.toml`'s `target-cpu=x86-64-v3`.
- **Hot-path evidence available:**
  - **Sampling profile: NOT REACHABLE.** `samply record` returns
    `Failed to start profiling: Operation not permitted (os error 1)`
    in this sandboxed environment even when run as root with
    `kernel.perf_event_paranoid = 1`. Per skill guidance, **all
    code-level findings in §5 are capped at Likely** — no Hot-path is
    reachable without profile evidence.
  - Criterion baseline (verbatim, this session, host build of bench
    binary):
    ```
    psp_reader/snp_typical_3_3M
                            time:   [284.63 ms 289.83 ms 291.92 ms]
                            thrpt:  [11.304 Melem/s 11.386 Melem/s 11.594 Melem/s]
    psp_reader/phase_chain_heavy_1M
                            time:   [165.82 ms 170.58 ms 174.60 ms]
                            thrpt:  [5.7274 Melem/s 5.8623 Melem/s 6.0308 Melem/s]
    psp_reader/multi_allele_500k
                            time:   [96.910 ms 100.04 ms 101.51 ms]
                            thrpt:  [4.9258 Melem/s 4.9982 Melem/s 5.1594 Melem/s]
    psp_reader_region/region_window_chr1_mid_100k
                            time:   [75.006 ms 78.817 ms 80.505 ms]
                            thrpt:  [1.2422 Melem/s 1.2688 Melem/s 1.3332 Melem/s]
    ```
    Saved as criterion baseline `reader_baseline` /
    `reader_region_baseline`. Re-run with
    `cargo bench --bench psp_reader_perf -- --baseline reader_baseline`
    after each candidate fix to gate.
  - **Indirect hot-path evidence via the psp::writer review** at
    [ia/reviews/perf_psp_writer_2026-05-13.md](perf_psp_writer_2026-05-13.md).
    The reader mirrors the writer's column structure on the decode
    side. The writer review *did* have a real samply profile, and its
    H1–H4 findings produced measured wins (commits 26f3d65, 47c1514,
    969de6c, 154c4eb). Findings here that name a reader analog of one
    of those sites cite the writer's profile as indirect evidence; the
    finding remains capped at **Likely** because the reader's actual
    hot path is unverified.
  - Workload-vs-writer analogy table:
    [tmp/perf_review_2026-05-13_psp_reader/_evidence.md](../../tmp/perf_review_2026-05-13_psp_reader/_evidence.md).
- **In-scope files:**
  - [src/per_sample_caller/psp/reader.rs](../../src/per_sample_caller/psp/reader.rs)
  - [src/per_sample_caller/psp/block.rs](../../src/per_sample_caller/psp/block.rs)
  - [src/per_sample_caller/psp/varint.rs](../../src/per_sample_caller/psp/varint.rs)
  - [src/per_sample_caller/psp/index.rs](../../src/per_sample_caller/psp/index.rs)
  - [src/per_sample_caller/psp/trailer.rs](../../src/per_sample_caller/psp/trailer.rs)
  - [src/per_sample_caller/psp/registry.rs](../../src/per_sample_caller/psp/registry.rs)
  - [benches/psp_reader_perf.rs](../../benches/psp_reader_perf.rs)
  - [Cargo.toml](../../Cargo.toml)
- **Deliberately out of scope:**
  - The writer at
    [src/per_sample_caller/psp/writer.rs](../../src/per_sample_caller/psp/writer.rs)
    — already reviewed.
  - The TOML header parser/serializer in
    [src/per_sample_caller/psp/header.rs](../../src/per_sample_caller/psp/header.rs)
    — once-per-file, cold.
  - Test code under `#[cfg(test)]`.
  - Pileup walker, BAQ engine, CRAM input — separate prior reviews.
  - The eventual `PileupRecord` / `AlleleObservation` API redesign to
    expose borrowed views — large enough to deserve its own design
    review. Filed under Speculative for cross-reference.
- **Categories dispatched:**
  - `methodology` — always; reader has no bench, no dhat example, no
    sampling profile available.
  - `allocations` — the writer review showed allocator pressure was
    35 % inclusive; the reader has the symmetric structure on every
    column.
  - `data_layout` — `DecodedBlock` carries the same `Vec<Vec<SlotId>>`
    ragged 2-D the writer's H2 attacked.
  - `hot_loops` — varint decode, scalar decode, list decode,
    per-record validation loops.
  - `io_and_syscalls` — `zstd::stream::decode_all` per column,
    `read_exact` per column, `read_block_header` grow loop, BufReader
    sizing on the production caller side.
  - `concurrency` — **skipped.** The reader is fully single-threaded;
    no `Arc`, `Mutex`, atomics, channels, `rayon`, or `async fn`. The
    pipeline as a whole is parallel at the per-sample level (one
    reader per worker), but each `PspReader` is owned by one thread.

## 2. Verdict

**Run experiments.** No Hot-path findings are reachable without a
sampling profile, so this review's primary deliverable is *gated
candidates*: each L-finding below carries a measurement plan against
the criterion baseline saved this session, and a complexity cost. The
top tier (L1–L4) is unusually well-supported because each one has a
direct structural analog on the writer side that produced a measured
win — applying the same fix on the symmetric site has a strong
plausibility argument plus a ready-to-run bench gate. The team should
land L1–L4 as separate small PRs, each gated by `cargo bench --bench
psp_reader_perf -- --baseline reader_baseline`.

The methodology-side gaps (M1–M5, mostly tooling) should land before
or alongside the code-level fixes — the **dhat-heap example** (M3)
is the primary instrument for L1/L2/L3/L4 in particular, and a real
sampling profile (M1's escape-hatch documentation) would upgrade
several of the Likely findings to Hot-path.

## 3. Measurement plan

Run before any code-level work; the listed wins below cite the same
benches.

1. **Use `benches/psp_reader_perf.rs` (added this session) as the
   gate.** Save the current numbers as criterion baseline
   `reader_baseline` (already done). After each candidate fix:
   `cargo bench --bench psp_reader_perf -- --baseline reader_baseline`.
   Threshold per finding stated below.

2. **Add `examples/dhat_psp_reader.rs`** mirroring
   [examples/dhat_psp_writer.rs](../../examples/dhat_psp_writer.rs) —
   the `dhat-heap` feature is already wired in
   [Cargo.toml](../../Cargo.toml#L72). Run inside the dev container
   (or on the host if the container is unavailable):
   `cargo run --release --example dhat_psp_reader --features dhat-heap`.
   Save the baseline `dhat-heap.json` under
   `tmp/dhat_psp_reader_baseline.json`. **This is the primary
   instrument for L1, L2, L3 — wall-clock alone will not show whether
   a per-column allocation site was eliminated.**

3. **Document a sampling-profile escape hatch.** The reader review
   could not collect a `samply` / `perf record` profile in this
   sandbox. Add `scripts/profile.sh` (or extend `scripts/dev.sh`) so
   the next reviewer has a documented recipe for either (a) running
   on a privileged host where `perf_event_paranoid <= 1` is set, or
   (b) extending the dev container with `--cap-add=PERFMON
   --security-opt seccomp=unconfined`. Until this lands every reader
   review will be capped at Likely.

4. **Add bench coverage gaps** flagged under M2 / M5:
   - `decode_block_only` sub-bench: pre-decompress one block, time
     `decode_block_payload` over a `Cursor<&[u8]>`. Isolates the
     per-column codec cost from the per-block I/O / open cost. Needs
     a `#[doc(hidden)] pub` re-export of `decode_block_payload`.
   - `region_open_only` + `region_records_only` split: the current
     `region_window_chr1_mid_100k` measures
     `PspReader::new` + iterate inside one `b.iter`; split so the
     binary-search + first-block-decode cost is separable from the
     in-window per-record cost.
   - `bufreader_file_1M_64KiB` + `bufreader_file_1M_8KiB` +
     `bufreader_file_1M_1MiB` workloads, mirroring the writer-side
     `bufwriter_file_1M_64KiB` from
     [benches/psp_writer_perf.rs](../../benches/psp_writer_perf.rs).
     Production callers will hit a `File`, not a `Cursor`; the I/O
     findings (L9, L10) only become measurable once the bench has a
     file-backed source.
   - Add `Throughput::Bytes` alongside `Throughput::Elements` so the
     records/byte ratio is comparable to disk-bandwidth bounds.
   - Add per-iteration `assert_eq!(n, NUM_RECORDS_X as u64)` so a
     future iterator-truncation bug is caught at bench time.

5. **Run an allocator A/B** once L1 + L2 have landed. The
   `alloc-mimalloc` feature is already wired in
   [Cargo.toml:79](../../Cargo.toml#L79) and the bench file has the
   `#[global_allocator]` declaration at lines 30-33. Threshold to
   merge a default change: ≥ 5 % wall-clock improvement on
   `snp_typical_3_3M`. Don't change the default if the win is workload-
   dependent.

## 4. Build / toolchain configuration

- **Add `[profile.profiling]` to
  [Cargo.toml](../../Cargo.toml).** `[profile.release]` uses
  `debug = "line-tables-only"` (good for production binary size); the
  bench profile correctly upgrades to `debug = true`. The gap is
  on the *production-binary* side: when the eventual `psp` CLI / merge
  driver is profiled in-situ on a customer dataset, line-tables-only
  gives callstacks but no inlined-frame resolution, which makes
  writer-style "self %" attribution hard to read. Trivial:
  ```toml
  [profile.profiling]
  inherits = "release"
  debug = true
  ```
  No effect on the production build (it stays `--release`); use
  `cargo build --profile profiling --bin <future-reader-cli>` when an
  in-situ profile is needed.
- **Already correct, no action needed (called out so they are not
  re-litigated):**
  - `[profile.release]`: `lto = "fat"`, `codegen-units = 1`,
    `debug = "line-tables-only"`, `panic = "abort"`.
  - `[profile.bench] inherits = "release"; debug = true`.
  - [.cargo/config.toml](../../.cargo/config.toml): `target-cpu =
    x86-64-v3` on Linux/x86_64 (AVX2 + FMA + BMI2 floor) and
    `apple-m1` on aarch64-macos.
  - [rust-toolchain.toml](../../rust-toolchain.toml) pinned at 1.95
    for bench reproducibility.
  - `alloc-mimalloc` feature wired via
    [Cargo.toml:79](../../Cargo.toml#L79); the new
    [benches/psp_reader_perf.rs](../../benches/psp_reader_perf.rs#L29-L33)
    declares the `#[global_allocator]` correctly.
- **PGO** is on the table for the reader once a CLI binary lands with
  a stable workload (per-sample WGS read end-to-end), but defer until
  L1–L4 are exhausted; PGO adds non-trivial build-pipeline complexity
  for ≤ 5–10 % typically. Same defer-call as the writer review.

## 5. Code-level findings

### Likely

- [src/per_sample_caller/psp/block.rs:531-533](../../src/per_sample_caller/psp/block.rs#L531-L533) + [src/per_sample_caller/psp/reader.rs:1184](../../src/per_sample_caller/psp/reader.rs#L1184) — **L1: Reuse one `zstd::bulk::Decompressor` across all columns and blocks (mirror of writer's H3)**
  - **Confidence:** High
  - **Hot-path evidence:** Indirect — the writer's H3 quoted samply
    inclusive shares `zstd::compress 17.43 % + cwksp 5.68 % +
    compressBegin 7.12 % + __brk 4.70 % + morecore 4.70 % + systrim
    4.62 %` for the symmetric site. The reader's
    `zstd::stream::decode_all` builds a fresh `DCtx` workspace per
    `decode_one_column` call, called 11 times per block × N blocks.
    The DCtx workspace is large enough at level-9 frames to push glibc
    past fastbins, so each call pays one `sbrk` up + one `systrim`
    down. The writer-side fix (commit 969de6c) drove `__brk` self
    from 4.70 % to noise.
  - **Pattern matched:** allocations checklist — *Allocations belong
    outside hot loops*; io_and_syscalls checklist — *Per-record /
    per-column decompression context creation is a hidden allocation*.
  - **Mechanism:** allocate one `zstd::bulk::Decompressor<'static>` on
    `RecordsIter` (or `PspReader` if shared across iterators); call
    `decompress_to_buffer(input, &mut out)` per column. The DCtx
    workspace is reused; the per-call cost drops to the actual
    decompression work plus a frame-header reset. Pairs with L2 (the
    same refactor wires both — `decompress_to_buffer` writes into a
    borrowed `Vec<u8>`).
  - **Measurement plan:** primary instrument is the dhat example
    (M3); threshold: `total_blocks` from `dhat-heap.json` drops by
    ≥ 30 % on `snp_typical_3_3M`, ≥ 20 % on `phase_chain_heavy_1M`.
    Secondary instrument: `cargo bench --bench psp_reader_perf --
    --baseline reader_baseline`; threshold ≥ 5 % wall-clock improvement
    on `snp_typical_3_3M`. Expect a larger lift on the region bench
    because per-block fixed cost dominates there. After both pass,
    rerun with `samply` (once M1 lands) and confirm `__brk` self < 1 %.
  - **Complexity cost:** one `Decompressor<'static>` field on
    `RecordsIter`, one reusable `Vec<u8>` for the decompressed bytes
    (paired with L2). `decode_one_column` gains `&mut Decompressor`
    + `&mut Vec<u8>` arguments — mechanical signature ripple. No
    `unsafe`, no new dependency. One regression test that decodes a
    known frame twice through the reused context and asserts the
    bytes match the single-shot path.
  - **Suggested experiment / fix:**
    ```rust
    // block.rs — symmetric to new_column_compressor() / zstd_compress_into()
    pub fn new_column_decompressor() -> io::Result<zstd::bulk::Decompressor<'static>> {
        zstd::bulk::Decompressor::new()
    }

    pub fn zstd_decompress_into(
        decompressor: &mut zstd::bulk::Decompressor<'static>,
        input: &[u8],
        uncompressed_len_hint: usize,
        out: &mut Vec<u8>,
    ) -> io::Result<()> {
        out.clear();
        out.reserve(uncompressed_len_hint);
        out.resize(uncompressed_len_hint, 0);
        let n = decompressor.decompress_to_buffer(input, out).map_err(io::Error::other)?;
        out.truncate(n);
        Ok(())
    }
    ```
    Hold one `Decompressor` and one reusable `decompressed_scratch:
    Vec<u8>` on `RecordsIter`; thread them through
    `decode_block_payload` → `decode_one_column`. Land together with
    L2 in one PR.

- [src/per_sample_caller/psp/reader.rs:1180,1143](../../src/per_sample_caller/psp/reader.rs#L1180) + [reader.rs:874](../../src/per_sample_caller/psp/reader.rs#L874) — **L2: Hoist per-block scratch buffers (`compressed`, `decompressed`, `block_header_buf`) onto `RecordsIter` (mirror of writer's L1)**
  - **Confidence:** High
  - **Hot-path evidence:** pattern-match plus the writer-side analog
    (L1, applied with H3 in commit 969de6c). Three sites:
    - `let mut compressed = vec![0u8; entry.compressed_len as usize];`
      ([reader.rs:1180](../../src/per_sample_caller/psp/reader.rs#L1180)) —
      every known column.
    - `let mut sink = vec![0u8; entry.compressed_len as usize];`
      ([reader.rs:1143](../../src/per_sample_caller/psp/reader.rs#L1143)) —
      every unknown optional column (cold today, future-extension path).
    - `Vec::with_capacity(BLOCK_HEADER_INITIAL_CHUNK)` (4 KiB) per
      `read_block_header` call ([reader.rs:874](../../src/per_sample_caller/psp/reader.rs#L874)) — once per block.
    Plus the `Vec<u8>` returned by `zstd_decompress` paired with L1.
    Per-block: ~22 short-lived `Vec<u8>` allocations.
  - **Pattern matched:** allocations checklist — *Allocations belong
    outside hot loops*.
  - **Mechanism:** three reusable `Vec<u8>` fields on `RecordsIter`
    (`compressed_scratch`, `decompressed_scratch`,
    `block_header_buf`). `Vec::clear()` between uses keeps capacity;
    the buffers converge to the largest column / header the reader
    has seen. Combined with L1: same allocation that paired the
    writer's L1 + H3.
  - **Measurement plan:** dhat (M3) is the primary instrument;
    threshold ≥ 30 % drop in `total_blocks` on `snp_typical_3_3M`.
    Wall-clock change likely small in isolation; bundled with L1 it
    will compound visibly.
  - **Complexity cost:** three new fields on `RecordsIter`; the
    accompanying `decode_one_column` and `read_block_header` signature
    changes are mechanical. No `unsafe`, no new dependency.
  - **Suggested experiment / fix:** sketch in
    [allocations.md](../../tmp/perf_review_2026-05-13_psp_reader/allocations.md).
    Land in the same PR as L1.

- [src/per_sample_caller/psp/reader.rs:417-434](../../src/per_sample_caller/psp/reader.rs#L417-L434) + [block.rs:320-378](../../src/per_sample_caller/psp/block.rs#L320-L378) — **L3: Collapse the three `Vec<Vec<SlotId>>` columns in `DecodedBlock` to flat CSR via a sibling `decode_list_column_csr` (mirror of writer's H2)**
  - **Confidence:** Medium
  - **Hot-path evidence:** Indirect — writer's H2 (commit 26f3d65)
    measured a 35.39 % realloc share for the symmetric encode-side
    `Vec<Vec<SlotId>>` ragged 2-D, and adopting CSR plus the encoder's
    `encode_list_column_csr` (block.rs:297-316) closed it. The reader
    rebuilds the same shape: `decode_list_column` allocates one
    `Vec<SlotId>` per row × three columns. Bench evidence:
    `phase_chain_heavy_1M` runs at 5.86 Mrec/s vs `snp_typical_3_3M`
    at 11.39 Mrec/s — the chain-slot-decoding workload is materially
    slower. Direct-cause attribution awaits a profile.
  - **Pattern matched:** data_layout checklist — *ragged `Vec<Vec<T>>`
    is the worst of both worlds (outer SoA spine + inner per-row heap
    allocation)*; allocations checklist — *Allocations belong outside
    hot loops*.
  - **Mechanism:** add a `decode_list_column_csr<T: WireScalar +
    bytemuck::Pod>(bytes, expected_count, name) -> CsrColumn<T>`
    primitive in `block.rs` (sibling of the writer's
    `encode_list_column_csr`). Replace `Vec<Vec<SlotId>>` in
    `DecodedBlock` with `CsrColumn<SlotId>` for each of the three
    list columns. The materialiser at
    [reader.rs:683-714,742-756,762-783](../../src/per_sample_caller/psp/reader.rs#L683)
    moves from `mem::take(&mut block.allele_chain_slots[j])` to
    `block.allele_chain_slots.row(j).to_vec()`. **Caveat:** today's
    `mem::take` hands ownership of the inner `Vec` to
    `PileupRecord` / `AlleleObservation` (which own their `Vec<SlotId>`
    by API contract). With CSR, the materialiser pays one allocation
    per row at emit time (row → owned `Vec<SlotId>`), in exchange for
    *eliminating* the per-row allocation at decode time. The trade is
    "many tiny allocations on decode, zero on emit (mem::take)" vs
    "two big allocations on decode, one tiny per row on emit
    (to_vec)". The allocator is materially faster on the small,
    fixed-size emit allocations than on the heterogeneous `Vec::with_capacity(k)`
    calls today, and the writer-side analog confirmed the direction.
    Net direction is toward CSR, but the emit-side allocation has to
    be measured — this is not a guaranteed wall-clock win.
  - **Measurement plan:** primary on `phase_chain_heavy_1M` from
    `reader_baseline` — threshold ≥ 5 % wall-clock improvement
    **and** SNP workload not slower at criterion's 95 % CI. dhat
    threshold: per-row inner-Vec allocations drop to zero on
    phase-heavy. If the wall-clock win evaporates because the emit-
    time allocations cancel it, the secondary fallback is to keep
    the public API's `Vec<SlotId>` ownership but pre-allocate all
    emit-time row Vecs in one bulk pass right after decode (one
    allocator hit per row, but better allocator locality than the
    interleaved per-row + per-column today).
  - **Complexity cost:** one new `CsrColumn<T>` type (~30 lines of
    helper struct in block.rs); three list-column field replacements
    in `DecodedBlock`; one new `decode_list_column_csr` primitive;
    materialiser-side rewrite of the three `mem::take` sites. No
    `unsafe`, no new dependency (`bytemuck` already in tree). The
    wire format is unchanged.
  - **Suggested experiment / fix:** sketch in
    [data_layout.md](../../tmp/perf_review_2026-05-13_psp_reader/data_layout.md).
    Land *after* L1+L2 so the bench gate is not muddled by other
    allocator-pressure changes.

- [src/per_sample_caller/psp/block.rs:417-470](../../src/per_sample_caller/psp/block.rs#L417-L470) — **L4: `decode_bytes_split` per-allele `to_vec()` — special-case `n == 0` and `n == 1`, or push the bytes column to CSR shape too**
  - **Confidence:** Medium
  - **Hot-path evidence:** pattern-match plus the writer-side L1
    analog. The `multi_allele_500k` bench is the workload that
    exercises this — 500 k records × 2-4 alleles × one
    `to_vec()` per allele = 1–2 M small `Vec<u8>` allocations per
    file. Bench landed at 5.0 Mrec/s, the slowest of the three
    streaming workloads, consistent with allele-bytes pressure.
    Direct attribution awaits a profile.
  - **Pattern matched:** allocations checklist — *Allocations belong
    outside hot loops*; hot_loops checklist — short-row inner-loop
    overhead dominates the actual memcpy.
  - **Mechanism:** two-stage fix.
    - **Stage A (in-category, minimal):** special-case `n == 0` →
      `Vec::new()` (no allocation). The SNP-typical workload has
      `n_alleles == 1` always, but multi-allele has lots of `n == 0`-
      length alleles only when REF is encoded as empty (rare today).
      Modest win.
    - **Stage B (cross-category, with L3):** decode the bytes column
      into a CSR pair `(data: Vec<u8>, offsets: Vec<u32>)` instead of
      `Vec<Vec<u8>>`; the materialiser still pays one allocation per
      allele at emit time (because `AlleleObservation::seq` is
      `Vec<u8>`), same as the L3 trade.
  - **Measurement plan:** dhat (M3) on `multi_allele_500k`;
    threshold: ≥ 15 % `total_blocks` drop. Wall-clock on
    `multi_allele_500k`; threshold ≥ 5 %. If Stage A alone produces
    no measurable change, defer Stage B until L3 has landed and the
    same-shape decoder primitive can be reused.
  - **Complexity cost:** Stage A is ~5 lines. Stage B is one new
    `decode_bytes_concat(bytes, lengths, name, max_entry_len) ->
    Result<(Vec<u8>, Vec<u32>), PspReadError>` primitive plus a
    `DecodedBlock::allele_seqs` shape change. Same shape as L3.
  - **Suggested experiment / fix:** see Stage A and Stage B sketches
    in [allocations.md](../../tmp/perf_review_2026-05-13_psp_reader/allocations.md)
    and [hot_loops.md](../../tmp/perf_review_2026-05-13_psp_reader/hot_loops.md).

- [src/per_sample_caller/psp/varint.rs:73-92](../../src/per_sample_caller/psp/varint.rs#L73-L92) — **L5: Split `decode_u64_leb128` into an inlined `< 0x80` fast path and a `#[cold] #[inline(never)]` multi-byte slow path (mirror of writer's L4)**
  - **Confidence:** Medium
  - **Hot-path evidence:** writer-side `encode_u64_leb128` measured
    1.22 % self in the writer profile and got the same fast/cold
    split (commit f95c6e5). Decoder is the symmetric site:
    `decode_varint_column` runs one varint per record per column
    (delta-pos, n-alleles), and `decode_list_column` runs one varint
    per row of three list columns. On the SNP-typical workload most
    decoded values are single-byte (`delta_pos = 1`, `n_alleles = 1`,
    list-counts = 0).
  - **Pattern matched:** hot_loops checklist — *Cold paths get
    `#[cold]`*; fast-path / slow-path split.
  - **Mechanism:** today the function loops with a per-iteration
    bounds-check + shift + OR + continuation-bit branch even on
    1-byte values. An inlined fast path returns in 3 instructions
    (load, compare, return); the multi-byte path moves to a
    `#[cold] #[inline(never)]` helper so the call-site icache stays
    tight.
  - **Measurement plan:** `cargo bench --bench psp_reader_perf --
    --baseline reader_baseline`. Threshold ≥ 3 % wall-time on
    `snp_typical_3_3M` or `phase_chain_heavy_1M`. Pair with
    `cargo asm` on the rewritten function to confirm the fast-path
    inlines.
  - **Complexity cost:** one extra fn (~25 lines, the cold path).
    No new types, no `unsafe`.
  - **Suggested experiment / fix:**
    ```rust
    #[inline]
    pub fn decode_u64_leb128(bytes: &[u8]) -> Result<(u64, usize), VarintError> {
        match bytes.first() {
            Some(&b) if b < 0x80 => Ok((b as u64, 1)),
            Some(_) => decode_u64_leb128_cold(bytes),
            None => Err(VarintError::Truncated),
        }
    }

    #[cold]
    #[inline(never)]
    fn decode_u64_leb128_cold(bytes: &[u8]) -> Result<(u64, usize), VarintError> {
        let mut value: u64 = 0;
        let mut shift: u32 = 0;
        for (i, &b) in bytes.iter().enumerate().take(MAX_VARINT_BYTES) {
            let data = u64::from(b & 0x7f);
            value |= data << shift;
            if b & 0x80 == 0 {
                return Ok((value, i + 1));
            }
            shift += 7;
        }
        if bytes.len() >= MAX_VARINT_BYTES {
            Err(VarintError::Overflow)
        } else {
            Err(VarintError::Truncated)
        }
    }
    ```

- [src/per_sample_caller/psp/block.rs:178-209](../../src/per_sample_caller/psp/block.rs#L178-L209) — **L6: Slab-cast `decode_scalar_column` on little-endian for `T: Pod` (mirror of writer's L6)**
  - **Confidence:** Medium
  - **Hot-path evidence:** writer-side L6 (commit 207b3e9) measured
    1.14 % self on the symmetric encode-side `encode_list_column`.
    Reader's `decode_scalar_column` is the mirror: one per fixed-width
    scalar column per block. v1.0 has five per-allele fixed-width
    scalar columns (`allele-obs-count` u32, `allele-q-sum-log` f64,
    `allele-fwd-count` u32, `allele-placed-left-count` u32,
    `allele-placed-start-count` u32) — five `decode_scalar_column`
    calls per block × `n_total_alleles` element decodes per call.
    `multi_allele_500k` is the per-allele-heavy workload (5.0 Mrec/s,
    the slowest streaming bench).
  - **Pattern matched:** hot_loops checklist — *Bounds checks: avoid
    by structure*; *Autovectorization needs the compiler's
    confidence*.
  - **Mechanism:** today's loop calls `T::decode_le(&bytes[cursor..])`
    per element (which does `bytes[..$width].try_into().unwrap()` +
    `<$t>::from_le_bytes(arr)` + per-iter bounds check + `Vec::push`).
    On `target_endian = "little"` plus `T: Pod`, the entire payload
    is `bytemuck::cast_slice::<u8, T>(bytes).to_vec()` — one length
    check, one allocation, one memcpy. The added `Pod` bound covers
    every v1.0 fixed-width scalar except `bool` (which keeps the
    per-element path because `0/1`-validation is not a memcpy).
  - **Measurement plan:** `cargo bench --bench psp_reader_perf --
    --baseline reader_baseline` on `multi_allele_500k`. Threshold
    ≥ 4 % wall-clock improvement; no regression on
    `snp_typical_3_3M`.
  - **Complexity cost:** one new `decode_scalar_column_pod<T:
    WireScalar + bytemuck::Pod>` primitive in `block.rs`; per-arm
    dispatch in `decode_one_column` to use it for the five `Pod` arms.
    No `unsafe` (bytemuck enforces it via trait). The existing
    `decode_scalar_column` keeps the `bool` callers.
  - **Suggested experiment / fix:** sketch in
    [hot_loops.md](../../tmp/perf_review_2026-05-13_psp_reader/hot_loops.md).

- [src/per_sample_caller/psp/reader.rs:1116-1276](../../src/per_sample_caller/psp/reader.rs#L1116-L1276) + [reader.rs:658-790](../../src/per_sample_caller/psp/reader.rs#L658-L790) — **L7: Mark error-construction sites `#[cold]` (mirror of writer's L5)**
  - **Confidence:** Medium
  - **Hot-path evidence:** symmetric with the writer's L5 (applied
    in commit 154c4eb). Every `return Err(PspReadError::...)` arm in
    `decode_one_column` and `materialise_next_record` builds the
    error inline. Without `#[cold]` LLVM may interleave the cold
    construction code with the hot decode loop, hurting branch
    layout and icache. Confirmed by the hot_loops sub-agent that
    `reader.rs` / `block.rs` / `index.rs` carry zero `#[cold]`
    markers today.
  - **Pattern matched:** hot_loops checklist — *Cold paths get
    `#[cold]`*.
  - **Mechanism:** extract each repeated error arm into a `#[cold]
    #[inline(never)] fn err_X(...) -> PspReadError` helper (or wrap
    in `#[cold]` closures for one-off sites). About 12–15 helpers
    across the two functions.
  - **Measurement plan:** `cargo bench --bench psp_reader_perf --
    --baseline reader_baseline`. Threshold ≥ 1 % wall-time on any
    workload (writer's L5 was modest in cycles but tightened the
    icache pattern). Optional: `cargo asm` to confirm the hot block
    is contiguous. Composes naturally with L5 (same files).
  - **Complexity cost:** mechanical refactor; no semantic change.
  - **Suggested experiment / fix:**
    ```rust
    #[cold] #[inline(never)]
    fn err_column_truncated(column: &str, decoded: usize, expected: usize) -> PspReadError {
        PspReadError::ColumnTruncated {
            column: column.to_string(),
            decoded,
            expected,
        }
    }
    // call sites: `return Err(err_column_truncated(name, 0, n));`
    ```

- [src/per_sample_caller/psp/reader.rs:108-121,343-356](../../src/per_sample_caller/psp/reader.rs#L108-L121) — **L8: Update the `PspReader::new` and `region_records` rustdoc examples to recommend `BufReader::with_capacity(64 * 1024, file)` (mirror of writer's L11)**
  - **Confidence:** Medium
  - **Hot-path evidence:** pattern-match. `PspReader::new` issues 5
    seeks + 4 `read_exact` calls just to open
    ([reader.rs:122-290](../../src/per_sample_caller/psp/reader.rs#L122-L290));
    each block then triggers ~12 `read_exact` per column. Default
    `BufReader::new(file)` allocates an 8 KiB buffer, which is too
    small to absorb the larger compressed columns in a single
    kernel `read(2)`. 64 KiB is the genomic-tooling default and what
    the writer review's L11 recommended for symmetric reasons.
  - **Pattern matched:** io_and_syscalls checklist — *Wrap in
    `BufReader::with_capacity(N, file)`. Default capacity is 8 KiB;
    for large sequential reads bump to 64 KiB or 1 MiB and benchmark*.
  - **Mechanism:** the library API doesn't change — the writer should
    not internally wrap its source (then the caller cannot recover
    the inner `File` for `mmap` / `lseek` later). The hook is at the
    call site; the writer docs advertise the requirement.
  - **Measurement plan:** validated by the `bufreader_file_*_*KiB`
    bench from M2/M5 (sweep 8 KiB / 64 KiB / 1 MiB).
  - **Complexity cost:** docstring-only.
  - **Suggested experiment / fix:**
    ```rust
    /// # Examples
    ///
    /// ```no_run
    /// use std::fs::File;
    /// use std::io::BufReader;
    /// use merge_per_sample_vcfs::per_sample_caller::psp::PspReader;
    ///
    /// let f = File::open("sample.psp")?;
    /// // 64 KiB matches the genomic-tooling default; per-block
    /// // column reads (up to ~tens of KiB compressed at level 9)
    /// // fit in a single buffer fill. For multi-GB files consider
    /// // 1 MiB or larger.
    /// let mut reader = PspReader::new(BufReader::with_capacity(64 * 1024, f))?;
    /// // ...
    /// ```
    ```
    Apply at both
    [reader.rs:108-121](../../src/per_sample_caller/psp/reader.rs#L108-L121)
    and the second example at
    [reader.rs:343-356](../../src/per_sample_caller/psp/reader.rs#L343-L356).

- (project root) — **M3 / L9: Add `examples/dhat_psp_reader.rs` mirroring `examples/dhat_psp_writer.rs`**
  - **Confidence:** High
  - **Hot-path evidence:** the `dhat-heap` feature is already wired
    in [Cargo.toml:72](../../Cargo.toml#L72); examples exist for
    `dhat_pileup`, `dhat_baq`, `dhat_psp_writer`, but not for the
    reader. The reader's allocator pressure is the primary driver
    of L1, L2, L3, L4 — DHAT is the right instrument.
  - **Pattern matched:** methodology checklist — *Use heap profilers
    for allocation hypotheses; wall time alone cannot tell you
    whether a fix removed allocator pressure or moved it*.
  - **Measurement plan:** N/A — this *is* the measurement plan for
    L1/L2/L3/L4.
  - **Complexity cost:** ~30 lines mirroring the existing example.
  - **Suggested experiment / fix:** copy
    [examples/dhat_psp_writer.rs](../../examples/dhat_psp_writer.rs)
    and adapt its harness to drive `PspReader::new` + the records
    walk over a 1 M-record fixture inside a
    `dhat::Profiler::new_heap()` block. Size the workload so ≥ 2
    block reads are triggered.

- [benches/psp_reader_perf.rs](../../benches/psp_reader_perf.rs) — **L10: Add `bufreader_file_1M_{8KiB,64KiB,1MiB}` workloads + a `decode_block_only` sub-bench + split the region bench into open-only / records-only**
  - **Confidence:** High
  - **Hot-path evidence:** the bench file currently has three
    streaming workloads + one bundled region bench. L1 / L8 / S1
    only become measurable with a real-file source; L3's payoff is on
    the per-block decode path which the streaming benches average
    against open-time cost; the region bench bundles open + decode +
    iterate so a single-axis fix is diluted in the headline number.
  - **Pattern matched:** methodology checklist — *Microbenchmarks lie
    about cache behavior; one change per measurement*.
  - **Measurement plan:** add workloads, save baselines as
    `bufreader_file_1M_64KiB`, `decode_block_only`,
    `region_open_only`, `region_records_only`. Per-workload pass
    criterion: each must clear 1 M records/s at the lower 95 % CI.
  - **Complexity cost:** ~80 lines of bench code; one
    `#[doc(hidden)] pub` re-export of `decode_block_payload` (or a
    `decode_block_payload_for_bench` shim) so the bench binary can
    reach the per-column codec in isolation. Bench runtime grows from
    ~5 × 30 s to ~9 × 30 s.
  - **Suggested experiment / fix:** see sketches in
    [io_and_syscalls.md](../../tmp/perf_review_2026-05-13_psp_reader/io_and_syscalls.md)
    and
    [methodology.md](../../tmp/perf_review_2026-05-13_psp_reader/methodology.md).

### Speculative

- [src/per_sample_caller/psp/reader.rs:521](../../src/per_sample_caller/psp/reader.rs#L521) — **S1: `active_chain_slots: Vec<SlotId>` → `SmallVec<[SlotId; 8]>`**
  - **Confidence:** Low
  - **Hot-path evidence:** pattern-match only — the typical active set
    is ~8 slots (per writer's S1 rationale). With `Vec<SlotId>` the
    24-byte header lives in `RecordsIter` but the data lives behind a
    pointer; `SmallVec<[u16; 8]>` keeps the typical case inline.
  - **Pattern matched:** allocations / data_layout — *Tiny, often-empty
    collections are SmallVec candidates*.
  - **Measurement plan:** `phase_chain_heavy_1M` is the workload.
    Threshold ≥ 5 % improvement on the heavy workload, ≤ 1 %
    regression on `snp_typical_3_3M`.
  - **Complexity cost:** new `smallvec` dep (widely used).
  - **Suggested experiment / fix:** two-line type swap; `cap_slots_for_error`
    already takes `&[SlotId]` (SmallVec derefs).

- [src/per_sample_caller/psp/reader.rs:683-714](../../src/per_sample_caller/psp/reader.rs#L683-L714) — **S2: Combine the four `binary_search` loops in `materialise_next_record` into one validate-and-apply pass**
  - **Confidence:** Low
  - **Hot-path evidence:** pattern-match. The function does four
    `binary_search` loops per record over `active_chain_slots`
    (validate-expired, validate-new, apply-expired, apply-new). When
    expired+new are both empty (SNP-typical), all four are no-ops.
    On phase-heavy workloads the 4× redundancy adds up but is dwarfed
    by L3's `decode_list_column` allocation pressure.
  - **Pattern matched:** hot_loops — *Avoid recomputing in tight
    loops*.
  - **Measurement plan:** revisit only if a profile (post-L1+L2+L3)
    still shows the validate/apply loops. Threshold ≥ 2 % on
    `phase_chain_heavy_1M`.
  - **Complexity cost:** loses the explicit fail-fast separation —
    an apply that mutated the active set before a later validation
    fails would need manual rollback. Tricky.

- [src/per_sample_caller/psp/index.rs:32-39](../../src/per_sample_caller/psp/index.rs#L32-L39) — **S3: Project `(chrom_id, first_pos)` into a parallel `Vec<(u32, u32)>` for `find_first_overlapping_block`'s binary search**
  - **Confidence:** Low
  - **Hot-path evidence:** pattern-match only — `find_first_overlapping_block`
    is `O(log n_blocks)` per region query, not per record. Today's
    `BlockIndexEntry` is 24 B; binary search compares only the
    leading 8 bytes. For 100 k+ block files (large WGS samples),
    the parallel-keys array fits 3× more probe targets per cache
    line.
  - **Pattern matched:** data_layout — *Project hot fields into a
    cache struct*.
  - **Measurement plan:** add `region_open_only_1M_blocks` micro-
    bench with a synthetic 1 M-entry index. Threshold ≥ 20 %
    wall-clock improvement on the binary-search step.
  - **Complexity cost:** one extra `Vec<(u32, u32)>` field on
    `PspReader`, kept in sync at construction. Dual-source-of-truth
    risk if the index ever becomes mutable.

- [src/per_sample_caller/psp/reader.rs:357-369](../../src/per_sample_caller/psp/reader.rs#L357-L369) — **S4: `mmap` experiment for region queries that touch only a few blocks of a large file**
  - **Confidence:** Low
  - **Hot-path evidence:** pattern-match only. Region queries
    `seek(block_offset)` per block boundary; against `BufReader<File>`
    each seek invalidates the buffer (one wasted-fill + one
    `lseek(2)` per block). `mmap` sidesteps both: random-access reads
    land in already-cached pages, kernel readahead does the heavy
    lifting.
  - **Pattern matched:** io_and_syscalls — *`mmap` for random-access
    reads of large files*.
  - **Measurement plan:** add `region_window_file_mmap` vs
    `region_window_file_bufreader` benches over a 5–10 block window;
    threshold ≥ 20 % faster on mmap.
  - **Complexity cost:** `memmap2` dependency + SIGBUS-on-truncation
    caveat. Document explicitly.

- [src/per_sample_caller/psp/reader.rs:417-434](../../src/per_sample_caller/psp/reader.rs#L417-L434) — **S5: Lending iterator with `BorrowedPileupRecord<'block>` to eliminate per-emit allocations**
  - **Confidence:** Low
  - **Hot-path evidence:** pattern-match only. The `mem::take` shape
    today is the right choice given `PileupRecord` / `AlleleObservation`
    own their `Vec<u8>` / `Vec<SlotId>`. A borrowed view tied to the
    block lifetime would let the reader keep column data in a per-
    block arena and emit slices, eliminating both the per-allele
    `to_vec()` and the per-row chain-slot Vec.
  - **Pattern matched:** allocations / data_layout — *Owning emission
    contract forces a per-emit allocation*.
  - **Measurement plan:** build a separate `consume_records` API
    yielding `BorrowedPileupRecord<'block>`; bench head-to-head
    against `records()`. Threshold ≥ 20 % to merit the API surface
    widening.
  - **Complexity cost:** **Significant.** Two parallel public types
    or a generic shape; new lifetime in the iterator (lending iterator
    pattern, GATs since Rust 1.65 but `Iterator::Item` is the wrong
    trait — needs `LendingIterator` or a hand-rolled consumer
    protocol). API churn for downstream callers.

- (perf-tooling) — **S6: Document a sampling-profile escape hatch in `ia/skills/performance_review/methodology.md` (or a sibling note)**
  - **Confidence:** High
  - **Hot-path evidence:** verbatim from `_evidence.md`: this review
    capped every code-level finding at Likely because no profile is
    reachable in the sandbox.
  - **Pattern matched:** methodology — *Profile first; pattern-match
    second*.
  - **Measurement plan:** N/A — process change.
  - **Complexity cost:** ~10 lines of bash + a paragraph in the skill
    doc.
  - **Suggested experiment / fix:** see M1 in
    [methodology.md](../../tmp/perf_review_2026-05-13_psp_reader/methodology.md).

### Note

- [src/per_sample_caller/psp/reader.rs:1085-1098](../../src/per_sample_caller/psp/reader.rs#L1085-L1098) — `DecodedColumn` enum is 32 bytes (24 B max variant + tag + padding); the 12-arm match in `decode_block_payload` runs once per column per block, not per record. Below the action threshold.
- [src/per_sample_caller/psp/block.rs:566-571](../../src/per_sample_caller/psp/block.rs#L566-L571) — `ColumnManifestEntry` is `repr(Rust)`; compiler reorders `(u32, u32, u16)` to 12 bytes (4+4+2+2 padding). Already optimal. Could pin with a `const _: () = assert!(std::mem::size_of::<ColumnManifestEntry>() == 12);` near the type.
- [src/per_sample_caller/psp/reader.rs:496](../../src/per_sample_caller/psp/reader.rs#L496) — `Option<DecodedBlock>` is niche-optimised (the `Vec` ptr is `NonNull`, compiler stuffs `None` discriminant there). Per-iteration branch is fully predictable.
- [src/per_sample_caller/psp/registry.rs](../../src/per_sample_caller/psp/registry.rs) — `ColumnDef` fields are `&'static str` + scalars; registry walks allocation-free. Linear-scan `lookup_by_tag` over 12 entries is correct (the function comment names this).
- [src/per_sample_caller/psp/trailer.rs:39-45](../../src/per_sample_caller/psp/trailer.rs#L39-L45), [src/per_sample_caller/psp/varint.rs](../../src/per_sample_caller/psp/varint.rs) — `decode_trailer` uses a stack `[u8; 32]`; varint codecs are pure / no allocation.
- [src/per_sample_caller/psp/reader.rs:870-917](../../src/per_sample_caller/psp/reader.rs#L870-L917) — `read_block_header`'s grow loop re-walks from byte 0 each iteration (quadratic in pathological case). Cold per-block; first chunk satisfies typical headers; bounded at 64 KiB. Not actionable today.
- [src/per_sample_caller/psp/reader.rs:457-481](../../src/per_sample_caller/psp/reader.rs#L457-L481) — `RangeClamp` `match` arms compile to cmovs. Don't rewrite by hand.
- [src/per_sample_caller/psp/reader.rs:1295-1301](../../src/per_sample_caller/psp/reader.rs#L1295-L1301) — `cap_slots_for_error` `.collect()` on the error path; reached at most once per `RecordsIter`. Cold by construction.
- [src/per_sample_caller/psp/reader.rs:1184-1187,1131-1132](../../src/per_sample_caller/psp/reader.rs#L1131-L1132) — `format!` / `to_string()` allocations in `decode_one_column` live inside `.map_err(...)` closures (only run on error) or inside the cold budget-check arm. Not on the success path.
- [src/per_sample_caller/psp/reader.rs:173,209](../../src/per_sample_caller/psp/reader.rs#L173) — open-sequence `vec![0u8; ...]` allocations are once-per-file, pre-sized correctly.
- [src/per_sample_caller/psp/index.rs:101-127](../../src/per_sample_caller/psp/index.rs#L101-L127) — `decode_index` reads the entire encoded index in one `read_exact` and pre-sizes its `Vec<BlockIndexEntry>`. Already correct.

## 6. Out-of-scope observations

- [benches/psp_reader_perf.rs](../../benches/psp_reader_perf.rs) — the bench was added in this session and is intentionally minimal. L9 / L10 trace the gaps; landing them is the next-PR menu before any §5 L code-level fix is merged with confidence.
- The eventual CLI / merge driver — when written — must wrap `File` sources in `BufReader::with_capacity(64 KiB, file)` per L8 (and consider `mmap` for region-heavy workloads per S4). Filed as a follow-up issue, not blocking these findings.
- The reader's symmetric `[profile.profiling]` is the once-correctly-set-up-stays-correct kind of build-config change; do it at the same time as the Cargo.toml-level mimalloc A/B (when it lands) so a single PR carries the build-config delta.
- The `samply`-blocked sandbox is the same constraint that hit the writer review (which solved it by running the container-built binary on the host). The host's perf_event_open is *also* blocked here, suggesting either (a) a deeper sandbox layer than the writer review faced, or (b) a regression in the reviewer's environment. Either way, every finding here that cites the writer profile as evidence is an *educated guess* until a real reader profile exists — file an environment ticket as a follow-up.

## 7. What's already good

- **Pre-sized scalar columns.** `decode_scalar_column`,
  `decode_varint_column`, and `decode_list_column` all use
  `Vec::with_capacity(expected_count)` from the cardinality-derived
  count ([block.rs:183](../../src/per_sample_caller/psp/block.rs#L183),
  [block.rs:233](../../src/per_sample_caller/psp/block.rs#L233),
  [block.rs:325](../../src/per_sample_caller/psp/block.rs#L325)) —
  no per-element grow on the column-payload path. The writer's H1
  fix (pre-size accumulator columns) does not have a reader analog
  because the reader was already doing the right thing.
- **Sorted-`Vec<SlotId>` active set.** `RecordsIter::active_chain_slots`
  ([reader.rs:521](../../src/per_sample_caller/psp/reader.rs#L521))
  matches the writer's `IngestState::active_slots` (`writer.rs:130`)
  pattern: one contiguous allocation, `binary_search` membership,
  block-start snapshot equality is a `Vec == Vec` memcmp.
- **`mem::take` on emission.** `materialise_next_record` already
  uses `std::mem::take(&mut block.allele_seqs[j])` to move the
  per-allele `Vec<u8>` out of the block instead of cloning
  ([reader.rs:765,773,780,781](../../src/per_sample_caller/psp/reader.rs#L765-L781)) —
  the right shape given today's `Vec<Vec<…>>` `DecodedBlock` layout.
  (CSR adoption per L3 will change the optimal shape; the current
  code is exactly the right thing for the current layout.)
- **Build configuration is tight.** Same as the writer review:
  `lto = "fat"`, `codegen-units = 1`, `debug = "line-tables-only"`,
  `panic = "abort"` on `[profile.release]`; `target-cpu = x86-64-v3`
  on Linux/x86_64; `rust-toolchain.toml` pinned at 1.95.

### Author response convention

Address each finding by its identifier (e.g., "L1", "L2") with one
of: `applied in <commit>` / `experiment shows no gain — closing` /
`disputed because …` / `deferred to <issue>` / `won't fix because …`.
The "experiment shows no gain" path is expected and welcome — that
is what the measurement plan is for.
