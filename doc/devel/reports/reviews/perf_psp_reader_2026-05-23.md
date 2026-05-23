# Performance Review: psp_reader (cohort hot path, post-H3)
**Date:** 2026-05-23
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** `.psp` reader at [src/per_sample_pileup/psp/](../../../../src/per_sample_pileup/psp/) on the cohort var-calling hot path
**Verdict:** **Apply the listed wins** — two Hot-path findings (CSR ragged-column collapse + `SeekFrom::Current` for per-block seeks). Three of six categories converged on the same `DecodedBlock` ragged-column site; the I/O fix is the 2026-05-20 H4 that was never applied.
**Hot-path evidence:** `perf record --call-graph=dwarf,16384 -F 997` against `examples/profile_cohort_e2e --threads 4` on real tomato N=10 cohort (post-H3 commit `503b7cc`); 142,374 samples; verbatim output in [tmp/perf_review_2026-05-23_ref_fetcher/psp_profile_summary.txt](../../../../tmp/perf_review_2026-05-23_ref_fetcher/psp_profile_summary.txt) and per-category sub-reports in [tmp/perf_review_2026-05-23_psp_reader/](../../../../tmp/perf_review_2026-05-23_psp_reader/).

---

## 1. Scope and constraints

- **What was reviewed:** the on-disk `.psp` *reader* — `PspReader`, `RecordsIter`, `DecodedBlock`, `region_records`, and the per-block decoders. Writer out of scope (separate review).
- **Reviewed against:** commit `503b7cc` (the H3 commit on `main`).
- **Throughput / latency targets, input sizes, target hardware:** the reader serves N=10 cohort samples × per-chrom rayon workers in `process_one_chromosome`. Tomato SL4.0, 13 chroms (~85 Mbp ch00 + ~50 Mbp average), ~1.4 M records emitted + dropped per cohort run. Production target: 27.7 s cohort var-calling on 13-thread server. Today's 4-thread host baseline (i7-1260P P-cores) at H3 = mean 40.34–40.60 s. PSP-related cpu_core share: ~18.8 % (the second-largest consumer after DUST).
- **Hot-path evidence available:** sampling profile (perf with dwarf call-graph). Both `cpu_atom` (E-cores) and `cpu_core` (P-cores) event types captured. PSP top symbols (cpu_atom / cpu_core, from [psp_profile_summary.txt](../../../../tmp/perf_review_2026-05-23_ref_fetcher/psp_profile_summary.txt)):
  ```
   2.38%   6.65%  ZSTD_decompressSequences_bmi2.constprop.0
   2.29%   6.55%  RecordsIter::next
   1.50%   1.94%  decode_list_column
   0.58%   1.66%  ZSTD_decompressMultiFrame
   0.51%   0.78%  decode_varint_column
   0.29%   0.74%  decode_bytes_split
   0.19%   0.50%  drop_in_place<Option<DecodedBlock>>
   0.12%   0.28%  ZSTD_buildFSETable_body_bmi2
  ```
  And the allocator-leaf attribution count (P-core samples one frame above `_int_malloc` / `cfree` / `malloc_consolidate` / etc):
  ```
      2099  pop_var_caller::per_sample_pileup::psp::block::decode_bytes_split
        10  pop_var_caller::per_sample_pileup::psp::block::decode_varint_column
         4  pop_var_caller::per_sample_pileup::psp::block::decode_list_column
  ```
  `decode_bytes_split` = **87 % of project-attributable P-core allocator samples in the entire cohort pipeline.** The next-biggest project-side allocator caller in the whole cohort (`PerGroupMerger::unify_alleles` at 142 samples) is 15× smaller.
- **In-scope files:**
  - [src/per_sample_pileup/psp/reader.rs](../../../../src/per_sample_pileup/psp/reader.rs) — `PspReader`, `RecordsIter`, `DecodedBlock`, `region_records`, `materialise_next_record` (2,248 LOC)
  - [src/per_sample_pileup/psp/block.rs](../../../../src/per_sample_pileup/psp/block.rs) — block decoders (`decode_scalar_column_pod`, `decode_varint_column`, `decode_list_column`, `decode_bytes_split`)
  - [src/per_sample_pileup/psp/varint.rs](../../../../src/per_sample_pileup/psp/varint.rs) — varint primitive
  - [src/per_sample_pileup/psp/index.rs](../../../../src/per_sample_pileup/psp/index.rs) — chrom-block index (used by `region_records`)
  - [src/per_sample_pileup/psp/header.rs](../../../../src/per_sample_pileup/psp/header.rs), [trailer.rs](../../../../src/per_sample_pileup/psp/trailer.rs), [registry.rs](../../../../src/per_sample_pileup/psp/registry.rs) — cold per record
  - [benches/psp_reader_perf.rs](../../../../benches/psp_reader_perf.rs) — maintained `Cursor`-backed bench
  - [src/pop_var_caller/cohort_driver.rs](../../../../src/pop_var_caller/cohort_driver.rs) — constructs the per-sample `PspReader<BufReader<File>>` chain in `process_one_chromosome`
- **Deliberately out of scope:**
  - PSP **writer** (`writer.rs`) — separate review surface, separate hot path.
  - **zstd library internals** — third-party (`ZSTD_decompressSequences_bmi2` at 6.65 % cpu_core is consequential but we can't tune the zstd crate's inner loop).
  - The **cohort orchestration** layer (already covered by the 2026-05-20 perf review).
- **Categories dispatched:**
  - `methodology` — Always.
  - `allocations` — `decode_bytes_split` is the project's #1 allocator caller.
  - `data_layout` — `drop_in_place<DecodedBlock>` 0.50 % cpu_core surfaces a ragged-column tear-down cost.
  - `concurrency` — Per-worker ownership per `process_one_chromosome`; expected `No findings`.
  - `hot_loops` — `RecordsIter::next` 6.55 % cpu_core, `decode_list_column` 1.94 %, `decode_varint_column` 0.78 %.
  - `io_and_syscalls` — `RecordsIter::next` does `SeekFrom::Start` per block; the 2026-05-20 H4 finding flagged this.

## 2. Verdict

**Apply the listed wins.** Two Hot-path findings, both profile-evidenced, both with contained complexity, both with clear measurement plans. The CSR refactor (H1) attacks the largest project allocator caller in the entire pipeline; the seek fix (H2) recovers a buffer that the cohort driver carefully constructs and the reader silently throws away on every block transition.

The Hot-path findings stack with the just-landed ref-fetcher H1 + H3 — there's no overlap. H3 of this review (per-allele bounds-check hoist) is gated by H1 of this review (CSR refactor); they should land together as a single PR.

## 3. Measurement plan

1. **H2 first (smaller blast radius, unblocks H1 measurement).** `seek_to_offset(target)` helper using `SeekFrom::Current(delta)` so a wrapped `BufReader` short-circuits when target is in-window; replace `SeekFrom::Start` at [reader.rs:587](../../../../src/per_sample_pileup/psp/reader.rs#L587) and [reader.rs:807](../../../../src/per_sample_pileup/psp/reader.rs#L807). Re-run cohort var-calling 3× back-to-back. Threshold: ≥ 1 s off the 40.34–40.60 s wall-time on the post-H3 baseline. Cross-check with `strace -c -e read,lseek ./profile_cohort_e2e ...` — expect a substantial drop in `read(2)` syscall count.
2. **H1: CSR refactor.** Replace `DecodedBlock.allele_seqs: Vec<Vec<u8>>` with `(allele_seq_data: Vec<u8>, allele_seq_offsets: Vec<u32>)`; same for `allele_chain_ids`. Add caller-owned `data` + `offsets` buffers on `RecordsIter`. Re-run cohort var-calling and capture a perf profile. Threshold: ≥ 5 % wall-time improvement *and* `decode_bytes_split` allocator-leaf samples drop from 2099 to <20 (mostly from cold paths). The `drop_in_place<Option<DecodedBlock>>` line should drop from 0.50 % cpu_core to noise.
3. **H3: bounds-check hoist.** Apply leading `assert!(allele_end <= block.<column>.len())` for the 9 ragged columns at the top of `materialise_next_record`. Cross-check via `cargo asm` — expect 8–9 `panic_bounds_check` branches in the per-allele loop body to collapse to one each, outside the loop.
4. **Bench-side: add a file-backed `psp_reader_perf` group** (L6 below). Current bench uses `Cursor<&[u8]>` which can't measure H2's effect. Write the SNP / phase / multi fixtures to a `tempfile::NamedTempFile` once outside the timed region, then time `PspReader::new(BufReader::with_capacity(64 KiB, File::open(&path)?))`. Without this, H2 has only the cohort-e2e wall-time as a signal — and the host has known ~2× variance across sessions.
5. **Bench-side: add `examples/dhat_psp_reader.rs`** (L7 below). DHAT directly counts allocations; after H1 lands, the per-block allocations should drop from "~N_alleles" to "2" per column.

## 4. Build / toolchain configuration

Already in place — do not touch:
- `[profile.release] lto = "fat"`, `codegen-units = 1`, `panic = "abort"`, `debug = "line-tables-only"`.
- `[profile.bench] inherits = "release"` with `debug = true`.
- `.cargo/config.toml` `target-cpu = x86-64-v3` (AVX2 + FMA + BMI2).
- 1.95.0 toolchain pin.
- `alloc-mimalloc` feature wired for benches.

The two cross-routed methodology asks from the previous review are still open: a `[profile.release-with-debug]` for reproducible production-binary profiling, and wiring `mimalloc::MiMalloc` as `#[global_allocator]` in `src/main.rs` behind `alloc-mimalloc`. Both already tracked under the 2026-05-23 ref-fetcher review.

## 5. Code-level findings

### Hot-path

#### H1: [src/per_sample_pileup/psp/reader.rs:432](../../../../src/per_sample_pileup/psp/reader.rs#L432), [reader.rs:440](../../../../src/per_sample_pileup/psp/reader.rs#L440) + [block.rs:521](../../../../src/per_sample_pileup/psp/block.rs#L521), [block.rs:412](../../../../src/per_sample_pileup/psp/block.rs#L412) — CSR collapse of `DecodedBlock`'s ragged `Vec<Vec<T>>` columns

- **Confidence:** High (three category sub-agents converge on this site with consistent profile evidence)
- **Hot-path evidence:** **decode_bytes_split holds 2099 / 2399 = 87 % of project-attributable P-core allocator samples** in the entire cohort var-calling profile. cpu_core self-time numbers:
  ```
   1.50%   1.94%  decode_list_column                      (allocates Vec::with_capacity(k) per allele entry, block.rs:412)
   0.29%   0.74%  decode_bytes_split                      (allocates Vec::with_capacity + per-allele to_vec, block.rs:521)
   0.19%   0.50%  drop_in_place<Option<DecodedBlock>>     (tears down the two Vec<Vec<T>> columns + their inner buffers)
  ```
  The 0.74 % self-time in `decode_bytes_split` understates the impact: the work happens in the allocator (`_int_malloc`, `cfree`, `_int_free_chunk`, `malloc_consolidate` — summing ~14–18 % combined of cpu_core in the post-H3 profile, with this site as the dominant project-side cause).
- **Pattern matched (3 categories):**
  - **Allocations:** "Allocations belong outside hot loops." Per-allele `Vec<u8>` (in `decode_bytes_split`) and `Vec<ChainId>` (in `decode_list_column`) inside the per-block decode loop.
  - **Data layout:** "Pointer chases fragment the cache. If the elements are uniform, store inline." `Vec<Vec<T>>` requires a pointer chase per row + per-row heap allocation.
  - **Hot loops:** "Cache-friendly memory layout enables vectorization." The per-allele assembly loop in `materialise_next_record` reads 9 separate `block.allele_<col>[j]` columns through 9 distinct `Vec` data pointers — AoS access pattern through SoA storage.
- **Mechanism:** `DecodedBlock` (reader.rs:424) holds two ragged columns:
  ```rust
  allele_seqs: Vec<Vec<u8>>,           // reader.rs:432
  allele_chain_ids: Vec<Vec<ChainId>>, // reader.rs:440
  ```
  `decode_bytes_split` (block.rs:472–525) does `out.push(bytes[cursor..cursor + n].to_vec())` per allele inside the per-block loop. `decode_list_column` does `Vec::with_capacity(k)` per row (block.rs:412). On a SNP-typical block with thousands of alleles, that's thousands of small `malloc` + `free` pairs per block per sample. `materialise_next_record` then `std::mem::take`s individual inner Vecs out into the emitted `AlleleObservation`, leaving zombie empty `Vec`s in `DecodedBlock` to be torn down at `drop` time — which the profile attributes at 0.50 % cpu_core.

  **CSR collapse:** replace each ragged column with `(data: Vec<T>, offsets: Vec<u32>)`. Both buffers live on `RecordsIter` (next to the existing `decompressed_scratch`) and are reused across blocks. `decode_bytes_split_csr` / `decode_list_column_csr` write into caller-owned buffers — symmetric with the writer's existing `encode_list_column_csr` (block.rs:352). The materialiser slices `data[offsets[j]..offsets[j+1]]` instead of `mem::take`ing inner Vecs. Per-block allocations collapse from `~N_alleles` to `2` for each column.

  Two-stage payoff:
  - **Stage A (no public API change):** `materialise_next_record` still ends with `seq: data[a..b].to_vec()` so `AlleleObservation.seq: Vec<u8>` is unchanged. Net allocation count downstream is the same (one Vec per emitted allele), but the decode side no longer round-trips through a per-row heap Vec. The 87 % allocator pressure share is the **decode side** — Stage A captures that.
  - **Stage B (cross-cutting):** Change `AlleleObservation.seq` from `Vec<u8>` to `Box<[u8]>` (or a borrow into a per-block arena). Pulls the remaining one-allocation-per-allele cost. Touches `pileup::AlleleObservation`'s public surface — defer until Stage A confirms the mechanism.
- **Measurement plan:** Step 2 of §3.
- **Complexity cost:** **Stage A:** `decode_bytes_split` / `decode_list_column` gain CSR-shaped twin entry points (signatures take `&mut Vec<u8>` data + `&mut Vec<u32>` offsets). `DecodedBlock` swaps two fields. `materialise_next_record` swaps `mem::take` for slice. `RecordsIter` gains 4 new scratch fields. No public API change. The existing non-CSR functions can remain as thin wrappers for tests / non-hot callers. **Stage B:** touches `AlleleObservation.seq` + downstream `unify_alleles` / `vcf_writer`. Gate behind dhat data confirming Stage A.
- **Suggested experiment / fix:**
  ```rust
  // block.rs — Stage A: CSR decode for bytes columns
  pub fn decode_bytes_split_csr(
      bytes: &[u8],
      lengths: &[u64],
      column_name: &str,
      max_entry_len: Option<u64>,
      slab_buf: &mut Vec<u8>,
      offsets_buf: &mut Vec<u32>,
  ) -> Result<(), PspReadError> {
      // same bounds / max_entry_len checks as today
      slab_buf.clear();
      slab_buf.extend_from_slice(bytes);       // one memcpy
      offsets_buf.clear();
      offsets_buf.reserve(lengths.len() + 1);
      offsets_buf.push(0);
      let mut cursor: u64 = 0;
      for &len in lengths {
          cursor = cursor.checked_add(len).ok_or(/* Truncated */)?;
          offsets_buf.push(cursor as u32);     // u32 fits — block-level cap enforces it
      }
      Ok(())
  }

  // reader.rs — DecodedBlock change
  struct DecodedBlock {
      // ...
      // Replace allele_seqs: Vec<Vec<u8>>
      allele_seq_data: Vec<u8>,
      allele_seq_offsets: Vec<u32>,         // len == n_total_alleles + 1
      // Replace allele_chain_ids: Vec<Vec<ChainId>>
      allele_chain_ids_data: Vec<ChainId>,
      allele_chain_ids_offsets: Vec<u32>,
  }

  // reader.rs — materialise_next_record
  for j in allele_start..allele_end {
      let s = block.allele_seq_offsets[j] as usize;
      let e = block.allele_seq_offsets[j + 1] as usize;
      let cs = block.allele_chain_ids_offsets[j] as usize;
      let ce = block.allele_chain_ids_offsets[j + 1] as usize;
      alleles.push(AlleleObservation {
          seq: block.allele_seq_data[s..e].to_vec(),               // Stage A: one alloc per emit
          support: AlleleSupportStats { /* unchanged */ },
          chain_ids: block.allele_chain_ids_data[cs..ce].to_vec(), // Stage A: one alloc per emit
      });
  }
  ```
  Sketches in [allocations.md](../../../../tmp/perf_review_2026-05-23_psp_reader/allocations.md) and [data_layout.md](../../../../tmp/perf_review_2026-05-23_psp_reader/data_layout.md).

#### H2: [src/per_sample_pileup/psp/reader.rs:587](../../../../src/per_sample_pileup/psp/reader.rs#L587) + [reader.rs:807](../../../../src/per_sample_pileup/psp/reader.rs#L807) — `SeekFrom::Start` per block discards the 64 KiB BufReader buffer between contiguous blocks (the 2026-05-20 H4 finding, still unapplied)

- **Confidence:** High (the H4 finding was profile-evidenced in 2026-05-20 and the current source still shows the seek; confirmed by reading the file)
- **Hot-path evidence:** `RecordsIter::next` is at 2.29 % / **6.55 % cpu_core** — the largest non-zstd PSP symbol. Two distinct per-block seek sites discard the BufReader cache that `process_one_chromosome` deliberately constructs at 64 KiB.
- **Pattern matched (io_and_syscalls):** "`File::read` / `File::write` without buffering are a per-call syscall ... wrap in `BufReader::with_capacity`." The wrap exists. The defect is that `SeekFrom::Start` against `BufReader` discards the buffer unconditionally; `BufReader::seek_relative` / `BufReader`'s special-cased `SeekFrom::Current` short-circuit when the target lies in the existing window.
- **Mechanism:** Blocks are written **contiguously** on disk by the writer (`writer.rs:507,525` — `block_offset = sink_offset`, `sink_offset += written`). The reader's own `block_byte_budget` (reader.rs:1218–1225) computes `next_offset - this_offset` as the block's geographic byte budget. After exhausting block N, the source position is already at `entry.block_offset` for block N+1 — but `load_next_block` ([reader.rs:585–588](../../../../src/per_sample_pileup/psp/reader.rs#L585-L588)) issues `SeekFrom::Start(entry.block_offset)`, a zero-distance jump that **`BufReader` still treats as a buffer-invalidating absolute seek**. Every block transition forces the next column read to re-issue `read(2)` for bytes that were in user-space microseconds earlier.

  Same anti-pattern at [reader.rs:806–808](../../../../src/per_sample_pileup/psp/reader.rs#L806-L808): `read_block_header` over-reads up to one 4 KiB chunk past the actual header, then rewinds to `header_start + consumed` via `SeekFrom::Start`. The rewind distance is at most ~4 KiB (well inside the 64 KiB buffer), but `SeekFrom::Start` invalidates the whole BufReader buffer anyway.

  The fix on both sites: a `seek_to_offset(source, target, context)` helper that uses `SeekFrom::Current(target - current)`. `BufReader`'s `Seek::seek` impl special-cases `SeekFrom::Current(n)` for `|n| < buffer.len()` to preserve the buffer (per std's source). For an unwrapped `File` this degenerates to a normal `lseek` — no slower than today.
- **Measurement plan:** Step 1 of §3. Threshold: ≥ 1 s off the 40 s wall, *and* `RecordsIter::next` cpu_core drops by ≥ 1 pp, *and* `strace -c` shows ≥ 50 % drop in `read(2)` count.
- **Complexity cost:** ~10 LOC helper. No new dependencies, no `unsafe`. The PSP reader's generic bound `R: Read + Seek` is preserved (the helper uses `SeekFrom::Current` which is on the bare `Seek` trait — does NOT require a `BufRead` bound).
- **Suggested experiment / fix:**
  ```rust
  // reader.rs — new helper
  /// Move `source` to absolute `target`. Uses `SeekFrom::Current(delta)`
  /// so a wrapped `BufReader` keeps its buffer when the target lies in
  /// the existing window. `BufReader::seek` with `SeekFrom::Start`
  /// discards the buffer unconditionally; with `SeekFrom::Current(n)`
  /// and `|n| < buffer.len()` it preserves it.
  fn seek_to_offset<R: Read + Seek>(
      source: &mut R,
      target: u64,
      context: &'static str,
  ) -> Result<(), PspReadError> {
      let current = source.stream_position().map_err(io_err(context))?;
      if current == target {
          return Ok(());
      }
      let delta = (target as i64) - (current as i64);
      source.seek(SeekFrom::Current(delta)).map_err(io_err(context))?;
      Ok(())
  }
  ```
  Use it at the two flagged sites. The L5 below (`stream_position` elimination) drops the `stream_position()` call by threading `entry.block_offset` into `read_block_header` — once that lands, `seek_to_offset` is only called from `load_next_block` and the helper's `stream_position()` overhead becomes the sole `lseek` of the per-block path.

#### H3 (gated by H1): [src/per_sample_pileup/psp/reader.rs:666](../../../../src/per_sample_pileup/psp/reader.rs#L666) — Per-allele inner loop hits 9 independent `Vec` bounds checks; one `assert!` per column hoists them all

- **Confidence:** Medium (mechanism is well-understood; the gain depends on H1 having flattened the columns first)
- **Hot-path evidence:** `RecordsIter::next` at 6.55 % cpu_core. `materialise_next_record` is single-call-site `#[inline]`-able and collapses into the `next` symbol after monomorphisation against `BufReader<File>`. The per-allele body reads from 9 separate `Vec`s through `&mut block`: `allele_seqs[j]`, `allele_obs_count[j]`, `allele_q_sum_log[j]`, `allele_fwd_count[j]`, `allele_placed_left_count[j]`, `allele_placed_start_count[j]`, `allele_mapq_sum[j]`, `allele_mapq_sum_sq[j]`, `allele_chain_ids[j]`. LLVM cannot share a bounds proof across disjoint Vecs.
- **Pattern matched (hot_loops):** "A leading `assert!(idx < slice.len())` hoists the bounds check out of the loop body — counterintuitively, more `assert!`s can be faster."
- **Mechanism:** 8–9 leading `assert!(allele_end <= block.<column>.len())` give LLVM a single dominator condition per column. The compiler then proves the in-loop `block.<column>[j]` accesses are safe and eliminates the `panic_bounds_check` calls from the loop body. The asserts are dead in practice (block decode already enforces `sum(n_alleles) == n_total_alleles`).
- **Measurement plan:** `cargo asm --rust ... materialise_next_record` to confirm `panic_bounds_check` edges in the `j` loop collapse to one each outside the loop. Microbench `materialise_next_record` over a 1k-record block; threshold ≥ 3 % drop in per-record cpu_core *or* visible asm change.
- **Complexity cost:** 9 `assert!` lines (or one helper `assert_block_lengths(block, allele_end)`). No new types, no `unsafe`. The asserts can never fire on well-formed input — `decode_block_payload` already validates the column lengths.
- **Suggested experiment / fix:** sketch in [hot_loops.md](../../../../tmp/perf_review_2026-05-23_psp_reader/hot_loops.md). After H1 (CSR), the column count to assert drops from 9 to 7 (the two ragged columns become `data` + `offsets` pairs); the offsets vec only needs one assert (`allele_end + 1 <= offsets.len()`).

### Likely

#### L1: [src/per_sample_pileup/psp/block.rs:412](../../../../src/per_sample_pileup/psp/block.rs#L412) — `decode_list_column` per-element loop blocks LE-slab cast even when `T: Pod`

- **Confidence:** Medium
- **Hot-path evidence:** `decode_list_column` at 1.94 % cpu_core. The writer's `encode_list_column_csr` (block.rs:352) already does the symmetric LE-slab cast on encode; the reader has no fast path on decode.
- **Pattern matched (hot_loops):** "Manual SIMD is rarely the right answer; fix what blocks autovec instead." The per-element `T::decode_le(&bytes[cursor..])` + `cursor += 8` + `list.push(v)` body looks autovectorizable but isn't due to the per-iteration cursor write + `Vec::push`. On `target_endian = "little"` with `T: Pod` (true for `ChainId = u64`, the only list-column instantiation today), the whole k-element run is one `bytemuck::cast_slice::<u8, T>(&bytes[cursor..cursor + k * width]).to_vec()`.
- **Mechanism:** Add a `decode_list_column_pod<T: WireScalar + Pod>` fast path gated on little-endian. After H1, this fix also writes directly into the CSR `data: Vec<ChainId>` buffer.
- **Measurement plan:** Microbench over `phase_chain_heavy_1M` (the workload with non-trivial list density). Threshold: ≥ 30 % self-time reduction on `decode_list_column` for lists-of-≥2 fixtures.
- **Complexity cost:** One new fast-path entry point, gated on `T: Pod` + `cfg(target_endian = "little")`. Big-endian fallback keeps the per-element loop.

#### L2: [src/per_sample_pileup/psp/reader.rs:780-783](../../../../src/per_sample_pileup/psp/reader.rs#L780-L783) — `stream_position()` inside `read_block_header` is an extra `lseek(2)` per block

- **Confidence:** Medium
- **Hot-path evidence:** Pattern-match. `stream_position()` on `BufReader<File>` calls `seek(SeekFrom::Current(0))`, which the std impl special-cases to NOT discard the buffer but still issues an `lseek(2)` to learn the file offset.
- **Mechanism:** `let header_start = source.stream_position()?;` exists only to compute the rewind target at reader.rs:807. `load_next_block` *just* seeked to `entry.block_offset` one line earlier — the caller already knows the offset. Threading `entry.block_offset` into `read_block_header` removes the call.
- **Measurement plan:** `strace -c -e lseek` before/after; expect a drop equal to the block count per cohort run.
- **Complexity cost:** One new `u64` parameter on `read_block_header`. Private function, local change.

#### L3: [src/per_sample_pileup/psp/reader.rs:374](../../../../src/per_sample_pileup/psp/reader.rs#L374) — `region_records` first-block jump

- **Confidence:** Medium
- **Hot-path evidence:** Subsumed by H2 — the cohort driver calls `region_records(chrom_id, 1, chrom_length)` 130 times per cohort run (10 samples × 13 chroms). The first block jump per region is a large absolute seek; subsequent block transitions inside the region are the same H2 problem.
- **Mechanism:** Once H2's `seek_to_offset` helper lands, the first region jump falls through to a real `SeekFrom::Current` because the displacement exceeds the buffer — that's correct (no in-buffer hit possible for the first jump). Subsequent transitions are recovered by H2.
- **Measurement plan:** Subsumed by H2.
- **Complexity cost:** Zero beyond H2.

#### L4: [src/per_sample_pileup/psp/reader.rs:621-693](../../../../src/per_sample_pileup/psp/reader.rs#L621-L693) — AoS-through-SoA scan in `materialise_next_record` (gated by H1)

- **Confidence:** Medium
- **Hot-path evidence:** `RecordsIter::next` at 6.55 % cpu_core. The per-allele loop reads from 7 fixed-width per-allele `Vec<T>` columns — that's 7 distinct stride streams for the hardware prefetcher, which can only track ~4–6 streams at once. Each `j++` crosses a fresh cache line in 7 distinct address ranges.
- **Pattern matched (data_layout):** "AoS vs SoA is a function of access pattern. AoS is correct when the hot loop reads several fields of each element together."
- **Mechanism:** Pack the 7 fixed-width per-allele attributes (`obs_count`, `q_sum_log`, `fwd_count`, `placed_left_count`, `placed_start_count`, `mapq_sum`, `mapq_sum_sq`) into one `Vec<DecodedAlleleStatsRow>` where `DecodedAlleleStatsRow` is `#[repr(C)]` and sized to one cache line (~40 B). The decode side gains a small SoA-to-AoS gather (paid once per block, not per allele).
- **Measurement plan:** Gated behind H1. `cargo bench` on `snp_typical_3_3M`. Threshold: ≥ 3 % wall-time improvement *after H1 has landed*. If <3 %, revert — `RecordsIter::next` may not be cache-miss-bound (a `samply record` with cache-miss counters would confirm).
- **Complexity cost:** One new private `repr(C)` struct. Internal repacking — `AlleleSupportStats` public type is unchanged.

#### L5: [src/per_sample_pileup/psp/pileup/mod.rs (AlleleObservation)](../../../../src/per_sample_pileup/pileup/mod.rs) — `AlleleObservation` is ~88 B; `chain_ids: Vec<ChainId>` is empty for most alleles

- **Confidence:** Low (cross-cutting; already on 2026-05-20 L2)
- **Hot-path evidence:** Pattern-match — `AlleleObservation` straddles 2 cache lines, and `chain_ids` is empty for most alleles (chain-id population is sparse). The 2026-05-20 cohort review filed this as L2; this review's PSP reader is one of the hot consumers that would migrate in lockstep.
- **Mechanism:** Migrate `chain_ids: Vec<ChainId>` to `SmallVec<[ChainId; 1]>` (or `Box<[ChainId]>` empty-as-null). PSP reader changes follow downstream of the AlleleObservation refactor; track as a coordinated PR with 2026-05-20 L2.
- **Measurement plan:** Downstream of 2026-05-20 L2 — same paired bench.
- **Complexity cost:** Cross-cutting (touches every reader of `AlleleObservation.chain_ids`).

#### L6: [benches/psp_reader_perf.rs:247](../../../../benches/psp_reader_perf.rs#L247) — Bench uses `Cursor<&[u8]>`; cannot measure H2 (seek-defeats-bufreader)

- **Confidence:** High
- **Hot-path evidence:** Bench feeds `PspReader::new(Cursor::new(bytes))`. `Cursor::seek` is a single `usize` write — no syscall, no BufReader buffer to discard. H2's fix is unmeasurable on the existing bench.
- **Pattern matched (methodology):** "Microbenchmarks lie about cache behavior" — extends to "in-memory `Cursor` lies about syscall + BufReader behavior."
- **Mechanism:** Add a `bench_reader_file_backed` group that writes the SNP / phase / multi fixtures to `tempfile::NamedTempFile` once outside the timed region, then times the reader against `BufReader<File>`. Without this, H2 / L1 / L2 / L3 all have only cohort-e2e as a signal — and the host has ~2× cross-session noise.
- **Measurement plan:** Add the bench group as a prerequisite to H2's measurement (Step 4 of §3).
- **Complexity cost:** ~50 LOC bench code. Reuses the existing fixture builders.

#### L7: missing `examples/dhat_psp_reader.rs`

- **Confidence:** Medium (deferred since 2026-05-13 as L9; the new profile is the justification)
- **Hot-path evidence:** Allocator-leaf attribution shows `decode_bytes_split` at 2099 P-core samples. DHAT directly counts allocation events; the perf sampling profile shows allocator self-time, not allocation count.
- **Mechanism:** Pattern from `examples/dhat_baq.rs` and `examples/dhat_pileup.rs`. Drive `PspReader::new(BufReader::with_capacity(64 KiB, File))` over a real `.psp` file with `dhat::Profiler::new_heap()` running. Compare pre/post-H1: expect per-block allocations to drop from ~N_alleles to 2 per ragged column.
- **Measurement plan:** Prerequisite for H1's Stage A / Stage B decision (Step 5 of §3).
- **Complexity cost:** ~80 LOC example, no new deps (dhat-heap feature exists).

#### L8: [benches/psp_reader_perf.rs](../../../../benches/psp_reader_perf.rs) — Bench bodies count records but never `assert_eq!(count, expected)`

- **Confidence:** High
- **Hot-path evidence:** Same as the 2026-05-23 ref_fetcher review L10 — without a drained-count assertion, a fix that silently truncates the iterator output would benchmark as faster.
- **Mechanism:** Add `assert_eq!(count, N_RECORDS)` after each bench inner loop. Trivial guardrail.
- **Measurement plan:** N/A — correctness guardrail.
- **Complexity cost:** 1 line per bench.

#### L9: [Cargo.toml](../../../../Cargo.toml) (cross-routed from 2026-05-23 ref_fetcher review) — wire `mimalloc::MiMalloc` as `#[global_allocator]` behind `alloc-mimalloc`

- Already tracked under the 2026-05-23 ref-fetcher review's L7. Listed here for completeness — would uniformly reduce allocator cost across PSP decode + the rest of the cohort pipeline.

#### L10: [Cargo.toml](../../../../Cargo.toml) (cross-routed) — add `[profile.release-with-debug]` for reproducible production-binary profiles

- Already tracked under the 2026-05-23 ref-fetcher review's L6. Listed for completeness.

### Speculative

#### S1: [common.rs:46](../../../../src/pop_var_caller/common.rs#L46) — `DEFAULT_BUFFERED_IO_CAPACITY = 64 * 1024` may be undersized vs. per-block read footprint

- **Confidence:** Low
- **Hot-path evidence:** Pattern-match only.
- **Mechanism:** If a typical compressed block exceeds 64 KiB, multiple `read(2)`s fire per block regardless of seek behavior. Sweep `{64 KiB, 256 KiB, 1 MiB, 4 MiB}` after H2 lands. Only meaningful with the buffer-preservation fix in place.
- **Measurement plan:** Gated behind H2.
- **Complexity cost:** One constant change. RSS grows by `N_samples × N_workers × capacity` — 4 MiB × 10 × 13 = ~500 MiB at the high end; check the production target's headroom.

#### S2: [reader.rs:553-554](../../../../src/per_sample_pileup/psp/reader.rs#L553-L554) — `Decompressor` per `RecordsIter` instead of per `PspReader`

- **Confidence:** Low
- **Hot-path evidence:** Pattern-match only — zstd `Decompressor::new()` allocates a DCtx workspace (~tens of KB). If multiple `RecordsIter` instances are constructed per `PspReader` lifetime, hoisting the decompressor onto `PspReader` saves N-1 setups.
- **Mechanism:** `cohort_driver.rs` constructs ONE `RecordsIter` per `PspReader` via `region_records`, so today's structure is fine. Re-evaluate if a workflow ever needs multiple region queries per `PspReader`.
- **Measurement plan:** None — wait for the use case.
- **Complexity cost:** Minor refactor of `RecordsIter` construction.

#### S3: [block.rs:283](../../../../src/per_sample_pileup/psp/block.rs#L283), [block.rs:231](../../../../src/per_sample_pileup/psp/block.rs#L231) — `decode_varint_column` / `decode_scalar_column_pod` allocate an outer `Vec` per call

- **Confidence:** Low (covered by H1's general scratch-buffer pattern)
- **Hot-path evidence:** 10 + 4 allocator-leaf samples respectively. Small compared to `decode_bytes_split`'s 2099.
- **Mechanism:** Generalize the scratch-buffer pattern from H1 — make all `decode_*_column` functions write into a caller-owned `&mut Vec<T>`.
- **Measurement plan:** Gated behind H1 — same scratch-management infrastructure.

#### S4: [reader.rs:1007](../../../../src/per_sample_pileup/psp/reader.rs#L1007) — `DecodedColumn` enum-of-Vecs

- **Confidence:** Low — pattern-match only.
- **Mechanism:** If a future refactor consolidates column storage, an enum-of-`Vec`s vs separate fields-on-DecodedBlock trade-off would matter. Today's flat-fields shape is fine.
- **Measurement plan:** None.

### Note

- **L5 (varint fast/cold split)** confirmed in place at [varint.rs:46-90](../../../../src/per_sample_pileup/psp/varint.rs#L46-L90). `#[cold] #[inline(never)]` on the multi-byte path.
- **L6 (LE-slab cast in `decode_scalar_column_pod`)** confirmed in place at [block.rs:231-264](../../../../src/per_sample_pileup/psp/block.rs#L231-L264). Used by every fixed-width per-allele column.
- **L8 (rustdoc `BufReader::with_capacity(64 * 1024, file)`)** confirmed in place at [reader.rs:107-115](../../../../src/per_sample_pileup/psp/reader.rs#L107-L115); cohort driver follows it at `cohort_driver.rs:493`.
- **L1+L2 from 2026-05-13** (persistent `Decompressor` + scratch buffers on `RecordsIter`) confirmed in place at [reader.rs:527-541](../../../../src/per_sample_pileup/psp/reader.rs#L527-L541).
- `RecordsIter` is **monomorphized** in production (always `PspReader<BufReader<File>>`) — no vtable on the inner reader.
- `PspReader` is `!Send` by design (`reader.rs:490`) — no `Sync`-only-for-the-trait-alias overhead.
- `V1_0_COLUMNS` is a flat compile-time `&[ColumnDef]` constant — no `Lazy` / `OnceLock` ([registry.rs:311](../../../../src/per_sample_pileup/psp/registry.rs#L311)).
- `PspReader::new` opens with 5 sequential seeks (trailer-first then header-first) — cold per file, doesn't merit a refactor.
- The five non-block `SeekFrom::Start` sites in `reader.rs` (lines 135, 149, 181, 200, 291) are all in `PspReader::new` — cold construction path, do not refactor.
- **No concurrency findings.** Per-worker ownership is clean; no `Mutex` / `RwLock` / atomics on the PSP reader hot path.

## 6. Out-of-scope observations

- **zstd at 8.3 % cpu_core combined.** `ZSTD_decompressSequences_bmi2` + `ZSTD_decompressMultiFrame` + `ZSTD_buildFSETable_body_bmi2`. The zstd crate is third-party; tuning would require switching compression level (writer-side decision) or switching algorithms entirely. Not actionable in this review.
- **PSP writer not in scope.** The writer has a separate review history ([perf_psp_writer_2026-05-13.md](perf_psp_writer_2026-05-13.md)). The cohort var-calling path doesn't exercise the writer, so the writer doesn't appear in this profile.
- **AlleleObservation refactor (2026-05-20 L2)** is the cross-cutting downstream change that would let H1's Stage B land cleanly. Not blocked on PSP reader — sequence it independently.

## 7. What's already good

- **L5 + L6 + L1/L2 + L8** from the 2026-05-13 review have stuck in place — no regression-via-rewrite.
- **`RecordsIter` is monomorphized over the inner reader type.** The cohort driver always constructs `PspReader<BufReader<File>>`, so the inner reader's `Read + Seek` calls inline through to `BufReader::read`/`seek` without a vtable. Static dispatch confirmed by reading the bench output — `materialise_next_record` inlines into the `next` symbol.
- **Per-worker ownership invariant** is clean (`PspReader` is `!Send` by design, documented at [reader.rs:490](../../../../src/per_sample_pileup/psp/reader.rs#L490)). Each per-chrom rayon worker constructs its own N readers in `process_one_chromosome` and never crosses thread boundaries — no `Sync`-only-for-trait-alias overhead like the one removed in the ref-fetcher H1.

### Author response convention

Address each finding by its identifier (H1, H2, H3, L1, …, S1, …) with one of:
- `applied in <commit>`
- `experiment shows no gain — closing`
- `disputed because …`
- `deferred to <issue>`
- `won't fix because …`

The "experiment shows no gain" path is expected and welcome — that's what the measurement plan is for. Per the 2026-05-23 ref-fetcher review's H3 outcome (mechanism worked but wall-time gain was in the noise), be ready to revert if H1's Stage A doesn't show ≥ 5 % wall improvement on the cohort fixture.
