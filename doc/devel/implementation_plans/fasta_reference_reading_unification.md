# FASTA reference reading — unification + `--regions` perf fix

Design + implementation plan to remove the per-region FASTA-reload
regression on the Stage 1 pileup path and unify the lifecycle of the
reference readers it uses.

Branch: `fasta-ref-reading-perf`.

## 1. The reported problem

A colleague observed that `pop_var_caller pileup` is much slower with
`--regions <bed>` than without it, even though both are meant to do the
same work. Reproduced on the `human_genome_bottle` fixture (HG002 CRAM
pre-restricted to a 1000-region BED, `--no-baq --threads 4`):

| Run | Wall | Peak RSS | `.psp` |
|---|---:|---:|---:|
| **with** `--regions` (1000 regions ≈ 5 Mb) | 253 s | **14.6 GB** | 95 MB |
| **without** `--regions` (whole genome) | 180 s | **0.78 GB** | 102 MB |

The `--regions` run produces *less* output (it sequences a strict subset
of the genome) yet is slower **and uses 18× the memory**. `sys` time
jumps 1.1 s → 17 s (a syscall storm). Peak RSS is already 5.3 GB at just
**20 regions** — ~265 MB per region, the size of a whole human
chromosome. That number is the tell.

## 2. Root cause (grounded in the code)

Both production paths run the **same loop** in
[`run_pileup`](../../src/pop_var_caller/cli.rs#L312): `pileup` only ever
calls [`AlignmentMergedReader::query`](../../src/bam/alignment_input.rs#L845)
(`new()` is test/diagnostics-only — see its call sites). The *only*
difference between the two runs is the region argument:

- **without** `--regions`: `RegionSet::whole_contigs` yields one
  whole-contig span per contig; the span is dropped to `None` at
  [alignment_input.rs:881](../../src/bam/alignment_input.rs#L881), so
  `query()` runs **once per contig**.
- **with** `--regions`: `query()` runs **once per BED region** —
  e.g. ~40 times for chr1.

Every `query()` call builds a **fresh** noodles `fasta::Repository`
([alignment_input.rs:893-900](../../src/bam/alignment_input.rs#L893)).
That repository is the CRAM reference resolver (CRAM is
reference-compressed, so decoding needs reference bases) and it also
backs the reader's F1 read-mismatch filter
([alignment_input.rs:1267](../../src/bam/alignment_input.rs#L1267)).

The noodles `Repository` (`noodles-fasta-0.61.0/src/repository.rs`) is a
**whole-contig cache keyed by name**: `get(name)` calls its
`IndexedReader` adapter, whose `get` issues
`reader.query(Region::new(name, ..))` — an **unbounded** range, i.e.
read the entire contig — and caches the result as an `Arc<Sequence>`
(`record.sequence().clone()`). The CRAM slice decoder pulls reference
bases through it per slice (`get_slice_reference_sequence` in
`noodles-cram-0.93.0/.../container/slice.rs`).

So **touching any base of chr1 loads all ~250 MB of chr1.** With
`--regions`, chr1 is loaded into ~40 short-lived per-region
repositories; the system allocator retains the churned buffers (RSS
plateaus ~15 GB). Without `--regions`, chr1 is loaded exactly once (the
single whole-contig `query()`), reused across all its slices, then
dropped before the next contig → one contig resident → 0.78 GB.

Arithmetic check: 20 regions × ~250 MB (chr1) ≈ 5 GB ≈ the measured
5.3 GB. The memory blowup *is* the repeated whole-contig load. The wall
penalty is the same redundant work: re-reading + re-decoding hundreds of
MB of reference per region.

A second, smaller per-region tax: the walker's reference fetcher
([`MultiChromStreamingRefFetcher`](../../src/fasta/fetcher.rs#L824)) is
also rebuilt per region in
[`with_stage1_chain`](../../src/pop_var_caller/stage1_pipeline.rs#L124),
re-parsing the 2,580-line `.fai` and re-opening the FASTA each time. It
*streams* (1 MiB sliding buffer), so it doesn't drive the memory blowup,
but it contributes to the `sys`-time storm.

## 3. Current architecture: two FASTA-reader families

There are **two independent reference-reading mechanisms** on this path,
each reading the same FASTA:

| # | Reader | Consumer | Strategy | Rebuilt |
|---|---|---|---|---|
| 1 | noodles `fasta::Repository` | CRAM decode + reader F1 mismatch filter | whole-contig cache by name (`Arc<Sequence>`) | **per region** (in `query`) |
| 2 | [`MultiChromStreamingRefFetcher`](../../src/fasta/fetcher.rs#L824) (wraps [`StreamingChromRefFetcher`](../../src/fasta/fetcher.rs#L82)) | Stage 1 walker | 1 MiB streaming sliding buffer, contig-keyed swap | **per region** (in `with_stage1_chain`) |
| 2′ | [`ManualEvictChromRefFetcher`](../../src/fasta/fetcher.rs#L988) | Stage 1 BAQ (per rayon worker) | contig buffer, caller-managed evict | **per chunk** (in `BaqStream`) |

Family 2 / 2′ (the `ChromRefFetcher` trait family) was already unified
by the [`unified_chrom_ref_fetcher`](unified_chrom_ref_fetcher.md) plan.
That plan **explicitly left the noodles `Repository` untouched** — it is
reader #1, a separate type the CRAM decoder mandates. Reader #1 is the
one that regresses.

### Why two readers, and why we keep both

They answer to two different masters:

- **noodles `fasta::Repository`** is mandated by the noodles **CRAM
  decoder API**: CRAM is reference-compressed, so reconstructing a read's
  `SEQ` needs the reference at the read's span, and noodles pulls it
  through a `Repository` (`slice.records(repository, …)`). We keep using
  noodles to read CRAM — these are complex files and reimplementing
  their reference resolution is not worth it.
- **`StreamingChromRefFetcher`** is *our* low-memory reader for the
  walker / cohort / BAQ, built for the memory-efficiency thesis (the
  RAM-sensitive cohort path). noodles' `Repository` can't provide a
  streaming, sub-range read, which is why this family exists at all.

### Do we really need the whole contig in memory for one region?

No, not fundamentally — the analysis only needs the region's span. The
whole-contig load is an artifact of the **`Repository` abstraction**, not
a real requirement: `Repository::get(name)` (and the underlying
`Adapter::get`) are *whole-sequence* operations — the `IndexedReader`
adapter reads `Region::new(name, ..)` (the unbounded contig) and caches
it as one `Arc<Sequence>`, which the CRAM decoder then indexes by
**absolute 1-based contig position**. You cannot substitute a sub-range
slice without corrupting those coordinates. Our own
`StreamingChromRefFetcher` proves sub-range reading is possible, but the
CRAM decoder can only consume a whole-contig `Repository`.

**Consequence:** for CRAM decode, whole-contig residency is effectively
the floor. The without-`--regions` run already sits at that floor (one
contig resident → 0.78 GB). So Phase 1 does **not** stop loading the
contig — it stops *reloading* it: load once per contig, reuse across that
contig's regions, evict on contig change, reaching the same floor the
whole-genome path already achieves. Making the reference truly
region-bounded for CRAM would require bypassing noodles' decode API and
is out of scope.

### The redundancy this leaves (Phase 2 territory)

For CRAM the `Repository` has already loaded the whole current contig
into RAM, yet the walker (and BAQ) re-read the same bytes through their
*own* streaming readers. That double-read is the genuine architectural
redundancy; Phase 2 collapses it (serve the walker/F1/BAQ from the
Repository's resident `Arc<Sequence>` for CRAM; stream for BAM).

### The remaining smell is lifecycle, not type

Every reference reader on this path is constructed *per region* (or per
chunk) when the access pattern is *per contig*. The region loop already
visits regions grouped by contig — `RegionSet` is sorted by
`(chrom_id, start)` and non-overlapping
([regions.rs:67](../../src/regions.rs#L67), enforced by `sort_and_merge`
and `whole_contigs`). So each contig is a single contiguous ascending
run, never revisited. Nothing about the workload requires per-region
construction; it is an artifact of `query()` being self-contained.

## 4. Proposed design: build-once, contig-scoped

The unifying principle: **reference reading is built once per run and
scoped to one resident contig at a time, advancing as the region loop
crosses contig boundaries.** This is exactly the memory/IO profile the
without-`--regions` path already achieves incidentally; we make the
`--regions` path do it deliberately, and give all readers the same
lifecycle.

Two phases. Phase 1 fixes the regression and unifies the *lifecycle*.
Phase 2 (optional, later) unifies the *reads* themselves for CRAM.

### Phase 1 — hoist construction out of the region loop

**1a. Share one `fasta::Repository` across all regions; bound it
per-contig.**

- `run_pileup` builds the `Repository` **once** (the code currently in
  `query()` lines 887-900 moves up; or `load_pileup_inputs` builds it
  and stores it in `PileupInputs`).
- `query()` gains a `repository: &fasta::Repository` parameter and uses
  the passed-in repository instead of building its own. Its
  `from_open_alignment_files(..., Some(repository.clone()))` tail is
  unchanged (the clone is an `Arc` bump).
- The region loop tracks the current `chrom_id`. When it changes, call
  `repository.clear()` **before** opening the first region of the new
  contig. `clear()` drops the cache `HashMap` (freeing the previous
  contig's `Arc<Sequence>`) but keeps the adapter's open file handle —
  no re-open, no `.fai` re-parse.

Result: each contig's sequence is loaded **once** and reused across all
its regions; only **one** contig is resident at a time. This is safe
because `RegionSet` is contig-grouped and ascending — once we leave a
contig we never return to it, so clearing its cache can't cause a
reload.

Expected: peak RSS drops from 14.6 GB → ≈ the 0.78 GB whole-genome
profile; wall drops to ≤ the whole-genome number (strictly less work).

**1b. Hoist the walker fetcher.**

This does **not** touch the whole-contig memory (the walker fetcher
streams and never loads a whole contig). It removes only the per-region
`.fai` re-parse + file re-open — a `sys`-time tax, not the RSS blowup.

- Build `MultiChromStreamingRefFetcher` **once** in `run_pileup`; pass
  `&MultiChromStreamingRefFetcher` into `with_stage1_chain` instead of
  constructing it per region
  ([stage1_pipeline.rs:124](../../src/pop_var_caller/stage1_pipeline.rs#L124)).
- It already rebuilds its inner `StreamingChromRefFetcher` on contig
  transition, so within a contig the inner streamer (and its 1 MiB
  buffer) is reused across regions, and the `.fai` is parsed once per
  contig instead of once per region.
- **Monotonicity holds across regions.** Regions on one contig are
  sorted and non-overlapping, so region *k+1*'s start exceeds region
  *k*'s end; the streamer's monotonic-forward `fetch` contract is
  satisfied across the region boundary (it simply refills forward). The
  walker uses `fetch` only (not `iter_bases`), so there is no
  phase-reset interaction to worry about.

### Phase 2 (optional, later) — one resident contig, read once

For CRAM input the `Repository` already holds the whole current contig
as an `Arc<Sequence>`. The walker (and the F1 filter, already on the
repository) re-read the same bytes through a *separate* streaming
reader. A deeper unification would expose one **`ContigReference`**
abstraction, built once per contig, that:

- is sourced from the `Repository`'s `Arc<Sequence>` for CRAM (zero
  extra bytes read), and from a streamed/loaded buffer for BAM (which
  has no repository);
- serves the walker (forward windows), the F1 filter (per-read slices),
  and BAQ (random-access-within-chunk views) from that one resident
  copy.

This removes reads #2/#2′ entirely for CRAM. It is deferred because it
touches the byte-identity-critical walker reference path and the BAM
path's memory profile (loading a whole contig for BAM is a memory
regression vs today's streaming; needs a per-input-format decision).
Phase 1 already eliminates the regression, so Phase 2 is a cleanliness
/ micro-optimization follow-up, not a fix.

BAQ's per-chunk `ManualEvictChromRefFetcher` rebuild (the `.fai`
re-parse per chunk) is folded into Phase 2 as well — or addressed
independently by handing BAQ a pre-parsed `.fai` / fetcher factory if a
BAQ-on profile shows it matters.

## 5. API changes (Phase 1)

- `AlignmentMergedReader::query(...)` — add a
  `repository: &fasta::Repository` parameter; delete the in-body
  `Repository` construction (lines 887-900). The `MissingFastaIndex`
  pre-check moves to the single build site.
- `run_pileup` — build the `Repository` once (or read it from
  `PileupInputs`); track `current_chrom_id` and `repository.clear()` on
  transition; build the `MultiChromStreamingRefFetcher` once.
- `with_stage1_chain(...)` — replace the `reference: &Path` →
  per-call-`MultiChromStreamingRefFetcher::new` with an injected
  `walker_fetcher: &MultiChromStreamingRefFetcher`. (The `reference`
  path is still needed by `BaqStream`; keep both, or build BaqStream's
  fetcher factory from the shared handle.)
- `new()` is left as-is (test/diagnostics; not in the hot loop).

## 6. Migration steps (each shippable, suite green)

```
Step 1 — 1a: share the Repository + clear-on-contig.
   query() takes &Repository; run_pileup builds it once and clears
   per contig. Smallest diff that kills the 18× memory + the reload.
        │
Step 2 — 1b: hoist the walker MultiChromStreamingRefFetcher to run scope.
   Removes the per-region .fai re-parse for the walker.
        │
Step 3 — (optional) Phase 2 ContigReference unification + BAQ.
```

Each step is independently measurable on `human_genome_bottle`.

## 6.1 Results (Phase 1, implemented)

`human_genome_bottle`, HG002 CRAM, `--no-baq --threads 4`, 1000-region
BED:

| Build | Wall | Peak RSS | `sys` |
|---|---:|---:|---:|
| baseline `--regions` (before) | 253 s | 14.6 GB | 17.1 s |
| **after Phase 1a+1b** | **202 s** | **0.75 GB** | 1.3 s |
| reference: without `--regions` (whole genome) | 180 s | 0.78 GB | 1.1 s |

- **Memory: −95 % (14.6 GB → 0.75 GB)** — `--regions` now sits at the
  one-contig floor the whole-genome path always had. This is the
  regression.
- **Wall: −20 % (253 s → 202 s)**, `sys` 17.1 s → 1.3 s. The remaining
  ~20 s over the whole-genome baseline is the genuinely-per-region CRAM
  indexed-query setup (one `cram::io::Reader<File>` open + seek per
  region); sharing CRAM readers per contig is a possible later
  optimization, not part of this regression fix.
- **Byte-identity:** the decoded pileup (`psp-to-pileup`) is identical
  before/after — same line count (5 050 273) and SHA-256. The raw
  `.psp` differs only in the header (the `created` timestamp and the
  recorded output path), which shifts block offsets; the cache changes
  are read-timing-only, never which bytes are returned.
- **1b alone** moved wall negligibly (the per-region `.fai` re-parse was
  small next to the CRAM query setup); it is kept for the consistent
  build-once lifecycle and to remove the redundant per-region work.

Code: [`build_fasta_repository`](../../src/bam/alignment_input.rs) +
shared-repository `query` param; `run_pileup`'s clear-on-contig loop and
hoisted `MultiChromStreamingRefFetcher`
([cli.rs](../../src/pop_var_caller/cli.rs)); injected `walker_fetcher`
in [stage1_pipeline.rs](../../src/pop_var_caller/stage1_pipeline.rs).

## 7. Correctness / byte-identity

- **Output must be byte-identical** to the current `--regions` run
  (and ideally to the without-`--regions` run on the same regions).
  Validate by diffing the `.psp` before/after Step 1 and Step 2 on the
  same BED. The repository is a *cache* — sharing or clearing it changes
  only when bytes are read, never which bytes are returned. The walker
  streamer returns identical bytes whether built per-region or per-run
  (same `.fai`, same uppercasing). So both steps are pure
  lifecycle/cache changes with no algorithmic effect.
- **Memory bound depends on the sorted-`RegionSet` invariant.** If that
  invariant ever weakened (interleaved contigs), clear-on-transition
  would thrash (reload), not miscompute — it degrades to today's
  behavior rather than breaking. The invariant is enforced in
  `regions.rs`; add a guard/assertion at the clear site documenting the
  dependency.

## 8. Risks / open items

1. **`query()` callers.** Several tests call `query()` directly
   (alignment_input.rs test module). Each needs a `&repository`
   argument — mechanical, but ~8 call sites.
2. **Repository ownership in `PileupInputs`.** Deciding whether the
   repository lives in `PileupInputs` (built in `load_pileup_inputs`) or
   as a local in `run_pileup`. The former centralizes the FASTA-repo
   construction the 2026-06-03 review flagged as duplicated ×3 (M4) —
   worth folding in.
3. **Phase 2 BAM memory.** Loading a whole contig for the BAM walker
   would regress the streaming memory profile; Phase 2 must keep BAM on
   a streamed source. Flagged for that phase, not Phase 1.
4. **BAQ-on path** is not in the measured bench (`--no-baq`); confirm
   the per-chunk fetcher rebuild isn't a separate regression with a
   BAQ-on profile before deciding Phase 2 priority.

## 9. Validation plan

For each step, on `human_genome_bottle` (`--no-baq --threads 4`,
1000-region BED):

1. `/usr/bin/time -l` wall + peak RSS, vs the 253 s / 14.6 GB baseline
   and the 180 s / 0.78 GB whole-genome target.
2. `cmp` the produced `.psp` against the pre-change `--regions` `.psp`
   (byte-identity).
3. `cargo test --lib --tests` green; `cargo clippy -D warnings`,
   `cargo fmt --check` clean.
