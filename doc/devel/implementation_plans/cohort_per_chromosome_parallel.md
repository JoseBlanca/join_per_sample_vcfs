# Cohort `var-calling` — per-chromosome parallelism

Implementation plan for the headline lever from the
[2026-05-20 perf review](../reports/reviews/perf_psp_to_vcf_2026-05-20.md):
parallelise the `.psp` → cohort-VCF pipeline by chromosome.

Today's reality (measured on real tomato `SRR7279725.small.psp × N=10`):

| Threads | Median wall |
|---:|---:|
| 1  | 12.7 s |
| 2  | 11.5 s |
| 4  | 12.0 s |
| 16 | **13.6 s** (slower than T=1) |

The pipeline is a single sequential pull-iterator chain
(`PerPositionMerger → DustFilter → VariantGrouper → PerGroupMerger →
PosteriorEngine → CohortVcfWriter`) — every stage runs on the caller
thread; adding rayon workers buys nothing because there is no
decomposition above `Iterator::next`. The 13-chromosome tomato workload
on a 16-core host suggests the natural decomposition: one rayon worker
per chromosome, end-to-end independent, fragments concatenated in
contig-table order at the writer.

This is **finding H1** in the perf review. Expected payoff: ~6–10×
wall reduction on a 13-chromosome input (theoretical max ~13×; the
~10× imbalance between SL4.0ch00 at 9.6 Mb and SL4.0ch01 at 90 Mb gates
the realistic ceiling).

## Spec / supporting documents

- Perf review (verdict + measurement plan):
  [perf_psp_to_vcf_2026-05-20.md](../reports/reviews/perf_psp_to_vcf_2026-05-20.md).
- Pipeline architecture spec:
  [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md)
  §"Stage 3 — low-complexity filter" through §"Stage 6 — posterior
  engine". Per-chrom independence is asserted in §"Cost and
  parallelism".
- Bench harness (used to verify the win):
  [benches/cohort_e2e_perf.rs](../../benches/cohort_e2e_perf.rs).
- Profiling driver (used during baseline + will be the
  before/after measurement tool):
  [examples/profile_cohort_e2e.rs](../../examples/profile_cohort_e2e.rs).

## Why concat, not sort

The perf review's recommendation was per-chrom workers + concat in
contig-table order. The user proposed a sort-after-write alternative
(spray records to one shared file in any order, sort at the end). We
went with concat for these reasons, recorded here so the choice is
not re-litigated:

- **Cost.** Concat is one streaming byte copy per fragment (linear
  in output size, I/O bound). External-merge sort on a 50 GB VCF is
  `O(n log n)` records, with a measurable serial tail at the end of
  an otherwise-parallel run. The architecture gives us
  sorted-within-chromosome for free; sorting again pays for
  something we already have.
- **No new runtime dependency.** A pure-Rust bgzf-aware concat is
  ~80 lines; pure-Rust VCF sort (external-merge) is several hundred
  + careful bgzf I/O.
- **Cleaner ordering invariant.** Each fragment is internally
  monotonic by construction (the per-chrom pipeline guarantees it
  + the writer's per-contig monotonicity check asserts it). Concat
  walks the contig table in declared order — final VCF is monotonic
  without a reorder buffer.

## Scope

### In scope

- New `process_one_chromosome` helper in
  [src/pop_var_caller/cohort_driver.rs](../../src/pop_var_caller/cohort_driver.rs)
  that runs the existing DUST→…→writer chain over **one chromosome's
  records** and writes to a caller-supplied fragment path.
- New `concat_vcf_fragments` helper in
  `src/var_calling/vcf_writer/concat.rs` (new file) — pure-Rust
  bgzf-aware concat that strips the VCF header from fragments 2..N.
- Re-shape [`run_var_calling`](../../src/pop_var_caller/var_calling.rs#L230)
  to: (1) build per-chromosome reader sets, (2) drive the per-chrom
  workers via `rayon::par_iter`, (3) call `concat_vcf_fragments`,
  (4) atomic-rename `<output>.tmp` → `<output>`.
- Drop the per-group merger's intra-batch `par_iter` at
  [per_group_merger.rs:594-608](../../src/var_calling/per_group_merger.rs#L594)
  (L1 from the perf review) **in the same PR** — nested rayon under
  the outer per-chrom decomposition is wasteful, and the inner
  par_iter is already net-negative on real data (<2 % useful work
  vs. 0.82 % rayon-bridge overhead).
- Tests: end-to-end output equivalence between the new parallel
  path and the existing sequential path on the synthetic cohort
  fixture from
  [tests/cohort_cli_integration.rs](../../tests/cohort_cli_integration.rs);
  edge cases (N_chroms=1, empty chromosomes, worker-error
  propagation).

### Out of scope

- **Streaming concat as workers finish.** Defer to v2. v1 waits for
  all workers to join, then concats serially. The concat step is
  I/O-bound; for v1 it sits at the end of the wall and is small
  relative to the per-worker compute.
- **Block-level bgzf concat with header surgery.** v1 decompresses +
  re-encodes during concat. The decompress + re-encode cost is
  bounded by total output size and is I/O-bound; if it ever shows
  up in a profile as a real bottleneck, the optimisation is "find
  the BGZF block boundary at or after the header's last byte and
  byte-copy the rest with EOF-marker stripping". Out of scope here.
- **Per-chromosome reader sharing.** v1 opens N readers per worker
  × N_chroms workers (130 file descriptors on the tomato N=10
  fixture). At realistic cohort sizes this could climb (N=256 × 13
  chroms = 3328 fds, which exceeds the default `RLIMIT_NOFILE` of
  1024 on most distros). A bump in the binary's startup is
  prudent — call out as a follow-up note, not blocking v1.
- **`var-calling-from-bam`**. The from-bam subcommand uses a
  walker-driven single-sample merger, not a multi-PSP merger. Out
  of scope; it stays serial in this slice.
- **`estimate-contamination`**. Same shape as `var-calling` upstream
  but emits a `ContaminationArtefact` (TOML) not a VCF. Worth a
  follow-up parallelisation pass under the same design, but
  separate slice.
- **Memory-pressure mitigations under N=1000+ cohorts.** The
  fragment-file approach is already mitigation-by-design (workers
  stream records to disk one at a time). Further pressure points
  (the per-chrom merger's `heads: Vec<Option<PileupRecord>>` growing
  with N) are tracked under H6 / L7 in the perf review and are
  independent.

## Design

### 1. Per-chromosome reader splitting

The `.psp` format has per-chromosome block boundaries (every block's
header carries a `chrom_id`; block index is sorted by
`(chrom_id, first_pos)`). The existing
[`PspReader::region_records(chrom_id, start, end)`](../../src/per_sample_pileup/psp/reader.rs#L364)
seeks to the first overlapping block and yields records clamped to
`[start, end]` on the requested chromosome.

For whole-chromosome iteration, the per-chrom worker calls
`reader.region_records(chrom_id, 1, chrom.length)`. Single seek per
reader, no PSP API changes.

**Ownership:** each worker opens its own `PspReader` set (one per
sample). No cross-worker sharing of readers. At N=10 / 13 chroms
that's 130 fds total; at N=256 / 13 chroms = 3328 fds — a future
slice may want to bump `RLIMIT_NOFILE` in the binary's startup, but
v1 only targets the tomato N=10 fixture (130 fds well within the
1024 default).

### 2. Worker helper — `process_one_chromosome`

New helper in [src/pop_var_caller/cohort_driver.rs](../../src/pop_var_caller/cohort_driver.rs):

```rust
pub(crate) fn process_one_chromosome(
    chrom_id: u32,
    psp_paths: &[PathBuf],         // shared with all workers; each opens its own readers
    sample_names: Vec<String>,
    all_chromosomes: Vec<ParsedChromosome>,  // full table — needed for the writer's contig header
    fragment_path: PathBuf,        // where this worker writes
    metadata: CohortMetadata,      // full cohort metadata (same for every worker)
    writer_cfg_template: WriterConfig,  // emit_gp + bgzf-vs-plain inherited
    pipeline_params: CohortPipelineParams,
) -> Result<(u32 /* chrom_id */, u64 /* records_written */), VarCallingCliError>
```

Body:

1. Open one `PspReader` per sample (re-opening the `.psp` files; cheap, single seek per open).
2. Build per-sample iterators via
   `reader.region_records(chrom_id, 1, all_chromosomes[chrom_id as usize].length)`.
3. Build a per-chrom `PerPositionMerger` over those iterators with
   the full `sample_names` + `all_chromosomes`. The merger's per-contig
   monotonicity check passes trivially (only one chrom appears in its
   input).
4. Build `WriterConfig` from the template with `output = fragment_path`.
5. Call the existing
   [`drive_cohort_pipeline`](../../src/pop_var_caller/cohort_driver.rs#L64)
   with the per-chrom merger + the supplied params + the fragment-path
   writer config + the full cohort metadata.
6. Return `(chrom_id, records_written)` on success; bubble the typed
   error otherwise.

Each fragment is a **complete self-contained `.vcf.gz`** (header +
records). No `CohortVcfWriter` API changes.

### 3. Tempdir + fragment paths

In `run_var_calling`, after the existing setup + validation:

```rust
let frags_dir = tempfile::TempDir::new_in(
    args.output.parent().unwrap_or(Path::new(".")),
)
.map_err(VarCallingCliError::Io)?;
```

Per-chrom path:

```rust
let extension = if final_output_is_bgzf(&args.output) { "vcf.gz" } else { "vcf" };
let fragment_path = |chrom_id: u32, name: &str| -> PathBuf {
    frags_dir.path().join(format!("chr_{chrom_id:03}_{name}.{extension}"))
};
```

The `TempDir`'s RAII `Drop` removes the whole subtree on success
**and** on panic. No working-directory pollution. Locating it next to
the output (not `/tmp`) follows the project's `no-system-tmp` rule from
[CLAUDE.md](../../CLAUDE.md) and keeps the final atomic rename within
one filesystem.

### 4. Rayon scope + error handling

```rust
let fragment_paths: Vec<PathBuf> = chromosomes
    .iter()
    .enumerate()
    .map(|(cid, c)| fragment_path(cid as u32, &c.name))
    .collect();

// Parallel per-chrom drive. Each worker owns its readers + writer.
let per_chrom_results: Vec<Result<(u32, u64), VarCallingCliError>> = chromosomes
    .par_iter()
    .enumerate()
    .map(|(cid, _chrom)| {
        process_one_chromosome(
            cid as u32,
            &args.psp_files,
            sample_names.clone(),
            chromosomes.clone(),
            fragment_paths[cid].clone(),
            metadata.clone(),
            writer_cfg.clone(),
            // CohortPipelineParams isn't Clone today — see "Open work" #2 below
            pipeline_params_for_worker(/* ... */),
        )
    })
    .collect();

// Fail-fast on the first error (option (a) from the design discussion).
let mut total_records: u64 = 0;
for result in per_chrom_results {
    let (_cid, n) = result?;
    total_records += n;
}
```

**Error semantics — option (a) chosen.** Rayon's `par_iter().collect()`
continues all workers even if some fail; we collect a
`Vec<Result<...>>` and return the first error after all workers join.
Matches the current sequential driver's typed-error API; simplest
shape. Latency cost on bad input is at most one chromosome's runtime
(the slowest still-running worker), acceptable for a 13-worker pool.

**Rayon thread pool.** Already wired through `configure_rayon_pool`
from `--threads`. The per-chrom `par_iter` automatically uses the
global pool sized by that helper. No new init code.

**`--threads` vs. chromosome count — soft cap, not an error.**
Effective parallelism is `min(threads, num_chromosomes)` because the
per-chrom `par_iter` only enqueues `num_chromosomes` jobs no matter
how big the pool is. Excess pool threads sit idle; no correctness
issue, no memory blowup. We do **not** error on `--threads >
num_chromosomes` because:

- Wrapper scripts running across multiple species
  (`--threads $(nproc)` for tomato + bacterium + human) would have to
  special-case the chromosome count per species; samtools / bcftools
  / GATK accept over-provisioned thread counts and just use what
  they can.
- The constraint is workload-shape-dependent. When finer-grained
  parallelism lands (rayon-over-records within a chromosome;
  sub-chromosome DUST chunking), "more threads than chromosomes"
  becomes useful — the error would have to be relaxed, breaking the
  same wrapper scripts.

Instead, the **run summary surfaces the cap**:

```
var-calling: n_samples=10 records_emitted=126499 \
  effective_threads=13 (requested 16; capped by 13 chromosomes)
```

The parenthetical only appears when the cap actually bit. When no
cap or no `--threads` flag, the line is the existing format plus an
`effective_threads=N` field. Implementation: in
[`print_run_summary`](../../src/pop_var_caller/var_calling.rs#L433),
add an `effective_threads` arg + an `Option<u32>` for the requested
count; compose the parenthetical only when `Some(requested) > n_chroms`.

### 5. Pure-Rust concat — `src/var_calling/vcf_writer/concat.rs`

New module exporting:

```rust
pub fn concat_fragments(
    output_tmp: &Path,        // <output>.tmp
    fragments: &[PathBuf],    // in contig-table order
    output_kind: SinkKind,    // Bgzf or Plain — must match each fragment's extension
) -> Result<(), VcfWriteError>
```

Algorithm:

1. Open `output_tmp` for write:
   - `Plain` → `BufWriter::with_capacity(64 * 1024, File::create(output_tmp)?)`.
   - `Bgzf`  → `bgzf::Writer::new(BufWriter::with_capacity(64 * 1024, File::create(output_tmp)?))`.
2. For each `fragment` in order:
   - Open with `bgzf::Reader::new(File::open(fragment)?)` if `Bgzf`,
     `BufReader::with_capacity(64 * 1024, File::open(fragment)?)` if `Plain`.
   - Walk lines (`BufRead::read_line` into a reused `String`).
   - **First fragment:** write every line.
   - **Fragments 2..N:** skip lines starting with `#` (the VCF
     header — both `##` meta-lines and the single `#CHROM\t...` line).
     Write the rest.
3. `writer.finish()?` (bgzf) or `writer.flush()?` (plain).

This decompresses + re-encodes each fragment (necessary to skip the
header in a bgzf stream — block-level surgery is the v2 optimisation).
Cost is bounded by total output size; I/O-bound in practice.

The `SinkKind` enum already exists or can be inferred from the output
path's extension via the same logic the writer uses
([vcf_writer/sink.rs](../../src/var_calling/vcf_writer/sink.rs)).

### 6. Tying it together in `run_var_calling`

The current sequential body of
[`run_var_calling`](../../src/pop_var_caller/var_calling.rs#L230) does
(steps numbered in the existing source):

1. configure rayon, 2. open readers, 3. ref basename check,
4. sample names + chrom agreement, 5. build per-stage configs,
7. build fetcher + FASTA MD5 verify, 8. cohort metadata + writer
config, 9. wire merger, 10. drive_cohort_pipeline, 11. run summary.

The new shape replaces steps 9-10 with:

```rust
// 9. Build per-chrom fragment paths under a TempDir next to the output.
let frags_dir = tempfile::TempDir::new_in(...)?;
let fragment_paths: Vec<_> = ...;

// 10. Drive per-chrom workers in parallel.
let per_chrom_results: Vec<Result<(u32, u64), _>> = chromosomes
    .par_iter().enumerate()
    .map(|(cid, _)| process_one_chromosome(...))
    .collect();
let mut total_records = 0;
for r in per_chrom_results { total_records += r?.1; }

// 10b. Concat fragments in contig-table order to <output>.tmp.
let output_tmp = tmp_path_for(&args.output);
let kind = SinkKind::from_path(&args.output);
concat_fragments(&output_tmp, &fragment_paths, kind)?;

// 10c. Atomic rename <output>.tmp → <output> + parent dir fsync.
fs::rename(&output_tmp, &args.output).map_err(VarCallingCliError::Io)?;
sync_parent_dir(&args.output).map_err(VarCallingCliError::Io)?;

// 11. Run summary (unchanged).
print_run_summary(&sample_names, total_records);
```

Steps 1-8 stay byte-for-byte the same.

The existing `drive_cohort_pipeline`'s tmp-cleanup-on-error discipline
(remove `<output>.tmp` if the driver fails mid-stream) moves into
`process_one_chromosome` for the per-chrom fragment paths; the
top-level `run_var_calling` only manages the final `<output>.tmp`.

### 7. Per-group merger `par_iter` removal (L1)

[per_group_merger.rs:594-608](../../src/var_calling/per_group_merger.rs#L594)
currently does `batch.into_par_iter().map(process_group).collect()`.
Under per-chrom outer parallelism this is nested rayon — wasteful even
when the inner work is non-trivial (and on real data the inner work is
<1 % of self-time per the perf review). One-line change:

```rust
let results: Vec<Result<Option<MergedRecord>, PerGroupMergerError>> = batch
    .into_iter()
    .map(|group| process_group(group, fetcher.as_ref(), &config, &genotype_tables))
    .collect();
```

Drops two `Arc::clone` lines (`&self.ref_fetcher`, `&self.genotype_tables`).
Drops the `use rayon::prelude::*` if not used elsewhere in the module.

Bundled into the same PR because:
- The two changes are interdependent (outer per-chrom parallelism +
  inner serial = correct; outer + inner par_iter = nested rayon).
- The test plan for H1 needs to measure parallel speedup; leaving the
  inner par_iter in place pollutes the measurement.

## Test plan

Integration tests in
[tests/cohort_cli_integration.rs](../../tests/cohort_cli_integration.rs):

1. **`var_calling_parallel_emits_same_records_as_serial`** —
   golden-comparison test. Run the existing `var_calling_happy_path_three_samples`
   fixture, dump both outputs to plain VCF, diff record-level (ignore
   header timestamps + `##source` lines). Output must be byte-identical
   in the body.
2. **`var_calling_one_chromosome_input_works`** — degenerate case:
   single-contig fixture, one worker, concat with one fragment ≡ identity.
3. **`var_calling_empty_chromosome_in_middle`** — three-contig fixture
   where the middle contig has zero records (worker produces a
   header-only fragment). Concat must handle the empty-body case.
4. **`var_calling_worker_error_surfaces_first_error`** — fixture
   crafted to make exactly one chromosome's worker fail (e.g. corrupt
   one PSP file's chr1 block); assert the typed error is the one we
   would have got from the serial driver.

Unit tests in `src/var_calling/vcf_writer/concat.rs`:

5. **`concat_plain_two_fragments_strips_second_header`** — synthesize
   two plain-text VCF fragments, concat, assert exactly one header
   block + both bodies in order.
6. **`concat_bgzf_two_fragments_round_trip`** — same with bgzf
   fragments; decompress the concat output and assert exactly one
   header + bodies.
7. **`concat_single_fragment_is_identity`** — N=1 fragment, output
   byte-identical to input.
8. **`concat_empty_fragments_list_writes_only_header`** — degenerate
   guard; should it be a hard error? Decision: yes, return
   `VcfWriteError::EmptyFragmentList` (new variant); a downstream
   pipeline producing zero fragments is a bug, not a quiet success.

Bench validation (not new benches, just re-runs):

- `cohort_e2e_full/scaling_samples` at N=10, 64, 256 before/after —
  expect comparable wall at N=10 (overhead of per-chrom setup), much
  faster at N=64+ where the parallel speedup amortises.
- `cohort_e2e_full/scaling_region` at L=1k, 5k, 20k before/after.
- The one-off `examples/profile_cohort_e2e.rs` script at T=1, 2, 4, 8,
  13, 16 against real tomato. **Acceptance threshold: T=8 wall ≥ 4×
  reduction vs. T=1** (justifies the per-chrom plumbing). Below 2×
  the design needs re-thinking before merging.

## Sequencing

A single PR. Split into reviewable commits in this order:

1. **Bench-side L1**: per-group merger `par_iter` → `into_iter`. Lands
   alone first so its perf delta is visible in isolation.
2. **Concat helper**: `src/var_calling/vcf_writer/concat.rs` + unit
   tests. Standalone, no callers yet.
3. **Worker helper**: `process_one_chromosome` in `cohort_driver.rs`.
   Unit-tested via the existing `drive_cohort_pipeline` tests — the
   new function is the same shape, just chrom-bounded.
4. **`run_var_calling` reshape**: wire everything together. Integration
   tests added here.
5. **Bench validation**: re-run `cohort_e2e_full` + the
   `profile_cohort_e2e` thread sweep on real data; record numbers in
   the impl report.

## Open work / non-goals

1. **`CohortPipelineParams` not `Clone`.** Today
   [`CohortPipelineParams`](../../src/pop_var_caller/cohort_driver.rs#L32)
   is constructed once and consumed by-value in
   `drive_cohort_pipeline`. For 13 parallel workers each needs its own
   instance. Two options:
   - (a) Add `#[derive(Clone)]` if every field is `Clone` (the configs
     are; `fetcher: SharedRefFetcher` is `Arc` so cheap; `chromosomes`
     is `Vec<ParsedChromosome>` which is `Clone`). Cheapest path.
   - (b) Build a `pipeline_params_for_worker(&shared_inputs) -> CohortPipelineParams`
     factory that takes the inputs by reference and builds fresh for
     each worker.

   Decision: (a) — derive `Clone`. The fetcher's `Arc::clone` is one
   atomic per worker, negligible.

2. **`SyncRefFetcher` contention under parallel workers (L5).** Not
   addressed in this slice; the perf review's L5 finding notes the
   `noodles::Repository::RwLock<HashMap>::read()` on every fetch will
   contend under 16 workers. Today (single-threaded) it's silent. The
   landed parallel implementation should be re-profiled; if `RwLock`
   contention is visible (look for `_raw_spin_lock` jumping in the
   per-worker call graph), schedule L5 as the next slice.

3. **`RLIMIT_NOFILE` bump.** At realistic cohort sizes (N=256 × 13
   chroms = 3328 fds) we may exceed default ulimit. Track as a note;
   the binary's startup could call `rlimit::increase_nofile_limit(8192)`
   defensively, but v1's tomato fixture is well within the default.

4. **Streaming concat (v2).** As soon as the chr_0 worker finishes,
   the main thread could start appending it to `<output>.tmp` while
   the slower chromosomes are still processing. End-to-end wall
   becomes `max(slowest_chrom, total_concat_time)`. Out of scope here.

5. **Block-level bgzf concat (v2).** Identify each fragment's
   header-end BGZF block boundary; byte-copy the body blocks raw with
   EOF-marker stripping. Avoids decompress + re-encode. Out of scope
   for v1; the v1 decompress-re-encode is bounded by total output
   bytes and is I/O-bound in practice.

6. **`var-calling-from-bam` parallelisation.** The from-bam subcommand
   runs the pileup walker for one sample's BAM(s) and feeds records
   directly into the cohort pipeline. The walker is per-sample (only
   one sample exists), so per-chromosome parallelism there is a
   different shape (walker-over-chrom + per-chrom var-calling). Not
   covered here; track as a separate slice if profiling indicates the
   single-sample path is worth parallelising.

7. **`estimate-contamination` parallelisation.** Same upstream shape
   as `var-calling` (cohort of `.psp`s → per-position records). Output
   is a contamination TOML, not a VCF. The per-chrom decomposition
   applies but the reduction step is different (one TOML, not many VCF
   fragments). Defer to a separate slice.

## Forward references

- Once H1 lands, **L1** (per-group `par_iter` removal) is included.
- **L5** (`SyncRefFetcher` pre-warm + drop `RwLock`) becomes the next
  ceiling under parallel workers; address as a separate slice if
  re-profiling shows it.
- **H7** (DUST criterion bench) is the gate for any subsequent DUST
  inner-loop tuning (H2/H3); independent of this slice.
- **H6** (PerPositionMerger `vec![None; n_samples]` per emit) is per-record
  and runs inside every per-chrom worker; the parallel speedup
  multiplies it but does not change the per-record cost. Address as a
  separate slice.

## Estimated effort

- Concat helper + unit tests: ~150 lines.
- Worker helper + reshape of `run_var_calling`: ~80 lines net diff.
- Per-group merger `par_iter` → `into_iter`: ~5 lines net diff.
- Integration tests: ~150 lines.

Total: ~400 net lines. The biggest risk surface is the concat helper
(bgzf header-skip correctness on real fragments); cover with both
synthetic unit tests and the real-tomato integration smoke before the
impl report is written.
