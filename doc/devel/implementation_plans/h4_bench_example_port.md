# H4 — Port broken benches/examples off the retired `RefSeqFetcher` trait

## Context

The Step-2 `ChromRefFetcher` migration plus the M12 (2026-05-23 code
review) legacy-trait removal deleted `RefSeqFetcher` from
`src/per_sample_pileup/pileup/`. Five files plus a sixth that's been
broken since the cohort refactor still reference the old trait or the
pre-refactor `CohortPipelineParams` / `drive_cohort_pipeline` shape.
The new `./scripts/precommit-check.sh` script surfaces them as the
gating compile failure between us and a working pre-commit hook.

PROJECT_STATUS labels this task **"Mechanical port; no new deps."**
This plan honors that — no new crates, no library surface changes.

## Scope (six files)

| File | Failure |
|---|---|
| `examples/dhat_pileup.rs` | E0432 unresolved `pileup::RefSeqFetcher` (line 27) |
| `examples/profile_posterior_engine.rs` | E0432 unresolved `pileup::RefSeqFetcher` (line 25) |
| `benches/pileup_walker_scaling.rs` | E0432 unresolved `pileup::RefSeqFetcher` (line 22) |
| `benches/baq_perf.rs` | E0432 unresolved import + `BaqEngine::process` / `BaqStream::new` signature drift |
| `benches/var_calling_perf.rs` | E0432 unresolved `pileup::RefSeqFetcher` (line 35) |
| `benches/cohort_e2e_perf.rs` | E0432 `SyncRefFetcher` (line 90); E0560 `CohortPipelineParams.fetcher`/`.chromosomes` (lines 272-273); E0061 `drive_cohort_pipeline` 7-vs-5 args (lines 616, 621) |

No other files in `src/`, `benches/`, `examples/`, or `tests/` reference
the removed surfaces.

## Trait situation — splits in two

Sealed-trait pattern recap (the user is new to Rust): a `pub trait` whose
supertrait is `pub(crate)` cannot be implemented from outside the crate,
because external code can't name (let alone satisfy) the supertrait. The
project uses this to keep the cross-call contract (`iter_bases` Drop-reset,
monotonic-forward `fetch`, etc.) under in-crate control.

- **Walker side** — `MultiChromRefFetcher` at [src/per_sample_pileup/pileup/mod.rs:515](../../src/per_sample_pileup/pileup/mod.rs#L515)
  is **not sealed**. The inline `ConstFasta`-style bench types stay; only
  the trait name and the `fetch` return-type change
  (`Result<Vec<u8>, ChromRefFetchError>` instead of
  `Result<Vec<u8>, io::Error>`).

- **Cohort side** — `ChromRefFetcher` at [src/per_sample_pileup/ref_fetcher.rs:359](../../src/per_sample_pileup/ref_fetcher.rs#L359)
  **is** sealed (M9 was an explicit design call). Inline `InMemRef` is
  impossible. The replacement is the existing on-disk path:
  `StreamingChromRefFetcher::for_contig(fasta_path, "chr0")` wrapped in
  `Arc<dyn ChromRefFetcher + Send>` (the `SharedRefFetcher` alias at
  [src/var_calling/per_group_merger.rs:489](../../src/var_calling/per_group_merger.rs#L489)).

- **BAQ side** — worst. [src/per_sample_pileup/baq/engine.rs:108](../../src/per_sample_pileup/baq/engine.rs#L108)
  `BaqEngine::process` takes a *concrete* `&mut ManualEvictChromRefFetcher`,
  not a trait, so the bench's `ConstRefFetcher` cannot satisfy it however
  defined. [src/per_sample_pileup/baq/stream.rs:124](../../src/per_sample_pileup/baq/stream.rs#L124)
  `BaqStream::new` now takes `(reads, cfg, fasta_path, contigs, chunk_size)`
  and builds its own per-worker fetcher inside `map_init` at stream.rs:234.

## Decisions baked in (Jose, 2026-05-24)

- **D1.** Accept the measurement-validity shift on the cohort benches.
  `StreamingChromRefFetcher` adds `RefCell::borrow_mut()` per fetch,
  page-cache effects, and a canonicalise byte loop. Strictest reading of
  "mechanical port, no new deps". Old absolute numbers stop being
  comparable; re-baseline after H4 lands.
- **D2.** For `baq_perf`'s engine bench, construct
  `ManualEvictChromRefFetcher` per-iter inside the `iter_batched` setup
  closure. Untimed. Matches the production constructor at
  `baq/stream.rs:234` more closely than `ConstRefFetcher` ever did.
- **D3.** All new bench/example FASTAs go under project-local `tmp/` via
  `tempfile::TempDir::new_in("tmp")`. Retrofit
  [benches/cohort_e2e_perf.rs:185](../../benches/cohort_e2e_perf.rs#L185)
  for consistency. `/tmp/` is already in `.gitignore`. Each fixture
  creates the `tmp/` dir if missing (`std::fs::create_dir_all("tmp")`).

## Ordered implementation steps

Run after each step inside the container:
`./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`.
Final gate: `./scripts/precommit-check.sh`.

### Step 0 — `tmp/` discipline

Ensure project-local `tmp/` exists at bench-fixture setup time. Add a
shared bench helper (or inline `std::fs::create_dir_all("tmp").expect(..)`
at the top of each ported fixture builder). One commit, no behaviour
change yet — sets the convention before Step 1.

### Step 1 — Walker trait rename

`examples/dhat_pileup.rs` and `benches/pileup_walker_scaling.rs`:
- Replace `use pop_var_caller::per_sample_pileup::pileup::RefSeqFetcher`
  with `MultiChromRefFetcher` and `ChromRefFetchError` (both re-exported
  via `pileup::mod`).
- Change `impl RefSeqFetcher for ConstFasta` → `impl MultiChromRefFetcher
  for ConstFasta`.
- Return `Result<Vec<u8>, ChromRefFetchError>`. The
  `io::Error::new(io::ErrorKind::UnexpectedEof, "off end")` path becomes
  `ChromRefFetchError::Io { chrom_name: "<bench-const>".to_string(),
  source: io::Error::new(io::ErrorKind::UnexpectedEof, "off end") }`
  (variant at [ref_fetcher.rs:316](../../src/per_sample_pileup/ref_fetcher.rs#L316)).
- Update the bench's doc comment to mention `MultiChromRefFetcher`.

Call sites (`run(reads, &ref_fetcher, ...)`) need no change — `run` is
generic over `F: MultiChromRefFetcher`.

### Step 2 — `var_calling_perf.rs` cohort rewrite

- Drop `InMemRef` struct + impl (lines 60-81) and the `RefSeqFetcher`
  import (line 35).
- Add `fn write_fasta(dir: &Path, ref_seq: &[u8]) -> PathBuf` helper.
  Writes 99 N's + `ref_seq` + newline + matching `.fai`. The 99 N's
  preserve the existing `base_offset = 100` offset math without shifting
  the synthetic groups (which would change *what's measured*). Model:
  `cohort_e2e_perf.rs:396-415`.
- Replace `fn shared_fetcher(seq, base_offset)` (lines 83-85) with
  `fn build_shared_fetcher(dir: &Path, ref_seq: &[u8]) -> SharedRefFetcher`:
  writes the FASTA then
  `Arc::new(StreamingChromRefFetcher::for_contig(&fasta_path, "chr0").expect(..))`.
- Each bench fn (`bench_per_group_merger`, `bench_posterior_engine`)
  opens its own `let dir = TempDir::new_in("tmp").expect("tempdir");` at
  the top so the FASTA outlives criterion measurements. Hold all
  `TempDir`s in a top-of-fn `Vec` if multiple fixtures need different
  FASTAs.
- Add `tempfile::TempDir` to imports.

### Step 3 — `profile_posterior_engine.rs` cohort rewrite

Same shape as Step 2. Replace the `InMemRef` block (lines 45-64) and the
`Arc::new(InMemRef { .. })` at line 184. `tempfile` is already a
dev-dep ([Cargo.toml:91-95](../../Cargo.toml#L91-L95)).

### Step 4 — `baq_perf.rs` reshape

- Drop `ConstRefFetcher` (lines 33-46) and the `RefSeqFetcher` import
  (line 28).
- Add `write_fasta` helper (one contig "chr0" of `vec![b'A'; span + 200]`).
- In `bench_engine_read_length`:
  - Build FASTA + `ContigList` once at fn top (untimed).
  - Change the `iter_batched` setup closure to construct
    `ManualEvictChromRefFetcher::for_contig(&fasta_path, "chr0")` per-iter
    alongside the engine (D2 decision). The fetcher is `!Sync` and
    `process` takes `&mut`, so each iter needs its own.
  - Update `process_all` signature to
    `fn process_all(engine: &mut BaqEngine, reads: Vec<MappedRead>,
    fetcher: &mut ManualEvictChromRefFetcher) -> u64`.
- In `bench_stream_chunk_size`:
  - Drop the `&F: RefSeqFetcher + Sync` generic on `drain_stream`.
  - New signature:
    `fn drain_stream(inputs, fasta_path: PathBuf, contigs: ContigList,
    chunk_size) -> u64`.
  - Call `BaqStream::new(inputs.into_iter(), BaqConfig::default(),
    fasta_path, contigs, chunk_size)`.
  - Setup closure clones `fasta_path` and `contigs` per iter
    (`ContigList` is one small `ContigEntry`; clone cost negligible
    against the 50 K × 30× workload).

### Step 5 — `cohort_e2e_perf.rs` drift fixes

Six discrete edits:

1. **Line 90:** `SyncRefFetcher` import → `StreamingChromRefFetcher`.
2. **Lines 165-167, 175-181:** update `CohortFixture` doc comment to
   mention `StreamingChromRefFetcher`.
3. **Lines 219-229:** drop the `contig_list` construction; drop the
   `ContigEntry`/`ContigList` import at line 81 (now unused); drop the
   `md5_hex_to_bytes` + `hex_digit` helpers at lines 376-394 (dead once
   `contig_list` is gone — `-D warnings` will reject otherwise). Replace
   lines 227-229 with:
   ```
   let fetcher_concrete = StreamingChromRefFetcher::for_contig(&fasta, &contig_name)
       .expect("StreamingChromRefFetcher::for_contig");
   let fetcher: SharedRefFetcher = Arc::new(fetcher_concrete);
   ```
4. **Lines 260-278:** delete `fetcher:` and `chromosomes:` from the
   `CohortPipelineParams` literal. The struct's current fields are the
   9 enumerated at [cohort_driver.rs:137-167](../../src/pop_var_caller/cohort_driver.rs#L137-L167).
5. **Lines 616-619, 621-624:** add `0u32` chrom_id and
   `fixture.fetcher.clone()` to both `drive_cohort_pipeline` call sites.
   New shape: `drive_cohort_pipeline::<_, VarCallingCliError>(0u32,
   merger, params, fixture.fetcher.clone(), &vcf_out, metadata,
   writer_cfg)`. `0u32` is correct because the synthetic fixture has
   exactly one contig ("chr1") by construction (line 188).
6. **Line 185:** `TempDir::new()` → `TempDir::new_in("tmp")` (D3
   retrofit). Add `std::fs::create_dir_all("tmp").expect("mkdir tmp")`
   above it.

### Step 6 — Second-wave compile sweep

`./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`.
The first-wave-only output we captured stopped at the unresolved-import
line for every target; this sweep flushes any cascading errors (dead
code, `clippy::needless_borrow` on new `&fasta_path` chains, etc.). Fix
in place.

### Step 7 — `cargo bench --no-run` end-to-end

`./scripts/dev.sh cargo bench --no-run`. Catches any criterion-version
drift not visible to clippy.

### Step 8 — Final gate + smoke

`./scripts/precommit-check.sh` end-to-end (fmt, clippy, full tests,
bench compile). Then one smoke run per ported bench shape to confirm
the new fetcher *runs*, not just compiles:
- `./scripts/dev.sh cargo bench --bench cohort_e2e_perf -- --test cohort_e2e_core/scaling_samples/10`
- `./scripts/dev.sh cargo bench --bench baq_perf -- --test baq_engine_read_length/150`

## Risks

- **R1 (D1 accepted).** Cohort benches re-baseline; old absolute numbers
  drop. Flag for PROJECT_STATUS update.
- **R2 (D2 accepted).** BAQ engine bench includes per-iter fetcher
  construction in the (untimed) setup. Adds ~µs of setup cost; closer
  to production reality.
- **R3.** `cohort_e2e_perf`'s `scaling_threads` sub-group measures
  *intra-chromosome* rayon work, not the per-chromosome parallel path.
  That was already the design; the H4 port doesn't change it. No
  blocker.
- **R4.** `tempfile::TempDir::new_in("tmp")` requires `tmp/` exist
  beforehand. Each fixture must `create_dir_all("tmp")` first. Cheap
  and idempotent.
- **R5.** With `-D warnings`, helpers will be scrutinised
  (`clippy::unnecessary_cast`, `clippy::format_in_format_args`, etc.).
  Fix incrementally as Step 6 surfaces them.

## Critical files

Edit:
- `examples/dhat_pileup.rs`
- `examples/profile_posterior_engine.rs`
- `benches/pileup_walker_scaling.rs`
- `benches/baq_perf.rs`
- `benches/var_calling_perf.rs`
- `benches/cohort_e2e_perf.rs`

Read-only reference:
- `src/per_sample_pileup/ref_fetcher.rs` (sealed-trait + constructor surface)
- `src/per_sample_pileup/pileup/mod.rs` (`MultiChromRefFetcher` trait)
- `src/per_sample_pileup/baq/{engine,stream}.rs` (new BAQ API)
- `src/pop_var_caller/cohort_driver.rs` (new `CohortPipelineParams` + `drive_cohort_pipeline`)
- `examples/dhat_baq.rs`, `examples/profile_cohort_e2e.rs` (already-ported precedents)
