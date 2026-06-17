# Code Review: ssr-pileup (Mark-2)
**Date:** 2026-06-17
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** the Mark-2 rebuild of the `ssr-pileup` Stage-1 module — `src/ssr/pileup/`
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the 6-file `src/ssr/pileup/` module on branch `ssr-pileup-mark2` (the Mark-2 empirical-candidate rebuild that deleted the Mark-1 rung model).
- **Reviewed against:** branch `ssr-pileup-mark2` @ `5cce21f` (diff vs `main`).
- **In-scope files:**
  - [src/ssr/pileup/mod.rs](../../../../src/ssr/pileup/mod.rs)
  - [src/ssr/pileup/alignment.rs](../../../../src/ssr/pileup/alignment.rs) (delimiter + quality gate)
  - [src/ssr/pileup/footprint.rs](../../../../src/ssr/pileup/footprint.rs) (CIGAR geometry, reach gate, region extraction)
  - [src/ssr/pileup/fetch_reads.rs](../../../../src/ssr/pileup/fetch_reads.rs) (reservoir + per-locus fetch)
  - [src/ssr/pileup/locus_tally.rs](../../../../src/ssr/pileup/locus_tally.rs) (per-locus fold)
  - [src/ssr/pileup/driver.rs](../../../../src/ssr/pileup/driver.rs) (error type, config, header build, rayon run loop, adapter)
- **Deliberately out of scope** (read only to judge the in-scope surface, no findings filed against them): [src/ssr/types.rs](../../../../src/ssr/types.rs) (Locus/Motif, reviewed 2026-06-12), [src/psp/registry_ssr.rs](../../../../src/psp/registry_ssr.rs) + `src/psp/*` (reviewed 2026-06-15), [src/bam/segment_reader.rs](../../../../src/bam/segment_reader.rs) + `src/bam/alignment_input.rs` (reviewed 2026-06-16), `src/ssr/catalog/*`.
- **Categories dispatched (11):** reliability (always), errors (always), naming (always), defaults (constants + provenance header), idiomatic (always), refactor_safety (always), module_structure (multi-file), unsafe_concurrency (rayon + the byte-identity claim), smells (always), tooling (crate gate), extras (parser/untrusted-input + hot path + stable output + diff-matches-intent).

## 2. Verdict

**Request-changes.** The architecture is sound, the determinism *mechanism* is correct by inspection, and the diff faithfully implements the stated empirical-candidate model. But three issues block merge: a broken `cargo doc` CI gate (verified), an unguarded slice that panics the whole run on a legal-but-edge-case record (empty `QUAL`), and the headline byte-identity-for-any-thread-count contract is untested on the only two paths that could break it.

## 3. Execution status

Commands run by the orchestrator in the dev container (per `CLAUDE.md`), output quoted verbatim where it drove a finding:

| Command | Exit | Result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | `Finished` — no warnings (so `--all-targets` compiles; no dangling Mark-1 bench target) |
| `cargo test --lib ssr::pileup` | 0 | `test result: ok. 40 passed; 0 failed; 0 ignored` |
| `cargo doc --no-deps` | **non-zero** | **FAIL** — `error: unresolved link to super::locus_record::aggregate` at `src/ssr/pileup/fetch_reads.rs:17:24`; `could not document pop_var_caller` |
| `cargo audit` | n/a | not run — cargo-audit not installed in the container (no dependency changes this branch) |

- Findings labeled "Needs verification": **1** (M2 — zero-length-flank reachability depends on Stage-0 catalog behavior not exercised here).

## 4. Open questions and assumptions

1. **Can the catalog emit a locus with a zero-length flank?** `Locus::new` permits `ref_bytes_start == start` (empty left flank) and `end == ref_bytes_end` (empty right flank). If Stage-0 (`ssr-catalog`) ever produces such a locus (a tract abutting a contig edge), then **M2** becomes a Blocker (wrong extracted bytes, no panic). Resolve before responding to M2.
2. **Is byte-identity required across machines/toolchains, or only across `--threads` on one build?** The module doc (driver.rs:9-11) states it unscoped. This determines whether **M5** (float-bit reproducibility) needs a golden file + fp-contract pinning or just a doc scoping sentence.
3. **Is `MAX_READS_PER_LOCUS = 1000` the production `--cap` default, or only a test value?** The runtime default is set in the CLI layer (`src/pop_var_caller/ssr_pileup.rs`, out of scope); the in-module const is referenced only by tests. Affects Mi2.
4. **Should a length-inconsistent record (empty/short `QUAL`, or a CIGAR that over-consumes the read) be a counted filter-drop or a hard error?** Determines the shape of the **B2** fix (a new `FilterCounts` bucket vs `MalformedRecord`).

## 5. Top 3 priorities

1. **B1** — fix the broken intra-doc link at [fetch_reads.rs:17](../../../../src/ssr/pileup/fetch_reads.rs#L17); the `cargo doc` gate is red. One-line change; also clears the stale module docs.
2. **B2** — guard the seq/qual slicing in [process_locus](../../../../src/ssr/pileup/driver.rs#L121-L125): an aligned record with empty `QUAL` (a legal SAM/BAM construct the SNP path already defends against) panics the entire Stage-1 run.
3. **B3** — add the two tests that actually exercise the byte-identity contract (reservoir eviction biting + a multi-chunk parallel split); today the only thread-invariance test runs 3 single-read loci in a single chunk, so the determinism-critical paths are never executed.

## 6. Findings

### Blocker

**B1: `src/ssr/pileup/fetch_reads.rs:17` — Broken intra-doc link fails the `cargo doc` CI gate**
- **Confidence:** High (verified by the orchestrator's own `cargo doc --no-deps` run, quoted in §3).
- **Categories:** tooling, smells, refactor_safety, module_structure, naming, idiomatic, defaults (convergent — surfaced by 7 sub-agents).
- **Problem:** The module doc links ``[`super::locus_record::aggregate`]``. The rebuild renamed `locus_record` → `locus_tally` and `aggregate` → `tally`, so the path does not exist. `Cargo.toml` sets `[lints.rustdoc] broken_intra_doc_links = "deny"`, so this is a hard error: `unresolved link to super::locus_record::aggregate` → `could not document pop_var_caller`.
- **Why it matters:** A required CI gate is red; the whole crate's docs cannot build. This project has historically treated a broken `cargo doc` gate as a Blocker (e.g. the 2026-06-17 ssr-pileup review B1, the 2026-06-12 ssr_types B1).
- **Suggested fix:** Retarget to the symbol that now exists, or de-link the historical name. The analyze→tally stage *did* land (`driver::process_locus` → `locus_tally::tally`):
  ```rust
  //! the `analyze_read` → [`super::locus_tally::tally`] worker stage, the
  ```
  This is part of a larger stale-doc fix (see Mi-docs); re-run `cargo doc --no-deps` to confirm green.

**B2: `src/ssr/pileup/driver.rs:121-125` — A length-inconsistent record panics the whole run; nothing validates `qual.len() == seq.len()` or the CIGAR's read-consumption**
- **Confidence:** High (reachability of the empty-`QUAL` trigger verified — see below).
- **Category:** extras.
- **Problem:** `process_locus` computes `region` purely from `read.seq.len()` via `extract_region`, then indexes **both** buffers with it:
  ```rust
  let region = extract_region(&read.cigar, fp, read.seq.len(), locus);
  let region_seq = &read.seq[region.clone()];
  let region_qual = &read.qual[region];
  ```
  `record_buf_to_mapped_read` ([alignment_input.rs:827](../../../../src/bam/alignment_input.rs#L827)) sets `qual = rb.quality_scores().as_ref().to_vec()` with **no reconciliation to `seq.len()`**, and the SSR fetch path `classify_segment_record` ([segment_reader.rs:407-448](../../../../src/bam/segment_reader.rs#L407-L448)) applies only flag/MAPQ/length filters — it never calls `cigar_is_bad` nor checks seq/qual/CIGAR consistency. Two reachable triggers:
  - **Empty `QUAL` (verified reachable, the common case).** A record with quality `*` decodes to an empty `qual` buffer while `seq` is full. Then `&read.qual[region]` slices `[..n)` of a zero-length `Vec` → index-out-of-bounds panic, killing the worker and aborting `run`. This is **proven to be a real input class the codebase already handles elsewhere**: the SNP BAQ engine explicitly guards `if read.qual.is_empty() || read.qual.len() != read.seq.len()` ([baq_engine.rs:117](../../../../src/pileup/per_sample/baq_engine.rs#L117)). The SSR path inherited none of that protection.
  - **CIGAR over-consuming the read.** `extract_region` clamps only the end: `r_start..r_end.min(read_len).max(r_start)` ([footprint.rs:135](../../../../src/ssr/pileup/footprint.rs#L135)). `ref_to_read` advances `read_cur` by the CIGAR's read-consuming ops with no bound; a CIGAR consuming more read bases than `seq` holds yields `r_start > read_len`, and `&read.seq[r_start..r_start]` panics (a `start > len` empty range still bounds-checks `start`).
  - The `debug_assert_eq!(region_seq.len(), region_qual.len(), …)` inside `delimit_read` (alignment.rs:154) does not help: it is compiled out in release, and the panicking slice happens *before* `delimit_read` is called.
- **Why it matters:** BAM/CRAM is untrusted input produced by other tools. One read with missing `QUAL` — entirely legal per the SAM spec — crashes the whole Stage-1 run with no typed error. Denial-of-service on a malformed or merely `*`-quality sample.
- **Suggested fix:** Reject (and count) length-inconsistent records once at the fetch boundary so every downstream consumer is protected, **and** clamp `r_start` as defense-in-depth. Minimal driver guard:
  ```rust
  for read in &fetched.reads {
      if read.qual.len() != read.seq.len() {
          // count as a malformed/filtered drop, not a panic
          continue;
      }
      // …
  }
  ```
  and in `extract_region`: `let r_start = r_start.min(read_len); r_start..r_end.min(read_len).max(r_start)`. Better: enforce `qual.len() == seq.len()` (and `cigar_is_bad` / read-consumption == `seq.len()`) in the fetch path, counted into `FilterCounts`. Add a malformed-input test (a spanning read with empty quality scores) asserting a counted drop, not a panic. (See Open Question 4 for the drop-vs-error decision.)

**B3: `src/ssr/pileup/fetch_reads.rs:101` + `src/ssr/pileup/driver.rs:629` — The byte-identity-for-any-thread-count contract is untested on the only paths that can break it**
- **Confidence:** High that the paths are untested; the determinism *mechanism* is verified sound by inspection (see §9), so this is test-coverage debt on a wrong-results-without-panic path, not a known bug.
- **Categories:** reliability (filed Blocker), unsafe_concurrency (filed Major) — convergent.
- **Problem:** The headline guarantee (driver.rs:9-11; the CLI promise "Output is identical for any value", driver.rs:191) rests on two runtime paths no test reaches:
  - **Reservoir eviction.** The only place the PRNG stream decides *which* reads survive is `Reservoir::offer`'s eviction branch (fetch_reads.rs:105-111). Every reservoir unit test and the end-to-end `run_output_is_thread_count_invariant` (driver.rs:629-641) run with `cap = MAX_READS_PER_LOCUS = 1000` against fixtures of **one read per locus**, so the eviction branch is never taken on the determinism path.
  - **Multi-chunk split.** `chunk_size = batch.len().div_ceil(n_chunks).max(MIN_FETCH_CHUNK)` floors at 64 (driver.rs:359), so a 3-locus batch is always a *single* chunk regardless of `--threads`. The cross-worker ordered-collect + indexed write loop (driver.rs:379-383) is never exercised at >1 chunk.
  - A regression in the chunk-boundary math, the collect/write ordering, or the reservoir offer order would produce divergent `.ssr.psp` bytes at different `--threads` and pass all 40 tests.
- **Why it matters:** Per the reliability rubric, the absence of tests for a code path that, if broken, would produce wrong results without panicking is Blocker-class — and this is the stage's load-bearing correctness contract.
- **Suggested fix:** Add (1) a `Reservoir` test offering ≫ capacity asserting the held set is byte-identical across two runs, `len() == capacity`, and a subset of the stream; (2) an end-to-end test with a locus carrying `> cap` spanning reads and a small explicit `cap` (e.g. 4 or 8), asserting `run(…, 1)` and `run(…, 4)` produce byte-identical records (eviction on the determinism path); (3) a fixture with `> 2·MIN_FETCH_CHUNK` (≥130) ascending loci asserting records at T=1 equal records at T=8 and are in ascending-`start` order (multi-chunk). Test bodies are specified in §8. (This needs `run_config` to take an explicit `cap` and ideally a `#[cfg(test)]` knob to lower `MIN_FETCH_CHUNK`.)

### Major

**M1: `src/ssr/pileup/driver.rs:64-65` — `Io(#[from] std::io::Error)` collapses five distinct I/O sites into one context-free variant**
- **Confidence:** High. **Category:** errors.
- **Problem:** `SsrPileupError::Io(#[from] io::Error)` is reached from `File::open(&cfg.catalog)` (327), `File::create(&tmp)` (338), the writer flush `into_inner` (389), `file.sync_all()` (390), and `std::fs::rename(&tmp, &cfg.output)` (391). The `#[from]` flattens all five into one variant whose `Display` is the bare `"I/O error"`, with no path and no operation. The rename is the worst loss — a cross-device `EXDEV` failure is otherwise unattributable.
- **Why it matters:** The most common production failures (missing catalog, unwritable output dir, full disk at sync, cross-device rename) become indistinguishable from logs alone.
- **Suggested fix:** Replace the blanket `#[from]` with operation-named variants carrying paths and `#[source]`:
  ```rust
  #[error("failed to open the catalog {path:?}")]
  OpenCatalog { path: PathBuf, #[source] source: std::io::Error },
  #[error("failed to create the temporary output {path:?}")]
  CreateOutput { path: PathBuf, #[source] source: std::io::Error },
  #[error("failed to flush the .ssr.psp writer")]
  FlushOutput { #[source] source: std::io::Error },
  #[error("failed to sync the output {path:?}")]
  SyncOutput { path: PathBuf, #[source] source: std::io::Error },
  #[error("failed to rename {tmp:?} into {dest:?}")]
  RenameOutput { tmp: PathBuf, dest: PathBuf, #[source] source: std::io::Error },
  ```
  then attach context at each `?` site via `.map_err(…)`.

**M2: `src/ssr/pileup/alignment.rs:283` — Zero-length-flank delimitation is untested**
- **Confidence:** Medium (Blocker-class if Open Question 1 resolves "yes"). **Category:** reliability. **Needs verification.**
- **Problem:** `delimit_read` decides off-end via `left_off = left_len > 0 && tract_start == 0` (line 283), with traceback branches keyed on `left_junction = left_len` and `right_len > 0` guards (lines 257, 268). Every test uses `ca_locus`/`locus6` with non-empty flanks on **both** sides. With an empty left flank `left_junction == 0`, the `k == 0` M/D branch sets `tract_start` while `left_off` is forced false — an uncovered path; the empty-right-flank case (`right_junction == n`, `tract_end` must stay `m`) is likewise uncovered. A bug there silently corrupts the extracted repeat bytes for a contig-edge locus.
- **Why it matters:** Wrong genotype evidence, no panic, for loci `Locus::new` permits today.
- **Suggested fix:** Add tests with `ref_bytes = CACACATTT` (empty left flank) and `GGGCACACA` (empty right flank) asserting the extracted region equals the tract for a clean read (bodies in §8). **First** confirm whether Stage-0 can emit such a locus (Open Question 1); if yes, raise to Blocker.

**M3: `src/ssr/pileup/footprint.rs:26` — `MIN_FLANK_BP` shapes the output tally but is absent from the `.ssr.psp` provenance header**
- **Confidence:** High. **Category:** defaults.
- **Problem:** `MIN_FLANK_BP = 5` is the reach-gate threshold: an unclipped read enters the reservoir only if its footprint brackets ≥ 5 ref bp past each tract end, so it determines *which* reads become evidence — i.e. it changes output bytes. `build_ssr_writer_header` (driver.rs:230-267) records `quality_q1_threshold`, `reservoir_cap`, `flank_bp` (the *catalog's* flank — a different quantity), `catalog_reference_md5`, and the filter knobs, but never `MIN_FLANK_BP`.
- **Why it matters:** A `.ssr.psp` cannot answer "what reach threshold produced this evidence?"; a future recalibration silently invalidates old outputs with no recorded discriminator.
- **Suggested fix:** Add it to the parameter map:
  ```rust
  parameters.insert("reach_min_flank_bp".to_string(),
      ParameterValue::Integer(i64::from(footprint::MIN_FLANK_BP)));
  ```

**M4: `src/ssr/pileup/alignment.rs:37,40,46,59` — The pair-HMM model constants shape delimitation output but are neither configurable nor recorded**
- **Confidence:** Medium (assumes the model is fixed for Mark-2). **Category:** defaults.
- **Problem:** `GAP_OPEN_PROB`, `GAP_EXTEND_PROB`, the `EMISSION_LN[256]` table, and `INS_EMIT_LN` drive the Viterbi delimitation (the core output-determining step), yet none is emitted into the header. Two `.ssr.psp` files from builds with a retuned gap model would be byte-different with nothing in provenance to explain it.
- **Suggested fix:** Record a `delimiter_model` tag (a string version, e.g. `"dindel-v1"`, plus the gap-open value is enough to discriminate) in `build_ssr_writer_header`; or add a `# Model provenance` doc note stating the model is frozen for the `format_version` and pointing at the spec. Prefer the recorded tag.

**M5: `src/ssr/pileup/alignment.rs:40,46-59` — Cross-platform float-bit reproducibility is asserted as determinism but not scoped or pinned**
- **Confidence:** Medium. **Category:** extras. (See Open Question 2.)
- **Problem:** DP scores are `f64` sums of `LazyLock`-computed log-probabilities (`powf`/`ln`/`exp`); `pick` resolves on strict `>`. On one target+toolchain the result is deterministic, so the *thread-count* invariance claim holds. But transcendental seeds and sum reassociation/fma can differ by 1 ULP across platforms, and a 1-ULP difference at a near-tie flips the traceback and the delimited bytes. The module doc presents determinism as an unscoped invariant.
- **Suggested fix:** Scope the guarantee in the doc ("byte-identical across runs and thread counts *on a fixed target + toolchain*; cross-platform identity is not guaranteed"). If cross-platform identity is actually required, add a golden `.ssr.psp` regenerated on the CI target plus fp-contract pinning.

**M6: `src/ssr/pileup/driver.rs:418` — `..FilterCounts::default()` in a test literal silently absorbs new filter buckets**
- **Confidence:** High. **Categories:** refactor_safety (convergent with the established `FilterCounts::merge` exhaustive-destructure convention).
- **Problem:** `qc_counts_excludes_dups_from_coverage_but_keeps_them_filtered` builds `FilterCounts { …, ..FilterCounts::default() }`. `FilterCounts` ([alignment_input.rs:109](../../../../src/bam/alignment_input.rs#L109)) is not `#[non_exhaustive]` and carries `high_mismatch_fraction`, `bad_cigar`, `baq_rejected` that the literal does not name; the `..` fills them silently. If a future bucket should feed `n_filtered`/`mapped_reads`, this regression guard keeps passing against the old arithmetic. Notably `FilterCounts::merge` already uses an exhaustive destructure for exactly this reason — the test is inconsistent with the convention.
- **Suggested fix:** Name all ten fields in the literal (drop the `..`), or use an exhaustive destructure of `FilterCounts::default()` so a new field compile-breaks the test.

**M7: `src/ssr/pileup/driver.rs:73,223-225` — `TimestampFormat` discards the parse cause and the offending value**
- **Confidence:** High. **Category:** errors.
- **Problem:** `rfc3339_now().parse().map_err(|_| SsrPileupError::TimestampFormat)?` throws away the `DatetimeParseError` and the input string; the chain dead-ends at `"could not format the run timestamp"`. Low-probability, but if a clock/format regression fires, there is nothing to debug.
- **Suggested fix:** `TimestampFormat { value: String, #[source] source: toml::value::DatetimeParseError }`, populated from the parsed string.

**M8: `src/ssr/pileup/footprint.rs:86 — `ref_to_read`'s indel branches (the reads the aligner exists for) are untested in isolation**
- **Confidence:** High. **Category:** reliability.
- **Problem:** `extract_region`'s two tests use pure-`M` / clip+`M` CIGARs where `ref_to_read` is trivial. The dual-cursor branches that matter — `target` inside a `Deletion`/`Skip` (returns `read_cur`, line 100), an `Insertion` advancing only `read_cur` (line 104), and `target == ref_end` falling through to `read_cur` (line 108) — are never exercised. A SSR read with a real indel in the window (the common case) routes through these, so a coordinate bug shifts the slice handed to the pair-HMM for *every* indel-bearing read.
- **Suggested fix:** Add `extract_region`/`ref_to_read` tests with `M D M` and `M I M` CIGARs spanning the window (bodies in §8).

**M9: `src/ssr/pileup/driver.rs:112 — `process_locus`'s multi-outcome routing is untested**
- **Confidence:** Medium. **Category:** reliability.
- **Problem:** The only coverage runs one clean spanning read per locus. No test puts a `Sequence`, a `LowQuality`, and a `BorderOffEnd` outcome at one locus, so the arm routing (driver.rs:128-140) and the `region_qual[r]` vs `region_seq[r]` slicing are unverified end-to-end; a swapped arm would shift counts between `n_low_quality`, `n_border_off_end`, and `observed` and still pass.
- **Suggested fix:** A fixture BAM mixing a clean read, a Q<15 spanning read, and a read whose right flank runs off the end at one locus; assert `depth`, `n_low_quality`, `n_border_off_end`, `observed` independently.

**M10: `src/ssr/pileup/driver.rs:197 — `build_ssr_writer_header`'s error variants and conditional parameter map are untested**
- **Confidence:** High. **Category:** reliability.
- **Problem:** Only reached via the happy-path `run` tests (always a 32-char md5, a 200 bp contig). The `ContigLengthOverflow` (line 205), `MissingMd5` (line 213), and present-vs-absent `min_mapq`/`min_read_length` insertions (247-258) are uncovered.
- **Suggested fix:** Unit tests passing a hand-built `ContigList` with (a) a `>u32::MAX` length → `ContigLengthOverflow`, (b) `None` md5 → `MissingMd5`, (c) a filter with/without `min_mapq` asserting the parameter key is present/absent.

**M11: `src/ssr/pileup/alignment.rs:294 — `passes_quality_gate` is untested at the threshold boundary and for tiny regions**
- **Confidence:** Medium. **Category:** reliability.
- **Problem:** The two gate tests use Q40 / Q10 / Q5 (far from the Phred-15 threshold) at lengths 4 and 8. The boundary (`>=` vs `>` at line 307), `len == 1` (`k = 0`), and `len == 2` are untested; a `>=`→`>` flip or a change to the `k = (len-1)/4` index would pass.
- **Suggested fix:** Add `[15]`→pass, `[14]`→fail, `[15,15]`→pass, and a region whose first-quartile element is exactly 15.

**M12: `src/ssr/pileup/alignment.rs:24-27 — DP state is a `usize`/`u8` primitive with a `_` catch-all traceback arm, not an enum**
- **Confidence:** Medium. **Category:** idiomatic.
- **Problem:** `const M/I/D: usize` are cast `as u8` at ~14 sites to build the backpointer cells and cast back on the traceback, whose `match state { M => …, D => …, _ => /* insertion */ }` (line 273) uses a wildcard for the I case. The backpointer matrix is the determinism-critical structure; an out-of-domain byte would be a silent miscompare rather than a compile error, and the `_` arm would swallow a future fourth state.
- **Suggested fix:** `#[repr(u8)] enum DpState { Match = 0, Insertion = 1, Deletion = 2 }` (1 byte, `Copy`); store `[DpState; 3]` in `back`, drop every cast, make the traceback exhaustive. `back` is scratch, never serialized, so the layout is free to change.

### Minor

- **Mi1 (convergent: extras/errors/idiomatic/refactor_safety/naming/smells): `alignment.rs:132-141` — `pick` indexes `cands[0]`/`cands[1..]` with no empty-slice guard.** A latent panic; all current callers pass non-empty 2/3-element literals so it cannot fire today. Add a `# Panics` note + `debug_assert!(!cands.is_empty())`, or use `split_first().expect(…)`.
- **Mi2 (defaults): `fetch_reads.rs:34` — `MAX_READS_PER_LOCUS` is documented as *the* cap but is only a test default in scope** (the runtime default lives in the CLI; the header records `cfg.cap`). Have the CLI default reference the const by name, or reword the doc (Open Question 3).
- **Mi3 (defaults): `alignment.rs:305` — the Q1 index `k = (len-1)/4` is undocumented.** Add a one-line comment naming the convention (nearest-rank lower quartile).
- **Mi4 (errors): `driver.rs:97-104` — unchecked `as u32` truncation + `u64` adds in `qc_counts`.** Implausible at real depth (Medium); use `saturating_add` + `u32::try_from(…).unwrap_or(u32::MAX)` so an overflow clamps visibly rather than silently corrupting QC scalars.
- **Mi5 (errors): `driver.rs:237` — `cfg.cap as i64` unchecked cast into the provenance header.** Low; `i64::try_from(cfg.cap).unwrap_or(i64::MAX)`.
- **Mi6 (extras): `fetch_reads.rs:101 — the reservoir uses `% self.seen` (modulo bias) but the doc claims "unbiased".** Determinism-preserving and intentional; soften the doc, do not change the code (would re-baseline output).
- **Mi7 (naming): `driver.rs:227,278` — `input_crams` is wrong for BAM input** and contradicts its sibling `alignment_files`. Rename the local to `input_alignment_files`; the `WriterProvenance.input_crams` field rename spans `src/psp/` (out of scope) — see §7.
- **Mi8 (smells): `alignment.rs:147-290 — `delimit_read` runs six non-interleaving phases in ~140 lines.** Extract `init_row0` / `fill_dp` / `traceback` helpers over the scratch borrows; the determinism invariants are easier to audit in isolation.
- **Mi9 (smells): `driver.rs:197-284 — `build_ssr_writer_header` (~90 lines) mixes contig validation, parameter assembly, and construction.** Extract `build_chromosome_entries` + `build_parameters`; drop the needless `input_fasta.clone()` (line 271).
- **Mi10 (smells): five near-identical `Locus` test-fixture builders** (`ca_locus`/`locus6` in alignment/footprint/fetch_reads/locus_tally + inlined at driver.rs:564). Hoist one `#[cfg(test)]` helper. Its natural home is `src/ssr/types.rs` (out of scope to edit) — see §7.
- **Mi11 (idiomatic): `alignment.rs:39-40,59 — `GAP_EXTEND_PROB` and `INS_EMIT_LN` are `LazyLock<f64>` for compile-time-knowable scalars.** Make them `const` literals (with the closed form in a comment); only `EMISSION_LN` needs `LazyLock`. (If exact `.exp()`/`.ln()` bits matter for byte-identity, keep them `static` but add a round-trip test.)
- **Mi12 (idiomatic): `alignment.rs:294-308 — `passes_quality_gate` takes `&mut ViterbiScratch` only to borrow `qual_buf`.** A `&mut Vec<u8>` would be the most-general bound.
- **Mi13 (reliability): `footprint.rs:50 — `read_footprint`'s `pos == 0` saturating-subtract branch is untested.** Add a test asserting `ref_start == 0` for `pos == 0` so a change to a plain `- 1` is caught.
- **Mi14 (reliability): `driver.rs:97 — `qc_counts` does not assert that `secondary`/`supplementary`/`unmapped` are excluded from every output.** Add a test with only those fields set asserting `(depth, n_filtered, mapped_reads) == (yielded, 0, yielded)`.
- **Mi-docs (convergent: module_structure/refactor_safety/smells/naming/tooling/defaults): `mod.rs:4-7` and `fetch_reads.rs:1-20` describe the module as "still to land"** and the fetcher-thread/bounded-queue topology that was never built (the driver is bulk-synchronous rayon `par_chunks`). Rewrite both module docs to the landed shape. Fixing fetch_reads.rs also fixes B1.

### Nits

Grouped — apply mechanically:
- `alignment.rs:24-27` — `M`/`I`/`D` cast `as u8` ~15 times; subsumed by M12 if the enum lands, otherwise store as `u8`.
- `alignment.rs:159,170,88,92` — `hap` (→ `ref_frame`), `nn` (→ `row_stride`), `prev`/`cur` (→ `prev_row`/`cur_row`): half-names/abbreviations in a long function.
- `driver.rs:99` — one-letter `d` alias for `filtered`, reused on three lines.
- `driver.rs:122-125` — add a one-line comment that `r` (the `Delimited::Region` range) is relative to the already-extracted `region`, not the full read.
- `alignment.rs:174-176` — the `ViterbiScratch { .. }` destructure drops `qual_buf`; a named-ignore `qual_buf: _` would compile-break on a rename.
- `module_structure` Nit — `LocusScratch`/`qc_counts`/`process_locus`/`to_container_record` are `pub(crate)` but driver-internal; demote to private (optional; the `ssr` tree is under `#![allow(dead_code)]`).

## 7. Out of scope observations

- **`WriterProvenance.input_crams` field name** (`src/psp/header.rs`) is misleading for BAM input crate-wide; the in-scope local rename (Mi7) is a stopgap. Suggest a separate crate-wide rename to `input_alignment_files`.
- **Shared `Locus` test-fixture helper** belongs in `src/ssr/types.rs` (Mi10) — out of scope to edit here; flag for the types owner.
- **Fetch-boundary validation** (`classify_segment_record`, `src/bam/segment_reader.rs`) is the root enabler of B2: the SSR path applies no `cigar_is_bad` / seq-qual-length consistency check that the SNP walker relies on. Decide whether the guard belongs at the fetch boundary (protecting all future consumers) or in `record_buf_to_mapped_read`. Not filed as a finding against the out-of-scope file, but the B2 fix should land there for the most leverage.

## 8. Missing tests to add now

Grouped by function. (Full bodies/specs live in the audit trail; the high-value ones:)

- `Reservoir` — **`reservoir_eviction_is_deterministic_and_subset_when_offered_far_past_capacity`**: offer `1..=10_000` into `Reservoir::new(8, seed)`; assert `seen == 10_000`, `held.len() == 8`, all held ∈ stream, and `run() == run()`. Catches a kept-set that depends on anything but (seed, order) — e.g. `held.len()` instead of `seen` in the modulus. (B3)
- `run` — **`run_output_is_thread_invariant_when_the_reservoir_cap_bites`**: a locus with `> cap` identical spanning reads, small explicit `cap` (4), assert `read_records(out@1) == read_records(out@4)`. (B3)
- `run` — **`run_output_is_thread_invariant_across_multiple_chunks`**: ≥130 ascending loci, assert records at T=1 equal records at T=8 and are `start`-ascending, and that `n_chunks > 1` at T=8. (B3)
- `delimit_read` — **`..._extracts_tract_with_empty_left_flank`** / **`..._empty_right_flank`**: `ref_bytes = CACACATTT` / `GGGCACACA`, assert the extracted region equals `CACACA`. (M2)
- `delimit_read` — **`..._returns_border_off_end_for_empty_read_region`** (`m == 0`). (M2-adjacent)
- `delimit_read` — **`..._handles_junction_indels_without_inverting_the_span`**: large insertions at both junctions; assert `Region(r)` with `r.start <= r.end` or `BorderOffEnd`, deterministic across two calls (exercises the untested `tract_start > tract_end` guard at line 285).
- `extract_region`/`ref_to_read` — **`..._maps_window_across_an_internal_deletion`** (`M6 D2 M10`) / **`...insertion`** (`M6 I2 M10`). (M8)
- `passes_quality_gate` — **`..._is_inclusive_at_the_threshold_and_handles_singletons`**: `[15]`→pass, `[14]`→fail, `[15,15]`→pass. (M11)
- `build_ssr_writer_header` — **`..._errors_on_overflow_and_missing_md5`**. (M10)
- `process_locus` — **`..._routes_low_quality_border_and_sequence_outcomes_independently`**. (M9)
- A **malformed-input test** for B2: a spanning read with empty quality scores asserts a counted drop / typed error, not a panic.
- `read_footprint` — **`..._saturates_ref_start_at_pos_zero`**. (Mi13)
- `qc_counts` — **`..._excludes_secondary_supplementary_unmapped`**. (Mi14)

## 9. What's good

- **The determinism mechanism is correct by inspection** (verified by `unsafe_concurrency`): per-locus self-contained work, reservoir seeded by `(chrom,start)`, fixed offer order, `tally` sorted by bytes, `par_chunks` indexed + ordered collect/write, `pool.install` correctly scoped, worker panics re-raised (not swallowed), and the reader `Mutex` recovers from poison so a worker panic cannot cascade or deadlock. (`driver.rs:342-383`, `fetch_reads.rs:174`, `locus_tally.rs:73-74`.) B3 is about *testing* this, not a defect in it.
- **The diff faithfully implements the stated empirical-candidate model** (verified by `extras`): the repeat region is extracted verbatim with no ladder snapping (`driver.rs:134`), the reference is only a coordinate frame, no Stage-1 likelihood, and all four determinism knobs (seed, M>D>I tie-break, 5′-junction indel rule, byte-sorted output) are present.
- **Reservoir Algorithm R and the query-window conversion are correct** (verified by `reliability`): `seen` is incremented before the modulus so `j ∈ [0, seen)`, and `seg_start`/`seg_end` = `win_start+1 ..= win_start+len` is the right 0-based-half-open → 1-based-inclusive shift.
- **`to_container_record` spells every field explicitly** (no `..`), so adding a field to the intra-crate `SsrLocusRecord` compile-breaks the adapter — the right call for a pre-alpha schema (`driver.rs:159-169`).
- **Naming is strong and domain-anchored** throughout (`Delimited`, `Footprint`, `Reservoir`, `reaches_locus`, `brackets`, `tally`); nouns for types, verbs for functions, units on the constants — `MIN_REGION_Q1` is the model the other constants should follow into the provenance header.

## 10. Commands to re-verify

Run in the dev container (`./scripts/dev.sh …`):
- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib ssr::pileup` (and the new tests from §8)
- `cargo doc --no-deps` — **currently failing**; must go green after B1.
- `cargo audit` — install + wire (deferred; no dep changes this branch).

### Author response convention
Address each finding by id (B1, M3, Mi7, …) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the §4 open questions first — Q1 determines whether M2 is a Blocker, and Q4 the shape of the B2 fix.
