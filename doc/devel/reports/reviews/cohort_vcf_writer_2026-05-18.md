# Code Review: cohort_vcf_writer
**Date:** 2026-05-18
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** new `src/var_calling/vcf_writer/` module plus its integration test
**Status:** Request-changes

---

## 1. Scope

- What was reviewed: the new in-tree module landed at commit `cd1977a`
  (cohort VCF writer) plus the matching integration test.
- Reviewed against: branch `review/vcf-writer` at HEAD `cd1977a`.
- In-scope files:
  - [src/var_calling/vcf_writer/mod.rs](../../../src/var_calling/vcf_writer/mod.rs)
  - [src/var_calling/vcf_writer/errors.rs](../../../src/var_calling/vcf_writer/errors.rs)
  - [src/var_calling/vcf_writer/sink.rs](../../../src/var_calling/vcf_writer/sink.rs)
  - [src/var_calling/vcf_writer/header.rs](../../../src/var_calling/vcf_writer/header.rs)
  - [src/var_calling/vcf_writer/record_encode.rs](../../../src/var_calling/vcf_writer/record_encode.rs)
  - [src/var_calling/vcf_writer/writer.rs](../../../src/var_calling/vcf_writer/writer.rs)
  - [tests/cohort_vcf_writer_integration.rs](../../../tests/cohort_vcf_writer_integration.rs)
- Deliberately out of scope:
  - `src/vcf_writer.rs` — legacy gVCF-merger code, scheduled deletion.
  - `src/var_calling/posterior_engine/*` — parallel SIMD WIP on another branch (also the source of the 10 fmt failures noted in §3).
  - The plan / impl-report Markdown.
- Categories dispatched (one-line reason each):
  - `reliability` — always.
  - `errors` — always; module's whole value proposition is typed loud-failure errors.
  - `naming` — always.
  - `defaults` — public API surface (`WriterConfig`, `CohortMetadata`, `CohortVcfWriter`).
  - `idiomatic` — always.
  - `refactor_safety` — always; module is expected to grow.
  - `smells` — always.
  - `extras` — stable wire format + hot-path + parser/validator at the upstream boundary.
- Categories deliberately skipped:
  - `unsafe_concurrency` — no `unsafe`, no `Arc`/`Mutex`, no channels, no async, no threads.
  - `tooling` — scope is a module, not a crate.

## 2. Verdict

**Request-changes.** One Blocker (silent allele-index drop in the
`AC` tally), an `errors.rs` redesign-class cluster (4 Majors), a
malformed-input panic cluster in `record_encode.rs`, and a
durability gap in `sink.rs::finish` that undermines the
documented atomic-rename contract. The code is clean otherwise:
clippy clean, fmt clean for in-scope files, all 23 unit tests +
4 integration tests pass.

## 3. Execution status

| Command | Exit | Result |
|---------|------|--------|
| `cargo fmt --all -- --check` | non-zero | 10 hunks need formatting — **all in OUT-OF-SCOPE files** (`benches/var_calling_perf.rs`, `examples/profile_posterior_engine.rs`, `src/var_calling/posterior_engine.rs`, `src/var_calling/posterior_engine/{backends,interp,shape}.rs`). vcf_writer files pass. |
| `cargo clippy --lib --tests --all-features -- -D warnings` | 0 | clean. |
| `cargo test --lib var_calling::vcf_writer` | 0 | `23 passed; 0 failed`. |
| `cargo test --test cohort_vcf_writer_integration` | 0 | `4 passed; 0 failed`. |

Not run:

- `cargo doc --no-deps` — the parallel `posterior_engine` WIP has
  broken intra-doc links unrelated to vcf_writer; re-run once
  that branch lands.
- `cargo audit` — not part of the standard pre-merge gate in this
  repo.

**"Needs verification" findings:** none. All Blocker/Major issues
are High confidence; a single Major (`M8`) on the `num_obs as i32`
wrap is Medium confidence because the trigger requires
`num_obs > i32::MAX` per allele — implausible at typical SNP
coverage but mechanical to fix.

## 4. Open questions and assumptions

1. **Should `Default for WriterConfig` exist at all?** Defaults,
   idiomatic, errors, and extras all flagged this independently.
   The current `Default` synthesises `output: PathBuf::new()` —
   an empty path that fails at `File::create` time. The choices
   are (a) drop the impl and force `output` to be an explicit
   constructor argument; (b) keep the impl but document that
   `output` is a required override. The fix proposed under M11
   assumes (a). Affects M11, M12, Mi11.
2. **Should the writer's `Encode(String)` collapse become typed
   variants with `#[source]` wrapping noodles error types?** This
   is the cleanest fix for M2 but it exposes `noodles_vcf` /
   `noodles_core` error types in the public `VcfWriteError`. If
   keeping noodles types out of the public surface is a hard
   constraint, fall back to `Encode(Box<dyn Error + Send +
   Sync>)` per the idiomatic agent's suggestion. Affects M2, M4.
3. **Is a parent-directory fsync required for this writer's
   durability contract?** M10 raises this. The module doc-comment
   advertises the tmp+rename dance as the loud-failure mode; a
   parent-dir fsync closes the only remaining crash window. The
   cost is one extra `open` + `sync_all` on `finish`. Confirm
   the contract.
4. **`VcfWriteError` `#[non_exhaustive]`?** This module's
   `VcfWriteError` is `pub use`-exported. The implementation
   report enumerates planned new variants (depth overflow,
   structural mismatch, allele-index out of bounds, etc.). Adding
   any of those is a semver break today. M3 proposes adding
   `#[non_exhaustive]` now. Affects M3.

## 5. Top 3 priorities

1. **`B1` — silent allele-index drop in `tally_called_alleles`.**
   Produces VCFs where `sum(AC) + REF != AN`, with no error
   signal. Highest-impact, smallest fix.
2. **`errors.rs` redesign (M1 + M2 + M3 + M4 + M5).** Five
   Majors compound onto one error type: opaque `Io` collapse,
   stringly `Encode`, missing `#[non_exhaustive]`, mechanism-
   prefixed `Display` messages, and a field-inverted
   `SampleCountMismatch` message. One coordinated pass through
   `errors.rs` and the seven `?`/`map_err` sites that feed it
   resolves all five with one set of tests.
3. **Panic-safety + DP-overflow cluster in `record_encode.rs`
   (M6 + M7 + M8).** All three convert "malformed upstream"
   into either a panic at the Stage 6 sink or a silently wrong
   VCF column. Same shape of fix — length/range validation
   at the top of `encode` (or per-site `?` with a typed error)
   plus regression tests.

## 6. Findings

### Blocker

#### B1: `src/var_calling/vcf_writer/record_encode.rs:256` — `tally_called_alleles` silently drops out-of-bounds ALT allele indices
**Categories:** reliability (primary) — convergent cross-category notes from errors, naming, defaults, idiomatic, smells, extras.
**Confidence:** High.

**Problem.** When walking each sample's argmax genotype, the loop
decodes allele indices from the canonical `genotype_order`
table and tallies `AC`. If `allele_idx > n_alts`, the guard
silently drops the index:

```rust
let alt_pos = (allele_idx as usize) - 1;
// Defensive bounds check; allele_idx > n_alts would mean
// the upstream produced a genotype indexing an allele
// beyond the record's allele set.
if alt_pos < ac_per_alt.len() {
    ac_per_alt[alt_pos] += 1;
}
```

`an_total` was already incremented one line above (line 248), so
the emitted VCF has `sum(AC) + (called REF) != AN` — a VCF-spec
invariant violation — without any error returned. `VcfWriteError`
has no variant for the situation. The project's "No logs — promote
to typed errors" memory rules out the silent-drop shape.

**Why it matters.** Downstream stats (allele-frequency QC, sample
counts) read wrong numbers from a VCF that parses cleanly. This
is "swallowed errors that hide failures" (Blocker per the rubric).

**Suggested fix.** Promote to a typed error and refuse to emit
the record. Reuse the existing five-field shape:

```rust
// errors.rs — new variant (or extend GenotypeIndexOutOfBounds)
#[error(
    "record at {chrom_id}:{pos}: sample {sample_idx} genotype decodes \
     allele index {allele_idx} but record carries only {n_alleles} alleles"
)]
AlleleIndexOutOfBounds {
    chrom_id: u32, pos: u32, sample_idx: usize,
    allele_idx: u8, n_alleles: usize,
},

// record_encode.rs — replace the `if alt_pos < …` guard:
if alt_pos >= ac_per_alt.len() {
    return Err(VcfWriteError::AlleleIndexOutOfBounds {
        chrom_id: record.locus.chrom_id,
        pos: record.locus.start,
        sample_idx,
        allele_idx,
        n_alleles: ac_per_alt.len() + 1,
    });
}
ac_per_alt[alt_pos] += 1;
```

Add the regression test
`tally_called_alleles_errors_on_out_of_bounds_allele_index`
described in §8.

### Major

#### M1: `src/var_calling/vcf_writer/errors.rs:20` — `Io(#[from] io::Error)` collapses 8 distinct call sites
**Confidence:** High.

**Problem.** `Io(#[from] io::Error)` auto-derives `From<io::Error>`
for `VcfWriteError`, which fires at *eight* in-scope sites:
`File::create` (`sink.rs:40`), `BufWriter::into_inner` (`sink.rs:58`),
`bgzf.finish` (`sink.rs:59`), `file.sync_all` (`sink.rs:61`),
`fs::rename` (`sink.rs:64`), `noodles_vcf::io::Writer::write_header`
(`writer.rs:71`), `noodles_vcf::io::Writer::write_variant_record`
(`writer.rs:109`), and any record-encode I/O reached via `?`. A
rename failure is reported with the same shape, message, and error
chain as a tmp-create failure or a bgzf-flush failure; operators
cannot tell which guarantee was violated.

**Why it matters.** The writer's whole point is that a failure
during `finish` is distinguishable from a failure during
construction or per-record write. The error type currently
hides exactly that distinction.

**Suggested fix.** Adopt operation-named variants carrying the
relevant context (`tmp_path`, `final_path`, `record_number`) plus a
`#[source] io::Error`. Drop `#[from]`; attach the variant with
`.map_err` at each `?`. Sketch:

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum VcfWriteError {
    #[error("failed to create tmp output {tmp_path}")]
    CreateTmp { tmp_path: PathBuf, #[source] source: io::Error },
    #[error("failed to write VCF header to {tmp_path}")]
    WriteHeader { tmp_path: PathBuf, #[source] source: io::Error },
    #[error("failed to write record at {chrom_id}:{pos}")]
    WriteRecord { chrom_id: u32, pos: u32, #[source] source: io::Error },
    #[error("failed to flush bgzf sink for {tmp_path}")]
    FinishBgzf { tmp_path: PathBuf, #[source] source: io::Error },
    #[error("failed to fsync {tmp_path}")]
    Fsync { tmp_path: PathBuf, #[source] source: io::Error },
    #[error("failed to rename {tmp_path} -> {final_path}")]
    Rename { tmp_path: PathBuf, final_path: PathBuf, #[source] source: io::Error },
    // ... existing typed variants kept
}
```

#### M2: `src/var_calling/vcf_writer/errors.rs:22` — `Encode(String)` is a catch-all stringly-typed variant that flattens the source
**Categories:** errors (primary), idiomatic (convergent).
**Confidence:** High.

**Problem.** `Encode(String)` collapses six distinct failure modes
into one variant: `##source` key parse (`header.rs:61`), `##source`
insert (`header.rs:67`), `##commandline` key parse (`header.rs:72`),
`##commandline` insert (`header.rs:77`), allele UTF-8
(`record_encode.rs:140-148`), `usize`-to-`Position` failure
(`record_encode.rs:84-95`). Each site does
`.map_err(|e| VcfWriteError::Encode(format!("…: {e}")))`,
embedding the noodles cause into the parent's `Display` string —
breaking error-chain rendering (`{:#}`, `eyre`, `tracing`) and
making programmatic recovery (`matches!` on the underlying parse
error) impossible.

**Suggested fix.** Replace with operation-named variants that hold
typed causes. Two viable shapes:

```rust
// Shape A — keep noodles types out of the public surface:
#[error("vcf encode")]
Encode(#[source] Box<dyn std::error::Error + Send + Sync>),

// Shape B — expose typed sources (depends on Open question 2):
#[error("record at {chrom_id}:{pos}: invalid VCF position")]
InvalidPosition { chrom_id: u32, pos: u32,
                  #[source] source: noodles_core::position::TryFromIntError },
#[error("record at {chrom_id}:{pos}: allele bytes are not valid UTF-8")]
InvalidAlleleBytes { chrom_id: u32, pos: u32,
                     #[source] source: std::str::Utf8Error },
```

Either way the call sites collapse to
`.map_err(|e| VcfWriteError::Encode(Box::new(e)))?` (Shape A) or
the per-variant builder (Shape B) — no `format!` per failure.

#### M3: `src/var_calling/vcf_writer/errors.rs:17` — public `VcfWriteError` is not `#[non_exhaustive]`
**Confidence:** High.

**Problem.** `VcfWriteError` is `pub use`-exported from
[mod.rs:35](../../../src/var_calling/vcf_writer/mod.rs#L35) and so
is part of the crate's public API. It is a plain `pub enum` with no
`#[non_exhaustive]`. The implementation report's "Open" list
already enumerates several planned variants (`PL`-related,
contamination-fraction passthrough, depth-overflow). Each new one
is a semver break against today's enum.

**Suggested fix.** Mark the enum `#[non_exhaustive]` now, while the
crate is internal-only and the cost is zero.

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum VcfWriteError { ... }
```

#### M4: `src/var_calling/vcf_writer/errors.rs:19-27` — `Display` messages prepend mechanism, flatten source
**Confidence:** High.

**Problem.** Three variants use `"<mechanism>: {0}"`: `"io: {0}"`,
`"vcf encode: {0}"`, `"invalid metadata: {0}"`. This concatenates
the inner cause into the parent's message and uses mechanism-keyed
prefixes rather than describing *what was attempted*. The rule
sheet is explicit: messages should describe operation + input;
let `source()` carry the cause. Tools that walk the error chain
print `"io: foo: foo"` (the cause appears twice).

**Suggested fix.** Drop the `: {0}` from the formats; let the
inner cause be reached via `#[source]`. After M1 + M2 land this
is mechanical:

```rust
#[error("failed to write VCF record to {tmp_path}")]
WriteRecord { tmp_path: PathBuf, #[source] source: io::Error },
```

not `#[error("io: {source}")]`.

#### M5: `src/var_calling/vcf_writer/errors.rs:62` — `SampleCountMismatch` message inverts the meaning of its two count fields
**Confidence:** High.

**Problem.** The format string reads `"posterior arrays were sized
for {expected_samples} samples but the cohort metadata names
{got_samples}"`. The call site at
[record_encode.rs:66-73](../../../src/var_calling/vcf_writer/record_encode.rs#L66-L73)
populates `expected_samples = <cohort metadata size>` and
`got_samples = record.n_samples` — the opposite of what the
message says. The existing test
(`sample_count_mismatch_is_an_error`) only asserts the variant
matches, not the rendered text, so the bug is silent under
`cargo test`.

**Why it matters.** A wrong-direction diagnostic actively misleads
operators investigating an upstream mismatch — the natural fix
("resize the posterior arrays") would be exactly backwards.

**Suggested fix.** Swap the field roles in the message and lock
the wording with a test:

```rust
#[error(
    "record at {chrom_id}:{pos}: cohort metadata names {expected_samples} \
     samples but the posterior arrays carry {got_samples}"
)]
SampleCountMismatch { ... }
```

Add `sample_count_mismatch_message_names_cohort_first` per §8.

#### M6: `src/var_calling/vcf_writer/record_encode.rs:97,174,237,274,287,293,304` — multiple unchecked slice indexes panic on malformed `PosteriorRecord`
**Categories:** errors (primary), reliability (convergent), refactor_safety (convergent), extras (convergent).
**Confidence:** High.

**Problem.** `encode` validates `record.n_samples ==
expected_samples` but does not check that any of the dependent
vectors are sized accordingly. Each of these sites panics on a
malformed upstream:

- `record.alleles[0].seq` (line 97) — panics on empty `alleles`.
- `record.alleles[1..]` (line 105) — same shape.
- `record.allele_frequencies[1..]` (line 174) — panics if
  `allele_frequencies.len() < 1` once `n_alts > 0`.
- `record.best_genotype.iter()` + `record.best_genotype[sample_idx]`
  (lines 237, 274) — `best_genotype` length not checked.
- `record.gq_phred[sample_idx]` (line 287).
- `record.scalars_row(sample_idx)` (line 289) and `scalars[a]`
  for `a in 0..n_alleles` (line 294) — `scalars.len() == n_samples
  * n_alleles` never verified.
- `record.posteriors_row(sample_idx)` (line 305) — same shape.

Defensive guards elsewhere in `encode` (`SampleCountMismatch`,
`UnknownChromId`, `GenotypeIndexOutOfBounds`) establish the
"typed error, never panic" pattern; these vectors are the
missing pieces. A panic at the Stage 6 sink crashes the entire
pipeline with a stack trace and leaves `<output>.tmp` on disk
with no diagnosis, bypassing the loud-failure contract.

**Suggested fix.** Validate every dependent length once at the
top of `encode` and return a typed `VcfWriteError::InconsistentRecord
{ chrom_id, pos, field, expected, actual }` on mismatch:

```rust
let n_alleles = record.alleles.len();
if n_alleles == 0 {
    return Err(VcfWriteError::InconsistentRecord {
        chrom_id: record.locus.chrom_id, pos: record.locus.start,
        field: "alleles", expected: 1, actual: 0,
    });
}
for (field, expected, actual) in [
    ("best_genotype", record.n_samples, record.best_genotype.len()),
    ("gq_phred",      record.n_samples, record.gq_phred.len()),
    ("allele_frequencies", n_alleles,   record.allele_frequencies.len()),
    ("scalars",       record.n_samples * n_alleles, record.scalars.len()),
    ("posteriors",    record.n_samples * record.n_genotypes,
                                      record.posteriors.len()),
    ("chain_anchor_flags", record.n_samples * n_alleles,
                                      record.chain_anchor_flags.len()),
] {
    if expected != actual {
        return Err(VcfWriteError::InconsistentRecord {
            chrom_id: record.locus.chrom_id, pos: record.locus.start,
            field, expected, actual,
        });
    }
}
```

Add the regression tests `encode_errors_on_empty_alleles` and one
per dependent vector (or one parameterised test covering all six).

#### M7: `src/var_calling/vcf_writer/record_encode.rs:213,291` — DP silently saturates at `i32::MAX` on overflow
**Categories:** reliability (primary), errors (convergent), naming (convergent), defaults (convergent), idiomatic (convergent), smells (convergent), extras (convergent). **Convergent finding** — six categories independently flagged this.
**Confidence:** High.

**Problem.** Two sites cast `u64` depths down to `i32` with
`unwrap_or(i32::MAX)`:

```rust
i32::try_from(dp_total).unwrap_or(i32::MAX)         // INFO=DP
i32::try_from(dp_sample).unwrap_or(i32::MAX)        // FORMAT=DP
```

On overflow (sum of `num_obs` > ~2.1 billion across the cohort, or
per-sample), the VCF reports `2147483647` with no signal — a
silent data-corruption pattern the project's "No logs — promote to
typed errors" memory explicitly bans. Filed at Major (rather than
Blocker) only because the trigger is implausible at typical SNP
coverage; the silent shape itself is Blocker-class.

**Suggested fix.** Promote to a typed error.

```rust
// errors.rs:
#[error("record at {chrom_id}:{pos}: \
         depth {depth} overflows i32 (sample {sample_idx:?})")]
DepthOverflow {
    chrom_id: u32, pos: u32, sample_idx: Option<usize>, depth: u64,
},

// record_encode.rs — replace each `unwrap_or(i32::MAX)` with:
let dp_total_i32 = i32::try_from(dp_total).map_err(|_| {
    VcfWriteError::DepthOverflow {
        chrom_id: record.locus.chrom_id, pos: record.locus.start,
        sample_idx: None, depth: dp_total,
    }
})?;
```

Add `encode_errors_when_per_sample_depth_overflows_i32`.

#### M8: `src/var_calling/vcf_writer/record_encode.rs:294` — `num_obs as i32` is a silent wrapping cast (negative AD on overflow)
**Confidence:** Medium (Major because the symptom is Blocker-class but the trigger requires `num_obs > i32::MAX`).

**Problem.** Adjacent to the saturating cast in M7, this line uses
a plain `as` cast:

```rust
let ad_values: Vec<Option<i32>> = (0..n_alleles)
    .map(|a| Some(scalars[a].num_obs as i32))
    .collect();
```

`AlleleSupportStats::num_obs` is `u32`
([per_sample_pileup/pileup/mod.rs:348](../../../src/per_sample_pileup/pileup/mod.rs#L348)).
Any `num_obs > i32::MAX` wraps to a negative value silently — data
corruption in the `AD` cell. Contrast with the DP-total path two
lines down that correctly uses `try_from(...).unwrap_or(i32::MAX)`
(which is itself wrong per M7, but at least non-negative).

**Suggested fix.** Use the same shape as the DP fix in M7:

```rust
let ad_values: Vec<Option<i32>> = (0..n_alleles)
    .map(|a| {
        i32::try_from(scalars[a].num_obs).map(Some).map_err(|_| {
            VcfWriteError::DepthOverflow {
                chrom_id: record.locus.chrom_id, pos: record.locus.start,
                sample_idx: Some(sample_idx),
                depth: u64::from(scalars[a].num_obs),
            }
        })
    })
    .collect::<Result<_, _>>()?;
```

#### M9: `src/var_calling/vcf_writer/record_encode.rs:287` — GQ cap `99.0` is an unnamed magic number that silently re-clamps below the engine's `max_gq_phred`
**Categories:** reliability (primary), naming (convergent), defaults (convergent). **Convergent finding.**
**Confidence:** High (Medium for the drift symptom; High for the missing constant).

**Problem.** Per-sample GQ:

```rust
let gq = record.gq_phred[sample_idx].round().clamp(0.0, 99.0) as i32;
```

Three issues compound:

1. `99.0` is an unnamed literal — there is no `pub const` to
   reference, while the sister `QUAL_MAX = 9999.0` immediately
   above is fully named and documented.
2. `PosteriorRecord.gq_phred` is documented
   ([posterior_engine.rs:428-430](../../../src/var_calling/posterior_engine.rs#L428-L430))
   as already clamped to `PosteriorEngineConfig::max_gq_phred`.
   Two independent ceilings on the same value will drift the
   moment one is bumped.
3. No test asserts the writer cap is `>= engine cap`.

**Suggested fix.** One named constant + a drift-detecting test.

```rust
/// Hard cap on per-sample `GQ`. GATK and bcftools convention is to
/// cap GQ at 99 (single-byte Phred). MUST be >= the engine's
/// `max_gq_phred`; the test below latches drift.
pub(super) const GQ_MAX: f32 = 99.0;

let gq = record.gq_phred[sample_idx]
    .round()
    .clamp(0.0, GQ_MAX as f64) as i32;
```

```rust
#[test]
fn gq_writer_cap_is_at_least_engine_cap() {
    use crate::var_calling::posterior_engine::PosteriorEngineConfig;
    let engine_cap = PosteriorEngineConfig::default().max_gq_phred;
    assert!((GQ_MAX as f64) >= engine_cap,
            "writer GQ cap {GQ_MAX} < engine cap {engine_cap}");
}
```

#### M10: `src/var_calling/vcf_writer/sink.rs:54` — `finish` does not fsync the parent directory after the atomic rename
**Confidence:** High.

**Problem.** `finish` flushes the file, calls `file.sync_all()`,
drops the file, then `fs::rename(&tmp_path, final_path)`. On
ext4/xfs/btrfs the rename metadata-update lives in the parent
directory's journal entry. A power loss between the rename and
the parent-dir journal flush can leave the file durable on disk
but the directory entry pointing at `<final_path>` lost,
surfacing as either the old name (`.tmp`) or nothing. The module
docs promise "atomic-renames `<output>.tmp` → `<output>`" as the
durability contract — this is the gap.

**Suggested fix.** After the rename, fsync the parent directory:

```rust
fs::rename(&tmp_path, final_path)?;
if let Some(parent) = final_path.parent() {
    let dir_path = if parent.as_os_str().is_empty() {
        std::path::Path::new(".")
    } else {
        parent
    };
    let dir = std::fs::File::open(dir_path)?;
    dir.sync_all()?;
}
Ok(())
```

Smoke-test it via `finish_renames_and_fsyncs_parent_dir` (the
success-path check that the new code is exercised).

#### M11: `src/var_calling/vcf_writer/mod.rs:62` — `Default for WriterConfig` produces a runtime-behavioural + structurally-broken default
**Categories:** defaults (primary), idiomatic (convergent), errors (convergent), extras (convergent). **Convergent finding.**
**Confidence:** High.

**Problem.** `WriterConfig::default()` returns:

```rust
Self {
    output: PathBuf::new(),
    default_filter_pass: true,
    emit_gp: false,
}
```

`output: PathBuf::new()` is an empty path that, used as-is, causes
`SinkKind::open_tmp` to write a `.tmp` file in the process cwd
and `fs::rename` to error confusingly. The two boolean fields
encode behaviour (always-PASS, no GP) via prose docs only — there
are no named `pub const`s to reference, so the docs and the
literals will drift. The in-module tests use `WriterConfig::default()`
at 9+ call sites; each silently inherits the broken `output`.

**Suggested fix.** Drop `impl Default` and add explicit constants
+ a constructor. Update affected tests to call `WriterConfig::new(out_path)`.

```rust
/// `WriterConfig::default_filter_pass` default. v1 always emits PASS.
pub const DEFAULT_FILTER_PASS: bool = true;

/// `WriterConfig::emit_gp` default. Off — `GP` is `Number=G` and
/// most consumers do not read it.
pub const DEFAULT_EMIT_GP: bool = false;

impl WriterConfig {
    pub fn new(output: PathBuf) -> Self {
        Self {
            output,
            default_filter_pass: DEFAULT_FILTER_PASS,
            emit_gp: DEFAULT_EMIT_GP,
        }
    }
}
```

See Open question 1.

#### M12: `src/var_calling/vcf_writer/writer.rs:33` — running `CohortVcfWriter` config is not inspectable
**Confidence:** High.

**Problem.** `CohortVcfWriter` holds `config: WriterConfig`
privately, exposes no accessor, has no `Debug` impl, and prints no
startup summary. The defaults rule "Inspectable at runtime"
requires every effective default be recoverable from a running
instance. Today, "are you emitting GP? what output path?
always-PASS?" can only be answered by re-reading the source.

**Suggested fix.** Add an accessor:

```rust
impl CohortVcfWriter {
    /// The config this writer was constructed with (frozen at `new`).
    pub fn config(&self) -> &WriterConfig {
        &self.config
    }
}
```

A manual `Debug` impl is a nice-to-have on top.

#### M13: `src/var_calling/vcf_writer/writer.rs:62` — `<output>.tmp` is leaked on header-write failure during `new`
**Categories:** reliability (primary), errors (convergent), idiomatic (convergent).
**Confidence:** High.

**Problem.** `new` calls `SinkKind::open_tmp(&final_path)?`
(creating `<final>.tmp` on disk), then constructs
`noodles_vcf::io::Writer` and calls `inner.write_header(&header)?`.
If `write_header` fails, the tmp file stays on disk with partial
bytes, `Self` was never constructed, and no `Drop` impl exists to
clean it up. There is no test for this path. The behaviour is
indistinguishable from "user forgot to call `finish`", which
muddles the forensic value of the leftover `.tmp`.

**Suggested fix.** Wrap the header-write so a failure removes
the tmp file before bubbling:

```rust
let sink = SinkKind::open_tmp(&final_path)?;
let mut inner = noodles_vcf::io::Writer::new(sink);
if let Err(e) = inner.write_header(&header) {
    let _ = std::fs::remove_file(super::sink::tmp_path_for(&final_path));
    return Err(e.into());
}
```

Regression test: `new_cleans_up_tmp_on_header_write_failure`
(may need a test-only seam that injects a sink that errors on
the first write).

#### M14: `src/var_calling/vcf_writer/writer.rs:33` — `CohortVcfWriter` is missing `#[must_use]`
**Confidence:** High.

**Problem.** `CohortVcfWriter` is a transaction handle: it opens
`<output>.tmp`, accepts records, and only `finish()` flushes the
bgzf EOF, fsyncs, and atomically renames. Module docs treat a
forgotten `finish` as the loud-failure mode, but the type has
neither `#[must_use]` nor a `Drop`-based warning. A caller writing
`let _ = CohortVcfWriter::new(...)?;` (or letting it drop without
chained `finish`) gets no signal.

**Suggested fix.** Annotate the type:

```rust
#[must_use = "CohortVcfWriter must be finalised by calling \
              `.finish()`; otherwise the output is left at \
              `<output>.tmp` and no `<output>` is produced"]
pub struct CohortVcfWriter { ... }
```

#### M15: `benches/` — no criterion benchmark for the cohort VCF writer despite the hot-path designation
**Confidence:** High.

**Problem.** The implementation report explicitly flags the
writer as on a potential hot path for cohort calls
(1000+ samples × millions of records). The extras rule requires
`criterion` benches with a `// REGRESSION THRESHOLD: N%` comment
for any hot-path module. The `benches/` directory holds six
benches today, none of which target `vcf_writer`. There is no
regression guard against a future change that drops the per-writer
format-key reuse and rebuilds keys per record (or drops the
encode short-circuits), and the existing per-record allocation
hotspots (genotype-table rebuild, `Vec<Vec<Option<SampleValue>>>`,
`gt_buf.clone()` per sample) have no baseline.

**Suggested fix.** Add `benches/vcf_writer_perf.rs` measuring:

1. `CohortVcfWriter::write_record` throughput at 1000 samples ×
   biallelic SNPs, in-memory sink.
2. The same with `emit_gp = true`.
3. `record_encode::encode` alone (isolates noodles encode cost).

Add `// REGRESSION THRESHOLD: 5%` on each bench. Wire into the
existing `[[bench]]` entries pattern in `Cargo.toml`.

### Minor

#### Mi1: `src/var_calling/vcf_writer/writer.rs:89` — missing regression test that `last_locus` is not advanced when `write_record` errors
**Confidence:** High.

The doc comment explicitly promises that a `RecordOutOfOrder`
keeps the writer running with the previous accepted record as the
new comparison anchor. The implementation honours this (the `?`
on `encode` short-circuits before `self.last_locus = Some(locus)`),
but no test pins it. A future refactor moving the update earlier
would regress silently. Add `out_of_order_does_not_advance_last_locus`
per §8.

#### Mi2: `src/var_calling/vcf_writer/sink.rs:95` — `path_is_bgzf` is case-sensitive but the convention is undocumented
**Confidence:** High.

A user writing to `out.VCF.GZ` gets the plain-text sink, producing
a file with a `.GZ` suffix that isn't gzipped — a quiet usability
trap. Either compare lowercased, or doc-comment the case-sensitivity
on `WriterConfig::output`. Pick consciously.

#### Mi3: `src/var_calling/vcf_writer/header.rs:198` — empty `contigs` is silently accepted, neither documented nor tested
**Confidence:** Medium.

`validate_metadata` rejects empty `sample_names` but accepts empty
`contigs`. A VCF with zero `##contig=` lines is structurally valid
but operationally meaningless. Either reject or document — and
add `empty_contigs_accepted_produces_header_with_no_contig_lines`
(or `empty_contigs_rejected`).

#### Mi4: `src/var_calling/vcf_writer/header.rs:82-97` — duplicated overflow checks; `ContigLengthOverflow` has no test
**Categories:** reliability + errors + naming + refactor_safety + smells (convergent — five categories surfaced it).
**Confidence:** High.

`usize::try_from(u32)` is infallible on both 32-bit and 64-bit
targets (the project deploys to x86_64 and aarch64), so the first
arm is dead and the `entry.length > i32::MAX as u32` arm is the
only one that ever fires. `ContigLengthOverflow` has zero test
coverage. Fix:

```rust
if entry.length > i32::MAX as u32 {
    return Err(VcfWriteError::ContigLengthOverflow {
        name: entry.name.clone(), length: entry.length,
    });
}
let length_usize = entry.length as usize;
```

Add `contig_length_above_i32_max_rejected`.

#### Mi5: `src/var_calling/vcf_writer/header.rs:60,71` — no test exercises empty `tool_string` / empty `command_line` skip paths
**Confidence:** High.

`build_vcf_header` skips the `##source` line when `tool_string` is
empty (line 60) and `##commandline` when `command_line` is empty
(line 71). The fixture used in every test sets both to non-empty,
so the skip paths are dead under `cargo test`. The doc on
`CohortMetadata::command_line` explicitly contracts the skip. Add
`header_omits_source_when_tool_string_empty` and
`header_omits_commandline_when_command_line_empty`.

#### Mi6: `src/var_calling/vcf_writer/header.rs:60-80` — duplicate `##source` / `##commandline` insert blocks
**Confidence:** High.

Two 9-line copies differing only in key name and field. The
extension point is obvious (next `##KEY=value` line will copy
the same 10 lines). Extract a helper:

```rust
fn insert_unstructured(
    builder: HeaderBuilder, key: &str, value: &str,
) -> Result<HeaderBuilder, VcfWriteError> {
    if value.is_empty() { return Ok(builder); }
    let parsed = key.parse::<…::Other>()
        .map_err(|e| /* per M2 shape */)?;
    builder.insert(parsed, HeaderValue::from(value))
        .map_err(|e| /* per M2 shape */)
}
```

Then both sites collapse to one line each.

#### Mi7: `src/var_calling/vcf_writer/record_encode.rs:237,275` — duplicated genotype-table lookup + identical 5-field error
**Confidence:** High.

`tally_called_alleles` and `build_samples` both do
`table.get(gt_idx).ok_or(VcfWriteError::GenotypeIndexOutOfBounds {
chrom_id, pos, sample_idx, got, n_genotypes })?`. Extract a
private `lookup_genotype(table, record, sample_idx, gt_idx)`
helper; both call sites become one line. Adding a sixth field
to the error variant later becomes a one-site edit.

#### Mi8: `src/var_calling/vcf_writer/record_encode.rs:272` — `gt_buf` clear/reuse + `.clone()` per row defeats the amortisation
**Confidence:** High.

The reuse pattern is meant to retain the buffer's capacity across
iterations, but `gt_buf.clone()` per row allocates a fresh
`String` anyway — same allocation count as the obvious version,
with worse readability. Either drop the buffer (one fresh
`String` per iteration, identical count) or `std::mem::take(&mut
gt_buf)` to actually amortise.

#### Mi9: `src/var_calling/vcf_writer/writer.rs:65` — `metadata.contigs.clone()` immediately before `metadata` is dropped
**Confidence:** High.

`metadata` is taken by value and not used after the contigs are
cloned out. Move instead:

```rust
let header = build_vcf_header(&metadata, &config)?;
let expected_samples = metadata.sample_names.len();
let contigs = metadata.contigs; // move
let final_path = config.output.clone();
```

#### Mi10: `src/var_calling/vcf_writer/record_encode.rs:630` — `..Default::default()` in `format_keys_track_emit_gp` test breaks the compile-fence-on-new-field discipline
**Confidence:** High.

Every other `WriterConfig` literal in scope (12 sites) spells the
three fields explicitly; this one elides them. When the planned
new flags land, the compiler will flag every other literal but
this one — and the test will keep passing under whatever default
the new field takes. Fix is mechanical:

```rust
let on_cfg = WriterConfig {
    output: PathBuf::new(),
    default_filter_pass: true,
    emit_gp: true,
};
```

(Matches `header_includes_gp_when_enabled` at
[header.rs:290-294](../../../src/var_calling/vcf_writer/header.rs#L290-L294).)

#### Mi11: `src/var_calling/vcf_writer/mod.rs:43` — `WriterConfig` is under-qualified vs `CohortVcfWriter` / `CohortMetadata`
**Confidence:** Medium (depends on Open question 1's outcome).

`pub struct WriterConfig` reads asymmetrically alongside
`CohortVcfWriter` and `CohortMetadata` — both prefixed `Cohort*`.
At a use site like
[tests/cohort_vcf_writer_integration.rs:21-23](../../../tests/cohort_vcf_writer_integration.rs#L21-L23),
the three `Cohort*`-prefixed types and the bare `WriterConfig`
are jarring. Rename to `CohortVcfWriterConfig`. Mechanical
rename; no behaviour change.

#### Mi12: `src/var_calling/vcf_writer/header.rs:45` — `tool_string` field name encodes type, not domain meaning
**Confidence:** High.

The field doc says "goes into `##source=...`". A name like
`source_label` or `tool_name` matches the domain. The sibling
`command_line` is correctly named after what it represents.
Public-API field — renaming later is a breaking change, easier
now.

#### Mi13: `tests/cohort_vcf_writer_integration.rs:25,61,89,127` — `ref_allele` factory used to build ALT alleles
**Confidence:** High.

The integration test defines `fn ref_allele(seq) -> MergedAllele`
and then calls it for both REF and ALT positions. A reader
skimming `vec![ref_allele(b"A"), ref_allele(b"T"), ref_allele(b"C")]`
has to stop and check whether the test is constructing a malformed
record. The in-module tests handle this with `fn alt_allele(seq)
{ ref_allele(seq) }`; same fix here.

#### Mi14: `record_encode.rs:340`, `writer.rs:132`, `tests/cohort_vcf_writer_integration.rs:25` — three-way fixture duplication (`ref_allele`, `support`, biallelic record factory)
**Confidence:** High.

Three test modules independently redefine the same fixtures. The
integration test already imports `AlleleSupportStats::new` because
the type is `#[non_exhaustive]` for out-of-crate callers, so the
fixtures are *already* spelled differently across the boundary
even though they mean the same thing. Hoist to a `#[cfg(test)]
pub(crate) mod test_fixtures` (or `tests/common/mod.rs`) and import
from all three sites. Next field addition becomes a one-site edit.

#### Mi15: `src/var_calling/vcf_writer/writer.rs:62, 88, 116` — `pub fn` methods on `CohortVcfWriter` lack `# Errors` rustdoc sections
**Confidence:** High.

`new`, `write_record`, and `finish` all return `Result<_,
VcfWriteError>` but enumerate no error variants in their docs.
The writer is the boundary the CLI talks to — callers should
not have to read the enum source to know what they may need to
handle. Add `# Errors` blocks to each. (Pair this work with M3:
once the enum is `#[non_exhaustive]`, the `# Errors` block can
still list the *known* variants.)

#### Mi16: BGZF EOF magic duplicated across `sink.rs:152` and `tests/cohort_vcf_writer_integration.rs:270`
**Categories:** smells (primary), extras (convergent).
**Confidence:** High.

The 28-byte htslib EOF marker is hand-written twice. Both copies
match today. A future edit that touches one would silently
diverge. Expose once as `pub(crate) const BGZF_EOF: &[u8; 28]`
from `sink.rs` (already there as `const BGZF_EOF` inside the test
module — promote scope) and import from the integration test.

#### Mi17: `src/var_calling/vcf_writer/record_encode.rs:117` — `genotype_order(record.ploidy, record.alleles.len())` is rebuilt per record
**Categories:** extras (perf, cross-cat).
**Confidence:** Medium.

`build_format_keys` is cached per writer (`writer.rs:50-52`,
`record_encode.rs:46-52`), but `genotype_order` is called once per
record inside `encode`. Memoise per `(ploidy, n_alleles)` pair on
`CohortVcfWriter` to avoid the rebuild — same pattern Stage 5's
per-group merger already uses
([per_group_merger.rs:443-447](../../../src/var_calling/per_group_merger.rs#L443-L447)).
Pair with M15's bench so the change can be defended.

#### Mi18: missing malformed-input tests for `allele_to_string`'s UTF-8 path and `GenotypeIndexOutOfBounds`
**Confidence:** High.

The extras rule asks for a test asserting the *specific* error
variant per malformed-input path. Two paths are uncovered:

```rust
#[test]
fn invalid_utf8_allele_surfaces_encode_error() {
    let mut record = biallelic_two_samples();
    record.alleles[1].seq = vec![0xff, 0xfe]; // not UTF-8
    let err = encode(&record, &fixture_contigs(),
                     &WriterConfig::default(),
                     &build_format_keys(&WriterConfig::default()), 2).unwrap_err();
    assert!(matches!(err, VcfWriteError::Encode(msg) if msg.contains("not valid UTF-8")));
}

#[test]
fn out_of_range_best_genotype_surfaces_typed_error() {
    let mut record = biallelic_two_samples();
    record.best_genotype = vec![0, 99]; // n_genotypes = 3
    let err = encode(/* ... */, 2).unwrap_err();
    assert!(matches!(err, VcfWriteError::GenotypeIndexOutOfBounds {
        sample_idx: 1, got: 99, .. }));
}
```

(`ContigLengthOverflow` is covered by Mi4.)

#### Mi19: `src/var_calling/vcf_writer/mod.rs:25` — plan rustdoc link points at `https://example.invalid`
**Categories:** smells (primary; flagged at Major), reliability (Nit), extras (Nit). Synthesizer's call: Minor.
**Confidence:** High.

The module doc-comment ends with `Plan: [doc/devel/implementation_plans/cohort_vcf_writer.md](https://example.invalid).`.
The label looks like a repo-relative path but the resolved URL is
the IETF reserved placeholder `example.invalid`. Drop the URL or
fix the target:

```rust
//! Plan: `doc/devel/implementation_plans/cohort_vcf_writer.md`.
```

### Nits

A handful of items grouped here per the rubric:

- `record_encode.rs:329` — `let _ = write!(out, "{a}");` silently
  discards the result; either `unwrap()` with a comment naming
  the discarded `fmt::Error` (impossible against a `String` sink)
  or use the explicit `out.write_fmt(format_args!(...))` shape.
- `record_encode.rs:46-52` — `build_format_keys` allocates an
  intermediate `Vec<String>` before `.collect()`-ing into `Keys`;
  chain `std::iter::once` over the conditional `GP` for one
  less allocation.
- `writer.rs:140` + `record_encode.rs:348-354` — `fn alt_allele
  { ref_allele(seq) }` test aliases that are literal pass-throughs.
  Either delete or have them differ.
- `header.rs:61-65, 72-76` — closure-parameter type annotation
  on `parse()` is needed today; turbofish on `parse::<Other>()`
  would be more idiomatic.
- `refactor_safety` agent flagged three `matches!` patterns in
  tests using bare `..` instead of `field: _` (`record_encode.rs:605`,
  `record_encode.rs:619-623`, `writer.rs:253`); a rename sweep on
  `VcfWriteError` field names makes these caught at compile time.
  Grouped here.
- `errors.rs:14-16` — variant doc comment names
  `GenotypeIndexOutOfBounds` and `ContigLengthOverflow` but omits
  the equally-defensive `UnknownChromId` and `SampleCountMismatch`.

## 7. Out of scope observations

- `src/vcf_writer.rs` (legacy) — pending deletion in a separate
  cleanup pass, no new findings.
- `src/var_calling/posterior_engine/{,backends,interp,shape}.rs`,
  `benches/var_calling_perf.rs`, `examples/profile_posterior_engine.rs`
  — 10 `cargo fmt` failures all live here; this is parallel SIMD
  WIP on `perf/posterior-samply`. Out of scope for the writer
  review.
- Spec gap (not for this review): the contamination side-pass
  produces `ContaminationEstimates` but no INFO/FORMAT field on
  the cohort VCF currently exposes per-sample `c_s` — surfaces
  when the cohort CLI slice lands; the writer is contract-ready.

## 8. Missing tests to add now

Grouped by function under test. Names follow `function_returns_expected_on_condition`.

### `tally_called_alleles` (Blocker B1)

- `tally_called_alleles_errors_on_out_of_bounds_allele_index` —
  drives a `PosteriorRecord` whose `best_genotype[s]` decodes
  through `genotype_order(2, 3)[3] == [0, 2]` while
  `record.alleles.len() == 2`. After the B1 fix, asserts
  `VcfWriteError::AlleleIndexOutOfBounds`.

### `encode` / structural validation (M6)

- `encode_errors_on_empty_alleles` — `record.alleles = vec![]` →
  typed error (after M6 fix).
- `encode_errors_on_undersized_best_genotype` — same shape for
  `best_genotype.len() != n_samples`.
- `encode_errors_on_undersized_scalars` — `scalars.len() !=
  n_samples * n_alleles`.

### Depth overflow (M7, M8)

- `encode_errors_when_per_sample_depth_overflows_i32` — `support(u32::MAX)`
  twice on the same sample, expect `VcfWriteError::DepthOverflow`.

### Error messages (M5)

- `sample_count_mismatch_message_names_cohort_first` — renders
  the error and asserts `contains("cohort metadata names …")` /
  `contains("posterior arrays carry …")`.

### Writer durability (M10, M13, M14)

- `finish_renames_and_fsyncs_parent_dir` — smoke that the new
  parent-dir-fsync code is exercised on the happy path.
- `new_cleans_up_tmp_on_header_write_failure` — sink injection
  to fail the header write; asserts `<output>.tmp` absent after.
- `dropped_without_finish_warns_via_must_use` — compile-time
  warning verification via `compile_fail` doctest is the closest
  to programmatic. Acceptable to skip; the type-level annotation
  is the actual guarantee.

### Writer state (Mi1)

- `out_of_order_does_not_advance_last_locus` — write at 200,
  fail at 100, then fail again at 150 (against 200), then accept
  201 — asserts the latch shape.

### GQ ceiling drift (M9)

- `gq_writer_cap_is_at_least_engine_cap` — compile-time-ish
  assert against `PosteriorEngineConfig::default().max_gq_phred`.

### Header validation (Mi3, Mi4, Mi5)

- `contig_length_above_i32_max_rejected`.
- `header_omits_source_when_tool_string_empty`.
- `header_omits_commandline_when_command_line_empty`.
- `empty_contigs_accepted_produces_header_with_no_contig_lines`
  (or `_rejected`, per Mi3 decision).

### Malformed inputs (Mi18)

- `invalid_utf8_allele_surfaces_encode_error`.
- `out_of_range_best_genotype_surfaces_typed_error`.

### Sink suffix (Mi2)

- `uppercase_suffixes_select_bgzf_sink` (only if the
  case-insensitive direction is chosen).

### Property test (reliability serializer rule)

- `proptest_encode_record_round_trips_through_noodles_reader` —
  random valid `PosteriorRecord` (ploidy ∈ {1, 2}, n_alleles ∈
  {1..=4}, n_samples ∈ {1..=8}, posteriors rows summing to 1 ± eps),
  emit through `CohortVcfWriter`, read back with
  `noodles_vcf::io::Reader`, assert REF/ALT/AC/AN/per-sample-GT
  round-trip. Required by the reliability rule for serializers.

## 9. What's good

- `tests/cohort_vcf_writer_integration.rs:268-282` — the bgzf
  EOF-byte-pattern assertion is the right level of test for a
  wire-format contract; it caught the htslib EOF requirement
  precisely.
- [src/var_calling/vcf_writer/sink.rs:25-29](../../../src/var_calling/vcf_writer/sink.rs#L25-L29)
  — `SinkKind` as a two-variant enum with a `Write` impl that
  dispatches the branches keeps the writer's per-record path
  branch-free at the type level. Clean.
- [src/var_calling/vcf_writer/record_encode.rs:34-41](../../../src/var_calling/vcf_writer/record_encode.rs#L34-L41)
  — `QUAL_MAX` is the textbook "named constant + doc explaining
  the rationale + grep anchor" pattern; the GQ cap (M9) failing
  to match this is the bug, the QUAL_MAX shape is the model.
- [src/var_calling/vcf_writer/writer.rs:84-100](../../../src/var_calling/vcf_writer/writer.rs#L84-L100)
  — non-decreasing-locus check as the writer's last line of
  defence is a clean defensive layer that the upstream merger's
  out-of-order check already catches; cheap belt-and-braces.
- [src/var_calling/vcf_writer/writer.rs:116-120](../../../src/var_calling/vcf_writer/writer.rs#L116-L120)
  — `finish(self)` consumes `Self` so a forgotten `finish` shows
  as "no output file", not a half-written file or a `Drop`-side-
  effect surprise. Pair with M14 to lock the compile-time guard.

## 10. Commands to re-verify

After fixes, re-run the same set inside `./scripts/dev.sh`:

- `cargo fmt --all -- --check` — must pass on the in-scope files
  (out-of-scope hunks are pre-existing parallel WIP).
- `cargo clippy --lib --tests --all-features -- -D warnings`.
- `cargo test --lib var_calling::vcf_writer`.
- `cargo test --test cohort_vcf_writer_integration`.
- `cargo test --all-targets` — every other test target green.

New commands the review introduces:

- `cargo test --test cohort_vcf_writer_proptest` — if the
  proptest from §8 is added as a separate test file.
- `cargo bench --bench vcf_writer_perf -- --save-baseline pre-M15`
  / re-run after — once M15 lands.
