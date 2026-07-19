# ng — reference info: types & interfaces

*Status: architecture draft (2026-07-18), companion to the spec
[`../spec/reference_info.md`](../spec/reference_info.md) (the design and every *why*) and to the
shared arch docs [`module_layout.md`](module_layout.md) (the `src/ng/` tree) and
[`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary). This module is **foundational
infra, not a pipeline step**, so `ng_step_interfaces.md` does not sketch it — this doc is its
interface. Naming follows [`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md):
domain nouns for types, verbs for functions; the module is marker-agnostic, so `ssr` does not
appear. **Signatures are illustrative; the contract is the deliverable.** See the spec for the
reasoning behind every decision here.*

## Module home

One file, `src/ng/reference_info.rs`, with its `#[cfg(test)]` block beside it. It is
**non-step infrastructure** — like `ref_seq.rs` and `pileup/`, shared by many consumers and owned
by none — so it carries **no trait and no bake-off**: it is a set of functions over plain data
types, not a swappable interface (`module_layout.md` principle 1a: a unit with no competing
implementations is a file). It splits into `reference_info/` when the pass, the cache, and the
writer grow enough to warrant it, the same file→folder rule `ref_seq.rs` follows. Deps: it imports
`crate::fasta::{ContigEntry, ContigList}` read-only (for `contig_list()` and the deferred
reconciliation), `noodles_fasta`, `md-5`, `std` — and **no ng-step type** (spec §3.6). Why it
exists, and its adoption mandate (the one builder of a contig table from a reference, and the one
place ng parses a `.fai`): spec §1.

## 1. Types

### 1.1 The data — `ContigInfo`, `ReferenceInfo`

ng-native, plain data, `pub` fields (unconstrained). `ContigInfo` is the everything-about-a-contig
type production never had in one place — a superset of `ContigEntry` (which has no geometry) and of
a `.fai` row (which has no md5). why: spec §3.9.

```rust
/// One contig, everything in one place: name and length (as a BAM `@SQ` carries), the `.fai`
/// geometry that locates it in the file (§3.8 of the spec), and — once the bases are read — its
/// MD5 (the `@SQ M5`). Named for what it *is*, not the `.fai` some fields come from.
pub struct ContigInfo {
    pub name: String,
    pub length: u64,
    pub offset: u64,          // faidx.5 OFFSET — byte offset of the first base
    pub line_bases: u64,      // faidx.5 LINEBASES — bases per line
    pub line_width: u64,      // faidx.5 LINEWIDTH — bytes per line, incl. the newline
    pub md5: Option<[u8; 16]>, // the @SQ M5; None until the FASTA is read (a .fai has no MD5 column)
}

/// What reading a reference yields: the whole-assembly digest and every contig's info.
pub struct ReferenceInfo {
    /// Whole-reference MD5 — every contig's uppercased bases *concatenated* in file order
    /// (the `ssr-catalog` convention; NOT a SAM field). Distinct from the per-contig `@SQ M5`
    /// on each `ContigInfo`. `None` from a `.fai`-only read. why: spec §3.4.
    pub md5: Option<[u8; 16]>,
    /// Every contig, **in file order** — the order that defines `ContigId` (`ContigId(i)` is
    /// `contigs[i]`). One source of truth for name/length/geometry/md5. why: spec §3.9.
    pub contigs: Vec<ContigInfo>,
}

impl ReferenceInfo {
    /// **Transitional** projection to the production `ContigList` (one `ContigEntry {name,
    /// length, md5}` per contig) for ng consumers that still speak it (`RefSeq`/`ContigTable`,
    /// `GenomeRegions`); geometry is dropped. Removed when those migrate to `ContigInfo`.
    /// The only place `ContigList` surfaces from this module. why: spec §3.6, §3.9.
    pub fn contig_list(&self) -> ContigList;
}
```

**Contract.** `contigs` is in FASTA/`.fai` file order and that order is the `ContigId` contract —
contigs are never reordered. `reference_md5` and each `ContigInfo.md5` are `Some` iff the FASTA
bases were read (both FASTA sources); the `.fai`-only source leaves them `None` — an absent MD5 is
a wildcard downstream (`ContigEntry` `PartialEq`, spec §5 T6). `contig_list()` is a pure projection,
so its `ContigEntry`s cannot disagree with `contigs`.

### 1.2 The input — `ReferenceSource`

Two ways to name a reference; **owns its paths** (`PathBuf`, no lifetime — the module stays
lifetime-free through the cache). why: spec §3.1.

```rust
/// What the caller hands in. The variant names the cost: `Fai` is a small file read, `Fasta`
/// reads the whole reference. why: spec §3.1, §3.2.
pub enum ReferenceSource {
    /// The index alone → names, order, lengths, geometry; `md5: None`. Cheap.
    Fai(PathBuf),
    /// The FASTA, optionally cross-checked against a `.fai`. `fai: None` → read alone (all MD5s).
    /// `fai: Some` → read the same **and** verify against that index (§3.3); disagreement errors.
    /// Both cost a whole-genome read — the index only adds a comparison.
    Fasta { fasta: PathBuf, fai: Option<PathBuf> },
}
```

**Contract.** A supplied `.fai` is *always* checked — there is no trust-the-index fast path
(spec §3.2). The module never probes the filesystem for a path it was not handed (spec §3.6); the
sibling naming convention lives only in `sibling_fai_path` (§2) and the orchestrator (§2, §3.11).

### 1.3 The error — `ReferenceInfoError`

`#[non_exhaustive]` `thiserror` enum (the house pattern, mirroring `RefSeqError`); each variant
carries the path. why: spec §2, §5.

```rust
#[non_exhaustive]
pub enum ReferenceInfoError {
    FaiRead { .. },            // noodles `fai::fs::read` failed (missing/unreadable .fai)
    FastaRead { .. },          // FASTA I/O failed mid-pass
    MalformedFasta { .. },     // bad geometry / bases-before-`>` / empty contig / non-uniform line (T3)
    DuplicateContigName { .. },// a name appears twice — stricter than htslib, which warns+drops (T2)
    FastaFaiMismatch { .. },   // fasta-vs-.fai disagreement, naming the field + contig (§3.3)
    FastqIndex { .. },         // a 6-column FASTQ index handed where a FASTA one was expected (§3.8)
    CompressedReference { .. },// a bgzip .fa.gz — out of scope, deferred (§3.8)
    FaiWrite { .. },           // the write step of `read_reference_verifying_or_creating_fai` failed (§3.11)
}
```

**Contract.** Nothing in a supplied reference may panic — its bytes are attacker-influenced input,
so every malformed shape is a named error, never an `assert!`/`debug_assert!` (spec §5 T3).

### 1.4 Background verification — `VerificationHandle`

```rust
/// The pending result of a background FASTA verification (§2's background entry point).
/// `#[must_use]`, and **must be joined before the caller commits any output** — an unverified
/// reference is a possibly-wrong run, and this is the only place that error surfaces. why: spec §3.10.
#[must_use]
pub struct VerificationHandle { /* crossbeam result channel + the worker's join handle */ }

impl VerificationHandle {
    /// Non-blocking: `true` once verification has finished; lets a caller defer `join`.
    pub fn is_finished(&self) -> bool;
    /// Block, take the result: the **verified** `ReferenceInfo` (MD5s filled) or the error.
    pub fn join(self) -> Result<Arc<ReferenceInfo>, ReferenceInfoError>;
}
// Drop without join() logs a warning to stderr (house convention) rather than swallow a
// possibly-pending failure; it does not block. why: spec §3.10.
```

**Contract.** The md5s arrive by **return** (`join`), never by in-place mutation — `ContigInfo`
stays a plain immutable struct (spec §3.10). Before `join`: `md5: None`; after: the returned value
has `md5: Some`. On failure the handle yields the error and the cache slot is left empty
(errors are not cached — §1.5, spec §3.7/§3.10).

### 1.5 The cache — `ReferenceInfoCache`

```rust
/// A thread-safe, compute-once, single-flight cache over `read_reference_info`. Caller-held
/// (not a global static): the caller constructs one per run, shares it, and every thread asking
/// for a reference already being read **waits** rather than re-reading. `Send + Sync`. why: spec §3.7.
pub struct ReferenceInfoCache { /* Mutex<HashMap<Key, Arc<Mutex<Option<Arc<ReferenceInfo>>>>>> */ }

impl ReferenceInfoCache {
    pub fn new() -> Self;
    /// `Arc`-shared result (a hit is a pointer clone). Keyed on `(source discriminant, per-file
    /// (path, size, mtime))` — `Fai`, `Fasta{None}`, `Fasta{Some}` are distinct keys.
    pub fn get_or_read(&self, source: ReferenceSource)
        -> Result<Arc<ReferenceInfo>, ReferenceInfoError>;
}
```

**Contract.** Two-level lock: the map lock is held only to get-or-create the per-key slot, never
across a read, so a read on one key never blocks a read on another; the slot lock is held across
the read, which *is* the single-flight (same-key callers wait, then hit). **Successes only** are
cached — a transient I/O error is not, so it can retry (spec §3.7, chosen over `OnceLock::get_or_init`
for exactly this). Key is stat-based, which is sound because a reference is write-once during a run
(spec §5 T8). If mtime is unavailable, `get_or_read` bypasses the cache rather than failing — the
cache is an optimisation, never a correctness gate.

## 2. Interfaces — the functions

Two layers: **pure primitives** (explicit paths, no filesystem probing, no side effects) and, on
top, one **orchestrator** that adopts the sibling convention and may write. why: spec §3.6.

```rust
// ---- pure primitives ----

/// Read a reference. Pure: no shared state, no probing beyond opening the paths given; every
/// call does the full work the source names. why: spec §3.1.
pub fn read_reference_info(source: ReferenceSource)
    -> Result<ReferenceInfo, ReferenceInfoError>;

/// Write a `.fai` from contig info — the five-column `faidx.5` form (`md5` ignored), **byte-
/// identical to `samtools faidx`**, written **atomically** (`.tmp` + rename). why: spec §3.9.
pub fn write_fai(contigs: &[ContigInfo], out: &Path) -> io::Result<()>;

/// `<fasta>` + ".fai" — the reference naming convention, as a path, **no I/O**. A caller wanting
/// sibling behaviour composes this with its own existence check. why: spec §2, §3.6.
pub fn sibling_fai_path(fasta: &Path) -> PathBuf;

// ---- background verification (through the cache) ----

/// Read the `.fai` **now** (`cache.get_or_read(Fai)` — fast, `md5: None`, returned immediately)
/// and verify the FASTA on a **background thread** (`cache.get_or_read(Fasta{Some})`). Both go
/// through the cache, so the verify **populates** it and a concurrent reader of the same key
/// coordinates on the single-flight slot instead of re-reading. Takes `&Arc<ReferenceInfoCache>`
/// because the detached thread outlives the call. why: spec §3.10.
pub fn read_fai_verify_in_background(
    cache: &Arc<ReferenceInfoCache>, fasta: PathBuf, fai: PathBuf,
) -> Result<(Arc<ReferenceInfo>, VerificationHandle), ReferenceInfoError>;

// ---- the batteries-included orchestrator ----

/// Read a reference the convenient way: derive `<fasta>.fai`; if present → verify in the
/// background (returns immediately, `Some(handle)` to join); if absent → scan now (verified,
/// MD5s), **write** the `.fai`, return `(info, None)`. The *one* place the module adopts the
/// sibling convention and the *one* place a read writes a file; a `.fai`-write failure is
/// **fatal** (read-only-dir escape hatch: call `get_or_read(Fasta{None})` directly). why: spec §3.11.
pub fn read_reference_verifying_or_creating_fai(
    cache: &Arc<ReferenceInfoCache>, fasta: PathBuf,
) -> Result<(Arc<ReferenceInfo>, Option<VerificationHandle>), ReferenceInfoError>;
```

**Contract on the orchestrator's asymmetric return.** `Some(handle)` = verification pending (fai
present, join to verify); `None` = nothing to await (scan path, already verified). A caller wanting
uniformly-verified info writes `if let Some(h) = handle { h.join()?; }`. why: spec §3.11.

## 3. The FASTA pass (contract only — body is the implementer's)

The `Fasta` sources run one **from-byte-zero, one-buffer streaming** loop that produces every
`ContigInfo` field and both MD5s in a single pass; the geometry it reconstructs is then compared to
the `.fai` (`Fasta{Some}`). The single predicate: a byte is a **base** iff it is in `[0x21, 0x7E]`
(printable, non-space) — the same predicate defines the `@SQ M5`, `line_bases`, and `length`. Full
loop, the CR-LF handling, and the four rejected readers (including *why not reuse ng's own
fai-driven streaming readers* — circular, spec §3.3): spec §4. **Do not pull the loop body into
this doc.**

## 4. Design decisions — decided

Crisp records; the *why* is the spec, cross-referenced. Open items marked `OPEN:`.

- **ng-native types, not reused production ones — decided.** `ContigInfo`/`ReferenceInfo` are ng's;
  ng codes what is clearest for ng and does not carry production types for port-back's sake.
  `contig_list()` is a *transitional* bridge for consumers still on `ContigList`, not a
  compatibility promise. Rejected: storing a production `ContigList` + a parallel `Vec<FaiRecord>`
  (the earlier draft) — name/length duplicated, a projection to keep in step; the merge removes it
  (spec §3.9, §3.6).
- **Two sources; a supplied `.fai` is always checked — decided.** No `verify: bool`, no
  trust-the-index path; passing a `.fai` *is* asking for the check. The cost is legible in the
  variant the caller chose (spec §3.1, §3.2).
- **"They match" means the whole `.fai`, geometry included — decided.** The comparison
  reconstructs and checks all five `faidx.5` columns, not just name+length, because a re-wrap
  keeps name+length and breaks the geometry the readers seek by. Rejected: `first_disagreement`
  (no geometry) and hashing-through-the-`.fai`-offsets (circular) (spec §3.3).
- **M5 rule = `samtools dict`'s: keep `[0x21,0x7E]`, uppercase, MD5 — decided.** One predicate for
  the M5 and the `.fai` lengths. This *diverges from* production's `compute_contig_md5_streaming`
  (which hashes spaces/tabs) deliberately — ours must equal a BAM's `@SQ M5` (spec §3.4).
- **Reject duplicate contig names — decided.** Stricter than htslib (warns + drops): ng resolves
  contigs by position, so a dropped duplicate renumbers every `ContigId` (spec §5 T2).
- **Caller-held, single-flight cache — decided.** Not a global static (hidden state, test leakage,
  process-wide stale hits). Errors not cached. Rejected: `OnceLock::get_or_init` (caches the
  failure) (spec §3.7).
- **`ContigInfo` carries the index; ship a `.fai` writer — decided.** Reverses spec §8's earlier
  "no", whose objection was a *noodles* type in the output; a plain `ContigInfo` removes it and the
  data is already computed. Real use: index a FASTA that lacks a `.fai` (spec §3.9).
- **Background verify: crossbeam + return, through the cache — decided.** Not async (no runtime in
  the tree); not `rayon::spawn` (a long blocking read starves the pool — rayon's place is *inside*,
  the deferred parallel hash). md5s arrive by **return via `join`**, not in-place mutation
  (shared-mutable-across-threads → locked fields → nondeterministic-until-join). Reads **through**
  the cache so the verify populates it and a concurrent reader coordinates (spec §3.10).
- **`read_reference_verifying_or_creating_fai` is the one sanctioned discovery+write point —
  decided.** A named orchestrator on the far side of §3.6's modularity line, so the pure primitives
  stay pure. Write failure is **fatal** (owner: fully-set-up-or-nothing); read-only-dir callers use
  the primitives (spec §3.11).

### Deferred, with a home (spec §7)

- **`@SQ` ↔ `.fai` reconciliation** (closes `ReadFilter::new`'s permutation hole) → read
  ingestion's spec; reuses `ContigList::first_disagreement`, compares *order* not just membership.
- **Consolidating ng's fai-driven readers onto `ContigInfo`** → `RawChromReader` takes a
  `&ContigInfo` instead of re-parsing the `.fai` per contig transition; retires the `contig_list()`
  bridge for the `RefSeq` impls. ng-internal.
- **bgzip-compressed reference** (`.fa.gz` + `.fai` + `.gzi`) → error today; a bgzf reader over §3's
  pass when a consumer needs it.
- **Parallelising a single read** (the two-phase hash) → measure-first; the cache is the cheaper win.

## 5. Reconciliation with existing code

Convergence points, not new types to invent alongside the old. Vendored oracles
(`samtools`/`htslib`) are gitignored, cited in plain text (no link). Verified against the code.

| ng name / piece | existing code | action |
|---|---|---|
| `ContigInfo` | — (supersets `ContigEntry` + a `.fai` row) | **new** ng type |
| `ReferenceInfo`, `ReferenceSource`, `ReferenceInfoError`, `VerificationHandle`, `ReferenceInfoCache` | — | **new** |
| `ContigList` / `ContigEntry` (+ wildcard-md5 `PartialEq`) | [`src/fasta/mod.rs:37-62`](../../../../src/fasta/mod.rs) (`:43-55`) | **reuse read-only** — `contig_list()` builds them |
| `first_disagreement` | [`src/fasta/mod.rs:69-100`](../../../../src/fasta/mod.rs) (`pub(crate)`) | reuse later for `@SQ` reconciliation (§4 deferred), not here |
| `.fai` parse (`Fai` arm) | `noodles_fasta::fai::fs::read` ([`Cargo.toml:93`](../../../../Cargo.toml)) | call as-is |
| FASTA-pass loop (§3) | `noodles_fasta::io::Indexer::index_record` | **copy the logic + fold MD5**; do *not* call it (it hides the bases) |
| `.fai` field guards (`line_bases>0`, `line_width≥line_bases`) | `ContigFai::validate` ([`src/ng/raw_chrom_reader.rs:77`](../../../../src/ng/raw_chrom_reader.rs)) | copy the guards |
| streaming MD5 (64 KiB buffer) | `compute_contig_md5_streaming` ([`src/pop_var_caller/common.rs:250`](../../../../src/pop_var_caller/common.rs)) | copy the streaming shape — but predicate is `isgraph`, not skip-`\n\r` (§3.4) |
| MD5 → hex | `format_md5_hex` ([`src/pop_var_caller/common.rs:60`](../../../../src/pop_var_caller/common.rs), `pub(crate)`) | reuse |
| M5 / `.fai` oracle | `samtools/dict.c:78-86`, `htslib/faidx.c:48-50`, `htslib/faidx.5` (vendored) | match byte-for-byte; test via `samtools dict`/`faidx` |
| readers to consolidate onto `ContigInfo` | `open_contig`/`ContigFai` ([`src/ng/raw_chrom_reader.rs:58`,`:111`](../../../../src/ng/raw_chrom_reader.rs)) | consolidation target (§4 deferred) |
| consumers migrating off `ContigList` | `ContigTable` / `WindowedRefSeq::new` ([`src/ng/ref_seq.rs:212`,`:524`](../../../../src/ng/ref_seq.rs)) | `contig_list()` bridges until they take `ContigInfo` |
| deps | `md-5 = "0.11"` ([`Cargo.toml:81`](../../../../Cargo.toml)), `noodles-fasta = "0.61.0"` (`:93`) | in tree already |

## 6. Open items

- `OPEN:` **Does noodles' buffered reader expose a clean per-record byte offset?** §3's pass needs
  the offset reached after each definition line; a counting `BufRead` gives it *if* the reader
  never reads ahead past a record boundary. Confirm at implementation time; fallback is to parse
  the definition line with `BufRead::read_until` directly (spec §8). Not a design decision.
- Concrete `ReferenceInfoError` variant fields/`#[from]` wiring — pin when coding; the *set* of
  variants is decided (§1.3).

## 7. Test & bench shape (spec §6)

- Tests in the module's own `#[cfg(test)]`. **First fixture to build:** a FASTA writer that wraps
  at a chosen width (the existing `ref_seq.rs::build_fasta` is one-line-per-contig, so it cannot
  exercise geometry) — everything else stands on it.
- **Oracle is `samtools`, resolved to committed constants** (no `samtools` at test time): run
  `samtools dict` + `samtools faidx` once on the golden reference
  (`tests/data/tandem_repeat/synthetic_ref.fa`), commit the `M5`/`LN`/`.fai` values, assert we
  reproduce them; `write_fai` output byte-identical to the committed `.fai`. The `faidx.5` LF/CR-LF
  worked example is a second committed vector.
- **Regression anchors** — the tests the feature exists for: the `fold -w`-defeats-it re-wrap must
  error on `Fasta{Some}` naming `line_bases` (mutation-verified against the cheaper names-only
  check); the background-verify single-flight-through-cache (real read count = 1, concurrent reader
  waits); the orchestrator's fatal write + escape hatch.
- **No `bench/`** — no bake-off (infra, single implementation).

---

*Cross-links: spec [`../spec/reference_info.md`](../spec/reference_info.md); build order
[`../impl_plan/reference_info.md`](../impl_plan/reference_info.md); tree
[`module_layout.md`](module_layout.md); indexed in [`../README.md`](../README.md).*
