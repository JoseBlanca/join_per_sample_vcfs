# ng step 1 — read filtering: types & interfaces

*Status: architecture draft (2026-07-14), companion to the spec
[`../spec/read_filtering.md`](../spec/read_filtering.md) (the design and its rationale) and
to the shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary + step
traits) and [`module_layout.md`](module_layout.md) (the `src/ng/` tree). `ng_step_interfaces.md`
§3 sketches step 1 only briefly and defers "the real interface"; this doc is it. Naming
follows [`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): domain
nouns for types, verbs for functions, newtypes for domain scalars, **STR** in prose ↔ `ssr`
in code (read filtering is generic, so `ssr` does not appear). Signatures are illustrative;
the **contract** is the deliverable. See the spec for the "why" behind every decision here.*

## Module home

One file, `src/ng/read/filtering.rs`, with its `#[cfg(test)]` block beside it — read
filtering is a **fixed prelude with no bake-off**, so it is a file, not a folder
(`module_layout.md` principle 1a). It shares the `read/` module with step 2 (`ReadPrep`),
which turns filtering's output into locus evidence and keeps its swappable impls side by side
(principle 1b). The scalar newtypes seed `src/ng/types.rs` (§2.2 below); everything else is
step-1-local.

## 1. The filter cascade (contract)

Nine filters in **hit-rate order** (cheapest, most-often-firing first), ported verbatim from
production (`per_sample_caller_cram_input.md`). A read is charged to the **first** filter it
fails; the keep/drop outcome itself is order-independent. Full rationale per filter: spec §3.

| # | filter | drops when | default | toggle | phase |
|---|---|---|---|---|---|
| 1 | duplicate | `flag & FLAG_DUPLICATE` | on | `drop_duplicate` | pre-decode |
| 2 | low MAPQ | `mapq < min_mapq` (0xFF → 0) | `Some(20)` | `min_mapq` | pre-decode |
| 3 | supplementary | `flag & FLAG_SUPPLEMENTARY` | on | — | pre-decode |
| 4 | secondary | `flag & FLAG_SECONDARY` | on | — | pre-decode |
| 5 | unmapped | `flag & FLAG_UNMAPPED` | on | — | pre-decode |
| 6 | QC fail | `flag & FLAG_QC_FAIL` | on | `drop_qc_fail` | pre-decode |
| 7 | too short | `len < min_read_length` | `Some(30)` | `min_read_length` | post-decode |
| 8 | high mismatch fraction | quality-clearing mismatches / comparable bases `> max` (reference-dependent) | `Some(0.10)` | `max_read_mismatch_fraction` | post-decode |
| 9 | bad CIGAR | adjacent indel pair, or boundary deletion | on | — | post-decode |

**The decode boundary sits between #6 and #7** (spec §3): #1–#6 read only flag/MAPQ, both
available before decode; #7–#9 need the decoded read. Filters #3–#5 have **no toggle** by
design — admitting them would double-count alleles or break per-read invariants. Filter #8 is
the only reference-dependent one; it turns off cleanly when `max_read_mismatch_fraction` is
`None`.

## 2. Types

### 2.1 Scalar newtypes — seed `types.rs`

Only the newtypes read filtering touches; names match `ng_step_interfaces.md` §1 so nothing
renames later. Conventions (spec §2.3): unconstrained newtypes keep `pub` fields and a
`.get()`; a constrained newtype hides its field behind a checked constructor.

```rust
pub struct MapQual(pub u8);          // SAM MAPQ; 0xFF (unavailable) → treated as 0
pub struct BaseQual(pub u8);         // Phred base quality, 0–93
pub struct Bp(pub u32);              // length in base pairs (generic length currency)
pub struct MismatchFraction(f32);    // constrained to [0, 1] — checked constructor

impl MismatchFraction {              // reject an out-of-range threshold loudly (untrusted config)
    pub fn try_new(x: f32) -> Result<Self, DomainError>;
    pub fn get(self) -> f32;
}
pub enum DomainError { MismatchFraction(f32), /* later constrained types add variants */ }
```

`DomainError` is introduced here (ng-wide domain-invariant error) with its first variant.

### 2.2 Step-1-local types

```rust
/// The filtering policy — one field per active filter, no dormant levers. `Default` is the
/// production policy. `Option<T>` = "no threshold" (never a sentinel; `Some(0)` ≠ `None`).
pub struct ReadFilterConfig {
    pub min_mapq: Option<MapQual>,
    pub min_read_length: Option<Bp>,
    pub drop_qc_fail: bool,
    pub drop_duplicate: bool,
    pub max_read_mismatch_fraction: Option<MismatchFraction>,   // None = filter #8 off
    pub mismatch_bq_floor: BaseQual,                            // only meaningful when the above is Some
}

/// One read's verdict. Drop variants line up 1:1 with the ReadFilterCounts fields.
pub enum FilterVerdict { Keep, Drop(DropReason) }
pub enum DropReason {
    Duplicate, LowMapq, Supplementary, Secondary, Unmapped,
    QcFail, TooShort, HighMismatchFraction, BadCigar,
}

/// Running per-sample tally — one counter per drop reason plus the kept count. The ng port
/// of `FilterCounts`. "No silent caps": every vanished read is accounted for.
pub struct ReadFilterCounts {
    pub kept: u64,
    pub duplicate: u64, pub low_mapq: u64, pub supplementary: u64, pub secondary: u64,
    pub unmapped: u64, pub qc_fail: u64, pub too_short: u64,
    pub high_mismatch_fraction: u64, pub bad_cigar: u64,
}
```

No new *read* type: a kept read is a bare `MappedRead` (step 2's input), not a `FilteredRead`
wrapper (§4).

## 3. Interfaces

Two traits describe the **input edge**; the filter is a std iterator on the **output edge**.
The input is a `RecordSource` that fills one reused buffer — deliberately *not* an
`Iterator<Item = RawRecord>` (§4). Both traits are ng's; their production impls are ng-owned
adapters over the noodles reader/record (§4).

```rust
/// A borrowed view of one alignment record — the seam that runs the flag/MAPQ cascade
/// BEFORE decode. `decode` borrows `&self`, so the underlying buffer stays reusable.
pub trait RawRecord {
    fn flag(&self) -> Flags;         // #1, #3–#6   (Flags = whatever MappedRead.flag carries)
    fn mapq(&self) -> MapQual;       // #2
    fn decode(&self) -> MappedRead;  // expensive phase (#7–#9 read the result); copies what MappedRead keeps
}

/// The filter's input: fills a caller-owned buffer with the next record, reusing its
/// allocations. Replaces `Iterator<Item = RawRecord>` precisely to get that reuse.
pub trait RecordSource {
    type Record: RawRecord + Default;
    fn read_next(&mut self, buf: &mut Self::Record) -> io::Result<bool>;   // Ok(false) = EOF; Err is fatal
}

/// A lazy iterator of kept reads, carrying its running counts. One RecordBuf is reused
/// across the whole pass; `next()` reads a record, runs #1–#6 on flag/MAPQ, decodes only on
/// survival, runs #7–#9, tallies every drop, and yields the first Keep.
pub struct ReadFilter<S: RecordSource, R /* : RawRefSeq */> {
    source: S,
    buf: S::Record,               // the single reused record buffer
    reference: R,                 // RawRefSeq — raw bytes for #8 (ref_seq.md)
    config: ReadFilterConfig,
    counts: ReadFilterCounts,
}

impl<S: RecordSource, R: RawRefSeq> Iterator for ReadFilter<S, R> {
    type Item = MappedRead;       // survivors, already decoded
    fn next(&mut self) -> Option<MappedRead>;
}

impl<S: RecordSource, R: RawRefSeq> ReadFilter<S, R> {
    /// Fail-fast setup: validates every source-header contig resolves in the reference.
    pub fn new(source: S, reference: R, config: ReadFilterConfig) -> Result<Self, RefSeqError>;
    pub fn counts(&self) -> &ReadFilterCounts;   // running tally; final once iteration is exhausted
}

// The decision, split on the decode boundary so each half unit-tests in isolation:
fn verdict_pre_decode(flag: Flags, mapq: MapQual, config: &ReadFilterConfig) -> FilterVerdict;
fn verdict_post_decode(read: &MappedRead, reference: &impl RawRefSeq,
                       config: &ReadFilterConfig) -> FilterVerdict;
```

**Contract:** content-preserving, per-sample, lazy; every dropped read charged to exactly one
reason. Decode runs only on pre-decode survivors. One buffer is reused; because a kept
`MappedRead` outlives it, `decode` **copies** the bytes it keeps (per-kept-read — reads
dropped pre-decode pay neither copy nor allocation). `reference` is consulted only by #8, only
when `max_read_mismatch_fraction` is `Some`. Any `RefSeqError` (#8) or `read_next` `Err` is
**fatal to the run**, not a per-read drop — which keeps the surface `Item = MappedRead` (no
per-item `Result`, no error bucket in the counts).

## 4. Design decisions — decided

Distilled from the spec; see it for the reasoning. Add new open items with `OPEN:`.

- **Filtering owns its entire policy — decided.** The reader is asked only to hand over
  records and decode on demand; the full cascade lives in one visible `ReadFilterConfig`, not
  inherited from a hidden reader default. A real, inspectable step means its policy is in one
  place (spec §2.5).
- **Pre-decode gate lives in the filter, decode only on survivors — decided.** Production's
  `classify_pre_decode` optimisation, relocated from the reader into read filtering so the
  whole policy and tally stay in one place. Result-preserving (same drops, byte-identical
  output); a read dropped on flag/MAPQ never pays decode cost (spec §3).
- **Input is a reused-buffer `RecordSource`, not an iterator of records — decided.** A std
  `Iterator`'s owned `Item` can't borrow a reused buffer (the lending-iterator problem), so an
  `Iterator<Item = RawRecord>` would force a fresh `RecordBuf` per read. A source that fills
  one buffer keeps reuse; the *output* stays `Iterator<Item = MappedRead>` so nothing
  downstream changes (spec §5).
- **`decode(&self)`, copy-on-keep — decided.** The output outlives the reused input buffer, so
  `decode` copies the bytes `MappedRead` keeps. Reuse-buffer and move-out-of-buffer are
  mutually exclusive (moving drains the buffer, defeating reuse); copy-on-keep is the
  deliberate choice, and it only costs on *kept* reads (spec §5).
- **No `FilteredRead` typestate — decided.** Kept reads stay bare `MappedRead`s. A wrapper
  would make "reached step 2 unfiltered" a compile error, but that per-read cost is not worth
  a guarantee the test suite already secures cheaply; composability with step 2's
  `prepare_read(&MappedRead, …)` stays adapter-free (spec §4, §7).
- **Error model is fatal / run-level — decided.** `ReadFilter::new` validates the source
  header against the reference up front (returns `Result`); in-loop fetch/read errors abort
  the run rather than degrade to a per-read drop or poison every `next()` with a `Result`
  (spec §5, §7).
- **`RecordSource`/`RawRecord` impl = ng-owned adapter — decided.** The traits are ng's; the
  production impl is an ng-owned adapter wrapping the noodles reader/record, so the dependency
  points ng → existing code and production never learns about ng. Rejected: adding the impl to
  the reader (would couple production to an ng trait for no gain) (spec §7).
- **No `bench/` this step — decided.** Read filtering has no competing implementations, so no
  frontier to plot; its measurement is the fixture-driven regression test (spec §2.6). First
  real `bench/` use is candidate generation (`ng_proposal.md` §3).

### Deferred, with a recommended home (spec §6)

- **Adaptor-clip *application*** and **mate-overlap quality reconciliation** → `pileup/`
  (generic) / pair-HMM prep (STR). Per-base operations needing pileup context — not per-read
  decisions. Read filtering *preserves* the `MappedRead.adaptor_boundary` annotation but does
  not apply it.
- **Coverage downsampling / read pooling** → future `ReadFilterConfig` knobs, added only when
  they enter the pipeline (not dormant fields now).
- **Output-`MappedRead` allocation reuse** → locus-stream memory design (a consumer-side
  free-list or per-window arena), not step 1. Buys lower allocator churn, not a lower peak;
  measure-first.

## 5. Reconciliation with existing code

The reuse map (spec §2.5) in the `ng_step_interfaces.md` §6 format — consolidation points, not
new types to invent alongside the old. Verify against the code when implementing.

| ng name | existing code | action |
|---|---|---|
| `ReadFilterConfig` | filtering subset of `AlignmentMergedReaderConfig` ([bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)) | new; mirrors the subset |
| `ReadFilterCounts` | `FilterCounts` ([bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)) | port / rename |
| `RawRecord`, `RecordSource` (traits) | the reader's record-yield + `classify_pre_decode`'s `RecordBuf` input | **new**; ng-owned adapters over the noodles reader/record |
| `ReadFilter` (driver) | the reader's filter-application plumbing | **new** driver; reuses the pure predicates below |
| filter #8 test | `read_exceeds_mismatch_fraction(...)` | call directly — gated on `RawRefSeq` (raw bytes, matching `RawContigRefCache`) |
| filter #9 test | `cigar_is_bad(...)` | call directly (pure — no reference) |
| flag constants | `FLAG_DUPLICATE`, `FLAG_SECONDARY`, `FLAG_SUPPLEMENTARY`, `FLAG_UNMAPPED`, `FLAG_QC_FAIL` | import as-is |
| default thresholds | `DEFAULT_MIN_MAPQ` (20), `DEFAULT_MIN_READ_LENGTH` (30), `DEFAULT_MAX_READ_MISMATCH_FRACTION` (0.10), `DEFAULT_MISMATCH_BQ_FLOOR` (10) | import as-is |
| `MapQual` / `BaseQual` / `Bp` / `MismatchFraction` | scalars scattered as raw `u8`/`u32`/`f32` | seed `types.rs` (match `ng_step_interfaces.md` §1) |
| `MappedRead` | `MappedRead` ([bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)) | reuse as-is (also §6 there) — the step-2 input |
| `RawRefSeq` | see [`../spec/ref_seq.md`](../spec/ref_seq.md) / `ng_step_interfaces.md` §6 | reuse (#8's `fetch_raw` capability) |

## 6. Open items

- **`Flags` concrete type** — the `RawRecord::flag` return resolves to whatever
  `MappedRead.flag` carries; confirm at implementation time. Not a design decision.
- **`SetupError` vs `RefSeqError`** — `new` returns `Result<Self, RefSeqError>` on a
  header/reference mismatch; if setup grows other failure modes, a dedicated setup-error type
  may be cleaner. Revisit only if a second failure mode appears.

## Test & bench shape (spec §2.6)

- Unit tests beside the code: each filter's boundary (a read exactly at `min_mapq` kept, one
  below dropped) and the cascade order; the `verdict_pre_decode` / `verdict_post_decode` split
  makes each half testable in isolation (post-decode half against the in-memory `RawRefSeq`).
- One fixture-driven integration test: filter a small known BAM/CRAM, assert the exact
  `ReadFilterCounts` against hand-counted expectations — the regression anchor for the port.
- No `bench/` (no bake-off) — see §4.
