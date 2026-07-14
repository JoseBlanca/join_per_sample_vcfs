# ng read filtering (step 1) — implementation plan

**Status:** draft, 2026-07-14. The build order for **step 1, read filtering**: the `read/`
module, the vocabulary it adds to `types.rs`, the two-phase filter cascade, the
`RecordSource`/`RawRecord` seam, and the `ReadFilter` iterator. Design is settled in
[`read_filtering.md`](../spec/read_filtering.md) (spec) and
[`../arch/read_filtering.md`](../arch/read_filtering.md) (types & interfaces), under the shared
arch docs ([step interfaces](../arch/ng_step_interfaces.md),
[module layout](../arch/module_layout.md)). This roadmap turns that design into build order; it
is **not** a place for new design — every open question is resolved in the spec §7.

This is the plan [`foundations.md`](foundations.md) deferred ("read filtering follows in its own
plan once its design settles").

---

## Scope

**In:** `src/ng/read/mod.rs` + `src/ng/read/filtering.rs`; the scalar newtypes read filtering
adds to `types.rs`; `ReadFilterConfig` / `FilterVerdict` / `DropReason` / `ReadFilterCounts`;
the `verdict_pre_decode` / `verdict_post_decode` cascade; the `RawRecord` + `RecordSource`
traits with a test fake **and** the ng-owned noodles adapter; the `ReadFilter` iterator.

**Out (later plans):**

- **Step 2, read preparation (`ReadPrep`)** — the `read/` sibling; its own plan.
- **Adaptor-clip *application*, mate-overlap reconciliation** — per-base work that lands in
  `pileup/` (spec §6), not here. Filtering only *preserves* `MappedRead.adaptor_boundary`.
- **Output-`MappedRead` allocation reuse** — a locus-stream memory decision (spec §6).
- **Downsampling / read pooling** — future `ReadFilterConfig` knobs, added when they enter the
  pipeline.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The algorithmic heart before the plumbing.** The two-phase `verdict_*` cascade (Milestone
  B) is built and unit-tested in isolation — against the in-memory `RawRefSeq` — *before* the
  `RecordSource` seam and the iterator wire it up. The decision is the thing that must match
  production; prove it first, plumb it second.
- **Reuse over rewrite.** Read filtering is a **port**: it calls the existing pure predicates
  (`read_exceeds_mismatch_fraction`, `cigar_is_bad`) and the `FLAG_*` / `DEFAULT_*` constants
  as-is, and reuses the `RecordBuf → MappedRead` decode path. It supplies only its own driver
  and config (spec §2.5). No filter logic is re-derived.
- **Verify against ground truth.** The north-star test is **drop-parity with the production
  filters** on a BAM/CRAM fixture (a hand-counted `ReadFilterCounts`), not self-consistency —
  a port is correct when it drops exactly what production drops.
- **Incremental, with pauses.** Land one milestone, stop for review, then the next.
- **Ungated / container builds.** `ng` stays a plain module; all `cargo` via `./scripts/dev.sh`
  (CLAUDE.md), a native host build at completion.

## Preconditions (already in place)

- **Foundations Milestones A + B are done** — so `RawRefSeq` (all filter #8 needs) is available:
  `InMemoryRefSeq` as the test oracle, `ResidentRefSeq` for real FASTA. `WindowedRefSeq`
  (foundations C) is **not** a dependency of read filtering. `types.rs` exists, seeded with
  `ContigId` + `RefSeqError` (foundations A1); `DomainError` is added here.
- **The production filter stack** in
  [`bam/alignment_input.rs`](../../../../src/bam/alignment_input.rs) —
  `read_exceeds_mismatch_fraction`, `cigar_is_bad`, `classify_pre_decode`, the `FLAG_*` and
  `DEFAULT_*` constants, `FilterCounts`, `compute_adaptor_boundary`, and `MappedRead` — is the
  reuse target **and** the parity oracle.

---

## The steps

### Milestone A — module scaffold + vocabulary + local types (types, no logic)

**A1. Scaffold the `read/` module.**  ✅
`src/ng/read/mod.rs` (declares `filtering`; the `ReadPrep` trait lands later, step 2) and an
empty `src/ng/read/filtering.rs` with its `#[cfg(test)]` block; wire `pub mod read;` into
`ng/mod.rs`. One file, no folder — no bake-off (module_layout 1a). *Source:* arch §Module home.

**A2. Extend `types.rs` with the scalar newtypes.**  ✅
`MapQual`, `BaseQual`, `Bp` (unconstrained — `pub` field, ergonomic derives, `.get()`) and
`MismatchFraction` (constrained — private field, `try_new` in `[0,1]`); add the
`MismatchFraction` variant to `DomainError`. Names match ng_step_interfaces §1. Unit test:
`MismatchFraction::try_new` accepts boundary values, rejects out-of-range. *Source:* spec
§2.2/§2.3, arch §2.1.

**A3. The step-1-local types.**  ✅
In `filtering.rs`: `ReadFilterConfig` (+ `Default` = the production policy, its thresholds the
reused `DEFAULT_*` constants), `FilterVerdict`, `DropReason` (variants 1:1 with the counts
fields), `ReadFilterCounts`. No logic. Test: `Default` reproduces the production defaults.
*Source:* spec §4, arch §2.2.

> **Checkpoint A:** types compile; `MismatchFraction` validation and the `Default` config are
> tested. Pause for review.

### Milestone B — the two-phase cascade (the decision, pure)

**B1. `verdict_pre_decode` (#1–#6).**  ✅
The flag/MAPQ bit-tests over the reused `FLAG_*` constants and `min_mapq`, in hit-rate order,
charging a read to the first filter it fails; reference-free and decode-free. Unit tests: each
filter's boundary (mapq exactly at `min_mapq` kept, one below dropped; each flag bit), the
`drop_qc_fail` / `drop_duplicate` toggles, and attribution order when several would fire.
*Source:* spec §3, §5.

**B2. `verdict_post_decode` (#7–#9).**  ✅
Too-short (`min_read_length` on decoded length), high-mismatch (call
`read_exceeds_mismatch_fraction` gated on `&impl RawRefSeq` + `mismatch_bq_floor`), bad-CIGAR
(call `cigar_is_bad`). Unit tests use `InMemoryRefSeq` as the reference; cover #8 on and off
(`None` ⇒ no reference access), a read exactly at the mismatch threshold, and the two bad-CIGAR
shapes. *Source:* spec §3, arch §3.

> **Checkpoint B:** the full cascade is decided and unit-tested in isolation — both halves,
> against the in-memory reference — with per-filter boundary coverage. Pause for review.

### Milestone C — the `RecordSource` / `RawRecord` seam

**C1. The traits + a test fake.**  ✅
`RawRecord` (`flag`, `mapq`, `decode(&self)`) and `RecordSource` (`read_next(&mut buf)`,
reused buffer) in `filtering.rs`. A trivial in-memory fake `RawRecord` and a fake
`RecordSource` yielding a fixed slice — enough to drive the cascade without a BAM. *Source:*
arch §3, spec §5. **Deviation (recorded):** `decode` returns `io::Result<MappedRead>`, not the
spec's illustrative infallible `MappedRead` — the reused `record_buf_to_mapped_read` is
fallible, and a decode failure is fatal like the #8 fetch (spec §7).

**C2. The ng-owned noodles adapters (BAM + CRAM).**  ✅
`NoodlesRawRecord` wraps a noodles `RecordBuf`: `flag`/`mapq` are cheap field reads,
`decode(&self)` runs the existing `RecordBuf → MappedRead` path (reusing
`compute_adaptor_boundary` via `record_buf_to_mapped_read`, whose visibility was lifted
`pub(super)`→`pub(crate)`). Two `RecordSource` impls: **`BamRecordSource<R>`** fills the reused
buffer in place via noodles `read_record_buf`; **`CramRecordSource<R>`** buffers one CRAM
*container*'s decoded records (CRAM decodes at container granularity into owned records and
consults a reference at decode time) and yields them one per read, refilling across containers,
passing its own `fasta::Repository` to the slice decoder. Both ng-owned — the dependency points
ng → existing code (spec §7, decision a). *(CRAM was initially deferred at Checkpoint C, then
added at the owner's request — a sibling `RecordSource`, no change to the seam.)*

> **Checkpoint C:** records flow through the seam; the fake source unit-tests it with no BAM;
> the BAM adapter reads flag/MAPQ pre-decode and decodes survivors (in-memory-BAM round-trip);
> the CRAM adapter matches noodles' own record iterator across multiple containers. Pause for
> review.

### Milestone D — the `ReadFilter` iterator (end-to-end)

**D1. `ReadFilter` + `new`.**  ☐
The struct (owns `source`, the reused `buf`, `reference`, `config`, `counts`) and
`new(...) -> Result<Self, RefSeqError>`: validates every source-header contig resolves in the
reference (fail-fast), seeds `buf = Default`. *Depends:* B, C. *Source:* spec §5, arch §3.

**D2. `Iterator::next` + `counts`.**  ☐
`read_next` into the reused `buf` → `verdict_pre_decode` → on survival `buf.decode()` →
`verdict_post_decode` → tally every drop → yield first `Keep`; an `Err` from `read_next` or #8
aborts the run (fatal, no per-item `Result`). `counts()` exposes the running tally. *Source:*
spec §5.

**D3. Fixture-driven integration test — the port anchor.**  ☐
Filter the reads of a small known **BAM** (and, since the CRAM source now exists, a CRAM too)
and assert the exact `ReadFilterCounts` against hand-counted expectations. Add a buffer-reuse
assertion (one `RecordBuf` across the pass) and a header/reference-mismatch case that fails in
`new`. **Note (from B):** hand-count against the **#7→#9→#8** order — per-bucket
`bad_cigar`/`high_mismatch_fraction` differ from production for both-failing reads by the
settled reorder, so *drop-parity vs the production filters is per-bucket-order-adjusted*, not
byte-identical. Also enforce the `DropReason`↔`ReadFilterCounts` 1:1 mapping (exhaustive
`match` at the tally site) — carried from A/B. *Source:* spec §2.6, arch §Test & bench shape.

> **Checkpoint D:** read filtering runs end-to-end on a real BAM; counts match the hand-count.
> **Step 1 is complete** (CRAM `RecordSource` is a tracked follow-up). Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | type-level tests — `MismatchFraction` range; `Default` config = production defaults |
| B | unit tests per filter boundary + cascade order/attribution, **both halves**, against `InMemoryRefSeq` |
| C | fake `RecordSource`/`RawRecord` unit tests (no BAM); the noodles adapter compiles and gates pre-decode |
| D | **fixture-driven integration:** exact `ReadFilterCounts` vs hand-count **and drop-parity vs the production filters**; buffer-reuse + fail-fast-`new` cases |

## Out of scope (next plans)

- **Step 2, read preparation (`ReadPrep`)** — the `read/` sibling; its own plan.
- **`pileup/`** (adaptor-clip application, mate-overlap reconciliation, the non-STR
  loci+evidence generator), the **locus router** (step 3), and the **research probes** —
  later (`ng_proposal.md` §3).
- **Output-`MappedRead` allocation reuse** — a locus-stream memory decision, measure-first
  (spec §6).
- **Pre-decode-gate magnitude** — the gate is result-preserving, so measuring its wall-clock
  win is a perf follow-up, never a correctness gate.
