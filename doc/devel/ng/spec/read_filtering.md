# ng step 1 — read filtering (and the ng foundations)

*Status: design spec (2026-07-10; renamed admission → filtering 2026-07-11). First
spec under [`ng_proposal.md`](ng_proposal.md); companion to the arch docs
[`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) (shared types +
step traits) and [`../arch/module_layout.md`](../arch/module_layout.md) (the
`src/ng/` tree). **No code yet** — this settles the design.*

*Naming follows the project convention: **STR** in prose, `ssr` in code
identifiers (`str` is a Rust primitive). Read filtering is generic — it has no
STR-specific surface — so `ssr` does not appear here.*

*The taxonomy step (`ng_proposal.md` §1) is named "Read admission & QC", the survey
label matching all five callers. Our module realises it as `read/filtering.rs`
because "filtering" names the mechanism plainly — and pairs with step 11's **site
filtering** (`locus_filter/`): filtering happens at two levels, reads and sites.*

**Reference dependency — now specced.** One filter (the mismatch-fraction filter, #8)
needs reference bases, as do the pileup, BAQ, and realignment. The reference-sequence
accessor is defined in its own spec, [`ref_seq.md`](ref_seq.md) (the `RefSeq` trait) —
written first because it is shared and #8 depends on it. Filter #8 takes a `&impl RawRefSeq`
and reads **raw** bytes (matching production, for a clean port-back); it stays in read
filtering (§3, §7).

---

## 1. What read filtering *is* — and what it is not

Read filtering is the **whole-read keep/drop prelude**. It takes the reads decoded
from one sample's alignment file (BAM/CRAM) and yields the subset worth carrying
forward, plus a tally of what was dropped and why. Three properties fix its scope:

- **Per-read.** Every decision is about one read in isolation: its flags, its
  mapping quality, its length, its per-base mismatches against the reference. No
  decision needs another read, a pair, or a locus.
- **Locus-independent.** Filtering runs *before* the caller knows which loci it will
  genotype (step 3, the router). A read is kept or dropped for the whole run, not per
  candidate site.
- **Content-preserving.** Filtering **selects** reads; it never rewrites them. A kept
  read's bytes, CIGAR, and qualities are exactly what was decoded. Every
  transformation of a read's *content* — realignment, per-base masking, quality
  reconciliation — belongs downstream.

That last property draws the one boundary line worth stating explicitly, because the
step map in `ng_proposal.md` §1 lists "overlapping-mate quality reconciliation" under
step 1, yet in the production code that work is **per-base and pairwise** and happens
far downstream, at pileup-walker time. See §6.

**Read filtering does *not*:**

- realign or re-map reads (that is **step 2, read preparation**);
- apply the adaptor-readthrough mask or reconcile overlapping mate-pair qualities
  (per-base operations that need the locus/pileup context — the **`pileup/`** module on
  the generic path, the pair-HMM prep on the STR path; see §6);
- assign reads to loci or decide STR-ness (**step 3, the router**);
- estimate any model parameter (**step 4, the pre-pass**).

Its output is simply *the reads that earned their place in the pipeline.*

---

## 2. The ng foundations settled here

### 2.1 Module structure

Read filtering and preparation are two stages of the same job — turning a decoded
alignment record into locus evidence — and they share the same input type
(`MappedRead`) and the same reference accessor.

```
src/ng/read/
├── mod.rs         – module declarations + the ReadPrep trait (step 2, added later)
├── filtering.rs   – step 1: the fixed filtering prelude  (this spec)  (+ #[cfg(test)])
└── …              – step 2 impls land here as siblings: trust_mapper.rs, pair_hmm.rs, …
```

This is a deliberate, documented deviation from the "one folder per step" rule in
`module_layout.md` principle 1, and it *keeps* that rule's real intent: a step's
competing implementations sit side by side. Read filtering has no competitors (it is a
fixed prelude), so it is a single file; step 2's `ReadPrep` implementations still sit
side by side within `read/`.

### 2.2 The shared vocabulary — seeding `types.rs`

Per `module_layout.md` principle 3, the cross-step domain newtypes begin in one
`src/ng/types.rs` and split into concept modules (`units`, `locus`, …) as they grow.
This spec seeds `types.rs` with **only the newtypes read filtering actually touches** —
not the whole vocabulary from `ng_step_interfaces.md` §1.

Read filtering touches three unconstrained scalar newtypes and one constrained one:

```rust
// Unconstrained — any value of the primitive is a valid instance, so the field
// is `pub` and there is no checked constructor (a constructor would be ceremony).
pub struct MapQual(pub u8);   // SAM mapping quality (MAPQ): the aligner's Phred-scaled
                              //   confidence that the read is placed at the right locus.
                              //   0 = "could be anywhere"; 60 = "as sure as this aligner gets".
pub struct BaseQual(pub u8);  // a single base call's Phred quality (0–93).
pub struct Bp(pub u32);       // a length in base pairs — here, a read's decoded length.

// Constrained — has a domain range, so the field is PRIVATE and construction goes
// through a checked constructor. This is the one place read filtering exercises the
// validated-newtype half of the convention (§2.3).
pub struct MismatchFraction(f32);   // in [0, 1]
```

`MapQual`, `BaseQual`, `Bp` match the target names in `ng_step_interfaces.md` §1
exactly, so no later renaming is needed when the rest of the vocabulary lands. (`Bp`
stays **unprefixed**: base pairs are the generic length currency both the SNP/indel
and STR paths speak — only *repeat-unit* quantities carry the `Ssr` prefix.)

### 2.3 The newtype and validation conventions

Two rules, both demonstrated above, that every later ng type will follow:

- **Unconstrained newtypes keep `pub` fields.** `MapQual`, `BaseQual`, `Bp` — there is
  no invariant to protect (every `u8`/`u32` is a legal value), so a checked
  constructor would be ceremony. They derive the ergonomic set
  (`Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug`) and expose `.get()`.

- **Constrained newtypes hide their field and construct through a checked
  constructor**, so an illegal value is *unrepresentable* rather than merely
  discouraged. `MismatchFraction` is our worked example. Its source is an untrusted
  CLI/config value, so the policy from `ng_step_interfaces.md` §1 applies: **fail
  loudly, never silently coerce.** (This is the *threshold* value — validated
  regardless of how the fraction itself ends up being computed, §3.)

  ```rust
  impl MismatchFraction {
      /// The only constructor. A fraction outside [0, 1] is a user error — reject it.
      pub fn try_new(x: f32) -> Result<Self, DomainError> {
          (0.0..=1.0).contains(&x).then_some(Self(x)).ok_or(DomainError::MismatchFraction(x))
      }
      pub fn get(self) -> f32 { self.0 }
  }
  ```

  `DomainError` is introduced here (the ng-wide error for a domain-invariant
  violation) with its first variant; later constrained types (`AlleleFreq`,
  `InbreedingF`, `Theta`, …) add their own variants as they arrive.

### 2.4 The config convention

`ReadFilterConfig` (§4) is the template every later step's config follows. Three
conventions it fixes:

- **`Option<T>` is the "no threshold" state, not a sentinel.** `min_mapq: Option<MapQual>`
  — `None` means *no minimum*, `Some(q)` means *drop reads below `q`*. Using `0` as a
  disable value would make `Some(0)` and `None` two spellings of the same behaviour;
  `Option` makes "absent" structurally distinct from "zero". (This mirrors the existing
  `AlignmentMergedReaderConfig`.)
- **Every default is a named `pub const` with a doc comment** giving its units and its
  source (a spec, a measurement, or a vendor/tool convention) — no magic numbers. Read
  filtering reuses the existing constants where the value is genuinely the same (§3).
- **`Default` is implemented** and is the production default the lab runs with. The
  config is **minimal**: it holds exactly the knobs of the filters that are actually
  active (§3). Dormant levers (downsampling, read pooling) are *not* added as disabled
  fields — they enter the config only when they enter the pipeline (§6).

### 2.5 Reuse over rewrite — the map to `alignment_input`

Read filtering is a **port**, and the port reuses the existing *pure predicates*
rather than re-deriving them. The existing filter code separates cleanly into "the
decision" (pure functions) and "where it runs" (the CRAM/segment reader plumbing); ng
reuses the former and supplies its own driver.

| what | existing code | ng reuse |
|---|---|---|
| flag bit constants | `FLAG_DUPLICATE`, `FLAG_SECONDARY`, `FLAG_SUPPLEMENTARY`, `FLAG_UNMAPPED`, `FLAG_QC_FAIL`, … | import as-is |
| default thresholds | `DEFAULT_MIN_MAPQ` (20), `DEFAULT_MIN_READ_LENGTH` (30), `DEFAULT_MAX_READ_MISMATCH_FRACTION` (0.10), `DEFAULT_MISMATCH_BQ_FLOOR` (10) | import as-is |
| mismatch-fraction test | `read_exceeds_mismatch_fraction(...)` | call directly — **but gated on the reference accessor** (§3, §7) |
| bad-CIGAR test | `cigar_is_bad(...)` | call directly (pure — no reference) |
| adaptor boundary (annotation) | `compute_adaptor_boundary(...)` → stored on `MappedRead.adaptor_boundary` | carried through, not applied (§6) |
| the read itself | `MappedRead` | **reuse as-is** — it is the step-2 input (`ng_step_interfaces.md` §6) |

The one predicate ng does **not** call directly is `classify_pre_decode`: it takes a
noodles `RecordBuf` (a not-yet-decoded record), because in production the flag and MAPQ
tests run *before* decoding to save work. ng operates on already-decoded `MappedRead`s
(the lab decodes everything; it is not chasing that optimisation), so it re-expresses
the same six flag/MAPQ tests over `MappedRead.flag` and `MappedRead.mapq` — a handful
of bit-tests reusing the same `FLAG_*` constants and `DEFAULT_MIN_MAPQ`. Same logic,
same constants, different input type.

**Where the filter config lives — a deliberate choice.** In production the flag/MAPQ
gate is applied by the reader (via its `AlignmentMergedReaderConfig`), so a read may
arrive already filtered by a policy set elsewhere. ng instead makes read filtering own
its **entire** policy in one visible `ReadFilterConfig` and applies the full cascade
itself: the reader is asked only to *decode* reads, not to filter them. Filtering being
a real, inspectable step — the whole point of the decomposition — means its policy must
be in one place, not inherited from a hidden reader default.

> **Boy-scout note.** Two decisions are left open to *improve if it helps readability*,
> not frozen: (a) whether kept reads deserve a distinct type (`FilteredRead`) rather
> than a bare `MappedRead` (§4), and (b) exactly how much to reuse the old predicates
> versus tidy them as we lift them into ng. The spec records the leaning; the
> implementation is free to leave the code a little better than it found it where a
> clear win appears, and to write that down when it does.

### 2.6 The test and bench shape

- **Unit tests beside the code** — the standard `#[cfg(test)] mod tests` block in
  `filtering.rs`, exercising each filter's boundary (e.g. a read exactly at `min_mapq`
  is kept; one below is dropped) and the cascade order. The existing `alignment_input`
  tests are the template and, for the ported predicates, the behaviour to match.
- **One fixture-driven integration test** — filter the reads of a small known BAM/CRAM
  and assert the resulting `ReadFilterStats` (§4) matches hand-counted expectations.
  This is the regression anchor for the port.
- **`bench/` is deferred for this step, on purpose.** The lab's `bench/` module scores
  *competing* implementations against gold/silver/synthetic standards
  (`module_layout.md` principle 4). Read filtering has no competitors, so it has no
  frontier to plot; its "measurement" is the regression test above. The `bench/` seam
  is described here only so later steps inherit the expectation — the first real use of
  `bench/` is the candidate-generation probe (`ng_proposal.md` §3).

---

## 3. The filters — the port

The filtering cascade, in application order. The order is **hit-rate-ordered** (the
cheapest, most-often-firing tests first) so the average read is dropped as early as
possible — ported verbatim from the production rationale in
`per_sample_caller_cram_input.md`. Order matters only for the drop *attribution* (a
read is charged to the first filter it fails), not for the keep/drop outcome, which is
order-independent.

| # | filter | drops a read when… | default | source const | toggle |
|---|---|---|---|---|---|
| 1 | **duplicate** | `flag & FLAG_DUPLICATE` — the read is a PCR/optical duplicate (another copy of the same original molecule; keeping both would double-count the evidence) | on | — | `drop_duplicate` |
| 2 | **low MAPQ** | `mapq < min_mapq` — the aligner is not confident the read is placed at the right locus. MAPQ unavailable (SAM `0xFF`) is treated as 0, so any non-zero minimum drops it | `Some(20)` | `DEFAULT_MIN_MAPQ` | `min_mapq` |
| 3 | **supplementary** | `flag & FLAG_SUPPLEMENTARY` — a chunk of a chimeric read whose other chunks are tracked separately | on (unconditional) | — | none |
| 4 | **secondary** | `flag & FLAG_SECONDARY` — a duplicate projection of a read already represented by its primary alignment | on (unconditional) | — | none |
| 5 | **unmapped** | `flag & FLAG_UNMAPPED` — no alignment position, so no allele evidence | on (unconditional) | — | none |
| 6 | **QC fail** | `flag & FLAG_QC_FAIL` — the sequencer/pipeline flagged the read as failing quality control | on | — | `drop_qc_fail` |
| 7 | **too short** | decoded read length `< min_read_length` — very short reads rarely contribute a reliable alignment | `Some(30)` | `DEFAULT_MIN_READ_LENGTH` | `min_read_length` |
| 8 | **high mismatch fraction** ⚠ | the read's fraction of quality-clearing mismatches against the reference exceeds the threshold — a defence against contamination, adaptor runthrough, and chimeric tails the aligner placed anyway. **Reference-dependent — pending the reference-accessor spec** | `Some(0.10)` | `DEFAULT_MAX_READ_MISMATCH_FRACTION` | `max_read_mismatch_fraction` |
| 9 | **bad CIGAR** | the alignment is ill-formed: an adjacent insertion/deletion pair, or a deletion at a read boundary — signals the aligner was genuinely confused | on | — | none (always applied) |

**Supplementary / secondary / unmapped have no toggle** — this matches the production
policy and its rationale (`AlignmentMergedReaderConfig`'s "why no
`drop_secondary`/`drop_supplementary` toggles" note): admitting either would silently
double-count alleles or break downstream per-read invariants. There is no correct
configuration that lets them through.

**Filter #8 is the reference-dependent one — now settled via `RefSeq`.** Counting
mismatches needs to compare each aligned read base to the reference base under it. With a
plain `M` CIGAR op (the common bwa/bowtie case), match-vs-mismatch is *not* encoded, so
this needs reference bases. Since the pileup, BAQ, and realignment all need real
reference bases regardless, the shared accessor is defined in [`ref_seq.md`](ref_seq.md).
Filter #8 therefore:

- **stays in read filtering**, taking a `&impl RawRefSeq` threaded into `ReadFilter::new`;
- reads **raw** (un-canonicalised) bytes via `RawRefSeq::fetch_raw`, matching production's
  `RawContigRefCache` → `read_exceeds_mismatch_fraction` path exactly, so the ported
  filter behaves identically (the point of matching the old byte convention);
- turns off cleanly when `max_read_mismatch_fraction` is `None`, in which case no
  reference access happens at all.

(The `MD` tag would also let this run reference-free, but `MappedRead` doesn't carry it
and the accessor is needed anyway — so raw `RawRefSeq` bytes it is.) Filters #1–#7 and #9
have no reference dependency and are settled independently.

**Jargon, once:** a **mismatch fraction** is `(mismatches / comparable bases)` over the
read's aligned (`M`-op) positions, counting only positions where both the read base and
the reference base are a real nucleotide (`A/C/G/T`) and the base quality clears a floor
(`mismatch_bq_floor`, default 10) — so genuinely low-quality bases, which the likelihood
model already de-emphasises, do not trip a whole-read drop.

---

## 4. The types read filtering introduces

Three small types, all in `read/filtering.rs` (they are step-1-local, not shared
vocabulary — only the *scalar newtypes* of §2.2 go in `types.rs`).

```rust
/// The filtering policy: which filters are active and their thresholds. Minimal by
/// design — one field per active filter, no dormant levers. `Default` is the
/// production policy the lab runs with. Mirrors the filtering subset of the existing
/// `AlignmentMergedReaderConfig`.
pub struct ReadFilterConfig {
    pub min_mapq: Option<MapQual>,               // None = no minimum
    pub min_read_length: Option<Bp>,             // None = no minimum
    pub drop_qc_fail: bool,
    pub drop_duplicate: bool,
    pub max_read_mismatch_fraction: Option<MismatchFraction>,  // None = filter off; see §3 (#8)
    pub mismatch_bq_floor: BaseQual,             // only meaningful when the above is Some
}

/// The verdict for one read. `Keep` carries the read on; `Drop` records which filter
/// fired (the first one, per the hit-rate order) for the tally.
pub enum FilterVerdict {
    Keep,
    Drop(DropReason),   // enum: Duplicate | LowMapq | Supplementary | Secondary |
                        //       Unmapped | QcFail | TooShort | HighMismatch | BadCigar
}

/// A per-sample tally of the filtering pass — one counter per drop reason, plus the
/// kept count. The ng port of `FilterCounts`. Surfacing every drop is the "no silent
/// caps" discipline: a read that vanished must be accounted for. It is a **running**
/// tally (§5): readable at any point, final once the input is exhausted.
pub struct ReadFilterStats {
    pub kept: u64,
    pub duplicate: u64,
    pub low_mapq: u64,
    pub supplementary: u64,
    pub secondary: u64,
    pub unmapped: u64,
    pub qc_fail: u64,
    pub too_short: u64,
    pub high_mismatch_fraction: u64,
    pub bad_cigar: u64,
}
```

**No new *read* type.** A kept read is a `MappedRead`, unchanged — filtering selects, it
does not transform (§1). `MappedRead` is already the documented input to step 2
(`prepare_read(&self, read: &MappedRead, …)` in `ng_step_interfaces.md` §3), so making
filtering yield `MappedRead` keeps the two steps composable with no adapter. *(Boy-scout
open item: if a typestate `FilteredRead(MappedRead)` wrapper buys enough "you cannot
skip filtering" safety to be worth the ceremony, add it when that value is demonstrated
— see §2.5. The leaning is to keep the bare `MappedRead`.)*

---

## 5. The public surface — a filtering iterator that carries its stats

Read filtering is exposed as an **iterator adapter**: construct a `ReadFilter` from the
input read iterator and the config, then iterate it to pull the reads that pass. It
holds the running `ReadFilterStats` as reads flow through. Lazy and streaming — no
upfront `Vec` — so it drops straight into the locus-stream spine (`ng_proposal.md` §1)
as its first stage.

```rust
/// Filters one sample's decoded reads, lazily. `next()` pulls from the input, applies
/// the §3 cascade to each read, tallies every drop, and returns the first read that
/// passes (or `None` when the input is exhausted).
///
/// `stats()` is a RUNNING tally: readable at any point — peek mid-stream to watch how
/// filtering is going — and complete once the iterator is exhausted. Reaching the final
/// totals is the caller's responsibility (iterate to the end); nothing forces full
/// consumption.
pub struct ReadFilter<I, R /*: I: Iterator<Item = MappedRead>, R: RawRefSeq */> {
    reads: I,
    reference: R,                 // RawRefSeq — raw bytes for filter #8; see §3 + ref_seq.md
    config: ReadFilterConfig,
    stats: ReadFilterStats,
}

impl<I: Iterator<Item = MappedRead>, R: RawRefSeq> Iterator for ReadFilter<I, R> {
    type Item = MappedRead;
    fn next(&mut self) -> Option<MappedRead> { /* pull → cascade → tally → first Keep */ }
}

impl<I, R> ReadFilter<I, R> {
    pub fn new(reads: I, reference: R, config: ReadFilterConfig) -> Self { /* … */ }
    /// The running tally — current counts, final once iteration is exhausted.
    pub fn stats(&self) -> &ReadFilterStats { &self.stats }
}
```

Usage — iterate by `&mut` so the filter (and its stats) outlives the loop:

```rust
let mut filter = ReadFilter::new(decoded_reads, reference, config);
for read in &mut filter {
    // … feed the locus stream
}
let stats = filter.stats();   // running totals; complete here because the loop drained it
```

The single-read decision is factored into a helper so the cascade is unit-testable in
isolation (with the in-memory `RawRefSeq` impl standing in for a FASTA, `ref_seq.md` §2):

```rust
fn verdict(read: &MappedRead, reference: &impl RawRefSeq, config: &ReadFilterConfig) -> FilterVerdict;
```

*Contract:* content-preserving, per-sample, lazy; every dropped read is charged to
exactly one reason in the running stats. `reference` is consulted only by filter #8, and
only when `max_read_mismatch_fraction` is `Some` (§3). (The `R: RawRefSeq` bound holds even
when #8 is disabled at runtime; a caller wanting *no* reference at all would be a future
`Option<R>` refinement — not needed now, since every real run has a reference anyway.)

---

## 6. Deferred, with a recommended home

Everything listed under "read admission" in `ng_proposal.md` §1 that this step does
**not** implement, and where it should live instead. Recording the home is the point —
nothing is silently dropped.

Adaptor masking and mate-overlap reconciliation are per-base operations that need the
evidence-gathering context, so on the **generic path they live in the `pileup/`
module** (the reused genome-walking pileup — `ng_proposal.md` §1, *The locus stream*),
exactly where production does them today; on the **STR path** they fold into the
pair-HMM read preparation. Neither is a per-read filtering decision.

- **Adaptor-clip *application* → `pileup/` (generic) / pair-HMM prep (STR).** The
  adaptor boundary (the reference position past which a read has sequenced into its
  mate-pair adaptor) is a cheap *per-read annotation*, already computed at decode by
  `compute_adaptor_boundary` and carried on `MappedRead.adaptor_boundary`. Filtering
  preserves it. But *applying* it — masking the bases past the boundary — is a per-base
  operation done at pileup-walker time in production
  ([`decompose.rs`](../../../../src/pileup/walker/decompose.rs)), so in ng it lands in
  `pileup/`, where a read's content is turned into locus evidence.

- **Mate-overlap quality reconciliation → `pileup/` (generic) / pair-HMM prep (STR).**
  Where the two reads of a pair overlap on the reference, their bases are not
  independent evidence (they come from one molecule), so agreeing bases have their
  qualities combined and disagreeing ones are down-weighted — the samtools-style
  reconciliation. This is inherently **pairwise and per-base**, needs the pileup-column
  context, and cannot be a per-read filtering decision. It is catalogued under step 1 in
  the step map, but the code has always done it downstream, and ng follows the code.

- **Coverage downsampling and read pooling → future config, off by default.** Neither
  exists in our production caller today ("No pooling/downsampling" in the step
  catalogue), and neither is needed yet. They are genuine future levers (downsampling
  caps reads per start position; pooling runs the likelihood once per unique read
  sequence). They enter `ReadFilterConfig` as real knobs only when they enter the
  pipeline — not as dormant disabled fields now (§2.4).

This is consistent with `ng_proposal.md` §1: the step map catalogues overlap
reconciliation under step 1 (as all five surveyed callers do), but the *locus-stream*
subsection there places the per-base evidence work in `pileup/` — read filtering is the
whole-read prelude, not the per-base gatherer.

---

## 7. Open questions

- **The reference-sequence accessor — resolved: `RefSeq`** ([`ref_seq.md`](ref_seq.md)),
  a dedicated spec written first because it is shared (pileup, BAQ, realignment, DUST) and
  #8 depends on it. Filter #8 stays in read filtering, takes a `&impl RawRefSeq`, and reads
  raw bytes (§3). (Coordinate base is settled: ng is uniformly **1-based**, matching
  production and VCF/SAM.)
- **`FilteredRead` typestate?** Keep the bare `MappedRead` (current leaning), or wrap it
  so "reached step 2 without filtering" is a compile error? Decide if/when the safety
  proves worth the ceremony (§4, §2.5).
- **Does the lab ever want the flag/MAPQ gate at decode time?** ng decodes everything
  and filters after, for simplicity. If a future large-cohort ng experiment makes decode
  cost bite, the pre-decode `classify_pre_decode` optimisation can be lifted in — but not
  before it measurably matters.
```
