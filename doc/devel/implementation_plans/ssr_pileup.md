# Stage 1 — `ssr-pileup` implementation sketch

**Status:** first draft, 2026-06-15. Files + main structs + central function
signatures — the *shape* we'll implement, not full pseudocode. Built on the
architecture ([ssr_pileup.md](../architecture/ssr_pileup.md), every structural
question decided in its §14), the shared types
([ssr_shared_types.md](../architecture/ssr_shared_types.md)), Stage 0
([ssr_catalog.md](ssr_catalog.md) — `Locus`/`Motif`/`CatalogReader`), and the spec
(§4). We iterate on this until we agree, then code it.

Grounding facts already checked in the tree:
- [`src/ssr/types.rs`](../../src/ssr/types.rs) today carries **only** `Motif` +
  `Locus`. `Locus` already exposes everything the SSR math reads off the catalog:
  `left_flank()` / `right_flank()` / `ref_tract()` / `motif()` / `ref_bytes()`,
  all from embedded bytes — **no FASTA opened by the algorithm** (arch §1). The
  `Allele` / `NormalizedSeq` types are *designed* (shared-types §2–§5) but **not
  built**; building them is task 1 (arch §13).
- [`AlignmentMergedReader`](../../src/bam/alignment_input.rs) already does the
  k-way coordinate merge, the CRAM `Repository` sharing, the `classify_pre_decode`
  MAPQ/flag/length admission gate, and a region `query()`. Stage 1 **reuses it
  verbatim** (arch §3.1/§8.1) — same reader, same `--min-mapq`.
- `normalize_alleles` lives **private** in
  [pileup/walker/indel_norm.rs](../../src/pileup/walker/indel_norm.rs); its body
  already operates on an abstract `(seqs, bounds)` pair. The lift to a shared
  module is task 2 (arch §7), behaviour-preserving + green before the SSR adapter.
- The Dindel/BAQ per-Q emission pattern (a process-wide lookup) is the precedent
  for our emission table — [baq/scratch.rs](../../src/baq/scratch.rs) `Q2P`,
  [baq/probaln.rs](../../src/baq/probaln.rs). We borrow the **pattern + the banded
  loop shape + the grow-and-keep scratch discipline**, not the code (arch §5.6).
- The generic columnar container + `.ssr.psp` SSR schema (`registry_ssr` /
  `SsrLocusRecord`, arch §10.4) is a writer prerequisite — task 3, independent of
  the SSR math.

---

## 0. Build order (arch §13 — strict prerequisites first)

```
1. types.rs    : Allele / NormalizedSeq + to_sequence / repeat_count   (§1 below)
2. shared lift : normalize_alleles → shared module, SNP tests green     (§2 below)
3. container   : registry_ssr + SsrLocusRecord schema (if not landed)   (§8 below)
4. stage       : candidate_generation → pair_hmm → triage(coverage+centre) → fetch_reads → mod  (§3–§9)
                 (count_repeats parked as the deferred, measured fast-path shortcut — §2/§4)
```
Each is landed and green before the next. The stage modules themselves go
**bottom-up** (pure scorers + candidate builders first, the I/O fetcher and the
driver last) so every layer can be unit-tested before it has a consumer.

---

## 1. Task 1 — extend `types.rs` with the allele representation (shared-types §2–§5)

The first consumer of the designed-but-unbuilt allele model. No new files; this
grows `src/ssr/types.rs`.

```rust
/// Canonical normalized off-ladder sequence — identity is its bytes (rule 3).
/// Left-aligned + minimally trimmed by `normalize_alleles` before construction,
/// so derived Eq/Hash give cross-sample identity for free.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct NormalizedSeq(Box<[u8]>);

/// A locus-relative allele (shared-types §3: does NOT carry its locus).
/// Two encodings of one tract sequence.
#[derive(Clone, PartialEq, Eq, Hash)]
pub enum Allele {
    OnLadder { units: u16 },          // clean rung; sequence = flank+motif×units+flank
    OffLadder(NormalizedSeq),         // non-rung; literal canonical bytes
}

impl Allele {
    /// Reconstruct the full tract+flank bytes against this locus (pure fn of
    /// (allele, locus); on-ladder tiles the motif, off-ladder returns its bytes).
    pub fn to_sequence(&self, locus: &Locus) -> Vec<u8>;
    /// Repeat count (integer for OnLadder; fractional ref-relative for OffLadder).
    pub fn repeat_count(&self, locus: &Locus) -> f64;
}
```
Round-trip invariant tested here (independent of the rest of the stage):
`to_sequence` of an `OnLadder { units }` re-parses back to the same `units`;
`NormalizedSeq` built two ways from the same allele is byte-identical.

---

## 2. Task 2 — lift `normalize_alleles` to a shared module (arch §7)

A real refactor with a regression gate, **not** a rename:

```rust
// new: src/norm_seqs/mod.rs  (representation-neutral, public — it normalizes
// sequences; neither "reads" nor "alleles" leaks into the shared name)
pub fn normalize_alleles(seqs: &mut [Vec<u8>], bounds: Range<usize>) -> Range<usize>;
```
- Move the existing private fn + its `Range`/helpers out of
  [pileup/walker/indel_norm.rs](../../src/pileup/walker/indel_norm.rs); the SNP
  CIGAR path keeps its `left_align_cigar` wrapper and now *calls* the shared fn.
- **Regression gate = SNP end-to-end tests pass unchanged** (it is GATK/freebayes
  cross-validated; must not drift).
- The SSR side gets a thin adapter (lives in `candidate_generation.rs`, §5) that builds
  `(seqs, bounds)` from `(off-ladder candidate tract, ref tract)` — not a
  CIGAR-faking shim.

*(The kernel keeps its `normalize_alleles` name inside `norm_seqs/`; only the
module path is new. Contract: "one kernel, two users".)*

---

## 3. File layout for the stage proper (arch §13)

```
src/ssr/pileup/
├── mod.rs           # SsrPileupArgs, run(): wire fetcher→queue→pool→collector→writer (§9)
├── fetch_reads.rs   # single I/O thread: forward index walk, depth-cap reservoir, bundles (§7)
├── triage.rs        # coverage classification + region extract + window centre (pre-probe) (§4)
├── count_repeats.rs # parked: confident pure-tiling count — the deferred fast-path shortcut (§4)
├── pair_hmm.rs      # NEW banded 3-state log-space forward + PairHmmScratch (pure scorer) (§6)
├── candidate_generation.rs  # candidate gen: rungs + read-derived off-ladder + normalize (§5)
└── locus_record.rs  # per-locus aggregation → SsrLocusRecord (hist_*/amb_*/offl_*) (§8)
```
Subcommand: an `SsrPileupArgs` (clap) dispatched from the existing router beside
`pileup` / `var-calling` / `ssr-catalog`, reusing the shared
`--reference` / `--threads` / `--regions` / `--block-window-bp` parsers.
`--reference` is **CRAM-decode-only** here (spec §3.2).

---

## 4. `triage.rs` — coverage classification, then extract + centre (arch §2/§3)

> **Revised (2026-06-15, realign-everything — §2 banner).** Triage no longer runs
> a fast/slow gate or trusts the CIGAR's indel placement. It classifies *coverage*
> (CIGAR-lite: mapping position + ref-consumed length + clip lengths) and, for a
> spanning read, extracts the locus-region bases and centres the pair-HMM window;
> **every spanning read is then realigned** (§5/§6). The `Algorithm` /
> `SlowReason` / `FastHit` gate types are the **deferred, measured optimization**,
> not v1.

Takes the (gated, capped) reads of one bundle and decides, per read: does its
footprint span the locus, and if so what bases to score + where to centre the
window? The subtle stage — everything downstream trusts the coverage call.

```rust
/// Outcome of triaging one read against one locus (arch §3.3) — how many of the
/// tract's two ends the read's footprint brackets (2/1/0). Only `Spanning` feeds
/// the likelihood; the other two are tallied for QC, then dropped.
enum TriageResult {
    Spanning(SpanningRead),  // footprint brackets the tract + MIN_FLANK_BP, both ends
    Flanking,                // brackets one end only → n_flanking, not used (v1)
    InRepeat,                // brackets neither → n_frr, not used (v1)
}

/// A spanning read ready for the pair-HMM: the read-coord span covering the
/// locus region (flanks + tract, soft-clips included) and the observed repeat
/// count that centres the `count ± W` window. No fast/slow field — v1 realigns
/// every spanning read (§2 revision); the direct-count shortcut adds one back later.
struct SpanningRead {
    region: Range<usize>,    // read-coord span covering [locus ± FLANK_BP], clips included
    observed_count: u16,     // pre-probe longest-run → window centre
    // base-quality slice handle for the emission model
}

/// Classify by footprint coverage; for a spanning read, extract the region and
/// centre via the content pre-probe. Uses the CIGAR only for the footprint
/// (start + ref-consumed length + clip lengths), never its indel placement.
fn triage_read(read: &MappedRead, locus: &Locus, p: &SsrParams) -> TriageResult;
```

`find_longest_stretch` (the content pre-probe) is built in this module — longest
contiguous motif run (window centre) + total copies. The coverage test: the
read's footprint (aligned ref-span plus any soft-clip on the missing side)
reaches ≥ `MIN_FLANK_BP` past the tract on a side; both → `Spanning`, one →
`Flanking`, neither → `InRepeat`. **No flank-byte match at triage** — the
pair-HMM's flank emission (§5.3) judges flank quality and scores adapter/junk
clips poorly on its own (triage dumb, HMM smart). Per-locus QC counts (`depth`,
`n_spanning`, `n_flanking`, `n_frr`, `n_filtered`, `mapped_reads`) are tallied in
the fetcher's full pass (§7).

### `count_repeats.rs` — the deferred fast-path shortcut (arch §4, §2 revision)

Built and parked. `count_pure_tiling(tract, quals, motif) -> Option<(u16, f32)>`
returns a confident `(units, weight)` when a read's tract is a clean tiling — the
**measured optimization** to skip the pair-HMM for the easy majority, added back
once the fast-path fraction is measured (arch §14). Its `pure_tiling_units` is
already in use by the off-ladder degenerate check (§5).

---

## 5. `candidate_generation.rs` — candidate generation + off-ladder normalization (arch §6)

The translation layer between raw read bytes/counts and typed `Allele`s.
**Generates** the candidate haplotype sequences the pair-HMM scores; the DP never
constructs a sequence (arch §5.7).

```rust
/// One candidate allele the read is scored against: its full local DNA sequence
/// (what the forward aligns the read to) + which Allele that sequence encodes
/// (the key the Qᵣ result is recorded under). `candidate_seq` = upper-cased bytes,
/// left_flank + middle + right_flank, where middle = motif×L (rung) or the
/// normalized tract (off-ladder).
struct CandidateAllele { candidate_seq: Vec<u8>, allele: Allele }

/// Job 1 — build the observed_count±W rungs (left_flank + motif×L + right_flank)
/// for the read's window, from the catalog locus (no FASTA). Byte-identical across samples.
fn build_rungs(locus: &Locus, observed_count: u16, w: u16, out: &mut Vec<CandidateAllele>);

/// Job 1b — at most one read-derived off-ladder candidate, GATED on SlowReason
/// (only Impure / InteriorIndel; soft-clip-but-pure + low-Q stay on-ladder, §5.8):
///   left_flank + normalize(observed_tract) + right_flank.
/// Dropped if normalize(observed_tract) is itself a pure tiling (it's already a rung).
fn build_offladder(locus: &Locus, observed_tract: &[u8], reason: SlowReason,
                   out: &mut Vec<CandidateAllele>);

/// Job 2 — canonicalize an off-ladder tract to ONE spelling via the shared kernel
/// (§2 lift), so the same allele in any two samples → identical bytes (cohort union).
fn normalize_offladder(observed_tract: &[u8], ref_tract: &[u8]) -> NormalizedSeq;
```
Scratch-buffer discipline: `out` is a reused per-worker buffer (cleared per read),
haplotype bytes built into a reused scratch `Vec`.

---

## 6. `pair_hmm.rs` — the banded 3-state forward (arch §5, the net-new risk)

HipSTR's flank pair-HMM **minus its stutter marginalization**: a forward
(sum-over-alignments, log-space), pure scorer `read × haplotype → Qᵣ`, **no
traceback** ⇒ rolling two-row scratch. Reuses the *pattern*, not the BAQ code.

```rust
/// Per-worker scratch — grow-and-keep, zero per-read allocation on the hot path
/// (BAQ scratch.rs ethos). Two DP rows (prev/cur) × 3 states, sized to band.
struct PairHmmScratch { prev: Vec<f64>, cur: Vec<f64> /* M/I/D interleaved */ }
impl PairHmmScratch {
    fn resize_for(&mut self, read_len: usize, band: usize);
}

/// Forward log-likelihood of `read` (with per-base quals) against one candidate's
/// `candidate_seq`, banded to PAIR_HMM_BAND_BP around the diagonal. log-sum-exp
/// with a rescale guard (BAQ RESCALE_THRESHOLD pattern). Returns Qᵣ for this candidate.
fn forward(read: &[u8], quals: &[u8], candidate_seq: &[u8],
           scratch: &mut PairHmmScratch, model: &HmmModel) -> f64;

/// Process-wide emission table (Dindel/BAQ per-Q model, a LazyLock 256-entry
/// lookup): match = ln(1−10^(−Q/10)), mismatch = ln(10^(−Q/10)/3).
struct HmmModel {
    emit: &'static [Emit; 256],         // per-Q (match_ln, mismatch_ln)
    gap_open: &'static [f64],           // DINDEL_GAP_OPEN[], homopolymer-run-indexed (1..=10, extrap, cap 15)
    gap_extend: f64,                    // e^-1; gap→M = 1 − e^-1
}

/// Score one read against its whole candidate set → a dense Qᵣ over candidates,
/// independent narrow-band forwards (arch §5.5 — the recommended first layout;
/// shared-flank lattice is a measured optimization, deferred).
fn score_candidates(read: &SpanningRead, cands: &[CandidateAllele],
                    scratch: &mut PairHmmScratch, model: &HmmModel) -> Vec<(Allele, f64)>;
```
On- and off-ladder candidates go through the *same* `forward` (arch §5.7); base
quality (the emission model) arbitrates off-ladder vs nearest rung (arch §5.8).

---

## 7. `fetch_reads.rs` — single I/O thread: index walk + reservoir + bundles (arch §8)

```rust
/// A self-contained unit of work handed to a worker: one locus + its (capped) reads
/// + the per-locus QC tallies (computed over the WHOLE pass, not just kept reads).
struct LocusBundle { locus: Locus, reads: Vec<MergedRecord>, qc: LocusQc }

struct LocusQc { depth, n_spanning, n_flanking, n_frr, n_filtered,
                 n_flank_indel, mapped_reads: u32 }   // all per-locus

/// Forward-only walk of the sorted catalog. For each locus (or close cluster, the
/// one free fetch optimization): index-jump to [start−FLANK_BP, end+FLANK_BP],
/// read forward to the first record past the window, apply the cheap coordinate-
/// reach gate, reservoir-sample admitted reads to MAX_READS_PER_LOCUS, tally QC,
/// push a LocusBundle onto the bounded queue (blocks when full = back-pressure).
fn run_fetcher(reader: AlignmentMergedReader, catalog: CatalogReader,
               p: &SsrParams, tx: BoundedSender<LocusBundle>) -> Result<()>;

/// Algorithm R reservoir, seeded DETERMINISTICALLY from (chrom, start) — NOT
/// wall-clock/thread-id (arch §8.3). Keeps first K admitted reads; i-th admitted
/// (i>K) kept w.p. K/i, evicting one uniformly. Reads arrive in the merged
/// reader's total order (ref_id,pos → source idx → record order), so the i-th
/// read is well-defined → byte-identical across --threads.
struct Reservoir { k: usize, rng: DeterministicRng, held: Vec<MergedRecord>, seen: u64 }
```

- **Cheap coordinate-reach gate** (in the fetcher, no sequence scan): skip
  reservoir-admission for reads whose aligned footprint clearly can't span and are
  not soft-clipped on the missing side. Soft-clips are **always admitted**
  (conservative — only §4's content scan can tell).
- **Bounded queue keyed by reads/bytes in flight**, not bundle count (a deep
  locus is a big bundle) — the memory bound + back-pressure.
- The cap threshold + seed scheme go into the `.ssr.psp` header
  `extraction_params` so the subsample is reproducible (arch §8.3/§10).

**Hard requirement:** input must be indexed (`.bai`/`.csi`/`.crai`); error out
cleanly if missing — no whole-file-scan fallback (arch §8.1).

---

## 8. `locus_record.rs` + container schema — aggregation & write (arch §8.2/§11)

```rust
/// Per-locus reduction over one bundle's worker outputs → the columnar record.
/// confident → histogram; ambiguous slow → sparse CSR; off-ladder → offl_*.
fn aggregate(locus: &Locus, qc: LocusQc, outcomes: Vec<ReadOutcome>) -> SsrLocusRecord;

enum ReadOutcome {
    ConfidentOnLadder { units: u16, weight: f32 },          // → hist_*
    AmbiguousOnLadder { profile: Vec<(u16, f32)> },         // → amb_* CSR (after prune)
    ConfidentOffLadder { seq: NormalizedSeq, weight: f32 }, // → offl_*
    AmbiguousOffLadder { /* off-ladder leg + rungs */ },    // → offl_amb_*
}
```
Sparse storage rule (arch §11): a slow read's forward is evaluated over `2W+1`
lengths but stored sparse — drop candidates `> AMB_LL_DROP` below the per-read
max, renormalize over survivors (typically 1–3), write CSR
(`amb_read_offsets` / `amb_lengths` / `amb_logliks`). All stored log-liks are
**stutter-free** (the invariant); lengths are `uint16`.

Container side (task 3 if not already landed): `SsrLocusRecord` + the
`registry_ssr` schema over the generic columnar container (arch §10.4) — the
`hist_*` / `amb_*` / `offl_*` columns of spec §4.3. `--block-window-bp` is the
output RSS lever, as on the SNP path.

---

## 9. `mod.rs` — driver & CLI (arch §8.2)

```rust
pub struct SsrPileupArgs {
    inputs: Vec<PathBuf>,       // one or more BAM/CRAM — all the SAME sample (lanes/reps)
    catalog: PathBuf,           // Stage 0 .bed.gz + index
    reference: Option<PathBuf>, // CRAM-decode ONLY
    output: PathBuf,            // sampleN.ssr.psp
    threads: usize,             // WITHIN-sample pool (cross-sample = user orchestration)
    regions: Option<PathBuf>,
    block_window_bp: u32,
    // extraction params surfaced as overrides (recorded in the header):
    // min_flank_bp, flank_bp, w, max_reads_per_locus, min_mapq, ...
}

pub fn run(args: SsrPileupArgs) -> Result<(), SsrPileupError>;
```

`run()` flow (the §8.2 pipeline):
```
1. open AlignmentMergedReader over inputs (shared CRAM Repository, k-way merge);
   reject inputs whose read groups disagree on sample; require the index.
2. open CatalogReader (or .query() when --regions); read SsrParams + header.
3. spawn 1 fetcher thread (fetch_reads.rs) → bounded queue.
4. worker pool (threads): per LocusBundle, no BAM access — for each spanning read:
      triage(coverage+centre) → candidate_generation(rungs+offladder) → score_candidates
      → ReadOutcome[] → aggregate → SsrLocusRecord
   (v1 realigns every spanning read; the count_repeats shortcut, when measured-in,
    would handle confident reads before candidate_generation.)
   each worker owns a PairHmmScratch + candidate scratch buffers (reused per read).
5. ordered collector: reorder finished records to (chrom, start) before the
      registry_ssr writer (determinism + the block grid).
6. write header (extraction_params: all named consts + seed scheme) + finalize.
```

**Determinism invariant** (arch §8.4, a regression gate): byte-identical
`.ssr.psp` across `--threads` — fetcher fixed coord order + total read order +
ordered collector + deterministic forward + per-locus reservoir seed.

```rust
pub enum SsrPileupError {
    MissingIndex(PathBuf), SampleMismatch{..}, Catalog(..),
    Bam(..), Fasta(..), Io(std::io::Error),
}
```

---

## 10. Parameters (arch §10 — named consts, recorded in `extraction_params`)

`MIN_FLANK_BP=5`, `FLANK_BP=30 (=bundle_threshold)`, `MIN_BASE_QUAL`,
`STUTTER_WINDOW_UNITS (W)`, `AMB_LL_DROP`, `MIN_MAPQ=20` (shared SNP reader),
`MIN_MOTIF_RUN_BY_PERIOD[]` (~5/4/3 for period 2/3/≥4), `MAX_READS_PER_LOCUS (K)`,
`PAIR_HMM_BAND_BP (~7)`, `DINDEL_GAP_OPEN[]` / `GAP_EXTEND_PROB (e^-1)`. All become
`const`s with CLI overrides, written into the header. **Values are validation
targets, not design** (arch §14): the fast/slow fraction is the one load-bearing
*assumption* to measure on real target genomes.

---

## 11. Tests (arch §12 — the safety net built alongside)

- **Unit — pair-HMM in isolation:** a known read × known haplotype with a
  hand-computable forward score; gap-open homopolymer indexing; the rescale guard.
- **Unit — fast/slow gate** on boundary cases (each `SlowReason` reachable).
- **Unit — ladder:** `to_sequence ∘ repeat_count` round-trip; off-ladder
  normalization gives identical `NormalizedSeq` for the same allele built two ways;
  a normalized tract that is a pure tiling collapses to a rung (no off-ladder).
- **Read/BAM-level (crate-internal simulator, not a subcommand):** synthesize
  reads from known genotypes with injected sequencing error **and deliberately
  soft-clipped long alleles** (the §3.2 failure mode); run `ssr-pileup`; assert
  recovered `hist_*` / `amb_*` / `offl_*` match the injected lengths. The only way
  to exercise soft-clip recovery + the forward together.
- **Determinism gate:** byte-identical `.ssr.psp` across `--threads ∈ {1, N}`,
  *including* the depth-cap subsample (per-locus seed + total read order).
- **Regression gate (task 2):** SNP end-to-end tests stay green after the
  `normalize_alleles` lift.

> **Deferred (arch §12/§14):** external truth-set accuracy / HipSTR–GangSTR
> cross-tool concordance — simulator + unit tests prove the stage recovers what it
> was *given*, not real-data accuracy; the named benchmark is pinned before any
> precision claim, tracked as a known gap, not a build blocker.

---

## 12. What this sketch deliberately leaves vague

- **`MergedRecord` accessor surface** the worker needs (sequence incl. soft-clip,
  per-base quals, CIGAR ends) — pinned against `alignment_input.rs` when coding.
- **Exact shared-`normalize` module path** and the SSR adapter's `(seqs, bounds)`
  construction (§2) — the contract is fixed, the wiring lands at code time.
- **Cluster-fetch bookkeeping** (one reservoir per locus over a shared cluster
  query, arch §8.2) — built only if close-locus clusters are common enough to
  matter; the common far-apart case is one-locus-per-query.
- **Shared-flank lattice** pair-HMM optimization (arch §5.5) — built only if
  profiling shows the slow path binds; independent narrow-band forwards first.
- **Const *values*** beyond the decided `MIN_FLANK_BP=5` / `FLANK_BP=30` — defaults
  in place, calibrated on real data (arch §14), not chosen here.
