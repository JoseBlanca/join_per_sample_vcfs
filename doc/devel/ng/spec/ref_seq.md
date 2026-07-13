# ng — reference-sequence access (`RefSeq`)

*Status: design spec (2026-07-11). Companion arch: the `RefSeq` trait +
coordinate contract land in [`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md)
§1; the module home in [`../arch/module_layout.md`](../arch/module_layout.md). **No code
yet.** Grounded in a source read of the production reference stack — `src/fasta/fetcher.rs`,
`src/fasta/mod.rs`, `src/bam/alignment_input.rs`, and
[`../../implementation_plans/fasta_reference_reading_unification.md`](../../implementation_plans/fasta_reference_reading_unification.md).*

*Naming: **STR** in prose, `ssr` in code. This module is generic (reference bases are
marker-agnostic), so `ssr` does not appear.*

---

## Why this exists

Several ng steps need reference-genome bases for a `(contig, range)`: the read-filtering
**mismatch-fraction filter** (#8), the **pileup** (per-locus REF windows and widening),
**BAQ**, indel **left-alignment**, and **DUST**. Rather than each reaching for noodles
directly, they share one abstraction — **`RefSeq`**. It is specced and built *first*
because filter #8 depends on it (that is why `read_filtering.md` §3 flags #8 as blocked
on this doc) and because it is the first shared piece the locus stream leans on.

The production reference stack was iterated and perf-tuned hard — the `--regions` memory
blow-up, the whole-contig-resident vs streaming split, the raw-vs-canonical byte
asymmetry. This spec **reuses those hard-won lessons** rather than relearning them.

**Jargon, once.** The *reference* is the genome assembly the reads were aligned to; a
*contig* is one reference sequence (a chromosome or scaffold); *canonical* bases are the
sequence folded to uppercase `{A,C,G,T,N}` (any ambiguity code or lowercase soft-mask
becomes `N`/uppercase); *raw* bases are the verbatim FASTA bytes before that folding.

---

## Decisions (all confirmed with the owner, 2026-07-11)

### 1. Carry **both** access patterns — resident *and* streaming

Production has two irreducible ways to read reference bases, and ng keeps both:

- **Whole-contig-resident** — load an entire contig into memory once (`Arc<Sequence>`),
  then slice it by absolute position. Random-access; mandated by CRAM decoding (which
  indexes by absolute contig position); simplest; fast for small references and tests.
- **Streaming / sub-range** — hold only a sliding sub-range of the current contig in a
  small buffer, advancing as the walk moves forward. Lower memory *and* — the point that
  settled this — **faster**: for the forward genome walk it avoids materialising a
  250 MB contig, so first-access latency and footprint both drop, which speeds the test
  suite too.

We do **not** cut streaming for the lab. It is not only a memory optimisation, and —
decisively — **building both in ng means the reference code ports back into the
production var-caller without a rewrite.** ng's reference layer is meant to *become* (or
directly seed) the production one, so it carries both patterns behind one trait from the
start.

### 2. One trait, `RefSeq`

The three implementations sit behind a single trait:

- **resident** — reuses the noodles `fasta::Repository` machinery (as `RepositoryRefFetcher` does);
- **streaming** — a sub-range buffer with **caller-driven eviction** (Decision 6), reusing
  the `StreamingChromRefFetcher` / `ManualEvictChromRefFetcher` buffer logic;
- **in-memory synthetic** — a trivial `HashMap<ContigId, Vec<u8>>` for tests, so the
  mismatch filter and pileup are testable with a hand-built reference, no FASTA on disk.

**Name — `RefSeq`.** `Sequence` alone collides with noodles' `fasta::record::Sequence`;
`RefSeq` is unambiguous in the code. It also overlaps NCBI's reference-sequence
**database** name, but the owner has accepted that — it causes no real confusion in code —
so `RefSeq` is settled (not `RefBases`). The module is `ref_seq.rs` (§5); the trait matches it.

### 3. Reuse production's coordinate & contig contract verbatim

- A contig is identified by **`chrom_id: u32`** — an index into the canonical
  **`ContigList`** (`@SQ` / `.fai` order), never a name at the call boundary.
- Positions are **1-based** (ng's uniform convention — `ng_step_interfaces.md` §1); a
  fetch is `(chrom_id, start_1based, length)`, i.e. the `length` bases from `start_1based`
  onward.

Reusing `ContigList` / `ContigEntry` and this exact signature means the port-back is a
**no-op at this seam** — production consumers already speak it.

### 4. Byte handling — do exactly what the old caller does

Production exposes **both** raw and canonical bases, and different consumers want
different ones. ng matches that split, so the ported filters behave identically:

- **canonical `{A,C,G,T,N}` uppercase** for the pileup, BAQ, and DUST (the walker
  canonicalises today);
- **raw verbatim bytes** for the indel left-alignment + **mismatch-fraction** path —
  which is exactly what read-filtering #8 consumes (production's `RawContigRefCache`
  feeds the raw slice to `read_exceeds_mismatch_fraction`).

Matching production here — rather than simplifying to canonical-only — is a deliberate
choice **to facilitate the port back**: the mismatch filter reads the same bytes it will
read in production, so its behaviour is portable, not merely similar. `canonicalise` is
reused from `src/fasta`.

Consequence for the shape: **raw is naturally a resident-impl capability** (as in
production, where the raw path uses the resident `RawContigRefCache`, not the streaming
reader whose buffer is already uppercased). The streaming impl serves canonical bases for
the forward walk; a consumer needing raw bytes (the mismatch filter) uses the resident
impl. (Raw from the streaming/windowed impl is YAGNI — see *Decisions & deferrals*.)

### 5. Home — `src/ng/ref_seq.rs`

A single file to start; it splits into `ref_seq/` (trait in `mod.rs`, one file per impl)
when the three implementations grow enough to warrant it — the same file-vs-folder
convention `types.rs` and `read/` follow. It is **non-step infrastructure**, like
`pileup/`, `pipeline.rs`, and `bench/` — shared by many steps, owned by none.

### 6. Buffer lifecycle — caller-driven eviction (the production win)

Beyond "build once", production's other memory lesson is that **the consumer, not the
fetcher, usually knows when reference bytes are no longer needed** — and giving the
consumer an explicit way to free them was a real improvement
([`ManualEvictChromRefFetcher::evict_before`](../../../../src/fasta/fetcher.rs)). ng keeps
it first-class. There are **three** lifecycle shapes — my earlier draft flattened the
third into the second:

- **Whole-contig `clear()`** (resident) — drop the entire resident contig at a **contig
  transition**; the caller (the locus-stream driver) picks the point. This is the
  `--regions` fix.
- **Auto-slide** (streaming, the rigid form) — the buffer advances forward on its own and
  frees behind the window, but **forces monotonic-forward access** (a large backward jump
  is an error, `OutOfPattern`).
- **Caller-driven `evict_before(pos)`** (windowed) — the buffer extends in *either*
  direction on demand and shrinks **only** when the consumer says "nothing before `pos` is
  needed" (production `drain`s `[.., pos)`, keeping capacity). It tolerates **any access
  pattern within the resident range** and trusts the caller to bound memory. Production's
  Stage-1 BAQ uses it per worker: coordinate-sorted reads → `evict_before(read.pos)` after
  each read keeps just "current window + the next read's forward reach."

**ng's streaming impl uses caller-driven eviction, because it generalises auto-slide.** A
forward walk is simply `evict_before(current_pos)` after each step; random access *within*
a window (BAQ, pileup widening) just defers the evict. So the rigid monotonic auto-slide
is redundant for the lab — one caller-evictable windowed impl covers **both** the pileup's
forward walk and BAQ's within-window random access, without the monotonic straitjacket.
That is the shape the owner remembers as the big production improvement, and it earns its
place over the stricter reader.

So ng's impls, refining Decisions 1–2: **resident** (whole-contig, `clear()` at a contig
transition), **streaming/windowed** (caller-evictable via `evict_before`), and **in-memory
synthetic** (tests; holds everything, no eviction).

**`evict_before` is a method on the windowed impl, not a `RefSeq` trait method.** A trait
method that silently did nothing on the resident / in-memory impls would let a caller
*think* it freed memory when it did not — a silent no-op hides the fact. So eviction lives
where it means something: an inherent `WindowedRefSeq::evict_before`, and a consumer that
manages memory holds that concrete type. (If a second evictable impl ever appears, promote
it to an `EvictableRefSeq` capability trait — cheap, and still no silent no-op.)

---

## The trait — a design sketch (contract, not final API)

```rust
/// Access to reference-genome bases by (contig, 1-based range). The shared reference
/// source for the pileup, BAQ, DUST, left-alignment, and read filtering (#8).
///
/// The trait carries only the **universal** surface: canonical fetch, which every impl
/// provides. Capabilities only some impls have — raw bytes, buffer eviction — live on
/// those impls (a `RawRefSeq` sub-trait; an inherent `evict_before`), NOT here. A trait
/// method that silently no-ops on half the impls would hide that nothing happened, so the
/// universal surface stays honest (see *Implementations* + Decision 6).
pub trait RefSeq {
    /// Canonical bases {A,C,G,T,N} for the `length` bases from `start_1based`, written
    /// into `dst` — the alloc-free hot path (mirrors production's `fetch_into`).
    fn fetch_into(&self, chrom_id: ContigId, start_1based: u32, length: u32,
                  dst: &mut Vec<u8>) -> Result<(), RefSeqError>;

    /// Owned convenience over `fetch_into`. Canonical fetch returns *owned* bytes because
    /// canonicalisation produces fresh bytes (the raw read below can borrow instead).
    fn fetch(&self, chrom_id: ContigId, start_1based: u32, length: u32)
        -> Result<Vec<u8>, RefSeqError> { /* default: fetch_into a fresh Vec */ }
}

/// Raw, un-canonicalised bytes (borrowed) — the left-align / mismatch-fraction path (#8).
/// A capability of impls that keep bytes resident (`ResidentRefSeq`, `InMemoryRefSeq`);
/// the windowed impl simply does NOT implement it, so "no raw here" is a **compile-time**
/// fact, not a runtime `Err`. Read filtering #8 binds `R: RawRefSeq`.
pub trait RawRefSeq: RefSeq {
    fn fetch_raw(&self, chrom_id: ContigId, start_1based: u32, length: u32)
        -> Result<&[u8], RefSeqError>;
}
```

Two contracts to nail in the doc, both lifted from production:

- **Memory is bounded by eviction, not a monotonic-forward constraint** (Decision 6). The
  windowed impl extends its buffer either direction and stays bounded because the consumer
  calls its inherent `evict_before` — so it does *not* impose the rigid auto-slide's
  monotonic-forward precondition, and within-window random access (BAQ, pileup widening) is
  fine; the locus-stream walk just evicts as it advances. (`RefSeqError` keeps a reserved
  `OutOfPattern` variant in case a strict auto-slide impl is ever added.)
- **Build once, one resident contig, clear on transition.** See the caching section.

**On `&self` vs `&mut self`.** The trait is uniformly `&self`, so a resident impl stays
shareable read-only (the property production's rayon workers rely on, and the port-back
wants). The windowed impl's `fetch_into` still mutates its buffer, so it uses interior
mutability (a `RefCell`, as production's streaming reader does) to keep `&self`; its
`evict_before` is an inherent `&mut self` method. Nothing forces `&mut self` into the trait.

`RefSeqError` is `#[non_exhaustive]`, mirroring production's `ChromRefFetchError`:
out-of-bounds, invalid start (0), the reserved monotonic-contract violation, and I/O. (No
`RawUnavailable` — raw absence is now a compile-time capability fact, not a runtime error.)

---

## Implementations — the three distinct structs

Each impl owns a different memory/access trade-off. All provide the universal `RefSeq`
(canonical fetch); `RawRefSeq` and `evict_before` are held only by the impls that support
them — no silent no-ops.

### `ResidentRefSeq` — whole-contig, random-access, raw-capable

```rust
pub struct ResidentRefSeq {
    repository: fasta::Repository,   // reuse: noodles whole-contig cache (Arc<Sequence>)
    contigs: ContigList,             // chrom_id -> name
    raw_contig: RefCell<Option<(ContigId, Arc<Sequence>)>>,  // last contig, for borrowed raw
}
impl RefSeq    for ResidentRefSeq { /* fetch_into: slice the resident contig + canonicalise */ }
impl RawRefSeq for ResidentRefSeq { /* fetch_raw: borrow the resident bytes, no folding */ }
impl ResidentRefSeq {
    pub fn clear(&self);             // drop the resident contig at a contig transition (--regions fix)
}
```
Random access anywhere in the resident contig; `Send + Sync` (shareable read-only).
Bounded by `clear()` at contig transitions. **The raw-capable impl**, so read filtering #8
uses it. Reuse target: `RepositoryRefFetcher` + `RawContigRefCache`.

### `WindowedRefSeq` — sub-range buffer, caller-evictable (canonical only)

```rust
pub struct WindowedRefSeq {
    file: File, contigs: ContigList,
    current_window: RefCell<Option<ContigWindow>>,   // rebuilt when chrom_id changes
}
struct ContigWindow { chrom_id: ContigId, fai: ContigFai, buf: Vec<u8>, buf_start_1based: u32 } // uppercased
impl RefSeq for WindowedRefSeq { /* fetch_into: extend buf either direction, then copy out */ }
impl WindowedRefSeq {
    /// THE capability (inherent, not a trait method): drain [.., pos), keep capacity.
    pub fn evict_before(&mut self, pos: u32);
}
```
Canonical only (the buffer is uppercased → no `RawRefSeq`; Decision 4). Any access within
the resident window; memory bounded by `evict_before` (Decision 6), not a monotonic
constraint. Used by the pileup forward walk and BAQ within-window access. `Send` (per-
worker; `!Sync` via `RefCell` — fine, ng is single-thread). Reuse target:
`StreamingChromRefFetcher` / `ManualEvictChromRefFetcher`.

### `InMemoryRefSeq` — synthetic, tests

```rust
pub struct InMemoryRefSeq { contigs: Vec<Vec<u8>> }   // indexed by chrom_id; hand-built bases
impl RefSeq    for InMemoryRefSeq { /* fetch_into: slice + canonicalise */ }
impl RawRefSeq for InMemoryRefSeq { /* fetch_raw: slice the stored bytes */ }
impl InMemoryRefSeq { pub fn from_contigs(contigs: Vec<Vec<u8>>) -> Self; }
```
Holds everything (no eviction, no `clear` needed); random access; raw-capable. Lets the
mismatch filter (#8) and the pileup be unit-tested with **no FASTA on disk**.

### Capability matrix

| impl | canonical `fetch` (`RefSeq`) | `fetch_raw` (`RawRefSeq`) | random access | `evict_before` | `clear()` | thread |
|---|---|---|---|---|---|---|
| `ResidentRefSeq` | ✓ | ✓ | whole contig | — | ✓ (contig transition) | `Send + Sync` |
| `WindowedRefSeq` | ✓ | — (not impl'd) | within window | ✓ (inherent) | — (rebuild on chrom change) | `Send` |
| `InMemoryRefSeq` | ✓ | ✓ | whole (in memory) | — | — | `Send + Sync` |

A blank cell is a **compile-time** absence (the impl doesn't provide that method), never a
silent no-op. A consumer requiring a capability binds on it (`R: RawRefSeq`) or holds the
concrete type (`WindowedRefSeq` for eviction).

---

## Reuse map — to `src/fasta`

| ng piece | production code | reuse |
|---|---|---|
| resident `RefSeq` impl | `RepositoryRefFetcher` ([fasta/fetcher.rs](../../../../src/fasta/fetcher.rs)) + `build_fasta_repository` ([bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)) | wrap / lift — it already does chrom_id → name → `Repository::get` → slice → canonicalise |
| raw resident access | `RawContigRefCache` ([pileup/per_sample/read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)) | the `fetch_raw` model |
| streaming `RefSeq` impl | `StreamingChromRefFetcher` ([fasta/fetcher.rs](../../../../src/fasta/fetcher.rs)) | the sliding-buffer / both-direction-extend logic |
| caller-driven eviction | `ManualEvictChromRefFetcher::evict_before` ([fasta/fetcher.rs](../../../../src/fasta/fetcher.rs)) | the windowed impl's `evict_before` (drain `[.., pos)`, keep capacity) — Decision 6 |
| coordinate / contig contract | `ContigList`, `ContigEntry` ([fasta/mod.rs](../../../../src/fasta/mod.rs)); the `(chrom_id, 1-based, length)` signature | reuse as-is |
| byte folding | `canonicalise` ([fasta/fetcher.rs](../../../../src/fasta/fetcher.rs)) | reuse as-is |
| the caching discipline | the `--regions` fix (build once, `clear()` on contig transition) | adopt the *rule*, not the plumbing |

The production traits themselves — `ChromRefFetcher` (sealed, single-contig) and
`MultiChromRefFetcher` (multi-contig) — are what `RefSeq` is meant to **consolidate** into
one clean name. That consolidation is the §6-reconciliation payoff of doing this in ng
and porting back (deferred to port-back; see *Decisions & deferrals*).

---

## Caching & the `--regions` lesson (do not relearn it)

Production once built a fresh `fasta::Repository` **per region**, each rebuild reloading
the whole ~250 MB contig → 14.6 GB peak RSS on a 1000-region run. The fix: **build the
reader once per run, keep one contig resident, and `clear()` on contig transition**,
relying on the sorted-`RegionSet` invariant (contigs grouped, ascending, never revisited)
— which dropped peak RSS 95%.

ng's single-phase **forward locus walk** makes this discipline natural: one contig
resident at a time, advanced as the walk crosses contigs. The spec's one hard rule:
**do not reintroduce the per-region rebuild.** ng does not need production's full
machinery, but it inherits the rule.

---

## Decisions & deferrals

Nothing is open. Records of what was settled and what is intentionally left for later:

- **Trait vs capabilities — decided.** `RefSeq` carries only the universal canonical
  fetch; **raw** is a `RawRefSeq` sub-trait, and **eviction** is an inherent
  `WindowedRefSeq::evict_before` — capabilities on the impls that have them, so a missing
  capability is a compile-time fact, never a silent no-op. An `EvictableRefSeq` trait is
  deferred until a second evictable impl needs it.
- **Coordinate base — decided: ng is 1-based**, matching production, VCF/SAM/IGV, and
  minimising conversion seams. `ng_step_interfaces.md` §1's `Position` was revised from
  0-based to 1-based to match; only BED input converts, at the boundary.
- **Trait name — decided: `RefSeq`** (the NCBI-database overlap is consciously accepted;
  Decision 2).
- **Raw from the windowed impl — YAGNI.** `WindowedRefSeq` is canonical-only by design; add
  a `RawRefSeq` impl (with an un-canonicalised buffer) only if a windowed consumer ever
  actually needs raw bytes.
- **Production adopting `RefSeq` at port-back — deferred.** The ambition is that `RefSeq`
  *replaces* the two production reference traits (an `ng_step_interfaces.md` §6
  consolidation) rather than just wrapping them; settle it when the port-back happens.
```
