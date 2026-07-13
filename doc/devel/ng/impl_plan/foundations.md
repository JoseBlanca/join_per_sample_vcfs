# ng foundations — implementation plan (`RefSeq` reference accessor)

**Status:** draft, 2026-07-13. The build order for the **first ng code**: the `src/ng/`
skeleton, the `types.rs` seed, and the `RefSeq` reference-sequence accessor with all three
implementations. Design is settled in [`ref_seq.md`](../spec/ref_seq.md) and the arch docs
([step interfaces](../arch/ng_step_interfaces.md), [module layout](../arch/module_layout.md)).
This roadmap turns that design into build order; it is **not** a place for new design.

**Read filtering (step 1) is deliberately *not* in this plan** — its spec and architecture are
not yet finalized. The `RefSeq` accessor is the clean, fully-specified foundation to build
first; read filtering follows in its own plan once its design settles.

---

## Scope

**In:** the `src/ng/` skeleton; `types.rs` (seeded with what `ref_seq` needs); `ref_seq.rs`
with the `RefSeq` + `RawRefSeq` traits and **all three** impls:

- `InMemoryRefSeq` (synthetic, for tests),
- `ResidentRefSeq` (whole-contig, real FASTA, reuse),
- `WindowedRefSeq` (streaming sub-range, caller-evictable).

**Out (later plans):**

- **Read filtering (step 1)** — spec/arch not finalized; a separate plan once they are.
- **`pileup/`, the locus router, and the research probes** — later.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every step (project rule).
- **Simplest impl first, as the test oracle for the next.** `InMemoryRefSeq` (no FASTA) →
  `ResidentRefSeq` (reuse `fasta::Repository`) → `WindowedRefSeq` (the most complex — buffered
  and evictable). Each earlier impl is a byte-parity oracle for the next.
- **Reuse over rewrite.** `ResidentRefSeq` stands on `fasta::Repository` /
  `RepositoryRefFetcher` / `RawContigRefCache`; `WindowedRefSeq` on the
  `StreamingChromRefFetcher` sliding-buffer + `ManualEvictChromRefFetcher` evict logic. See
  ref_seq.md *Reuse map*.
- **Foundations set the conventions.** First ng code: the newtype rules (`ContigId`), the
  `#[non_exhaustive]` error pattern (`RefSeqError`), the capability-not-silent-no-op split
  (`RawRefSeq` sub-trait / inherent `evict_before`), the file-vs-folder rule, and the test
  shape all land here correctly (ng_step_interfaces §1; ref_seq.md). Getting them right on this
  well-specified module de-risks the harder steps.
- **Verify against ground truth.** Every real impl is proven by **byte-parity with the
  production fetcher** (and with the earlier ng impl) on a FASTA fixture, not just
  self-consistency.
- **Incremental, with pauses.** Land one milestone, stop for review, then the next.
- **Ungated.** `ng` compiles as a plain module (no `cargo` feature) for now; gate only if
  compile time later bites (module_layout *Open items*).
- **Container builds.** All `cargo` runs via `./scripts/dev.sh` (CLAUDE.md); a final native
  host build on completion.

## Preconditions (already in place)

The production reference stack — `src/fasta/fetcher.rs` (`RepositoryRefFetcher`,
`StreamingChromRefFetcher`, `ManualEvictChromRefFetcher`, `canonicalise`, `ContigFai`),
`build_fasta_repository`, and `ContigList` — exists and is both the reuse target and the parity
oracle. `MappedRead` / `ContigList` are reused as-is.

---

## The steps

### Milestone A — skeleton + `RefSeq` traits + `InMemoryRefSeq`

**A1. Scaffold + `types.rs` seed.**  ✅
`src/ng/mod.rs` (declarations + re-exports, minimal) and `pub mod ng;` in `lib.rs`. `types.rs`
seeded with **only** what `ref_seq` needs: `ContigId` (unconstrained newtype) and `RefSeqError`
(`#[non_exhaustive]`). Nouns/errors only, no logic. *Source:* ng_step_interfaces §1; ref_seq.md
§The trait.

**A2. The traits.**  ✅
`RefSeq` (universal canonical `fetch_into` + owned `fetch`) and `RawRefSeq: RefSeq` (`fetch_raw`)
in `ref_seq.rs`. Definitions + docs only. *Source:* ref_seq.md §The trait. *Note:* capabilities
are separate traits by decision — no silent no-ops.

**A3. `InMemoryRefSeq` + tests.**  ✅
The synthetic impl (`Vec<Vec<u8>>` by `chrom_id`), implementing `RefSeq` + `RawRefSeq`. Unit
tests: canonical fetch folds to `{A,C,G,T,N}`, raw returns stored bytes, bounds / start-0 error
paths. *Source:* ref_seq.md §Implementations.

> **Checkpoint A:** `cargo test` green; the traits are usable and testable with a hand-built
> reference, no FASTA. Pause for review.

### Milestone B — `ResidentRefSeq` (whole-contig, real FASTA, reuse)

**B1. Canonical resident fetch + `clear`.**  ☐
`ResidentRefSeq` wrapping `fasta::Repository` + `ContigList` (reuse `build_fasta_repository`);
`fetch_into` / `fetch` slice the resident contig and reuse `canonicalise`; `clear()` drops the
resident contig at a contig transition. *Depends:* A2. *Source:* ref_seq.md §Implementations +
Reuse map; the `--regions` caching rule.

**B2. `RawRefSeq` for `ResidentRefSeq`.**  ☐
Borrowed raw bytes (the `RawContigRefCache` model — a cached `Arc<Sequence>` reachable via
`&self`). *Source:* ref_seq.md Decision 4.

**B3. Tests — byte-parity with production.**  ☐
Against a small FASTA fixture, sweep `(chrom_id, start, length)` windows and assert
`ResidentRefSeq` canonical bytes match `RepositoryRefFetcher` and raw bytes match
`RawContigRefCache`. Also assert agreement with `InMemoryRefSeq` on a hand-built reference.
*Source:* ref_seq.md §Implementations.

> **Checkpoint B:** real references read correctly, byte-identical to production. Pause for review.

### Milestone C — `WindowedRefSeq` (streaming sub-range, caller-evictable)

**C1. The windowed buffer + canonical fetch.**  ☐
`WindowedRefSeq { file, contigs, current_window: RefCell<Option<ContigWindow>> }` with
`ContigWindow { chrom_id, fai: ContigFai, buf, buf_start_1based }`. `impl RefSeq`: `fetch_into`
**extends the buffer in either direction** on demand and rebuilds it on a `chrom_id` change,
canonicalising into `dst`; interior mutability (`RefCell`) keeps `fetch` `&self`. **Canonical
only — no `RawRefSeq` impl** (the buffer is uppercased). Reuse the `StreamingChromRefFetcher`
sliding-buffer + `ManualEvictChromRefFetcher` both-direction-extend logic. *Depends:* A2.
*Source:* ref_seq.md §Implementations, Decision 6, Reuse map.

**C2. Caller-driven eviction.**  ☐
Inherent `evict_before(&mut self, pos)` — drain `[.., pos)`, retain capacity (the
`ManualEvictChromRefFetcher::evict_before` model); no-op when `pos ≤ buf_start_1based` or the
buffer is empty; clear when `pos` is past the buffer end. *Source:* ref_seq.md Decision 6.

**C3. Tests.**  ☐
- **Byte-parity:** canonical bytes match `ResidentRefSeq` (B) and the production
  `StreamingChromRefFetcher` across `(chrom_id, start, length)` sweeps, including a contig
  transition (buffer rebuild).
- **Access patterns:** a forward walk; within-window random access (backward within the
  resident range); both-direction extend.
- **Eviction:** after `evict_before(pos)`, `buf_start_1based` advances and capacity is retained;
  re-fetching `≥ pos` still works; memory stays bounded across a simulated forward walk with
  eviction.

*Source:* ref_seq.md §Implementations, Decision 6.

> **Checkpoint C:** `WindowedRefSeq` reads byte-identically to Resident/production; eviction
> bounds memory; any access pattern within the resident range works. The `RefSeq` accessor is
> complete (all three impls). Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | unit tests on the synthetic impl (fold, raw, error paths) |
| B | **byte-parity** vs `RepositoryRefFetcher` / `RawContigRefCache` (and `InMemoryRefSeq`) on a FASTA fixture |
| C | **byte-parity** vs `ResidentRefSeq` + `StreamingChromRefFetcher`; access-pattern + eviction-bounds tests |

## Out of scope (next plans)

- **Read filtering (step 1)** — once its spec + architecture are finalized; a separate plan.
  (`RawRefSeq`, which #8 will consume, is fully built here — read filtering only adds its own
  module on top.)
- **`pileup/`** (the non-STR loci+evidence generator — `WindowedRefSeq`'s first real consumer),
  the **locus router** (step 3), and the **research probes** (candidate generation + realignment
  first — `ng_proposal.md` §3).
- **`bench/`** — first used by the candidate-generation probe, not the foundations.
