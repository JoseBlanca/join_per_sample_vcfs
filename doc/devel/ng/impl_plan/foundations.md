# ng foundations — implementation plan (`ref_seq` + read filtering)

**Status:** draft, 2026-07-13. The build order for the **first ng code**: the `src/ng/`
skeleton, the `types.rs` seed, the `RefSeq` reference accessor, and step 1 (read
filtering). Design is settled in the specs — [`ref_seq.md`](../spec/ref_seq.md) and
[`read_filtering.md`](../spec/read_filtering.md) — and the arch docs
([step interfaces](../arch/ng_step_interfaces.md), [module layout](../arch/module_layout.md)).
This roadmap turns that design into build order; it is **not** a place for new design.

---

## Scope

**In:** the `src/ng/` skeleton; `types.rs` (seeded incrementally); `ref_seq.rs` with the
`RefSeq` + `RawRefSeq` traits and **two** impls — `InMemoryRefSeq` (synthetic, tests) and
`ResidentRefSeq` (real FASTA, reuse) — plus `read/filtering.rs` (step 1).

**Deferred, with homes (owner decision 2a — build what the foundations need, not
everything the context makes cheap now):**

- **`WindowedRefSeq`** — the caller-evictable streaming impl (`ref_seq.md` §6). Nothing in
  the foundations needs it: read filtering #8 uses a *raw-capable resident* ref. It lands
  with the **pileup** step, which is its first consumer.
- **`pileup/`, the locus router, and the research steps** — later plans.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every step (project rule).
- **Dependency order: `ref_seq` before read filtering.** Read filtering's mismatch filter
  (#8) consumes a `RawRefSeq`, so the accessor is built first. Within `ref_seq`, the
  synthetic `InMemoryRefSeq` comes before `ResidentRefSeq` so the traits get tests with no
  FASTA on disk.
- **Reuse over rewrite.** `ResidentRefSeq` stands on `fasta::Repository` /
  `RepositoryRefFetcher` / `RawContigRefCache`; read filtering calls the existing pure
  predicates (`read_exceeds_mismatch_fraction`, `cigar_is_bad`) and `FLAG_*` / `DEFAULT_*`
  consts. See the reuse maps in both specs.
- **Foundations set the conventions.** This is the first ng code, so the newtype rules,
  the config pattern, the file-vs-folder rule, and the test shape land here correctly
  (read_filtering.md §2) — getting them right on easy code de-risks the hard steps.
- **Verify against ground truth, not just self-consistency.** `ResidentRefSeq` is proven
  by **byte-parity with the production fetcher** on a FASTA fixture; read filtering by
  **parity with the production predicates** it ports. (read_filtering.md §2.6.)
- **Incremental, with pauses.** Land one milestone, stop for review, then the next — not
  the whole module at once.
- **Ungated.** `ng` compiles as a plain module (no `cargo` feature) for now; gate only if
  compile time later bites (module_layout *Open items*).
- **Container builds.** All `cargo` runs via `./scripts/dev.sh` (CLAUDE.md).

## Preconditions (already in place)

The production reference stack (`src/fasta/fetcher.rs`, `build_fasta_repository`) and the
per-read filter predicates (`src/bam/alignment_input.rs`) exist and are the reuse targets.
`MappedRead` and `ContigList` are reused as-is.

---

## The steps

### Milestone A — skeleton + `RefSeq` traits + `InMemoryRefSeq`

**A1. Scaffold + `types.rs` seed.**  ☐
`src/ng/mod.rs` (declarations + re-exports, minimal) and `pub mod ng;` in `lib.rs`.
`types.rs` seeded with **only** what `ref_seq` needs first: `ContigId` and `RefSeqError`
(`#[non_exhaustive]`). Nouns/errors only, no logic. *Source:* ng_step_interfaces §1;
ref_seq.md §The trait.

**A2. The traits.**  ☐
`RefSeq` (universal canonical `fetch_into` + owned `fetch`) and `RawRefSeq: RefSeq`
(`fetch_raw`) in `ref_seq.rs`. Trait definitions + docs only. *Source:* ref_seq.md
§The trait. *Note:* capabilities are separate traits by decision — no silent no-ops.

**A3. `InMemoryRefSeq` + tests.**  ☐
The synthetic impl (`Vec<Vec<u8>>` by `chrom_id`), implementing `RefSeq` + `RawRefSeq`.
Unit tests: canonical fetch folds to `{A,C,G,T,N}`, raw returns stored bytes, bounds/start-0
error paths. *Source:* ref_seq.md §Implementations.

> **Checkpoint A:** `cargo test` green; the traits are usable and testable with a hand-built
> reference, no FASTA. Pause for review.

### Milestone B — `ResidentRefSeq` (real FASTA, reuse)

**B1. Canonical resident fetch + `clear`.**  ☐
`ResidentRefSeq` wrapping `fasta::Repository` + `ContigList` (reuse `build_fasta_repository`);
`fetch_into` / `fetch` slice the resident contig and reuse `canonicalise`; `clear()` drops
the resident contig at a contig transition. *Depends:* A2. *Source:* ref_seq.md
§Implementations + Reuse map; the `--regions` caching rule.

**B2. `RawRefSeq` for `ResidentRefSeq`.**  ☐
Borrowed raw bytes (the `RawContigRefCache` model — a cached `Arc<Sequence>` reachable via
`&self`). *Source:* ref_seq.md Decision 4.

**B3. Tests — byte-parity with production.**  ☐
Against a small FASTA fixture, sweep `(chrom_id, start, length)` windows and assert
`ResidentRefSeq` canonical bytes match the production `RepositoryRefFetcher`, and raw bytes
match `RawContigRefCache`. Confirms the reuse is faithful. *Source:* read_filtering.md §2.6.

> **Checkpoint B:** real references read correctly, byte-identical to production. `ref_seq`
> is complete for the foundations (Windowed deferred). Pause for review.

### Milestone C — read filtering (step 1)

**C1. `types.rs` grows.**  ☐
Add `MapQual`, `BaseQual`, `Bp` (unconstrained, `pub` field) and `MismatchFraction`
(constrained, checked `try_new`) + `DomainError`. *Source:* read_filtering.md §2.2–2.3.

**C2. Filter config + result types.**  ☐
`ReadFilterConfig` (minimal — one field per active filter; `Default`; named `pub const`
defaults reused from `alignment_input`), `FilterVerdict { Keep, Drop(DropReason) }`,
`ReadFilterStats` (running tally). Types only. *Source:* read_filtering.md §4.

**C3. The `ReadFilter` iterator.**  ☐
`ReadFilter<I, R: RawRefSeq>` implementing `Iterator<Item = MappedRead>`: the hit-rate
cascade (filters #1–#7, #9 pure over `MappedRead`; #8 via `RawRefSeq::fetch_raw`, reusing
`read_exceeds_mismatch_fraction`), tallying every drop into the running `ReadFilterStats`.
`stats()` accessor. *Depends:* B (a `RawRefSeq` impl), C1–C2. *Source:* read_filtering.md §3, §5.

**C4. Tests.**  ☐
Unit: each filter's boundary (a read exactly at `min_mapq` kept, one below dropped) and the
cascade order/attribution; `InMemoryRefSeq` as the reference for #8. Integration: filter a
small known BAM/CRAM fixture and assert the `ReadFilterStats` matches hand-counted
expectations. *Source:* read_filtering.md §2.6, §3.

> **Checkpoint C:** read filtering runs end-to-end over a fixture; drop counts hand-verified;
> #8 parity with the production predicate. The ng foundations are in. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | unit tests on the synthetic impl (fold, raw, error paths) |
| B | **byte-parity** vs `RepositoryRefFetcher` / `RawContigRefCache` on a FASTA fixture |
| C | per-filter boundary + cascade unit tests; **fixture integration** with hand-counted stats; #8 parity with `read_exceeds_mismatch_fraction` |

## Out of scope (next plans)

- `WindowedRefSeq` — with the pileup step (its first consumer).
- `pileup/` (the non-STR loci+evidence generator), the locus router (step 3), and the
  research probes (candidate generation + realignment first — `ng_proposal.md` §3).
- `bench/` — first used by the candidate-generation probe, not the foundations.
