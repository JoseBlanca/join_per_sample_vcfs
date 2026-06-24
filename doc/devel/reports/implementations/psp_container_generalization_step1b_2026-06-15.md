# `.psp` container generalization — step 1b (reader's typed path on the schema)

**Date:** 2026-06-15
**Branch:** `ssr-architecture`
**Plan:** [psp_container_generalization.md](../../../doc/devel/implementation_plans/psp_container_generalization.md)
(architecture §10 of [ssr_genotyping_architecture.md](../../../doc/devel/architecture/ssr_genotyping_architecture.md))
**Builds on:** [step 1a](psp_container_generalization_step1a_2026-06-15.md)

## Goal

Step 1b: make the reader's **typed record path** (`RecordsIter`) generic over
the container schema `S`, mirroring the writer (1a). Behaviour-preserving; SNP
the only impl. Acceptance was raised beyond "e2e green" at the user's request:
the produced `.psp` must stay **byte-identical** to the pre-refactor baseline,
and records must **round-trip** through the generalized reader.

## Design — symmetric with the writer

The discussion settled that reader and writer are symmetric, and the earlier
worry ("the typed reader is too perf-tuned to generalize") was misplaced: the
tuning is *SNP decode scratch*, which simply moves behind a schema hook, exactly
as the SNP column `Vec`s live in `SnpBlock` on the write side.

| | Writer (1a) | Reader (1b) |
|---|---|---|
| Shared core | block index/header/sink; **compress + write** columns | block index/header/source; **read + decompress** columns |
| Schema state | `S::Block` (`BlockAccumulator`) | `S::Decoder` (`BlockDecoder`) |
| Record bridge | `BlockAccumulator::append`: record → columns | `BlockDecoder::next_record`: columns → record |

This satisfies arch Q3 ("each schema keeps its own closed decode") — `SnpDecoder`
is its own closed decode; a future `SsrDecoder` is another. The shared driver
only shares decompression, never merges the two decodes.

## What changed

### [src/psp/kind.rs](../../../src/psp/kind.rs)
- `PspKind` gains `type Decoder` (**unbounded** in the `pub` trait — see
  visibility note) and `fn record_coord(&Record) -> (u32, u32)` (for the region
  clamp; the §10.5 interval step widens the right edge later).
- New `pub(crate) trait BlockDecoder` (mirror of `BlockAccumulator`):
  `new_decoder`, `decode_block` (shared decompression scratch handed in),
  `next_record`.

### [src/psp/reader.rs](../../../src/psp/reader.rs)
- New `pub struct SnpDecoder` (`SnpKind::Decoder`) — owns the `DecodedBlock`,
  the two cross-block CSR slabs (the H1 allocation collapse, **moved off**
  `RecordsIter`), and the per-block record cursor. `impl BlockDecoder`:
  `decode_block` wraps the **unchanged** `decode_block_payload`;
  `next_record`/`materialise_next_record` is the unchanged SNP materialization.
- `RecordsIter<'r, R, S = SnpKind>` is now generic: it holds the generic
  decompression scratch + `S::Decoder` and drives block framing + the clamp via
  `decoder.next_record()` + `S::record_coord`. `records()` / `region_records()`
  keep returning `RecordsIter<'_, R>` (default `SnpKind`) — call sites unchanged.
- `RangeClamp` predicates are now coordinate-based (`coord_past_window` /
  `coord_before_window`) instead of `&PileupRecord`-typed.
- **`BlockColumnReader` / `BlockColumns` / `decode_block_payload` / the
  two-phase decode are untouched** — the cohort columnar hot path keeps its exact
  behaviour. `decode_block_payload` is now shared between the (unchanged)
  columnar reader and `SnpDecoder`.

### [src/psp/writer.rs](../../../src/psp/writer.rs)
- `impl PspKind for SnpKind` gains `type Decoder = super::reader::SnpDecoder` and
  `record_coord` (returns `(chrom_id, pos)`).

## Decisions / notes

- **Scope = typed path only.** `decode_block_payload`/`DecodedBlock` are shared
  by both `RecordsIter` and the cohort `BlockColumnReader`. Rather than risk the
  perf-critical cohort path, 1b wraps that shared decode behind `SnpDecoder` and
  leaves the columnar reader byte-for-byte as-is. Generalizing the columnar view
  to a second schema is future work (and arch Q2 doubts SSR needs a full typed
  reader anyway — SSR Stage-2 rides the columnar layer).
- **Visibility:** `PspKind::Decoder` is left **unbounded** in the `pub` trait so
  the trait surface doesn't leak the `pub(crate)` `BlockDecoder` (whose method
  signatures name internal wire types — `BlockHeader`). The
  `Decoder: BlockDecoder` bound lives only on `RecordsIter`'s impls, carrying a
  documented `#[allow(private_bounds)]` (sealed-style idiom: the bound is an
  implementation detail callers can't observe — the public surface is just
  `Iterator<Item = Result<Record>>`). This keeps the wire types encapsulated
  rather than `pub`-re-exporting them. (`SnpDecoder` itself is `pub` because it
  is the `SnpKind::Decoder` associated value — same reason `SnpBlock` is `pub`.)

## Verification

Built a deterministic public-API `.psp` producer (multi-block via window cuts,
multi-allele, chain-id lists, a chromosome boundary) and ran it on the
pre-refactor baseline (`aa6a105`, via a worktree) and on HEAD:

- **Byte-identity:** `cmp` of the produced `.psp` — **identical** (writer is
  untouched across 1a+1b, so this is the strong confirmation that the refactor
  changed implementation only).
- **Round-trip:** the exact produced bytes read back through the generalized
  `RecordsIter<SnpKind>` — **11/11 records equal** the written input.

The producer + baseline worktree were throwaway and have been removed; step 5
will add a committed round-trip test.

## Gates (dev container)

- `cargo build --lib` / `clippy --all-targets -- -D warnings` / `fmt --check` — clean
- `cargo test --lib` — **1133 passed, 0 failed, 1 ignored**
- SNP e2e: `psp_to_pileup_integration` (4), `pileup_cli_integration` (10) — green

## Next (remaining §10.6 steps)

- **2** — `kind` header tag; reader selects the registry by `kind`;
  `cross_check_against_registry` validates against it (`SnpKind::KIND` already wired).
- **3** — interval-overlap block index (`last_pos` = block max-`end`).
- **4** — `registry_ssr` + `SsrLocusRecord` (all-CSR: `amb_*` + QC scalars).
- **5** — synthetic `.ssr.psp` round-trip (write → columnar read → typed read).
