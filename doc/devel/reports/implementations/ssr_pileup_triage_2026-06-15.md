# SSR Stage-1 — `triage` read classification (realign-everything)

**Date:** 2026-06-15 · **Skill:** rust-feature-implementation · **Branch:** `ssr-architecture`

The trust-critical read-classification stage, built to the **realign-everything**
design ([arch §2 revision](../../../doc/devel/architecture/ssr_pileup.md), plan §4):
coverage classification → region extract → window centre, with no CIGAR fast/slow
gate and no trust in the mapper's indel placement.

## 1. What was built

In [src/ssr/pileup/triage.rs](../../../src/ssr/pileup/triage.rs) (atop the
existing `find_longest_stretch`):

- `Footprint` + `read_footprint(cigar, pos)` — the read's reference footprint
  from CIGAR op-lengths + soft-clips only (`ref_start`, `ref_end`, `leading_clip`,
  `trailing_clip`); never the indel placement.
- `brackets(fp, locus) -> (bool, bool)` — does the footprint, extended by each
  side's soft-clip (the optimistic long-allele reach), reach ≥ `MIN_FLANK_BP` past
  the tract on each side? Signed arithmetic (no underflow near a contig end).
- `ref_to_read` + `extract_region` — the read-coordinate span covering the locus's
  embedded reference window; where the window extends beyond the alignment, the
  whole soft-clip on that side is grabbed (that's where a long allele's excess +
  far flank live), and the pair-HMM realigns within.
- `TriageResult { Spanning(SpanningRead) | Flanking | InRepeat }` + `SpanningRead
  { region, observed_count }`.
- `triage_read(read: &MappedRead, locus) -> TriageResult` — footprint → coverage →
  (for spanning) extract region + centre via `find_longest_stretch`.
- `MIN_FLANK_BP = 5` (arch §9, calibration placeholder).

## 2. Decisions applied (from the design discussion)

- **Coverage, not flank-match** (the agreed policy): spanning is decided by
  footprint position + clip lengths; no flank-byte comparison at triage. The
  pair-HMM's flank emission judges flank quality ("triage dumb, HMM smart").
- **Clip = optimistic reach**: a soft-clip on a side extends that side's reach by
  its length, so a read whose far flank is in the clip still classifies Spanning;
  extraction then grabs the whole clip.
- **Read seam = `MappedRead`** (the sanctioned SNP-reader reuse, arch §3.1).
  `triage_read` takes `&MappedRead`; the inner helpers take primitives
  (`cigar`/`pos`) so they're unit-tested without building full records.
- **`observed_count` = pre-probe longest contiguous run** (the window centre).

## 3. The realign-everything payoff, demonstrated

`triage_recovers_a_soft_clipped_long_allele_via_the_clip`: a read the mapper
aligned as **3 units + a 12 bp soft-clip** is classified Spanning and its
`observed_count` recovered as the true **6 units** — from the pre-probe over the
clip-included region, never trusting the CIGAR's 3. This is exactly the
mapper-mis-alignment case the realign-everything pivot targets.

## 4. Tests added (14; 22 in the module total)

Footprint (plain match, both-end soft-clips, deletion extends ref-span-only,
hard-clip-then-soft-clip); `brackets` (fully-aligned spanning, one-side flanking,
buried-in-tract neither, clips-extend-reach); `extract_region` (full alignment,
window-opens-left-of-alignment grabs the clip); `triage_read` (clean spanning
centres on the tract, in-repeat dropped, flanking dropped, soft-clipped long
allele recovered).

## 5. Validation results

Dev container (`./scripts/dev.sh`):
- `cargo fmt -- --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ssr::pileup::triage` — **22 passed, 0 failed**.

## 6. Tradeoffs and follow-ups

- Coverage is optimistic for clipped reads by design; a clipped read that doesn't
  truly span will score poorly / off-centre at the pair-HMM — acceptable, and the
  alternative (flank-byte matching at triage) was deliberately rejected.
- `triage_read` returns a `region` `Range` into `seq`/`qual`; the worker slices
  both. No per-read allocation beyond the candidate set.
- **Next:** the worker wiring — `triage_read` → `candidate_generation::build_rungs`
  (+ `build_offladder`) → `pair_hmm::score_candidates` → `ReadOutcome`, then the
  aggregation (`locus_record`, which needs the container schema) and `fetch_reads`.
- Deferred + measured: the `count_repeats` fast-path shortcut.
