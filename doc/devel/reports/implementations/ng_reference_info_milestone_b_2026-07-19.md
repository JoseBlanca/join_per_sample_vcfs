# ng reference_info — Milestone B (the FASTA streaming pass)

*Implementation report, 2026-07-19. Plan-driven (Milestone B: B1+B2 one commit, B3 its own),
continuing the same implement → review → apply → commit loop. Design authority: spec
[`../../ng/spec/reference_info.md`](../../ng/spec/reference_info.md) (§4, §3.3, §3.4, §3.8)
and arch [`../../ng/arch/reference_info.md`](../../ng/arch/reference_info.md). All code in
[`src/ng/reference_info.rs`](../../../../src/ng/reference_info.rs).*

## Plan

Build the algorithmic heart: the from-byte-zero FASTA streaming pass. B1 reads
`Fasta { fai: None }` in one buffer (geometry + per-contig and whole-reference MD5); B2 is
its `samtools`/`.cat`/`faidx.5` oracle, landing in the same commit as its guard; B3 adds the
`Fasta { fai: Some }` fasta-vs-`.fai` cross-check. The silent-failure steps (a moved offset
or wrong M5 is a wrong number, not a crash) each land as their own commit with the oracle
green.

## Recorded deviations / absorbed additions

- **Compression detection landed in B1** (not named in the B1 bullet): a first-window gzip
  magic check → `CompressedReference`. Spec §3.8/§6 want a compressed reference to be a clean
  error, and the `CompressedReference` variant (declared in A1) needed a producer. Three
  lines, within the Fasta reader's scope.
- **B1 review fixes**: added the missing empty-contig-name guard (T3), and converted
  `finalize_contig`'s provably-unreachable `line_bases0.expect` to a graceful error — the
  module never panics on a supplied reference.
- **B3 review fix**: read the cheap `.fai` before the expensive genome pass, so a
  missing/bad index fails fast.
- **The pass is a hand-rolled byte-level state machine**, not a call to any existing reader —
  spec §4 rejects each one (the `.fai`-driven fetchers are circular; `sequence_reader` drops
  geometry; `Indexer` hides the bases). It copies `compute_contig_md5_streaming`'s streaming
  *shape* (fixed buffer, batched MD5) with the `isgraph` predicate.

## Changes made

- **B1** (`7ba63e5`) — `read_fasta` + `FastaPass` (Start/Header/Sequence state machine): from
  byte zero, one 64 KiB buffer, never a whole contig. One predicate — a base is `[0x21,0x7E]`
  — defines the M5, `line_bases`, and `length`; `\r`/space/tab fall out of the count and the
  MD5, so CR-LF and trailing whitespace need no special case (the deliberate divergence from
  production's space/tab-hashing rule, §3.4). Geometry reconstructed htslib-style (first line
  fixes it, later lines must match except a shorter last). Guards, all named errors:
  bases-before-`>`, empty contig, empty contig name, non-uniform interior line, duplicate
  name (T2), compressed `.fa.gz`. Wires `Fasta { fai: None }`.
- **B2** (`7ba63e5`, with B1) — the oracle: per-contig M5 == committed `samtools dict`; LN +
  reconstructed `.fai` == committed `samtools faidx`; whole-reference digest == the golden
  `.cat` header (`8a23ad24…`); the `faidx.5` man-page vector reproduced byte-for-byte under
  LF **and** CR-LF; per-contig MD5 == one-shot `Md5::digest`; the space/tab predicate edge.
- **B3** (`8766061`) — `read_fasta_verifying` + `first_fai_field_disagreement`: field-for-field
  comparison of the pass's reconstruction against the supplied `.fai`, first-difference
  ordered name → geometry → `FastaFaiMismatch`. Wires `Fasta { fai: Some }`.

## Tests added (18 new; 31 total in the module)

B1/B2: the golden index + per-contig M5; the reference digest vs the golden `.cat`; the
`faidx.5` worked example under both line endings; streaming-vs-one-shot MD5; the space/tab
predicate; and the guards (bases-before-header, non-uniform line, empty contig, empty name,
compressed, duplicate, last-line-without-newline). B3: a matching `.fai`; a re-wrap naming
`line_bases` (mutation-verified against a names-only check); a single-contig re-wrap
(offset unchanged); a reordering caught on `name` + the digest differing (T1); a count
mismatch; a missing supplied `.fai` (T4).

## Validation

Every step in the container: `cargo fmt --check` clean; `cargo clippy --lib --tests
--all-features -D warnings` clean; `cargo test` green (full lib suite **1922 passed** after
B3). Cross-checked on the host (cargo 1.95.0): the 31 module tests build and pass natively.

**Pre-existing, unrelated** (unchanged from Milestone A): `cargo clippy --all-targets
-D warnings` fails on `examples/ssr_psp_seqdump.rs:41` (`unnecessary_sort_by`, toolchain 1.95
bump; outside `src/ng`; left untouched per the freeze). The library and its tests are clean.

## Follow-ups (Milestone C onward)

- `write_fai` (C) — the writer, byte-identical to `samtools faidx`, atomic.
- The single-flight cache (D), the two entry points (E).
- Perf: the pass is byte-at-a-time over a match; correctness-first and cached, so untuned.
  If a single uncached read ever proves too slow, the deferred two-phase parallel hash
  (spec §7) is the lever — not needed now.
