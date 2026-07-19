# ng reference info — implementation plan (`src/ng/reference_info.rs`)

**Status:** draft, 2026-07-18. The build order for the reference-info reader: the types, the
`.fai` and FASTA readers, the `.fai` writer, the single-flight cache, and the two composing entry
points. Design is settled in [`../spec/reference_info.md`](../spec/reference_info.md) (the *why*)
and [`../arch/reference_info.md`](../arch/reference_info.md) (the types & interfaces). This
roadmap turns that design into build order; it is **not** a place for new design — a step that
seems to need a decision goes back to the spec/arch, not decided here.

Builds on [`foundations.md`](foundations.md) (the `src/ng/` skeleton and `types.rs` already
exist). It does **not** touch production; ng owns divergent copies where it copies.

---

## Scope

**In:** `src/ng/reference_info.rs` — the types (`ContigInfo`, `ReferenceInfo`, `ReferenceSource`,
`ReferenceInfoError`, `VerificationHandle`, `ReferenceInfoCache`), the readers
(`read_reference_info` for both sources), `write_fai`, `sibling_fai_path`, the cache, and the two
entry points (`read_fai_verify_in_background`, `read_reference_verifying_or_creating_fai`), each
with its tests.

**Out (later plans / homes):**

- **`@SQ` ↔ `.fai` reconciliation** (closes `ReadFilter::new`'s permutation hole) → read
  ingestion's plan (arch §4 deferred).
- **Consolidating ng's fai-driven readers onto `ContigInfo`** (`RawChromReader` takes a
  `&ContigInfo`; retires the `contig_list()` bridge) → its own ng refactor (arch §4).
- **bgzip-compressed reference** and the **two-phase parallel hash** → error / deferred today
  (arch §4).

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The algorithmic heart before the plumbing.** The FASTA streaming pass (geometry + MD5, spec §4)
  is the heart; the cache, the background verify, and the orchestrator are plumbing that composes
  it. Build and prove the pass in isolation first.
- **Isolate the silent-failure steps.** An off-by-one in the geometry/offset reconstruction or a
  wrong M5 predicate is a *wrong `.fai`/M5*, not a panic. Those steps (B1, B2, C1) land as **their
  own commit with the `samtools` oracle green before and after** — never bundled — so `git bisect`
  can find a moved number.
- **Verify against ground truth, not self-consistency.** The north-star is **byte-parity with
  `samtools`** (`dict` for M5/LN, `faidx` for the `.fai`), resolved to committed constants so
  `samtools` is not a test-time dependency (arch §7).
- **Reuse over rewrite.** The pass copies `compute_contig_md5_streaming`'s streaming shape and
  `ContigFai::validate`'s guards; the `Fai` arm calls `noodles fai::fs::read`; hex via
  `format_md5_hex` (arch §5). ng owns the copies; production is untouched.
- **Incremental, with pauses.** One milestone, then stop for review.
- **Ungated, container builds.** `ng` compiles as a plain module; all `cargo` via `./scripts/dev.sh`
  (CLAUDE.md); a native host build on completion.

## Preconditions (already in place)

- The `src/ng/` skeleton and `types.rs` exist ([`foundations.md`](foundations.md)); adding the
  module is one `pub mod reference_info;` line.
- **Reuse targets** exist and are read-only: `compute_contig_md5_streaming`
  ([`common.rs:250`](../../../../src/pop_var_caller/common.rs)) and `format_md5_hex` (`:60`);
  `ContigFai::validate` ([`raw_chrom_reader.rs:77`](../../../../src/ng/raw_chrom_reader.rs));
  `ContigList`/`ContigEntry` ([`fasta/mod.rs:37-62`](../../../../src/fasta/mod.rs)); `noodles-fasta`
  0.61 (`fai::fs::read`, `io::Indexer`) and `md-5` ([`Cargo.toml:81`,`:93`](../../../../Cargo.toml)).
- **Oracle**: the vendored `samtools`/`htslib` (gitignored) and the golden reference
  `tests/data/tandem_repeat/synthetic_ref.fa` exist (arch §7).

---

## The steps

### Milestone A — types + the cheap (`.fai`) reader

**A1. Types + module scaffold.**  ✅
`pub mod reference_info;` in `ng/mod.rs`; the data types `ContigInfo`, `ReferenceInfo` (+ the
`contig_list()` projection), `ReferenceSource`, and the `#[non_exhaustive]` `ReferenceInfoError`
(all variants, per-variant docs). Nouns/errors only, no logic. *Source:* arch §1.

**A2. `sibling_fai_path`.**  ✅
Pure path helper (`<fasta>` + `.fai`, no I/O — copy `with_fai_extension`'s five lines). Test: it
appends `.fai` and touches no filesystem. *Depends:* A1. *Source:* arch §2, spec §2/§3.6.

**A3. Test fixtures + committed oracle constants.**  ✅
A width-configurable FASTA writer (the existing `ref_seq.rs::build_fasta` is one-line-per-contig,
so it cannot exercise geometry); run `samtools dict` + `samtools faidx` **once** on the golden
reference and **commit** the `M5`/`LN`/`.fai` values as constants, plus the `faidx.5` LF/CR-LF
worked example as a second vector. Test infra — pairs with A4 and guards all of B/C. *Source:*
arch §7, spec §6/§3.8.

**A4. `read_reference_info` — the `Fai` arm.**  ✅
`ReferenceSource::Fai` → `noodles fai::fs::read` → `Vec<ContigInfo>` (`md5: None`); reject duplicate
names (T2), a 6-column FASTQ index (§3.8), and a `.fai` failing the field guards (`line_bases > 0`,
`line_width ≥ line_bases` — copy `ContigFai::validate`). Tests: parse the committed `.fai`, assert
fields; the three rejections error (named). *Depends:* A1, A3. *Source:* spec §4 (Fai arm), §3.8,
§5 T2/T3; arch §5.

> **Checkpoint A:** the cheap path reads a `.fai` into `ContigInfo`s; types are usable; the oracle
> constants are committed. Pause for review.

### Milestone B — the FASTA streaming pass (the heart; silent-failure isolated)

**B1. The one-buffer, byte-zero pass — `Fasta { fai: None }`.**  ✅  **Own commit, do not bundle.**
The from-byte-zero streaming loop (one ~64 KiB buffer, never a whole contig): per contig produce
name, length, `offset`, `line_bases`, `line_width`, the per-contig MD5 and the running
whole-reference MD5, in one pass. The single predicate — a base is `[0x21,0x7E]` (`isgraph`); CR-LF
falls out of `line_width − line_bases`; uniform-line-width guard (T3). Copy
`compute_contig_md5_streaming`'s streaming shape but with the `isgraph` predicate (not skip-`\n\r`).
**Silent** (a wrong number is a wrong `.fai`/M5): lands with its oracle green. *Depends:* A1, A3.
*Source:* spec §4, §3.4, §3.8.

**B2. Oracle test for the pass.**  ✅  *(lands with B1 — it is B1's guard.)*
Against the committed constants: per-contig M5 = `samtools dict`; `LN` and the reconstructed `.fai`
fields = `samtools faidx`; the LF/CR-LF vector; per-contig MD5 vs `Md5::digest` one-shot; the
whole-reference digest vs the golden `.cat` header; the space/tab predicate edge (hashed-as-absent).
*Depends:* B1. *Source:* spec §6, §3.4.

**B3. The fasta-vs-`.fai` check — `Fasta { fai: Some }`.**  ✅  **Own commit, do not bundle.**
Add the field-for-field comparison of the pass's reconstruction against the supplied `.fai`;
disagreement → `FastaFaiMismatch` naming the field + contig. **Silent** (a wrong comparison is a
false pass/fail). Tests: a matching `.fai` passes; a **re-wrapped (stale)** `.fai` errors naming
`line_bases` — **mutation-verify** against a names-only check; a single-contig re-wrap (offset
unchanged); a reordering caught by the digest (T1). *Depends:* B1. *Source:* spec §3.3, §5 T1.

> **Checkpoint B:** the FASTA pass reconstructs the index byte-for-byte vs `samtools`, and a stale
> `.fai` cannot survive. Pause for review.

### Milestone C — the `.fai` writer

**C1. `write_fai`.**  ✅  **Own commit, do not bundle.**
Five-column `faidx.5` from `&[ContigInfo]` (`md5` ignored), written **atomically** (`.tmp` + rename).
**Silent** (a wrong byte is a wrong index). Tests: **byte-identical to the committed `samtools faidx`
`.fai`**; round-trip (`write_fai` → read back as `Fasta { Some }` verifies clean — a `read→write→read`
fixpoint); atomicity (a simulated mid-write failure leaves no partial file). *Depends:* B1 (for a
`ContigInfo` to write) and B3 (for the round-trip's verify). *Source:* spec §3.9.

> **Checkpoint C:** ng can index a reference itself, indistinguishably from `samtools faidx`. Pause
> for review.

### Milestone D — the single-flight cache

**D1. `ReferenceInfoCache` + `get_or_read`.**  ✅
Two-level `Mutex<HashMap<Key, Arc<Mutex<Option<Arc<ReferenceInfo>>>>>>`; key = `(source
discriminant, per-file (path, size, mtime))`; hold the map lock only for the slot lookup, the slot
lock across the read (the single-flight); **successes only** cached; mtime-unavailable → bypass;
`Send + Sync`. Tests: a hit computes **once** (read counter + `Arc::ptr_eq`); **single-flight** under
threads (barrier — the real read runs once, all threads get the same `Arc`); different keys do not
serialize; `Fai` / `Fasta{None}` / `Fasta{Some}` are distinct keys (no poisoning); an error is not
cached (retry re-runs); mtime-unavailable returns the same value uncached. *Depends:* A4 (a source
to read). *Source:* spec §3.7, §5 T8.

> **Checkpoint D:** N parallel readers of one reference pay one genome read; the failure and
> bypass paths are honest. Pause for review.

### Milestone E — the two entry points (compose everything)

**E1. `read_fai_verify_in_background`.**  ✅
`&Arc<ReferenceInfoCache>` in. Foreground `cache.get_or_read(Fai)` (returns now, `md5: None`);
spawn a **crossbeam** thread (not `rayon::spawn`) running `cache.get_or_read(Fasta { Some })`,
delivering its `Result` to the `VerificationHandle` (`#[must_use]`; `is_finished`; `join` returns
the **verified** info by **return**, not mutation; `Drop` without `join` warns to stderr). Tests:
info available before `join` and `join` upgrades it (`md5` None→Some); a stale `.fai` errors **on
`join`, not before** (mutation-verify a swallowed error fails the test); the verified result is
**cached** (later `get_or_read` = hit, read count 1); a concurrent `get_or_read(Fasta{Some})` **waits
on the slot** (single-flight-through-cache, same `Arc`); a failed verify does not poison the key;
abandoning the handle warns. *Depends:* B3, D1. *Source:* spec §3.10.

**E2. `read_reference_verifying_or_creating_fai`.**  ✅
`sibling_fai_path`, then branch: `.fai` present → E1 (return `(info, Some(handle))`); absent →
`cache.get_or_read(Fasta { None })` + `write_fai` (a write failure is **fatal** — `FaiWrite`), return
`(info, None)`. Tests: fai-present → background verify, **no write**, handle present; fai-absent →
scan + `.fai` written **byte-identical to `samtools faidx`**, no handle; a `.fai`-write failure
(read-only dir) is fatal, **and** the escape hatch works (`get_or_read(Fasta{None})` succeeds,
writes nothing). *Depends:* A2, C1, E1. *Source:* spec §3.11.

> **Checkpoint E:** the batteries-included entry point does the right thing on both branches, and
> the module is complete. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | unit tests: `Fai` parse vs committed `.fai`; dup-name / FASTQ-index / malformed rejections; `sibling_fai_path` no-I/O |
| B | **byte-parity vs `samtools dict`/`faidx`** (committed constants + LF/CR-LF vector); per-contig & reference MD5 vs one-shot/golden `.cat`; **stale-`.fai` errors** (mutation-verified) |
| C | **`write_fai` byte-identical to committed `samtools faidx` `.fai`**; `read→write→read` fixpoint; atomicity |
| D | read-count-once hit; **single-flight** under threads (same `Arc`); source-in-key; error-not-cached; mtime-bypass |
| E | info-before-`join` + `join`-upgrade; **stale-`.fai` errors on `join`** (mutation-verified); cache-populated + concurrent-reader-waits; orchestrator both branches + fatal-write + escape hatch |

## Out of scope (next plans)

- **`@SQ` ↔ `.fai` reconciliation** (the `ReadFilter::new` permutation hole) → **read ingestion's
  plan**; reuses `ContigList::first_disagreement` ([`fasta/mod.rs:69-100`](../../../../src/fasta/mod.rs)),
  compares *order* not just membership (arch §4).
- **Reader consolidation onto `ContigInfo`** (`RawChromReader`/`WindowedRefSeq` off `ContigList`,
  retiring `contig_list()`) → an ng refactor once this module lands (arch §4).
- **bgzip-compressed reference** and the **two-phase parallel hash** → error/deferred today; built
  when a consumer needs them (arch §4).
- **No `bench/`** — infra, single implementation, no bake-off (arch §7).
