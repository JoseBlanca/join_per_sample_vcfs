# Code Review: psp_container_generalization

**Date:** 2026-06-15
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** the `.psp` container generalization (architecture §10), `git diff aa6a105 f0dfffc -- src/psp/` on branch `ssr-architecture`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** PR diff — the 5 commits `1eae1e9..f0dfffc` that generalize the `.psp` per-sample artefact format from SNP-only to a generic columnar container hosting two schemas (`snp`, `ssr`) via three traits (`PspKind` / `BlockAccumulator` / `BlockDecoder`).
- **Reviewed against:** `aa6a105` (pre-feature parent) → `f0dfffc` (tip), branch `ssr-architecture`.
- **In-scope files:**
  - [src/psp/kind.rs](../../../../src/psp/kind.rs) (new — the trait surface)
  - [src/psp/registry_ssr.rs](../../../../src/psp/registry_ssr.rs) (new — the SSR schema)
  - [src/psp/writer.rs](../../../../src/psp/writer.rs) (`PspWriter<W, S>`, `SnpKind`, `SnpBlock`, SSR writer)
  - [src/psp/reader.rs](../../../../src/psp/reader.rs) (`RecordsIter<'r, R, S>`, `SnpDecoder`, `read_and_inflate_column`)
  - [src/psp/registry.rs](../../../../src/psp/registry.rs) (`ColumnKey::tag/from_tag`, `columns_for_kind`)
  - [src/psp/header.rs](../../../../src/psp/header.rs) (`kind` tag)
  - [src/psp/block.rs](../../../../src/psp/block.rs) (dropped invariant)
  - [src/psp/errors.rs](../../../../src/psp/errors.rs) (`UnknownKind`)
  - [src/psp/index.rs](../../../../src/psp/index.rs) (doc-only)
  - [src/psp/mod.rs](../../../../src/psp/mod.rs) (visibility + re-exports)
  Direct callers/callees of changed items were in-scope surface.
- **Deliberately out of scope:** `src/ssr/*` (the SSR Stage-1 compute path, separately unit-tested); the cohort columnar reader path (`BlockColumnReader`/two-phase decode — deliberately left untouched, confirmed); docs commit `a0644f8` and all `doc/`+`ia/reports/` files; pre-existing untouched psp code (pre-existing Blockers would still be raised — none found).
- **Categories dispatched:** reliability (always), errors (always), naming (always), defaults (public API + the `kind`/serde defaults), idiomatic (always), refactor_safety (always — primary, byte-identity-critical refactor), module_structure (multi-file), smells (always), tooling (crate), extras (parser/stable-output/hot-path/diff-matches-intent). **Skipped:** `unsafe_concurrency` — the diff introduces no `unsafe`, `Arc`, `Mutex`, atomics, channels, `async`, or thread spawning (crate-level `#![forbid(unsafe_code)]` holds; the zstd decompressor is passed single-threaded).

## 2. Verdict

**Approve-with-changes.**

The SNP production path is sound. `refactor_safety` confirmed all five byte-identity / behaviour-preservation claims by reasoning over the diff (no Blockers, no Majors on the SNP path): the `ColumnDef.key` removal preserves the dispatch bijection for all twelve v1.0 tags, the interval-clamp rewrite is bit-for-bit equivalent for the SNP degenerate point, and the `SnpDecoder` extraction preserved every cursor/reset. Gates are green (fmt, clippy `-D warnings`, 1138 lib tests, doc).

Every Major below is in the **new `ssr` schema**, which has no production consumer yet (the Stage-1 driver is deferred). None can cause wrong results in production today because nothing writes `.ssr.psp` outside the round-trip test. They are must-fix-before-the-SSR-driver-ships items: a self-produced file the reader rejects (M1), a silent-truncation parser gap (M2), mislabeled corruption errors (M3), a silent cross-schema misdecode trap (M4), and a weakened compile-time guarantee the design explicitly leans on (M5).

## 3. Execution status

Commands run in the dev container (`./scripts/dev.sh cargo ...`), output quoted verbatim:

- `cargo fmt --check` — **exit 0**, clean (no diff emitted).
- `cargo clippy --all-targets --all-features -- -D warnings` — **exit 0**: `Checking pop_var_caller v0.1.0` … `Finished \`dev\` profile [unoptimized + debuginfo] target(s) in 3.56s` (0 warnings).
- `cargo test --lib psp` — `test result: ok. 200 passed; 0 failed; 0 ignored; 0 measured; 939 filtered out`.
- `cargo test --lib` (full) — `test result: ok. 1138 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 31.79s`.
- `cargo doc --no-deps` — **exit 0**: `Documenting pop_var_caller v0.1.0` … `Finished` (no warnings).

**Commands not run:**
- `cargo audit` — **unavailable in the container** (`error: no such command: \`audit\``; `cargo-audit` is not installed). Advisory-DB scan not performed; the diff adds no dependencies (`Cargo.toml` byte-identical across the range), so the dependency surface is unchanged.

**Findings labeled "Needs verification":** 0. (One sub-agent class of inflated line numbers for `registry.rs`/`registry_ssr.rs` was detected and corrected against the actual files during synthesis; all cited locations below were re-verified by direct read.)

## 4. Open questions and assumptions

1. **`last_pos` convention for interval schemas (gates M1).** Should `SsrBlock` store the *inclusive* last-touched position (`max(end - 1)`, matching the SNP point and the new index.rs doc) — or should the writer/reader bound be reconciled to accept exclusive `end == chrom.length + 1`? The current code does neither consistently, so a legal end-of-contig locus writes but won't reopen. Recommended: store inclusive (smaller surface, uniform field meaning).
2. **Missing-`kind` policy (gates the serde-default Minor).** The prompt states *no on-disk back-compat is required*, yet `#[serde(default = "default_kind")]` silently infers `"snp"` for a tag-less header. Is the back-compat default wanted for files written during steps 1a/1b (which predate the tag), or should a missing `kind` be a hard parse error now? If kept, should the inference be surfaced (a flag on `ParsedHeader`) rather than silent?
3. **`from_tag` exhaustiveness (judgment call #2, gates M5).** Adopt the `column_key!` declarative macro to restore the compile-time bijection the `ColumnDef.key` field used to give — or keep the hand-listed array + drift-guard test and *soften the doc comments* that currently claim a compile-time guarantee `from_tag` does not provide? Either resolves M5; the docs must not over-claim.
4. **`records_of::<S>` kind enforcement (judgment call resolution, gates M4).** Make it fallible (`try_records_of -> Result<_, KindMismatch>` comparing `S::KIND` to `header().kind`), or accept caller-discipline-by-doc? The header already holds the kind, so the guard is one comparison.
5. **Judgment call #4 (`#[allow(private_bounds)]`) — resolved, no change.** Three reviewers (idiomatic, module_structure, tooling) independently concluded the `#[allow(private_bounds)]` + unbounded `PspKind::Decoder` + re-bound-at-use-site pattern is the cleanest option that keeps the `pub(crate)` wire types out of the public API; a sealed trait or public marker supertrait would add surface for no caller benefit. Optional enhancement: `#[expect(private_bounds)]` (stable on the pinned 1.95 toolchain) auto-warns if `BlockDecoder` ever goes public.
6. **Judgment call #1 (two `SsrLocusRecord` types) — see Mi5.** The duplication (chrom-name-keyed Stage-1 record vs chrom_id-keyed container record) is defensible as a deliberate mirror; the only ask is to disambiguate the names or add a glossary note before the Stage-1 adapter makes both coexist.
7. **Judgment call #3 (dropped block-header invariant) — resolved.** Confirmed: SNP integrity is fully enforced elsewhere (writer `validate_record` rejects zero-allele records; reader B3 `Σ n_alleles == n_total_alleles` + per-allele column-count checks). The relaxation is correctly SNP→container-scoped. The only residue is the dead `AllelesLessThanRecords` variant (Mi6) and a narrowing of foreign-SNP-file strictness (acceptable, see M5/Findings).

## 5. Top 3 priorities

1. **M1 — the SSR writer can emit a file its own reader refuses.** A legal end-of-contig locus (`end == chrom.length + 1`, explicitly permitted by `validate_locus`) stores block-index `last_pos = chrom.length + 1`, which `PspReader::new` rejects as `BlockIndexPosOutOfRange`. Round-trip break on legal input; the existing test never reaches the contig end. Fix before any SSR data is produced.
2. **M2 — SSR decode silently truncates evidence on a malformed file.** No `sum(n_spanning) == n_total_alleles` cross-column check; an under-counting foreign/corrupt block decodes "successfully" with trailing per-locus profiles dropped, no error. Add the symmetric check the SNP path already has (B3).
3. **M4 — `records_of::<S>()` silently misdecodes on a schema mismatch.** No runtime `S::KIND == header.kind` check; the `S = SnpKind` default routes a future SSR caller through the SNP decoder by construction. One comparison turns a silent garbage-decode into a typed error.

## 6. Findings

### Blocker

None. The two Blocker-class behaviours (M1 round-trip break; M2 silent truncation) are filed Major because they live entirely in the unshipped SSR schema (no production consumer; only the round-trip test writes `.ssr.psp`), so they cannot produce wrong results in production today. Both are must-fix before the SSR Stage-1 driver lands.

### Major

#### M1: src/psp/registry_ssr.rs:379 — SSR writer emits a file its own reader rejects (exclusive `last_pos` vs reader's `last_pos <= chrom.length`)
- **Confidence:** High
- **Categories:** smells, errors, extras (convergent); confirmed directly by the orchestrator.
- **Problem:** `validate_locus` explicitly accepts `record.end == chrom.length + 1` (the half-open interval may sit one past the last 1-based position — comment at [writer.rs:584-587](../../../../src/psp/writer.rs#L584-L587), check at [writer.rs:596](../../../../src/psp/writer.rs#L596)). `SsrBlock::append` stores the exclusive end verbatim: `self.last_pos = self.last_pos.max(record.end)` ([registry_ssr.rs:379](../../../../src/psp/registry_ssr.rs#L379)). On read, `PspReader::new`'s block-index validation rejects any `last_pos > chrom.length` ([reader.rs:255-263](../../../../src/psp/reader.rs#L255-L263)). So a locus whose tract ends at the contig's final base (`end = length + 1`) is written successfully but fails to reopen with `BlockIndexPosOutOfRange`. The new index.rs doc compounds the confusion by calling `last_pos` the "maximum (**inclusive**) end coordinate" ([index.rs:26](../../../../src/psp/index.rs#L26)) while SSR stores an exclusive end — different semantics in one bare `u32` across the two kinds.
- **Why it matters:** A serializer that produces output its own deserializer rejects is data loss on legal input. The existing round-trip test (contig length 1,000,000, loci ≤ 276) never approaches the contig end, so the gate is green while the bug is live.
- **Suggested fix:** Store the inclusive last-touched position for SSR so the field's meaning is uniform across kinds and within the reader's `[1, chrom.length]` bound (`end > start` is guaranteed, so `end - 1 >= start >= 1`):
  ```rust
  // registry_ssr.rs, SsrBlock::append
  self.last_pos = self.last_pos.max(record.end - 1);
  ```
  Add `ssr_locus_at_contig_end_round_trips` (a locus with `end == chrom.length + 1`, write then `PspReader::new` + `records_of`). The reader's overlap arithmetic (`find_first_overlapping_block`) already treats `last_pos` as an inclusive max, so this also aligns it with the index.rs doc.

#### M2: src/psp/registry_ssr.rs:581 — SSR decode lacks `sum(n_spanning) == n_total_alleles`; an under-counting file silently drops trailing profiles
- **Confidence:** High
- **Categories:** reliability, extras (convergent).
- **Problem:** The SNP reader enforces B3 (`Σ n_alleles[i] == n_total_alleles`) on both the eager and two-phase paths and rejects interior/trailing zero-allele records. `SsrDecoder` has no equivalent. It decodes the two CSR columns sized by `n_profiles = header.n_total_alleles` and decodes `n_spanning` independently, but never reconciles `Σ n_spanning` against `n_profiles`. `next_record` ([registry_ssr.rs:581](../../../../src/psp/registry_ssr.rs#L581)) advances `next_profile` by `n_spanning[i]` per record and bounds-checks only the *over-run* case ([registry_ssr.rs:611](../../../../src/psp/registry_ssr.rs#L611)). On a foreign/corrupt file where `Σ n_spanning < n_profiles`, the trailing decoded profiles are simply never visited — a clean decode with silently dropped per-locus evidence (the swallowed-corruption class).
- **Why it matters:** This is a parser of an on-disk format that "may be malformed or foreign." A producer/transport desync between `n_spanning` and the CSR column count yields a sample whose evidence is truncated with no failure signal → wrong genotype likelihoods downstream.
- **Suggested fix:** After the offsets-agreement check in `decode_block`, add the symmetric reconciliation (and prefer a typed variant per M3 over the `Io` strings):
  ```rust
  let declared: u64 = self.n_spanning.iter().map(|&n| n).sum();
  if declared != header.n_total_alleles as u64 {
      return Err(PspReadError::SsrProfileCountMismatch {
          expected: header.n_total_alleles, got: declared,
      });
  }
  ```
  Add `decode_block_rejects_n_spanning_sum_below_profile_count` (synthetic decoded vectors: `n_total_alleles = 3`, `n_spanning = [1, 1]`).

#### M3: src/psp/registry_ssr.rs:564,611 — SSR structural-corruption failures mislabeled as `PspReadError::Io` with a synthesized `io::Error::other`
- **Confidence:** High
- **Categories:** errors, idiomatic, extras, smells, defaults, tooling (strongly convergent — 6 categories).
- **Problem:** Both SSR integrity checks report a malformed/foreign-file condition as an I/O failure:
  ```rust
  // registry_ssr.rs:564
  return Err(PspReadError::Io {
      context: "amb-lengths / amb-logliks CSR offsets disagree",
      source: std::io::Error::other("ssr profile column mismatch"),
  });
  // registry_ssr.rs:611 (and the M2 check, if added)
  return Some(Err(PspReadError::Io {
      context: "ssr profile index past CSR offsets",
      source: std::io::Error::other("n-spanning exceeds decoded profiles"),
  }));
  ```
  Neither is an I/O fault. The rest of `psp` surfaces corruption through dedicated variants (`BlockHeaderInvariant`, `ColumnTruncated`, `UncompressedLenMismatch`, `MissingRequiredColumnInManifest`, and the new `UnknownKind`). The synthesized `io::Error::other("…")` carries no real `source()` chain — it is a freshly-minted string smuggled through the `Io` variant. This is the mechanism-named-error anti-pattern and the stringly-typed-error smell.
- **Why it matters:** A caller matching on `PspReadError` cannot distinguish a torn disk read (retryable) from deterministic structural corruption (fatal); a retry-on-`Io` loop would spin forever on a corrupt file. Malformed-input tests can only assert `Io`, not the specific fault.
- **Suggested fix:** Add typed variants and route both sites (plus M2) through them:
  ```rust
  /// A schema-specific per-block structural invariant was violated on
  /// decode (corrupt/foreign block, not an I/O fault).
  #[error("ssr block structural invariant violated: {context}")]
  SsrBlockStructureInvalid { context: &'static str },
  ```
  (or finer-grained `SsrProfileOffsetsDisagree` / `SsrProfileCountMismatch { expected, got }`). Keep `Io` for real `read_exact`/source failures. Note: the pre-existing SNP `materialise_next_record` "internal invariant" `Io` ([reader.rs:772](../../../../src/psp/reader.rs)) is the same anti-pattern in untouched code — see Out-of-scope.

#### M4: src/psp/reader.rs:359 — `records_of::<S>()` does not check `S::KIND == header.kind`; cross-schema misdecode is silent and the `SnpKind` default invites it
- **Confidence:** High
- **Categories:** defaults, reliability (convergent).
- **Problem:** `records_of` ([reader.rs:359](../../../../src/psp/reader.rs#L359)) just calls `RecordsIter::new`; the decoder is chosen purely from `S::Decoder` with no consultation of the `kind` the reader already parsed into `header().kind`. Two hazards compound: (1) `RecordsIter<'r, R, S: PspKind = SnpKind>` means the obvious `records()` is silently `S = SnpKind`, so a future SSR caller reaching for it gets the SNP decoder against an SSR file — wrong path by construction, no compile or runtime error; (2) even an explicit `records_of::<SnpKind>()` on an SSR file (or vice-versa) is accepted and decodes the wrong columns. The mismatch is caught *incidentally* today only because the B1 required-column check rejects a manifest lacking the other schema's tags — that is disjoint-tag luck, not a guard, and there is no test for it.
- **Why it matters:** The entire point of the `kind` tag + `UnknownKind` refusal is to prevent decoding columns under the wrong schema. `records_of` re-opens that hole on the typed read path.
- **Suggested fix:** Make it fallible and compare against the parsed kind:
  ```rust
  pub fn records_of<S: PspKind>(&mut self) -> Result<RecordsIter<'_, R, S>, PspReadError>
  where S::Decoder: BlockDecoder<Record = S::Record> {
      if S::KIND != self.header().kind {
          return Err(PspReadError::KindMismatch { expected: S::KIND, found: self.header().kind.clone() });
      }
      Ok(RecordsIter::new(self, RangeClamp::None))
  }
  ```
  Add `records_of_rejects_wrong_kind`. (Also generalize `region_records` over `S` with the same guard — see Out-of-scope: SSR region queries are currently unreachable.)

#### M5: src/psp/registry.rs:229 — `from_tag`'s hand-listed array weakens the M4 compile-time exhaustiveness the design claims (and the docs over-state it)
- **Confidence:** High
- **Categories:** idiomatic, smells (convergent). Directly answers judgment-call #2.
- **Problem:** Removing `ColumnDef.key` split the key↔tag pairing into `tag()` (exhaustive `match Self`, [registry.rs:209](../../../../src/psp/registry.rs#L209)) plus `from_tag(tag)` ([registry.rs:229](../../../../src/psp/registry.rs#L229)), a linear scan over a **hand-written array literal** of the variants (mirrored by `SsrColumnKey::from_tag`, [registry_ssr.rs:78](../../../../src/psp/registry_ssr.rs#L78)). Adding a `ColumnKey` variant forces an arm in `tag()` and in every dispatch `match key` — but **not** an entry in the `from_tag` array. A variant omitted from the array compiles green; `from_tag` returns `None` for its tag and the column silently no-ops at encode/decode. The only backstop is the runtime test (`column_keys_are_unique_and_cover_v1_0` / `ssr_columns_are_well_formed`). Yet the doc comments on `ColumnKey`/`from_tag` claim the M4 "compile-time hole" property the original `key` field actually provided — the tag→key direction is now a runtime-test guarantee, not a compile-time one. The prompt asked whether the M4 guarantee still holds: **for tag→key, no.**
- **Why it matters:** A contributor adding a column, building green, and not running the psp suite ships a column whose codec silently no-ops — the "wrong results without panicking" class, for the very mechanism (M4) the design cites as compile-safe.
- **Suggested fix:** Restore the compile-time bijection with a declarative macro that generates the enum + both directions from one variant↔tag list (justified per the macro rule — it generates *variants*, which a generic fn cannot):
  ```rust
  macro_rules! column_key {
      ($name:ident { $($variant:ident = $tag:expr),+ $(,)? }) => {
          #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
          pub enum $name { $($variant),+ }
          impl $name {
              pub const fn tag(self) -> u16 { match self { $(Self::$variant => $tag),+ } }
              pub fn from_tag(tag: u16) -> Option<Self> { match tag { $($tag => Some(Self::$variant),)+ _ => None } }
          }
      };
  }
  ```
  This also collapses the registry.rs/registry_ssr.rs duplication. If deferred, at minimum rewrite both `from_tag`s as `match tag { … _ => None }` (no array, no scan) **and** soften the doc claims to state the tag→key direction is checked by the named test, not the compiler.

### Minor

#### Mi1: src/psp/block.rs:821 — `BlockHeader.n_total_alleles` is a SNP-only name on a now-generic container field
- **Confidence:** High · **Categories:** naming (+ refactor_safety, smells touch the same field via M1).
- **Problem:** The write/read trait method was deliberately renamed to the schema-neutral `n_entries()` ([kind.rs:117](../../../../src/psp/kind.rs#L117)), but the value it carries is still the field `n_total_alleles` on the generic `BlockHeader` ([block.rs:821](../../../../src/psp/block.rs#L821)), shuttled under that name through `assemble_block_header::<S>` and read back by the SSR decoder as `n_profiles`. `validate_block_header_invariants` even adds a prose NOTE ([block.rs:937](../../../../src/psp/block.rs#L937)) explaining the field is *not* allele-counted for SSR — a comment patching a misleading name.
- **Why it matters:** Re-introduces the SNP assumption the abstraction removed; a reader of the SSR path meets a field named "alleles" holding profile counts.
- **Suggested fix:** Rename the in-memory wire-struct field to `n_entries` (the field is varint-positional on disk, not name-keyed, so it is a pure rename): update `block.rs`, `writer.rs` flush, `reader.rs` SNP decode, `registry_ssr.rs:489`. If judged too broad for this PR, add a documented alias note on the field.

#### Mi2: src/psp/reader.rs:905 — eager block-buffer drop on exhaust was lost in the `SnpDecoder` extraction (memory-residency regression)
- **Confidence:** High · **Categories:** refactor_safety.
- **Problem:** The old `Iterator::next` freed an exhausted block's columns immediately (`self.cur_block = None`); the new code only flips `self.cur_block_loaded = false` ([reader.rs:905](../../../../src/psp/reader.rs#L905)) and never clears `SnpDecoder::cur_block` until the next `decode_block` overwrites it. `DecodedBlock`'s per-record/per-allele fixed-width `Vec`s are read by index/copy (not `mem::take`), so the last-yielded block stays fully resident — transiently between blocks, and permanently for the final block while the iterator lives. Correctness is unaffected (`next_record` is gated on `cur_block_loaded`). The `materialise_next_record` doc still claims `mem::take` ownership-move + "dropped wholesale once exhausted," now both false.
- **Why it matters:** The project's headline thesis is RAM-for-sample-count scaling; this is a silent behaviour change vs the documented contract that the next memory-measurement run would surprise on. Bounded to one block's fixed-width columns (the heavy CSR slabs are reused either way), hence Minor.
- **Suggested fix:** Add `BlockDecoder::unload()` (SnpDecoder sets `cur_block = None`; SsrDecoder clears its column vecs) and call it from `Iterator::next` when `next_record()` returns `None`. Update the stale `materialise_next_record` doc.

#### Mi3: src/psp/header.rs:99 — missing `kind` silently defaults to `"snp"` under a stated no-back-compat constraint, invisibly
- **Confidence:** High · **Categories:** defaults. See Open Question 2.
- **Problem:** `#[serde(default = "default_kind")]` ([header.rs:99](../../../../src/psp/header.rs#L99)) makes a tag-less header deserialize to `"snp"` and sail through with no diagnostic — an invisible inferred default right beside the `UnknownKind` refusal that exists to close exactly this misdecode hole. Under the prompt's *no on-disk back-compat* constraint, the back-compat justification is moot.
- **Why it matters:** A hand-edited/truncated header (or an SSR-content file that lost its tag) is silently decoded as SNP rather than rejected.
- **Suggested fix:** Either drop the default and make `kind` mandatory (a missing field becomes a hard parse error; invert the `missing_kind_defaults_to_snp` test), or keep it but surface the inference (e.g. `ParsedHeader::kind_was_defaulted: bool`). Resolve via Open Question 2.

#### Mi4: src/psp/registry.rs:528 — `columns_for_kind` puts the cross-kind dispatch in the SNP registry, coupling it to `registry_ssr`
- **Confidence:** High · **Categories:** module_structure.
- **Problem:** `columns_for_kind` ([registry.rs:528](../../../../src/psp/registry.rs#L528)) lives in the SNP registry but reaches into `super::registry_ssr::{SSR_KIND, SSR_COLUMNS}`, while `registry_ssr` imports `ColumnDef`/`Cardinality`/… back from `registry` — a mutual `use` edge, with `registry.rs` accreting "the SNP schema *and* the global kind table." `SNP_KIND`/`KNOWN_KINDS`/`columns_for_kind` are conceptually a peer of *both* schema modules.
- **Why it matters:** Obscures the dependency graph; every future kind adds an arm in the SNP module; `registry_ssr` can't be understood without `registry`.
- **Suggested fix:** Move `SNP_KIND` (or a neutral kind table), `KNOWN_KINDS`, and `columns_for_kind` into `kind.rs` and dispatch via `SnpKind::columns()`/`SsrKind::columns()` through the trait both implement, so neither registry imports the other. Updates the `header.rs` call site to `kind::columns_for_kind`.

#### Mi5: src/psp/mod.rs:35 — `pub mod registry_ssr` over-exposes a consumer-less `#![allow(dead_code)]` schema module
- **Confidence:** High · **Categories:** module_structure.
- **Problem:** `registry_ssr` is `pub mod` ([mod.rs:35](../../../../src/psp/mod.rs#L35)) while its sibling `registry` is `pub(crate) mod`; the module carries `#![allow(dead_code)]` and (grep-confirmed) has no consumer outside `src/psp/`. The `pub` publishes `SsrKind`/`SsrBlock`/`SsrLocusRecord`/`SSR_COLUMNS`/`SsrColumnKey` to crate rustdoc before any caller exists; the SNP equivalents are reached through `pub writer`/`pub reader`, not a `pub` registry.
- **Why it matters:** Freezes a pre-alpha surface (with `pub(crate)`-typed signatures) for no consumer, and the asymmetry with `registry` misleads on which modules are entry points.
- **Suggested fix:** `pub(crate) mod registry_ssr;` to match `registry`; re-export the specific SSR types through `mod.rs` (mirroring `pub use registry::ColumnDef`) when an external API is actually wanted.

#### Mi6: src/psp/errors.rs:92 — dead `BlockHeaderInvariantKind::AllelesLessThanRecords`; retained-for-compat rationale lives only in block.rs
- **Confidence:** High · **Categories:** errors, tooling, naming, refactor_safety (convergent).
- **Problem:** The §10 change dropped the only construction site; grep confirms only the definition ([errors.rs:92](../../../../src/psp/errors.rs#L92)) + a block.rs comment remain. Clippy can't flag it (`pub` + re-exported variant). The "retained for compatibility but no longer raised" note sits at the removed call site ([block.rs:937](../../../../src/psp/block.rs#L937)), invisible at the declaration.
- **Why it matters:** A dead public error variant implies a guarantee (`n_total_alleles >= n_records`) the container no longer makes — the opposite of what it once signalled; the next maintainer may delete it (format-compat risk) or rewire it (re-SNP-ifying the generic validator).
- **Suggested fix:** Either remove the variant (`#[non_exhaustive]` already makes that non-breaking for external matchers) or co-locate a doc note on the variant stating it is intentionally unraised and where SNP integrity now lives.

#### Mi7: src/psp/writer.rs:596 — `PosOutOfRange` reused for the SSR `end` field reports a wrong range and conflates two failures
- **Confidence:** High · **Categories:** errors.
- **Problem:** `validate_locus` rejects a bad `end` with `PspWriteError::PosOutOfRange { pos: record.end, chrom_length }` ([writer.rs:596](../../../../src/psp/writer.rs#L596)), but `end` is legally `[start+1, chrom.length + 1]` while the variant's `Display` says `out of [1, {chrom_length}]` — off by one for the legal `length + 1`, and the `end <= start` (zero/negative span) vs `end` past-contig failures are conflated into one positional variant.
- **Why it matters:** An operator reading the log is told a false range and can't tell which rule fired.
- **Suggested fix:** Add an SSR-meaningful variant (`LocusEndOutOfRange { start, end, chrom_length }` with an honest `({start}, {chrom_length}+1]` message) and split the `end <= start` case into its own arm.

#### Mi8: src/psp/registry_ssr.rs:217 — `amb-logliks` `finite_constraint: false` + no SSR finite sweep lets NaN logliks round-trip silently
- **Confidence:** Medium · **Categories:** reliability, extras (convergent).
- **Problem:** The finite-float sweep lives only on the SNP scalar path ([reader.rs:1369](../../../../src/psp/reader.rs), matching `DecodedColumn::AlleleQSumLog`); the SSR decoder calls `decode_list_column_csr::<f32>` directly and never consults `finite_constraint`, and `append` pushes loglik with no `is_finite` check. `amb-logliks` is documented as renormalized log-probabilities feeding a downstream EM; `-inf` is legitimate (`log(0)`), but `NaN` is never legitimate and currently round-trips undetected. Even setting `finite_constraint: true` would not help today — the sweep does not cover CSR-shaped f32 lists.
- **Why it matters:** A NaN loglik surviving to the deferred Stage-1 caller poisons the EM sum with no decode-time signal.
- **Suggested fix:** Decide explicitly: add a writer-side NaN-only reject in `SsrBlock::append` (permit `-inf`), or document on the `finite_constraint: false` line why NaN-rejection is the consumer's responsibility. If finiteness is required, extend the generic sweep to list f32 columns. Add `append_rejects_nan_loglik` if guarding.

#### Mi9: src/psp/registry_ssr.rs:603 — `span[i] as u32` truncates silently, inconsistent with the adjacent `delta-start` `u32::try_from`
- **Confidence:** High · **Categories:** extras, idiomatic.
- **Problem:** `span` is decoded as `Vec<u64>` (varint); `let end = start.saturating_add(self.span[i] as u32)` does an unchecked narrowing cast. The sibling `delta-start` cast a few lines up is correctly guarded with `u32::try_from` → `ColumnElementDecode { VarintOverflow }`. On a hand-built/corrupt block a span ≥ 2³² wraps to a small value, yielding a wrong `end` with no error (the writer's `validate_locus` prevents it for well-formed files).
- **Why it matters:** Same untrusted-input class as the guard right beside it; cheap to make symmetric.
- **Suggested fix:** `let span = u32::try_from(self.span[i]).map_err(|_| PspReadError::ColumnElementDecode { column: "span".into(), entry: i, source: ScalarDecodeError::VarintOverflow })?;` then `start.saturating_add(span)`.

#### Mi10: src/psp/kind.rs:9 & src/psp/mod.rs:1 — stale module docs after the generalization
- **Confidence:** High · **Categories:** module_structure.
- **Problem:** `kind.rs` module doc still says "this module is the **writer-side surface only**" ([kind.rs:9-11](../../../../src/psp/kind.rs#L9-L11)) though step 1b added the read-side `BlockDecoder` trait in the same file. `mod.rs` still frames `.psp` as the SNP-only "per-sample pileup byte format" with `registry` as "the single source of truth" ([mod.rs:1-22](../../../../src/psp/mod.rs#L1-L22)), omitting the new `kind`/`registry_ssr` modules and the second registry.
- **Why it matters:** Top-level module docs are the entry point; both now misdescribe the format's scope and structure.
- **Suggested fix:** Update `kind.rs` to mention the read-side `BlockDecoder` mirror; reframe `mod.rs` as a "generic `.psp` columnar container hosting the `snp` + `ssr` schemas," qualify the `registry` line as SNP-specific, and add `kind` + `registry_ssr` to the layout list.

#### Mi11: src/psp/registry_ssr.rs:39 — `INITIAL_LOCI/PROFILES/PAIRS_HINT` magic capacities decoupled from the block-size knob
- **Confidence:** High · **Categories:** defaults.
- **Problem:** The three hints (256 / 4096 / 8192, [registry_ssr.rs:39-41](../../../../src/psp/registry_ssr.rs#L39-L41)) are fixed literals, unlike the SNP hints which are *derived from* `TARGET_BLOCK_BYTES` precisely so retuning the block size retunes the hint. If a user dials block target/window up, the SSR block under-reserves and reallocs on the write path.
- **Why it matters:** Perf-only and no consumer yet, but the values set SSR write-path allocation behaviour the moment the driver lands, and they don't track the knob the SNP path deliberately tracks.
- **Suggested fix:** Either derive them from the block target like the SNP hints, or document the per-value basis (locus density → profiles → pairs) so a future tuner can re-derive, and reference them by name in `SsrBlock::new_block`.

### Nits

Grouped (per the rubric, not individually blocking):

- **`record_interval` returns a bare `(u32, u32, u32)`** ([kind.rs:76](../../../../src/psp/kind.rs#L76), both impls) — a transposed `start`/`end` would typecheck silently. A `struct RecordInterval { chrom_id, start, end }` (end exclusive, documented once) names the fields at the single destructure site. Low blast radius (one consumer); idiomatic + smells agree it is acceptable as-is but cleaner named.
- **`KNOWN_KINDS = "snp, ssr"`** ([registry.rs:523](../../../../src/psp/registry.rs#L523)) is a hand-maintained second source of truth vs `columns_for_kind`; derive both from one `&[(&str, &[ColumnDef])]` table (diagnostic-text drift only).
- **`projected_bytes` magic arithmetic** `2 + 1 + 4 * 5 + …` ([registry_ssr.rs:384](../../../../src/psp/registry_ssr.rs#L384)) — name the constituents like the SNP sibling; `4 * 5` silently couples to the per-locus u32-column count.
- **`spanning: Vec<Vec<(u16, f32)>>`** ([registry_ssr.rs:244](../../../../src/psp/registry_ssr.rs#L244)) — nested-Vec + tuple obsession; a `struct LadderEntry { length: u16, loglik: f32 }` + `type Profile = Vec<LadderEntry>` names the shape (in-memory only; wire layout unaffected). Deferrable until SSR has consumers.
- **Four `from_tag(...).expect("…")` sites** ([writer.rs](../../../../src/psp/writer.rs), [reader.rs:1595](../../../../src/psp/reader.rs#L1595), [registry_ssr.rs:271,512](../../../../src/psp/registry_ssr.rs)) are genuinely unreachable (verified — the tag is matched against the registry before `from_tag` runs) but lack the canonical `// UNREACHABLE:` / `// PANIC-FREE:` marker the convention asks for.
- **`#![allow(dead_code)]`** ([registry_ssr.rs:5](../../../../src/psp/registry_ssr.rs#L5)) — on the pinned 1.95 toolchain, `#[expect(dead_code, reason = "…")]` auto-warns when the driver makes items live, converting the "remove once the driver lands" intent into a compiler tripwire. Matches existing house style, so optional. (The inner `#![allow]` also sits *above* the `//!` module doc — conventionally `//!` comes first.)
- **`#[allow(private_bounds)]`** (three sites, reader.rs) — correct tool, tightly scoped; `#[expect(private_bounds)]` would auto-warn if `BlockDecoder` ever goes public.
- **SSR `decode_block` unknown-column skip** allocates a throwaway `Vec` via `read_compressed_blob` + `let _ =` ([registry_ssr.rs:497](../../../../src/psp/registry_ssr.rs#L497)) instead of reusing scratch / seeking; latent until optional columns land.
- **Stale `///` in the SSR test** referencing `key_ssr` (the concept is `SsrColumnKey::from_tag`).
- **`frr` / `amb` abbreviations** in `SsrColumnKey` variants and fields ([registry_ssr.rs:46](../../../../src/psp/registry_ssr.rs#L46)) have meaning only in the `ColumnDef::description`; add a one-line glossary to the module doc (the no-abbreviation rule requires it for non-universal terms).

## 7. Out of scope observations

- **`region_records` is not generalized over `S`** ([reader.rs:398](../../../../src/psp/reader.rs#L398) returns `RecordsIter<'_, R>` = SNP-only). The §10.5 interval-clamp was generalized but its only public window-query entry point was not, so SSR region/window queries are unreachable through the public reader even though `RecordsIter` + the clamp are generic. Follow-up: add `region_records_of::<S>` alongside the M4 kind-guard (same PR).
- **Pre-existing `reader.rs` "internal invariant" `Io` error** (`materialise_next_record without a loaded block`) is the same mechanism-named/synthetic-source anti-pattern as M3, in untouched code. Fold into the M3 typed-variant cleanup.
- **No golden `.psp` fixture exists** (grep for `include_bytes`/golden: none). The SNP byte-identity claim rests on the round-trip tests + reasoning, not a committed regression artefact. If a hard guarantee is wanted, capture a golden at the step-2 baseline and diff (separate task; consistent with the project's out-of-tree byte-identity practice).
- **The two `SsrLocusRecord` types** (chrom-name-keyed `src/ssr/pileup/locus_record.rs` vs chrom_id-keyed `src/psp/registry_ssr.rs`) — judgment-call #1. Defensible as a deliberate mirror; if kept, rename the container form (`SsrLocusColumns` / `PspSsrLocus`) or add a glossary note before the Stage-1 adapter makes both coexist unqualified.

## 8. Missing tests to add now

Grouped by the function under test (the reliability challenge-tests pass + extras feed this section):

- **`SsrBlock` / `PspReader` round-trip — `ssr_locus_at_contig_end_round_trips`** *(catches M1)*: write a locus with `end == chrom.length + 1`, then `PspReader::new` + `records_of::<SsrKind>()`. Currently the reader rejects the writer's own output with `BlockIndexPosOutOfRange`.
- **`SsrDecoder::decode_block` — `decode_block_rejects_n_spanning_sum_below_profile_count`** *(catches M2)*: synthetic decoded vectors `n_total_alleles = 3`, `n_spanning = [1, 1]`; assert a structural error rather than two records silently materialised.
- **`SsrDecoder::decode_block` — `decode_block_returns_error_on_csr_offset_disagreement`** *(locks the registry_ssr.rs:564 guard)*: encode `amb-lengths` rows `[[1,2],[3]]` and `amb-logliks` rows `[[0.0],[0.0,0.0]]` (same total, different boundaries), `n_total_alleles = 2`; assert the typed structural variant (per M3).
- **`RecordsIter<SsrKind>` region clamp — `ssr_region_query_keeps_locus_straddling_left_edge`** *(catches an off-by-one in `interval_before_window` on intervals)*: write `[40,60)` and a strictly-left locus on chr0, query `[50, 80]`; expect `[40,60)` kept (end 60 > 50), strictly-left dropped. Requires the generic region entry point (Out-of-scope item).
- **`records_of` — `records_of_rejects_wrong_kind`** *(catches M4)*: open a `snp` file, call `records_of::<SsrKind>()`, assert `KindMismatch` (not an incidental B1 failure or garbage records).
- **`validate_locus` — `validate_locus_rejects_end_equal_to_start` + `validate_locus_rejects_out_of_order_locus`**: `end == start` would underflow `span = end - start` in `append`; out-of-order would corrupt delta-start. Pin both writer guards.
- **`SsrBlock`/`SsrDecoder` — `ssr_round_trips_under_proptest`**: random `Vec<SsrLocusRecord>` (ascending starts per chrom, profiles incl. empty / single-pair / many-pair / `u16::MAX` lengths), a *small* byte cap to exercise the `projected_bytes` auto-flush, write→read equality (compare f32 by `to_bits`, or exclude NaN from the generator). The shared-offsets CSR pairing is the one structurally-new wire shape and is covered only by five hand-picked fixtures.
- **`SsrBlock::append` — `append_rejects_nan_loglik`** *(only if Mi8 adopts the writer guard)*.

## 9. What's good

- **The `BlockAccumulator` / `BlockDecoder` write/read mirror** ([kind.rs](../../../../src/psp/kind.rs)) cleanly separates schema-specific column buffers from the shared block-framing + zstd machinery, so the SSR schema is genuinely "a table + a record mapping" with no new wire code.
- **`read_and_inflate_column`** ([reader.rs:1700+](../../../../src/psp/reader.rs)) extracts the budget-check + read + inflate + length-verify as a shared scratch-reusing helper, letting `SsrDecoder` ride the exact decompression path `SnpDecoder` uses.
- **The shared CSR codec constructs offsets from decoded per-row counts** (offsets are never trusted from the wire; monotonic + terminal-at-`data.len()` by construction with a per-row DoS guard) — this is what makes the SSR within-column slice indexing panic-safe and is the right design for an untrusted parser.
- **`#[allow(private_bounds)]` + unbounded `PspKind::Decoder`** keeps the `pub(crate)` wire types out of the public API rather than widening it to satisfy the lint — three reviewers independently judged it the cleanest available option.
- **The `kind`-tag + `UnknownKind` hard refusal** ([header.rs](../../../../src/psp/header.rs), [errors.rs:294](../../../../src/psp/errors.rs)) is the right shape for selecting the column registry on read (let down only by the `records_of` gap in M4 and the silent default in Mi3).

## 10. Commands to re-verify

Run in the dev container (`./scripts/dev.sh cargo ...`):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib` (1138 expected; add the §8 tests, which should raise the count and — for M1/M2/M4 — fail until fixed)
- `cargo doc --no-deps`
- `cargo audit` — **not run** (cargo-audit absent in the container); install or run on the host if a dependency advisory scan is wanted (the diff adds no dependencies).

### Author response convention

Address each finding by identifier (`M1`…`M5`, `Mi1`…`Mi11`) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer Open Questions 1–7 first — several findings (M1, M3, M4, M5, Mi3) are gated on those decisions.

---

*Per-category audit trail: `tmp/review_2026-06-15_psp-container/`.*
