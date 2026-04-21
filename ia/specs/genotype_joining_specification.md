# Genotype joining
The overall objective is to replicate GATK's GenotypeGVCFs that joins the variants found in different gVCFs per sample files into a VCF with all samples.

We want to join the variants found in an OverlappingVariantGroup.
The result of processing an OverlappingVariantGroup could be:

- An error if the variants can not be joined.
- A vector of Variants with one or more joined variants.

In general the phase should be maintained.
If there is per sample data with coverages these data should be also maintained in the final joined genotypes.

The output FORMAT fields should be: GT, GQ, PL, DP, AD.

# Nomenclature
We are going to create two new terms in order to clarify the code:

- vallele: an allele found in a variant. Example: In this variant there are three valleles, AC,A,AG
- gallele: an allele found in a genotype. Example: This sample is diploid and homozygous, so it has two galleles (0/0), but only one vallele.

# Examples of Genotype joining

## Example: simple SNP
gVCF sample1
pos REF ALT GT_sample1
1   A   T   1/1

gVCF sample2
pos REF ALT GT_sample2
1   A   C   0/1

joined VCF
pos REF ALT GT_sample1 GT_sample2
1   A   T,C 1/1        0/2

## Example: Mixed:Deletion + SNP
gVCF sample1
pos REF ALT GT_sample1
1   AC  A   1/1
2   C   .   0/0

gVCF sample2
pos REF ALT GT_sample2
1   A   C   0/1
2   C   .   0/0

joined VCF
pos REF ALT  GT_sample1 GT_sample2
1   AC  A,CC 1/1        0/2
2   C   .    0/0        0/0

## Example: Complex:Deletion + SNP in region covered by deletion
gVCF sample1
pos REF ALT GT_sample1
1   AC  A   1/1
2   C   .   0/0

gVCF sample2
pos REF ALT GT_sample2
1   A   .   0/0
2   C   G   1/1

By default we want to remove these kinds of OverlappingVariantGroup. This is not what GATK does. We will have the option of doing what GATK does, but only if the user requests to do it that way explicitly. The parameter should be remove_complex and by default should be true.

joined VCF GATK style
pos REF ALT  GT_sample1 GT_sample2
1   AC  A    1/1        0/0
2   C   *,G  1/1        2/2

# Algorithm

## 1. Group variants by position

Group the variants within the OverlappingVariantGroup by their start position, because each position should be processed almost independently.
Track which samples have active deletions from earlier positions (deletion state carries forward across positions).

## 2. Complex variant detection

If remove_complex: scan all positions. If the OverlappingVariantGroup covers more than one position and any non-first position has a real (non-*, non-ref) ALT allele, return an Error of kind complex_variant_removed.

The rationale: the `*` (spanning deletion) allele is expected in GATK output for deletions and is not inherently complex. The complexity arises when a real variant (e.g., a SNP) co-occurs at a position that is also covered by a deletion from another sample.

## 3. Process each position

### 3.a Build the unified vallele list

Collect all alleles from all samples at this position:
- Strip `<NON_REF>` from the output allele list (but retain the `<NON_REF>` PL values for use in step 3.b).
- For samples inside an active deletion (tracked from step 1), their vallele is `*`.
- Build the unified vallele list: [ref, alt1, alt2, ...].
  The order of the gVCFs given in the input should be respected when adding new alleles.
- At the first position of a multi-position group (where deletions start), if the ref alleles have different lengths, extend shorter alleles by appending trailing ref bases. For example, if REF=AC and a sample has A->C (SNP), the SNP allele becomes CC (extended to match the 2-base ref).
- Multi-allelic variants within a single sample's gVCF are supported: they just contribute more alleles to the unified set.

### 3.b Remap each sample's PLs to the unified vallele set

This is a critical step. Each sample's original PLs cover only the genotypes for its own allele set. The EM needs PLs for all genotypes of the unified allele set.

For each sample:
- Map the sample's original allele indices to unified allele indices.
- For genotypes involving alleles the sample didn't have, use the `<NON_REF>` likelihood values if available. If no `<NON_REF>` data exists, use a high PL penalty (e.g., the max PL from that sample or a fixed high value).
- Samples with no variant at this position (they are in a reference block or absent from the group): use synthetic hom-ref PLs, or remap from the reference block's PLs if the reference block carries PL data.
- Samples with GT-only (no PL field): generate synthetic PLs from the GT call (as already implemented: 99% confidence on the called genotype, remainder spread uniformly).

### 3.c Assign tentative galleles

Once we have the unified valleles, change the galleles for each sample accordingly (mapping original allele indices to unified indices). These galleles are tentative and will be used as fallback; the final galleles come from the EM posteriors.

### 3.d Run the EM algorithm

Run the EM algorithm on the remapped PLs to estimate allele frequencies and genotype posteriors. This uses the existing `estimate_posteriors()` function.

### 3.e Derive output fields from posteriors

- Normalize the posteriors to phred scale: PP = -10 * log10(posterior). Normalize so that the best (most likely genotype) = 0.
- GT = genotype with PP = 0 (highest posterior).
- GQ = second-lowest PP value (confidence in the GT call).
- PL = the remapped input likelihoods from step 3.b (this is what GATK outputs as PL; it preserves the original evidence).
- QUAL = phred-scaled probability that the site is non-variable (all samples hom-ref).

### 3.f Prune unused valleles

After the EM assigns final genotypes, some valleles might no longer be used by any sample (e.g., the EM reassigned all carriers of a rare allele to a different genotype). These unused valleles should be removed from the allele list, and all gallele indices should be renumbered accordingly. PL values should also be recomputed for the reduced allele set.

### 3.g Skip invariant positions

If the site is non-variable after pruning (no sample has a non-ref GT, i.e., all samples are hom-ref), do not emit a variant for this position.

### 3.h Update deletion tracking

Update the deletion state for each sample/haplotype so the next position knows which samples are inside an active deletion.

# Implementation Plan

The new joining algorithm will live in `src/genotype_joining.rs`. The old `genotype_merging.rs` stays untouched until the new module is fully wired in.

## Data structures

### DeletionState
Tracks, per global-sample-index and per haplotype, how many positions remain in an active deletion. A flat `Vec<i32>` of size `total_samples * ploidy`. A value > 0 means that haplotype is inside a deletion and its allele at the current position is `*`.

### PositionGroup
All input variants that share the same start position within an OverlappingVariantGroup, together with their source iterator indices. Fields: `pos: u32`, `variants: Vec<&Variant>`, `source_iter_idxs: Vec<usize>`.

### UnifiedAlleles
The merged allele list for a single position. Fields: `alleles: Vec<String>` (index 0 = ref), plus `per_sample_allele_maps: Vec<Vec<usize>>` — for each global sample, maps its original allele index to the unified allele index.

### PositionOutput
The output for one position after EM. Fields: `pos`, `alleles`, `genotypes`, `phase`, `pls`, `gqs`, `qual`.

## Function decomposition and TDD order

Each step below is a single implementation unit: write the test first (red), then implement (green). Functions marked *(trivial)* are small enough to not need their own dedicated test — they will be covered by the tests of the function that uses them.

### Step 1: `group_variants_by_position`

```
fn group_variants_by_position<'a>(
    group: &'a OverlappingVariantGroup,
) -> Vec<PositionGroup<'a>>
```

Takes an OverlappingVariantGroup and returns a list of PositionGroups sorted by position. Each PositionGroup gathers all variants that start at that position, together with their `source_var_iter_idx`.

Tests:
- Single-position group → one PositionGroup
- Multi-position group (e.g., deletion + SNPs) → multiple PositionGroups sorted by pos
- Variants from different iterators at the same position are grouped together

### Step 2: `detect_complex_variant`

```
fn detect_complex_variant(
    position_groups: &[PositionGroup],
) -> bool
```

Returns `true` if the group covers more than one position AND any non-first position has a real ALT allele (not `*`, not ref-only). Used to decide whether to return `ComplexVariantRemoved` error.

Tests:
- Single position → false
- Multi-position, all non-first are ref-only → false
- Multi-position, a non-first position has a real SNP alt → true
- Multi-position, a non-first position has only `*` → false

### Step 3: `build_unified_alleles`

```
fn build_unified_alleles(
    pos_group: &PositionGroup,
    iter_infos: &[VariantIteratorInfo],
    deletion_state: &DeletionState,
    ref_context: Option<&str>,  // trailing ref bases for extending short alleles at first position
) -> UnifiedAlleles
```

For one position, builds the unified allele list and per-sample allele index mapping.

- Reference allele comes from the longest ref at this position (or is extended with `ref_context`).
- Each sample's alleles are mapped to the unified set. New alleles are appended in input-order.
- Samples inside an active deletion get `*` as their allele.
- Samples not present at this position (no variant in any iterator) get a synthetic ref-only mapping.
- `<NON_REF>` is excluded from the unified allele list.

Tests:
- Simple SNP: two samples with different ALTs → unified list [ref, alt1, alt2], correct mappings
- Deletion at first position: ref alleles of different lengths → shorter alleles extended
- Sample inside active deletion → its allele is `*`
- Multi-allelic sample → all its alleles added to unified list
- Allele deduplication: two samples with same ALT → single entry, both map to it
- `<NON_REF>` stripped from output but accounted for in mapping

### Step 4: `remap_pls`

```
fn remap_pls(
    sample_pls: &[f64],
    sample_num_alleles: usize,   // includes <NON_REF> if present
    unified_num_alleles: usize,
    allele_mapping: &[usize],    // sample allele idx -> unified allele idx
    non_ref_index: Option<usize>, // index of <NON_REF> in sample's alleles, if present
    ploidy: usize,
) -> Vec<f64>
```

Maps one sample's PL values from its local genotype space to the unified genotype space.

- For genotypes that the sample had: copy the PL value directly.
- For genotypes involving alleles the sample didn't have: use the `<NON_REF>` PL value if available (by finding the corresponding genotype in the sample's PL array where the unknown allele maps to `<NON_REF>`). If no `<NON_REF>` data exists, use a high penalty PL.
- This is the most algorithmically complex function in the module.

Tests:
- Identity mapping (sample has same alleles as unified) → PLs unchanged
- Sample has 2 alleles, unified has 3 → new genotype slots filled from `<NON_REF>` PLs
- Sample has no PLs (empty) → return empty (handled upstream with synthetic PLs)
- Diploid with `<NON_REF>`: genotype involving unknown allele uses `<NON_REF>` likelihood

### Step 5: `build_position_pls`

```
fn build_position_pls(
    pos_group: &PositionGroup,
    unified: &UnifiedAlleles,
    iter_infos: &[VariantIteratorInfo],
    deletion_state: &DeletionState,
    ploidy: usize,
) -> Vec<f64>
```

Builds the full PL matrix for all samples at one position by calling `remap_pls` for each sample. Handles:
- Samples with real PLs → remap
- Samples with GT only (no PL) → generate synthetic PLs then remap
- Samples not present at this position → synthetic hom-ref PLs
- Samples inside a deletion → synthetic hom-`*` PLs

Tests:
- Two samples, one with PLs, one GT-only → combined flat PL vector of correct size
- Sample absent from position → gets hom-ref synthetic PLs
- Sample inside deletion → gets hom-`*` synthetic PLs

### Step 6: `posteriors_to_calls`

```
fn posteriors_to_calls(
    posteriors: &SitePosteriors,
    num_samples: usize,
    ploidy: usize,
    num_alleles: usize,
) -> (Vec<i8>, Vec<u16>, f64)  // (genotypes, gqs, qual)
```

Converts EM posteriors into final GT calls and GQ values.

- GT = genotype with highest posterior for each sample.
- GQ = phred-scaled difference between best and second-best posterior, capped at 99.
- QUAL = from `SitePosteriors.qual`.

Tests:
- Clear het (one genotype has posterior ~1.0) → correct GT, GQ ≈ 99
- Ambiguous call (two genotypes close) → lower GQ
- All hom-ref → QUAL near 0

### Step 7: `prune_unused_alleles`

```
fn prune_unused_alleles(
    alleles: &[String],
    genotypes: &[i8],
    pls: &[f64],
    ploidy: usize,
    num_samples: usize,
) -> (Vec<String>, Vec<i8>, Vec<f64>)
```

Removes alleles not referenced by any sample's GT. Renumbers allele indices in genotypes and recomputes PLs for the reduced allele set.

Tests:
- No unused alleles → unchanged
- One unused alt allele → removed, indices shifted, PL columns dropped
- Ref always kept even if all samples are alt

### Step 8: `update_deletion_state`

```
fn update_deletion_state(
    pos_group: &PositionGroup,
    unified: &UnifiedAlleles,
    genotypes: &[i8],
    deletion_state: &mut DeletionState,
    iter_infos: &[VariantIteratorInfo],
    ploidy: usize,
)
```

After processing a position, updates `DeletionState`. For each sample/haplotype: if the called genotype is a deletion allele (ref_len > alt_len), set remaining positions. Decrement existing counters. *(trivial — tested via the orchestrator integration test)*

### Step 9: `join_genotypes` (orchestrator)

```
pub fn join_genotypes(
    group: &OverlappingVariantGroup,
    iter_infos: &[VariantIteratorInfo],
    prior: &PriorConfig,
    remove_complex: bool,
) -> VcfResult<Vec<Variant>>
```

Top-level function. Orchestrates steps 1-8:

1. `group_variants_by_position`
2. If `remove_complex` and `detect_complex_variant` → return `ComplexVariantRemoved` error
3. Initialize `DeletionState`
4. For each position:
   a. `build_unified_alleles`
   b. `build_position_pls`
   c. Build tentative genotypes from `UnifiedAlleles` mappings *(inline, trivial)*
   d. `estimate_posteriors` (existing)
   e. `posteriors_to_calls`
   f. `prune_unused_alleles`
   g. Skip if invariant
   h. `update_deletion_state`
   i. Build output `Variant`
5. Return collected variants

Tests (integration-level within the module):
- Simple SNP joining (the spec example)
- Deletion + SNP (the spec example)
- Complex variant detection + removal
- Multi-position group producing multiple output variants
- Phase preservation through the full pipeline

## Implementation sequence summary

| Order | Function                    | Test? | Depends on        |
|-------|-----------------------------|-------|--------------------|
| 1     | `group_variants_by_position`| Yes   | —                  |
| 2     | `detect_complex_variant`    | Yes   | Step 1             |
| 3     | `build_unified_alleles`     | Yes   | —                  |
| 4     | `remap_pls`                 | Yes   | —                  |
| 5     | `build_position_pls`        | Yes   | Steps 3, 4         |
| 6     | `posteriors_to_calls`       | Yes   | —                  |
| 7     | `prune_unused_alleles`      | Yes   | —                  |
| 8     | `update_deletion_state`     | No    | Step 3             |
| 9     | `join_genotypes`            | Yes   | Steps 1-8          |

## Deferred work (not part of this plan)

- Adding DP, AD fields to `Variant` and parsing them from gVCF
- Writing GQ, PL, DP, AD in `VcfWriter` (currently GT-only)
- Wiring `join_genotypes` into the pipeline (replacing `genotype_merging`)
- Removing `genotype_merging.rs`