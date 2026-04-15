# Genotype Joining Implementation Plan

## Goal

Implement the GATK-style genotype joining algorithm as described in `specs/genotype_joining_specification.md`. This produces a new module `genotype_joining.rs` (and test file `genotype_joining_test.rs`) that processes `OverlappingVarGroup`s by splitting them into per-position slices, building unified allele sets, remapping PLs, running the EM algorithm, and emitting one output `Variant` per variable position.

The existing `genotype_merging.rs` stays untouched. Wiring the new module into `pipeline.rs` and `main.rs` is deferred to a later step.

## Key differences from the current merging approach

| Aspect | Current merging | New joining |
|---|---|---|
| Output granularity | One merged `Variant` per `OverlappingVarGroup` | One `Variant` per *position* within the group |
| Deletion handling | Haplotype-level string concatenation | Per-position `*` allele + `DeletionState` tracker |
| PL data | Ignored; synthetic PLs from GT | Parsed from input when available; synthetic as fallback |
| Complex variants | Not detected | Detected and optionally removed (`remove_complex`) |
| Allele pruning | Not performed | Unused alleles pruned after EM |
| Output FORMAT fields | GT only | GT, GQ, PL (DP, AD deferred) |

## Data structures

### DeletionState

Tracks which sample-haplotypes are inside an active deletion. Flat `Vec<i32>` of size `total_samples * ploidy`. A value > 0 means that haplotype's allele at the current position is `*`. Decremented after each position. Set when a deletion genotype is called.

```rust
struct DeletionState {
    // remaining[sample_idx * ploidy + haplotype] = positions left in deletion
    remaining: Vec<i32>,
    ploidy: usize,
}

impl DeletionState {
    fn new(total_samples: usize, ploidy: usize) -> Self { ... }
    fn is_in_deletion(&self, sample_idx: usize, haplotype: usize) -> bool { ... }
    fn any_haplotype_in_deletion(&self, sample_idx: usize) -> bool { ... }
    fn decrement_all(&mut self) { ... }
    fn set_deletion(&mut self, sample_idx: usize, haplotype: usize, remaining_positions: i32) { ... }
}
```

### PositionGroup

All input variants sharing the same start position within a group, with their source iterator indices.

```rust
struct PositionGroup<'a> {
    pos: u32,
    variants: Vec<&'a Variant>,
    source_iter_idxs: Vec<usize>,
}
```

### UnifiedAlleles

The merged allele list for one position, plus the mapping from each sample's original allele indices to unified indices.

```rust
struct UnifiedAlleles {
    // alleles[0] = ref, alleles[1..] = alts. No <NON_REF>.
    alleles: Vec<String>,
    // per_sample_maps[global_sample_idx][original_allele_idx] = unified_allele_idx
    // For samples not present at this position: empty vec.
    per_sample_maps: Vec<Vec<usize>>,
    // For each sample, the index of <NON_REF> in its original allele list (if present).
    per_sample_non_ref_idx: Vec<Option<usize>>,
}
```

### JoinedVariant (intermediate, not exported)

Collects the per-position output before it is converted to a `Variant`.

```rust
struct JoinedVariant {
    pos: u32,
    alleles: Vec<String>,
    genotypes: Vec<i8>,    // flat: [sample * ploidy .. (sample+1) * ploidy]
    phases: Vec<bool>,     // per sample
    gqs: Vec<u16>,         // per sample
    pls: Vec<f64>,         // flat: [sample * num_genotypes .. (sample+1) * num_genotypes]
    qual: f64,
}
```

## Function decomposition, pseudocode, and TDD order

Each step is one implementation unit: write the test(s) first, then implement.

---

### Step 1: `group_variants_by_position`

```rust
fn group_variants_by_position<'a>(
    group: &'a OverlappingVarGroup,
) -> Vec<PositionGroup<'a>>
```

**Pseudocode:**
```
fn group_variants_by_position(group):
    pos_map: BTreeMap<u32, PositionGroup> = empty

    for (variant, source_idx) in zip(group.variants, group.source_var_iter_idxs):
        entry = pos_map.entry(variant.pos).or_insert(PositionGroup { pos: variant.pos, variants: [], source_iter_idxs: [] })
        entry.variants.push(&variant)
        entry.source_iter_idxs.push(source_idx)

    return pos_map.into_values().collect()   // already sorted because BTreeMap
```

**Tests:**
- Single-position group -> one `PositionGroup`
- Multi-position group (deletion + SNPs) -> multiple `PositionGroup`s sorted by pos
- Variants from different iterators at same position are grouped together

---

### Step 2: `detect_complex_variant`

```rust
fn detect_complex_variant(position_groups: &[PositionGroup]) -> bool
```

**Pseudocode:**
```
fn detect_complex_variant(position_groups):
    if position_groups.len() <= 1:
        return false

    for pos_group in position_groups[1..]:
        for variant in pos_group.variants:
            for allele in variant.alleles[1..]:     // skip ref
                if allele != "*" and allele != "." and allele != "<NON_REF>":
                    return true
    return false
```

**Tests:**
- Single position -> false
- Multi-position, all non-first are ref-only -> false
- Multi-position, non-first has real SNP alt -> true
- Multi-position, non-first has only `*` -> false

---

### Step 3: `build_unified_alleles`

```rust
fn build_unified_alleles(
    pos_group: &PositionGroup,
    iter_infos: &[VarIteratorInfo],
    deletion_state: &DeletionState,
    is_first_position: bool,
    max_ref_len: usize,       // longest ref at first position, for extending
    total_samples: usize,
    sample_offsets: &[usize],  // sample_offsets[iter_idx] = first global sample index
) -> UnifiedAlleles
```

**Pseudocode:**
```
fn build_unified_alleles(...):
    // Determine the reference allele
    // At first position of multi-position group: use longest ref
    // At other positions: use the ref from any variant (they should all agree on ref)
    ref_allele = longest ref allele among variants in pos_group
    if is_first_position and max_ref_len > ref_allele.len():
        // Shouldn't happen since max_ref_len comes from this position's variants
        // But handle anyway
        pass

    // Initialize unified allele list with ref
    alleles = [ref_allele.clone()]
    allele_to_idx: HashMap<String, usize> = { ref_allele -> 0 }

    // Initialize per-sample maps: one empty vec per global sample
    per_sample_maps = vec![vec![]; total_samples]
    per_sample_non_ref_idx = vec![None; total_samples]

    // Track which samples have data at this position
    samples_with_data: HashSet<usize> = empty

    // Process each variant in the position group
    for (variant, &iter_idx) in zip(pos_group.variants, pos_group.source_iter_idxs):
        sample_offset = sample_offsets[iter_idx]
        n_samples_in_iter = iter_infos[iter_idx].samples.len()

        for local_sample in 0..n_samples_in_iter:
            global_sample = sample_offset + local_sample
            samples_with_data.insert(global_sample)

            // Build the allele mapping for this sample
            sample_map = vec![0; variant.alleles.len()]

            for (orig_idx, allele) in variant.alleles.iter().enumerate():
                if allele == "<NON_REF>":
                    per_sample_non_ref_idx[global_sample] = Some(orig_idx)
                    // <NON_REF> is NOT added to unified list
                    // Map it to a sentinel; we'll handle it in remap_pls
                    sample_map[orig_idx] = usize::MAX  // sentinel for <NON_REF>
                    continue

                // At first position, extend short alleles to match ref length
                actual_allele = allele
                if is_first_position and orig_idx > 0 and allele.len() < ref_allele.len():
                    // Extend with trailing ref bases
                    // e.g., REF=AC, ALT=C (SNP), ref_allele=AC -> ALT becomes CC
                    trailing = ref_allele[allele.len()..]
                    actual_allele = allele + trailing

                if orig_idx == 0:
                    // Ref allele: always maps to 0
                    sample_map[orig_idx] = 0
                else:
                    // Alt allele: add to unified set if not already present
                    unified_idx = allele_to_idx.get(actual_allele)
                    if unified_idx is None:
                        unified_idx = alleles.len()
                        alleles.push(actual_allele)
                        allele_to_idx.insert(actual_allele, unified_idx)
                    sample_map[orig_idx] = unified_idx

            per_sample_maps[global_sample] = sample_map

    // Handle samples inside active deletions: their allele is *
    // Add * to unified alleles if any sample is in a deletion
    for global_sample in 0..total_samples:
        if deletion_state.any_haplotype_in_deletion(global_sample):
            if "*" not in allele_to_idx:
                star_idx = alleles.len()
                alleles.push("*")
                allele_to_idx.insert("*", star_idx)

    // Samples not present at this position get a ref-only mapping (identity: [0])
    // (handled downstream in build_position_pls; per_sample_maps stays empty for them)

    return UnifiedAlleles { alleles, per_sample_maps, per_sample_non_ref_idx }
```

**Tests:**
- Simple SNP: two samples with different ALTs -> unified list [ref, alt1, alt2], correct mappings
- Allele deduplication: two samples with same ALT -> single entry, both map to it
- Deletion at first position: short alleles extended to match ref length
- Sample inside active deletion -> `*` added to allele list
- Multi-allelic sample -> all its alleles in unified list
- `<NON_REF>` stripped from output but tracked per sample

---

### Step 4: `remap_pls`

```rust
fn remap_pls(
    sample_pls: &[f64],         // PL values in sample's local genotype space
    sample_num_alleles: usize,  // number of alleles for this sample (including <NON_REF>)
    unified_num_alleles: usize, // number of alleles in unified set (no <NON_REF>)
    allele_mapping: &[usize],   // sample allele idx -> unified allele idx (usize::MAX for <NON_REF>)
    non_ref_index: Option<usize>, // index of <NON_REF> in sample's alleles
    ploidy: usize,
) -> Vec<f64>
```

This is the most algorithmically complex function. It maps PLs from a sample's local genotype space to the unified genotype space.

**Pseudocode:**
```
fn remap_pls(...):
    // Enumerate genotypes in unified space
    unified_genotypes = enumerate_genotypes(unified_num_alleles, ploidy)
    num_unified_gts = unified_genotypes.len()

    // Enumerate genotypes in sample's local space
    sample_genotypes = enumerate_genotypes(sample_num_alleles, ploidy)

    // Find max PL in sample (used as penalty for completely unknown genotypes)
    max_pl = max(sample_pls) or 255.0

    // Build reverse mapping: unified_allele_idx -> sample_allele_idx (if exists)
    // For alleles the sample doesn't have, check if <NON_REF> can serve as proxy
    unified_to_sample: Vec<Option<usize>> = vec![None; unified_num_alleles]
    for (sample_idx, &unified_idx) in allele_mapping.iter().enumerate():
        if unified_idx != usize::MAX and unified_idx < unified_num_alleles:
            unified_to_sample[unified_idx] = Some(sample_idx)

    result = vec![0.0; num_unified_gts]

    for (g, unified_gt) in unified_genotypes.iter().enumerate():
        // Try to map this unified genotype back to a sample genotype
        // Convert unified allele_counts to sample allele_counts
        sample_allele_counts = vec![0; sample_num_alleles]
        can_map_exactly = true
        non_ref_count = 0

        for (unified_allele, &count) in unified_gt.allele_counts.iter().enumerate():
            if count == 0:
                continue
            match unified_to_sample[unified_allele]:
                Some(sample_allele):
                    sample_allele_counts[sample_allele] += count
                None:
                    // This allele doesn't exist in the sample
                    // Use <NON_REF> as proxy if available
                    if non_ref_index is Some(nri):
                        sample_allele_counts[nri] += count
                        non_ref_count += count
                    else:
                        can_map_exactly = false

        if can_map_exactly:
            // Find the PL index for these allele counts in the sample's genotype list
            for (sg, sample_gt) in sample_genotypes.iter().enumerate():
                if sample_gt.allele_counts == sample_allele_counts:
                    result[g] = sample_pls[sg]
                    break
        else:
            // No way to represent this genotype in sample's space
            result[g] = max_pl

    return result
```

**Tests:**
- Identity mapping (sample has same alleles as unified) -> PLs unchanged
- Sample has 2 alleles, unified has 3 -> new genotype slots filled from `<NON_REF>` PLs
- Diploid with `<NON_REF>`: genotype involving unknown allele uses `<NON_REF>` likelihood
- No `<NON_REF>` and unknown allele -> max penalty PL used
- Haploid case (ploidy=1)

---

### Step 5: `build_position_pls`

```rust
fn build_position_pls(
    pos_group: &PositionGroup,
    unified: &UnifiedAlleles,
    iter_infos: &[VarIteratorInfo],
    deletion_state: &DeletionState,
    sample_offsets: &[usize],
    total_samples: usize,
    ploidy: usize,
) -> Vec<f64>
```

Builds the full PL matrix for all samples at one position.

**Pseudocode:**
```
fn build_position_pls(...):
    unified_num_alleles = unified.alleles.len()
    num_unified_gts = num_genotypes(unified_num_alleles, ploidy)

    all_pls = vec![0.0; total_samples * num_unified_gts]

    // Build set of which global samples have a variant at this position
    sample_to_variant: HashMap<usize, &Variant> = {}
    for (variant, &iter_idx) in zip(pos_group.variants, pos_group.source_iter_idxs):
        offset = sample_offsets[iter_idx]
        for local in 0..iter_infos[iter_idx].samples.len():
            sample_to_variant[offset + local] = variant

    for global_sample in 0..total_samples:
        pl_offset = global_sample * num_unified_gts

        if deletion_state.any_haplotype_in_deletion(global_sample):
            // Sample is inside a deletion: synthetic hom-* PLs
            // Find index of * in unified alleles
            star_idx = unified.alleles.iter().position(|a| a == "*")
            if star_idx is Some:
                // Generate synthetic PLs: hom-* is the called genotype
                synthetic_gt = vec![star_idx as i8; ploidy]
                synthetic = synthetic_pls_from_gt(&synthetic_gt, 1, unified_num_alleles, ploidy, 0.99)
                all_pls[pl_offset..pl_offset+num_unified_gts].copy_from_slice(&synthetic)
            continue

        variant = sample_to_variant.get(global_sample)
        if variant is None:
            // Sample not present: synthetic hom-ref PLs
            synthetic_gt = vec![0i8; ploidy]
            synthetic = synthetic_pls_from_gt(&synthetic_gt, 1, unified_num_alleles, ploidy, 0.99)
            all_pls[pl_offset..pl_offset+num_unified_gts].copy_from_slice(&synthetic)
            continue

        // Sample has a variant at this position
        allele_mapping = &unified.per_sample_maps[global_sample]
        non_ref_idx = unified.per_sample_non_ref_idx[global_sample]

        // Try to get real PL values from the variant
        pl_field_idx = variant.gt_field_index("PL")
        if pl_field_idx is Some(idx):
            // Parse PL values for this sample
            raw_pl_strings = variant.get_gt_field_by_index(idx)
            // Find which local sample index this global sample corresponds to
            // (determined by iter_idx and sample_offset)
            local_sample_idx = global_sample - sample_offsets[iter_idx_for_this_sample]
            pl_str = raw_pl_strings[local_sample_idx]
            sample_pls = parse pl_str as Vec<f64> (comma-separated)

            // Remap to unified space
            remapped = remap_pls(
                &sample_pls,
                variant.alleles.len(),
                unified_num_alleles,
                allele_mapping,
                non_ref_idx,
                ploidy,
            )
            all_pls[pl_offset..pl_offset+num_unified_gts].copy_from_slice(&remapped)
        else:
            // No PL field: generate synthetic PLs from GT, then remap
            // Extract this sample's GT from the variant
            local_sample_idx = ...
            gt_start = local_sample_idx * ploidy
            sample_gt = variant.genotypes[gt_start..gt_start+ploidy]

            // Generate synthetic PLs in sample's local space
            sample_num_alleles = variant.alleles.len()
            synthetic = synthetic_pls_from_gt(&sample_gt, 1, sample_num_alleles, ploidy, 0.99)

            // Remap to unified space
            remapped = remap_pls(
                &synthetic,
                sample_num_alleles,
                unified_num_alleles,
                allele_mapping,
                non_ref_idx,
                ploidy,
            )
            all_pls[pl_offset..pl_offset+num_unified_gts].copy_from_slice(&remapped)

    return all_pls
```

**Tests:**
- Two samples, one with PLs, one GT-only -> combined flat PL vector of correct size
- Sample absent from position -> gets hom-ref synthetic PLs in unified space
- Sample inside deletion -> gets hom-`*` synthetic PLs
- Real PL values correctly remapped to wider allele set

---

### Step 6: `posteriors_to_calls`

```rust
fn posteriors_to_calls(
    posteriors: &SitePosteriors,
    num_samples: usize,
    ploidy: usize,
    num_alleles: usize,
) -> (Vec<i8>, Vec<u16>, f64)  // (genotypes, gqs, qual)
```

**Pseudocode:**
```
fn posteriors_to_calls(posteriors, num_samples, ploidy, num_alleles):
    genotypes_list = enumerate_genotypes(num_alleles, ploidy)
    num_gts = genotypes_list.len()

    genotypes = vec![0i8; num_samples * ploidy]
    gqs = vec![0u16; num_samples]

    for sample in 0..num_samples:
        post_offset = sample * num_gts

        // Find best and second-best posterior
        best_idx = 0
        best_post = posteriors.genotype_posteriors[post_offset]
        second_best_post = f64::NEG_INFINITY

        for g in 1..num_gts:
            p = posteriors.genotype_posteriors[post_offset + g]
            if p > best_post:
                second_best_post = best_post
                best_post = p
                best_idx = g
            elif p > second_best_post:
                second_best_post = p

        // GT: expand allele_counts of best genotype into sorted allele indices
        gt_offset = sample * ploidy
        pos = 0
        for (allele, &count) in genotypes_list[best_idx].allele_counts.iter().enumerate():
            for _ in 0..count:
                genotypes[gt_offset + pos] = allele as i8
                pos += 1

        // GQ: phred-scaled difference between best and second-best
        // GQ = -10 * log10(1 - best_post) but more practically:
        // GQ = -10 * log10(second_best_post / best_post) ... but since posteriors
        // are already probabilities, GQ = -10 * log10(sum_of_non_best / total)
        // Simplification: GQ = min(99, round(-10 * log10(1 - best_post)))
        if best_post >= 1.0 - 1e-15:
            gqs[sample] = 99
        else:
            raw_gq = -10.0 * (1.0 - best_post).log10()
            gqs[sample] = min(99, round(raw_gq)) as u16

    qual = posteriors.qual
    return (genotypes, gqs, qual)
```

**Tests:**
- Clear hom-alt (posterior ~1.0) -> correct GT, GQ near 99
- Ambiguous call (two genotypes close) -> lower GQ
- All hom-ref -> QUAL near 0

---

### Step 7: `prune_unused_alleles`

```rust
fn prune_unused_alleles(
    alleles: &[String],
    genotypes: &[i8],
    pls: &[f64],
    ploidy: usize,
    num_samples: usize,
) -> (Vec<String>, Vec<i8>, Vec<f64>)
```

**Pseudocode:**
```
fn prune_unused_alleles(alleles, genotypes, pls, ploidy, num_samples):
    // Find which allele indices are actually used
    used_alleles: HashSet<usize> = {0}  // ref always kept
    for &a in genotypes:
        if a >= 0:
            used_alleles.insert(a as usize)

    if used_alleles.len() == alleles.len():
        return (alleles.to_vec(), genotypes.to_vec(), pls.to_vec())  // nothing to prune

    // Build old->new index mapping
    new_idx = 0
    old_to_new: Vec<Option<usize>> = vec![None; alleles.len()]
    new_alleles = []
    for (old_idx, allele) in alleles.iter().enumerate():
        if old_idx in used_alleles:
            old_to_new[old_idx] = Some(new_idx)
            new_alleles.push(allele.clone())
            new_idx += 1

    new_num_alleles = new_alleles.len()

    // Remap genotypes
    new_genotypes = genotypes.iter().map(|&a|
        if a < 0: a
        else: old_to_new[a as usize].unwrap() as i8
    ).collect()

    // Recompute PLs: keep only genotype columns for the new allele set
    old_gts = enumerate_genotypes(alleles.len(), ploidy)
    new_gts = enumerate_genotypes(new_num_alleles, ploidy)
    num_old_gts = old_gts.len()
    num_new_gts = new_gts.len()

    new_pls = vec![0.0; num_samples * num_new_gts]
    for sample in 0..num_samples:
        for (new_g, new_gt) in new_gts.iter().enumerate():
            // Find the corresponding old genotype
            // Convert new allele counts to old allele counts
            old_allele_counts = vec![0; alleles.len()]
            for (new_allele, &count) in new_gt.allele_counts.iter().enumerate():
                // new_allele -> which old allele?
                // We need reverse of old_to_new. Since we know which alleles were kept:
                old_allele = kept_alleles_old_indices[new_allele]
                old_allele_counts[old_allele] = count

            // Find old genotype with these counts
            for (old_g, old_gt) in old_gts.iter().enumerate():
                if old_gt.allele_counts == old_allele_counts:
                    new_pls[sample * num_new_gts + new_g] = pls[sample * num_old_gts + old_g]
                    break

    return (new_alleles, new_genotypes, new_pls)
```

**Tests:**
- No unused alleles -> unchanged
- One unused alt allele -> removed, indices shifted, PL columns dropped
- Ref always kept even if all samples are alt

---

### Step 8: `update_deletion_state`

```rust
fn update_deletion_state(
    genotypes: &[i8],
    alleles: &[String],
    deletion_state: &mut DeletionState,
    total_samples: usize,
    ploidy: usize,
)
```

**Pseudocode:**
```
fn update_deletion_state(genotypes, alleles, deletion_state, total_samples, ploidy):
    ref_len = alleles[0].len()

    // For each sample+haplotype, check if the called allele starts a new deletion
    for sample in 0..total_samples:
        for h in 0..ploidy:
            allele_idx = genotypes[sample * ploidy + h]
            if allele_idx < 0:
                continue

            allele = &alleles[allele_idx as usize]
            if allele == "*":
                // Already in deletion, state was already decremented
                continue

            del_len = ref_len as i32 - allele.len() as i32
            if del_len > 0:
                // This allele is a deletion: it spans del_len more positions
                deletion_state.set_deletion(sample, h, del_len)
            else:
                // Not a deletion (SNP, insertion, ref): clear any leftover
                // (shouldn't happen if logic is correct, but defensive)
                deletion_state.set_deletion(sample, h, 0)

    // Already tracked; next position's build_unified_alleles will use the state
```

No dedicated tests (tested via orchestrator integration tests).

---

### Step 9: `assign_tentative_genotypes`

```rust
fn assign_tentative_genotypes(
    pos_group: &PositionGroup,
    unified: &UnifiedAlleles,
    deletion_state: &DeletionState,
    iter_infos: &[VarIteratorInfo],
    sample_offsets: &[usize],
    total_samples: usize,
    ploidy: usize,
) -> (Vec<i8>, Vec<bool>)  // (genotypes, phases)
```

Maps each sample's original GT to the unified allele indices. Also determines phase.

**Pseudocode:**
```
fn assign_tentative_genotypes(...):
    genotypes = vec![-1i8; total_samples * ploidy]
    phases = vec![true; total_samples]

    // Samples with a variant at this position
    for (variant, &iter_idx) in zip(pos_group.variants, pos_group.source_iter_idxs):
        offset = sample_offsets[iter_idx]
        n_local = iter_infos[iter_idx].samples.len()

        for local in 0..n_local:
            global = offset + local
            allele_map = &unified.per_sample_maps[global]
            gt_start = local * ploidy

            for h in 0..ploidy:
                orig_allele = variant.genotypes[gt_start + h]
                if orig_allele < 0:
                    genotypes[global * ploidy + h] = -1
                else:
                    mapped = allele_map[orig_allele as usize]
                    if mapped == usize::MAX:
                        // Was <NON_REF>, treat as ref
                        genotypes[global * ploidy + h] = 0
                    else:
                        genotypes[global * ploidy + h] = mapped as i8

            phases[global] = variant.phases[local]

    // Samples inside a deletion (not present at this position via a variant)
    for global in 0..total_samples:
        if genotypes[global * ploidy] == -1 and deletion_state.any_haplotype_in_deletion(global):
            // Find * index in unified alleles
            star_idx = unified.alleles.iter().position(|a| a == "*")
            if star_idx is Some(si):
                for h in 0..ploidy:
                    if deletion_state.is_in_deletion(global, h):
                        genotypes[global * ploidy + h] = si as i8
                    else:
                        genotypes[global * ploidy + h] = 0  // ref

    // Samples completely absent: leave as hom-ref (0)
    for global in 0..total_samples:
        if genotypes[global * ploidy] == -1:
            for h in 0..ploidy:
                genotypes[global * ploidy + h] = 0
            phases[global] = false  // unknown phase

    return (genotypes, phases)
```

No dedicated tests (tested via orchestrator).

---

### Step 10: `join_genotypes` (orchestrator)

```rust
pub fn join_genotypes(
    group: &OverlappingVarGroup,
    iter_infos: &[VarIteratorInfo],
    prior: &PriorConfig,
    ploidy: usize,
    remove_complex: bool,
) -> VcfResult<Vec<Variant>>
```

**Pseudocode:**
```
fn join_genotypes(group, iter_infos, prior, ploidy, remove_complex):
    // Compute sample offsets
    sample_offsets = compute_sample_offsets(iter_infos)
    total_samples = sum of all iter_infos[i].samples.len()

    // Step 1: group by position
    position_groups = group_variants_by_position(group)

    // Step 2: detect complex variants
    if remove_complex and detect_complex_variant(&position_groups):
        return Err(VcfParseError::ComplexVariantRemoved)

    // Initialize deletion state
    deletion_state = DeletionState::new(total_samples, ploidy)

    // Determine max ref length at first position (for allele extension)
    max_ref_len = max(v.alleles[0].len() for v in position_groups[0].variants)

    output_variants = []

    for (i, pos_group) in position_groups.iter().enumerate():
        is_first = (i == 0)

        // Step 3a: build unified alleles
        unified = build_unified_alleles(
            pos_group, iter_infos, &deletion_state,
            is_first, max_ref_len, total_samples, &sample_offsets,
        )

        // Step 3b: build PLs for all samples
        position_pls = build_position_pls(
            pos_group, &unified, iter_infos, &deletion_state,
            &sample_offsets, total_samples, ploidy,
        )

        // Step 3c: assign tentative genotypes (for phase preservation)
        (tentative_gts, phases) = assign_tentative_genotypes(
            pos_group, &unified, &deletion_state,
            iter_infos, &sample_offsets, total_samples, ploidy,
        )

        num_alleles = unified.alleles.len()

        // Step 3d: run EM
        posteriors = estimate_posteriors(num_alleles, ploidy, total_samples, &position_pls, prior)

        // Step 3e: derive GT, GQ, QUAL from posteriors
        (em_genotypes, gqs, qual) = posteriors_to_calls(&posteriors, total_samples, ploidy, num_alleles)

        // Prefer tentative GTs when EM agrees (preserves phase)
        final_genotypes = vec![0i8; total_samples * ploidy]
        for sample in 0..total_samples:
            tentative_counts = count_alleles(tentative_gts, sample, ploidy, num_alleles)
            em_counts = count_alleles(em_genotypes, sample, ploidy, num_alleles)
            if tentative_counts == em_counts:
                // EM agrees: keep tentative (preserves phase/order)
                copy tentative_gts[sample*ploidy..(sample+1)*ploidy] -> final_genotypes
            else:
                // EM changed the call: use EM genotype
                copy em_genotypes[sample*ploidy..(sample+1)*ploidy] -> final_genotypes
                phases[sample] = false  // phase no longer reliable

        // Step 3f: prune unused alleles
        (pruned_alleles, pruned_gts, pruned_pls) = prune_unused_alleles(
            &unified.alleles, &final_genotypes, &position_pls, ploidy, total_samples,
        )

        // Step 3g: skip invariant
        if pruned_alleles.len() <= 1:
            // Update deletion state before skipping
            // (need to decrement counters even if we don't emit)
            deletion_state.decrement_all()
            update_deletion_state(&final_genotypes, &unified.alleles, &mut deletion_state, total_samples, ploidy)
            continue

        // Step 3h: update deletion state
        deletion_state.decrement_all()
        update_deletion_state(&pruned_gts, &pruned_alleles, &mut deletion_state, total_samples, ploidy)

        // Step 3i: build output Variant
        variant = Variant::new(
            group.chrom.clone(),
            pos_group.pos,
            pruned_alleles,
            pruned_gts,
            phases,
            total_samples,
        )
        variant.qual = qual as f32
        // TODO: store GQ, PL in variant (requires extending Variant or output format)
        // For now, we could store them in sample_gt_fields or add new fields

        output_variants.push(variant)

    return Ok(output_variants)
```

**Integration tests:**
- Simple SNP joining (spec example): two samples, different ALTs at same position
- Deletion + SNP (spec example): one sample has deletion, other has SNP at same position
- Complex variant detection + removal with `remove_complex=true`
- Complex variant passes through with `remove_complex=false`
- Multi-position group producing multiple output variants
- Phase preservation through the full pipeline
- All-homref group -> empty output (no variants emitted)
- Single sample, single SNP -> trivial join

---

## New error variant needed

Add to `VcfParseError`:

```rust
#[error("Complex variant removed at {chrom}:{pos}")]
ComplexVariantRemoved { chrom: String, pos: u32 },
```

## Implementation sequence summary

| Order | Function                      | Dedicated tests? | Depends on     |
|-------|-------------------------------|-------------------|----------------|
| 1     | `DeletionState`               | Yes               | -              |
| 2     | `group_variants_by_position`  | Yes               | -              |
| 3     | `detect_complex_variant`      | Yes               | Step 2         |
| 4     | `build_unified_alleles`       | Yes               | Step 1         |
| 5     | `remap_pls`                   | Yes               | -              |
| 6     | `build_position_pls`          | Yes               | Steps 4, 5    |
| 7     | `posteriors_to_calls`         | Yes               | -              |
| 8     | `prune_unused_alleles`        | Yes               | -              |
| 9     | `update_deletion_state`       | No (via Step 11)  | Step 1         |
| 10    | `assign_tentative_genotypes`  | No (via Step 11)  | Step 4         |
| 11    | `join_genotypes`              | Yes (integration) | Steps 1-10     |

## Deferred work (not part of this plan)

- Extending `Variant` struct to carry GQ/PL data natively (currently GT-only fields).
- Extending `VcfWriter` to output GQ, PL, DP, AD FORMAT fields.
- Wiring `join_genotypes` into `pipeline.rs` (replacing `genotype_merging`).
- Adding `--remove-complex` / `--join-mode` CLI flags to `main.rs`.
- Parsing DP and AD fields from input gVCFs.
- Removing `genotype_merging.rs` once joining is fully validated.

## Files affected

- **New:** `src/genotype_joining.rs` - all joining logic
- **New:** `tests/genotype_joining_test.rs` - unit and integration tests
- **Modified:** `src/errors.rs` - add `ComplexVariantRemoved` variant
- **Modified:** `src/lib.rs` - add `pub mod genotype_joining;`
