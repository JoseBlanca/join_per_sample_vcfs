use rayon::prelude::*;

use crate::gvcf_parser::{Variant, VcfResult};
use crate::variant_group::{OverlappingVariantGroup, VariantIteratorInfo};

/// Analyzes a single `OverlappingVariantGroup` and produces the merged variant(s).
///
/// Currently a mock: copies the first variant in the group.
/// The real implementation will combine genotypes across samples.
fn analyze_group(
    group: &OverlappingVariantGroup,
    _iter_info: &[VariantIteratorInfo],
) -> VcfResult<Vec<Variant>> {
    let first = &group.variants[0];
    let variant = Variant {
        chrom: first.chrom.clone(),
        pos: first.pos,
        alleles: first.alleles.clone(),
        ref_allele_len: first.ref_allele_len,
        qual: first.qual,
        gt_index: 0,
        pl_index: None,
        genotypes: first.genotypes.clone(),
        phase: first.phase.clone(),
        n_samples: first.n_samples,
    };
    Ok(vec![variant])
}

/// Processes all `OverlappingVariantGroup`s in parallel, returning an ordered
/// sequence of merged variants.
///
/// Groups are analyzed in parallel using rayon's thread pool, and results
/// are collected in the original group order.
pub fn analyze_groups(
    groups: Vec<OverlappingVariantGroup>,
    iter_info: &[VariantIteratorInfo],
) -> VcfResult<Vec<Variant>> {
    let results: Vec<VcfResult<Vec<Variant>>> = groups
        .par_iter()
        .map(|group| analyze_group(group, iter_info))
        .collect();

    let mut variants = Vec::new();
    for result in results {
        variants.extend(result?);
    }
    Ok(variants)
}
