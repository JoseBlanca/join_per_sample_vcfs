use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::errors::VcfParseError;
use crate::gvcf_parser::{Variant, VcfResult};
use crate::variant_grouping::{OverlappingVariantGroup, VariantIteratorInfo};

/// Default number of groups to process in each parallel batch.
const DEFAULT_BATCH_SIZE: usize = 1000;

/// Creates a merged variant for the region covered by a variant group.
///
/// This is the core merging algorithm: it takes all overlapping variants from
/// multiple samples and produces a single merged variant with unified alleles
/// and genotypes across all samples.
fn create_variant_for_region(
    var_group: &OverlappingVariantGroup,
    var_iter_infos: &[VariantIteratorInfo],
) -> VcfResult<Variant> {
    // Compute sample offsets: sample_idx = sample_offsets[var_iter_idx] + sample_idx_in_var_iter
    let n_samples_per_var_iter: Vec<usize> = var_iter_infos
        .iter()
        .map(|info| info.samples.len())
        .collect();
    let sample_offsets: Vec<usize> = n_samples_per_var_iter
        .iter()
        .scan(0usize, |acc, &n| {
            let offset = *acc;
            *acc += n;
            Some(offset)
        })
        .collect();
    let total_samples: usize = n_samples_per_var_iter.iter().sum();

    // Per-sample state
    // alleles_for_samples[sample_idx] = per-haplotype list of allele string fragments
    let mut alleles_for_samples: Vec<Option<Vec<Vec<String>>>> = vec![None; total_samples];
    let mut positions_left_in_del: Vec<Vec<i32>> = vec![Vec::new(); total_samples];
    let mut first_het_seen = vec![false; total_samples];
    let mut phase_broken_since_het = vec![false; total_samples];
    let mut missing_samples: HashSet<usize> = HashSet::new();

    // Process each variant in the group
    for (variant, &var_iter_idx) in var_group
        .variants
        .iter()
        .zip(var_group.source_var_iter_idxs.iter())
    {
        let ploidy = variant.genotypes.len() / variant.n_samples;
        let ref_allele = &variant.alleles[0];

        for sample_idx_in_var_iter in 0..variant.n_samples {
            let sample_idx = sample_offsets[var_iter_idx] + sample_idx_in_var_iter;
            let gt_start = sample_idx_in_var_iter * ploidy;
            let sample_gt = &variant.genotypes[gt_start..gt_start + ploidy];
            let sample_phase = variant.phase[sample_idx_in_var_iter];

            let is_first = alleles_for_samples[sample_idx].is_none();
            if is_first {
                alleles_for_samples[sample_idx] = Some(vec![Vec::new(); ploidy]);
                positions_left_in_del[sample_idx] = vec![0i32; ploidy];
            }

            // Save previous deletion state before potential update
            let prev_pos_left = positions_left_in_del[sample_idx].clone();
            let deletion_created = if ref_allele.len() > 1 {
                let pos_left = &mut positions_left_in_del[sample_idx];
                for (h, &allele_int) in sample_gt.iter().enumerate() {
                    // Missing alleles (-1) are treated as ref for deletion tracking
                    let allele_idx = if allele_int < 0 {
                        0usize
                    } else {
                        allele_int as usize
                    };
                    let sample_allele = &variant.alleles[allele_idx];
                    let del_len = ref_allele.len() as i32 - sample_allele.len() as i32;
                    pos_left[h] = del_len + prev_pos_left[h];
                }
                true
            } else {
                false
            };

            // Phase tracking: detect when haplotype can't be built due to broken phase
            let is_het = sample_gt.iter().any(|&a| a != sample_gt[0]);
            if !is_first && first_het_seen[sample_idx] && !sample_phase {
                phase_broken_since_het[sample_idx] = true;
            }
            if is_het {
                if first_het_seen[sample_idx] && phase_broken_since_het[sample_idx] {
                    missing_samples.insert(sample_idx);
                }
                first_het_seen[sample_idx] = true;
                phase_broken_since_het[sample_idx] = false;
            }

            // Build alleles per haplotype
            let haplo_alleles = alleles_for_samples[sample_idx].as_mut().unwrap();
            let pos_left = &mut positions_left_in_del[sample_idx];

            for (h, &allele_int) in sample_gt.iter().enumerate() {
                // Missing alleles (-1) are treated as ref for allele building
                let allele_idx = if allele_int < 0 {
                    0usize
                } else {
                    allele_int as usize
                };
                let full_allele = &variant.alleles[allele_idx];

                // For multi-base ref (deletion context), keep only first character
                // For single-base ref (SNP/insertion), keep the full allele
                let fragment = if ref_allele.len() > 1 {
                    full_allele[..1].to_string()
                } else {
                    full_allele.clone()
                };

                // Skip positions consumed by a prior deletion
                if pos_left[h] > 0 && (!deletion_created || prev_pos_left[h] > 0) {
                    haplo_alleles[h].push(String::new());
                    pos_left[h] -= 1;
                } else {
                    haplo_alleles[h].push(fragment);
                }
            }
        }
    }

    // Build reference allele from first character of each variant's ref allele
    let mut ref_parts_per_iter: Vec<String> = vec![String::new(); var_iter_infos.len()];
    for (variant, &var_iter_idx) in var_group
        .variants
        .iter()
        .zip(var_group.source_var_iter_idxs.iter())
    {
        ref_parts_per_iter[var_iter_idx].push_str(&variant.alleles[0][..1]);
    }

    // Use the longest ref allele across var_iters (different var_iters may
    // contribute different numbers of variants to the group)
    let ref_allele = ref_parts_per_iter
        .iter()
        .filter(|s| !s.is_empty())
        .max_by_key(|s| s.len())
        .ok_or_else(|| VcfParseError::RuntimeError {
            message: "No variants in group".into(),
        })?
        .clone();

    // Build allele-to-ID mapping (ref allele is always ID 0)
    let mut allele_id_map: HashMap<String, usize> = HashMap::new();
    allele_id_map.insert(ref_allele.clone(), 0);

    let mut alleles_str_per_sample: Vec<Vec<String>> = Vec::with_capacity(total_samples);
    for sample_idx in 0..total_samples {
        if missing_samples.contains(&sample_idx) || alleles_for_samples[sample_idx].is_none() {
            alleles_str_per_sample.push(Vec::new());
            if alleles_for_samples[sample_idx].is_none() {
                missing_samples.insert(sample_idx);
            }
            continue;
        }
        let haplo_parts = alleles_for_samples[sample_idx].as_ref().unwrap();
        let mut sample_alleles = Vec::with_capacity(haplo_parts.len());
        for parts in haplo_parts {
            let allele: String = parts.concat();
            let next_id = allele_id_map.len();
            allele_id_map.entry(allele.clone()).or_insert(next_id);
            sample_alleles.push(allele);
        }
        alleles_str_per_sample.push(sample_alleles);
    }

    // Build ordered alleles vector (index = allele ID)
    let mut alleles_vec = vec![String::new(); allele_id_map.len()];
    for (allele, &id) in &allele_id_map {
        alleles_vec[id] = allele.clone();
    }

    // Determine ploidy from first variant
    let ploidy = if let Some(v) = var_group.variants.first() {
        v.genotypes.len() / v.n_samples
    } else {
        2
    };

    // Build flat genotypes vector
    let mut genotypes = Vec::with_capacity(total_samples * ploidy);
    for sample_idx in 0..total_samples {
        if missing_samples.contains(&sample_idx) {
            genotypes.extend(std::iter::repeat(-1i8).take(ploidy));
        } else {
            for allele in &alleles_str_per_sample[sample_idx] {
                genotypes.push(*allele_id_map.get(allele).unwrap() as i8);
            }
        }
    }

    let phases = vec![false; total_samples];

    Ok(Variant {
        chrom: var_group.chrom.clone(),
        pos: var_group.start,
        alleles: alleles_vec,
        ref_allele_len: ref_allele.len() as u8,
        qual: f32::NAN,
        gt_index: 0,
        pl_index: None,
        genotypes,
        phase: phases,
        n_samples: total_samples,
    })
}

/// Merges a single `OverlappingVariantGroup` into one or more output variants.
pub fn merge_variant_group(
    group: &OverlappingVariantGroup,
    iter_info: &[VariantIteratorInfo],
) -> VcfResult<Vec<Variant>> {
    let variant = create_variant_for_region(group, iter_info)?;

    let has_missing = variant.genotypes.iter().any(|&g| g < 0);
    if variant.alleles.len() <= 1 && !has_missing {
        return Err(VcfParseError::NotVariable);
    }
    Ok(vec![variant])
}

/// Iterator that pulls groups in batches, processes each batch in parallel
/// with rayon, and yields merged variants one at a time.
///
/// This gives both streaming output (bounded memory) and multi-core throughput.
struct BatchedParallelIter<'a, I> {
    groups: I,
    iter_info: &'a [VariantIteratorInfo],
    batch_size: usize,
    buffer: std::vec::IntoIter<VcfResult<Variant>>,
}

impl<'a, I> Iterator for BatchedParallelIter<'a, I>
where
    I: Iterator<Item = VcfResult<OverlappingVariantGroup>>,
{
    type Item = VcfResult<Variant>;

    fn next(&mut self) -> Option<Self::Item> {
        // Drain the current batch buffer first
        if let Some(item) = self.buffer.next() {
            return Some(item);
        }

        // Pull the next batch from the source iterator
        let mut batch = Vec::with_capacity(self.batch_size);
        for group_result in self.groups.by_ref() {
            batch.push(group_result);
            if batch.len() >= self.batch_size {
                break;
            }
        }

        if batch.is_empty() {
            return None;
        }

        // Process the batch in parallel (rayon preserves input order)
        let iter_info = self.iter_info;
        let results: Vec<VcfResult<Variant>> = batch
            .into_par_iter()
            .flat_map(|group_result| match group_result {
                Ok(group) => match merge_variant_group(&group, iter_info) {
                    Ok(variants) => variants.into_iter().map(Ok).collect(),
                    Err(VcfParseError::NotVariable) => vec![],
                    Err(e) => vec![Err(e)],
                },
                Err(e) => vec![Err(e)],
            })
            .collect();

        self.buffer = results.into_iter();
        self.buffer.next()
    }
}

/// Processes `OverlappingVariantGroup`s from an iterator, returning a streaming
/// iterator of merged variants.
///
/// Groups are pulled in batches and processed in parallel with rayon, then
/// yielded one at a time. This bounds memory to one batch while using all
/// available cores for the merge computation.
pub fn merge_vars_in_groups<'a, I>(
    groups: I,
    iter_info: &'a [VariantIteratorInfo],
) -> impl Iterator<Item = VcfResult<Variant>> + 'a
where
    I: Iterator<Item = VcfResult<OverlappingVariantGroup>> + 'a,
{
    BatchedParallelIter {
        groups,
        iter_info,
        batch_size: DEFAULT_BATCH_SIZE,
        buffer: Vec::new().into_iter(),
    }
}
