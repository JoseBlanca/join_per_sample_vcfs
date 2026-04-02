use std::collections::{HashMap, HashSet};
use std::sync::mpsc;
use std::thread;

use rayon::prelude::*;

use crate::errors::VcfParseError;
use crate::genotype_posteriors::{self, PriorConfig};
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
    // Track per-sample, per-haplotype missing status
    let mut missing_haplotypes: Vec<Vec<bool>> = vec![Vec::new(); total_samples];
    // Phase output: true unless first variant is unphased, any het is unphased, or sample is missing
    let mut sample_phase = vec![true; total_samples];

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
            let sample_phase_input = variant.phase[sample_idx_in_var_iter];

            let is_first = alleles_for_samples[sample_idx].is_none();
            if is_first {
                alleles_for_samples[sample_idx] = Some(vec![Vec::new(); ploidy]);
                positions_left_in_del[sample_idx] = vec![0i32; ploidy];
                missing_haplotypes[sample_idx] = vec![false; ploidy];
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

            // Phase tracking: detect when haplotype can't be built due to broken phase.
            // If phase is false after the first het was seen, we can't determine
            // haplotype assignment for the next het.
            let is_het = sample_gt.iter().any(|&a| a != sample_gt[0]);
            if !is_first && first_het_seen[sample_idx] && !sample_phase_input {
                phase_broken_since_het[sample_idx] = true;
            }
            if is_het {
                if first_het_seen[sample_idx] && phase_broken_since_het[sample_idx] {
                    missing_samples.insert(sample_idx);
                }
                first_het_seen[sample_idx] = true;
                phase_broken_since_het[sample_idx] = false;
            }

            // Output phase: false if first variant unphased or any het unphased
            if is_first && !sample_phase_input {
                sample_phase[sample_idx] = false;
            }
            if is_het && !sample_phase_input {
                sample_phase[sample_idx] = false;
            }

            // Build alleles per haplotype
            let haplo_alleles = alleles_for_samples[sample_idx].as_mut().unwrap();
            let pos_left = &mut positions_left_in_del[sample_idx];

            for (h, &allele_int) in sample_gt.iter().enumerate() {
                // Missing alleles (-1) are treated as ref for allele building,
                // but the haplotype is flagged as missing for the final genotype
                let allele_idx = if allele_int < 0 {
                    missing_haplotypes[sample_idx][h] = true;
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
            for (h, allele) in alleles_str_per_sample[sample_idx].iter().enumerate() {
                if missing_haplotypes[sample_idx][h] {
                    genotypes.push(-1i8);
                } else {
                    genotypes.push(*allele_id_map.get(allele).unwrap() as i8);
                }
            }
        }
    }

    // Missing samples get phase=false
    for &si in &missing_samples {
        sample_phase[si] = false;
    }
    let phases = sample_phase;

    Ok(Variant {
        chrom: var_group.chrom.clone(),
        pos: var_group.start,
        alleles: alleles_vec,
        ref_allele_len: ref_allele.len() as u8,
        qual: f32::NAN,
        genotypes,
        phase: phases,
        gt_format_fields: Vec::new(),
        sample_gt_fields: Vec::new(),
        n_samples: total_samples,
    })
}

/// Merges a single `OverlappingVariantGroup` into one or more output variants.
///
/// After merging alleles and genotypes, computes genotype posterior
/// probabilities using the EM algorithm.  If the input variants have PL data
/// it is used directly; otherwise synthetic PLs are generated from the GT
/// calls (99% confidence on the called genotype).
pub fn merge_variant_group(
    group: &OverlappingVariantGroup,
    iter_info: &[VariantIteratorInfo],
    prior: &PriorConfig,
) -> VcfResult<Vec<Variant>> {
    let mut variant = create_variant_for_region(group, iter_info)?;

    // A variant is non-variable if it has no alt alleles and no input variants
    // had alt alleles (i.e. missing genotypes from broken phase don't count,
    // but missing from ref-only inputs should be filtered).
    let any_input_has_alt = group
        .variants
        .iter()
        .any(|v| v.alleles.len() > 1 && v.genotypes.iter().any(|&g| g > 0));
    if variant.alleles.len() <= 1 && !any_input_has_alt {
        return Err(VcfParseError::NotVariable);
    }

    compute_posteriors_for_variant(&mut variant, prior);

    Ok(vec![variant])
}

/// Computes genotype posteriors for a merged variant and updates its QUAL and
/// GT fields.
fn compute_posteriors_for_variant(variant: &mut Variant, prior: &PriorConfig) {
    let num_alleles = variant.alleles.len();
    let num_samples = variant.n_samples;
    if num_alleles == 0 || num_samples == 0 {
        return;
    }

    // Determine ploidy from the genotypes vector
    let ploidy = variant.genotypes.len() / num_samples;
    if ploidy == 0 {
        return;
    }

    // Generate synthetic PLs from GT calls (99% confidence on the called genotype)
    let sample_pls = genotype_posteriors::synthetic_pls_from_gt(
        &variant.genotypes,
        num_samples,
        num_alleles,
        ploidy,
        0.99,
    );

    let posteriors = genotype_posteriors::estimate_posteriors(
        num_alleles,
        ploidy,
        num_samples,
        &sample_pls,
        prior,
    );

    // Update QUAL score
    variant.qual = posteriors.qual as f32;

    // Update genotypes: assign each sample the genotype with the highest
    // posterior probability
    let genotypes_list = genotype_posteriors::enumerate_genotypes(num_alleles, ploidy);
    let num_genotypes = genotypes_list.len();

    for sample in 0..num_samples {
        let gt_offset = sample * ploidy;

        // Skip samples with any missing allele — the EM cannot reliably
        // reassign partial genotypes, so we preserve the original call.
        let has_missing = (0..ploidy).any(|i| variant.genotypes[gt_offset + i] < 0);
        if has_missing {
            continue;
        }

        let post_offset = sample * num_genotypes;

        // Find the genotype with the highest posterior
        let mut best_gt_index = 0;
        let mut best_posterior = posteriors.genotype_posteriors[post_offset];
        for g in 1..num_genotypes {
            let p = posteriors.genotype_posteriors[post_offset + g];
            if p > best_posterior {
                best_posterior = p;
                best_gt_index = g;
            }
        }

        // Check if the EM changed the genotype.  If the allele counts match
        // the original GT, keep the original ordering (preserves phase like
        // 1|0 vs 0|1).  Only overwrite when the EM picked a different genotype.
        let best_genotype = &genotypes_list[best_gt_index];
        let mut original_counts = vec![0usize; num_alleles];
        for i in 0..ploidy {
            let a = variant.genotypes[gt_offset + i];
            if a >= 0 {
                original_counts[a as usize] += 1;
            }
        }

        if original_counts != best_genotype.allele_counts {
            // The EM changed the call — write the new genotype in sorted order
            let mut pos = 0;
            for (allele, &count) in best_genotype.allele_counts.iter().enumerate() {
                for _ in 0..count {
                    variant.genotypes[gt_offset + pos] = allele as i8;
                    pos += 1;
                }
            }
        }
        // Otherwise: keep the original GT (preserving allele order / phase)
    }

}

/// Collect up to `batch_size` groups from the source iterator.
fn collect_batch<I>(groups: &mut I, batch_size: usize) -> Vec<VcfResult<OverlappingVariantGroup>>
where
    I: Iterator<Item = VcfResult<OverlappingVariantGroup>>,
{
    let mut batch = Vec::with_capacity(batch_size);
    for group_result in groups.by_ref() {
        batch.push(group_result);
        if batch.len() >= batch_size {
            break;
        }
    }
    batch
}

/// Process a batch of groups in parallel with rayon.
fn process_batch(
    batch: Vec<VcfResult<OverlappingVariantGroup>>,
    iter_info: &[VariantIteratorInfo],
    prior: &PriorConfig,
) -> Vec<VcfResult<Variant>> {
    batch
        .into_par_iter()
        .flat_map(|group_result| match group_result {
            Ok(group) => match merge_variant_group(&group, iter_info, prior) {
                Ok(variants) => variants.into_iter().map(Ok).collect(),
                Err(VcfParseError::NotVariable) => vec![],
                Err(e) => vec![Err(e)],
            },
            Err(e) => vec![Err(e)],
        })
        .collect()
}

/// Processes `OverlappingVariantGroup`s from an iterator and writes merged
/// variants to a callback, one batch at a time.
///
/// A background thread continuously collects groups into batches while the main
/// thread processes the previous batch in parallel with rayon. This overlaps
/// the sequential group-collection work with the parallel merge computation,
/// keeping all cores busy.
///
/// Uses `thread::scope` so the collector thread can borrow the group iterator
/// without requiring `'static`.
pub fn merge_vars_in_groups<I, F>(
    groups: I,
    iter_info: &[VariantIteratorInfo],
    prior: &PriorConfig,
    mut on_variant: F,
) -> VcfResult<()>
where
    I: Iterator<Item = VcfResult<OverlappingVariantGroup>> + Send,
    F: FnMut(VcfResult<Variant>) -> VcfResult<()>,
{
    let batch_size = DEFAULT_BATCH_SIZE;

    thread::scope(|scope| {
        // Bounded channel with capacity 1: the collector can be at most one
        // batch ahead, bounding memory to ~2 batches.
        let (batch_sender, batch_receiver) =
            mpsc::sync_channel::<Vec<VcfResult<OverlappingVariantGroup>>>(1);

        // Background collector: pulls groups into batches.
        scope.spawn(move || {
            let mut groups = groups;
            loop {
                let batch = collect_batch(&mut groups, batch_size);
                if batch.is_empty() {
                    break;
                }
                if batch_sender.send(batch).is_err() {
                    break;
                }
            }
        });

        // Main thread: processes batches in parallel with rayon, then
        // streams each result to the callback.
        while let Ok(batch) = batch_receiver.recv() {
            if batch.is_empty() {
                break;
            }
            let results = process_batch(batch, iter_info, prior);
            for result in results {
                on_variant(result)?;
            }
        }

        Ok(())
    })
}
