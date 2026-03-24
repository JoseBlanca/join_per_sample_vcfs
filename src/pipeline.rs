use std::io::{BufRead, Write};

use crate::genotype_merging::merge_vars_in_groups;
use crate::gvcf_parser::{VariantIterator, VcfResult};
use crate::variant_grouping::VariantGroupIterator;
use crate::vcf_writer::VcfWriter;

/// Runs the full gVCF merging pipeline: groups overlapping variants across
/// samples, merges their alleles and genotypes, and writes the result as a
/// multi-sample VCF.
pub fn merge_alleles_and_genotypes<B: BufRead + Send>(
    var_iters: Vec<VariantIterator<B>>,
    sorted_chromosomes: Vec<String>,
    writer: Box<dyn Write>,
) -> VcfResult<()> {
    let all_samples: Vec<String> = var_iters
        .iter()
        .flat_map(|iter| iter.samples().iter().cloned())
        .collect();

    let grouper = VariantGroupIterator::new(var_iters, sorted_chromosomes)?;
    let iter_info = grouper.iter_info().to_vec();

    let mut vcf_writer =
        VcfWriter::from_writer(writer, &all_samples).map_err(crate::errors::VcfParseError::from)?;

    merge_vars_in_groups(grouper, &iter_info, |result| {
        let variant = result?;
        vcf_writer
            .write_variant(&variant)
            .map_err(crate::errors::VcfParseError::from)
    })?;

    vcf_writer
        .flush()
        .map_err(crate::errors::VcfParseError::from)?;

    Ok(())
}
