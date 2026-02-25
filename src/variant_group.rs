use std::collections::HashSet;
use std::io::BufRead;

use crate::errors::VcfParseError;
use crate::gvcf_parser::{GVcfRecord, GVcfRecordIterator, VcfResult};

/// A bin of overlapping variants collected from multiple per-sample VCFs.
///
/// The bin covers a contiguous genomic region where all contained variants
/// have reference-allele spans that overlap with the bin's span.
#[derive(Debug)]
pub struct OverlappingVariantGroup {
    pub chrom: String,
    /// Start position of the bin (1-based, inclusive).
    pub start: u32,
    /// End position of the bin (1-based, inclusive).
    pub end: u32,
    /// All variant records in this bin, in consumption order.
    pub records: Vec<GVcfRecord>,
}

impl OverlappingVariantGroup {
    pub fn span(&self) -> (&str, u32, u32) {
        (&self.chrom, self.start, self.end)
    }
}

/// Computes the reference-allele span end for a variant (1-based, inclusive).
#[inline]
fn ref_span_end(record: &GVcfRecord) -> u32 {
    record.pos + record.ref_allele_len as u32 - 1
}

/// Iterator that merges multiple per-sample VCF streams, grouping
/// overlapping variants into bins based on reference-allele span overlap.
///
/// Variants from all input VCFs are consumed in genomic order. Two variants
/// overlap when they are on the same chromosome and one's start position
/// falls within the other's reference-allele span.
pub struct VariantGroupIterator<B: BufRead> {
    vcf_iters: Vec<GVcfRecordIterator<B>>,
    sorted_chromosomes: Vec<String>,
    current_chrom_idx: usize,
    chroms_seen: HashSet<String>,
    last_bin_start: Option<u32>,
    done: bool,
}

impl<B: BufRead> VariantGroupIterator<B> {
    /// Creates a new `VariantGroupIterator`.
    ///
    /// # Arguments
    /// * `vcf_iters` - Per-sample VCF iterators to merge.
    /// * `sorted_chromosomes` - The expected chromosome processing order.
    pub fn new(
        vcf_iters: Vec<GVcfRecordIterator<B>>,
        sorted_chromosomes: Vec<String>,
    ) -> VcfResult<Self> {
        let mut seen_samples = HashSet::new();
        for iter in &vcf_iters {
            for sample in iter.samples() {
                if !seen_samples.insert(sample) {
                    return Err(VcfParseError::RuntimeError {
                        message: format!("Sample '{}' appears in multiple VCF files", sample),
                    });
                }
            }
        }

        Ok(Self {
            vcf_iters,
            sorted_chromosomes,
            current_chrom_idx: 0,
            chroms_seen: HashSet::new(),
            last_bin_start: None,
            done: false,
        })
    }

    fn fail(&mut self, error: VcfParseError) -> Option<VcfResult<OverlappingVariantGroup>> {
        self.done = true;
        Some(Err(error))
    }

    fn get_next_variant_group(&mut self) -> Option<VcfResult<OverlappingVariantGroup>> {
        if self.done {
            return None;
        }

        // Phase A: find the next chromosome with variants and the earliest position
        let (current_chrom, bin_start, mut bin_end) = loop {
            if self.current_chrom_idx >= self.sorted_chromosomes.len() {
                self.done = true;
                return None;
            }

            let chrom = self.sorted_chromosomes[self.current_chrom_idx].clone();
            let mut earliest_pos: Option<u32> = None;
            let mut earliest_end: u32 = 0;

            for i in 0..self.vcf_iters.len() {
                let peek_info = match self.vcf_iters[i].peek_variant() {
                    Ok(Some(r)) => Some((r.chrom.clone(), r.pos, r.ref_allele_len)),
                    Ok(None) => None,
                    Err(e) => return self.fail(e),
                };

                let (p_chrom, p_pos, p_ref_len) = match peek_info {
                    Some(info) => info,
                    None => continue,
                };

                // Validate: chromosome not already fully processed
                if self.chroms_seen.contains(&p_chrom) {
                    return self.fail(VcfParseError::RuntimeError {
                        message: format!(
                            "A chromosome already seen has appeared: {}:{}, VCF seems not to be ordered",
                            p_chrom, p_pos
                        ),
                    });
                }

                if p_chrom != chrom {
                    continue;
                }

                // Validate: position ordering within chromosome
                if let Some(lbs) = self.last_bin_start {
                    if p_pos < lbs {
                        return self.fail(VcfParseError::RuntimeError {
                            message: format!(
                                "The VCF seems not to be ordered: {}:{}",
                                p_chrom, p_pos
                            ),
                        });
                    }
                }

                let p_end = p_pos + p_ref_len as u32 - 1;
                match earliest_pos {
                    None => {
                        earliest_pos = Some(p_pos);
                        earliest_end = p_end;
                    }
                    Some(ep) if p_pos < ep => {
                        earliest_pos = Some(p_pos);
                        earliest_end = p_end;
                    }
                    Some(ep) if p_pos == ep && p_end > earliest_end => {
                        earliest_end = p_end;
                    }
                    _ => {}
                }
            }

            match earliest_pos {
                Some(pos) => break (chrom, pos, earliest_end),
                None => {
                    self.chroms_seen.insert(chrom);
                    self.current_chrom_idx += 1;
                    self.last_bin_start = None;
                    continue;
                }
            }
        };

        // Phase B: build the bin by consuming all overlapping variants
        let mut records: Vec<GVcfRecord> = Vec::new();

        loop {
            let mut added_any = false;

            for i in 0..self.vcf_iters.len() {
                loop {
                    let peek_info = match self.vcf_iters[i].peek_variant() {
                        Ok(Some(r)) => Some((r.chrom.clone(), r.pos, r.ref_allele_len)),
                        Ok(None) => None,
                        Err(e) => return self.fail(e),
                    };

                    let (p_chrom, p_pos, p_ref_len) = match peek_info {
                        Some(info) => info,
                        None => break,
                    };

                    // Validate ordering
                    if self.chroms_seen.contains(&p_chrom) {
                        return self.fail(VcfParseError::RuntimeError {
                            message: format!(
                                "A chromosome already seen has appeared: {}:{}, VCF seems not to be ordered",
                                p_chrom, p_pos
                            ),
                        });
                    }
                    if p_chrom == current_chrom {
                        if let Some(lbs) = self.last_bin_start {
                            if p_pos < lbs {
                                return self.fail(VcfParseError::RuntimeError {
                                    message: format!(
                                        "The VCF seems not to be ordered: {}:{}",
                                        p_chrom, p_pos
                                    ),
                                });
                            }
                        }
                    }

                    // Check overlap: same chromosome and start within bin span
                    if p_chrom != current_chrom || p_pos > bin_end {
                        break;
                    }

                    // Consume the variant
                    let record = match self.vcf_iters[i].next() {
                        Some(Ok(r)) => r,
                        Some(Err(e)) => return self.fail(e),
                        None => unreachable!("peek returned Some but next returned None"),
                    };

                    let new_end = record.pos + p_ref_len as u32 - 1;
                    if new_end > bin_end {
                        bin_end = new_end;
                    }
                    records.push(record);
                    added_any = true;
                }
            }

            if !added_any {
                break;
            }
        }

        // Phase C: yield the bin and update state
        self.last_bin_start = Some(bin_start);

        // Check if current chromosome is exhausted across all iterators
        let mut any_remaining = false;
        for i in 0..self.vcf_iters.len() {
            match self.vcf_iters[i].peek_variant() {
                Ok(Some(r)) if r.chrom == current_chrom => {
                    any_remaining = true;
                    break;
                }
                Ok(_) => {}
                Err(e) => return self.fail(e),
            }
        }
        if !any_remaining {
            self.chroms_seen.insert(current_chrom.clone());
            self.current_chrom_idx += 1;
            self.last_bin_start = None;
        }

        Some(Ok(OverlappingVariantGroup {
            chrom: current_chrom,
            start: bin_start,
            end: bin_end,
            records,
        }))
    }
}

impl<B: BufRead> Iterator for VariantGroupIterator<B> {
    type Item = VcfResult<OverlappingVariantGroup>;

    fn next(&mut self) -> Option<Self::Item> {
        self.get_next_variant_group()
    }
}
