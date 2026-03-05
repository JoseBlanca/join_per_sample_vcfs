use std::collections::HashSet;
use std::io::BufRead;

use crate::errors::VcfParseError;
use crate::gvcf_parser::{Variant, VariantIterator, VcfResult};

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
    pub variants: Vec<Variant>,
}

impl OverlappingVariantGroup {
    pub fn span(&self) -> (&str, u32, u32) {
        (&self.chrom, self.start, self.end)
    }
}

/// Computes the reference-allele span end for a variant (1-based, inclusive).
#[inline]
fn ref_span_end(record: &Variant) -> u32 {
    record.pos + record.ref_allele_len as u32 - 1
}

/// Iterator that merges multiple per-sample VCF streams, grouping
/// overlapping variants into bins based on reference-allele span overlap.
///
/// Variants from all input VCFs are consumed in genomic order. Two variants
/// overlap when they are on the same chromosome and one's start position
/// falls within the other's reference-allele span.
pub struct VariantGroupIterator<B: BufRead> {
    vcf_iters: Vec<VariantIterator<B>>,
    sorted_chromosomes: Vec<String>,
    current_chrom_idx: usize,
    chroms_seen: HashSet<String>,
    last_group_start: Option<u32>,
    done: bool,
}

impl<B: BufRead> VariantGroupIterator<B> {
    /// Creates a new `VariantGroupIterator`.
    ///
    /// # Arguments
    /// * `vcf_iters` - Per-sample VCF iterators to merge.
    /// * `sorted_chromosomes` - The expected chromosome processing order.
    pub fn new(
        vcf_iters: Vec<VariantIterator<B>>,
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
            last_group_start: None,
            done: false,
        })
    }

    fn fail(&mut self, error: VcfParseError) -> Option<VcfResult<OverlappingVariantGroup>> {
        self.done = true;
        Some(Err(error))
    }

    fn validate_ordering(&self, p_chrom: &str, p_pos: u32, expected_chrom: &str) -> VcfResult<()> {
        if self.chroms_seen.contains(p_chrom) {
            return Err(VcfParseError::RuntimeError {
                message: format!(
                    "Out-of-order chromosome '{}' at position {}: chromosome already processed",
                    p_chrom, p_pos
                ),
            });
        }

        if p_chrom == expected_chrom {
            if let Some(lbs) = self.last_group_start {
                if p_pos < lbs {
                    return Err(VcfParseError::RuntimeError {
                        message: format!(
                            "Out-of-order position on '{}': {} < last bin start {}",
                            p_chrom, p_pos, lbs
                        ),
                    });
                }
            }
        }

        Ok(())
    }

    /// Phase A: find the next chromosome span seed:
    /// (chromosome, earliest start, current end).
    fn compute_next_span_seed(&mut self) -> VcfResult<Option<(String, u32, u32)>> {
        loop {
            if self.current_chrom_idx >= self.sorted_chromosomes.len() {
                self.done = true;
                return Ok(None);
            }

            let chrom = self.sorted_chromosomes[self.current_chrom_idx].clone();
            let mut earliest_pos: Option<u32> = None;
            let mut earliest_end: u32 = 0;

            for i in 0..self.vcf_iters.len() {
                let (peek_chrom, peek_pos, p_end) = match self.vcf_iters[i].peek_variant()? {
                    Some(r) => (r.chrom.clone(), r.pos, ref_span_end(r)),
                    None => continue,
                };

                self.validate_ordering(&peek_chrom, peek_pos, &chrom)?;

                if peek_chrom != chrom {
                    continue;
                }

                match earliest_pos {
                    None => {
                        earliest_pos = Some(peek_pos);
                        earliest_end = p_end;
                    }
                    Some(ep) if peek_pos < ep => {
                        earliest_pos = Some(peek_pos);
                        earliest_end = p_end;
                    }
                    Some(ep) if peek_pos == ep && p_end > earliest_end => {
                        earliest_end = p_end;
                    }
                    _ => {}
                }
            }

            match earliest_pos {
                Some(pos) => return Ok(Some((chrom, pos, earliest_end))),
                None => {
                    self.chroms_seen.insert(chrom);
                    self.current_chrom_idx += 1;
                    self.last_group_start = None;
                }
            }
        }
    }

    /// Phase B: consume all records overlapping the current bin span.
    fn group_overlapping_vars(
        &mut self,
        current_chrom: &str,
        mut group_end: u32,
    ) -> VcfResult<(Vec<Variant>, u32)> {
        let mut records: Vec<Variant> = Vec::new();

        loop {
            let mut added_any = false;

            for i in 0..self.vcf_iters.len() {
                loop {
                    let (peek_chrom, peek_pos) = match self.vcf_iters[i].peek_variant()? {
                        Some(r) => (r.chrom.clone(), r.pos),
                        None => break,
                    };

                    self.validate_ordering(&peek_chrom, peek_pos, current_chrom)?;

                    // overlap rule: same chromosome and start <= group_end
                    if peek_chrom != current_chrom || peek_pos > group_end {
                        break;
                    }

                    let record =
                        self.vcf_iters[i]
                            .next()
                            .ok_or_else(|| VcfParseError::RuntimeError {
                                message:
                                    "Iterator ended unexpectedly while consuming overlapping record"
                                        .to_string(),
                            })??;

                    let new_end = ref_span_end(&record);
                    if new_end > group_end {
                        group_end = new_end;
                    }

                    records.push(record);
                    added_any = true;
                }
            }

            if !added_any {
                break;
            }
        }

        Ok((records, group_end))
    }

    fn chromosome_has_remaining_variants(&mut self, current_chrom: &str) -> VcfResult<bool> {
        for i in 0..self.vcf_iters.len() {
            match self.vcf_iters[i].peek_variant()? {
                Some(r) if r.chrom == current_chrom => return Ok(true),
                _ => {}
            }
        }
        Ok(false)
    }

    fn get_next_variant_group(&mut self) -> Option<VcfResult<OverlappingVariantGroup>> {
        if self.done {
            return None;
        }

        let (current_chrom, group_start, group_end) = match self.compute_next_span_seed() {
            Ok(Some(v)) => v,
            Ok(None) => return None,
            Err(e) => return self.fail(e),
        };

        let (overlapping_vars, group_end) =
            match self.group_overlapping_vars(&current_chrom, group_end) {
                Ok(v) => v,
                Err(e) => return self.fail(e),
            };

        self.last_group_start = Some(group_start);

        let any_remaining = match self.chromosome_has_remaining_variants(&current_chrom) {
            Ok(v) => v,
            Err(e) => return self.fail(e),
        };

        if !any_remaining {
            self.chroms_seen.insert(current_chrom.clone());
            self.current_chrom_idx += 1;
            self.last_group_start = None;
        }

        Some(Ok(OverlappingVariantGroup {
            chrom: current_chrom,
            start: group_start,
            end: group_end,
            variants: overlapping_vars,
        }))
    }
}

impl<B: BufRead> Iterator for VariantGroupIterator<B> {
    type Item = VcfResult<OverlappingVariantGroup>;

    fn next(&mut self) -> Option<Self::Item> {
        self.get_next_variant_group()
    }
}
