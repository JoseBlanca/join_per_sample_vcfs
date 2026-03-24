use std::io::BufReader;

use join_per_sample_vcfs::gvcf_parser::VariantIterator;
use join_per_sample_vcfs::pipeline::merge_alleles_and_genotypes;

fn make_iter(vcf_data: &str) -> VariantIterator<BufReader<BufReader<&[u8]>>> {
    let reader = BufReader::new(vcf_data.as_bytes());
    VariantIterator::from_reader(reader).expect("Failed to create parser")
}

/// Run the full pipeline: parse gVCFs -> group -> merge -> write VCF.
/// Returns the output VCF as a string.
fn run_pipeline(gvcf_inputs: &[&str], chromosomes: Vec<String>) -> String {
    let iters: Vec<_> = gvcf_inputs.iter().map(|input| make_iter(input)).collect();

    let shared_buf = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let shared_clone = shared_buf.clone();

    struct SharedWriter(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
    impl std::io::Write for SharedWriter {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    merge_alleles_and_genotypes(iters, chromosomes, Box::new(SharedWriter(shared_clone)))
        .expect("Pipeline failed");

    let buf = shared_buf.lock().unwrap();
    String::from_utf8(buf.clone()).unwrap()
}

// GATK-style gVCF for Sample1 covering chr1:1-5
// - Position 1: ref (hom-ref with <NON_REF>)
// - Position 2: ref (hom-ref with <NON_REF>)
// - Position 3: het SNP A->T (with <NON_REF>)
// - Position 4: ref (hom-ref with <NON_REF>)
// - Position 5: ref (hom-ref with <NON_REF>)
const GVCF_SAMPLE1: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n\
    chr1\t1\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t2\t.\tC\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t3\t.\tG\tT,<NON_REF>\t30\tPASS\t.\tGT\t0|1\n\
    chr1\t4\t.\tT\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t5\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0";

// GATK-style gVCF for Sample2 covering chr1:1-5
// - Position 1: ref (hom-ref with <NON_REF>)
// - Position 2: het SNP C->A (with <NON_REF>)
// - Position 3: hom-alt SNP G->C (with <NON_REF>)
// - Position 4: ref (hom-ref with <NON_REF>)
// - Position 5: het SNP A->G (with <NON_REF>)
const GVCF_SAMPLE2: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample2\n\
    chr1\t1\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t2\t.\tC\tA,<NON_REF>\t25\tPASS\t.\tGT\t0/1\n\
    chr1\t3\t.\tG\tC,<NON_REF>\t40\tPASS\t.\tGT\t1|1\n\
    chr1\t4\t.\tT\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t5\t.\tA\tG,<NON_REF>\t35\tPASS\t.\tGT\t0/1";

#[test]
fn test_full_pipeline_two_samples() {
    let output = run_pipeline(
        &[GVCF_SAMPLE1, GVCF_SAMPLE2],
        vec!["chr1".to_string()],
    );

    let lines: Vec<&str> = output.lines().collect();

    // Header
    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert_eq!(
        lines[1],
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2"
    );

    // Only variable positions should be emitted (positions 2, 3, 5).
    // Position 1 and 4 are hom-ref in both samples => filtered as non-variable.

    // Position 2: Sample1=0|0 (ref), Sample2=0/1 (het C->A)
    assert_eq!(
        lines[2],
        "chr1\t2\t.\tC\tA\t.\t.\t.\tGT\t0|0\t0/1"
    );

    // Position 3: Sample1=0|1 (het G->T), Sample2=1|1 (hom-alt G->C)
    // Both have different ALT alleles, so merged alleles: G (ref), T, C
    assert_eq!(
        lines[3],
        "chr1\t3\t.\tG\tT,C\t.\t.\t.\tGT\t0|1\t2|2"
    );

    // Position 5: Sample1=0|0 (ref), Sample2=0/1 (het A->G)
    assert_eq!(
        lines[4],
        "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0|0\t0/1"
    );

    // Exactly 5 lines total: 2 header + 3 data
    assert_eq!(lines.len(), 5);
}

// GATK-style gVCF for Sample3 covering chr1:1-5
// All positions are hom-ref
const GVCF_ALL_REF: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample3\n\
    chr1\t1\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t2\t.\tC\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t3\t.\tG\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t4\t.\tT\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t5\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0";

#[test]
fn test_all_ref_samples_produce_no_variants() {
    let output = run_pipeline(
        &[GVCF_ALL_REF],
        vec!["chr1".to_string()],
    );

    let lines: Vec<&str> = output.lines().collect();
    // Only header lines, no data lines (everything is non-variable)
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0], "##fileformat=VCFv4.2");
}

// gVCF with a deletion spanning positions 2-4
const GVCF_WITH_DELETION: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample4\n\
    chr1\t1\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t2\t.\tCGT\tC,<NON_REF>\t30\tPASS\t.\tGT\t0|1\n\
    chr1\t5\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0";

// gVCF covering the same region with SNPs at each position
const GVCF_WITH_SNPS: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample5\n\
    chr1\t1\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t2\t.\tC\tA,<NON_REF>\t20\tPASS\t.\tGT\t0|1\n\
    chr1\t3\t.\tG\tT,<NON_REF>\t20\tPASS\t.\tGT\t1|0\n\
    chr1\t4\t.\tT\tC,<NON_REF>\t20\tPASS\t.\tGT\t0|1\n\
    chr1\t5\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0";

#[test]
fn test_deletion_overlapping_with_snps() {
    let output = run_pipeline(
        &[GVCF_WITH_DELETION, GVCF_WITH_SNPS],
        vec!["chr1".to_string()],
    );

    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert_eq!(
        lines[1],
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample4\tSample5"
    );

    // The deletion at pos 2 (CGT->C) overlaps with SNPs at pos 2, 3, 4 in Sample5
    // This creates a merged group spanning positions 2-4.
    // Ref allele is CGT (3 bases from the deletion context).
    // Sample4 haplotypes: 0|1 for CGT->C, so alleles: CGT (ref) and C (del)
    // Sample5 haplotypes at 2,3,4: (C,T,C) and (A,G,T) per haplotype
    //   hap1: 0|1, 1|0, 0|1 => C+G+T = CGT (ref), A+T+C = ATC
    //   wait — let me think through this more carefully.
    // Actually, Sample5 at pos2: 0|1 means hap1=C(ref), hap2=A
    //   at pos3: 1|0 means hap1=T, hap2=G(ref)
    //   at pos4: 0|1 means hap1=T(ref), hap2=C
    // So Sample5 hap1 = C+T+T = CTT, hap2 = A+G+C = AGC
    // Merged alleles: CGT(ref=0), C(del=1), CTT(=2), AGC(=3)
    // Sample4: 0|1 => CGT|C
    // Sample5: 2|3 => CTT|AGC

    // There should be data lines for the merged group
    assert!(lines.len() > 2, "Expected variant output lines");

    // Verify the merged variant at pos 2
    let fields: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(fields[0], "chr1"); // CHROM
    assert_eq!(fields[1], "2"); // POS
    assert_eq!(fields[3], "CGT"); // REF (3 bases from deletion context)
    // ALT should include the deletion allele and the SNP-derived haplotype alleles
    assert!(fields[4].contains("C"), "ALT should contain deletion allele 'C'");
}

// Two samples, one with a missing genotype
const GVCF_WITH_MISSING: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample6\n\
    chr1\t1\t.\tA\tT,<NON_REF>\t30\tPASS\t.\tGT\t0|.\n\
    chr1\t2\t.\tC\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t3\t.\tG\t<NON_REF>\t.\t.\t.\tGT\t0|0";

const GVCF_SIMPLE_VARIANT: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample7\n\
    chr1\t1\t.\tA\tG,<NON_REF>\t40\tPASS\t.\tGT\t1|1\n\
    chr1\t2\t.\tC\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t3\t.\tG\t<NON_REF>\t.\t.\t.\tGT\t0|0";

#[test]
fn test_missing_genotype_propagation() {
    let output = run_pipeline(
        &[GVCF_WITH_MISSING, GVCF_SIMPLE_VARIANT],
        vec!["chr1".to_string()],
    );

    let lines: Vec<&str> = output.lines().collect();

    // Position 1 has variants in both samples
    let fields: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(fields[0], "chr1");
    assert_eq!(fields[1], "1");
    assert_eq!(fields[3], "A"); // REF

    // Sample6 has a missing allele (0|.), Sample7 is 1|1
    // The missing allele should be preserved as '.'
    let sample6_gt = fields[9];
    assert!(
        sample6_gt.contains('.'),
        "Sample6 genotype should contain missing allele '.', got: {}",
        sample6_gt
    );
    let sample7_gt = fields[10];
    assert!(
        !sample7_gt.contains('.'),
        "Sample7 genotype should not have missing alleles, got: {}",
        sample7_gt
    );
}

// Three samples across two chromosomes
const GVCF_MULTI_CHROM_S1: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample8\n\
    chr1\t1\t.\tA\tT,<NON_REF>\t30\tPASS\t.\tGT\t0|1\n\
    chr1\t2\t.\tC\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr2\t1\t.\tG\tA,<NON_REF>\t25\tPASS\t.\tGT\t1|0\n\
    chr2\t2\t.\tT\t<NON_REF>\t.\t.\t.\tGT\t0|0";

const GVCF_MULTI_CHROM_S2: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample9\n\
    chr1\t1\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr1\t2\t.\tC\tG,<NON_REF>\t20\tPASS\t.\tGT\t1|0\n\
    chr2\t1\t.\tG\t<NON_REF>\t.\t.\t.\tGT\t0|0\n\
    chr2\t2\t.\tT\tC,<NON_REF>\t35\tPASS\t.\tGT\t0|1";

#[test]
fn test_multi_chromosome_pipeline() {
    let output = run_pipeline(
        &[GVCF_MULTI_CHROM_S1, GVCF_MULTI_CHROM_S2],
        vec!["chr1".to_string(), "chr2".to_string()],
    );

    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(
        lines[1],
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample8\tSample9"
    );

    // Data lines: chr1:1, chr1:2, chr2:1, chr2:2 (all have at least one non-ref sample)
    // Collect (chrom, pos) from data lines
    let data_positions: Vec<(&str, &str)> = lines[2..]
        .iter()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (f[0], f[1])
        })
        .collect();

    // Chromosomes should appear in order: chr1 first, then chr2
    for (chrom, _) in &data_positions {
        assert!(
            *chrom == "chr1" || *chrom == "chr2",
            "Unexpected chromosome: {}",
            chrom
        );
    }
    let chr1_count = data_positions.iter().filter(|(c, _)| *c == "chr1").count();
    let chr2_count = data_positions.iter().filter(|(c, _)| *c == "chr2").count();
    assert!(chr1_count > 0, "Should have chr1 variants");
    assert!(chr2_count > 0, "Should have chr2 variants");

    // chr1 variants should come before chr2 variants
    let first_chr2_idx = data_positions
        .iter()
        .position(|(c, _)| *c == "chr2")
        .unwrap();
    let last_chr1_idx = data_positions
        .iter()
        .rposition(|(c, _)| *c == "chr1")
        .unwrap();
    assert!(last_chr1_idx < first_chr2_idx, "chr1 variants should precede chr2");

    // Check specific genotypes
    // chr1:1 — Sample8=0|1 (het A->T), Sample9=0|0 (ref)
    assert_eq!(lines[2], "chr1\t1\t.\tA\tT\t.\t.\t.\tGT\t0|1\t0|0");

    // chr1:2 — Sample8=0|0 (ref), Sample9=1|0 (het C->G)
    assert_eq!(lines[3], "chr1\t2\t.\tC\tG\t.\t.\t.\tGT\t0|0\t1|0");

    // chr2:1 — Sample8=1|0 (het G->A), Sample9=0|0 (ref)
    assert_eq!(lines[4], "chr2\t1\t.\tG\tA\t.\t.\t.\tGT\t1|0\t0|0");

    // chr2:2 — Sample8=0|0 (ref), Sample9=0|1 (het T->C)
    assert_eq!(lines[5], "chr2\t2\t.\tT\tC\t.\t.\t.\tGT\t0|0\t0|1");

    assert_eq!(lines.len(), 6); // 2 header + 4 data
}
