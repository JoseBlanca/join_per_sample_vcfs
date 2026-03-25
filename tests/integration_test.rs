use std::io::BufReader;

use join_per_sample_vcfs::gvcf_parser::VariantIterator;
use join_per_sample_vcfs::pipeline::merge_alleles_and_genotypes;

/// Convert a human-readable VCF string (spaces between fields) into proper
/// tab-delimited VCF.  For each line: trim leading/trailing whitespace, then
/// replace every run of 2+ spaces with a single tab.  Empty lines are dropped.
fn spaces_to_tabs(vcf: &str) -> String {
    vcf.lines()
        .map(|line| {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                return String::new();
            }
            let mut result = String::new();
            let mut space_count = 0usize;
            for ch in trimmed.chars() {
                if ch == ' ' {
                    space_count += 1;
                } else {
                    if space_count >= 2 {
                        result.push('\t');
                    } else if space_count == 1 {
                        result.push(' ');
                    }
                    space_count = 0;
                    result.push(ch);
                }
            }
            result
        })
        .filter(|l| !l.is_empty())
        .collect::<Vec<_>>()
        .join("\n")
}

fn make_iter(vcf_data: &str) -> VariantIterator<BufReader<BufReader<&[u8]>>> {
    let reader = BufReader::new(vcf_data.as_bytes());
    VariantIterator::from_reader(reader).expect("Failed to create parser")
}

/// Run the full pipeline and assert the output matches the expected VCF exactly.
///
/// Both the expected VCF and the actual output are normalized (trimmed, tabs)
/// before comparison, so constants can use spaces for readability.
fn assert_pipeline(gvcf_inputs: &[&str], chromosomes: Vec<String>, expected_vcf: &str) {
    let output = run_pipeline(gvcf_inputs, chromosomes);
    let got: Vec<&str> = output.lines().map(|l| l.trim()).collect();

    let normalized_expected = spaces_to_tabs(expected_vcf);
    let want: Vec<&str> = normalized_expected.lines().collect();

    assert_eq!(
        got.len(),
        want.len(),
        "Line count mismatch: got {} lines, expected {}\n--- got ---\n{}\n--- expected ---\n{}",
        got.len(),
        want.len(),
        got.join("\n"),
        want.join("\n"),
    );

    for (i, (g, w)) in got.iter().zip(want.iter()).enumerate() {
        assert_eq!(
            g,
            w,
            "Line {} differs:\n  got:      {}\n  expected: {}",
            i + 1,
            g,
            w
        );
    }
}

/// Run the full pipeline: parse gVCFs -> group -> merge -> write VCF.
/// Returns the output VCF as a string.
fn run_pipeline(gvcf_inputs: &[&str], chromosomes: Vec<String>) -> String {
    let inputs: Vec<String> = gvcf_inputs.iter().map(|s| spaces_to_tabs(s)).collect();
    let iters: Vec<_> = inputs.iter().map(|s| make_iter(s)).collect();

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

// ============================================================================
// Test data
// ============================================================================

// GATK-style gVCF for Sample1 covering chr1:1-5
// - Position 1: hom-ref
// - Position 2: hom-ref
// - Position 3: het SNP G->T
// - Position 4: hom-ref
// - Position 5: hom-ref
const GVCF_SAMPLE1: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample1
    chr1    1    .   A    <NON_REF>     .     .       .     GT      0|0
    chr1    2    .   C    <NON_REF>     .     .       .     GT      0|0
    chr1    3    .   G    T,<NON_REF>   30    PASS    .     GT      0|1
    chr1    4    .   T    <NON_REF>     .     .       .     GT      0|0
    chr1    5    .   A    <NON_REF>     .     .       .     GT      0|0
";

// GATK-style gVCF for Sample2 covering chr1:1-5
// - Position 1: hom-ref
// - Position 2: het SNP C->A
// - Position 3: hom-alt SNP G->C (different allele than Sample1)
// - Position 4: hom-ref
// - Position 5: het SNP A->G
const GVCF_SAMPLE2: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample2
    chr1    1    .   A    <NON_REF>     .     .       .     GT      0|0
    chr1    2    .   C    A,<NON_REF>   25    PASS    .     GT      0/1
    chr1    3    .   G    C,<NON_REF>   40    PASS    .     GT      1|1
    chr1    4    .   T    <NON_REF>     .     .       .     GT      0|0
    chr1    5    .   A    G,<NON_REF>   35    PASS    .     GT      0/1
";

const EXPECTED_SAMPLE1_SAMPLE2: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
    chr1    2    .   C    A    .     .       .     GT      0|0      0/1
    chr1    3    .   G    T,C  .     .       .     GT      0|1      2|2
    chr1    5    .   A    G    .     .       .     GT      0|0      0/1
";

#[test]
fn test_full_pipeline_two_samples() {
    assert_pipeline(
        &[GVCF_SAMPLE1, GVCF_SAMPLE2],
        vec!["chr1".to_string()],
        EXPECTED_SAMPLE1_SAMPLE2,
    );
}

// All positions hom-ref
const GVCF_ALL_REF: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT        QUAL  FILTER  INFO  FORMAT  Sample3
    chr1    1    .   A    <NON_REF>  .     .       .     GT      0|0
    chr1    2    .   C    <NON_REF>  .     .       .     GT      0|0
    chr1    3    .   G    <NON_REF>  .     .       .     GT      0|0
    chr1    4    .   T    <NON_REF>  .     .       .     GT      0|0
    chr1    5    .   A    <NON_REF>  .     .       .     GT      0|0
";

const EXPECTED_ALL_REF: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample3
";

#[test]
fn test_all_ref_samples_produce_no_variants() {
    assert_pipeline(&[GVCF_ALL_REF], vec!["chr1".to_string()], EXPECTED_ALL_REF);
}

// Deletion spanning positions 2-4
const GVCF_WITH_DELETION: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample4
    chr1    1    .   A    <NON_REF>     .     .       .     GT      0|0
    chr1    2    .   CGT  C,<NON_REF>   30    PASS    .     GT      0|1
    chr1    5    .   A    <NON_REF>     .     .       .     GT      0|0
";

// SNPs at each position overlapping with the deletion
const GVCF_WITH_SNPS: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample5
    chr1    1    .   A    <NON_REF>     .     .       .     GT      0|0
    chr1    2    .   C    A,<NON_REF>   20    PASS    .     GT      0|1
    chr1    3    .   G    T,<NON_REF>   20    PASS    .     GT      1|0
    chr1    4    .   T    C,<NON_REF>   20    PASS    .     GT      0|1
    chr1    5    .   A    <NON_REF>     .     .       .     GT      0|0
";

// Deletion at pos 2 (CGT->C) overlaps with SNPs at pos 2, 3, 4 in Sample5
// Merged ref = CGT (3 bases from deletion context)
// Sample4: both haplotypes produce "C" (ref first-char = "C", alt = "C")
// Sample5: hap1=C+T+T=CTT, hap2=A+G+C=AGC
const EXPECTED_DELETION_SNPS: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT        QUAL  FILTER  INFO  FORMAT  Sample4  Sample5
    chr1    2    .   CGT  C,CTT,AGC  .     .       .     GT      1|1      2|3
";

#[test]
fn test_deletion_overlapping_with_snps() {
    assert_pipeline(
        &[GVCF_WITH_DELETION, GVCF_WITH_SNPS],
        vec!["chr1".to_string()],
        EXPECTED_DELETION_SNPS,
    );
}

// One sample with a missing genotype allele
const GVCF_WITH_MISSING: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample6
    chr1    1    .   A    T,<NON_REF>   30    PASS    .     GT      0|.
    chr1    2    .   C    <NON_REF>     .     .       .     GT      0|0
    chr1    3    .   G    <NON_REF>     .     .       .     GT      0|0
";

const GVCF_SIMPLE_VARIANT: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample7
    chr1    1    .   A    G,<NON_REF>   40    PASS    .     GT      1|1
    chr1    2    .   C    <NON_REF>     .     .       .     GT      0|0
    chr1    3    .   G    <NON_REF>     .     .       .     GT      0|0
";

// Sample6 has 0|. (missing haplotype), so T is not contributed as an alt allele.
// Only Sample7's G survives as alt. Sample7: 1|1 -> G|G.
const EXPECTED_MISSING: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample6  Sample7
    chr1    1    .   A    G    .     .       .     GT      0|.      1|1
";

#[test]
fn test_missing_genotype_propagation() {
    assert_pipeline(
        &[GVCF_WITH_MISSING, GVCF_SIMPLE_VARIANT],
        vec!["chr1".to_string()],
        EXPECTED_MISSING,
    );
}

// Two samples across two chromosomes
const GVCF_MULTI_CHROM_S1: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample8
    chr1    1    .   A    T,<NON_REF>   30    PASS    .     GT      0|1
    chr1    2    .   C    <NON_REF>     .     .       .     GT      0|0
    chr2    1    .   G    A,<NON_REF>   25    PASS    .     GT      1|0
    chr2    2    .   T    <NON_REF>     .     .       .     GT      0|0
";

const GVCF_MULTI_CHROM_S2: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample9
    chr1    1    .   A    <NON_REF>     .     .       .     GT      0|0
    chr1    2    .   C    G,<NON_REF>   20    PASS    .     GT      1|0
    chr2    1    .   G    <NON_REF>     .     .       .     GT      0|0
    chr2    2    .   T    C,<NON_REF>   35    PASS    .     GT      0|1
";

const EXPECTED_MULTI_CHROM: &str = "
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample8  Sample9
    chr1    1    .   A    T    .     .       .     GT      0|1      0|0
    chr1    2    .   C    G    .     .       .     GT      0|0      1|0
    chr2    1    .   G    A    .     .       .     GT      1|0      0|0
    chr2    2    .   T    C    .     .       .     GT      0|0      0|1
";

#[test]
fn test_multi_chromosome_pipeline() {
    assert_pipeline(
        &[GVCF_MULTI_CHROM_S1, GVCF_MULTI_CHROM_S2],
        vec!["chr1".to_string(), "chr2".to_string()],
        EXPECTED_MULTI_CHROM,
    );
}

// Sample1: het SNP at pos 1, hom-ref at pos 2
const GVCF_SIMPLE_SNPS_S1: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample1
chr1    1    .   A    T,<NON_REF>   .     .       .     GT      0|1
chr1    2    .   A    <NON_REF>     .     .       .     GT      0|0
";

// Sample2: hom-alt SNP at pos 1 (different allele), hom-ref at pos 2
const GVCF_SIMPLE_SNPS_S2: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample2
chr1    1    .   A    C,<NON_REF>   .     .       .     GT      1/1
chr1    2    .   A    <NON_REF>     25    .       .     GT      0/0
";

// Expected: pos 1 merges T and C as alt alleles; pos 2 filtered (both hom-ref)
const EXPECTED_SIMPLE_SNPS: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
chr1    1    .   A    T,C  .     .       .     GT      0|1      2/2
";

#[test]
fn test_template_two_samples_snps() {
    assert_pipeline(
        &[GVCF_SIMPLE_SNPS_S1, GVCF_SIMPLE_SNPS_S2],
        vec!["chr1".to_string()],
        EXPECTED_SIMPLE_SNPS,
    );
}

// Two non reference SNPs
const GVCF_SIMPLE_SNP_NON_REF_S1: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample1
chr1    1    .   A    T,<NON_REF>   .     .       .     GT      1|1
chr1    2    .   A    <NON_REF>     .     .       .     GT      0|0
";
const GVCF_SIMPLE_SNP_NON_REF_S2: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample2
chr1    1    .   A    C,<NON_REF>   .     .       .     GT      1/1
chr1    2    .   A    <NON_REF>     25    .       .     GT      0/0
";
const EXPECTED_SIMPLE_SNP_NON_REF: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
chr1    1    .   A    T,C  .     .       .     GT      1|1      2/2
";

#[test]
fn test_template_two_samples_non_ref_snps() {
    assert_pipeline(
        &[GVCF_SIMPLE_SNP_NON_REF_S1, GVCF_SIMPLE_SNP_NON_REF_S2],
        vec!["chr1".to_string()],
        EXPECTED_SIMPLE_SNP_NON_REF,
    );
}

// Deletion in sample 1 and snp in sample 2
const GVCF3_S1: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample1
chr1    1    .   AA   A,<NON_REF>   .     .       .     GT      0|1
chr1    2    .   A    <NON_REF>     .     .       .     GT      0|0
";
const GVCF3_S2: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample2
chr1    1    .   A    C,<NON_REF>   .     .       .     GT      1/1
chr1    2    .   A    <NON_REF>     25    .       .     GT      0/0
";
const EXPECTED3: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
chr1    1    .   AA   A,CA  .     .       .     GT      0|1      2/2
";

#[test]
fn test_del1() {
    assert_pipeline(&[GVCF3_S1, GVCF3_S2], vec!["chr1".to_string()], EXPECTED3);
}

// Deletion in sample 1 and snp in sample 2
const GVCF4_S1: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample1
chr1    1    .   AA   A,<NON_REF>   .     .       .     GT      0|1
chr1    2    .   A    <NON_REF>     .     .       .     GT      0|0
chr1    3    .   A    <NON_REF>     .     .       .     GT      0|0
";
const GVCF4_S2: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT           QUAL  FILTER  INFO  FORMAT  Sample2
chr1    1    .   A    C,<NON_REF>   .     .       .     GT      1/1
chr1    2    .   A    T,<NON_REF>   .     .       .     GT      1/1
chr1    3    .   A    G,<NON_REF>   .     .       .     GT      1|0
";
const EXPECTED4: &str = "##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
chr1    1    .   AA   A,CT  .     .       .     GT      0|1      2/2
chr1    3    .   A    G     .     .       .     GT      0|0      1|0
";

#[test]
fn test_del2() {
    assert_pipeline(&[GVCF4_S1, GVCF4_S2], vec!["chr1".to_string()], EXPECTED4);
}
