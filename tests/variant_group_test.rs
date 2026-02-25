use std::io::BufReader;

use join_per_sample_vcfs::gvcf_parser::GVcfRecordIterator;
use join_per_sample_vcfs::variant_group::VariantGroupIterator;

fn make_iter(vcf_data: &str) -> GVcfRecordIterator<BufReader<BufReader<&[u8]>>> {
    let reader = BufReader::new(vcf_data.as_bytes());
    GVcfRecordIterator::from_reader(reader).expect("Failed to create parser")
}

fn collect_spans(
    vcf_data: &[&str],
    sorted_chromosomes: Vec<String>,
) -> Vec<(String, u32, u32)> {
    let iters: Vec<_> = vcf_data.iter().map(|d| make_iter(d)).collect();
    let grouper =
        VariantGroupIterator::new(iters, sorted_chromosomes).expect("Failed to create grouper");
    grouper
        .map(|r| {
            let bin = r.expect("Unexpected error");
            (bin.chrom, bin.start, bin.end)
        })
        .collect()
}

// VCF1
//    GTATGG
// 20 123456
//    AAGTAA
const VCF1: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n\
    20\t1\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t.|0\n\
    20\t3\t.\tA\tG\t20\tPASS\t.\tGT\t.|0\n\
    20\t4\t.\tT\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    20\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0/1\n\
    20\t6\t.\tG\tA\t20\tPASS\t.\tGT\t0/1";

// VCF2
//    GTATGG
// 20 123456
//    AAGTAA
const VCF2: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00002\n\
    20\t1\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t.|0\n\
    20\t3\t.\tA\tG\t20\tPASS\t.\tGT\t.|0\n\
    20\t4\t.\tT\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    20\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0/1\n\
    20\t6\t.\tG\tA\t20\tPASS\t.\tGT\t0/1";

// VCF3 — deletion spanning positions 1-7
//    GATCGAT
// 20 1234567
//     ------
const VCF3: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00003\n\
    20\t1\t.\tGATCGAT\tG\t20\tPASS\t.\tGT\t0|0";

// VCF4 — deletions on three different chromosomes
//       CGATGAT
// 20 1234567890
//        ------
const VCF4: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00004\n\
    1\t4\t.\tGATCGAT\tG\t20\tPASS\t.\tGT\t0|0\n\
    20\t4\t.\tCGATGAT\tC\t20\tPASS\t.\tGT\t0|0\n\
    21\t4\t.\tGATCGAT\tG\t20\tPASS\t.\tGT\t0|0";

// VCF5 — SNPs at positions 3, 8, and 20
const VCF5: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00005\n\
    20\t3\t.\tT\tA\t20\tPASS\t.\tGT\t0|0\n\
    20\t8\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    20\t20\t.\tG\tA\t20\tPASS\t.\tGT\t0|0";

const VCF_WITH_WRONG_CHROMOSOME_ORDER: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00005\n\
    1\t3\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    2\t8\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    1\t20\t.\tG\tA\t20\tPASS\t.\tGT\t0|0";

const VCF_WITH_WRONG_ORDER: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00005\n\
    1\t3\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    1\t4\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    1\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    1\t8\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    1\t2\t.\tG\tA\t20\tPASS\t.\tGT\t0|0";

#[test]
fn test_simple_binning() {
    // VCF1 + VCF2: same SNP positions → one bin per position
    let spans = collect_spans(&[VCF1, VCF2], vec!["1".into(), "20".into()]);
    assert_eq!(
        spans,
        vec![
            ("20".into(), 1, 1),
            ("20".into(), 2, 2),
            ("20".into(), 3, 3),
            ("20".into(), 4, 4),
            ("20".into(), 5, 5),
            ("20".into(), 6, 6),
        ]
    );

    // Single VCF
    let spans = collect_spans(&[VCF1], vec!["20".into()]);
    assert_eq!(
        spans,
        vec![
            ("20".into(), 1, 1),
            ("20".into(), 2, 2),
            ("20".into(), 3, 3),
            ("20".into(), 4, 4),
            ("20".into(), 5, 5),
            ("20".into(), 6, 6),
        ]
    );
}

#[test]
fn test_binning_with_deletion_spanning_several_snps() {
    // VCF1 (6 SNPs at 1-6) + VCF3 (deletion at 1 spanning 1-7)
    // → all merge into one bin
    let spans = collect_spans(&[VCF1, VCF3], vec!["1".into(), "20".into()]);
    assert_eq!(spans, vec![("20".into(), 1, 7)]);

    // VCF1 + VCF3 + VCF4 + VCF5
    // chr 1: VCF4 has deletion at 4-10 → bin (1, 4, 10)
    // chr 20: VCF1 SNPs 1-6, VCF3 deletion 1-7, VCF4 deletion 4-10, VCF5 SNPs 3,8
    //   all overlap transitively → bin (20, 1, 10)
    //   VCF5 SNP at 20 is separate → bin (20, 20, 20)
    // chr 21: VCF4 has deletion at 4-10 → bin (21, 4, 10)
    let spans = collect_spans(
        &[VCF1, VCF3, VCF4, VCF5],
        vec!["1".into(), "20".into(), "21".into()],
    );
    assert_eq!(
        spans,
        vec![
            ("1".into(), 4, 10),
            ("20".into(), 1, 10),
            ("20".into(), 20, 20),
            ("21".into(), 4, 10),
        ]
    );
}

#[test]
fn test_wrong_chrom_order() {
    let iter = make_iter(VCF_WITH_WRONG_CHROMOSOME_ORDER);
    let grouper =
        VariantGroupIterator::new(vec![iter], vec!["1".into(), "2".into()]).unwrap();
    let results: Vec<_> = grouper.collect();
    assert!(
        results.iter().any(|r| r.is_err()),
        "Expected an error for wrong chromosome order"
    );
}

#[test]
fn test_wrong_order() {
    let iter = make_iter(VCF_WITH_WRONG_ORDER);
    let grouper = VariantGroupIterator::new(vec![iter], vec!["1".into()]).unwrap();
    let results: Vec<_> = grouper.collect();
    assert!(
        results.iter().any(|r| r.is_err()),
        "Expected an error for wrong position order"
    );
}

#[test]
fn test_duplicate_samples() {
    let iter1 = make_iter(VCF1);
    let iter2 = make_iter(VCF1); // same sample name NA00001
    let result = VariantGroupIterator::new(vec![iter1, iter2], vec!["20".into()]);
    assert!(result.is_err(), "Expected error for duplicate sample names");
}
