use std::io::BufReader;

use merge_per_sample_vcfs::gvcf_parser::VarIterator;
use merge_per_sample_vcfs::variant_grouping::VarGroupIterator;

fn make_iter(vcf_data: &str) -> VarIterator<BufReader<BufReader<&[u8]>>> {
    let reader = BufReader::new(vcf_data.as_bytes());
    VarIterator::from_reader(reader, 2).expect("Failed to create parser")
}

fn collect_spans(vcf_data: &[&str], sorted_chromosomes: Vec<String>) -> Vec<(String, u32, u32)> {
    let iters: Vec<_> = vcf_data.iter().map(|d| make_iter(d)).collect();
    let grouper =
        VarGroupIterator::new(iters, sorted_chromosomes).expect("Failed to create grouper");
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
    20\t1\t.\tGTATGG\tG\t20\tPASS\t.\tGT\t0|0";

// VCF4 — deletions on three different chromosomes
//       CGATGAT
// 20 1234567890
//        ------
const VCF4: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00004\n\
    1\t4\t.\tGATCGAT\tG\t20\tPASS\t.\tGT\t0|0\n\
    20\t4\t.\tTGGAGCAT\tT\t20\tPASS\t.\tGT\t0|0\n\
    21\t4\t.\tGATCGAT\tG\t20\tPASS\t.\tGT\t0|0";

// VCF5 — SNPs at positions 3, 8, and 12
const VCF5: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00005\n\
    20\t3\t.\tA\tT\t20\tPASS\t.\tGT\t0|0\n\
    20\t8\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    20\t12\t.\tG\tA\t20\tPASS\t.\tGT\t0|0";

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
    //    GTATGG
    // 20 123456
    // VCF1
    //    AAG AA
    // VCF2
    //    AAG AA

    let spans = collect_spans(&[VCF1, VCF2], vec!["1".into(), "20".into()]);
    assert_eq!(
        spans,
        vec![
            ("20".into(), 1, 1),
            ("20".into(), 2, 2),
            ("20".into(), 3, 3),
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
            ("20".into(), 5, 5),
            ("20".into(), 6, 6),
        ]
    );
}

#[test]
fn test_binning_with_deletion_spanning_several_snps() {
    //         GTATGG
    //      20 123456
    // VCF1    AAG.AA
    // VCF3    .-----
    // → all merge into one bin
    let spans = collect_spans(&[VCF1, VCF3], vec!["1".into(), "20".into()]);
    assert_eq!(spans, vec![("20".into(), 1, 6)]);

    // VCF1 + VCF3 + VCF4 + VCF5
    //            GATCGAT     GTATGGAGCATG     CTAGATCGAT
    //       1 1234567890  20 123456789012  21 1234567890
    // VCF1    ..........     AAG.AA......     ..........
    // VCF3    ..........     .-----......     ..........
    // VCF4    ....------     ....-------.     ....------
    // VCF5    ..........     ..T....A...A     ..........

    let spans = collect_spans(
        &[VCF1, VCF3, VCF4, VCF5],
        vec!["1".into(), "20".into(), "21".into()],
    );
    assert_eq!(
        spans,
        vec![
            ("1".into(), 4, 10),
            ("20".into(), 1, 11),
            ("20".into(), 12, 12),
            ("21".into(), 4, 10),
        ]
    );
}

#[test]
fn test_bin_contains_expected_variants() {
    // VCF1 (6 SNPs at 1-6) + VCF3 (deletion GTATGG→G at pos 1, spanning 1-6)
    //         GTATGG
    //      20 123456
    // VCF1    AAG.AA
    // VCF3    .-----
    // All variants merge into one bin: ("20", 1, 6)
    let iters = vec![make_iter(VCF1), make_iter(VCF3)];
    let grouper = VarGroupIterator::new(iters, vec!["1".into(), "20".into()]).unwrap();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();

    assert_eq!(groups.len(), 1);
    let group = &groups[0];
    assert_eq!(group.span(), ("20", 1, 6));

    // 6 SNPs from VCF1 + 1 deletion from VCF3 = 7 variants
    assert_eq!(group.variants.len(), 7);

    let positions: Vec<u32> = group.variants.iter().map(|v| v.pos).collect();
    assert_eq!(positions, vec![1, 2, 3, 4, 5, 6, 1]);

    let ref_alleles: Vec<&str> = group
        .variants
        .iter()
        .map(|v| v.alleles[0].as_str())
        .collect();
    assert_eq!(ref_alleles, vec!["G", "T", "A", "T", "G", "G", "GTATGG"]);

    // The deletion is the last variant (from VCF3, consumed after VCF1's variants)
    let deletion = &group.variants[6];
    assert_eq!(deletion.alleles[0], "GTATGG");
    assert_eq!(deletion.alleles[1], "G");
    assert_eq!(deletion.alleles[0].len(), 6);
}

#[test]
fn test_wrong_chrom_order() {
    let iter = make_iter(VCF_WITH_WRONG_CHROMOSOME_ORDER);
    let grouper = VarGroupIterator::new(vec![iter], vec!["1".into(), "2".into()]).unwrap();
    let results: Vec<_> = grouper.collect();
    assert!(
        results.iter().any(|r| r.is_err()),
        "Expected an error for wrong chromosome order"
    );
}

#[test]
fn test_wrong_order() {
    let iter = make_iter(VCF_WITH_WRONG_ORDER);
    let grouper = VarGroupIterator::new(vec![iter], vec!["1".into()]).unwrap();
    let results: Vec<_> = grouper.collect();
    assert!(
        results.iter().any(|r| r.is_err()),
        "Expected an error for wrong position order"
    );
}

#[test]
fn test_source_iter_samples() {
    // VCF1 has sample NA00001, VCF3 has sample NA00003
    // All variants merge into one bin: ("20", 1, 6)
    // First 6 variants come from VCF1 (iter 0), last from VCF3 (iter 1)
    let iters = vec![make_iter(VCF1), make_iter(VCF3)];
    let grouper = VarGroupIterator::new(iters, vec!["1".into(), "20".into()]).unwrap();
    let iter_info = grouper.iter_info().to_vec();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();

    assert_eq!(iter_info.len(), 2);
    assert_eq!(iter_info[0].samples, vec!["NA00001"]);
    assert_eq!(iter_info[1].samples, vec!["NA00003"]);

    let group = &groups[0];
    assert_eq!(group.variants.len(), 7);
    assert_eq!(group.source_var_iter_idxs.len(), 7);

    // First 6 variants from VCF1 (iter 0), last from VCF3 (iter 1)
    assert_eq!(group.source_var_iter_idxs, vec![0, 0, 0, 0, 0, 0, 1]);

    // Verify we can look up sample names for each variant via iter_info
    let sample_for_last_variant = &iter_info[group.source_var_iter_idxs[6]].samples;
    assert_eq!(sample_for_last_variant, &vec!["NA00003".to_string()]);
}

#[test]
fn test_duplicate_samples() {
    let iter1 = make_iter(VCF1);
    let iter2 = make_iter(VCF1); // same sample name NA00001
    let result = VarGroupIterator::new(vec![iter1, iter2], vec!["20".into()]);
    assert!(result.is_err(), "Expected error for duplicate sample names");
}

#[test]
fn test_empty_input() {
    let grouper: VarGroupIterator<BufReader<BufReader<&[u8]>>> =
        VarGroupIterator::new(vec![], vec!["1".into(), "20".into()]).unwrap();
    let groups: Vec<_> = grouper.collect();
    assert!(groups.is_empty());
}

#[test]
fn test_no_matching_chromosomes() {
    // All iterators have variants, but none on chromosomes in sorted_chromosomes
    let spans = collect_spans(&[VCF1, VCF2], vec!["X".into(), "Y".into()]);
    assert!(spans.is_empty());
}

#[test]
fn test_single_variant_group() {
    // A single SNP from a single iterator produces a single-variant group
    let vcf = "##\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
        1\t100\t.\tA\tT\t30\tPASS\t.\tGT\t0/1";
    let iters = vec![make_iter(vcf)];
    let grouper = VarGroupIterator::new(iters, vec!["1".into()]).unwrap();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();

    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].span(), ("1", 100, 100));
    assert_eq!(groups[0].variants.len(), 1);
    assert_eq!(groups[0].source_var_iter_idxs, vec![0]);
}

// VCF with only non-variable records (all hom-ref with single allele after <NON_REF> stripping)
const VCF_ALL_HOMREF: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00010\n\
    1\t1\t.\tA\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    1\t2\t.\tC\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    1\t3\t.\tG\t<NON_REF>\t20\tPASS\t.\tGT\t0|0";

const VCF_ALL_HOMREF2: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00011\n\
    1\t1\t.\tA\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    1\t2\t.\tC\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    1\t3\t.\tG\t<NON_REF>\t20\tPASS\t.\tGT\t0|0";

#[test]
fn test_non_variable_groups_are_skipped() {
    // When all iterators have non-variable records at the same position,
    // the skip path should consume them and produce no groups.
    let spans = collect_spans(&[VCF_ALL_HOMREF, VCF_ALL_HOMREF2], vec!["1".into()]);
    assert!(
        spans.is_empty(),
        "Expected no groups when all variants are non-variable, got: {:?}",
        spans
    );
}

#[test]
fn test_skip_path_produces_fewer_groups() {
    // VCF_ALL_HOMREF has 3 non-variable positions; VCF1 has 5 variable positions on chr 20.
    // Without skipping, positions from VCF_ALL_HOMREF would create groups on chr 1.
    // With skipping, only chr 20 groups appear.
    let spans = collect_spans(&[VCF_ALL_HOMREF, VCF1], vec!["1".into(), "20".into()]);
    // chr 1 positions should be skipped (all non-variable), only chr 20 groups remain
    assert!(
        spans.iter().all(|(chrom, _, _)| chrom == "20"),
        "Expected only chr 20 groups, got: {:?}",
        spans
    );
    assert_eq!(spans.len(), 5); // 5 variable positions from VCF1
}
