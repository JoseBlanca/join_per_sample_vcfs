use std::collections::HashSet;
use std::io::BufReader;

use merge_per_sample_vcfs::genotype_merging::{merge_variant_group, merge_vars_in_groups};
use merge_per_sample_vcfs::genotype_posteriors::PriorConfig;
use merge_per_sample_vcfs::gvcf_parser::{Variant, VarIterator};
use merge_per_sample_vcfs::variant_grouping::VariantGroupIterator;
use merge_per_sample_vcfs::variant_grouping::{OverlappingVariantGroup, VariantIteratorInfo};

fn default_prior() -> PriorConfig {
    PriorConfig::default()
}

fn make_iter(vcf_data: &str) -> VarIterator<BufReader<BufReader<&[u8]>>> {
    let reader = BufReader::new(vcf_data.as_bytes());
    VarIterator::from_reader(reader).expect("Failed to create parser")
}

/// Collect merged variants from groups into a Vec (test helper).
fn collect_merged_variants<I>(groups: I, iter_info: &[VariantIteratorInfo]) -> Vec<Variant>
where
    I: Iterator<Item = merge_per_sample_vcfs::gvcf_parser::VcfResult<OverlappingVariantGroup>>
        + Send,
{
    let mut vars = Vec::new();
    merge_vars_in_groups(groups, iter_info, &default_prior(), |result| {
        vars.push(result?);
        Ok(())
    })
    .unwrap();
    vars
}

const VCF1: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n\
    20\t1\t.\tG\tA\t20\tPASS\t.\tGT\t1|1\n\
    20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t.|0\n\
    20\t3\t.\tA\tG\t20\tPASS\t.\tGT\t.|0\n\
    20\t4\t.\tT\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    20\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0/1\n\
    20\t6\t.\tG\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t7\t.\tA\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t8\t.\tG\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t9\t.\tC\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t10\t.\tA\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t11\t.\tT\t.\t20\tPASS\t.\tGT\t0/0";

const VCF2: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00002\n\
    20\t1\t.\tG\tA\t20\tPASS\t.\tGT\t1|1\n\
    20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t.|0\n\
    20\t3\t.\tA\tG\t20\tPASS\t.\tGT\t.|0\n\
    20\t4\t.\tT\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    20\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0/1\n\
    20\t6\t.\tG\tA\t20\tPASS\t.\tGT\t0/1\n\
    20\t7\t.\tA\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t8\t.\tG\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t9\t.\tC\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t10\t.\tA\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t11\t.\tT\t.\t20\tPASS\t.\tGT\t0/0";

const VCF3: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00003\n\
    20\t1\t.\tGTATGG\tG\t20\tPASS\t.\tGT\t1|1\n\
    20\t7\t.\tA\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t8\t.\tG\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t9\t.\tC\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t10\t.\tA\t.\t20\tPASS\t.\tGT\t0/0\n\
    20\t11\t.\tT\t.\t20\tPASS\t.\tGT\t0/0";

const VCF4: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00004\n\
    1\t4\t.\tGATCGAT\tG\t20\tPASS\t.\tGT\t1|1\n\
    20\t4\t.\tTGGAGCAT\tT\t20\tPASS\t.\tGT\t1|1\n\
    21\t4\t.\tGATCGAT\tG\t20\tPASS\t.\tGT\t1|1";

const VCF5: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00005\n\
    20\t3\t.\tA\tT\t20\tPASS\t.\tGT\t1|1\n\
    20\t8\t.\tG\tA\t20\tPASS\t.\tGT\t1|1\n\
    20\t12\t.\tG\tA\t20\tPASS\t.\tGT\t1|1";

const VCF6: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n\
    20\t1\t.\tG\tA,<NON_REF>\t20\tPASS\t.\tGT\t0|0";

const VCF7: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00002\n\
    20\t1\t.\tG\tA,<NON_REF>\t20\tPASS\t.\tGT\t0|0";

const VCF8: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n\
    20\t1\t.\tT\tA\t20\tPASS\t.\tGT\t0|.\n\
    20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t0|0\n";

const VCF9: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00002\n\
    20\t1\t.\tT\tA\t20\tPASS\t.\tGT\t0|.\n\
    20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t0|0\n";

#[test]
fn test_group_merging_creates_correct_number_of_vars() {
    let iters = vec![make_iter(VCF1), make_iter(VCF2)];
    let grouper = VariantGroupIterator::new(iters, vec!["1".into(), "20".into()]).unwrap();
    let iter_info = grouper.iter_info().to_vec();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();

    let vars = collect_merged_variants(groups.into_iter().map(Ok), &iter_info);

    // Non-variable vars are filtered out
    assert!(vars.len() == 3);
}

#[test]
fn test_non_variant_vars_are_removed() {
    let iters = vec![make_iter(VCF8), make_iter(VCF9)];
    let grouper = VariantGroupIterator::new(iters, vec!["1".into(), "20".into()]).unwrap();
    let iter_info = grouper.iter_info().to_vec();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();
    assert!(groups.len() == 2);

    let vars = collect_merged_variants(groups.into_iter().map(Ok), &iter_info);

    // Non-variable vars are filtered out
    assert!(vars.len() == 0);
}

#[test]
fn test_grouper_preserves_genomic_order() {
    //            chrom 1   chrom 20     chrom 21
    //            4567890   123456789012  4567890
    // REF        GATCGAT   GTATGGAGCATG  GATCGAT
    // VCF1
    // NA00001-1            ATATGGAGCAT
    // NA00001-2            A..TAAAGCAT
    // NA00003-1            G-----AGCAT
    // NA00003-2            G-----AGCAT
    // NA00004-1  G------   GTAT-------   G------
    // NA00004-2  G------   GTAT-------   G------
    // NA00005-1              T     A  A
    // NA00005-1              T     A  A

    let iters = vec![
        make_iter(VCF1),
        make_iter(VCF3),
        make_iter(VCF4),
        make_iter(VCF5),
    ];
    let grouper =
        VariantGroupIterator::new(iters, vec!["1".into(), "20".into(), "21".into()]).unwrap();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();

    let group_positions: Vec<(String, u32)> =
        groups.iter().map(|g| (g.chrom.clone(), g.start)).collect();
    assert_eq!(
        group_positions,
        vec![
            ("1".to_string(), 4),
            ("20".to_string(), 1),
            ("20".to_string(), 12),
            ("21".to_string(), 4),
        ]
    );
}

#[test]
fn test_analyze_single_group_from_merged_bin() {
    // VCF1 + VCF3: all variants merge into one group ("20", 1, 6)
    let iters = vec![make_iter(VCF1), make_iter(VCF3)];
    let grouper = VariantGroupIterator::new(iters, vec!["1".into(), "20".into()]).unwrap();
    let iter_info = grouper.iter_info().to_vec();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();

    assert_eq!(groups.len(), 1);

    let variants = collect_merged_variants(groups.into_iter().map(Ok), &iter_info);

    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].chrom, "20");
    assert_eq!(variants[0].pos, 1);
}

// ============================================================================
// Test helpers for merge algorithm tests
// ============================================================================

/// Creates a Variant from test parameters (mirrors Python's variant_from_dict).
fn make_variant(
    chrom: &str,
    pos: u32,
    alleles: &[&str],
    genotypes_per_sample: &[&[i8]],
    phases: &[bool],
) -> Variant {
    let flat_gts: Vec<i8> = genotypes_per_sample
        .iter()
        .flat_map(|g| g.iter().copied())
        .collect();
    Variant::new(
        chrom.to_string(),
        pos,
        alleles.iter().map(|a| a.to_string()).collect(),
        flat_gts,
        phases.to_vec(),
        genotypes_per_sample.len(),
    )
}

/// Creates an OverlappingVariantGroup from variants and their source iterator indices.
fn make_group(variants: Vec<Variant>, source_idxs: Vec<usize>) -> OverlappingVariantGroup {
    let chrom = variants[0].chrom.clone();
    let start = variants.iter().map(|v| v.pos).min().unwrap();
    let end = variants
        .iter()
        .map(|v| v.pos + v.ref_allele_len as u32 - 1)
        .max()
        .unwrap();
    OverlappingVariantGroup {
        chrom,
        start,
        end,
        variants,
        source_var_iter_idxs: source_idxs,
    }
}

/// Validates a merged variant against expected values (mirrors Python's check_expected_result).
fn check_result(
    variant: &Variant,
    expected_chrom: &str,
    expected_pos: u32,
    expected_ref: &str,
    expected_alts: &[&str],
    expected_gts: &[&[&str]],
    expected_phases: &[bool],
) {
    assert_eq!(variant.chrom, expected_chrom);
    assert_eq!(variant.pos, expected_pos);
    assert_eq!(&variant.alleles[0], expected_ref);

    let actual_alts: HashSet<&str> = variant.alleles[1..].iter().map(|s| s.as_str()).collect();
    let expected_alt_set: HashSet<&str> = expected_alts.iter().copied().collect();
    assert_eq!(actual_alts, expected_alt_set, "Alt alleles mismatch");

    let ploidy = variant.genotypes.len() / variant.n_samples;
    for (si, exp_gt) in expected_gts.iter().enumerate() {
        let start = si * ploidy;
        let gt = &variant.genotypes[start..start + ploidy];
        let actual: Vec<&str> = gt
            .iter()
            .map(|&idx| variant.alleles[idx as usize].as_str())
            .collect();
        assert_eq!(actual, *exp_gt, "Sample {} genotype mismatch", si);
    }

    assert_eq!(variant.phase, expected_phases);
}

// ============================================================================
// Merge algorithm tests (translated from Python test_group_analyzer.py)
// ============================================================================

#[test]
fn test_simple_merge() {
    // ref     AG
    // chrom 1 01
    // s1-1    TC
    // s1-2    AC
    // phase   00
    // s2-1    AT
    // s2-2    TT
    // phase   00
    // merged
    //               s1  s2
    // chrom 1 0 A T T/T T/A
    let var1_s1 = make_variant("1", 10, &["A", "T"], &[&[1, 0]], &[false]);
    let var2_s1 = make_variant("1", 11, &["G", "C"], &[&[1, 1]], &[false]);
    let var1_s2 = make_variant("1", 10, &["A", "T"], &[&[0, 1]], &[false]);
    let var2_s2 = make_variant("1", 11, &["G", "T"], &[&[1, 1]], &[false]);

    let group = make_group(vec![var1_s1, var2_s1, var1_s2, var2_s2], vec![0, 0, 1, 1]);
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        10,
        "AG",
        &["TC", "AC", "AT", "TT"],
        &[&["TC", "AC"], &["AT", "TT"]],
        &[false, false],
    );
}

#[test]
fn test_simple_insertion() {
    // ref    A--G
    // s1-1   A--C
    // s1-2   ATTC
    // phase  0  0
    // s2-1   A--T
    // s2-2   T--T
    // phase  0  0
    let var10_s1 = make_variant("1", 10, &["A", "ATT"], &[&[0, 1]], &[false]);
    let var11_s1 = make_variant("1", 11, &["G", "C"], &[&[1, 1]], &[false]);
    let var10_s2 = make_variant("1", 10, &["A", "T"], &[&[0, 1]], &[false]);
    let var11_s2 = make_variant("1", 11, &["G", "T"], &[&[1, 1]], &[false]);

    let group = make_group(
        vec![var10_s1, var11_s1, var10_s2, var11_s2],
        vec![0, 0, 1, 1],
    );
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        10,
        "AG",
        &["AC", "ATTC", "TT", "AT"],
        &[&["AC", "ATTC"], &["AT", "TT"]],
        &[false, false],
    );
}

#[test]
fn test_simple_deletion() {
    // ref  AT
    // s1-1 A-
    // s1-2 A-
    // s2-1 AT
    // s2-2 TT
    let var1_s1 = make_variant("1", 1, &["AT", "A"], &[&[1, 1]], &[false]);
    let var2_s1 = make_variant("1", 2, &["T"], &[&[0, 0]], &[false]);
    let var1_s2 = make_variant("1", 1, &["A", "T"], &[&[0, 1]], &[false]);
    let var2_s2 = make_variant("1", 2, &["T"], &[&[0, 0]], &[false]);

    let group = make_group(vec![var1_s1, var2_s1, var1_s2, var2_s2], vec![0, 0, 1, 1]);
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "AT",
        &["A", "TT"],
        &[&["A", "A"], &["AT", "TT"]],
        &[false, false],
    );
}

#[test]
fn test_deletion_len_2() {
    // ref  ATT
    // s1-1 A--
    // s1-2 A--
    // s2-1 ATA
    // s2-2 TTA
    let var1_s1 = make_variant("1", 1, &["ATT", "A"], &[&[1, 1]], &[false]);
    let var2_s1 = make_variant("1", 2, &["T"], &[&[0, 0]], &[false]);
    let var3_s1 = make_variant("1", 3, &["T"], &[&[0, 0]], &[false]);
    let var1_s2 = make_variant("1", 1, &["A", "T"], &[&[0, 1]], &[false]);
    let var2_s2 = make_variant("1", 2, &["T"], &[&[0, 0]], &[false]);
    let var3_s2 = make_variant("1", 3, &["T", "A"], &[&[1, 1]], &[false]);

    let group = make_group(
        vec![var1_s1, var2_s1, var3_s1, var1_s2, var2_s2, var3_s2],
        vec![0, 0, 0, 1, 1, 1],
    );
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "ATT",
        &["A", "ATA", "TTA"],
        &[&["A", "A"], &["ATA", "TTA"]],
        &[false, false],
    );
}

#[test]
fn test_overlapping_deletions() {
    // ref  ATT
    // s1-1 A-T
    // s1-2 A-C
    // s2-1 GT-
    // s2-2 GT-
    let var1_s1 = make_variant("1", 1, &["AT", "A"], &[&[1, 1]], &[false]);
    let var2_s1 = make_variant("1", 2, &["T"], &[&[0, 0]], &[false]);
    let var3_s1 = make_variant("1", 3, &["T", "C"], &[&[0, 1]], &[false]);
    let var1_s2 = make_variant("1", 1, &["A", "G"], &[&[1, 1]], &[false]);
    let var2_s2 = make_variant("1", 2, &["TT", "T"], &[&[1, 1]], &[false]);
    let var3_s2 = make_variant("1", 3, &["T"], &[&[0, 0]], &[false]);

    let group = make_group(
        vec![var1_s1, var2_s1, var3_s1, var1_s2, var2_s2, var3_s2],
        vec![0, 0, 0, 1, 1, 1],
    );
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "ATT",
        &["AT", "AC", "GT"],
        &[&["AT", "AC"], &["GT", "GT"]],
        &[false, false],
    );
}

#[test]
fn test_het_deletion() {
    // ref  AT
    // s1-1 AT
    // s1-2 A-
    // s2-1 AT
    // s2-2 TT
    let var1_s1 = make_variant("1", 1, &["AT", "A"], &[&[0, 1]], &[false]);
    let var2_s1 = make_variant("1", 2, &["T"], &[&[0, 0]], &[false]);
    let var1_s2 = make_variant("1", 1, &["A", "T"], &[&[0, 1]], &[false]);
    let var2_s2 = make_variant("1", 2, &["T"], &[&[0, 0]], &[false]);

    let group = make_group(vec![var1_s1, var2_s1, var1_s2, var2_s2], vec![0, 0, 1, 1]);
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "AT",
        &["A", "TT"],
        &[&["AT", "A"], &["AT", "TT"]],
        &[false, false],
    );
}

#[test]
fn test_two_deletions_in_same_sample() {
    // ref  ATG
    // s1-1 A--
    // s1-2 A-G
    // s2-1 GTG
    // s2-2 GT-
    let var1_s1 = make_variant("1", 1, &["ATG", "A", "AG"], &[&[1, 2]], &[false]);
    let var2_s1 = make_variant("1", 2, &["T"], &[&[0, 0]], &[false]);
    let var3_s1 = make_variant("1", 3, &["G"], &[&[0, 0]], &[false]);
    let var1_s2 = make_variant("1", 1, &["A", "G"], &[&[1, 1]], &[false]);
    let var2_s2 = make_variant("1", 2, &["TG", "T"], &[&[0, 1]], &[false]);
    let var3_s2 = make_variant("1", 3, &["G"], &[&[0, 0]], &[false]);

    let group = make_group(
        vec![var1_s1, var2_s1, var3_s1, var1_s2, var2_s2, var3_s2],
        vec![0, 0, 0, 1, 1, 1],
    );
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "ATG",
        &["A", "AG", "GTG", "GT"],
        &[&["A", "AG"], &["GTG", "GT"]],
        &[false, false],
    );
}

#[test]
fn test_two_overlapping_deletions_in_same_sample() {
    // ref  ATG
    // s1-1 AT-
    // s1-2 A-G
    // s2-1 GTG
    // s2-2 GT-
    let var1_s1 = make_variant("1", 1, &["AT", "A"], &[&[0, 1]], &[true]);
    let var2_s1 = make_variant("1", 2, &["TG", "T"], &[&[1, 0]], &[true]);
    let var3_s1 = make_variant("1", 3, &["G"], &[&[0, 0]], &[true]);
    let var1_s2 = make_variant("1", 1, &["A", "G"], &[&[1, 1]], &[true]);
    let var2_s2 = make_variant("1", 2, &["TG", "T"], &[&[0, 1]], &[true]);
    let var3_s2 = make_variant("1", 3, &["G"], &[&[0, 0]], &[true]);

    let group = make_group(
        vec![var1_s1, var2_s1, var3_s1, var1_s2, var2_s2, var3_s2],
        vec![0, 0, 0, 1, 1, 1],
    );
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "ATG",
        &["AT", "AG", "GTG", "GT"],
        &[&["AT", "AG"], &["GTG", "GT"]],
        &[true, true],
    );
}

// ============================================================================
// Phase tracking tests
// ============================================================================

#[test]
fn test_phase_single_het_no_phase_needed() {
    // ref  TA
    // s1-1 TA
    // s1-2 AA
    let var1 = make_variant("1", 1, &["T", "A"], &[&[0, 1]], &[false]);
    let var2 = make_variant("1", 2, &["A"], &[&[0, 0]], &[false]);

    let group = make_group(vec![var1, var2], vec![0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(variant, "1", 1, "TA", &["AA"], &[&["TA", "AA"]], &[false]);
}

#[test]
fn test_phase_kept_between_hets() {
    // ref  TAG
    // s1-1 TAG
    // s1-2 AAC
    let var1 = make_variant("1", 1, &["T", "A"], &[&[0, 1]], &[true]);
    let var2 = make_variant("1", 2, &["A"], &[&[0, 0]], &[true]);
    let var3 = make_variant("1", 3, &["G", "C"], &[&[0, 1]], &[true]);

    let group = make_group(vec![var1, var2, var3], vec![0, 0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "TAG",
        &["AAC"],
        &[&["TAG", "AAC"]],
        &[true],
    );
}

#[test]
fn test_phase_false_on_first_position_still_solvable() {
    // ref  TAG
    // s1-1 TAG
    // s1-2 AAC
    let var1 = make_variant("1", 1, &["T", "A"], &[&[0, 1]], &[false]);
    let var2 = make_variant("1", 2, &["A"], &[&[0, 0]], &[true]);
    let var3 = make_variant("1", 3, &["G", "C"], &[&[0, 1]], &[true]);

    let group = make_group(vec![var1, var2, var3], vec![0, 0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "TAG",
        &["AAC"],
        &[&["TAG", "AAC"]],
        &[false],
    );
}

#[test]
fn test_phase_first_het_not_at_first_position() {
    // ref  ATAG
    // s1-1 ATAG
    // s1-2 AAAC
    let var1 = make_variant("1", 1, &["A"], &[&[0, 0]], &[false]);
    let var2 = make_variant("1", 2, &["T", "A"], &[&[0, 1]], &[false]);
    let var3 = make_variant("1", 3, &["A"], &[&[0, 0]], &[true]);
    let var4 = make_variant("1", 4, &["G", "C"], &[&[0, 1]], &[true]);

    let group = make_group(vec![var1, var2, var3, var4], vec![0, 0, 0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "ATAG",
        &["AAAC"],
        &[&["ATAG", "AAAC"]],
        &[false],
    );
}

#[test]
fn test_phase_broken_between_hets_missing_genotype() {
    // ref  TGAG
    // s1-1 TGAG -> missing because of phasing
    // s1-2 AGAC -> missing because of phasing
    let var1 = make_variant("1", 1, &["T", "A"], &[&[0, 1]], &[true]);
    let var2 = make_variant("1", 2, &["G"], &[&[0, 0]], &[false]); // phase broken!
    let var3 = make_variant("1", 3, &["A"], &[&[0, 0]], &[true]);
    let var4 = make_variant("1", 4, &["G", "C"], &[&[0, 1]], &[true]);

    let group = make_group(vec![var1, var2, var3, var4], vec![0, 0, 0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    // Genotype should be missing (-1) for this sample due to broken phase
    let ploidy = variant.genotypes.len() / variant.n_samples;
    let gt = &variant.genotypes[0..ploidy];
    assert_eq!(gt, &[-1, -1]);
    assert_eq!(variant.phase, &[false]);
}

#[test]
fn test_allele_merge_with_non_ref_filtered() {
    // VCF6 and VCF7 both have ALT=A,<NON_REF> at chr20:1 with GT=0|0.
    // <NON_REF> is stripped during parsing, so alleles should be ["G", "A"].
    let iters = vec![make_iter(VCF6), make_iter(VCF7)];
    let grouper = VariantGroupIterator::new(iters, vec!["20".into()]).unwrap();
    let iter_info = grouper.iter_info().to_vec();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();

    assert_eq!(groups.len(), 1);

    // Both samples are hom-ref, so no alt alleles survive the merge.
    // Non-variable groups are filtered out by analyze_groups.
    let variants = collect_merged_variants(groups.into_iter().map(Ok), &iter_info);

    assert_eq!(variants.len(), 0);
}

#[test]
fn test_phase_broken_in_one_sample() {
    // ref  TG
    // s1-1 TA -> missing
    // s1-2 GC
    // s2-1 TG
    // s2-2 AC
    // Sample 1: phase broken between two hets
    let var1_s1 = make_variant("1", 1, &["T", "A"], &[&[0, 1]], &[true]);
    let var2_s1 = make_variant("1", 2, &["G", "C"], &[&[0, 1]], &[false]);
    // Sample 2: phase maintained
    let var1_s2 = make_variant("1", 1, &["T", "A"], &[&[0, 1]], &[true]);
    let var2_s2 = make_variant("1", 2, &["G", "C"], &[&[0, 1]], &[true]);

    let group = make_group(vec![var1_s1, var2_s1, var1_s2, var2_s2], vec![0, 0, 1, 1]);
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    let ploidy = variant.genotypes.len() / variant.n_samples;

    // Sample 1: missing genotype
    let gt_s1 = &variant.genotypes[0..ploidy];
    assert_eq!(gt_s1, &[-1, -1]);

    // Sample 2: valid genotype (0/1 with ref="TG" and alt containing "AC")
    let gt_s2 = &variant.genotypes[ploidy..2 * ploidy];
    assert_eq!(gt_s2, &[0, 1]);

    let phases = &variant.phase[0..2];
    assert_eq!(phases, &[false, true]);
}

#[test]
fn test_missing_allele_in_merged_haplotype() {
    //           12
    // ref       TC
    // sample1-1  A
    // sample1-2  .
    // sample2-1  -
    // sample2-2  -
    // vcf sample1 chrom1 2 C A 1/.
    // vcf sample2 chrom1 1 TC T 1/1
    // merged vcf chrom1 1 TC TA,T 1/. 2/2

    // sample1: ref T at pos 1, SNP C>A with one missing allele at pos 2
    let var1_s1 = make_variant("chrom1", 1, &["T"], &[&[0, 0]], &[false]);
    let var2_s1 = make_variant("chrom1", 2, &["C", "A"], &[&[1, -1]], &[false]);
    // sample2: deletion TC>T at pos 1, ref C at pos 2
    let var1_s2 = make_variant("chrom1", 1, &["TC", "T"], &[&[1, 1]], &[false]);
    let var2_s2 = make_variant("chrom1", 2, &["C"], &[&[0, 0]], &[false]);

    let group = make_group(vec![var1_s1, var2_s1, var1_s2, var2_s2], vec![0, 0, 1, 1]);
    let iter_infos = vec![
        VariantIteratorInfo {
            samples: vec!["sample1".into()],
        },
        VariantIteratorInfo {
            samples: vec!["sample2".into()],
        },
    ];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    assert_eq!(variant.chrom, "chrom1");
    assert_eq!(variant.pos, 1);
    assert_eq!(&variant.alleles[0], "TC");

    let actual_alts: HashSet<&str> = variant.alleles[1..].iter().map(|s| s.as_str()).collect();
    let expected_alts: HashSet<&str> = ["TA", "T"].iter().copied().collect();
    assert_eq!(actual_alts, expected_alts, "Alt alleles mismatch");

    let ploidy = variant.genotypes.len() / variant.n_samples;

    // Sample 1: one allele should be "TA" and the other missing (-1)
    let gt_s1 = &variant.genotypes[0..ploidy];
    assert!(gt_s1.contains(&-1), "Sample 1 should have a missing allele");
    let non_missing_alleles: Vec<&str> = gt_s1
        .iter()
        .filter(|&&g| g != -1)
        .map(|&g| variant.alleles[g as usize].as_str())
        .collect();
    assert_eq!(
        non_missing_alleles,
        vec!["TA"],
        "Sample 1 non-missing allele should be TA"
    );

    // Sample 2: both alleles should be "T" (the deletion)
    let gt_s2 = &variant.genotypes[ploidy..2 * ploidy];
    let gt_s2_alleles: Vec<&str> = gt_s2
        .iter()
        .map(|&g| variant.alleles[g as usize].as_str())
        .collect();
    assert_eq!(
        gt_s2_alleles,
        vec!["T", "T"],
        "Sample 2 should be homozygous deletion"
    );
}

#[test]
fn test_phase_lost_in_het() {
    // ref   TA
    // s1-1  TA
    // s1-2  TT
    // phase 10 -> phase result false
    let var1_s1 = make_variant("1", 1, &["T"], &[&[0, 0]], &[true]);
    let var2_s1 = make_variant("1", 2, &["A", "T"], &[&[0, 1]], &[false]);

    let group = make_group(vec![var1_s1, var2_s1], vec![0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(variant, "1", 1, "TA", &["TT"], &[&["TA", "TT"]], &[false]);
}

#[test]
fn test_phase_lost_in_het_and_not_recovered() {
    // ref   TAA
    // s1-1  TAA
    // s1-2  TTT
    // phase 101 -> phase result false
    let var1_s1 = make_variant("1", 1, &["T"], &[&[0, 0]], &[true]);
    let var2_s1 = make_variant("1", 2, &["A", "T"], &[&[0, 1]], &[false]);
    let var3_s1 = make_variant("1", 3, &["A", "T"], &[&[0, 1]], &[true]);

    let group = make_group(vec![var1_s1, var2_s1, var3_s1], vec![0, 0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "TAA",
        &["TTT"],
        &[&["TAA", "TTT"]],
        &[false],
    );
}

#[test]
fn test_phase_lost_in_het2() {
    // ref   TAA
    // s1-1  TAA -> missing because of phasing
    // s1-2  TTT -> missing because of phasing
    // phase 110 -> two hets without phase between them -> missing genotype
    let var1_s1 = make_variant("1", 1, &["T"], &[&[0, 0]], &[true]);
    let var2_s1 = make_variant("1", 2, &["A", "T"], &[&[0, 1]], &[true]);
    let var3_s1 = make_variant("1", 3, &["A", "T"], &[&[0, 1]], &[false]);

    let group = make_group(vec![var1_s1, var2_s1, var3_s1], vec![0, 0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    // Genotype should be missing (-1) for this sample due to broken phase
    let ploidy = variant.genotypes.len() / variant.n_samples;
    let gt = &variant.genotypes[0..ploidy];
    assert_eq!(gt, &[-1, -1]);
    assert_eq!(variant.phase, &[false]);
}

#[test]
fn test_phase_maintained_despited_not_phased_in_hom_position() {
    // ref   TAA
    // s1-1  TAA
    // s1-2  TAT
    // phase 101 -> phase result true
    let var1_s1 = make_variant("1", 1, &["T"], &[&[0, 0]], &[true]);
    let var2_s1 = make_variant("1", 2, &["A"], &[&[0, 0]], &[false]);
    let var3_s1 = make_variant("1", 3, &["A", "T"], &[&[0, 1]], &[true]);

    let group = make_group(vec![var1_s1, var2_s1, var3_s1], vec![0, 0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(
        variant,
        "1",
        1,
        "TAA",
        &["TAT"],
        &[&["TAA", "TAT"]],
        &[true],
    );
}

#[test]
fn test_phase_conserved_in_het() {
    // ref   TA
    // s1-1  TA
    // s1-2  TT
    // phase 11 -> phase result true
    let var1_s1 = make_variant("1", 1, &["T"], &[&[0, 0]], &[true]);
    let var2_s1 = make_variant("1", 2, &["A", "T"], &[&[0, 1]], &[true]);

    let group = make_group(vec![var1_s1, var2_s1], vec![0, 0]);
    let iter_infos = vec![VariantIteratorInfo {
        samples: vec!["sample1".into()],
    }];

    let result = merge_variant_group(&group, &iter_infos, &default_prior()).unwrap();
    let variant = &result[0];

    check_result(variant, "1", 1, "TA", &["TT"], &[&["TA", "TT"]], &[true]);
}

// TODO If you don't get data for one position in the merger it should return an error
// TODO if the reference allele is different in one sample it should fail
