use std::io::BufReader;

use join_per_sample_vcfs::gvcf_parser::VariantIterator;
use join_per_sample_vcfs::variant_group::VariantGroupIterator;
use join_per_sample_vcfs::variant_group_analyzer::analyze_groups;

fn make_iter(vcf_data: &str) -> VariantIterator<BufReader<BufReader<&[u8]>>> {
    let reader = BufReader::new(vcf_data.as_bytes());
    VariantIterator::from_reader(reader).expect("Failed to create parser")
}

const VCF1: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n\
    20\t1\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t.|0\n\
    20\t3\t.\tA\tG\t20\tPASS\t.\tGT\t.|0\n\
    20\t4\t.\tT\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    20\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0/1\n\
    20\t6\t.\tG\tA\t20\tPASS\t.\tGT\t0/1";

const VCF2: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00002\n\
    20\t1\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t.|0\n\
    20\t3\t.\tA\tG\t20\tPASS\t.\tGT\t.|0\n\
    20\t4\t.\tT\t<NON_REF>\t20\tPASS\t.\tGT\t0|0\n\
    20\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0/1\n\
    20\t6\t.\tG\tA\t20\tPASS\t.\tGT\t0/1";

const VCF3: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00003\n\
    20\t1\t.\tGTATGG\tG\t20\tPASS\t.\tGT\t0|0";

const VCF4: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00004\n\
    1\t4\t.\tGATCGAT\tG\t20\tPASS\t.\tGT\t0|0\n\
    20\t4\t.\tTGGAGCAT\tT\t20\tPASS\t.\tGT\t0|0\n\
    21\t4\t.\tGATCGAT\tG\t20\tPASS\t.\tGT\t0|0";

const VCF5: &str = "##\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00005\n\
    20\t3\t.\tA\tT\t20\tPASS\t.\tGT\t0|0\n\
    20\t8\t.\tG\tA\t20\tPASS\t.\tGT\t0|0\n\
    20\t12\t.\tG\tA\t20\tPASS\t.\tGT\t0|0";

#[test]
fn test_analyze_produces_one_variant_per_group() {
    let iters = vec![make_iter(VCF1), make_iter(VCF2)];
    let grouper = VariantGroupIterator::new(iters, vec!["1".into(), "20".into()]).unwrap();
    let iter_info = grouper.iter_info().to_vec();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();
    let n_groups = groups.len();

    let variants = analyze_groups(groups, &iter_info).unwrap();

    assert_eq!(variants.len(), n_groups);
}

#[test]
fn test_analyze_preserves_genomic_order() {
    let iters = vec![make_iter(VCF1), make_iter(VCF3), make_iter(VCF4), make_iter(VCF5)];
    let grouper = VariantGroupIterator::new(
        iters,
        vec!["1".into(), "20".into(), "21".into()],
    )
    .unwrap();
    let iter_info = grouper.iter_info().to_vec();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();

    let group_positions: Vec<(String, u32)> = groups
        .iter()
        .map(|g| (g.chrom.clone(), g.start))
        .collect();

    let variants = analyze_groups(groups, &iter_info).unwrap();

    let variant_positions: Vec<(String, u32)> = variants
        .iter()
        .map(|v| (v.chrom.clone(), v.pos))
        .collect();

    // Output variants must match the order of the input groups
    assert_eq!(variant_positions, group_positions);
}

#[test]
fn test_analyze_single_group_from_merged_bin() {
    // VCF1 + VCF3: all variants merge into one group ("20", 1, 6)
    let iters = vec![make_iter(VCF1), make_iter(VCF3)];
    let grouper = VariantGroupIterator::new(iters, vec!["1".into(), "20".into()]).unwrap();
    let iter_info = grouper.iter_info().to_vec();
    let groups: Vec<_> = grouper.map(|r| r.unwrap()).collect();

    assert_eq!(groups.len(), 1);

    let variants = analyze_groups(groups, &iter_info).unwrap();

    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].chrom, "20");
    assert_eq!(variants[0].pos, 1);
}
