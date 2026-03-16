from variant_group_analyzer import (
    VariantGroupIterable,
    VariantIteratorInfo,
    OverlappingVariantGroup,
    variant_from_dict,
    merge_variant_group,
    Variant,
)
from typing import Iterator


def variant_generator(variants: list[dict]) -> Iterator[Variant]:
    for variant in variants:
        yield variant_from_dict(variant)


def create_variant_group_iterator(
    vars_in_vcfs: list[list[dict]], samples_in_vcfs: list[list[str]]
):
    # for every VCF there should be a list of samples
    assert len(vars_in_vcfs) == len(samples_in_vcfs)
    iter_infos = []
    vars_iters = []
    for vars_in_vcf, samples in zip(vars_in_vcfs, samples_in_vcfs):
        iter_infos.append(VariantIteratorInfo(samples))
        vars_iters.append(variant_generator(vars_in_vcf))

    return VariantGroupIterable(variant_iterators=vars_iters, var_iter_infos=iter_infos)


def test_simple_merge():
    var1_sample1 = {
        "chrom": 1,
        "pos": 10,
        "alleles": ["A", "T"],
        "genotypes": [[1, 1]],
        "phases": [False],
    }
    var2_sample1 = {
        "chrom": 1,
        "pos": 11,
        "alleles": ["G", "C"],
        "genotypes": [[0, 0]],
        "phases": [False],
    }
    vars_in_vfc1 = [var1_sample1, var2_sample1]
    samples_in_vcf1 = ["sample1"]
    var1_sample2 = {
        "chrom": 1,
        "pos": 10,
        "alleles": ["A", "T"],
        "genotypes": [[0, 1]],
        "phases": [False],
    }
    var2_sample2 = {
        "chrom": 1,
        "pos": 11,
        "alleles": ["A", "C"],
        "genotypes": [[1, 1]],
        "phases": [False],
    }
    vars_in_vfc2 = [var1_sample2, var2_sample2]
    samples_in_vcf2 = ["sample2"]

    var_group_iter = create_variant_group_iterator(
        vars_in_vcfs=[vars_in_vfc1, vars_in_vfc2],
        samples_in_vcfs=[samples_in_vcf1, samples_in_vcf2],
    )
    var_group: OverlappingVariantGroup = list(var_group_iter)[0]
    merge_variant_group(
        var_group,
        samples_have_one_var_per_position=True,
        var_iter_infos=var_group_iter.var_iter_infos,
    )
