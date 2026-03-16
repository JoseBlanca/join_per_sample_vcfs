from variant_group_analyzer import (
    VariantIteratorInfo,
    OverlappingVariantGroup,
    variant_from_dict,
    merge_variant_group,
    Variant,
    var_end,
    Ok,
    Err,
)
from typing import Iterator
import pytest


class VariantGroupIterator:
    def __init__(
        self,
        variant_iterators: list[Iterator[Variant]],
        iter_infos: list[VariantIteratorInfo],
    ):
        self.variant_iterators = variant_iterators
        self.var_iter_infos = iter_infos

    def __next__(self) -> OverlappingVariantGroup:
        vars_in_group: list[Variant] = []
        var_iter_of_origin: list[int] = []
        chrom = None
        var_starts = []
        var_ends = []
        for var_iter_idx, var_iter in enumerate(self.variant_iterators):
            for var in var_iter:
                if chrom is None:
                    chrom = var.chrom
                else:
                    assert chrom == var.chrom
                var_starts.append(var.pos)
                var_ends.append(var_end(var))

                vars_in_group.append(var)
                var_iter_of_origin.append(var_iter_idx)

        if not vars_in_group:
            raise StopIteration

        start = min(var_starts)
        end = max(var_ends)
        var_group = OverlappingVariantGroup(
            variants=vars_in_group,
            chrom=chrom,
            start=start,
            end=end,
            var_iter_of_origin=var_iter_of_origin,
        )
        return var_group

    def __iter__(self):
        return self


class VariantGroupIterable:
    def __init__(
        self,
        variant_iterators: list[Iterator[Variant]],
        var_iter_infos: list[VariantIteratorInfo],
    ):
        assert len(variant_iterators) == len(var_iter_infos)

        samples_seen = set()
        for iter_info in var_iter_infos:
            for sample in iter_info.samples:
                if sample in samples_seen:
                    raise RuntimeError(
                        f"One sample should not be found in more than one VariantIterator: {sample}"
                    )
                samples_seen.add(sample)

        self.variant_iterators = variant_iterators
        self.var_iter_infos = var_iter_infos

    def __iter__(self) -> VariantGroupIterator:
        return VariantGroupIterator(self.variant_iterators, self.var_iter_infos)


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


def check_expected_result(result, expected):
    if isinstance(result, Ok):
        variants, sample_ids = result.value
        variant = variants[0]
    else:
        pytest.fail("Got error, but Ok expected")

    assert expected["chrom"] == variant.chrom
    assert expected["pos"] == variant.pos
    assert expected["ref_allele"] == variant.alleles[0]
    assert expected["alt_alleles"] == set(variant.alleles[1:])

    alleles = dict(enumerate(variant.alleles))
    for var_sample_gts, exp_sample_gts in zip(variant.genotypes, expected["genotypes"]):
        var_sample_gts = [alleles[allele] for allele in var_sample_gts]
        assert var_sample_gts == exp_sample_gts

    assert variant.phases == expected["phases"]


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
        "genotypes": [[1, 1]],
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
        "alleles": ["G", "T"],
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
    result = merge_variant_group(
        var_group,
        samples_have_one_var_per_position=True,
        var_iter_infos=var_group_iter.var_iter_infos,
    )

    expected = {
        "chrom": 1,
        "pos": 10,
        "ref_allele": "AG",
        "alt_alleles": {"TC", "AT", "TT"},
        "genotypes": [["TC", "TC"], ["AT", "TT"]],
        "samples": ["sample1", "sample2"],
        "phases": [False, False],
    }
    check_expected_result(result, expected)


def test_simple_deletion():
    var1_sample1 = {
        "chrom": 1,
        "pos": 10,
        "alleles": ["AT", "A"],
        "genotypes": [[1, 1]],
        "phases": [False],
    }
    var2_sample1 = {
        "chrom": 1,
        "pos": 11,
        "alleles": ["T"],
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
        "alleles": ["T"],
        "genotypes": [[0, 0]],
        "phases": [False],
    }
    vars_in_vfc2 = [var1_sample2, var2_sample2]
    samples_in_vcf2 = ["sample2"]

    var_group_iter = create_variant_group_iterator(
        vars_in_vcfs=[vars_in_vfc1, vars_in_vfc2],
        samples_in_vcfs=[samples_in_vcf1, samples_in_vcf2],
    )
    var_group: OverlappingVariantGroup = list(var_group_iter)[0]
    result = merge_variant_group(
        var_group,
        samples_have_one_var_per_position=True,
        var_iter_infos=var_group_iter.var_iter_infos,
    )

    expected = {
        "chrom": 1,
        "pos": 10,
        "ref_allele": "AT",
        "alt_alleles": {"A", "TT"},
        "genotypes": [["A", "A"], ["AT", "TT"]],
        "samples": ["sample1", "sample2"],
        "phases": [False, False],
    }
    check_expected_result(result, expected)
