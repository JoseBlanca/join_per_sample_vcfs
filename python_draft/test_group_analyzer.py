from variant_group_analyzer import (
    VariantGroupIterable,
    VariantIteratorInfo,
    variant_from_dict,
    variant_generator,
    merge_variant_group,
)


def create_variant_group_iterator(variant_dicts: list[dict], samples: list[str]):
    assert len(variant_dicts) == len(samples)
    variants = [variant_from_dict(variant) for variant in variant_dicts]
    variant_generators = [variant_generator([variant]) for variant in variants]
    iter_infos = [VariantIteratorInfo(samples=[sample]) for sample in samples]

    variant_group_iterable = VariantGroupIterable(variant_generators, iter_infos)
    return iter(variant_group_iterable)


def test_boiler_plate_for_var_group_creation():
    variant_sample1 = {
        "chrom": 1,
        "pos": 10,
        "alleles": ["A", "T"],
        "genotypes": [[1, 1]],
        "phase": [False],
    }
    variant_sample2 = {
        "chrom": 1,
        "pos": 10,
        "alleles": ["A", "T"],
        "genotypes": [[0, 1]],
        "phase": [False],
    }
    var_group_iter = create_variant_group_iterator(
        [variant_sample1, variant_sample2], ["sample1", "sample2"]
    )
    var_group, samples = list(var_group_iter)[0]
    assert var_group.start == 10
    assert var_group.end == 10
    assert samples == [["sample1"], ["sample2"]]


def test_simple_merge():
    variant_sample1 = {
        "chrom": 1,
        "pos": 10,
        "alleles": ["A", "T"],
        "genotypes": [[1, 1]],
        "phase": [False],
    }
    variant_sample2 = {
        "chrom": 1,
        "pos": 10,
        "alleles": ["A", "T"],
        "genotypes": [[0, 1]],
        "phase": [False],
    }
    var_group_iter = create_variant_group_iterator(
        [variant_sample1, variant_sample2], ["sample1", "sample2"]
    )
    var_group, samples = list(var_group_iter)[0]
    merge_variant_group(
        var_group,
        samples_have_one_var_per_position=True,
        samples_for_variants=samples,
    )
