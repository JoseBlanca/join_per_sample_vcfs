from dataclasses import dataclass
from typing import Iterator


@dataclass
class VariantIteratorInfo:
    samples: list[str]


@dataclass
class Variant:
    chrom: str
    pos: int
    alleles: list[str]
    genotypes: list[int]
    phase: list[bool]


def variant_new(
    chrom: str, pos: int, alleles: list[str], genotypes: list[int], phase: list[bool]
) -> Variant:
    return Variant(
        chrom=chrom, pos=pos, alleles=alleles, genotypes=genotypes, phase=phase
    )


def variant_from_dict(variant: dict) -> Variant:
    return variant_new(**variant)


def variant_generator(variants: list[Variant]) -> Iterator[Variant]:
    for variant in variants:
        yield variant


class VariantGroupIterator:
    def __init__(
        self,
        variant_iterators: list[Iterator[Variant]],
        iter_infos: list[VariantIteratorInfo],
    ):
        self.variant_iterators = variant_iterators
        self.iter_infos = iter_infos

    def __next__(self):
        var_group = [next(var_iter) for var_iter in self.variant_iterators]
        return var_group

        return next(self.variant_iterators)

    def __iter__(self):
        return self


class VariantGroupIterable:
    def __init__(
        self,
        variant_iterators: list[Iterator[Variant]],
        iter_infos: list[VariantIteratorInfo],
    ):
        assert len(variant_iterators) == len(iter_infos)
        self.variant_iterators = variant_iterators
        self.iter_infos = iter_infos

    def __iter__(self):
        return VariantGroupIterator(self.variant_iterators, self.iter_infos)


def create_variant_group_iterator(variant_dicts: list[dict], samples: list[str]):
    assert len(variant_dicts) == len(samples)
    variants = [variant_from_dict(variant) for variant in variant_dicts]
    variant_generators = [variant_generator([variant]) for variant in variants]
    iter_infos = [VariantIteratorInfo(samples=[sample]) for sample in samples]

    variant_group_iterable = VariantGroupIterable(variant_generators, iter_infos)
    return iter(variant_group_iterable)


def test():
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
    print(list(var_group_iter))


def main():
    test()


if __name__ == "__main__":
    main()
