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
