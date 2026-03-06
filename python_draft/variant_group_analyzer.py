from dataclasses import dataclass
from typing import Iterator, Iterable


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


@dataclass
class OverlappingVariantGroup:
    chrom: str
    start: int
    end: int
    variants: list[Variant]
    source_var_iter_idxs: list[int]

    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end} ({len(self.variants)})"


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

    def __next__(self) -> OverlappingVariantGroup:
        variants = [next(var_iter) for var_iter in self.variant_iterators]
        start = min([variant.pos for variant in variants])
        end = max([len(variant.alleles[0]) + variant.pos - 1 for variant in variants])
        var_group = OverlappingVariantGroup(
            variants=variants,
            chrom=variants[0].chrom,
            start=start,
            end=end,
            source_var_iter_idxs=list(range(len(self.variant_iterators))),
        )
        return var_group

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

    def __iter__(self) -> Iterable[OverlappingVariantGroup]:
        return VariantGroupIterator(self.variant_iterators, self.iter_infos)
