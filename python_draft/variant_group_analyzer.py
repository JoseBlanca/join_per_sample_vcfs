from dataclasses import dataclass
from typing import Iterator, Iterable
from typing import Generic, TypeVar, Union

T = TypeVar("T")
E = TypeVar("E")


@dataclass(frozen=True)
class Ok(Generic[T]):
    value: T


@dataclass(frozen=True)
class Err(Generic[E]):
    error: E


Result = Union[Ok[T], Err[E]]


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
    sample_origins: list[tuple[int, int]]  # list[(iter_idx, sample_idx_in_vcf)]

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

    def __next__(self) -> tuple[OverlappingVariantGroup, list[list[str]]]:
        variants: list[Variant] = []
        sample_origins: list[tuple[int, int]] = []
        for iter_idx, (var_iter, iter_info) in enumerate(
            zip(self.variant_iterators, self.iter_infos)
        ):
            for sample_idx_in_iter, sample in enumerate(iter_info.samples):
                sample_origins.append((iter_idx, sample_idx_in_iter))
                variants.append(next(var_iter))

        variant_samples: list[list[str]] = [
            iter_info.samples for iter_info in self.iter_infos
        ]

        start = min([variant.pos for variant in variants])
        end = max([len(variant.alleles[0]) + variant.pos - 1 for variant in variants])
        var_group = OverlappingVariantGroup(
            variants=variants,
            chrom=variants[0].chrom,
            start=start,
            end=end,
            sample_origins=sample_origins,
        )
        return var_group, variant_samples

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

    def __iter__(self) -> Iterator[tuple[OverlappingVariantGroup, list[list[str]]]]:
        return VariantGroupIterator(self.variant_iterators, self.iter_infos)


def create_variant_for_region(
    variant_group: OverlappingVariantGroup,
    region: tuple[int, int],
    samples_for_variants: list[list[str]],
):

    samples = []
    sample_origins = []
    for iter_idx, samples_in_var in enumerate(samples_for_variants):
        for sample_idx_in_iter, sample in enumerate(samples_in_var):
            samples.append(sample)
            sample_origins.append((iter_idx, sample_idx_in_iter))

    samples = [
        sample for samples_in_var in samples_for_variants for sample in samples_in_var
    ]
    print(samples)

    for samples_in_var, variant in zip(samples_for_variants, variant_group.variants):
        print(variant)
        print(samples_in_var)


def merge_variant_group(
    variant_group: OverlappingVariantGroup,
    samples_for_variants: list[list[str]],
    samples_have_one_var_per_position: bool,
) -> Result[list[Variant], str]:

    if not samples_have_one_var_per_position:
        return Err(
            "NotImplemented, only implemented for VCFs with one line per position"
        )

    # in some cases all the region covered by the group might be solved in a single variant
    # the problem happens when we have to heterozygotes and the phase between them is unknown
    # right now we a going to set the genotype for these heterozyous samples as missing and, thus
    # we are going to workaround the problem
    # so the region covered by the resulting variant will the same as the region covered by the
    # group and only one variant will be constructed for a given group.
    variant_regions: list[tuple[int, int]] = [(variant_group.start, variant_group.end)]

    variants: list[Variant] = []
    for region in variant_regions:
        variant = create_variant_for_region(variant_group, region, samples_for_variants)
        variants.append(variant)

    return Ok(variants)
