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
    phases: list[bool]


def variant_new(
    chrom: str, pos: int, alleles: list[str], genotypes: list[int], phases: list[bool]
) -> Variant:
    return Variant(
        chrom=chrom, pos=pos, alleles=alleles, genotypes=genotypes, phases=phases
    )


def var_len(variant: Variant):
    return len(variant.alleles[0])


def var_end(variant: Variant):
    return var_len(variant) + variant.pos - 1


def variant_from_dict(variant: dict) -> Variant:
    return variant_new(**variant)


@dataclass
class OverlappingVariantGroup:
    chrom: str
    start: int
    end: int
    variants: list[Variant]
    var_iter_of_origin: list[
        int
    ]  # from which VCF variant iterator did each of the variants originated

    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end} ({len(self.variants)})"


def create_variant_for_region(
    var_group: OverlappingVariantGroup,
    region: tuple[int, int],
    var_iter_infos: list[VariantIteratorInfo],
) -> Variant:

    alleles_for_samples = {}
    for variant, var_iter_of_origin in zip(
        var_group.variants, var_group.var_iter_of_origin
    ):
        alleles = variant.alleles
        for sample_idx_in_var_iter, (sample_gt, sample_phase) in enumerate(
            zip(variant.genotypes, variant.phases)
        ):
            sample_id = (var_iter_of_origin, sample_idx_in_var_iter)
            if sample_id not in alleles_for_samples:
                alleles_for_samples[sample_id] = [[] for _ in range(len(sample_gt))]
            print(sample_id, sample_gt, sample_phase)
            for haplo_chrom_idx, sample_allele_haplo_int in enumerate(sample_gt):
                sample_allele_haplo = alleles[sample_allele_haplo_int]
                alleles_for_samples[sample_id][haplo_chrom_idx].append(
                    sample_allele_haplo
                )
    print(alleles_for_samples)


def merge_variant_group(
    variant_group: OverlappingVariantGroup,
    var_iter_infos: list[VariantIteratorInfo],
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
        variant = create_variant_for_region(
            variant_group, region, var_iter_infos=var_iter_infos
        )
        variants.append(variant)

    return Ok(variants)
