from dataclasses import dataclass
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
) -> tuple[Variant, list[tuple[int, int]]]:

    # create the alleles for each sample by adding the alleles for each var
    # alleles_for_samples dict
    # key: sample_id
    #   sample_id, tuple with two items:
    #       - VariantIterator in which the variant is located (int index)
    #       - int index of that sample in that VariantIterator
    # value is a list of alleles for the corresponding sample, one per haploid chromosome.
    # each item in the list is a list of str alleles, one per variant
    alleles_for_samples = {}
    positions_left_in_del = {}
    for variant, var_iter_of_origin in zip(
        var_group.variants, var_group.var_iter_of_origin
    ):
        alleles = variant.alleles
        for sample_idx_in_var_iter, (sample_gt, sample_phase) in enumerate(
            zip(variant.genotypes, variant.phases)
        ):
            sample_id = (var_iter_of_origin, sample_idx_in_var_iter)

            ploidy = len(sample_gt)

            if sample_id not in alleles_for_samples:
                alleles_for_samples[sample_id] = [[] for _ in range(len(sample_gt))]
                positions_left_in_del[sample_id] = [0] * ploidy

            ref_allele = variant.alleles[0]
            if len(ref_allele) > 1:
                deletion_created = True
                positions_left_in_del_for_sample = []
                for sample_allele_haplo_int in sample_gt:
                    sample_allele = variant.alleles[sample_allele_haplo_int]
                    print(f"{sample_allele=}")
                    del_len = len(ref_allele) - len(sample_allele)
                    positions_left_in_del_for_sample.append(del_len)
                positions_left_in_del[sample_id] = positions_left_in_del_for_sample
                print(f"{positions_left_in_del=}")
            else:
                deletion_created = False

            for haplo_chrom_idx, sample_allele_haplo_int in enumerate(sample_gt):
                sample_allele_haplo = alleles[sample_allele_haplo_int]

                # the allele could have more than 1 nucleotides when is the first in a deletion
                # but we should only add the first nucleotide, the one corresponding to this position
                sample_allele_haplo = sample_allele_haplo[0]

                if (
                    positions_left_in_del[sample_id][haplo_chrom_idx] > 0
                    and not deletion_created
                ):
                    sample_allele_haplo = ""
                    positions_left_in_del[sample_id][haplo_chrom_idx] -= 1

                alleles_for_samples[sample_id][haplo_chrom_idx].append(
                    sample_allele_haplo
                )

    print(f"{alleles_for_samples=}")

    # create reference alleles
    # allele_ids dict:
    #   - value
    ref_allele_per_var_iter = []
    for _ in var_iter_infos:
        ref_allele_per_var_iter.append([])

    for variant, var_iter_of_origin in zip(
        var_group.variants, var_group.var_iter_of_origin
    ):
        var_ref_allele = variant.alleles[0]
        ref_allele_per_var_iter[var_iter_of_origin].append(var_ref_allele[0])
    ref_allele = "".join(ref_allele_per_var_iter[0])
    print(f"{ref_allele_per_var_iter=}")
    for other_ref_per_var_iter in ref_allele_per_var_iter[1:]:
        assert ref_allele == "".join(other_ref_per_var_iter)

    allele_ids = {ref_allele: 0}
    # create str alleles and modify alleles_for_samples to hold strs instead of lists
    for sample_id, alleles_for_sample in alleles_for_samples.items():
        alleles_for_sample_str = []
        for allele_list in alleles_for_sample:
            allele = "".join(allele_list)

            if allele in allele_ids:
                allele_id = allele_ids[allele]
            else:
                allele_id = len(allele_ids)
                allele_ids[allele] = allele_id

            alleles_for_sample_str.append(allele)
        alleles_for_samples[sample_id] = alleles_for_sample_str

    # create sample and genotypes for each sample
    sample_ids = []
    genotypes = []
    for var_iter_idx, var_iter_info in enumerate(var_iter_infos):
        for sample_idx, sample in enumerate(var_iter_info.samples):
            sample_id = (var_iter_idx, sample_idx)
            sample_ids.append(sample_id)
            alleles_for_sample = alleles_for_samples[sample_id]
            sample_gts = tuple(allele_ids[allele] for allele in alleles_for_sample)
            genotypes.append(sample_gts)

    for sample_id, alleles_for_sample in alleles_for_samples.items():
        print(f"{sample_id=}")
        var_iter_of_origin, sample_idx_in_var_iter = sample_id
        sample = var_iter_infos[var_iter_of_origin].samples[sample_idx_in_var_iter]
        print(f"{sample=}")
        print(f"{alleles_for_sample=}")

    alleles = list(allele_ids.keys())
    alleles = sorted(alleles, key=lambda x: allele_ids[x])

    print(f"{genotypes=}")
    print(f"{alleles=}")
    phases = [False] * len(sample_ids)
    variant = Variant(
        chrom=var_group.chrom,
        pos=var_group.start,
        alleles=alleles,
        genotypes=genotypes,
        phases=phases,
    )
    return variant, sample_ids


def merge_variant_group(
    variant_group: OverlappingVariantGroup,
    var_iter_infos: list[VariantIteratorInfo],
    samples_have_one_var_per_position: bool,
) -> Result[tuple[list[Variant], list[tuple[int, int]]], str]:

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
        variant, sample_ids = create_variant_for_region(
            variant_group, region, var_iter_infos=var_iter_infos
        )
        variants.append(variant)

    return Ok((variants, sample_ids))
