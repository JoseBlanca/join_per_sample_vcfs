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
) -> tuple[Variant, list[int]]:

    # how many samples are per vcf
    n_samples_per_var_iter = []
    var_iter_seen = set()
    for variant, var_iter_of_origin_idx in zip(
        var_group.variants, var_group.var_iter_of_origin
    ):
        if var_iter_of_origin_idx not in var_iter_seen:
            samples = var_iter_infos[var_iter_of_origin_idx].samples
            n_samples_per_var_iter.append(len(samples))
            var_iter_seen.add(var_iter_of_origin_idx)

    total_samples = sum(n_samples_per_var_iter)

    # create the alleles for each sample by adding the alleles for each var
    # alleles_for_samples list indexed by sample_idx
    # sample_idx = sum(n_samples_per_var_iter[:var_iter_of_origin_idx]) + sample_idx_in_var_iter
    # value is a list of alleles for the corresponding sample, one per haploid chromosome.
    # each item in the list is a list of str alleles, one per variant
    alleles_for_samples = [None] * total_samples
    positions_left_in_del = [None] * total_samples
    # Phase tracking: detect when haplotype can't be built due to broken phase
    first_het_seen = [False] * total_samples
    phase_broken_since_het = [False] * total_samples
    missing_samples = set()
    # Let's go through every variant of every variant iterator (VCF)
    # There should be one variant per genomic position and all genomic
    # positions should be covered by one variant
    for variant, var_iter_of_origin_idx in zip(
        var_group.variants, var_group.var_iter_of_origin
    ):
        print(f"{var_iter_of_origin_idx=}")
        print(f"{variant.pos=}")
        alleles = variant.alleles
        # Now we go through every genotype of every variant
        # One genotype corresponds to one sample, so this for can be interpreted as going
        # through each sample
        for sample_idx_in_var_iter, (sample_gt, sample_phase) in enumerate(
            zip(variant.genotypes, variant.phases)
        ):
            sample_idx: int = (
                sum(n_samples_per_var_iter[:var_iter_of_origin_idx])
                + sample_idx_in_var_iter
            )
            print(f"{var_iter_of_origin_idx=} {sample_idx_in_var_iter=} {sample_idx=}")

            ploidy = len(sample_gt)
            is_first_variant_for_sample = alleles_for_samples[sample_idx] is None

            if alleles_for_samples[sample_idx] is None:
                # we create the data structure to store the alleles that we are going to be building
                alleles_for_samples[sample_idx] = [[] for _ in range(ploidy)]
                # We need to keep track of the genomic positions to ignore after a deletion is found
                positions_left_in_del[sample_idx] = [0] * ploidy

            ref_allele = variant.alleles[0]
            prev_positions_left_in_del = positions_left_in_del[sample_idx]
            if (
                len(ref_allele) > 1
            ):  # This will be a deletion in the following positions
                deletion_created = True
                positions_left_in_del_for_sample = []
                # we update the deletion length
                for sample_allele_haplo_int, left_in_del in zip(
                    sample_gt, positions_left_in_del[sample_idx]
                ):
                    sample_allele = variant.alleles[sample_allele_haplo_int]
                    print(f"{sample_allele=}")
                    del_len = len(ref_allele) - len(sample_allele)
                    positions_left_in_del_for_sample.append(del_len + left_in_del)
                positions_left_in_del[sample_idx] = positions_left_in_del_for_sample
                print(f"{positions_left_in_del=}")
            else:
                deletion_created = False

            # Phase tracking: check if haplotype can be built for this sample
            is_het = len(set(sample_gt)) > 1
            if not is_first_variant_for_sample and first_het_seen[sample_idx]:
                if not sample_phase:
                    phase_broken_since_het[sample_idx] = True
            if is_het:
                if first_het_seen[sample_idx] and phase_broken_since_het[sample_idx]:
                    missing_samples.add(sample_idx)
                first_het_seen[sample_idx] = True
                phase_broken_since_het[sample_idx] = False

            # Now we add the corresponding nucleotide to each the alelle of each haplotype
            for haplo_chrom_idx, sample_allele_haplo_int in enumerate(sample_gt):
                sample_allele_haplo = alleles[sample_allele_haplo_int]

                # When the ref allele spans multiple positions (deletion), the allele may have
                # extra characters for subsequent positions. We only keep the first nucleotide.
                # But for insertions (allele longer than ref with single-base ref), keep the full allele.
                if len(ref_allele) > 1:
                    sample_allele_haplo = sample_allele_haplo[0]

                if positions_left_in_del[sample_idx][haplo_chrom_idx] > 0 and (
                    not deletion_created
                    or prev_positions_left_in_del[haplo_chrom_idx] > 0
                ):
                    sample_allele_haplo = ""
                    positions_left_in_del[sample_idx][haplo_chrom_idx] -= 1

                alleles_for_samples[sample_idx][haplo_chrom_idx].append(
                    sample_allele_haplo
                )
    print(f"{alleles_for_samples=}")

    # create reference alleles
    # allele_ids dict:
    #   - value
    ref_allele_per_var_iter = []
    for _ in var_iter_infos:
        ref_allele_per_var_iter.append([])

    for variant, var_iter_of_origin_idx in zip(
        var_group.variants, var_group.var_iter_of_origin
    ):
        var_ref_allele = variant.alleles[0]
        ref_allele_per_var_iter[var_iter_of_origin_idx].append(var_ref_allele[0])
    ref_allele = "".join(ref_allele_per_var_iter[0])
    print(f"{ref_allele_per_var_iter=}")
    for other_ref_per_var_iter in ref_allele_per_var_iter[1:]:
        assert ref_allele == "".join(other_ref_per_var_iter)

    allele_ids = {ref_allele: 0}
    # create str alleles and modify alleles_for_samples to hold strs instead of lists
    for sample_idx, alleles_for_sample in enumerate(alleles_for_samples):
        if sample_idx in missing_samples:
            continue
        alleles_for_sample_str = []
        for allele_list in alleles_for_sample:
            allele = "".join(allele_list)

            if allele in allele_ids:
                allele_id = allele_ids[allele]
            else:
                allele_id = len(allele_ids)
                allele_ids[allele] = allele_id

            alleles_for_sample_str.append(allele)
        alleles_for_samples[sample_idx] = alleles_for_sample_str

    # create sample and genotypes for each sample
    sample_ids = []
    genotypes = []
    for var_iter_idx, var_iter_info in enumerate(var_iter_infos):
        for sample_idx_in_var_iter, sample in enumerate(var_iter_info.samples):
            sample_idx = (
                sum(n_samples_per_var_iter[:var_iter_idx]) + sample_idx_in_var_iter
            )
            sample_ids.append(sample_idx)
            if sample_idx in missing_samples:
                ploidy = len(alleles_for_samples[sample_idx])
                genotypes.append(tuple(-1 for _ in range(ploidy)))
            else:
                alleles_for_sample = alleles_for_samples[sample_idx]
                sample_gts = tuple(allele_ids[allele] for allele in alleles_for_sample)
                genotypes.append(sample_gts)

    for sample_idx, alleles_for_sample in enumerate(alleles_for_samples):
        print(f"{sample_idx=}")
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
) -> Result[tuple[list[Variant], list[int]], str]:

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
