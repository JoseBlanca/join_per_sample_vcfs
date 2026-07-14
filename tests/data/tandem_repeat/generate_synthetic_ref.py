#!/usr/bin/env python3
"""Generate the synthetic STR-diversity reference used to validate the tandem-repeat
scanner against the trf-mod catalog (impl plan Milestone D).

Deterministic (fixed LCG seed), so re-running reproduces the committed
`synthetic_ref.fa` byte-for-byte. Two contigs of aperiodic filler with a variety of
embedded STR tracts — perfect and imperfect, periods 2-6, varied copy number — plus
edge cases (a homopolymer, a compound motif, boundary-hugging tracts). trf-mod is run
on the *output* to produce the golden catalog; the scanner must reproduce that catalog's
locus set.

Usage:  python3 generate_synthetic_ref.py > synthetic_ref.fa
"""

MOTIFS_PERFECT = [
    ("AT", 18),
    ("CAG", 14),
    ("AAAG", 9),
    ("GATA", 11),
    ("ATCG", 7),
    ("AAT", 12),
    ("CA", 24),
    ("GAA", 13),
    ("TTTA", 8),
    ("ACGGT", 8),
    ("AACGTT", 7),
    ("TG", 30),
]

# LCG for deterministic aperiodic filler.
_state = 0x1234_5678_9ABC_DEF0


def _rand():
    global _state
    _state = (_state * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
    return _state


def filler(n: int) -> str:
    """`n` bases of aperiodic ACGT (no long short-period tract)."""
    return "".join("ACGT"[(_rand() >> 40) & 3] for _ in range(n))


def tract(motif: str, copies: int) -> str:
    return motif * copies


def substitute(s: str, positions) -> str:
    """Flip the base at each 0-based position to a different one (introduce mismatches)."""
    b = list(s)
    for p in positions:
        b[p] = "ACGT"[("ACGT".index(b[p]) + 1) & 3]
    return "".join(b)


def build_contig() -> str:
    parts = []
    parts.append(filler(120))  # leading flank
    # perfect tracts, each isolated by aperiodic filler
    for motif, copies in MOTIFS_PERFECT:
        parts.append(tract(motif, copies))
        parts.append(filler(90))
    # imperfect tracts: substitutions
    parts.append(substitute(tract("CAG", 16), [18, 33]))  # 2 substitutions
    parts.append(filler(90))
    parts.append(substitute(tract("GATA", 14), [24]))  # 1 substitution
    parts.append(filler(90))
    # imperfect: a single-base insertion mid-tract
    at = tract("AT", 22)
    parts.append(at[:20] + "C" + at[20:])
    parts.append(filler(90))
    # a homopolymer (period-1; dropped downstream) and a compound motif
    parts.append("A" * 25)
    parts.append(filler(90))
    parts.append(tract("ATAT", 8))  # trf may call period 4; compound -> dropped
    parts.append(filler(120))  # trailing flank
    return "".join(parts)


def build_short_contig() -> str:
    # A second, shorter contig with a couple of tracts, incl. one near the start.
    parts = []
    parts.append(filler(70))
    parts.append(tract("TATC", 10))
    parts.append(filler(90))
    parts.append(tract("CCG", 12))
    parts.append(filler(110))
    return "".join(parts)


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main():
    import sys

    c1 = build_contig()
    c2 = build_short_contig()
    sys.stdout.write(">ctg1\n")
    sys.stdout.write(wrap(c1) + "\n")
    sys.stdout.write(">ctg2\n")
    sys.stdout.write(wrap(c2) + "\n")


if __name__ == "__main__":
    main()
