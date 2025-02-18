import sys
from itertools import groupby
from operator import itemgetter
import argparse
from typing import Generator

from Bio.SeqIO import parse, to_dict


def runs(data: list[int]) -> Generator[str, None, None]:
    """
    Take a list of integers and yield a generator of strings, where each string
    is either a single integer, or a range of integers, depending on whether
    the list of integers is a single element or a run of two or more elements.
    """
    for _, g in groupby(enumerate(data), lambda ix: ix[0] - ix[1]):
        group = list(map(itemgetter(1), g))
        yield str(group[0]) if len(group) == 1 else f"{group[0]}-{group[-1]}"


def main():

    parser = argparse.ArgumentParser(
        "ere_make_chimera",
        description="Make a chimeric sequence sequence characters that are missing in `incomplete` "
        "are filled by those in `donor`. Fasta file read from stdin must contain both incomplete "
        "and donor sequences. The result is printed to stdout, a summary of the sites used from "
        "the donor is printed to stderr.",
    )
    parser.add_argument(
        "--incomplete", help="ID of the incomplete sequence.", required=True
    )
    parser.add_argument("--donor", help="ID of the donor sequence.", required=True)
    parser.add_argument(
        "--chars", help="Characters to replace (default='n-').", default="n-"
    )
    parser.add_argument(
        "--name",
        help="Name to use for the chimera. (default={incomplete}_x_{donor})",
        required=False,
        default=None,
    )
    args = parser.parse_args()

    seqs = to_dict(parse(sys.stdin, "fasta"))

    chimeric_seq = []
    chimeric_sites = []

    chars = set(args.chars)
    if not chars:
        raise ValueError("chars must not be empty")

    seq_donor = seqs[args.donor].lower()

    try:
        seq_incomplete = seqs[args.incomplete].lower()
    except KeyError:
        # If the incomplete sequence is missing, then construct it entirely out of characters that
        # will get replaced by the donor
        seq_incomplete = args.chars[0] * len(seq_donor)

    for i, (a, b) in enumerate(zip(seq_incomplete, seq_donor), start=1):
        if a in chars:
            chimeric_seq.append(b)
            chimeric_sites.append(i)
        elif b in chars:
            raise ValueError(
                f"incomplete and donor both have character from {chars} at site {i}"
            )
        else:
            chimeric_seq.append(a)

    print(
        f"Sites taken from {args.donor}:",
        ",".join(runs(chimeric_sites)) if chimeric_sites else "<none>",
        file=sys.stderr,
    )

    name = f"{args.incomplete}_x_{args.donor}" if args.name is None else args.name
    print(f">{name}")
    print("".join(chimeric_seq))
