from sys import stdin, stdout
from functools import partial
import argparse
import Bio
from Bio.SeqIO import parse
from Bio.Align import PairwiseAligner


def alignment_score(a: Bio.SeqRecord.SeqRecord, b: Bio.SeqRecord.SeqRecord) -> float:
    """Alignment score between sequences a and b"""
    return PairwiseAligner().align(a, b).score


def main():
    parser = argparse.ArgumentParser(
        "ere_revcomp_match",
        description="Test if a given sequence or its reverse complement is the best match to a "
        "target. Prints the best match to stdout.",
    )
    parser.add_argument(
        "--target",
        required=True,
        help="Target sequence to match against. Uses the first record in the file.",
    )
    args = parser.parse_args()

    with open(args.target, "r") as handle:
        target = next(parse(handle, format="fasta"))

    score = partial(alignment_score, target)

    for record in parse(stdin, format="fasta"):
        revcomp = record.reverse_complement()
        revcomp.id = record.id  # .reverse_complement() does not copy id or description
        revcomp.description = record.description + " rc"

        best = max(record, revcomp, key=score)

        Bio.SeqIO.write(best, stdout, format="fasta")
