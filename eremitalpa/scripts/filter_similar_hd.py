from eremitalpa import filter_similar_hd
from Bio.SeqIO import parse, write
from sys import stdin, stdout
import argparse


def main():
    parser = argparse.ArgumentParser(
        "ere_filter_similar_hd",
        description="Filter sequences that have a hamming distance of less than n to any "
        "already seen. Reads from stdin, writes to stdout. Case sensitive.",
    )
    parser.add_argument("--n", required=True, type=float)
    parser.add_argument("--ignore", default="-X")
    parser.add_argument("--progress_bar", default=False, type=bool)
    args = parser.parse_args()

    records = parse(stdin, "fasta")
    subset = filter_similar_hd(
        records, n=args.n, progress_bar=args.progress_bar, ignore=set(args.ignore)
    )
    write(subset, stdout, "fasta")
