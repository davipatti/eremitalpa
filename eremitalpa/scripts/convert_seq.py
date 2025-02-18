import argparse
import sys

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(
        "ere_convert_seq", description="Convert between sequence formats"
    )
    parser.add_argument(
        "-i", "--in", help="Format of input sequences", required=True, dest="in_format"
    )
    parser.add_argument("-o", "--out", help="Format of output sequences", required=True)
    args = parser.parse_args()

    for record in SeqIO.parse(sys.stdin.buffer, args.in_format):
        SeqIO.write(record, sys.stdout, format=args.out)
