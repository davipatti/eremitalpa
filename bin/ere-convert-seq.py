#!/usr/bin/env python3

"""
ere-convert-seq.py

Convert between sequence formats.

Example usage:

    ere-convert-seq.py --in ab1 --out fasta < seq.ab1 > seq.fasta
"""

import argparse
import sys

from Bio import SeqIO

parser = argparse.ArgumentParser("ere-convert-seq.py")
parser.add_argument(
    "-i", "--in", help="Format of input sequences", required=True, dest="in_format"
)
parser.add_argument("-o", "--out", help="Format of output sequences", required=True)
args = parser.parse_args()

for record in SeqIO.parse(sys.stdin.buffer, args.in_format):
    SeqIO.write(record, sys.stdout, format=args.out)
