#!/usr/bin/env python3
from eremitalpa import filter_similar_hd
from Bio.SeqIO import parse, write
from sys import stdin, stdout
import argparse

parser = argparse.ArgumentParser(description=
    'Filter sequences that have a hamming distance of less than n to any '
    'already seen. Reads from stdin, writes to stdout.')
parser.add_argument('--n', required=True, type=float)
parser.add_argument('--ignore', default='-X')
parser.add_argument('--case_sensitive', default=False, type=bool)
parser.add_argument('--progress_bar', default=False, type=bool)
args = parser.parse_args()

records = parse(stdin, 'fasta')
subset = filter_similar_hd(records, **vars(args))
write(subset, stdout, 'fasta')
