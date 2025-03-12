import argparse
import os
import sys

from eremitalpa.bio import align_to_reference
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import parse, write
from tqdm import tqdm


def main():
    parser = argparse.ArgumentParser(
        "ere_align_to_reference",
        description="Align FASTA sequences to a reference. Results get printed to stdout. "
        "Sequences that introduce internal gaps in the reference sequence get printed on stderr. "
        "A progress bar gets written to ali-progress.tmp.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        help="Reference fasta sequence. First sequence in this file is taken.",
        required=True,
    )
    parser.add_argument("-s", "--sequences", help="Sequences to align.", required=True)
    args = parser.parse_args()

    with open(args.reference, "r") as fobj:
        reference = next(parse(fobj, "fasta"))

    naughty = []
    nice = []

    # Load everything into a list in order to use tqdm
    with open(args.sequences, "r") as fobj:
        records = list(parse(fobj, "fasta"))

    with open("ali-progress.tmp", "w") as progress:
        for record in tqdm(records, file=progress):

            aligned, internal_gap_in_ref = align_to_reference(reference, record)

            aligned_record = SeqRecord(
                Seq(aligned),
                id=record.id,
                description=record.description,
            )

            if internal_gap_in_ref:
                naughty.append(aligned_record)
            else:
                nice.append(aligned_record)

    write(nice, sys.stdout, "fasta")
    write(naughty, sys.stderr, "fasta")

    os.remove("ali-progress.tmp")
