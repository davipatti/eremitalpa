import argparse

from eremitalpa import read_iqtree_ancestral_states, write_fasta


def main():
    parser = argparse.ArgumentParser(
        "ere_write_iqtree_ancestral_seqs",
        description="Convert an IQTREE .state file containing ancestral sequences to a FASTA file.",
    )
    parser.add_argument("--states", help="Ancestral state file.", required=True)
    parser.add_argument("--fasta", help="Name of fasta file to write.")
    parser.add_argument(
        "--translate",
        help="Output amino acid sequences.",
        action="store_true",
        default=False,
        required=False,
    )
    parser.add_argument(
        "--partition_names",
        nargs="+",
        help=(
            "Optional names for partitions. Names of partitions will be used to name "
            "FASTA files that get dumped in the current directory."
        ),
    )
    args = parser.parse_args()

    seqs = read_iqtree_ancestral_states(
        state_file=args.states, translate_nt=args.translate
    )

    # Check if the states contain partitions
    key = tuple(seqs.keys())[0]

    if isinstance(seqs[key], dict):

        if args.partition_names:

            if len(args.partition_names) != len(seqs):
                raise ValueError(
                    "partition names were passed but a different number of partitions were "
                    "found in the ancestral state file"
                )
            else:
                names = args.partition_names

        else:
            names = range(len(seqs))

        for partition, name in zip(seqs, names):
            write_fasta(f"{name}.ancstate.fasta", seqs[partition])

    else:
        write_fasta(args.fasta, seqs)
