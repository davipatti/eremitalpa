import argparse

from eremitalpa import read_iqtree_ancestral_states, write_fasta


def main():
    parser = argparse.ArgumentParser(
        "ere_write_iqtree_ancestral_seqs",
        description="Convert an IQTREE .state file containing ancestral sequences to a FASTA file.",
    )
    parser.add_argument("--states", help="Ancestral state file.")
    parser.add_argument(
        "--translate",
        help="Output amino acid sequences.",
        action="store_true",
        default=False,
    )
    parser.add_argument("--fasta", help="Name of fasta file to write.")
    args = parser.parse_args()

    seqs = read_iqtree_ancestral_states(
        state_file=args.states, translate_nt=args.translate
    )

    # Check if the states contain partitions
    key = tuple(seqs.keys())[0]
    if isinstance(seqs[key], dict):
        raise NotImplementedError(
            "State file contains different sequences for different"
        )

    else:
        write_fasta(args.fasta, seqs)
