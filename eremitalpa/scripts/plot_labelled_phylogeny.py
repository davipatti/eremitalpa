import argparse

import matplotlib.pyplot as plt
import pandas as pd

import eremitalpa as ere


def main():
    parser = argparse.ArgumentParser(
        "ere_plot_labelled_phylogeny", description="Plot a labelled phylogeny"
    )
    parser.add_argument("-t", "--tree", help="Path to newick tree")
    parser.add_argument("--filename", help="Path of filename to save")
    parser.add_argument(
        "--metadata",
        help="Tab delimited metadata. First column must be node labels used in treefile. "
        "First row must be names of each column.",
        required=False,
    )
    parser.add_argument(
        "--leaf_labels", help="Column in metadata to use as tip labels.", required=False
    )
    parser.add_argument(
        "--append_original_labels",
        help="Append the original label in the tree to the labels when using --leaf_labels. "
        "Sometimes it's useful to keep the ID in the tree.",
        action="store_true",
    )
    parser.add_argument(
        "--leaf_colors",
        help="Column in metadata to use as colors. --leaf_size must be non zero to show "
        "leaves.",
        required=False,
    )
    parser.add_argument(
        "--leaf_size",
        help="Size of circles drawn on leaves. [Default=0].",
        required=False,
        type=int,
        default=0,
    )
    parser.add_argument("--outgroup", required=False, default=None)
    parser.add_argument("--width", required=False, default=None)
    parser.add_argument("--height", required=False, default=None)
    parser.add_argument("--fontsize", required=False, default=None)
    parser.add_argument(
        "--csv", required=False, action="store_true", help="Metadata is a CSV file."
    )
    args = parser.parse_args()

    tree = ere.Tree.from_disk(
        path=args.tree,
        schema="newick",
        get_kwds=dict(preserve_underscores=True),
        outgroup=args.outgroup,
    )
    tree.ladderize()

    if args.leaf_labels and not args.metadata:
        raise ValueError("Must pass metadata")

    if args.metadata:
        df = pd.read_table(args.metadata, sep="," if args.csv else "\t", index_col=0)
        leaves = [leaf.taxon.label for leaf in tree.leaf_nodes()]
    else:
        df = None

    # Labels
    if df is not None and args.leaf_labels:
        leaf_labels = (
            df.loc[leaves, args.leaf_labels]
            if df is not None and args.leaf_labels
            else True
        )
        if args.append_original_labels:
            leaf_labels = [f"{a} {b}" for a, b in zip(leaf_labels, leaves)]

        # Padding before label (space char easier than finding a suitable distance to move
        # text anchor for trees with different scales)
        leaf_labels = [f" {label}" for label in leaf_labels]
    else:
        leaf_labels = True  # Just shows original labels

    leaf_colors = (
        df.loc[leaves, args.leaf_colors]
        if df is not None and args.leaf_colors
        else "black"
    )

    tree = ere.compute_tree_layout(tree)
    height = len(tree.leaf_nodes()) / 25 if args.height is None else int(args.height)
    width = 10 if args.width is None else int(args.width)

    # Plot
    fig, ax = plt.subplots(figsize=(width, height))
    ere.plot_tree(
        tree,
        labels=leaf_labels,
        label_kws=dict(fontsize=2 if args.fontsize is None else int(args.fontsize)),
        leaf_kws=dict(s=args.leaf_size, c=leaf_colors),
        compute_layout=False,
    )
    filename = f"{args.tree}.pdf" if not args.filename else args.filename
    plt.savefig(filename, bbox_inches="tight")
