import argparse
from functools import reduce

import dendropy as dp


def main():
    parser = argparse.ArgumentParser(
        "ere_print_phylogeny_leaves",
        description="Print taxon labels in a newick tree. Passing multiple trees prints "
        "either the union (default) or the intersection of their taxon sets.",
    )
    parser.add_argument("newick", help="Path to newick tree(s)", nargs="+")
    parser.add_argument(
        "--intersection",
        help="Show the intersection of taxon sets.",
        default=set.union,
        dest="operate",
        action="store_const",
        const=set.intersection,
    )
    args = parser.parse_args()

    trees = (
        dp.Tree.get(path=path, schema="newick", preserve_underscores=True)
        for path in args.newick
    )

    taxon_labels = reduce(
        args.operate,
        (set(leaf.taxon.label for leaf in tree.leaf_node_iter()) for tree in trees),
    )

    for label in sorted(taxon_labels):
        print(label)
