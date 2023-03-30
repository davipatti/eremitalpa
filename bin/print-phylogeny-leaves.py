#!/usr/bin/env python3
import dendropy
import argparse

parser = argparse.ArgumentParser(
    description="Print labels of leaves in a newick tree to stdout"
)
parser.add_argument("--newick", help="Path to newick tree")
args = parser.parse_args()
tree = dendropy.Tree.get(path=args.newick, schema="newick")
for leaf in tree.leaf_node_iter():
    print(leaf.taxon.label)
