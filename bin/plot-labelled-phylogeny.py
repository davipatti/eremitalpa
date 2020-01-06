#!/usr/bin/env python3
import argparse

import dendropy
import matplotlib.pyplot as plt

import eremitalpa as ere

parser = argparse.ArgumentParser(description="Plot a labelled phylogeny")
parser.add_argument("--newick", help="Path to newick tree")
parser.add_argument("--filename", help="Path of filename to save")
args = parser.parse_args()
tree = dendropy.Tree.get(path=args.newick, schema="newick")
tree.ladderize()
# Compute a good figure size
tree = ere.compute_tree_layout(tree)
n_leaves = len(tuple(tree.leaf_node_iter()))
height = n_leaves / 50.0
width = 10
# Plot
fig, ax = plt.subplots(figsize=(width, height))
ere.plot_tree(tree, labels=True, label_kws=dict(fontsize=1),
              leaf_kws=dict(s=0), compute_layout=False)
plt.savefig(args.filename, bbox_inches="tight")
