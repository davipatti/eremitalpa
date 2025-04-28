from collections import namedtuple, Counter, defaultdict
from operator import attrgetter, itemgetter
from typing import Optional, Generator, Iterable, Union, Literal, Any, Mapping, Callable
import itertools
import logging
import random
import warnings

from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
import dendropy as dp
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bio import amino_acid_colors, sloppy_translate, find_substitutions

# MARKERS
_DEFAULT_MARKER_KWS = dict(
    color="black", s=0, marker="o", edgecolor="white", lw=0.1, clip_on=False
)
DEFAULT_LEAF_KWS = dict(zorder=15, **_DEFAULT_MARKER_KWS)
DEFAULT_INTERNAL_KWS = dict(zorder=12, **_DEFAULT_MARKER_KWS)

# EDGES
DEFAULT_EDGE_KWS = dict(
    color="black", linewidth=0.5, clip_on=False, capstyle="round", zorder=10
)

# LABELS
DEFAULT_LABEL_KWS = dict(
    horizontalalignment="left", verticalalignment="center", fontsize=8
)


class Tree(dp.Tree):
    def plot_tree_msa(
        self,
        msa_plot_kwds: Optional[dict] = None,
        axes: Optional[tuple[mp.axes.Axes, mp.axes.Axes]] = None,
    ) -> tuple[mp.axes.Axes, mp.axes.Axes]:
        """
        Plot the tree and multiple sequence alignment.
        """
        msa_plot_kwds = {} if msa_plot_kwds is None else msa_plot_kwds

        if axes is None:
            _, axes = plt.subplots(
                ncols=2, figsize=(12, 4), gridspec_kw=dict(wspace=0.7)
            )

        plot_tree(self, ax=axes[0], fill_dotted_lines=True)

        self.multiple_sequence_alignment.plot(
            variable_sites_kwds=msa_plot_kwds.pop(
                "variable_sites_kwds", dict(min_2nd_most_freq=2)
            ),
            rotate_xtick_labels=msa_plot_kwds.pop("rotate_xtick_labels", True),
            **msa_plot_kwds,
            ax=axes[1],
        )

        # Make the ylim of the tree align with the MSA plot
        tree_ylim = axes[0].get_ylim()
        axes[0].set_ylim(tree_ylim[0] + 0.5, tree_ylim[1] - 0.5)

        axes[1].invert_yaxis()

        return axes

    @property
    def multiple_sequence_alignment(self):
        """
        Generate an eremitalpa.MultipleSequence alignment object from a tree. Leaf nodes
        on the tree must have 'sequence' attributes and taxon labels.
        """
        return MultipleSequenceAlignment(
            [
                SeqRecord(node.sequence, description=node.taxon.label)
                for node in self.leaf_nodes()
            ]
        )

    @classmethod
    def from_disk(
        cls,
        path: str,
        schema: str = "newick",
        preserve_underscores: bool = True,
        outgroup: Optional[str] = None,
        msa_path: Optional[str] = None,
        get_kwds: Optional[dict] = None,
        **kwds,
    ) -> "Tree":
        """
        Load a tree from a file.

        Args:
            path: Path to file containing tree.
            schema: See dendropy.Tree.get
            preserve_underscores: Preserve underscores in taxon labels. (Overwrites
                'preserve_underscores' key if passed in get_kwds.)
            outgroup: Name of taxon to use as outgroup.
            msa_path: Path to fasta file containing leaf sequences.
            get_kwds: Passed to dendropy.Tree.get.
            kwds: Passed to add_sequences_to_tree
        """
        get_kwds = {} if get_kwds is None else get_kwds

        get_kwds["preserve_underscores"] = preserve_underscores

        tree = cls.get(path=path, schema=schema, **get_kwds)

        if msa_path:
            add_sequences_to_tree(tree, path=msa_path, **kwds)

        if outgroup is not None:
            og = tree.find_node_with_taxon_label(outgroup)
            tree.reroot_at_edge(og.edge)

        try:
            tree.ladderize(default_order=True)
        except TypeError:
            tree.ladderize()

        return tree


def compute_tree_layout(
    tree: dp.Tree, has_brlens: bool = True, copy: bool = False
) -> dp.Tree:
    """Compute layout parameters for a tree.

    Each node gets _x and _y values.
    The tree gets _xlim and _ylim values (tuples).

    Args:
        tree
        has_brlens: Does the tree have branch lengths?
        copy: Make a fresh copy of the tree.
    """
    if copy:
        tree = dp.Tree(tree)

    # Add branch lengths if necessary
    if not has_brlens:
        for node in tree.preorder_node_iter():
            node.edge.length = 1

    # Compute x for nodes
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            node._x = 0
        else:
            node._x = node.edge.length + node.parent_node._x

    # Compute y for leaf nodes
    for _y, node in enumerate(tree.leaf_node_iter()):
        node._y = _y

    # Compute y for internal nodes
    for node in tree.postorder_node_iter():
        if not hasattr(node, "_y"):
            child_y = tuple(child._y for child in node.child_node_iter())
            node._y = sum(child_y) / len(child_y)

    # X and Y limits
    tree._xlim = 0, max(node._x for node in tree.leaf_nodes())
    tree._ylim = 0, max(node._y for node in tree.leaf_nodes())

    return tree


def plot_tree(
    tree: dp.Tree,
    has_brlens: bool = True,
    edge_kws: dict = DEFAULT_EDGE_KWS,
    leaf_kws: dict = DEFAULT_LEAF_KWS,
    internal_kws: dict = DEFAULT_INTERNAL_KWS,
    ax: mp.axes.Axes = None,
    labels: Optional[Union[Iterable[str], Literal["all"]]] = None,
    label_kws: dict = DEFAULT_LABEL_KWS,
    compute_layout: bool = True,
    fill_dotted_lines: bool = False,
    color_leaves_by_site_aa: Optional[int] = None,
    color_internal_nodes_by_site_aa: Optional[int] = None,
    sequences: Optional[dict[str, str]] = None,
    jitter_x: Optional[float | str] = None,
    scale_bar: Optional[bool] = True,
) -> mp.axes.Axes:
    """
    Plot a dendropy tree object.

    Tree nodes are plotted in their current order. So, to ladderize, call tree.ladderize() before
    plotting.

    Args:
        tree
        has_brlens: Does the tree have branch lengths? If not, all branch lengths are plotted
            length 1.
        edge_kws: Keyword arguments for edges, passed to
            matplotlib.collections.LineCollection
        leaf_kws: Keyword arguments for leafs, passed to ax.scatter.
            For arguments that can be a vector, the order and length should
            match tree.leaf_node_iter().
        label_kwds: Passed to plt.text.
        internal_kws: Keyword arguments for internal nodes. Passed to
            ax.scatter. For arguments that can be a vector, the order and
            length should match tree.internal_nodes().
        ax: Matplotlib ax.
        labels: Taxon labels to annotate, or "all".
        compute_layout. Compute the layout or not. If the tree nodes
            already have _x and _y attributes, then just plot.
        fill_dotted_lines: Show dotted lines from leaves to the right hand edge of the tree.
        color_leaves_by_site_aa: Pass an integer to color the leaves by the amino acid at this site
            (1-based). This will overwrite the 'c' kwarg in leaf_kws. `sequences` must be passed.
        color_internal_nodes_by_site_aa: Same behaviour as color_leaves_by_site_aa but for internal
            nodes.
        sequences: A mapping of taxon labels and to sequences. Required for
            `color_leaves_by_site_aa`.
        jitter_x: Add a small amount of noise to the x value of the leaves to avoid over plotting.
            Either pass a float (the amount of noise) or 'auto' to try to automatically calculate a
            suitable value. 'auto' tries to calculate the fundamental 'unit' of branch length in the
            tree and then jitters x values by 1/2 of this value in either direction. See
            estimate_unit_branch_length for more information. Currently, positions of labels are
            not jittered.
        scale_bar: Show a scale bar at the bottom of the tree.

    Returns:
        tuple containing (Tree, ax). The tree and matplotlib ax. The tree has
            these additional attributes:

                _xlim (tuple) Min and max x value of nodes.
                _ylim (tuple) Min and max y value of nodes.

            Each node has these attributes:

                _x (number) X-coordinate of the nodes layout
                _y (number) Y-coordinate of the node's layout

    """
    ax = plt.gca() if ax is None else ax

    if labels == "all":
        labels = [node.taxon.label for node in tree.leaf_nodes()]

    elif labels is None:
        labels = []

    label_kws = {**DEFAULT_LABEL_KWS, **label_kws}
    leaf_kws = {**DEFAULT_LEAF_KWS, **leaf_kws}
    edge_kws = {**DEFAULT_EDGE_KWS, **edge_kws}
    internal_kws = {**DEFAULT_INTERNAL_KWS, **internal_kws}

    tree = compute_tree_layout(tree, has_brlens) if compute_layout else tree

    # Draw edges
    edges = []
    for node in tree.preorder_node_iter():
        # Horizontal
        if node.parent_node:
            edges.append(((node._x, node._y), (node.parent_node._x, node._y)))

        # Vertical
        if node.child_nodes():
            max_y = max(node._y for node in node.child_node_iter())
            min_y = min(node._y for node in node.child_node_iter())
            edges.append(((node._x, max_y), (node._x, min_y)))

    ax.add_artist(mp.collections.LineCollection(segments=edges, **edge_kws))

    # Dotted lines from the leaves to the right hand edge of the tree
    if fill_dotted_lines:
        max_x = max(node._x for node in tree.leaf_nodes())
        dotted_edges = [
            ((node._x, node._y), (max_x, node._y)) for node in tree.leaf_nodes()
        ]
        ax.add_artist(
            mp.collections.LineCollection(
                segments=dotted_edges, ls=(2, (1, 10)), color="black", linewidth=0.5
            )
        )

    # Infer suitable jitter_x value if need be
    if jitter_x == "auto":
        jitter_x = estimate_unit_branch_length(
            [edge.length for edge in tree.edges() if edge.length is not None]
        )
        logging.info(f"Auto jitter_x: {jitter_x}")

    # Draw leaves
    if color_leaves_by_site_aa is not None:

        # Group leaves by amino acid at site so that each group can be colored and a label passed
        # for the legend.

        def _get_aa(node):
            """
            Temporary helper function to get the amino acid at the site. Passed to sorted and
            itertools.groupby to group leaves regardless of input order.
            """
            return sequences[node.taxon.label][color_leaves_by_site_aa - 1]

        # Make groups of nodes that all have a particular amino acid at the site
        sorted_nodes = sorted(tree.leaf_node_iter(), key=_get_aa)
        aa_groups = {
            aa: list(nodes)
            for aa, nodes in itertools.groupby(sorted_nodes, key=_get_aa)
        }

        # Order the groups by size so that the smallest groups are plotted last to make them more
        # visible. Put unknown amino acids at the back
        for aa in reversed(
            sorted(
                aa_groups,
                key=lambda aa: len(aa_groups[aa]) if aa != "X" else 0,
            )
        ):
            nodes = aa_groups[aa]
            leaf_kws["color"] = amino_acid_colors[aa]
            x, y = node_x_y(nodes, jitter_x=jitter_x)
            ax.scatter(x, y, **leaf_kws, label=aa if aa != "X" else None)

    else:
        x, y = node_x_y(tree.leaf_node_iter(), jitter_x=jitter_x)
        ax.scatter(x, y, **leaf_kws)

    # Draw internal nodes
    if color_internal_nodes_by_site_aa is not None:
        internal_kws["color"] = [
            amino_acid_colors[
                sequences[node.label][color_internal_nodes_by_site_aa - 1]
            ]
            for node in tree.internal_nodes()
        ]

    if internal_kws:
        x, y = node_x_y(tree.internal_nodes())
        ax.scatter(x, y, **internal_kws)

    # Labels

    # If labels is True but not iterable, simply label all leaf nodes
    if not isinstance(labels, Iterable) and labels:
        for node in tree.leaf_node_iter():
            ax.text(node._x, node._y, node.taxon.label, **label_kws)

    # If labels is a mapping then look up the label for each node
    elif isinstance(labels, (Mapping, pd.Series)):
        for node in tree.leaf_node_iter():
            if label := labels.get(node.taxon.label):
                ax.text(node._x, node._y, label, **label_kws)

    # If all nodes are passed, plot all their labels
    elif all(isinstance(item, dp.Node) for item in labels):
        for node in labels:
            ax.text(node._x, node._y, node.taxon.label, **label_kws)

    elif all(isinstance(item, str) for item in labels):

        # If all strings are passed, and there is one per leaf, plot each on a leaf
        if len(labels) == len(tree.leaf_nodes()):
            for node, label in zip(tree.leaf_node_iter(), labels):
                ax.text(node._x, node._y, label, **label_kws)

        # If all strings are passed, and there are fewer than one per leaf, find
        # the nodes that have these taxon labels and label them
        elif len(labels) < len(tree.leaf_nodes()):
            for node in tree.find_nodes(lambda n: taxon_in_node_labels(labels, n)):
                ax.text(
                    node._x,
                    node._y,
                    node.taxon.label,
                    **label_kws,
                )

        else:
            raise ValueError("passed more labels than number of leaf nodes")

    else:
        raise ValueError("couldn't process labels")

    if scale_bar:
        length = tree._xlim[1] / 10
        length = float(f"{length:.1g}")  # round length to 1 significant figure
        bottom = tree._ylim[1]
        ax.plot((0, length), (bottom, bottom), c="black", lw=1, clip_on=False)
        ax.text(length / 2, bottom, str(length), ha="center", va="bottom")

    # Finalise
    ax.set_xlim(tree._xlim)
    ax.set_ylim(tree._ylim)
    ax.axis("off")
    ax.set_yticks([])
    ax.invert_yaxis()

    return tree, ax


def node_x_y(
    nodes: Iterable[dp.Node], jitter_x: Optional[float] = None
) -> tuple[tuple, tuple]:
    """
    x and y coordinates of nodes.

    Args:
        nodes (Iterable[dp.Node]): An iterable collection of dp.Node objects.
        jitter_x (Optional[float]): The amount of jitter to add to x coordinates. X is jittered
            by a quarter of this value above and below.

    Returns:
        tuple[tuple, tuple]: A tuple containing two tuples, the first with all x coordinates and the
            second with all y coordinates.
    """
    if jitter_x is None:
        return zip(*((node._x, node._y) for node in nodes))
    else:
        lo = -jitter_x / 4
        hi = jitter_x / 4
        return zip(*((node._x + random.uniform(lo, hi), node._y) for node in nodes))


def node_x_y_from_taxon_label(tree: Tree, taxon_label: str) -> tuple[float, float]:
    """
    Find the x and y attributes of a node in a tree from a taxon label.

    Args:
        tree: Tree
        taxon_label: str

    Returns:
        tuple[float, float]
    """
    node = tree.find_node_with_taxon_label(taxon_label)
    return node._x, node._y


def plot_leaves_with_labels(
    tree: dp.Tree, labels: list[str], ax: mp.axes.Axes = None, **kws
):
    """
    Plot leaves that have taxon labels in labels.

    Args:
        tree
        labels: Taxon labels to plot.
        ax: Matplotlib ax
        **kws: Passed to plt.scatter
    """
    ax = plt.gca() if ax is None else ax
    s = kws.pop("s", 20)
    c = kws.pop("c", "red")
    zorder = kws.pop("zorder", 19)
    linewidth = kws.pop("linewidth", 0.5)
    edgecolor = kws.pop("edgecolor", "white")
    nodes = tree.find_nodes(lambda n: taxon_in_node_labels(labels, n))

    if not nodes:
        raise ValueError("No node with taxon labels in labels found in tree.")

    try:
        x = [node._x for node in nodes]
    except AttributeError as err:
        print("Node(s) do not have _x attribute. Run compute_tree_layout.")
        raise (err)

    try:
        y = [node._y for node in nodes]
    except AttributeError as err:
        print("Node(s) do not have _y attribute. Run compute_tree_layout.")
        raise (err)

    ax.scatter(
        x, y, s=s, c=c, zorder=zorder, linewidth=linewidth, edgecolor=edgecolor, **kws
    )


def plot_subs_on_tree(
    tree: dp.Tree,
    sequences: dict[str, str],
    exclude_leaves: bool = True,
    site_offset: int = 0,
    ignore_chars: str = "X-",
    arrow_length: float = 40,
    arrow_facecolor: str = "black",
    fontsize: float = 6,
    **kwds,
) -> Counter:
    """
    Plot substitutions on a tree. This function plots substitutions on the tree by finding
    substitutions between each node and its parent node. The substitutions are then plotted at the
    midpoint of the edge between the node and its parent node.

    Args:
        tree (dendropy.Tree): The tree to annotate.
        sequences (dict[str, str]): A mapping of node labels to sequences.
        exclude_leaves (bool): If True, exclude leaves from getting substitutions.
        site_offset (int): Value added to substitution sites. E.g. if site '1' is actually at
            index 16 in the sequences, then pass 16.
        ignore_chars (str): Substitutions involving characters in this string will not be shown in
            substitutions.
        arrow_length (float): The length of the arrow pointing to the mutation.
        arrow_facecolor (str): The facecolor of the arrow pointing to the mutation.
        fontsize (float): The fontsize of the text.
        **kwds: Other keyword arguments to pass to plt.annotate.

    Returns:
        Counter containing the number of times each substitution appears in the tree.
    """
    ignore = set(ignore_chars)

    sub_counts = Counter()

    for node in tree.nodes():
        if (parent := node.parent_node) and not (exclude_leaves and node.is_leaf()):

            parent_seq = sequences[get_label(parent)]
            this_seq = sequences[get_label(node)]

            if parent_seq is None or this_seq is None:
                continue

            subs = [
                sub
                for sub in find_substitutions(parent_seq, this_seq, offset=site_offset)
                if all(char not in sub for char in ignore)
            ]

            if len(subs) == 0:
                continue

            sub_counts.update(subs)

            x = (node._x + parent._x) / 2

            plt.annotate(
                "\n".join(map(str, subs)),
                (x, node._y),
                xytext=(-arrow_length, arrow_length),
                va="bottom",
                ha="right",
                textcoords="offset pixels",
                arrowprops=dict(
                    facecolor=arrow_facecolor,
                    shrink=0,
                    linewidth=0,
                    width=0.3,
                    headwidth=2,
                    headlength=2,
                ),
                fontsize=fontsize,
                **kwds,
            )

    return sub_counts


def get_label(node: dp.Node):
    """
    Return the label of a node. If the node itself has a label, use that. Otherwise
    return the label of the node's taxon.
    """
    if node.label is not None:
        return node.label
    else:
        return node.taxon.label


def taxon_in_node_labels(labels, node):
    """True if node has taxon label in labels, else False"""
    try:
        return node.taxon.label in labels
    except AttributeError:
        return False


def taxon_in_node_label(label, node):
    """True if a node has a matching taxon label"""
    try:
        return node.taxon.label == label
    except AttributeError:
        return False


def get_trunk(tree, attr="_x"):
    """
    Ordered nodes in tree, from deepest leaf to root.

    Args:
        tree (dendropy Tree)
        attr (str)

    Returns:
        tuple containing dendropy Nodes
    """
    node = deepest_leaf(tree, attr)
    trunk = []
    while hasattr(node, "parent_node"):
        trunk.append(node)
        node = node.parent_node
    return tuple(trunk)


def deepest_leaf(tree, attr="_x"):
    """
    Find the deepest leaf node in the tree.

    Args:
        tree (dendropy Tree)
        attr (str): Either _x or _y. Gets node with max attribute.

    Returns:
        dendropy Node
    """
    try:
        return max(tree.leaf_node_iter(), key=attrgetter(attr))

    except AttributeError:
        tree = compute_tree_layout(tree)
        return max(tree.leaf_node_iter(), key=attrgetter(attr))


def read_iqtree_ancestral_states(
    state_file, partition_names: Optional[list[str]] = None, translate_nt: bool = False
) -> dict[str : dict[str, str]] | dict[str:str]:
    """
    Read an ancestral state file generated by IQ-TREE. If the file contains multiple partitions
    (i.e. a 'Part' column is present), then return a dict of dicts containing sequences accessed by
    [partition][node]. Otherwise return a dict of sequences accessed by node.

    Args:
        state_file: Path to .state file generated by iqtree --ancestral
        partition_names: Partitions are numbered from 1 in the .state file. Pass names for each
            segment (i.e. the order that partition_names appear in the partitions). Only takes
            effect if multiple partitions are present.
        translate_nt: If ancestral states are nucleotide sequences then translate them.

    Returns:
        dict of dicts that maps [node][partition] -> sequence, or dict that maps node -> sequence.
    """
    df = pd.read_table(state_file, comment="#")

    if "Part" in df:
        seqs = defaultdict(dict)
        for (node_name, part), sub_df in df.groupby(["Node", "Part"], sort=False):
            part_key = (
                partition_names[part - 1] if partition_names is not None else part
            )
            seq = "".join(sub_df["State"])
            seqs[part_key][node_name] = sloppy_translate(seq) if translate_nt else seq

    else:
        seqs = {}
        for node_name, sub_df in df.groupby("Node", sort=False):
            seq = "".join(sub_df["State"])
            seqs[node_name] = sloppy_translate(seq) if translate_nt else seq

    return seqs


def read_raxml_ancestral_sequences(
    tree, node_labelled_tree, ancestral_seqs, leaf_seqs=None
):
    """
    Read a tree and ancestral sequences estimated by RAxML.

    RAxML can estimate marginal ancestral sequences for internal nodes on a
    tree using a call like:

        raxmlHPC -f A -t {treeFile} -s {sequenceFile} -m {model} -n {name}

    The analysis outputs several files:

    - RAxML_nodeLabelledRootedTree.{name} contains a copy of the input tree
        where all internal nodes have a unique identifier {id}.
    - RAxML_marginalAncestralStates.{name} contains the ancestral sequence for
        each internal node. The format of each line is '{id} {sequence}'
    - RAxML_marginalAncestralProbabilities.{name} contains probabilities of
        each base at each site for each internal node. (Not used by this
        function.)

    Notes:
        Developed with output from RAxML version 8.2.12.

    Args:
        tree (str): Path to original input tree ({treeFile}).
        node_labelled_tree (str): Path to the tree with node labels.
            (RAxML_nodeLabelledRootedTree.{name})
        ancestral_seqs (str): Path to file containing the ancestral sequences.
            (RAxML_marginalAncestralStates.{name})
        leaf_seqs (str): (Optional) path to fasta file containing leaf
            sequences. ({sequenceFile}). If this is provided, also attach
            sequences to leaf nodes.

    Returns:
        (dendropy Tree) with sequences attached to nodes. Sequences are
            attached as 'sequence' attributes on Nodes.
    """
    tree = dp.Tree.get(path=tree, schema="newick", preserve_underscores=True)
    labelled_tree = dp.Tree.get(
        path=node_labelled_tree, schema="newick", preserve_underscores=True
    )

    # Dict mapping leaf labels -> node label
    leaves_to_labelled_node = {
        sorted_leaf_labels(node): node.label for node in labelled_tree.nodes()
    }

    internal_sequences = {}
    with open(ancestral_seqs, "r") as handle:
        for i, line in enumerate(handle.readlines()):
            try:
                key, sequence = line.strip().split()
            except ValueError as err:
                print(
                    f"Problem reading sequence on line {i + 1} in " f"{ancestral_seqs}."
                )
                raise err
            internal_sequences[key] = sequence

    for node in tree.internal_nodes():
        leaves = sorted_leaf_labels(node)
        key = leaves_to_labelled_node[leaves]
        node.sequence = internal_sequences[key]

    if leaf_seqs:
        add_sequences_to_tree(tree, leaf_seqs)

    return tree


def add_sequences_to_tree(
    tree: dp.Tree,
    path: str,
    labeller: Optional[callable] = None,
    seq_formatter: Optional[callable] = None,
) -> None:
    """
    Add sequences to leaves inplace.

    Args:
        tree: Dendropy tree.
        path: Path to multiple sequence alignment. Taxon labels in tree must match the
            fasta description.
        labeller: Function that takes a FASTA description and returns the name of
            associated taxon in the tree. For instance, a full fasta description might
            look like:

                >cdsUUB77424 A/Maryland/12786/2022 2022/04/07 HA

            RAxML would call this taxon just 'cdsUUB77424'. So the callable would have to
            be something like: lambda x: x.split()[0]
        seq_formatter
    """
    seqs = {}
    with open(path, "r") as handle:
        for r in SeqIO.parse(handle, format="fasta"):
            key = labeller(r.description) if labeller is not None else r.description
            seq = seq_formatter(r.seq) if seq_formatter is not None else r.seq
            seqs[key] = seq

    for node in tree.leaf_node_iter():
        node.sequence = seqs[node.taxon.label]


def sorted_leaf_labels(node):
    """Tuple containing the sorted leaf labels of a node."""
    return tuple(sorted(leaf.taxon.label for leaf in node.leaf_nodes()))


def compare_trees(
    left,
    right,
    gap=0.1,
    x0=0,
    connect_kws=dict(),
    extend_kws=dict(),
    extend_every=10,
    left_kws=dict(),
    right_kws=dict(),
    connect_colors=dict(),
    extend_colors=dict(),
):
    """Plot two phylogenies side by side, and join the same taxa in each tree.

    Args:
        left (dendropy Tree)
        right (dendropy Tree)
        gap (float): Space between the two trees.
        x0 (float): The x coordinate of the root of the left hand tree.
        connect_kws (dict): Keywords passed to matplotlib LineCollection.
            These are used for the lines that connect matching taxa.
        extend_kws (dict): Keywords passed to matplotlib LineCollection.
            These are used for lines that connect taxa to the connection lines.
        extend_every (n): Draw branch extension lines every n leaves.
        left_kws (dict): Passed to plot_tree for the left tree.
        right_kws (dict): Passed to plot_tree for the right tree.
        connect_colors (dict or Callable): Maps taxon labels to colors. Ignored if
            'colors' is used in connect_kws.
        extend_colors (dict or Callable): Maps taxon labels to colors. Ignored if
            'colors' is used in extend_kws.

    Returns:
        (2-tuple) containing dendropy Trees with _x and _y plot locations on
            nodes.
    """
    left = compute_tree_layout(left)
    right = compute_tree_layout(right)

    # Reflect the right tree
    constant = left._xlim[1] + right._xlim[1] + gap + x0
    for node in right.nodes():
        node._x *= -1
        node._x += constant

    # # Move the left tree by x0
    for node in left.nodes():
        node._x += x0

    # Criss-crossing lines that connect matching taxa in left and right
    if connect_kws:
        segments = []
        colors = [] if "colors" not in connect_kws and connect_colors else None

        for node in left.leaf_node_iter():
            other = right.find_node_with_taxon_label(node.taxon.label)

            if other:
                segments.append(
                    (
                        (left._xlim[1] + x0, node._y),
                        (left._xlim[1] + x0 + gap, other._y),
                    )
                )

                if colors is not None:
                    try:
                        c = connect_colors[node.taxon.label]
                    except TypeError:
                        c = connect_colors(node.taxon.label)
                    colors.append(c)

        if colors is not None:
            connect_kws["colors"] = colors

        plt.gca().add_artist(mp.collections.LineCollection(segments, **connect_kws))

    # Extend branches horizontally from the left and right trees to meet the
    # criss-crossing lines
    if extend_kws:
        segments = []
        colors = [] if "colors" not in extend_kws and extend_colors else None
        key = attrgetter("_y")

        for node in sorted(left.leaf_node_iter(), key=key)[::extend_every]:
            segments.append(((node._x, node._y), (left._xlim[1] + x0, node._y)))

            if colors is not None:
                try:
                    c = extend_colors[node.taxon.label]
                except TypeError:
                    c = extend_colors(node.taxon.label)
                colors.append(c)

        for node in sorted(right.leaf_node_iter(), key=key)[::extend_every]:
            segments.append(((left._xlim[1] + x0 + gap, node._y), (node._x, node._y)))

            if colors is not None:
                try:
                    c = extend_colors[node.taxon.label]
                except TypeError:
                    c = extend_colors(node.taxon.label)
                colors.append(c)

        if colors is not None:
            extend_kws["colors"] = colors

        plt.gca().add_artist(mp.collections.LineCollection(segments, **extend_kws))

    plot_tree(left, compute_layout=False, **left_kws)
    plot_tree(right, compute_layout=False, **right_kws)

    # plt.xlim(0, constant)

    return left, right


def prune_nodes_with_labels(tree, labels):
    """
    Prune nodes from tree that have a taxon label in labels.

    Args:
        tree (dendropy Tree)
        labels (iterable containing str)

    Returns:
        (dendropy Tree)
    """
    nodes = []
    not_found = []
    for label in labels:
        node = tree.find_node_with_taxon_label(label)
        if node is None:
            not_found.append(node)
        else:
            nodes.append(node)

    if not_found and len(nodes) == 0:
        raise ValueError("No taxa found with any of these labels.")
    elif not_found:
        warnings.warn(f"Couldn't find:\n{''.join(not_found)}")

    tree.prune_nodes(nodes)
    tree.prune_leaves_without_taxa(recursive=True)
    return tree


Column = namedtuple("Column", ["site", "aas"])


class MultipleSequenceAlignment(Align.MultipleSeqAlignment):
    def variable_sites(
        self, min_2nd_most_freq: int = 1
    ) -> Generator[Column, None, None]:
        """
        Generator for variable sites in the alignment.

        Args:
            min_2nd_most_freq: Used to filter out sites that have low variability. For
                instance if min_2nd_most_freq is 2 a column containing 'AAAAT' should be
                excluded because the second most frequent character (T) has a frequency
                of 1.
        """
        for i in range(self.get_alignment_length()):
            site = i + 1
            aas = self[:, i]
            if len(set(aas)) != 1:
                # Frequency of the second most common
                _, count = Counter(aas).most_common()[1]
                if count >= min_2nd_most_freq:
                    yield Column(site, aas)

    def plot(
        self,
        ax: Optional[mp.axes.Axes] = None,
        fontsize: int = 6,
        variable_sites_kwds: Optional[dict] = None,
        rotate_xtick_labels: bool = False,
        sites: Optional[Iterable[int]] = None,
    ) -> mp.axes.Axes:
        """
        Plot variable sites in the alignment.

        Args:
            ax: Matplotlib ax.
            fontsize: Fontsize of the character labels.
            variable_sites_kwds: Passed to MultipleSequenceAlignment.variable_sites.
            rotate_xtick_labels: Rotate the xtick labels 90 degrees.
            sites: Only plot these sites. (Note: Only variable sites are plotted, so if a
                site is passed in this argument but it is not variable it will not be
                displayed.)
        """
        ax = plt.gca() if ax is None else ax
        variable_sites_kwds = {} if variable_sites_kwds is None else variable_sites_kwds

        variable_sites = tuple(self.variable_sites(**variable_sites_kwds))

        if len(variable_sites) == 0:
            raise ValueError("No variable sites in alignment")

        if sites is not None:
            sites = set(sites)
            variable_sites = tuple(
                column for column in variable_sites if column.site in sites
            )

        for x, site in enumerate(variable_sites):
            for y, aa in enumerate(site.aas):
                rect = mp.patches.Rectangle(
                    (x, y), width=1, height=1, facecolor=amino_acid_colors[aa]
                )
                ax.add_artist(rect)
                ax.text(
                    x + 0.5,
                    y + 0.5,
                    aa,
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    color="black" if aa.upper() in "M" else "white",
                )

        max_x = rect.xy[0]

        ax.set(
            xlim=(0, max_x + 1),
            ylim=(0, len(self)),
            yticks=np.arange(0.5, len(self) + 0.5),
            yticklabels=[record.description for record in self],
            xticks=np.arange(0.5, max_x + 1.5),
            xticklabels=[column.site for column in variable_sites],
        )

        if rotate_xtick_labels:
            ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=90)

        for spine in "top", "bottom", "left", "right":
            ax.spines[spine].set_visible(False)

        return ax


def color_stack(
    tree: Tree,
    values: dict[str, Any],
    color_dict: dict[str, str],
    default_color: Optional[str] = None,
    x: float = 0,
    ax: Optional[mp.axes.Axes] = None,
    leg_kwds: Optional[dict] = None,
) -> tuple[mp.axes.Axes, mp.legend.Legend]:
    """
    A stack of colored patches that can be plotted adjacent to a tree to show how values
    vary on the tree leaves.

    Must have called eremitalpa.compute_layout on the tree in order to know y values for
    leaves (done anyway by eremitalpa.plot_tree).

    Args:
        tree: The tree to be plotted next to.
        values: Maps taxon labels to values to be plotted.
        color_dict: Maps values to colors.
        default_color: Color to use for values missing from color_dict.
        x: The x value to plot the stack at.
        ax: Matplotlib ax
    """
    ax = ax or plt.gca()
    leg_kwds = leg_kwds or dict()

    labels = [leaf.taxon.label for leaf in tree.leaf_nodes()]

    leaf_ys = [leaf._y for leaf in tree.leaf_nodes()]

    colors = [color_dict.get(values[label], default_color) for label in labels]

    # Group color and leaf y values by batches of the same color in order to make a
    # single larger patch if consecutive patches would be the same color.
    for _, grouped in itertools.groupby(zip(colors, leaf_ys), key=itemgetter(0)):

        # color will all be the same (it is what is being grouped by)
        # leaf_y will be the y values of each leaf in this group
        color, leaf_y = zip(*grouped)

        ax.add_patch(
            mp.patches.Rectangle(
                (x, leaf_y[0] - 0.5),  # bottom of the patch is the first y value
                width=1,
                height=len(color),  # height of the patch is just the size of the group
                color=color[0],
                linewidth=0,
            )
        )

    # Patches off the ax for easy legend
    handles_labels = [(mp.patches.Patch(color=v), k) for k, v in color_dict.items()]
    handles, labels = zip(*handles_labels)
    leg = ax.legend(handles, labels, **leg_kwds)
    ax.add_artist(leg)
    ax.axis(False)
    return ax, leg


def estimate_unit_branch_length(
    branch_lengths: list[float], min_diff: float = 1e-6, round_to: int = 6
) -> float:
    """
    Estimates the fundamental unit length in a set of phylogenetic branch lengths. Assumes that
    branch lengths occur in approximate integer multiples of a small unit. Algorithm:

        1. Compute all pairwise absolute differences between branch lengths.
        2. Construct a histogram of these differences with an adaptive bin size.
        3. Identify the most common small difference (mode of the histogram),
        which represents the estimated unit length.

    Args:
        branch_lengths (list or np.array): A list of branch lengths.
        min_diff (float): Minimum difference between branch lengths to consider. Branch length
            differences smaller than this value are considered to be zero.
        round_to (int): Round branch_lengths to this many decimal places.

    Returns:
        float: Estimated fundamental unit length.
    """
    branch_lengths = np.array([round(edge, round_to) for edge in branch_lengths])

    # Compute pairwise absolute differences
    diffs = np.abs(np.subtract.outer(branch_lengths, branch_lengths))
    diffs = diffs[np.triu_indices_from(diffs, k=1)]  # Extract unique values
    diffs = diffs[diffs > min_diff]  # Remove near-zero values

    # Define bin width adaptively based on small quantile
    bin_width = np.percentile(diffs, 1)
    hist, bin_edges = np.histogram(diffs, bins=np.arange(0, np.max(diffs), bin_width))

    # Find bin with the highest count (mode) and return it
    return bin_edges[np.argmax(hist)]


def plot_tree_with_subplots(
    tree: Tree,
    aa_seqs: dict,
    site: int,
    subplot_taxa_shifts: dict[str, tuple[float, float]],
    fun: Callable,
    fun_kwds: Optional[dict] = None,
    subplot_width: float = 0.2,
    subplot_height: float = 0.1,
    figsize: tuple[float, float] = (8, 12),
    sharex: bool = True,
    sharey: bool = True,
) -> None:
    """
    Plot a phylogeny tree with subplots for specified taxa.

    This function draws a phylogeny based on a given tree and amino acid sequences.
    It colors leaves (and internal nodes) according to their amino acid at a specified site,
    and attaches additional subplots at user-defined nodes for further custom visualization.

    Args:
        tree: eremitalpa.Tree
            The phylogenetic tree to be plotted.
        aa_seqs: dict
            A dictionary containing amino acid sequences for each taxon.
            Keys should match the node names in the tree, and values should be the sequences.
        site: int
            The site (1-based) to color the tree's leaves and internal nodes.
        subplot_taxa_shifts: dict of str -> tuple of float
            A mapping from taxon names to tuples (x_shift, y_shift).
            These values control the position of the subplot axes relative to their respective nodes.
        fun: Callable
            A callable function to generate each subplot.
            Must accept the current taxon as the first argument and an axes object as the second argument.
        fun_kwds: dict
            A dictionary of additional keyword arguments passed to the subplot function ``fun``.
        subplot_width: float, optional
            The width of each subplot in figure coordinates, by default 0.2.
        subplot_height: float, optional
            The height of each subplot in figure coordinates, by default 0.1.
        figsize: tuple of float, optional
            The overall size of the figure, by default (8, 12).
        sharex: bool, Have the sub axes share x-axes.
        sharey: bool, Have the sub axes share y-axes.

    Returns:
        None
            The function draws and displays a matplotlib figure with the main phylogeny
            and subplots at specified taxa. It does not return any objects.
    """

    fig, main_ax = plt.subplots(figsize=figsize)
    plot_tree(
        tree,
        color_leaves_by_site_aa=site,
        color_internal_nodes_by_site_aa=site,
        sequences=aa_seqs,
        leaf_kws=dict(s=8),
        internal_kws=dict(s=2),
        edge_kws=dict(lw=0.5, color="grey"),
        jitter_x="auto",
        ax=main_ax,
    )
    main_ax.legend(
        markerscale=5, loc="lower left", bbox_to_anchor=(1.05, 0), title=f"Site {site}"
    )

    data_to_fig = (main_ax.transData + fig.transFigure.inverted()).transform

    first_sub_ax = None

    for taxon in subplot_taxa_shifts:
        x, y = node_x_y_from_taxon_label(tree, taxon)

        x_shift, y_shift = subplot_taxa_shifts.get(taxon, (0.1, 0.1))
        x_shifted = x + tree._xlim[1] * x_shift
        y_shifted = y + tree._ylim[1] * -y_shift

        left, bottom = data_to_fig([x_shifted, y_shifted])

        # Lines connecting root viruses to their sub axes
        main_ax.plot(
            (x, x_shifted), (y, y_shifted), c="black", lw=0.5, clip_on=False, zorder=16
        )

        # Big marker showing where the root virus is
        main_ax.scatter(
            x, y, c="black", s=40, clip_on=False, marker="*", zorder=16, lw=0
        )

        # Add the sub axes
        position = left, bottom, subplot_width, subplot_height
        if first_sub_ax is None:
            sub_ax = plt.axes(position)
            first_sub_ax = sub_ax
        else:
            sub_ax = plt.axes(
                position,
                sharex=first_sub_ax if sharex else False,
                sharey=first_sub_ax if sharey else False,
            )

        fun_kwds = {} if fun_kwds is None else fun_kwds
        fun(taxon, ax=sub_ax, **fun_kwds)
