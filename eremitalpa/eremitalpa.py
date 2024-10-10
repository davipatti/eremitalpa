from collections import namedtuple, Counter, defaultdict
from operator import attrgetter, itemgetter
from typing import Optional, Generator, Iterable, Union, Literal, Any, Mapping
import itertools
import warnings

from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
import dendropy as dp
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bio import amino_acid_colors, sloppy_translate, find_substitutions


# Defaults
default_edge_kws = dict(color="black", linewidth=0.5, clip_on=False, capstyle="round")
default_leaf_kws = dict(s=0)
default_internal_kws = dict()
default_label_kws = dict(
    horizontalalignment="left", verticalalignment="center", fontsize=8
)


def load_fasta(path: str, translate_nt: bool = False) -> dict[str, str]:
    """
    Load fasta file sequences.

    Args:
        path: Path to fasta file.
        translate_nt: Translate nucleotide sequences.
    """
    with open(path) as fobj:
        return {
            record.description: (
                sloppy_translate(str(record.seq)) if translate_nt else str(record.seq)
            )
            for record in SeqIO.parse(fobj, format="fasta")
        }


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
            outgroup: Name of taxon to use as outgroup.
            msa_path: Path to fasta file containing leaf sequences.
            get_kwds: Passed to dendropy.Tree.get.
            kwds: Passed to add_sequences_to_tree
        """
        get_kwds = {} if get_kwds is None else get_kwds

        tree = cls.get(path=path, schema=schema, **get_kwds)

        if msa_path:
            add_sequences_to_tree(tree, path=msa_path, **kwds)

        if outgroup is not None:
            og = tree.find_node_with_taxon_label(outgroup)
            tree.reroot_at_node(og)

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
    edge_kws: dict = default_edge_kws,
    leaf_kws: dict = default_leaf_kws,
    internal_kws: dict = default_internal_kws,
    ax: mp.axes.Axes = None,
    labels: Optional[Union[Iterable[str], Literal["all"]]] = None,
    label_kws: dict = default_label_kws,
    compute_layout: bool = True,
    fill_dotted_lines: bool = False,
) -> mp.axes.Axes:
    """Plot a dendropy tree object.

    Tree nodes are plotted in their current order. So, to ladderize, call
    tree.ladderize() before plotting.

    Args:
        tree
        has_brlens: Does the tree have branch lengths? If not, all
            branch lengths are plotted length 1.
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
        fill_dotted_lines: Show dotted lines from leaves to the right hand edge of the
            tree.

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

    if isinstance(labels, str) and labels == "all":
        labels = [node.taxon.label for node in tree.leaf_nodes()]
    elif labels is None:
        labels = []

    label_kws = {**default_label_kws, **label_kws}
    leaf_kws = {**default_leaf_kws, **leaf_kws}
    edge_kws = {**default_edge_kws, **edge_kws}
    internal_kws = {**default_internal_kws, **internal_kws}

    tree = compute_tree_layout(tree, has_brlens) if compute_layout else tree

    # Draw edges
    edges = []
    for node in tree.preorder_node_iter():
        # Horizontal
        if node.parent_node:
            edges.append(((node._x, node._y), (node.parent_node._x, node._y)))

        # Vertical
        if node.child_nodes():
            maxy = max(node._y for node in node.child_node_iter())
            miny = min(node._y for node in node.child_node_iter())
            edges.append(((node._x, maxy), (node._x, miny)))

    lc = mp.collections.LineCollection(segments=edges, **edge_kws)
    ax.add_artist(lc)

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

    # Draw leaves
    ax.scatter(
        tuple(node._x for node in tree.leaf_node_iter()),
        tuple(node._y for node in tree.leaf_node_iter()),
        marker=[[-2, -1], [-2, 1], [0, 1], [0, -1]],
        clip_on=False,
        **leaf_kws,
    )

    # Draw internal nodes
    if internal_kws:
        ax.scatter(
            tuple(node._x for node in tree.internal_nodes()),
            tuple(node._y for node in tree.internal_nodes()),
            **internal_kws,
        )

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

    # Finalise
    ax.set_xlim(tree._xlim)
    ax.set_ylim(tree._ylim)
    ax.axis("off")
    ax.set_yticks([])
    ax.invert_yaxis()

    return tree, ax


def plot_leaves_with_labels(tree, labels, ax=None, **kws):
    """
    Plot leaves that have taxon labels in labels.

    Args:
        tree (dendropy.Tree)
        labels (iterable): Containing taxon labels to plot.
        ax (mpl ax)
        kws (dict): Passed to plt.scatter
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
    missing_seq_policy: Literal["ignore", "warn", "raise", "ancestor"] = "raise",
    **kwds,
):
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
    """
    ignore = set(ignore_chars)

    def get_seq(node):
        """Return the sequence for the given node."""
        label = get_label(node)
        try:
            return sequences[label]
        except KeyError as err:
            if missing_seq_policy == "raise":
                raise err
            elif missing_seq_policy == "warn":
                warnings.warn(f"no sequence for {label}")
            elif missing_seq_policy == "ignore":
                return None
            elif missing_seq_policy == "ancestor":
                try:
                    return get_seq(node.parent)
                except AttributeError:
                    # Nodes without parents will raise an AttributeError
                    warnings.warn(f"no parent with ancestral sequence for {label}")
                    return None
            else:
                raise ValueError(
                    "missing_seq_policy must be one of 'ignore', 'warn', 'raise' or 'ancestor'."
                )

    for node in tree.nodes():
        if node.parent_node and not (exclude_leaves and node.is_leaf()):
            parent = node.parent_node

            parent_seq = get_seq(parent)
            this_seq = get_seq(node)

            if parent_seq is None or this_seq is None:
                continue

            subs = [
                sub
                for sub in find_substitutions(parent_seq, this_seq, offset=site_offset)
                if all(char not in sub for char in ignore)
            ]

            if len(subs) == 0:
                continue

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


def get_label(node: dp.Node):
    """Return the label of a node. If the node itself has a label, use that. Otherwise
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
    """Ordered nodes in tree, from deepest leaf to root.

    Args:
        tree (dendropy Tree)
        attr (str)

    Returns:
        tuple containin dendropy Nodes
    """
    node = deepest_leaf(tree, attr)
    trunk = []
    while hasattr(node, "parent_node"):
        trunk.append(node)
        node = node.parent_node
    return tuple(trunk)


def deepest_leaf(tree, attr="_x"):
    """Find the deepest leaf node in the tree.

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
) -> dict[str : dict[str, str]]:
    """Read an ancestral state file generated by IQTREE.

    Args:
        state_file: Path to .state file generated by iqtree --ancestral
        partition_names: Partitions are numbered from 1 in the .state file. Pass names for each segment
            (i.e. the order that partition_names appear in the partitions).
        translate_nt: If ancestral states are nucleotide sequences then translate them.

    Returns:
        dict of dicts that maps [node name][partition] -> sequence
    """
    df = pd.read_table(state_file, comment="#")

    seqs = defaultdict(dict)

    for (node_name, part), sub_df in df.groupby(["Node", "Part"], sort=False):
        part_key = partition_names[part - 1] if partition_names is not None else part
        seq = "".join(sub_df["State"])
        seqs[part_key][node_name] = sloppy_translate(seq) if translate_nt else seq

    return seqs


def read_raxml_ancestral_sequences(
    tree, node_labelled_tree, ancestral_seqs, leaf_seqs=None
):
    """Read a tree and ancestral sequences estimated by RAxML.

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
    seqformatter: Optional[callable] = None,
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
        seqformatter
    """
    seqs = {}
    with open(path, "r") as handle:
        for r in SeqIO.parse(handle, format="fasta"):
            key = labeller(r.description) if labeller is not None else r.description
            seq = seqformatter(r.seq) if seqformatter is not None else r.seq
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
            These are used for lines that connect taxa to the conection lines.
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

    # Cris-crossing lines that connect matching taxa in left and right
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
    # cris-crossing lines
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
    """Prune nodes from tree that have a taxon label in labels.

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
        warnings.warn(f"Counldn't find:\n{''.join(not_found)}")

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
    default_color: Optional = None,
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
