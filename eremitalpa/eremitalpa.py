"""
Drawing phylogenetic trees (dendropy.Tree instances) using matplotlib.
"""

from operator import attrgetter
import warnings
import dendropy as dp
import matplotlib
import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.collections import LineCollection

# Defaults
default_edge_kws = dict(color="black", linewidth=0.5)
default_leaf_kws = dict(color="black", s=1)
default_internal_kws = dict()
default_label_kws = dict(
    horizontalalignment="left", verticalalignment="center", fontsize=8
)


def compute_tree_layout(tree, has_brlens=True):
    """Compute layout parameters for a tree.

    Each node gets _x and _y values.
    The tree gets _xlim and _ylim values (tuples).

    Args:
        tree (dendropy.Tree)
        has_brlens (bool): Does the tree have branch lengths?

    Returns:
        dendropy.Tree
    """
    # Fresh tree instance
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
    tree,
    has_brlens=True,
    edge_kws=default_edge_kws,
    leaf_kws=default_leaf_kws,
    internal_kws=default_internal_kws,
    ax=None,
    labels=(),
    label_kws=default_label_kws,
    compute_layout=True,
):
    """Plot a dendropy tree object.

    Tree nodes are plotted in their current order. So, to ladderize, call
    tree.ladderize() before plotting.

    Args:
        tree (dendropy.Tree)
        has_brlens (bool). Does the tree have branch lengths? If not, all
            branch lengths are plotted length 1.
        edge_kws (dict). Keyword arguments for edges, passed to
            matplotlib.collections.LineCollection
        leaf_kws (dict). Keyword arguments for leafs, passed to ax.scatter.
            For arguments that can be a vector, the order and length should
            match tree.leaf_node_iter().
        internal_kws (dict). Keyword arguments for internal nodes. Passed to
            ax.scatter. For arguments that can be a vector, the order and
            length should match tree.internal_nodes().
        ax (matplotlib.ax).
        labels (iterable). Names of taxa to add labels for.
        compute_layout (bool). Compute the layout or not. If the tree nodes
            already have _x and _y attributes, then just plot.

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

    lc = matplotlib.collections.LineCollection(segments=edges, **edge_kws)
    ax.add_artist(lc)

    # Draw leaves
    ax.scatter(
        tuple(node._x for node in tree.leaf_node_iter()),
        tuple(node._y for node in tree.leaf_node_iter()),
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
    # Test if labels is an iterable.
    # If it is, label nodes contained in iterable.
    # If it isn't, but is truthy, label all leaf nodes.
    try:
        iter(labels)
    except TypeError:
        if labels:
            for node in tree.leaf_node_iter():
                plt.text(node._x, node._y, node.taxon.label, **label_kws)
        else:
            pass
    else:
        for node in tree.find_nodes(lambda n: taxon_in_node_labels(labels, n)):
            plt.text(node._x, node._y, node.taxon.label, **label_kws)

    # Finalise
    ax.set_xlim(tree._xlim)
    ax.set_ylim(tree._ylim)
    ax.axis("off")
    ax.set_yticks([])
    ax.invert_yaxis()

    return tree, ax


def plot_leaves_with_labels(tree, labels, ax=None, **kws):
    """Plot leaves that have taxon labels in labels.

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
        leaf_sequences = {}
        with open(leaf_seqs, "r") as handle:
            for r in SeqIO.parse(handle, format="fasta"):
                leaf_sequences[r.description] = str(r.seq)

        for node in tree.leaf_node_iter():
            node.sequence = leaf_sequences[node.taxon.label]

    return tree


def sorted_leaf_labels(node):
    """Tuple containing the sorted leaf labels of a node."""
    return tuple(sorted(leaf.taxon.label for leaf in node.leaf_nodes()))


def compare_trees(
    left,
    right,
    gap=0.1,
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
        connect_kws (dict): Keywords passed to matplotlib LineCollection.
            These are used for the lines that connect matching taxa.
        extend_kws (dict): Keywords passed to matplotlib LineCollection.
            These are used for lines that connect taxa to the conection lines.
        extend_every (n): Draw branch extension lines every n leaves.
        left_kws (dict): Passed to plot_tree for the left tree.
        right_kws (dict): Passed to plot_tree for the right tree.
        connect_colors (dict): Maps taxon labels to colors. Ignored if 'colors'
            is used in connect_kws.

    Returns:
        (2-tuple) containing dendropy Trees with _x and _y plot locations on
            nodes.
    """
    left = compute_tree_layout(left)
    right = compute_tree_layout(right)

    # Reflect the right tree
    constant = left._xlim[1] + right._xlim[1] + gap
    for node in right.nodes():
        node._x *= -1
        node._x += constant

    # Cris-crossing lines that connect matching taxa in left and right
    if connect_kws:
        segments = []
        colors = [] if "colors" not in connect_kws and connect_colors else None

        for node in left.leaf_node_iter():
            other = right.find_node_with_taxon_label(node.taxon.label)

            if other:
                segments.append(
                    ((left._xlim[1], node._y), (left._xlim[1] + gap, other._y))
                )

                if colors is not None:
                    colors.append(connect_colors[node.taxon.label])

        if colors is not None:
            connect_kws["colors"] = colors

        plt.gca().add_artist(LineCollection(segments, **connect_kws))

    # Extend branches horizontally from the left and right trees to meet the
    # cris-crossing lines
    if extend_kws:
        segments = []
        colors = [] if "colors" not in extend_kws and extend_colors else None
        key = attrgetter("_y")

        for node in sorted(left.leaf_node_iter(), key=key)[::extend_every]:
            segments.append(((node._x, node._y), (left._xlim[1], node._y)))

            if colors is not None:
                colors.append(extend_colors[node.taxon.label])

        for node in sorted(right.leaf_node_iter(), key=key)[::extend_every]:
            segments.append(((left._xlim[1] + gap, node._y), (node._x, node._y)))

            if colors is not None:
                colors.append(extend_colors[node.taxon.label])

        if colors is not None:
            extend_kws["colors"] = colors

        plt.gca().add_artist(LineCollection(segments, **extend_kws))

    plot_tree(left, compute_layout=False, **left_kws)
    plot_tree(right, compute_layout=False, **right_kws)

    plt.xlim(0, constant)

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
