"""
Drawing phylogenetic trees (dendropy.Tree instances) using matplotlib.
"""

import dendropy
import matplotlib
import matplotlib.pyplot as plt

# Defaults
edge_kws = dict(
    color="black",
    linewidth=0.5
)
leaf_kws = dict(
    color="black",
    s=1
)
label_kws = dict(
    horizontalalignment="left",
    verticalalignment="center",
    fontsize=8
)


def plot_dendropy_tree(tree, has_brlens=True, edge_kws=edge_kws,
                       leaf_kws=leaf_kws, ax=None, labels=(),
                       label_kws=label_kws, compute_layout=True):
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

    # Fresh tree instance
    tree = dendropy.Tree(tree)

    if compute_layout:
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
        xs = tuple(node._x for node in tree.preorder_node_iter())
        ys = tuple(node._y for node in tree.preorder_node_iter())
        tree._xlim = min(xs), max(xs)
        tree._ylim = min(ys), max(ys)

    ############
    # Draw edges
    ############
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

    #############
    # Draw leaves
    #############
    ax.scatter(
        tuple(node._x for node in tree.leaf_node_iter()),
        tuple(node._y for node in tree.leaf_node_iter()),
        **leaf_kws)


    ax.set_xlim(tree._xlim)
    ax.set_ylim(tree._ylim)
    ax.axis("off")
    ax.set_yticks([])
    ax.invert_yaxis()

    return tree, ax



    # def highlight(self, labels, **kws):
    #     """kws passed to plt.scatter"""
    #     s = kws.pop("s", 20)
    #     c = kws.pop("c", "red")
    #     zorder = kws.pop("zorder", 19)
    #     clip_on = kws.pop("clip_on", False)
    #     linewidth = kws.pop("linewidth", 0.5)
    #     edgecolor = kws.pop("edgecolor", "white")
    #     nodes = tree.find_nodes(lambda n: Tree.nodeHasTaxonInLabels(labels, n))
    #     if not nodes:
    #         raise ValueError("No nodes labels in tree.")
    #     x = [node._x for node in nodes]
    #     y = [node._y for node in nodes]
    #     plt.scatter(x, y, s=s, c=c, zorder=zorder, clip_on=clip_on,
    #                 linewidth=linewidth, edgecolor=edgecolor, **kws)


    ############
    # Label taxa
    ############
    for node in tree.find_nodes(lambda n: taxon_in_node_labels(labels, n)):
        plt.text(node._x, node._y, node.taxon.label, **label_kws)



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
