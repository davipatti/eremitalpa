import eremitalpa

from .bio import (
    amino_acid_colors,
    AMINO_ACIDS,
    consensus_seq,
    filter_similar_hd,
    find_mutations,
    find_substitutions,
    grouped_sample,
    hamming_dist_lt,
    hamming_dist,
    pairwise_hamming_dists,
    sloppy_translate,
    variable_sites,
    Substitution,
    TiedCounter,
)
from .eremitalpa import (
    load_fasta,
    compare_trees,
    compute_tree_layout,
    color_stack,
    deepest_leaf,
    get_trunk,
    plot_leaves_with_labels,
    plot_tree,
    prune_nodes_with_labels,
    read_raxml_ancestral_sequences,
    taxon_in_node_label,
    taxon_in_node_labels,
    MultipleSequenceAlignment,
    Tree,
)
from .influenza import (
    aa_counts_thru_time,
    _cluster_motifs,
    cluster_colors,
    cluster_from_ha_2,
    cluster_from_ha,
    cluster_transitions,
    Cluster,
    clusters,
    ClusterTransition,
    guess_clusters_in_tree,
    hamming_to_all_clusters,
    hamming_to_cluster,
    plot_subs_on_tree,
    plot_aa_freq_thru_time,
    plot_tree_coloured_by_cluster,
    translate_trim_default_ha,
    NHSeason,
)
from .flu_wider import b7, wider5, wider10, wider15
from .lib import (
    annotate_points,
    compute_errorbars,
    find_runs,
    jupyter_code_cell_toggle,
    log_df_func,
    split_pairs,
    cal_months_diff,
)
from .spread_points import spread_points

__version__ = "1.1.4"

__all__ = [
    "__version__",
    "_cluster_motifs",
    "aa_counts_thru_time",
    "amino_acid_colors",
    "AMINO_ACIDS",
    "annotate_points",
    "b7",
    "cluster_colors",
    "cluster_from_ha_2",
    "cluster_from_ha",
    "cluster_transitions",
    "Cluster",
    "clusters",
    "ClusterTransition",
    "compare_trees",
    "compute_errorbars",
    "compute_tree_layout",
    "consensus_seq",
    "color_stack",
    "eremitalpa",
    "deepest_leaf",
    "filter_similar_hd",
    "find_mutations",
    "find_runs",
    "get_trunk",
    "grouped_sample",
    "guess_clusters_in_tree",
    "hamming_dist_lt",
    "hamming_dist",
    "hamming_to_all_clusters",
    "NHSeason",
    "hamming_to_cluster",
    "jupyter_code_cell_toggle",
    "log_df_func",
    "pairwise_hamming_dists",
    "plot_leaves_with_labels",
    "plot_subs_on_tree",
    "plot_tree_coloured_by_cluster",
    "plot_tree",
    "prune_nodes_with_labels",
    "read_raxml_ancestral_sequences",
    "sloppy_translate",
    "Substitution",
    "spread_points",
    "taxon_in_node_label",
    "taxon_in_node_labels",
    "translate_trim_default_ha",
    "TiedCounter",
    "wider10",
    "wider15",
    "wider5",
    "split_pairs",
    "cal_months_diff",
    "plot_aa_freq_thru_time",
    "find_substitutions",
    "MultipleSequenceAlignment",
    "Tree",
    "variable_sites",
    "load_fasta",
]
