b7 = 145, 155, 156, 158, 159, 189, 193

"""
Map cluster -> motifs
"""
_cluster_motifs = {
    "HK68": (
        "STKGSQS",
    ),
    "EN72": (
        "SYKGSQS",
        "SYKGNQS",
        "SYKGSQN",
    ),
    "VI75": (
        "NYKGSKD",
    ),
    "TX77": (
        "NYKESKN",
    ),
    "BK79": (
        "NYEESKN",
    ),
    "SI87": (
        "NHEEYRN",
    ),
    "BE89": (
        "KHEDYRS",
        "KHEEYRS"
    ),
    "BE92": (
        "NHKEYSS",
    ),
    "WU95": (
        "KHKEYSS",
    ),
    "SY97": (
        "KHQKYSS",
    ),
    "FU02": (
        "KTHKYSS",
        "KTHKFNS",
        "KTHKYNS",
    ),
    "CA04": (
        "NTHKFNS",
        "STHKFNS",
    ),
    "WI05": (
        "NTHKFNF",
    ),
    "PE09": (
        "NTHNFKF",
        "STHNFKF",
    ),
    "SW13": (
        "STHNSKF",
    ),
    "HK14": (
        "STHNYKF",
    ),
}

"""
Map motif -> cluster
"""
_motif_to_cluster = {}
for cluster, motifs in _cluster_motifs.items():
    for motif in motifs:
        _motif_to_cluster[motif] = cluster

"""Cluster colours."""
cluster_colors = {
    "HK68": "#A208BD",  # Purple
    "EN72": "#33CCCC",  # Dark cyan
    "VI75": "#F9D004",  # Yellow
    "TX77": "#AB4C00",  # Brown
    "BK79": "#3BBA30",  # Green
    "SI87": "#0000FF",  # Blue
    "BE89": "#FF0000",  # Red
    "BE92": "#F894F8",  # Pink
    "WU95": "#37802B",  # Dark green
    "SY97": "#00AFFF",  # Light blue
    "FU02": "#B3C261",  # Green
    "CA04": "#FC5A03",  # Lemon
    "WI05": "#3E809C",  # Blue
    "PE09": "#F0FC03",  # Orange
    "SW13": "#FFD700",  # Gold
    "HK14": "#9CA9B5",  # Grey
}


def cluster_from_ha(sequence, seq_type="long"):
    """Classify an amino acid sequence as an antignenic cluster.

    This function only tries to classify the one or two dominant sequences in
    any cluster.

    Args:
        sequence (str): HA amino acid sequence.
        seq_type (str): "long" or "b7". If long, sequence must contain at least
            the fist 193 positions of HA1. If b7, sequence should be the b7
            positions.

    Raises:
        ValueError: If the sequence can't be classified.

    Returns:
        (str): The name of the cluster.
    """
    if seq_type == "long":
        sequence = "".join(sequence[i - 1] for i in b7)
    elif seq_type == "b7":
        if len(sequence) != 7:
            raise ValueError("sequence should be len 7 if seq_type = b7.")
    else:
        raise ValueError("seq_type should be 'long' or 'b7'")

    try:
        return _motif_to_cluster[sequence.upper()]
    except KeyError:
        raise ValueError("Can't classify {}".format(sequence))
