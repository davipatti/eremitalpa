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
