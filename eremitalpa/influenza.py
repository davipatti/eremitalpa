from functools import partial
from operator import itemgetter

"""
Influenza related data.
"""

b7 = 145, 155, 156, 158, 159, 189, 193
clusters = ("BE89", "BE92", "BK79", "CA04", "EN72", "FU02", "HK14", "HK68",
            "KA17", "PE09", "SI87", "SW13", "SY97", "TX77", "VI75", "WI05", "WU95")

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
        "NYEEYKN",
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
Dict of dicts. For each cluster, map site -> amino acid, for all sites in
cluster transitions in and out.
"""
_cluster_key_residues = {
    "HK68": {
        155: "T"
    },
    "EN72": {
        155: "Y",
        189: "Q",
    },
    "VI75": {
        158: "G",
        189: "K",
        193: "D",
    },
    "TX77": {
        156: "K",
        158: "E",
        193: "N",
    },
    "BK79": {
        155: "Y",
        156: "E",
        159: "S",
        189: "K",
    },
    "SI87": {
        145: "N",
        155: "H",
        156: "E",
        159: "Y",
        189: "R",
    },
    "BE89": {
        145: "K",
    },
    "BE92": {
        156: "K",
        145: "N",
    },
    "WU95": {
        145: "K",
        156: "K",
        158: "E",
    },
    "SY97": {
        156: "Q",
        158: "K",
    },
    "FU02": {
        145: "K",
        156: "H",
    },
    "CA04": {
        145: "N",
        193: "S",
    },
    "WI05": {
        193: "F",
        158: "K",
        189: "N",
    },
    "PE09": {
        158: "N",
        189: "K",
        159: "F"
    },
    "SW13": {
        159: "S"
    },
    "HK14": {
        159: "Y"
    },
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
    "BE89": "#FF0000",  # Red
    "BE92": "#F894F8",  # Pink
    "BK79": "#3BBA30",  # Green
    "CA04": "#FC5A03",  # Lemon
    "EN72": "#33CCCC",  # Dark cyan
    "FU02": "#B3C261",  # Green
    "HK14": "#9CA9B5",  # Grey
    "HK68": "#A208BD",  # Purple
    "PE09": "#F0FC03",  # Orange
    "SI87": "#0000FF",  # Blue
    "SW13": "#B3DE69",  # Green
    "SY97": "#00AFFF",  # Light blue
    "TX77": "#AB4C00",  # Brown
    "VI75": "#F9D004",  # Yellow
    "WI05": "#3E809C",  # Blue
    "WU95": "#37802B",  # Dark green
}


def cluster_from_ha(sequence, seq_type="long"):
    """Classify an amino acid sequence as an antignenic cluster by checking
    whether the sequences Bjorn 7 sites match exactly sites that are known in
    a cluster.

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
        return Cluster(_motif_to_cluster[sequence.upper()])
    except KeyError:
        raise ValueError("Can't classify {}".format(sequence))


def cluster_from_ha_2(sequence):
    """Classify an amino acid sequence into an antigenic cluster.

    First, test if the B7 motif has been seen before (i.e. run
    cluster_from_ha). If the motif has been seen before, then return the
    corresponding cluster. If it hasn't been seen, then look for clusters
    with matching key residues. If there are no clusters that match all key
    residues, then give up. Else, check which cluster has the consensus sequence
    with the lowest hamming distance to the sequence. If the cluster with the
    lowest hammning distance also has matching key resiudes, return this
    cluster.

    Args:
        sequence (str)

    Returns:
        Cluster
    """
    try:
        return cluster_from_ha(sequence, seq_type="long")
    except ValueError:
        candidate_clusters = clusters_with_matching_key_residues(sequence)

        if not candidate_clusters:
            # Sequence doesn't have key residues that match any single cluster
            raise ValueError("Can't classify {}".format(sequence))

        # Find cluster with the smallest hamming distance to sequence
        nearest_seq_cluster, hd = min(
            hamming_to_all_clusters(sequence), key=itemgetter(1))

        if nearest_seq_cluster in candidate_clusters:
            return Cluster(nearest_seq_cluster)

        else:
            raise ValueError(
                "Cluster with lowest hamming distance is {} ({}), but key "
                "residues differ. Sequence:\n{}".format(
                    nearest_seq_cluster, hd, sequence))


def clusters_with_matching_key_residues(sequence, ignore="-X"):
    """H3N2 clusters that have matching cluster transition substitution key
    residues.

    Args:
        sequence (str): Amino acid sequence. At least 193 residues long.
        ignore (str): Ignore these characters.

    Returns:
        list of clusters with matching key residues.
    """
    sequence = sequence.upper()
    matches = []
    for cluster, site_aa in _cluster_key_residues.items():
        for site, aa in site_aa.items():
            i = site - 1
            if (sequence[i] != aa) and (sequence[i] not in ignore):
                break
        else:
            matches.append(Cluster(cluster))
    return matches


def hamming_to_all_clusters(sequence):
    """The hamming distance from sequence to all known clusters.

    Args:
        sequence (str)

    Returns:
        2-tuples containing (cluster, hamming distance)
    """
    return ((c, hamming_to_cluster(sequence, c)) for c in clusters)


def hamming_to_cluster(sequence, cluster, **kws):
    """The hamming distance from sequence to the sequence of a cluster.

    Args:
        sequence (str)
        cluster (str or Cluster)
        **kws passed to hamming_dist

    Returns
        int
    """
    return hamming_dist(sequence, Cluster(cluster).sequence, **kws)


def hamming_dist(a, b, ignore="-X", per_site=False):
    """The hamming distance between a and b.

    Args:
        a (str)
        b (str)
        ignore (str): String containing characters to ignore. If there is a
            mismatch where one string has a character in igonre, this does not
            contribute to the hamming distance.
        per_site (bool): Divide the hamming distance by the length of a and b,
            minus the number of sites with ignored characters.


    Returns:
        int
    """
    if len(a) != len(b):
        raise ValueError("a and b must be same len")
    a = a.upper()
    b = b.upper()
    d = 0
    if per_site:
        l = 0
        for m, n in zip(a, b):
            if (m not in ignore) and (n not in ignore):
                l += 1
                if m != n:
                    d += 1
        d /= l
    else:
        for m, n in zip(a, b):
            if (m != n) and (m not in ignore) and (n not in ignore):
                d += 1
    return d


class Cluster():

    def __init__(self, cluster):
        self._name = str(cluster)

    def __repr__(self):
        return "Cluster('{}')".format(self._name)

    def __str__(self):
        return self._name.upper()

    def __gt__(self, other):
        return self.year > other.year

    def __lt__(self, other):
        return self.year < other.year

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

    @property
    def year(self):
        digits = int(self._name[-2:])
        if digits < 67:
            return digits + 2000
        else:
            return digits + 1900

    @property
    def key_residues(self):
        return _cluster_key_residues[self._name]

    @property
    def b7_motifs(self):
        return _cluster_motifs[self._name]

    @property
    def color(self):
        return cluster_colors[self._name]

    @property
    def sequence(self):
        return _cluster_sequences[self._name]


_cluster_sequences = {
    "BE89": "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILDGKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFTNEDFNWTGVAQSGESYACKRGSVKSFFSRLNWLHESDYKYPALNVTMPNNGKFDKLYIWGVHHPSTDREQTSLYVRASGRVTVSTKRSQQTVIPNIGSRPWVRGLSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT",
    "BE92": "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILDGKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFINEDFNWTGVAQDGKSYACKRGSVNSFFSRLNWLHKLEYKYPALNVTMPNNGKFDKLYIWGVHHPSTDSDQTSLYVRASGRVTVSTKRSQQTVIPNIGSRPWVRGLSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRNGKSSIMRSDAPIGNCSSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT",
    "BK79": "QNLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILDGKNCTLIDALLGDPHCDGFQNEKWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFINEGFNWTGVTQSGGSYACKRGSDNSFFSRLNWLYESESKYPVLNVTMPNNGKFDKLYIWGVHHPSTDKEQTNLYVRASGRVTVSTKRSQQTVIPNIGSRPWVRGLSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNKITYGACPRYVKQNTLKLATGMRNVPEKQT",
    "CA04": "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT",
    "EN72": "QDFPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDGFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFINEGFTWTGVTQNGGSNACKRGPDSGFFSRLNWLYKSGSTYPVLNVTMPNNDNFDKLYIWGVHHPSTDQEQTSLYVQASGRVTVSTKRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDILVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIGTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT",
    "FU02": "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRRSNKSFFSRLNWLTHLKYKYPALNVTMPNNEKFDKLYIWGVHHPGTDSDQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRDVSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT",
    "HK14": "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSSSSFFSRLNWLTHLNYTYPALNVTMPNNEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKHSTLKLATGMRNVPEKQT",
    "HK68": "QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGINCTLIDALLGDPHCDVFQDETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGVHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT",
    "KA17": "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERNKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWAGVTQNGTSSSCIRGSKSSFFSRLNWLTHLNSKYPALNVTMPNNEQFDKLYIWGVHHPGTDKDQISLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPERQT",
    "PE09": "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSNSSFFSRLNWLTHLNFKYPALNVTMPNNEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPEKQT",
    "SI87": "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILDGKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFINEDFNWTGVAQSGGSYACKRGSVNSFFSRLNWLHESEYKYPALNVTMPNNGKFDKLYIWGVHHPSTDREQTNLYVRASGRVTVSTKRSQQTVIPNIGSRPWVRGLSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT",
    "SW13": "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWAGVTQNGTSSSCIRGSNSSFFSRLNWLTHLNSKYPALNVTMPNNEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPERQT",
    "SY97": "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHQILDGENCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVAQNGTSSACKRRSIKSFFSRLNWLHQLKYKYPALNVTMPNNEKFDKLYIWGVHHPSTDSDQISLYAQASGRVTVSTKRSQQTVIPNIGSRPWVRGVSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDASIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT",
    "TX77": "QDLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICNNPHRILDGINCTLIDALLGDPHCDGFQNKKWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFINEGFNWTGVTQNGGSYACKRGPDNGFFSRLNWLYKSESTYPVLNVTMPNNDNFDKLYIWGVHHPSTDKEQTNLYVQASGRVTVSTKRSQQTIIPNVGSRPWVRGLSSRISIYWTIVKPGDILVINSNGNLIAPRGYFKIRNGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT",
    "VI75": "QDLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICDNPHRILDGINCTLIDALLGDPHCDGFQNEKWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFINEGFNWTGVTQNGGSSACKRGPDNGFFSRLNWLYKSGSTYPVQNVTMPNNDNSDKLYIWGVHHPSTDKEQTDLYVQASGKVTVSTKRSQQTVIPNVGSRPWVRGLSSRVSIYWTIVKPGDILVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT",
    "WI05": "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRVRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT",
    "WU95": "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILDGKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFTNEGFNWTGVAQDGKSYACKRGSVKSFFSRLNWLHKLEYKYPALNVTMPNNDKFDKLYIWGVHHPSTDSDQTSLYVQASGRVTVSTKRSQQTVIPNIGSRPWVRGISSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRNGKSSIMRSDAPIGNCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT",
}
