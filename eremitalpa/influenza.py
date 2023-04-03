from typing import Union, Iterable
from operator import itemgetter

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from .bio import find_mutations, hamming_dist
from .eremitalpa import get_trunk, plot_tree

"""
Influenza related data.
"""

b7 = 145, 155, 156, 158, 159, 189, 193
_clusters = (
    "HK68",
    "EN72",
    "VI75",
    "TX77",
    "BK79",
    "SI87",
    "BE89",
    "BE92",
    "WU95",
    "SY97",
    "FU02",
    "CA04",
    "WI05",
    "PE09",
    "SW13",
    "HK14",
    "KA17",
    "SW17",
    "HK19",
    "CA20",
    "DA21",
)

"""
Map cluster -> motifs
"""
_cluster_motifs = {
    "HK68": ("STKGSQS",),
    "EN72": (
        "SYKGSQS",
        "SYKGNQS",
        "SYKGSQN",
    ),
    "VI75": ("NYKGSKD",),
    "TX77": ("NYKESKN",),
    "BK79": (
        "NYEESKN",
        "NYEEYKN",
    ),
    "SI87": ("NHEEYRN",),
    "BE89": ("KHEDYRS", "KHEEYRS"),
    "BE92": ("NHKEYSS",),
    "WU95": ("KHKEYSS",),
    "SY97": ("KHQKYSS",),
    "FU02": (
        "KTHKYSS",
        "KTHKFNS",
        "KTHKYNS",
    ),
    "CA04": (
        "NTHKFNS",
        "STHKFNS",
    ),
    "WI05": ("NTHKFNF",),
    "PE09": (
        "NTHNFKF",
        "STHNFKF",
    ),
    "SW13": ("STHNSKF",),
    "HK14": ("STHNYKF",),
}

"""
Dict of dicts. For each cluster, map site -> amino acid, for all sites in
cluster transitions in and out.
"""
_cluster_key_residues = {
    "HK68": {155: "T"},
    "EN72": {155: "Y", 189: "Q"},
    "VI75": {158: "G", 189: "K", 193: "D"},
    "TX77": {156: "K", 158: "E", 193: "N"},
    "BK79": {155: "Y", 156: "E", 159: "S", 189: "K"},
    "SI87": {145: "N", 155: "H", 156: "E", 159: "Y", 189: "R"},
    "BE89": {145: "K"},
    "BE92": {156: "K", 145: "N"},
    "WU95": {145: "K", 156: "K", 158: "E"},
    "SY97": {156: "Q", 158: "K"},
    "FU02": {145: "K", 156: "H"},
    "CA04": {145: "N", 193: "S"},
    "WI05": {193: "F", 158: "K", 189: "N"},
    "PE09": {158: "N", 189: "K", 159: "F"},
    "SW13": {159: "S", 193: "F"},
    "HK14": {142: {"R", "G"}, 159: "Y", 193: "F"},
    "KA17": {193: "S"},
    "SW17": {142: "K"},
    "HK19": {135: "K", 193: "S"},
    "CA20": {159: "Y", 193: "S"},
    "DA21": {159: "N"},
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
    "HK68": "#A208BD",  # Purple
    "BE92": "#F894F8",  # Pink
    "BK79": "#3BBA30",  # Green
    "CA04": "#FC5A03",  # Lemon
    "EN72": "#33CCCC",  # Dark cyan
    "FU02": "#B3C261",  # Green
    "HK14": "#9CA9B5",  # Grey
    "PE09": "#F0FC03",  # Orange
    "SI87": "#0000FF",  # Blue
    "SW13": "#B3DE69",  # Green
    "SY97": "#00AFFF",  # Light blue
    "TX77": "#AB4C00",  # Brown
    "VI75": "#F9D004",  # Yellow
    "WI05": "#3E809C",  # Blue
    "WU95": "#37802B",  # Dark green
    "KA17": "#DB073D",  # Hot pink
    "SW17": "#4c1273",  # Dark purple
    "HK19": "#806205",  # Middle Brown
    "CA20": "#ab47bc",  # Pink
    "DA21": "#07485B",  # Deep turquoise
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


def cluster_from_ha_2(sequence: str, strict_len: bool = True, max_hd: float = 10):
    """
    Classify an amino acid sequence into an antigenic cluster.

    First identify clusters that have matching key residues with the sequence.
    If multiple clusters are found, find the one with the lowest hamming
    distance to the sequence. If the resulting hamming distance is less than
    10, return the cluster.

    Args:
        sequence (str)
        strict_len (bool): See hamming_to_cluster.
        hd (int): Queries that have matching key residues to a cluster are not
            classified as a cluster if the hamming distance to the cluster
            consensus is > hd.

    Returns:
        Cluster
    """
    candidates = clusters_with_matching_key_residues(sequence)

    if len(candidates) == 0:
        raise NoMatchingKeyResidues(sequence)

    elif len(candidates) > 1:
        cluster, hd = min(
            hamming_to_clusters(sequence, candidates, strict_len), key=itemgetter(1)
        )

    else:  # len(candidates) == 1
        cluster = candidates[0]
        hd = hamming_to_cluster(sequence, cluster, strict_len)

    if hd <= max_hd:
        return cluster

    else:
        raise HammingDistTooLargeError(
            f"{sequence}\nmatches key residues with {cluster} "
            f"but hamming distance is >{max_hd} ({hd})"
        )


def clusters_with_matching_key_residues(
    sequence: str, ignore: str = "-X"
) -> list["Cluster"]:
    """
    H3N2 clusters that have matching cluster transition substitution key
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
            if (
                (sequence[i] != aa)  # aa can be a single character
                and (sequence[i] not in ignore)
                and (sequence[i] not in aa)  # aa can be a set of possible amino acids
            ):
                break
        else:
            matches.append(Cluster(cluster))
    return matches


def hamming_to_all_clusters(sequence: str, strict_len: bool = True) -> list[float]:
    """The hamming distance from sequence to all known clusters.

    Args:
        sequence (str)
        strict_len (bool): See hamming_to_cluster

    Returns:
        2-tuples containing (cluster, hamming distance)
    """
    return [(c, hamming_to_cluster(sequence, c, strict_len)) for c in clusters]


def hamming_to_clusters(
    sequence: str, clusters: Iterable[Union[str, "Cluster"]], strict_len: bool = True
) -> list[float]:
    """The hamming distance from sequence to given clusters.

    Args:
        sequence (str)
        clusters (iterable)
        strict_len (bool): See hamming_to_cluster

    Returns:
        2-tuples containing (cluster, hamming distance)
    """
    return [(c, hamming_to_cluster(sequence, c, strict_len)) for c in clusters]


def hamming_to_cluster(
    sequence: str, cluster: Union[str, "Cluster"], strict_len: bool = True
) -> float:
    """The hamming distance from sequence to the consensus sequence of a
    cluster.

    Args:
        sequence (str)
        cluster (str or Cluster)
        strict_len (bool): Cluster consensus sequences are for HA1 only, and are
            328 residues long. If strict_len is True, then don't check whether
            sequence matches this length. If False, the sequence is truncated
            to 328 residues to match. If a sequence is less than 328 residues
            then an error will still be raised.

    Returns
        int
    """
    if not strict_len:
        sequence = sequence[:328]

    return hamming_dist(
        sequence,
        Cluster(cluster).sequence,
        ignore="-X",
        case_sensitive=False,
        per_site=False,
    )


class Cluster:
    def __init__(self, cluster):
        if str(cluster).upper() not in _clusters:
            raise ValueError(f"unknown cluster: {cluster}")

        self._name = str(cluster)

    def __repr__(self):
        return "Cluster('{}')".format(self._name)

    def __str__(self):
        return self._name.upper()

    def __gt__(self, other):
        if isinstance(other, Cluster):
            return self.year > other.year
        else:
            return str(self) > other

    def __lt__(self, other):
        if isinstance(other, Cluster):
            return self.year < other.year
        else:
            return str(self) < other

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

    @property
    def year(self) -> int:
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
    "HK68": (  # Cluster consensus
        "QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILD"
        "GINCTLIDALLGDPHCDVFQDETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEF"
        "ITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIW"
        "GVHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPG"
        "DVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKI"
        "TYGACPKYVKQNTLKLATGMRNVPEKQT"
    ),
    "EN72": (  # Cluster consensus
        "QDFPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRILD"
        "GIDCTLIDALLGDPHCDGFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEF"
        "INEGFTWTGVTQNGGSNACKRGPDSGFFSRLNWLYKSGSTYPVLNVTMPNNDNFDKLYIW"
        "GVHHPSTDQEQTSLYVQASGRVTVSTKRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPG"
        "DILVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIGTCISECITPNGSIPNDKPFQNVNKI"
        "TYGACPKYVKQNTLKLATGMRNVPEKQT"
    ),
    "VI75": (  # Cluster consensus
        "QDLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICDNPHRILD"
        "GINCTLIDALLGDPHCDGFQNEKWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEF"
        "INEGFNWTGVTQNGGSSACKRGPDNGFFSRLNWLYKSGSTYPVQNVTMPNNDNSDKLYIW"
        "GVHHPSTDKEQTDLYVQASGKVTVSTKRSQQTVIPNVGSRPWVRGLSSRVSIYWTIVKPG"
        "DILVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNKI"
        "TYGACPKYVKQNTLKLATGMRNVPEKQT"
    ),
    "TX77": (  # Cluster consensus
        "QDLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICNNPHRILD"
        "GINCTLIDALLGDPHCDGFQNKKWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEF"
        "INEGFNWTGVTQNGGSYACKRGPDNGFFSRLNWLYKSESTYPVLNVTMPNNDNFDKLYIW"
        "GVHHPSTDKEQTNLYVQASGRVTVSTKRSQQTIIPNVGSRPWVRGLSSRISIYWTIVKPG"
        "DILVINSNGNLIAPRGYFKIRNGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNKI"
        "TYGACPKYVKQNTLKLATGMRNVPEKQT"
    ),
    "BK79": (  # Cluster consensus
        "QNLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILD"
        "GKNCTLIDALLGDPHCDGFQNEKWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEF"
        "INEGFNWTGVTQSGGSYACKRGSDNSFFSRLNWLYESESKYPVLNVTMPNNGKFDKLYIW"
        "GVHHPSTDKEQTNLYVRASGRVTVSTKRSQQTVIPNIGSRPWVRGLSSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNKI"
        "TYGACPRYVKQNTLKLATGMRNVPEKQT"
    ),
    "SI87": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILD"
        "GKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "INEDFNWTGVAQSGGSYACKRGSVNSFFSRLNWLHESEYKYPALNVTMPNNGKFDKLYIW"
        "GVHHPSTDREQTNLYVRASGRVTVSTKRSQQTVIPNIGSRPWVRGLSSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQNTLKLATGMRNVPEKQT"
    ),
    "BE89": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILD"
        "GKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "TNEDFNWTGVAQSGESYACKRGSVKSFFSRLNWLHESDYKYPALNVTMPNNGKFDKLYIW"
        "GVHHPSTDREQTSLYVRASGRVTVSTKRSQQTVIPNIGSRPWVRGLSSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQNTLKLATGMRNVPEKQT"
    ),
    "BE92": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILD"
        "GKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "INEDFNWTGVAQDGKSYACKRGSVNSFFSRLNWLHKLEYKYPALNVTMPNNGKFDKLYIW"
        "GVHHPSTDSDQTSLYVRASGRVTVSTKRSQQTVIPNIGSRPWVRGLSSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRNGKSSIMRSDAPIGNCSSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQNTLKLATGMRNVPEKQT"
    ),
    "WU95": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILD"
        "GKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "TNEGFNWTGVAQDGKSYACKRGSVKSFFSRLNWLHKLEYKYPALNVTMPNNDKFDKLYIW"
        "GVHHPSTDSDQTSLYVQASGRVTVSTKRSQQTVIPNIGSRPWVRGISSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRNGKSSIMRSDAPIGNCNSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQNTLKLATGMRNVPEKQT"
    ),
    "SY97": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHQILD"
        "GENCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "NNESFNWTGVAQNGTSSACKRRSIKSFFSRLNWLHQLKYKYPALNVTMPNNEKFDKLYIW"
        "GVHHPSTDSDQISLYAQASGRVTVSTKRSQQTVIPNIGSRPWVRGVSSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDASIGKCNSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQNTLKLATGMRNVPEKQT"
    ),
    "FU02": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILD"
        "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "NNESFNWTGVTQNGTSSACKRRSNKSFFSRLNWLTHLKYKYPALNVTMPNNEKFDKLYIW"
        "GVHHPGTDSDQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRDVSSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQNTLKLATGMRNVPEKQT"
    ),
    "CA04": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILD"
        "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "NNESFNWTGVTQNGTSSACKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIW"
        "GVHHPGTDNDQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRDIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQNTLKLATGMRNVPEKQT"
    ),
    "WI05": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGEICDSPHQILD"
        "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "NNESFNWTGVTQNGTSSACIRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIW"
        "GVHHPGTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRVRNIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQNTLKLATGMRNVPEKQT"
    ),
    "PE09": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGEICDSPHQILD"
        "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "NNESFNWTGVTQNGTSSACIRRSNSSFFSRLNWLTHLNFKYPALNVTMPNNEQFDKLYIW"
        "GVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRNIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQSTLKLATGMRNVPEKQT"
    ),
    "SW13": (  # Cluster consensus
        "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "NNESFNWAGVTQNGTSSSCIRGSNSSFFSRLNWLTHLNSKYPALNVTMPNNEQFDKLYIW"
        "GVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQSTLKLATGMRNVPERQT"
    ),
    "HK14": (  # Cluster consensus
        "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "NNESFNWTGVTQNGTSSACIRRSSSSFFSRLNWLTHLNYTYPALNVTMPNNEQFDKLYIW"
        "GVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKHSTLKLATGMRNVPEKQT"
    ),
    "KA17": (  # Vaccine strain sequence https://www.ncbi.nlm.nih.gov/protein/AVG71503
        "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GENCTLIDALLGDPQCDGFQNKKWDLFVERNKAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "NNESFNWAGVTQNGTSSSCIRGSKSSFFSRLNWLTHLNSKYPALNVTMPNNEQFDKLYIW"
        "GVHHPGTDKDQISLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQSTLKLATGMRNVPERQT"
    ),
    "SW17": (  # A/SWITZERLAND/8060/2017 https://www.ncbi.nlm.nih.gov/protein/WCF71421
        "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSSCYPYDVPDYASLRSLVASSGTLEF"
        "NNESFNWTGVKQNGTSSACIRKSSSSFFSRLNWLTHLNYKYPALNVTMPNNEQFDKLYIW"
        "GVHHPGTDKDQIFPYAQSSGRIIVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIQSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKHSTLKLATGMRNVPEKQT"
    ),
    "HK19": (  # A/HONGKONG/2671/2019 gisaid EPI1698489
        "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GGNCTLIDALLGDPQCDGFQNKKWDLFVERSRAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "KNESFNWAGVTQNGKSFSCIRGSSSSFFSRLNWLTHLNYIYPALNVTMPNKEQFDKLYIW"
        "GVHHPVTDKDQISLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRNIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQSTLKLATGMRNVPEKQT"
    ),
    "CA20": (  # A/CALIFORNIA/55/2020 https://www.ncbi.nlm.nih.gov/protein/QQO51865
        "QKIPVNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GGNCTLIDALLGDPQCDGFQNKEWDLFVERSRANSNCYPYDVPDYASLRSLVASSGTLEF"
        "KNESFNWTGVKQNGTSSACIRGSSSSFFSRLNWLTHLNYTYPALNVTMPNNEQFDKLYIW"
        "GVHHPSTDKDQISLFAQPSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQSTLKLATGMRNVPEKQT"
    ),
    "DA21": (  # A/DARWIN/45/2021  gisaid EPI1928850
        "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GGNCTLIDALLGDPQCDGFQNKEWDLFVERSRANSNCYPYDVPDYASLRSLVASSGTLEF"
        "KNESFNWTGVKQNGTSSACIRGSSSSFFSRLNWLTHLNNIYPAQNVTMPNKEQFDKLYIW"
        "GVHHPDTDKNQISLFAQSSGRITVFTKRSQQTVIPNIGSRPRIRDIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQSTLKLATGMRNVPEKQT"
    ),
}

# Temporally sorted clusters
clusters = tuple(map(Cluster, _clusters))


class NoMatchingKeyResidues(Exception):
    pass


class HammingDistTooLargeError(Exception):
    pass


def plot_tree_coloured_by_cluster(
    tree,
    legend=True,
    leg_kws=dict(),
    unknown_color="black",
    leaf_kws=dict(),
    internal_kws=dict(),
    **kws,
):
    """Plot a tree with nodes coloured according to cluster.

    Args:
        tree (dendropy Tree): Nodes that have 'cluster' attribute
            will be coloured.
        legend (bool): Add a legend showing the clusters.
        leg_kws (dict): Keyword arguments passed to plt.legend.
        unknown_color (mpl color): Color if cluster is not known.
        **kws: Keyword arguments passed to plot_tree.
    """
    leaf_color = [
        cluster_colors[l.cluster]
        if hasattr(l, "cluster") and l.cluster in cluster_colors
        else unknown_color
        for l in tree.leaf_node_iter()
    ]
    internal_color = [
        cluster_colors[n.cluster]
        if hasattr(n, "cluster") and n.cluster in cluster_colors
        else unknown_color
        for n in tree.internal_nodes()
    ]

    leaf_kws = {
        **dict(color=leaf_color, zorder=10, s=5, linewidth=0.2, edgecolor="white"),
        **leaf_kws,
    }

    internal_kws = {
        **dict(color=internal_color, s=5, zorder=5, linewidth=0.2, edgecolor="white"),
        **internal_kws,
    }

    plot_tree(
        tree,
        leaf_kws=leaf_kws,
        internal_kws=internal_kws,
        ax=plt.gca(),
        compute_layout=True,
        **kws,
    )

    if legend:
        # Find clusters in this tree
        leaf_clusters = set(
            l.cluster if hasattr(l, "cluster") else None for l in tree.leaf_node_iter()
        )
        internal_clusters = set(
            n.cluster if hasattr(n, "cluster") else None for n in tree.internal_nodes()
        )
        all_clusters = leaf_clusters.union(internal_clusters)

        # Remove anything not a known cluster
        all_clusters = set(cluster_colors) & all_clusters
        all_clusters = sorted(Cluster(c) for c in all_clusters)

        handles = [Patch(facecolor=c.color, label=c) for c in all_clusters]
        plt.legend(handles=handles, **leg_kws)


def has_different_cluster_descendent(node):
    """Test if node has a descendent in a cluster different to its own.

    Args:
        node (dendropy Node)

    Returns:
        (bool)
    """
    descendents = list(node.postorder_internal_node_iter()) + node.leaf_nodes()
    for d in descendents:
        if d.cluster and d.cluster != node.cluster:
            return True
    return False


def guess_clusters_in_tree(node):
    """
    If a node is in a known cluster, and all of it's descendents are in the
    same cluster, or an unknown cluster, then, update all the descendent nodes
    to the matching cluster.
    """
    if hasattr(node, "seed_node"):
        # This is a tree
        guess_clusters_in_tree(node.seed_node)

    elif (node.cluster is None) or has_different_cluster_descendent(node):
        return

    else:
        descendents = list(node.postorder_internal_node_iter()) + node.leaf_nodes()
        for d in descendents:
            if d.cluster is None:
                d.cluster = node.cluster


def plot_subs_on_tree(
    tree,
    seq_attr,
    cluster_change_only=None,
    length=30,
    exclude_leaves=True,
    find_mutation_offset=0,
    max_mutations=20,
    only_these_positions=None,
    exclude_characters="X",
    either_side_trunk=True,
    trunk_attr="_x",
    **kws,
):
    """Annotate a tree with substitutions.

    Args:
        tree (dendropy Tree)
        seq_attr (str): Name of the attribute on nodes that contain the
            sequence.
        cluster_change_only (bool): Only plot substitutions on nodes when a
            cluster has changed.
        length (scalar): Length of the line.
        exclude_leaves (bool): Don't label substitutions on branches leading to
            leaf nodes.
        find_mutation_offset (int): See ere.find_mutations.
        max_mutations (int): Annotate at most this number of mutations.
        exclude_characters (str): If a mutation contains a character in this
            string, don't annotate it.
        only_these_positions (iterable): Contains ints. Only show mutations that
            at these positions.
        either_side_trunk (bool): Plot labels both sides of the trunk.
        trunk_attr (str): _x or _y. Trunk is defined as root to deepest leaf.
            Deepest leaf is the leaf with maximum trunk_attr.

        **kws: Keyword arguments passed to plt.annotate.
    """
    trunk = get_trunk(tree, trunk_attr)

    length = (length**2 / 2) ** 0.5

    for node in tree.internal_nodes():
        if exclude_leaves and node.is_leaf():
            continue

        if node.parent_node:
            parent = node.parent_node

            # try:
            #     if cluster_change_only and node.cluster == parent.cluster:
            #         continue
            # except AttributeError:
            #     pass

            a = getattr(parent, seq_attr)
            b = getattr(node, seq_attr)

            mutations = find_mutations(a, b, offset=find_mutation_offset)

            # Apply filters to mutations
            if only_these_positions:
                mutations = filter(lambda m: m.pos in only_these_positions, mutations)

            if exclude_characters:

                def has_filtered(mutation):
                    for aa in mutation.a, mutation.b:
                        if aa in exclude_characters:
                            return True
                    return False

                mutations = filter(lambda m: not has_filtered(m), mutations)

            mutations = sorted(mutations)

            if len(mutations) == 0:
                continue

            elif len(mutations) > max_mutations:
                mutations = mutations[:max_mutations]
                mutations += "+"

            if either_side_trunk and node in trunk:
                xytext = -1 * length, -1 * length
                va = "top"
                ha = "right"
            else:
                xytext = length, length
                va = "bottom"
                ha = "left"

            label = "\n".join(map(str, mutations))
            xy = (node._x + parent._x) / 2, node._y

            plt.annotate(
                label,
                xy,
                xytext=xytext,
                va=va,
                ha=ha,
                textcoords="offset pixels",
                arrowprops=dict(
                    facecolor="darkgrey",
                    shrink=0,
                    linewidth=0,
                    width=0.3,
                    headwidth=2,
                    headlength=2,
                ),
                **kws,
            )
