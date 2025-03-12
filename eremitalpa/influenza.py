from functools import partial
from pathlib import Path
from typing import Union, Iterable, Generator
import json

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

from .lib import cal_months_diff
from .bio import find_mutations, hamming_dist, sloppy_translate, amino_acid_colors
from .eremitalpa import get_trunk, plot_tree

"""
Influenza related data.
"""

b7 = 145, 155, 156, 158, 159, 189, 193

# See also `clusters` defined below
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

# See also `cluster_transitions` defined below
_cluster_transitions = (
    ("HK68", "EN72"),
    ("EN72", "VI75"),
    ("VI75", "TX77"),
    ("TX77", "BK79"),
    ("BK79", "SI87"),
    ("SI87", "BE89"),
    ("SI87", "BE92"),
    ("BE92", "WU95"),
    ("WU95", "SY97"),
    ("SY97", "FU02"),
    ("FU02", "CA04"),
    ("CA04", "WI05"),
    ("WI05", "PE09"),
    ("PE09", "SW13"),
    ("PE09", "HK14"),
    ("SW13", "KA17"),
    ("HK14", "SW17"),
    ("HK14", "HK19"),
    ("HK14", "CA20"),
    ("CA20", "DA21"),
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
    """Classify an amino acid sequence as an antigenic cluster by checking
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
        cluster_hd = dict(
            hamming_to_clusters(sequence, candidates, strict_len=strict_len)
        )

        # Filter candidates by hamming distance, the lower the better.
        # If there are multiple candidate clusters that have the same lowest hamming
        # distance, then raise an error
        lowest_hd = min(cluster_hd.values())

        # All clusters that have the lowest HD
        lowest_hd_clusters = [
            cluster for cluster, hd in cluster_hd.items() if hd == lowest_hd
        ]

        if len(lowest_hd_clusters) > 1:
            raise TiedHammingDistances(lowest_hd_clusters)

        else:
            cluster = lowest_hd_clusters[0]
            hd = lowest_hd

    else:  # len(candidates) == 1
        cluster = candidates[0]
        hd = hamming_to_cluster(sequence, cluster, strict_len=strict_len)

    if hd <= max_hd:
        return cluster

    else:
        raise HammingDistTooLargeError(
            f"{sequence}\nmatches key residues with {cluster} "
            f"but hamming distance is >{max_hd} ({hd})"
        )


def clusters_with_matching_key_residues(sequence: str) -> list["Cluster"]:
    """
    List of H3N2 clusters that have matching key residues.

    Args:
        sequence (str): Amino acid sequence. At least 193 residues long (highest numeric
            position of  a key residue).
    """
    sequence = sequence.upper()
    matches = []
    for cluster, site_aa in _cluster_key_residues.items():
        for site, aa in site_aa.items():
            i = site - 1

            # aa can be a single char str or a set of single char str
            if sequence[i] not in aa:
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
        Cluster(cluster).aa_sequence,
        ignore="-X",
        case_sensitive=False,
        per_site=False,
    )


def load_cluster_nt_consensus() -> dict[str, str]:
    "Load cluster nt consensus seqs."
    path = Path(__file__).parent.parent.joinpath("data", "flu", "cluster_cons.json")
    with open(path, "r") as fp:
        return json.load(fp)


class Cluster:

    def __init__(self, cluster) -> None:
        if str(cluster).upper() not in _clusters:
            raise ValueError(f"unknown cluster: {cluster}")

        self._name = str(cluster)

    def __repr__(self) -> str:
        return "Cluster('{}')".format(self._name)

    def __str__(self) -> str:
        return self._name.upper()

    def __gt__(self, other: Union[str, "Cluster"]) -> bool:
        return self.year > Cluster(other).year

    def __lt__(self, other: Union[str, "Cluster"]) -> bool:
        return self.year < Cluster(other).year

    def __eq__(self, other: Union[str, "Cluster"]) -> bool:
        return str(self) == str(other)

    def __hash__(self) -> int:
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
    def aa_sequence(self: str):
        """Representative amino acid sequence."""
        return _cluster_aa_sequences[self._name]

    @property
    def nt_sequence(self) -> str:
        """Representative nucleotide sequence."""
        if not hasattr(self, "_cluster_nt_sequences"):
            self._cluster_nt_sequences = load_cluster_nt_consensus()
        return self._cluster_nt_sequences[self._name]

    def codon(self, n: int) -> str:
        """Codon at amino acid position n. 1-indexed."""
        return self.nt_sequence[(n - 1) * 3 : n * 3]


class ClusterTransition:
    def __init__(self, c0: Union[str, Cluster], c1: Union[str, Cluster]) -> None:
        """A cluster transition."""
        if (c0, c1) not in _cluster_transitions:
            raise ValueError(f"unrecognised cluster transition: {c0, c1}")

        self.c0 = Cluster(c0)
        self.c1 = Cluster(c1)
        self._clusters = self.c0, self.c1

        if self.c0.year >= self.c1.year:
            raise ValueError("c0 year must be less than c1 year")

    @classmethod
    def from_tuple(cls, c0c1) -> "ClusterTransition":
        """Make an instance from a tuple"""
        return cls(*c0c1)

    def __repr__(self) -> str:
        return f"ClusterTransition({self.c0}, {self.c1})"

    def __str__(self) -> str:
        return f"{self.c0} -> {self.c1}"

    def __eq__(self, other) -> bool:
        if isinstance(other, tuple):
            other = ClusterTransition.from_tuple(other)
        return self.c0 == other.c0 and self.c1 == other.c1

    def __lt__(self, other) -> bool:
        return self.c1 < other.c1

    def __gt__(self, other) -> bool:
        return self.c1 > other.c1

    def __getitem__(self, item) -> Cluster:
        return self._clusters[item]

    def __hash__(self) -> int:
        return hash(self._clusters)

    @property
    def preceding_transitions(self) -> Generator["ClusterTransition", None, None]:
        """All preceding cluster transitions"""
        for ct in cluster_transitions:
            if ct < self:
                yield ct


_cluster_aa_sequences = {
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
    "HK19": (  # A/HONGKONG/2671/2019 GISAID EPI1698489
        "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GGNCTLIDALLGDPQCDGFQNKKWDLFVERSRAYSNCYPYDVPDYASLRSLVASSGTLEF"
        "KNESFNWAGVTQNGKSFSCIRGSSSSFFSRLNWLTHLNYIYPALNVTMPNKEQFDKLYIW"
        "GVHHPVTDKDQISLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRNIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQSTLKLATGMRNVPEKQT"
    ),
    "CA20": (  # A/Cambodia/e0826360/2020 SIAT GISAID EPI1843589
        "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GGNCTLIDALLGDPQCDGFQNKEWDLFVERSRANSNCYPYDVPDYASLRSLVASSGTLEF"
        "KNESFNWTGVKQNGTSSACIRGSSSSFFSRLNWLTHLNYTYPALNVTMPNNEQFDKLYIW"
        "GVHHPSTDKDQISLFAQPSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQSTLKLATGMRNVPEKQT"
    ),
    "DA21": (  # A/DARWIN/45/2021  GISAID EPI1928850
        "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
        "GGNCTLIDALLGDPQCDGFQNKEWDLFVERSRANSNCYPYDVPDYASLRSLVASSGTLEF"
        "KNESFNWTGVKQNGTSSACIRGSSSSFFSRLNWLTHLNNIYPAQNVTMPNKEQFDKLYIW"
        "GVHHPDTDKNQISLFAQSSGRITVFTKRSQQTVIPNIGSRPRIRDIPSRISIYWTIVKPG"
        "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
        "TYGACPRYVKQSTLKLATGMRNVPEKQT"
    ),
}

# Temporally sorted clusters and cluster transitions
clusters = tuple(Cluster(c) for c in _clusters)
cluster_transitions = tuple(ClusterTransition(*pair) for pair in _cluster_transitions)


class NoMatchingKeyResidues(Exception):
    pass


class HammingDistTooLargeError(Exception):
    pass


class TiedHammingDistances(Exception):
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
        (
            cluster_colors[leaf.cluster]
            if hasattr(leaf, "cluster") and leaf.cluster in cluster_colors
            else unknown_color
        )
        for leaf in tree.leaf_node_iter()
    ]
    internal_color = [
        (
            cluster_colors[node.cluster]
            if hasattr(node, "cluster") and node.cluster in cluster_colors
            else unknown_color
        )
        for node in tree.internal_nodes()
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
            leaf.cluster if hasattr(leaf, "cluster") else None
            for leaf in tree.leaf_node_iter()
        )
        internal_clusters = set(
            node.cluster if hasattr(node, "cluster") else None
            for node in tree.internal_nodes()
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
    descendants = list(node.postorder_internal_node_iter()) + node.leaf_nodes()
    for d in descendants:
        if d.cluster and d.cluster != node.cluster:
            return True
    return False


def guess_clusters_in_tree(node):
    """
    If a node is in a known cluster, and all of it's descendants are in the
    same cluster, or an unknown cluster, then, update all the descendent nodes
    to the matching cluster.
    """
    if hasattr(node, "seed_node"):
        # This is a tree
        guess_clusters_in_tree(node.seed_node)

    elif (node.cluster is None) or has_different_cluster_descendent(node):
        return

    else:
        descendants = list(node.postorder_internal_node_iter()) + node.leaf_nodes()
        for d in descendants:
            if d.cluster is None:
                d.cluster = node.cluster


def plot_subs_on_tree(
    tree,
    seq_attr,
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


def translate_trim_default_ha(nt: str) -> str:
    """
    Take a default HA nucleotide sequence and return an HA1 sequence.
    """
    return sloppy_translate(nt)[16 : 328 + 16]


def aa_counts_thru_time(df_seq: pd.DataFrame, site: int, ignore="-X") -> pd.DataFrame:
    """
    Make a DataFrame containing counts of amino acids that were sampled in months at a
    particular sequence site.

    Columns in the returned DataFrame are dates, the index is amino acids.

    Args:
        df_seq: Must contain columns "aa" which contains amino acid sequences and "dt"
            which contains datetime objects of collection dates.
        site: Count amino acids at this site. Note, this is 1-indexed.
        ignore: Don't include these characters in the counts.
    """
    df = (
        pd.DataFrame(
            {
                month: df_grp["aa"].str[site - 1].value_counts().sort_index()
                for (month, df_grp) in df_seq.groupby(pd.Grouper(key="dt", freq="M"))
            }
        )
        .T.resample("M")  # Make sure we're not missing any months
        .asfreq()
        .fillna(0)
        .T
    )
    aas = sorted(set(df.index) - set(ignore))
    return df.loc[aas]


def plot_aa_freq_thru_time(
    t0: pd.Timestamp,
    t_end: pd.Timestamp,
    df_seq: pd.DataFrame,
    site: int,
    proportion=False,
    ax=None,
    ignore="X-",
    blank_xtick_labels=False,
):
    ax = plt.gca() if ax is None else ax

    df = aa_counts_thru_time(df_seq.query("dt > @t0"), site=site, ignore=ignore)

    if proportion:
        df /= df.sum(axis=0)

    x = [cal_months_diff(month, t0) for month in df.columns]

    bottom = 0
    for aa, row in df.iterrows():
        ax.bar(
            x=x,
            height=row,
            bottom=bottom,
            facecolor=amino_acid_colors[aa],
            width=1,
            align="edge",
            label=aa,
        )
        bottom += row

    xticks = df.columns[::12]

    ax.set_xticks(
        [cal_months_diff(x, t0) for x in xticks],
        labels=(
            ["" for _ in xticks]
            if blank_xtick_labels
            else xticks.astype(str).str.slice(0, 7)
        ),
    )
    ax.set_xlim(0, cal_months_diff(t_end, t0))


class NHSeason:
    def __init__(self, years: Union[tuple[int, int], "NHSeason"]) -> None:
        """
        A northern hemisphere flu season.
        """
        for y in years:
            if not isinstance(y, int):
                raise ValueError("years should be integers")

        self.years = years
        self.y0, self.y1 = years

        if self.y1 != self.y0 + 1:
            raise ValueError("years should be consecutive increasing integers")

    def __repr__(self) -> str:
        return f"NHSeason(({self.y0}, {self.y1}))"

    def __str__(self) -> str:
        a = str(self.y0)
        b = str(self.y1)
        return f"{a[-2:]}-{b[-2:]}"

    def __gt__(self, other: "NHSeason") -> bool:
        return self.y0 > other.y0

    def __lt__(self, other: "NHSeason") -> bool:
        return self.y0 < other.y0

    def __eq__(self, other: "NHSeason") -> bool:
        return self.y0 == other.y0

    def __le__(self, other: "NHSeason") -> bool:
        return self.y0 <= other.y0

    def __ge__(self, other: "NHSeason") -> bool:
        return self.y0 >= other.y0

    def __hash__(self) -> int:
        return hash((self.y0, self.y1))

    @classmethod
    def from_datetime(cls, dt):
        if dt.month > 8:
            return cls((dt.year, dt.year + 1))
        else:
            return cls((dt.year - 1, dt.year))

    def __len__(self) -> int:
        return 2

    def __getitem__(self, item) -> int:
        return self.years[item]


def splice(
    sequence: str, start_ends: tuple[tuple[int, int], ...], translate: bool = True
) -> str:
    """Splice a sequence.

    Args:
        start_ends: Contains 2-tuples that define the start and end of coding sequences. Values are
            1-indexed, and inclusive. E.g. passing (2, 4) for the sequence 'ACTGT' would return
            'CTG'.
    """
    spliced = "".join(sequence[start - 1 : end] for (start, end) in start_ends)
    return sloppy_translate(spliced) if translate else spliced


def translate_segment(sequence: str, segment: str) -> dict[str, str]:
    """
    Transcribe an influenza A segment. MP, PA, PB1 and NS all have splice variants, so simply
    transcribing the ORF of the segment would miss out proteins. This function returns a list
    containing coding sequences that are transcribed from a particular segment.

    Args:
        sequence (str): The RNA sequence of the segment.
        segment (str): The segment to translate. Must be one of 'HA', 'NA', 'NP', 'PB2', 'PA', 'MP' or 'PB1'.

    Returns:
        dict[str, str]: A dictionary mapping the segment name to the translated protein sequence.
    """
    splice_translate = partial(splice, sequence=sequence, translate=True)

    if segment == "PA":
        if len(sequence) != 2151:
            raise ValueError("expected PA to be 2151 nts")

        return {
            "PA": splice_translate(start_ends=((1, None),)),
            "PA-X": splice_translate(start_ends=((1, 570), (572, 760))),
        }

    elif segment == "MP":
        if len(sequence) != 982:
            raise ValueError("expected MP to be 982 nts")

        return {
            "M1": splice_translate(start_ends=((1, 819),)),
            "M2": splice_translate(start_ends=((1, 26), (715, None))),
        }

    elif segment == "PB1":
        if len(sequence) != 2274:
            raise ValueError("expected PB1 to be 2274 nts")

        return {
            "PB1": splice_translate(start_ends=((1, None),)),
            "PB1-F2": splice_translate(start_ends=((95, 367),)),
        }

    elif segment == "NS":
        return {
            "NS1": splice_translate(start_ends=ns1_splice_sites(sequence)),
            "NS2": splice_translate(start_ends=ns2_splice_sites(sequence)),
        }

    else:
        # For these remaining segments, just translate the full sequence.
        if segment not in {"HA", "NA", "NP", "PB2"}:
            raise ValueError(
                "segment must be one of 'HA', 'NA', 'NP', 'PB2', 'PB2', 'PA', 'MP' or 'PB1'."
            )

        return {segment: sloppy_translate(sequence)}


def ns1_splice_sites(sequence) -> tuple[tuple[int, int]]:
    """
    Return the start and end positions of the NS1 gene in a given sequence.

    Given the sequence, find the splice donor and acceptor sites, and then return
    the start and end positions of the NS1 gene in terms of 1-based indexing.

    Returns:
        tuple[tuple[int, int]]: A single tuple containing the start and end positions
            of the NS1 gene.
    """
    _, accept_loc = find_ns_splice_sites(sequence)

    # NS1 has 193 nts after splice acceptor (AG)
    ns1_end = accept_loc + 193

    # Testing against A/Texas/37/2024, the NS1 end location was adding a single additional NT.
    # So, correct that here. I'm not sure if this is a difference between the sequences that the
    # flu-ngs pipeline is typically run against (where the bulk of this NS splicing code was written
    # for).
    ns1_end -= 1

    # +1 for 1-based indexing
    return ((1, ns1_end + 1),)


def ns2_splice_sites(sequence):
    """
    Return start and end positions for the two exons of the NS2 gene.

    Args:
        sequence : str
            The sequence to extract the exons from.

    Returns:
        tuple[tuple[int, int]]
            A tuple of two tuples. The first tuple contains the start and end positions
            of the first exon, and the second tuple contains the start and end positions
            of the second exon. The positions are 1-based.
    """
    donor_loc, accept_loc = find_ns_splice_sites(sequence)

    # NS2 has 337 additional nts after splice acceptor (AG)
    ns2_end = accept_loc + 337

    # donor_loc corresponds to start of AGGT signal. 'GT' is trimmed, leaving 'AG'.
    # So, the position of the first G is the end of the first exon.
    ns2_exon1_end = donor_loc + 1

    # accept_loc is the start of the AG splice acceptor signal. In splicing the AG is
    # lost. So, need the position of the next nucleotide after the AG.
    ns2_exon2_start = accept_loc + 2

    # +1 for 1-based indexing
    return ((1, ns2_exon1_end + 1), (ns2_exon2_start + 1, ns2_end + 1))


def findall(sub: str, string: str) -> list[int]:
    """
    Return indexes of all substrings in a string
    """
    indexes = []
    length = len(sub)
    for i in range(len(string)):
        if string[i : i + length] == sub:
            indexes.append(i)
    return indexes


def find_ns_splice_donor(seq: str) -> int:
    """
    Find the AGGT splice donor signal. Returns an int which is the index of the 'A'

    Notes:
        In NS1 sequences Gabi has sent the AGGT is on average at position 38.2.
    """
    seq = seq.upper()
    candidates = findall("AGGT", seq)

    if not candidates:
        raise ValueError(f"No NS1 splice donor site ('AGGT') in {seq}")

    # Splice site should be in the first 100 nucs
    candidates = filter(lambda x: x < 100, candidates)

    if not candidates:
        raise ValueError(f"No NS1 splice donor sites ('AGGT') in first 100 nt of {seq}")

    # Choose the site that is closest to position 38.2
    return min(candidates, key=lambda x: abs(x - 38.2))


def find_ns_splice_acceptor(seq: str, donor_loc=None) -> int:
    """
    Find the 'AG' splice acceptor signal. Returns an int which is the index of the 'A'.

    Notes:
        Should be >350 nts downstream of the splice donor location.

    Args:
        seq (str)
        donor_loc (int): Location of the splice donor
    """
    seq = seq.upper()

    donor_loc = find_ns_splice_donor(seq) if donor_loc is None else donor_loc

    # Find all candidate acceptor sites 350 nts downstream of the donor site
    candidates = findall("AG", seq[donor_loc + 350 :])
    candidates = [c + donor_loc + 350 for c in candidates]  # Fix indexing

    # sequence around the splice site should be FQDI
    candidates = list(
        filter(
            lambda x: four_aas_around_splice_site(seq, donor_loc, x) == "FQDI",
            candidates,
        )
    )

    if not candidates:
        raise ValueError(f"No NS splice acceptor sites in {seq}")

    # Pick candidate that is closest to position 507
    return min(candidates, key=lambda x: abs(x - 507))


def splice_ns(seq: str, donor_loc: int, accept_loc: int) -> str:
    """
    Splice an NS sequence given a splice donor and acceptor locations.

    Args:
        seq (str):
        donor_loc (int): Location of the AGGT. The 'AG' remains in the
            transcript, the 'GT' is lost.
        acceptor_loc (int): Location of the 'AG'. The 'AG' is lost.
    """
    seq = seq.upper()

    if (accept_loc - donor_loc) < 350:
        raise ValueError(
            f"Splice acceptor signal location ({accept_loc}) should be at least 350 nts "
            f"downstream of the donor signal ({donor_loc}) location, but it is "
            f"{accept_loc - donor_loc}"
        )

    if seq[donor_loc : donor_loc + 4] != "AGGT":
        raise ValueError(f"No AGGT at position {donor_loc} in {seq}")

    if seq[accept_loc : accept_loc + 2] != "AG":
        raise ValueError(f"No AG at position {accept_loc} in {seq}")

    return seq[: donor_loc + 2] + seq[accept_loc + 2 :]


def four_aas_around_splice_site(seq: str, donor_loc: int, accept_loc: int) -> str:
    """
    What are the four amino acids either side of the splice site given a sequence,
    donor location and acceptor location?
    """
    spliced = splice_ns(seq, donor_loc, accept_loc)
    return sloppy_translate(extract_12_nts_around_splice_site(spliced, donor_loc))


def extract_12_nts_around_splice_site(seq: str, donor_loc: int) -> str:
    """
    Start 4 nts downstream of the donor location

            Start of splice donor 'AGGT' signal
            ⌄     Splice site
            |     ⌄
        XXXXTTTCAG|GA...
        ^
        Extracts 12 nts from here
    """
    start = donor_loc - 4
    end = start + 12
    return seq[start:end]


def find_ns_splice_sites(seq: str) -> tuple[int, int]:
    """
    Lookup the splice donor and acceptor locations for an NS1 transcript.
    """
    donor_loc = find_ns_splice_donor(seq)
    accept_loc = find_ns_splice_acceptor(seq, donor_loc)
    return donor_loc, accept_loc
