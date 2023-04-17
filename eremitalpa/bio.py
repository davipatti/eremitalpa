from typing import Iterable, TypeVar, Generator
import random
from collections import Counter
from itertools import combinations, groupby
import warnings

_T = TypeVar("_T")

from tqdm import tqdm

amino_acid_colors = {
    "A": "#F76A05",
    "C": "#dde8cf",
    "D": "#a020f0",
    "E": "#9e806e",
    "F": "#f1b066",
    "G": "#675b2c",
    "H": "#ffc808",
    "I": "#8b8989",
    "K": "#03569b",
    "L": "#9B84AD",
    "M": "#93EDC3",
    "N": "#a2b324",
    "P": "#e9a390",
    "Q": "#742f32",
    "R": "#75ada9",
    "S": "#e72f27",
    "T": "#049457",
    "V": "#00939f",
    "W": "#ED93BD",
    "X": "#777777",  # unknown AA
    "Y": "#a5b8c7",
}


FORWARD_CODON_TABLE = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAT": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGT": "S",
    "ATA": "I",
    "ATC": "I",
    "ATG": "M",
    "ATT": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAT": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAT": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "TAA": "*",
    "TAC": "Y",
    "TAG": "*",
    "TAT": "Y",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TGA": "*",
    "TGC": "C",
    "TGG": "W",
    "TGT": "C",
    "TTA": "L",
    "TTC": "F",
    "TTG": "L",
    "TTT": "F",
}

TRANSLATION_TABLE = dict(FORWARD_CODON_TABLE)
TRANSLATION_TABLE["---"] = "-"


def sloppy_translate(sequence):
    """Translate a nucleotide sequence.

    Don't check that the sequence length is a multiple of three. If any 'codon'
    contains any character not in [ATCG] then return X.

    Args:
        sequence (str): Lower or upper case.

    Returns:
        (str)
    """
    sequence = sequence.upper()
    peptide = ""
    for i in range(0, len(sequence), 3):
        j = i + 3
        codon = sequence[i:j]
        try:
            peptide += TRANSLATION_TABLE[codon]
        except KeyError:
            peptide += "X"
    return peptide


def find_mutations(a, b, offset=0):
    """Find mutations between strings a and b.

    Args:
        a (str)
        b (str)
        offset (int)

    Raises:
        ValueError if lengths of a an b differ.

    Returns:
        list of tuples. tuples are like: ("N", 145, "K") The number indicates
            the 1-indexed position of the mutation. The first element is the a
            character. The last element is the b character.
    """
    if len(a) != len(b):
        raise ValueError("a and b must have same length")

    return tuple(
        Mutation(_a, i + offset, _b)
        for i, (_a, _b) in enumerate(zip(a, b), start=1)
        if _a != _b
    )


class Mutation:
    def __init__(self, *args):
        """Change of a character at a site.

        Instantiate using either 1 or three arguments:
            Mutation("N145K") or Mutation("N", 145, "K")
        """
        if len(args) == 1:
            arg = args[0]
            self.a = arg[0]
            self.pos = int(arg[1:-1])
            self.b = arg[-1]
        elif len(args) == 3:
            self.a = args[0]
            self.pos = int(args[1])
            self.b = args[-1]
        else:
            raise ValueError(
                "Pass 1 or 3 arguments. E.g. Mutation('N145K') or "
                "Mutation('N', 145, 'K')"
            )
        self._elements = (self.a, self.pos, self.b)

    def __repr__(self):
        return "Mutation({}, {}, {})".format(self.a, self.pos, self.b)

    def __str__(self):
        return "{}{}{}".format(self.a, self.pos, self.b)

    def __gt__(self, other):
        return (self.pos, self.a, self.b) > (other.pos, other.a, other.b)

    def __lt__(self, other):
        return (self.pos, self.a, self.b) < (other.pos, other.a, other.b)

    def __eq__(self, other):
        return self._elements == other._elements

    def __getitem__(self, pos):
        return self._elements[pos]

    def __hash__(self):
        return hash(str(self))


def hamming_dist(
    a: str,
    b: str,
    ignore: Iterable[str] = "-X",
    case_sensitive: bool = True,
    per_site: bool = False,
) -> float:
    """
    The hamming distance between a and b.

    Args:
        a: Sequence.
        b: Sequence.
        ignore: String containing characters to ignore. If there is a
            mismatch where one string has a character in ignore, this does not
            contribute to the hamming distance.
        per_site: Divide the hamming distance by the length of a and b,
            minus the number of sites with ignored characters.

    Returns:
        float
    """
    if len(a) != len(b):
        raise ValueError(
            f"Length mismatch ({len(a)} vs. {len(b)}):\n" f"a: {a}\n" f"b: {b}"
        )
    if not case_sensitive:
        a = a.upper()
        b = b.upper()
        ignore = ignore.upper()
    ignore = set(ignore)
    d = 0
    if per_site:
        l = 0
        for m, n in zip(a, b):
            if (m not in ignore) and (n not in ignore):
                l += 1
                if m != n:
                    d += 1
        try:
            return d / l
        except ZeroDivisionError:
            return 0.0
    else:
        for m, n in zip(a, b):
            if (m != n) and (m not in ignore) and (n not in ignore):
                d += 1
        return float(d)


def hamming_dist_lt(a, b, n, ignore=None):
    """
    Test if hamming distance between a and b is less than n. This is case
    sensitive and does not check a and b have matching lengths.

    Args:
        a (iterable)
        b (iterable)
        n (scalar)
        ignore (set or None)

    Returns:
        bool
    """
    ignore = set() if ignore is None else ignore
    hd = 0
    for u, v in zip(a, b):
        if (u != v) and (u not in ignore) and (v not in ignore):
            hd += 1
            if hd >= n:
                return False
    return True


def pairwise_hamming_dists(collection, ignore="-X", per_site=False):
    """Compute all pairwise hamming distances between items in collection.

    Args:
        collection (iterable)

    Returns:
        list of hamming distances
    """
    return [
        hamming_dist(a, b, ignore=ignore, per_site=per_site)
        for a, b in combinations(collection, 2)
    ]


def grouped_sample(population, n, key=None):
    """Randomly sample a population taking at most n elements from each group.

    Args:
        population (iterable)
        n (int): Take at most n samples from each group.
        key (callable): Function by which to group elements. Default (None).

    Returns:
        list
    """
    sample = []
    population = sorted(population, key=key)
    for _, group in groupby(population, key=key):
        group = list(group)
        if len(group) <= n:
            sample += group
        else:
            sample += random.sample(group, n)
    return sample


def filter_similar_hd(sequences, n, progress_bar=False, ignore=None):
    """
    Iterate through sequences excluding those that have a hamming distance of
    less than n to a sequence already seen. Return the non-excluded sequences.

    Args:
        sequences (iterable of str / Bio.SeqRecord)
        progress_bar (bool)
        ignore (set or None)

    Returns:
        list
    """
    if n == 0:
        return list(sequences)

    subset = []
    append = subset.append
    sequences = tqdm(tuple(sequences)) if progress_bar else sequences
    for sequence in sequences:
        for included in subset:
            if hamming_dist_lt(sequence, included, n, ignore=ignore):
                break
        else:
            append(sequence)
    return subset


class TiedCounter(Counter):
    def most_common(self, n: int | None = None) -> list[tuple[_T, int]]:
        """
        If n=1 and there are more than one item that has the maximum count, return all of
        them, not just one. If n is not 1, do the same thing as normal
        Counter.most_common.
        """
        if n == 1:
            max_value = max(self.values())
            return [(k, v) for (k, v) in self.items() if v == max_value]
        else:
            return super().most_common(n)


def _generate_consensus_chars(
    seqs: tuple[str], error_without_strict_majority=True
) -> Generator[str, None, None]:
    """
    Args:
        seqs: Sequences. Must be the same length.
        error_without_strict_majority: Raise an error if a position has a tied most
            common character. If set to False, a warning is raised and a single value
            is chosen.
    """
    if len(set(map(len, seqs))) != 1:
        raise ValueError("seqs differ in length")

    for i in range(len(seqs[0])):
        counts = TiedCounter(seq[i] for seq in seqs)  # E.g. A:10, T:5, C: 2, G:10
        most_common = counts.most_common(1)

        if len(most_common) == 1:
            yield most_common[0][0]
        else:
            msg = f"no strict majority at index {i}: {most_common}"
            if error_without_strict_majority:
                raise ValueError(msg)
            else:
                warnings.warn(msg)
                yield most_common[0][0]


def consensus_seq(seqs: Iterable[str], case_sensitive: bool = True, **kwds) -> str:
    """
    Compute the consensus of sequences.

    Args:
        case_sensitive: If False, all seqs are converted to lowercase.
        **kwds passed to _generate_consensus_chars.
    """
    seqs = tuple(seq.lower() for seq in seqs) if not case_sensitive else tuple(seqs)
    return "".join(_generate_consensus_chars(seqs, **kwds))
