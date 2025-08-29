from collections import Counter, namedtuple
from itertools import combinations, groupby
from more_itertools import bucket, unzip
from typing import Iterable, Generator, Any
import random
import warnings

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from tqdm import tqdm
import pandas as pd
import matplotlib as mpl

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


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
    "X": "#777777",  # Unknown AA
    "Y": "#a5b8c7",
    "-": "#000000",  # Gap / insertion
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

    Doesn't check that the sequence length is a multiple of three. If any
    'codon' contains any character not in [ACTG] then return X.

    Args:
        sequence (str): Lower or upper case.

    Returns:
        str: The translated sequence.
    """
    sequence = sequence.upper()
    peptide = []
    append = peptide.append
    for i in range(0, len(sequence), 3):
        j = i + 3
        codon = sequence[i:j]
        try:
            amino_acid = TRANSLATION_TABLE[codon]
        except KeyError:
            amino_acid = "X"
        append(amino_acid)
    return "".join(peptide)


def find_mutations(*args, **kwargs):
    raise NotImplementedError("This function is now called 'find_substitutions'")


def find_substitutions(a, b, offset=0):
    """Find mutations between strings a and b.

    Args:
        a (str): The first string.
        b (str): The second string.
        offset (int): An offset to be added to the mutation position.

    Raises:
        ValueError: If lengths of a and b differ.

    Returns:
        tuple[Substitution, ...]: A tuple of Substitution objects.
    """
    if len(a) != len(b):
        raise ValueError("a and b must have same length")

    return tuple(
        Substitution(_a, i + offset, _b)
        for i, (_a, _b) in enumerate(zip(a, b), start=1)
        if _a != _b
    )


class Substitution:
    """A change of a character at a site."""

    def __init__(self, *args):
        """Initializes a Substitution object.

        Instantiate using either 1 or three arguments:
            Substitution("N145K") or Substitution("N", 145, "K")

        Args:
            *args: Either a single string like "N145K" or three arguments
                ("N", 145, "K").

        Raises:
            ValueError: If the number of arguments is not 1 or 3.
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
                "Pass 1 or 3 arguments. E.g. Substitution('N145K') or "
                "Substitution('N', 145, 'K')"
            )
        self._elements = self.a, self.pos, self.b

    def __repr__(self):
        return "Substitution({}, {}, {})".format(self.a, self.pos, self.b)

    def __str__(self):
        return "{}{}{}".format(self.a, self.pos, self.b)

    def __gt__(self, other):
        return (self.pos, self.a, self.b) > (other.pos, other.a, other.b)

    def __lt__(self, other):
        return (self.pos, self.a, self.b) < (other.pos, other.a, other.b)

    def __eq__(self, other):
        return str(self) == str(other)

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
    """Computes the Hamming distance between two sequences.

    Args:
        a (str): The first sequence.
        b (str): The second sequence.
        ignore (Iterable[str]): A string containing characters to ignore.
            Mismatches involving these characters will not contribute to the
            Hamming distance.
        case_sensitive (bool): If True, the comparison is case-sensitive.
        per_site (bool): If True, the Hamming distance is divided by the
            length of the sequences, excluding ignored sites.

    Returns:
        float: The Hamming distance.
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
        length = 0
        for m, n in zip(a, b):
            if (m not in ignore) and (n not in ignore):
                length += 1
                if m != n:
                    d += 1
        try:
            return d / length
        except ZeroDivisionError:
            return 0.0
    else:
        for m, n in zip(a, b):
            if (m != n) and (m not in ignore) and (n not in ignore):
                d += 1
        return float(d)


def hamming_dist_lt(a, b, n, ignore=None):
    """Checks if the Hamming distance between two iterables is less than n.

    This is case-sensitive and does not check if a and b have matching
    lengths.

    Args:
        a (iterable): The first iterable.
        b (iterable): The second iterable.
        n (scalar): The threshold value.
        ignore (set or None): A set of characters to ignore during comparison.

    Returns:
        bool: True if the Hamming distance is less than n, False otherwise.
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
    """Computes all pairwise Hamming distances between items in a collection.

    Args:
        collection (iterable): A collection of sequences.
        ignore (str): Characters to ignore during comparison.
        per_site (bool): Whether to calculate per-site distance.

    Returns:
        list[float]: A list of Hamming distances.
    """
    return [
        hamming_dist(a, b, ignore=ignore, per_site=per_site)
        for a, b in combinations(collection, 2)
    ]


def grouped_sample(population, n, key=None):
    """Randomly samples a population, taking at most n elements from each group.

    Args:
        population (iterable): The population to sample from.
        n (int): The maximum number of samples to take from each group.
        key (callable, optional): A function to group elements by. Defaults to None.

    Returns:
        list: The sampled elements.
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
    """Filters sequences based on Hamming distance.

    Iterates through sequences, excluding those that have a Hamming distance
    of less than n to a sequence already seen.

    Args:
        sequences (iterable[str | Bio.SeqRecord]): The sequences to filter.
        n (int): The Hamming distance threshold.
        progress_bar (bool): Whether to display a progress bar.
        ignore (set, optional): Characters to ignore during comparison. Defaults to None.

    Returns:
        list: The filtered sequences.
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
    """A Counter that handles ties in most_common(1)."""

    def most_common(self, n: int | None = None) -> list[tuple[Any, int]]:
        """Returns the most common elements.

        If n=1 and there is a tie for the most common element, all tied
        elements are returned. Otherwise, it behaves like Counter.most_common.

        Args:
            n (int, optional): The number of most common elements to return.
                Defaults to None.

        Returns:
            list[tuple[Any, int]]: A list of the most common elements and their counts.
        """
        if n == 1:
            max_value = max(self.values())
            return [(k, v) for (k, v) in self.items() if v == max_value]
        else:
            return super().most_common(n)


def _generate_consensus_chars(
    seqs: tuple[str], error_without_strict_majority=True
) -> Generator[str, None, None]:
    """Generates consensus characters from a collection of sequences.

    This is a generator called by consensus_seq.

    Args:
        seqs (tuple[str]): A tuple of sequences. Must be the same length.
        error_without_strict_majority (bool): If True, raises an error if a
            position has a tied most common character. If False, a warning is
            raised and one of the tied characters is chosen.

    Yields:
        str: The next consensus character.
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
    """Computes the consensus of a set of sequences.

    Args:
        seqs (Iterable[str]): The sequences to compute the consensus from.
        case_sensitive (bool): If False, all sequences are converted to
            lowercase.
        **kwds: Additional keyword arguments passed to
            _generate_consensus_chars.

    Returns:
        str: The consensus sequence.
    """
    seqs = tuple(seq.lower() for seq in seqs) if not case_sensitive else tuple(seqs)
    return "".join(_generate_consensus_chars(seqs, **kwds))


def variable_sites(
    seq: pd.Series, max_biggest_prop: float = 0.95, ignore: str = "-X"
) -> Generator[int, None, None]:
    """Finds variable sites among sequences.

    Args:
        seq (pd.Series): A pandas Series of sequences.
        max_biggest_prop (float): The maximum proportion for the most common
            character at a site for it to be considered variable.
        ignore (str): Characters to ignore when calculating proportions.

    Yields:
        int: The 1-indexed position of the next variable site.
    """
    ignore = set(ignore)
    df = seq.str.split("", expand=True)
    for site in df.columns[1:-1]:
        vc = df[site].value_counts().drop(index=ignore, errors="ignore").sort_values()
        biggest_prop = vc.max() / vc.sum()
        if biggest_prop < max_biggest_prop:
            yield site


def load_fasta(
    path: str,
    translate_nt: bool = False,
    convert_to_upper: bool = False,
    start: int = 0,
) -> dict[str, str]:
    """Loads sequences from a FASTA file.

    Args:
        path (str): The path to the FASTA file.
        translate_nt (bool): If True, translate nucleotide sequences to amino
            acids.
        convert_to_upper (bool): If True, convert sequences to uppercase.
        start (int): The 0-based index of the first character to include from
            each sequence. This is applied before translation.

    Returns:
        dict[str, str]: A dictionary mapping sequence descriptions to
            sequences.
    """
    with open(path) as fobj:
        seqs = {
            record.description: (
                sloppy_translate(str(record.seq)[start:])
                if translate_nt
                else str(record.seq)[start:]
            )
            for record in SeqIO.parse(fobj, format="fasta")
        }

    return {k: v.upper() for k, v in seqs.items()} if convert_to_upper else seqs


def write_fasta(path: str, records: dict[str, str]) -> None:
    """Writes sequences to a FASTA file.

    Args:
        path (str): The path to the output FASTA file.
        records (dict[str, str]): A dictionary where keys are sequence headers
            and values are the sequences.
    """
    with open(path, "w") as fobj:
        for header, sequence in records.items():
            fobj.write(f">{header}\n{sequence}\n")


ReferenceAlignment = namedtuple(
    "ReferenceAlignment", ("aligned", "internal_gap_in_ref")
)


def align_to_reference(reference_seq: str, input_seq: str) -> ReferenceAlignment:
    """Aligns an input sequence to a reference sequence.

    Returns the aligned input sequence trimmed to the region of the reference.

    Args:
        reference_seq (str): The reference sequence.
        input_seq (str): The input sequence to align.

    Returns:
        ReferenceAlignment: A named tuple containing the aligned sequence and a
            boolean indicating if an internal gap was introduced in the
            reference.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"

    # Bigger penalty to open gaps
    aligner.target_internal_open_gap_score = -5
    aligner.query_internal_open_gap_score = -5

    # Perform the alignment
    alignments = aligner.align(reference_seq, input_seq)
    best_alignment = alignments[0]

    # Extract the aligned sequences
    # Biopython doesn't provide a nice way of extracting the aligned sequences. This is the
    # 'cleanest' I could come up with.
    aligned_ref, _, aligned_input = best_alignment._format_unicode().strip().split("\n")

    # Check for internal gaps in the reference sequence
    first_non_gap, last_non_gap = idx_first_and_last_non_gap(aligned_ref)
    internal_gap_in_ref = "-" in aligned_ref[first_non_gap:last_non_gap]

    # Trim the input sequence to match the reference length
    # Get the positions corresponding to the non-gap parts of the reference
    trimmed_input_seq = aligned_input[first_non_gap : last_non_gap + 1]

    return ReferenceAlignment(trimmed_input_seq, internal_gap_in_ref)


def idx_first_and_last_non_gap(sequence: str) -> tuple[int, int]:
    """Finds the indices of the first and last non-gap characters.

    If all characters are gaps, the behavior is determined by the loop logic
    (likely resulting in an UnboundLocalError if not handled).

    Args:
        sequence (str): The input sequence, which may contain gaps ('-').

    Returns:
        tuple[int, int]: A tuple containing the indices of the first and last
            non-gap characters.
    """
    for i, char in enumerate(sequence):
        if char != "-":
            first_non_gap = i
            break

    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] != "-":
            last_non_gap = i
            break

    return first_non_gap, last_non_gap


def group_sequences_by_character_at_site(
    seqs: dict[str, str], site: int
) -> dict[str, str]:
    """Groups sequences by the character at a specific site.

    Args:
        seqs (dict[str, str]): A dictionary mapping sequence names to
            sequences.
        site (int): The 1-based site to group by.

    Returns:
        dict[str, str]: A dictionary where keys are characters at the given
            site and values are lists of sequence names.
    """
    bucketed = bucket(seqs.items(), key=lambda x: x[1][site - 1])
    grouped = {}
    for amino_acid in bucketed:
        libs, _ = unzip(bucketed[amino_acid])
        grouped[amino_acid] = list(libs)
    return grouped


def plot_amino_acid_colors(ax: "matplotlib.axes.Axes" = None) -> "matplotlib.axes.Axes":
    """Creates a simple plot to display amino acid colors.

    Args:
        ax (matplotlib.axes.Axes, optional): The matplotlib axes to plot on.
            If None, the current axes are used. Defaults to None.

    Returns:
        matplotlib.axes.Axes: The axes with the plot.
    """
    ax = ax or mpl.pyplot.gca()
    width = 0.5
    for i, (aa, color) in enumerate(reversed(amino_acid_colors.items())):
        patch = mpl.patches.Rectangle((0, i), width=width, height=1, facecolor=color)
        ax.add_artist(patch)
        ax.text(width + 0.01, i + 0.5, aa, va="center")
    ax.set(ylim=(0, i + 1), xlim=(0, 2))
    ax.axis(False)
    return ax
