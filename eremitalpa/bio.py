FORWARD_CODON_TABLE = {
    'AAA': 'K',
    'AAC': 'N',
    'AAG': 'K',
    'AAT': 'N',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'AGA': 'R',
    'AGC': 'S',
    'AGG': 'R',
    'AGT': 'S',
    'ATA': 'I',
    'ATC': 'I',
    'ATG': 'M',
    'ATT': 'I',
    'CAA': 'Q',
    'CAC': 'H',
    'CAG': 'Q',
    'CAT': 'H',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'GAA': 'E',
    'GAC': 'D',
    'GAG': 'E',
    'GAT': 'D',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    'TAA': '*',
    'TAC': 'Y',
    'TAG': '*',
    'TAT': 'Y',
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'TGA': '*',
    'TGC': 'C',
    'TGG': 'W',
    'TGT': 'C',
    'TTA': 'L',
    'TTC': 'F',
    'TTG': 'L',
    'TTT': 'F',
}


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
    i = 0
    peptide = ''
    while i <= len(sequence) - 3:
        j = i + 3
        codon = sequence[i:j]
        try:
            peptide += FORWARD_CODON_TABLE[codon]
        except KeyError:
            if codon == '---':
                peptide += '-'
            else:
                peptide += 'X'
        i = j
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
        list of tuples. tuples are like: ('N', 145, 'K') The number indicates
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


class Mutation():

    def __init__(self, a, pos, b):
        self.a = a
        self.pos = int(pos)
        self.b = b
        self._elements = (a, pos, b)

    def __repr__(self):
        return "Mutation(a={}, pos={}, b={})".format(self.a, self.pos, self.b)

    def __str__(self):
        return "{}{}{}".format(self.a, self.pos, self.b)

    def __getitem__(self, pos):
        return self._elements[pos]

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self._elements == other._elements
