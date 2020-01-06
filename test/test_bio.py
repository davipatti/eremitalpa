#!/usr/bin/env python3

import unittest

import eremitalpa as ere


class TestSloppyTranslate(unittest.TestCase):

    def test_atg(self):
        self.assertEqual('M', ere.sloppy_translate('ATG'))

    def test_lower_case(self):
        self.assertEqual('M', ere.sloppy_translate('atg'))

    def test_two_codon(self):
        self.assertEqual('FA', ere.sloppy_translate('TTTGCT'))

    def test_seven_bases(self):
        """Trailing nucleotides should generate trailing Xs."""
        self.assertEqual('FAX', ere.sloppy_translate('TTTGCTA'))

    def test_eight_bases(self):
        """Trailing nucleotides should generate trailing Xs."""
        self.assertEqual('FAX', ere.sloppy_translate('TTTGCTAC'))

    def test_question(self):
        self.assertEqual('X', ere.sloppy_translate('?TG'))


class TestFindMutations(unittest.TestCase):

    def test_no_mutations(self):
        self.assertEqual(0, len(ere.find_mutations('ABC', 'ABC')))

    def test_one_mutation(self):
        self.assertEqual((ere.Mutation("D", 3, "C"),),
                         ere.find_mutations("ABD", "ABC"))

    def test_raises_with_len_mismatch(self):
        with self.assertRaises(ValueError):
            ere.find_mutations("ABC", "AB")


class TestMutation(unittest.TestCase):

    def test_from_3_args(self):
        m = ere.Mutation("N", 145, "K")
        self.assertEqual("N", m.a)
        self.assertEqual("K", m.b)
        self.assertEqual(145, m.pos)

    def test_from_1_arg(self):
        m = ere.Mutation("N145K")
        self.assertEqual("N", m.a)
        self.assertEqual("K", m.b)
        self.assertEqual(145, m.pos)

    def test_two_instances_same_hash(self):
        a = ere.Mutation("N145K")
        b = ere.Mutation("N145K")
        self.assertEqual(hash(a), hash(b))


class TestHammingDistance(unittest.TestCase):

    def test_empty_str(self):
        self.assertEqual(0, ere.hamming_dist("", ""))

    def test_1_mismatch(self):
        self.assertEqual(1, ere.hamming_dist("A", "B"))

    def test_case_insensititve(self):
        self.assertEqual(0, ere.hamming_dist("A", "a", case_sensitive=False))

    def test_ignore_X(self):
        self.assertEqual(0, ere.hamming_dist("A", "X"))

    def test_ignore_gap(self):
        self.assertEqual(0, ere.hamming_dist("-", "A"))

    def test_ignore_argument(self):
        self.assertEqual(0, ere.hamming_dist("A", "B", ignore="A"))

    def test_len_mismatch_raises_valueerror(self):
        with self.assertRaises(ValueError):
            ere.hamming_dist("A", "AB")

    def test_longer(self):
        self.assertEqual(1, ere.hamming_dist(
            "D-VIDPATTINSON", "NAVIDPaTTIXSON", case_sensitive=False))

    def test_per_site_single(self):
        self.assertEqual(1, ere.hamming_dist("A", "B", per_site=True))

    def test_per_site_double(self):
        self.assertEqual(0.5, ere.hamming_dist("AB", "AA", per_site=True))

    def test_per_site_longer(self):
        self.assertEqual(1 / 12, ere.hamming_dist(
            "D-VIDPATTINSON", "NAVIDPaTTIXSON", per_site=True,
            case_sensitive=False))

    def test_per_site_zero(self):
        self.assertEqual(0, ere.hamming_dist("A", "A", per_site=True))


class TestPairwiseHammingDist(unittest.TestCase):

    def test_returns_list(self):
        collection = ("DAVID", "PAVID", "DAVID")
        self.assertIsInstance(ere.pairwise_hamming_dists(collection), list)


if __name__ == "__main__":
    unittest.main()
