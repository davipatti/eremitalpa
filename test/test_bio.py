#!/usr/bin/env python3

import unittest
from collections import Counter
import random
from operator import itemgetter
from string import ascii_lowercase

from Bio.Seq import Seq

import eremitalpa as ere


class TestSloppyTranslate(unittest.TestCase):
    def test_atg(self):
        self.assertEqual("M", ere.sloppy_translate("ATG"))

    def test_lower_case(self):
        self.assertEqual("M", ere.sloppy_translate("atg"))

    def test_two_codon(self):
        self.assertEqual("FA", ere.sloppy_translate("TTTGCT"))

    def test_seven_bases(self):
        """Trailing nucleotides should generate trailing Xs."""
        self.assertEqual("FAX", ere.sloppy_translate("TTTGCTA"))

    def test_eight_bases(self):
        """Trailing nucleotides should generate trailing Xs."""
        self.assertEqual("FAX", ere.sloppy_translate("TTTGCTAC"))

    def test_question(self):
        """? should generate X"""
        self.assertEqual("X", ere.sloppy_translate("?TG"))

    def test_gap(self):
        """- should generate X"""
        self.assertEqual("X", ere.sloppy_translate("AT-"))

    def test_gap_question(self):
        """Combination of gaps and question marks."""
        self.assertEqual("MXX", ere.sloppy_translate("ATG-TG?TG"))


class TestFindSubstitutions(unittest.TestCase):
    def test_no_substitutions(self):
        self.assertEqual(0, len(ere.find_substitutions("ABC", "ABC")))

    def test_one_substitution(self):
        self.assertEqual(
            (ere.Substitution("D", 3, "C"),), ere.find_substitutions("ABD", "ABC")
        )

    def test_raises_with_len_mismatch(self):
        with self.assertRaises(ValueError):
            ere.find_substitutions("ABC", "AB")


class TestSubstitution(unittest.TestCase):
    def test_from_3_args(self):
        m = ere.Substitution("N", 145, "K")
        self.assertEqual("N", m.a)
        self.assertEqual("K", m.b)
        self.assertEqual(145, m.pos)

    def test_from_1_arg(self):
        m = ere.Substitution("N145K")
        self.assertEqual("N", m.a)
        self.assertEqual("K", m.b)
        self.assertEqual(145, m.pos)

    def test_two_instances_same_hash(self):
        a = ere.Substitution("N145K")
        b = ere.Substitution("N145K")
        self.assertEqual(hash(a), hash(b))


class TestHammingDistance(unittest.TestCase):
    def test_empty_str(self):
        """
        Two empty strings have a HD of 0.
        """
        self.assertEqual(0, ere.hamming_dist("", ""))

    def test_passed_tuple(self):
        """
        Should handle passing an iterable that isn't a string.
        """
        ere.hamming_dist(tuple("abcd"), tuple("abcd"))

    def test_1_mismatch(self):
        """
        Simple test case of HD == 1.
        """
        self.assertEqual(1, ere.hamming_dist("A", "B"))

    def test_case_insensititve(self):
        """
        Test the case_sensitive flag.
        """
        self.assertEqual(0, ere.hamming_dist("A", "a", case_sensitive=False))

    def test_case_sensititve(self):
        """
        Test the case_sensitive flag.
        """
        self.assertEqual(1, ere.hamming_dist("A", "a", case_sensitive=True))

    def test_ignore_X(self):
        """
        X should be ignored by default.
        """
        self.assertEqual(0, ere.hamming_dist("A", "X"))

    def test_ignore_gap(self):
        """
        Gaps (-) should be ignored by default.
        """
        self.assertEqual(0, ere.hamming_dist("-", "A"))

    def test_ignore_argument(self):
        """
        Specify custom ignore characters.
        """
        self.assertEqual(0, ere.hamming_dist("A", "B", ignore="A"))

    def test_len_mismatch_raises_valueerror(self):
        """
        Should not be able to pass arguments of different length.
        """
        with self.assertRaises(ValueError):
            ere.hamming_dist("A", "AB")

    def test_longer(self):
        """
        More complex test case.
        """
        self.assertEqual(
            1,
            ere.hamming_dist("D-VIDPATTINSON", "NAVIDPaTTIXSON", case_sensitive=False),
        )

    def test_per_site_single(self):
        """
        Test per site flag with 100% mismatch.
        """
        self.assertEqual(1, ere.hamming_dist("A", "B", per_site=True))

    def test_per_site_double(self):
        """
        Test per site flag with 50% mismatch.
        """
        self.assertEqual(0.5, ere.hamming_dist("AB", "AA", per_site=True))

    def test_per_site_longer(self):
        """
        Test per site flag with longer test case.
        """
        self.assertEqual(
            1 / 12,
            ere.hamming_dist(
                "D-VIDPATTINSON", "NAVIDPaTTIXSON", per_site=True, case_sensitive=False
            ),
        )

    def test_per_site_zero(self):
        """
        Test per site flat with 0 mismatches.
        """
        self.assertEqual(0, ere.hamming_dist("A", "A", per_site=True))


class TestPairwiseHammingDist(unittest.TestCase):
    def test_returns_list(self):
        collection = ("DAVID", "PAVID", "DAVID")
        self.assertIsInstance(ere.pairwise_hamming_dists(collection), list)


class TestGroupedSample(unittest.TestCase):
    def test_returns_list(self):
        """
        Should return a list.
        """
        self.assertIsInstance(ere.grouped_sample("abcd", n=1), list)

    def test_asking_for_at_most_one_item(self):
        """
        Construct case where grouping by the second item, and asking for at
        most 1 item per group should return a smaller population.
        """
        items = (("a", 1), ("b", 1), ("c", 2))
        rv = ere.grouped_sample(items, n=1, key=itemgetter(1))
        self.assertEqual(2, len(rv))

    def test_ask_for_at_most_two(self):
        """
        Construct case where grouping by the second item, and asking for at
        most 2 item per group should return the same population.
        """
        items = (("a", 1), ("b", 1), ("c", 2))
        rv = ere.grouped_sample(items, n=2, key=itemgetter(1))
        self.assertEqual(3, len(rv))

    def test_no_oversampling(self):
        """
        Asking for more items than are in any one group should not cause
        oversampling.
        """
        items = (("a", 1), ("b", 1), ("c", 2))
        rv = ere.grouped_sample(items, n=10, key=itemgetter(1))
        self.assertEqual(3, len(rv))

    def test_asking_for_no_items(self):
        """
        Asking for at most 0 items per group should return no elements.
        """
        items = (("a", 1), ("b", 1), ("c", 2))
        rv = ere.grouped_sample(items, n=0, key=itemgetter(1))
        self.assertEqual(0, len(rv))

    def test_setting_seed_causes_repeatable(self):
        """
        Setting random.seed should cause the results to be repeatable.
        """
        items = [
            (random.choice(ascii_lowercase), random.choice(range(10)))
            for _ in range(100)
        ]

        random.seed(1234)
        a = ere.grouped_sample(items, n=2, key=itemgetter(1))

        random.seed(1234)
        b = ere.grouped_sample(items, n=2, key=itemgetter(1))

        self.assertEqual(a, b)

    def test_parity(self):
        """
        Test case sampling at most 2 odd or even numbers.
        """
        gs = ere.grouped_sample(range(10), n=2, key=lambda x: x % 2)
        self.assertEqual(4, len(gs))


class TestFilterSimilarHD(unittest.TestCase):
    def test_returns_list(self):
        """
        Should return a list.
        """
        sequences = "abc", "abc"
        result = ere.filter_similar_hd(sequences, 1)
        self.assertIsInstance(result, list)

    def test_n_zero(self):
        """
        n == 0 should mean that multiple identical sequences are returned.
        """
        sequences = "abc", "abc"
        result = ere.filter_similar_hd(sequences, 0)
        self.assertEqual(["abc", "abc"], result)

    def test_n_one(self):
        """
        n == 1 should mean that only a single copy of groups of identical
        sequences are returned.
        """
        sequences = "abc", "abc", "abc"
        result = ere.filter_similar_hd(sequences, 1)
        self.assertEqual(["abc"], result)

    def test_n_two(self):
        """
        n == 2 should mean that sequences with one mismatch to something already
        seen are not returned. Sequences with 2 or more mismatches should be
        returned.
        """
        sequences = (
            "abc",
            "abd",  # HD < 2 to 'abc' --> should not be returned
            "ade",  # HD ! <2 to 'abc' --> should be returned
        )
        result = ere.filter_similar_hd(sequences, 2)
        self.assertEqual(["abc", "ade"], result)

    def test_hamming_dist_kws_are_passed(self):
        """
        Test ignore keyword.
        """
        sequences = "abc", "ade"
        result = ere.filter_similar_hd(sequences, 1, ignore="D")
        # 'd' is ignored so these sequences have an HD of 1
        # therefore both seqs should be returned.
        self.assertEqual(["abc", "ade"], result)

    def test_biopython_seq_records(self):
        """
        Test using Bio.Seq.Seq.
        """
        sequences = Seq("ACTG"), Seq("ACTT")
        result = ere.filter_similar_hd(sequences, 2)
        self.assertEqual(1, len(result))


class TestHammingDistanceLt(unittest.TestCase):
    def test_returns_bool(self):
        """
        Should return a boolean.
        """
        self.assertIsInstance(ere.hamming_dist_lt("a", "b", 1), bool)

    def test_case_a(self):
        """
        Strings that differ by 1 have a hamming distance less than 2.
        """
        self.assertTrue(ere.hamming_dist_lt("a", "b", 2))

    def test_case_b(self):
        """
        Identical strings have a hamming distance less than 1.
        """
        result = ere.hamming_dist_lt("DavidPattinson", "DavidPattinson", 1)
        self.assertTrue(result)

    def test_ignore_argument(self):
        """
        Ignored characters shouldn't contribute to HD.
        """
        self.assertTrue(ere.hamming_dist_lt("abfgh", "abcde", 1, ignore="fgh"))

    def test_case_sensititve(self):
        """
        Test the function is case sensitive.
        """
        self.assertFalse(ere.hamming_dist_lt("abc", "ABC", 1))


class TestConsensusSeq(unittest.TestCase):
    """
    Tests for eremitalpa.bio.consensus_seq.
    """

    def test_returns_str(self):
        """
        Should return a str
        """
        self.assertIsInstance(ere.consensus_seq(["abc", "abc"]), str)

    def test_consensus_identical_seqs(self):
        """
        The consensus of identical sequences should be one of the sequences.
        """
        self.assertEqual("abc", ere.consensus_seq(["abc", "abc"]))

    def test_warns_no_unique_consensus(self):
        """
        Warning is raised with strict_majority=False if there is no unique consensus.

        (E.g. with matched character counts at one or more positions).
        """
        with self.assertWarnsRegex(Warning, "no strict majority at index"):
            ere.consensus_seq(["abc", "abd"], error_without_strict_majority=False)

    def test_value_error_no_unique_consensus(self):
        """
        ValueError is raised with strict_majority=True (default value) if there is no
        unique consensus.

        (E.g. with matched character counts at one or more positions).
        """
        with self.assertRaisesRegex(ValueError, "no strict majority at index"):
            ere.consensus_seq(["abc", "abd"])

    def test_case_insensitive(self):
        """
        Check the case_insensitive flag works correctly.
        """
        self.assertEqual("abc", ere.consensus_seq(["abc", "ABC"], case_sensitive=False))

    def test_length_mismatch(self):
        """
        A ValueError should be raised if sequences of different lengths are passed.
        """
        with self.assertRaisesRegex(ValueError, "seqs differ in length"):
            ere.consensus_seq(["abc", "ab"])

    def test_case_three_sequences(self):
        """
        Test case with more sequences.
        """
        seqs = ["ABCDEFG", "ABCXXXX", "ABCDEFX"]
        self.assertEqual("ABCDEFX", ere.consensus_seq(seqs))


class TestTiedCounter(unittest.TestCase):
    """Tests for eremitalpa.bio.TiedCounter."""

    def test_regular_counter_tied_count(self):
        """
        collections.Counter returns only a single most common item, even if there are
        multiple most common items.
        """
        self.assertEqual(
            [
                ("a", 2),
            ],
            Counter("ababc").most_common(1),
        )

    def test_tied_count_n_eq_1(self):
        """
        With multiple most common items, TiedCounter.most_common should return
        all most common items.
        """
        self.assertEqual(
            [
                ("a", 2),
                ("b", 2),
            ],
            ere.TiedCounter("ababc").most_common(1),
        )

    def test_tied_count_n_None(self):
        """
        With n as None, ere.TiedCounter should behave like collections.Counter.
        """
        self.assertEqual(
            Counter("ababc").most_common(),
            ere.TiedCounter("ababc").most_common(),
        )

    def test_tied_count_n_2(self):
        """
        With n as 2, ere.TiedCounter should behave like collections.Counter.
        """
        self.assertEqual(
            Counter("ababc").most_common(2),
            ere.TiedCounter("ababc").most_common(2),
        )


if __name__ == "__main__":
    unittest.main()
