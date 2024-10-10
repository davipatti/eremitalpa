#!/usr/bin/env python3

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import eremitalpa as ere


class TestMultipleSequenceAlignmentVariableSites(unittest.TestCase):
    """
    Tests for eremitalpa.MultipleSequenceAlignment.variable_sites
    """

    def test_no_variable_sites(self):
        """
        If there are no variable sites in a multiple sequence alignment, variable_sites
        should return an empty generator.
        """
        msa = ere.MultipleSequenceAlignment([SeqRecord(Seq("atcg")) for _ in range(5)])
        self.assertEqual((), tuple(msa.variable_sites()))

    def test_one_variable_site(self):
        """
        Test a single variable site.
        """
        msa = ere.MultipleSequenceAlignment(
            [
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcc")),
            ]
        )
        vs = tuple(msa.variable_sites())
        self.assertEqual(1, len(vs))

    def test_variable_sites_contain_Column_instances(self):
        """
        Variable sites should be Column instances.
        """
        msa = ere.MultipleSequenceAlignment(
            [
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcc")),
            ]
        )
        column = next(msa.variable_sites())
        self.assertIsInstance(column, ere.eremitalpa.Column)

    def test_column_site(self):
        """
        Site should be 1-indexed.
        """
        msa = ere.MultipleSequenceAlignment(
            [
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcc")),
            ]
        )
        column = next(msa.variable_sites())
        self.assertEqual(4, column.site)

    def test_min_2nd_most_freq_2_casea(self):
        """
        In this test case, passing min_2nd_most_freq=2 should exclude all columns.
        """
        msa = ere.MultipleSequenceAlignment(
            [
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcg")),
                SeqRecord(Seq("atcc")),
            ]
        )
        vs = tuple(msa.variable_sites(min_2nd_most_freq=2))
        self.assertEqual(0, len(vs))

    def test_min_2nd_most_freq_2_caseb(self):
        """
        Here, passing min_2nd_most_freq=2 should yield only a single column.
        """
        msa = ere.MultipleSequenceAlignment(
            [
                SeqRecord(Seq("atcga")),
                SeqRecord(Seq("atcga")),
                SeqRecord(Seq("atcgt")),
                SeqRecord(Seq("atcct")),
            ]
        )
        vs = tuple(msa.variable_sites(min_2nd_most_freq=2))
        self.assertEqual(1, len(vs))


if __name__ == "__main__":
    unittest.main()
