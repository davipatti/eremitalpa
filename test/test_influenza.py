#!/usr/bin/env python3

import unittest
from functools import reduce
from operator import add

import eremitalpa as ere


class TestClassifyCluster(unittest.TestCase):

    def test_cluster_motifs_all_unique(self):
        motifs = reduce(add, [list(v) for v in ere._cluster_motifs.values()])
        for motif in motifs:
            self.assertEqual(1, motifs.count(motif),
                             "{} occurs more than once in "
                             "_cluster_motifs".format(motif))

    def test_cluster_motifs_are_tuples(self):
        for v in ere._cluster_motifs.values():
            self.assertIsInstance(v, tuple)

    def test_cluster_motif_tuples_contain_str(self):
        for v in ere._cluster_motifs.values():
            for motif in v:
                self.assertIsInstance(motif, str)
                self.assertEqual(7, len(motif))

    def test_hk68(self):
        self.assertEqual("HK68", ere.cluster_from_ha("STKGSQS", seq_type="b7"))

    def test_raises_if_cant_classify(self):
        with self.assertRaises(ValueError):
            ere.cluster_from_ha("DAVIDPA", seq_type="b7")


class TestHammingDistance(unittest.TestCase):

    def test_empty_str(self):
        self.assertEqual(0, ere.hamming_dist("", ""))

    def test_1_mismatch(self):
        self.assertEqual(1, ere.hamming_dist("A", "B"))

    def test_case_insensititve(self):
        self.assertEqual(0, ere.hamming_dist("A", "a"))

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
            "D-VIDPATTINSON", "NAVIDPaTTIXSON"))

    def test_per_site_single(self):
        self.assertEqual(1, ere.hamming_dist("A", "B", per_site=True))

    def test_per_site_double(self):
        self.assertEqual(0.5, ere.hamming_dist("AB", "AA", per_site=True))

    def test_per_site_longer(self):
        self.assertEqual(1 / 12 , ere.hamming_dist(
            "D-VIDPATTINSON", "NAVIDPaTTIXSON", per_site=True))

    def test_hk68_year(self):
        self.assertEqual(1968, ere.Cluster("HK68").year)

    def test_fuo2_year(self):
        self.assertEqual(2002, ere.Cluster("FU02").year)

    def test_vi75_gt_hk68(self):
        self.assertGreater(ere.Cluster("VI75"), ere.Cluster("HK68"))

    def test_fu02_gt_sy97(self):
        self.assertGreater(ere.Cluster("FU02"), ere.Cluster("SY97"))


if __name__ == "__main__":
    unittest.main()
