#!/usr/bin/env python3

import unittest
import eremitalpa as ere
from functools import reduce
from operator import add


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


if __name__ == "__main__":
    unittest.main()
