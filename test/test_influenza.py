#!/usr/bin/env python3

import unittest
from functools import reduce
from operator import add

import eremitalpa as ere


class TestInfluenzaData(unittest.TestCase):
    def test_cluster_key_residues_expected_length(self):
        self.assertEqual(
            len(ere.influenza.clusters), len(ere.influenza._cluster_key_residues)
        )


class TestClassifyCluster(unittest.TestCase):
    def test_cluster_motifs_all_unique(self):
        motifs = reduce(add, [list(v) for v in ere._cluster_motifs.values()])
        for motif in motifs:
            self.assertEqual(
                1,
                motifs.count(motif),
                "{} occurs more than once in " "_cluster_motifs".format(motif),
            )

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


class TestClassifyCluster2(unittest.TestCase):
    def test_ca04_prototype(self):
        seq = "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSSCKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTNNDQISLYTQASGRITVSTKRSQQTVIPNIGSRPRVRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT"
        self.assertEqual("CA04", ere.cluster_from_ha_2(seq))

    def test_fu02_prototype(self):
        seq = "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRRSNKSFFSRLNWLTHLKYKYPALNVTMPNNEKFDKLYIWGVHHPGTDSDQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRDVSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT"
        self.assertEqual("FU02", ere.cluster_from_ha_2(seq))

    def test_wi05_prototype(self):
        """This sequence is from the WI05 prototype virus."""
        seq = "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNDESFNWTGVTQNGTSSSCKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPVTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRIRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT"
        self.assertEqual("WI05", ere.cluster_from_ha_2(seq))


class TestCluster(unittest.TestCase):
    def test_hk68_year(self):
        self.assertEqual(1968, ere.Cluster("HK68").year)

    def test_fuo2_year(self):
        self.assertEqual(2002, ere.Cluster("FU02").year)

    def test_vi75_gt_hk68(self):
        self.assertGreater(ere.Cluster("VI75"), ere.Cluster("HK68"))

    def test_fu02_gt_sy97(self):
        self.assertGreater(ere.Cluster("FU02"), ere.Cluster("SY97"))

    def test_key_residues_HK68(self):
        kr = ere.Cluster("HK68").key_residues
        self.assertIsInstance(kr, dict)
        self.assertEqual(1, len(kr))
        self.assertEqual("T", kr[155])

    def test_key_residues_EN72(self):
        kr = ere.Cluster("EN72").key_residues
        self.assertEqual(2, len(kr))
        self.assertEqual("Y", kr[155])
        self.assertEqual("Q", kr[189])

    def test_key_residues_TX77(self):
        kr = ere.Cluster("TX77").key_residues
        self.assertEqual(3, len(kr))
        self.assertEqual("E", kr[158])
        self.assertEqual("N", kr[193])
        self.assertEqual("K", kr[156])

    def test_key_residues_BK79(self):
        kr = ere.Cluster("BK79").key_residues
        self.assertEqual(4, len(kr))
        self.assertEqual("E", kr[156])
        self.assertEqual("Y", kr[155])
        self.assertEqual("S", kr[159])
        self.assertEqual("K", kr[189])

    def test_key_residues_SI87(self):
        kr = ere.Cluster("SI87").key_residues
        self.assertEqual(5, len(kr))
        self.assertEqual("H", kr[155])
        self.assertEqual("Y", kr[159])
        self.assertEqual("R", kr[189])
        self.assertEqual("N", kr[145])
        self.assertEqual("E", kr[156])

    def test_key_residues_BE89(self):
        kr = ere.Cluster("BE89").key_residues
        self.assertEqual(1, len(kr))
        self.assertEqual("K", kr[145])

    def test_key_residues_BE92(self):
        kr = ere.Cluster("BE92").key_residues
        self.assertEqual(2, len(kr))
        self.assertEqual("K", kr[156])
        self.assertEqual("N", kr[145])

    def test_key_residues_WU95(self):
        kr = ere.Cluster("WU95").key_residues
        self.assertEqual(3, len(kr))
        self.assertEqual("K", kr[145])
        self.assertEqual("K", kr[156])
        self.assertEqual("E", kr[158])

    def test_key_residues_SY97(self):
        kr = ere.Cluster("SY97").key_residues
        self.assertEqual(2, len(kr))
        self.assertEqual("Q", kr[156])
        self.assertEqual("K", kr[158])

    def test_key_residues_FU02(self):
        kr = ere.Cluster("FU02").key_residues
        self.assertEqual(2, len(kr))
        self.assertEqual("H", kr[156])
        self.assertEqual("K", kr[145])

    def test_key_residues_CA04(self):
        kr = ere.Cluster("CA04").key_residues
        self.assertEqual(2, len(kr))
        self.assertEqual("N", kr[145])
        self.assertEqual("S", kr[193])

    def test_key_residues_WI05(self):
        kr = ere.Cluster("WI05").key_residues
        self.assertEqual(3, len(kr))
        self.assertEqual("F", kr[193])
        self.assertEqual("K", kr[158])
        self.assertEqual("N", kr[189])

    def test_key_residues_PE09(self):
        kr = ere.Cluster("PE09").key_residues
        self.assertEqual(3, len(kr))
        self.assertEqual("N", kr[158])
        self.assertEqual("K", kr[189])
        self.assertEqual("F", kr[159])

    def test_key_residues_SW13(self):
        kr = ere.Cluster("SW13").key_residues
        self.assertEqual(2, len(kr))
        self.assertEqual("S", kr[159])
        self.assertEqual("F", kr[193])

    def test_key_residues_HK14(self):
        kr = ere.Cluster("HK14").key_residues
        self.assertEqual(4, len(kr))
        self.assertEqual("Y", kr[159])
        self.assertEqual({"R", "G"}, kr[142])
        self.assertEqual("F", kr[193])
        self.assertEqual("T", kr[131])


class TestHammingToCluster(unittest.TestCase):
    def setUp(self):
        self.seq = "QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGINCTLIDALLGDPHCDVFQDETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGVHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT"

    def test_exact_match(self):
        self.assertEqual(0, ere.hamming_to_cluster(self.seq, "HK68"))

    def test_raises_error_with_len_mismatch(self):
        with self.assertRaises(ValueError):
            ere.hamming_to_cluster("ABCD", "HK68")


class TestHammingToAllClusters(unittest.TestCase):
    def setUp(self):
        self.seq = "QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGINCTLIDALLGDPHCDVFQDETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGVHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT"

    def test_len_of_return_value(self):
        """Return value should be length 16 -- the number of clusters."""
        rv = ere.hamming_to_all_clusters(self.seq)
        self.assertEqual(16, len(rv))

    def test_returns_list(self):
        self.assertIsInstance(ere.hamming_to_all_clusters(self.seq), list)

    def test_hk68_hd_0(self):
        rv = ere.hamming_to_all_clusters(self.seq)
        self.assertEqual(0, dict(rv)["HK68"])


if __name__ == "__main__":
    unittest.main()
