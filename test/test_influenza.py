#!/usr/bin/env python3

import unittest
from functools import reduce
from operator import add
from pathlib import Path

import eremitalpa as ere


class TestInfluenzaData(unittest.TestCase):
    """
    Generic tests for top level objects.
    """

    def test_cluster_key_residues_expected_length(self):
        self.assertEqual(
            len(ere.influenza.clusters), len(ere.influenza._cluster_key_residues)
        )


class TestAllClusters(unittest.TestCase):
    def test_cluster_sequences_all_len_328(self):
        for cluster in ere.influenza.clusters:
            with self.subTest(cluster=cluster):
                self.assertEqual(328, len(ere.Cluster(cluster).aa_sequence))

    def test_all_clusters_have_color(self):
        for cluster in ere.influenza.clusters:
            cluster.color


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
        seq = (
            "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILD"
            "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
            "NNESFNWTGVTQNGTSSSCKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIW"
            "GVHHPGTNNDQISLYTQASGRITVSTKRSQQTVIPNIGSRPRVRDIPSRISIYWTIVKPG"
            "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRI"
            "TYGACPRYVKQNTLKLATGMRNVPEKQT"
        )
        self.assertEqual("CA04", ere.cluster_from_ha_2(seq))

    def test_fu02_prototype(self):
        seq = (
            "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILD"
            "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
            "NNESFNWTGVTQNGTSSACKRRSNKSFFSRLNWLTHLKYKYPALNVTMPNNEKFDKLYIW"
            "GVHHPGTDSDQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRDVSSRISIYWTIVKPG"
            "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRI"
            "TYGACPRYVKQNTLKLATGMRNVPEKQT"
        )
        self.assertEqual("FU02", ere.cluster_from_ha_2(seq))

    def test_wi05_prototype(self):
        """This sequence is from the WI05 prototype virus."""
        seq = (
            "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILD"
            "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
            "NDESFNWTGVTQNGTSSSCKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIW"
            "GVHHPVTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRIRNIPSRISIYWTIVKPG"
            "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRI"
            "TYGACPRYVKQNTLKLATGMRNVPEKQT"
        )
        self.assertEqual("WI05", ere.cluster_from_ha_2(seq))

    def test_hk14_prototype(self):
        seq = (
            "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILD"
            "GENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEF"
            "NNESFNWTGVTQNGTSSACIRRSSSSFFSRLNWLTHLNYTYPALNVTMPNNEQFDKLYIW"
            "GVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPG"
            "DILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRI"
            "TYGACPRYVKHSTLKLATGMRNVPEKQT"
        )
        self.assertEqual("HK14", ere.cluster_from_ha_2(seq))

    def test_all_gap(self):
        """
        A sequence with only gap characters should raise a NoMatchingKeyResidues error.
        """
        with self.assertRaises(ere.influenza.NoMatchingKeyResidues):
            ere.cluster_from_ha_2("-" * 328)

    def test_all_x(self):
        """
        A sequence of all X characters should raise a NoMatchingKeyResidues error.
        """
        with self.assertRaises(ere.influenza.NoMatchingKeyResidues):
            ere.cluster_from_ha_2("-" * 328)


class TestCluster(unittest.TestCase):
    def test_clusters_have_nt_sequence(self):
        for cluster in ere.influenza.clusters:
            with self.subTest(cluster=cluster):
                self.assertIsInstance(cluster.nt_sequence, str)

    def test_can_only_instantiate_known_clusters(self):
        with self.assertRaisesRegex(ValueError, "unknown cluster:"):
            ere.Cluster("XXXX")

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
        self.assertEqual(3, len(kr))
        self.assertEqual("Y", kr[159])
        self.assertEqual({"R", "G"}, kr[142])
        self.assertEqual("F", kr[193])


class TestHammingToCluster(unittest.TestCase):
    def setUp(self):
        self.seq = (
            "QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILD"
            "GINCTLIDALLGDPHCDVFQDETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEF"
            "ITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIW"
            "GVHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPG"
            "DVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKI"
            "TYGACPKYVKQNTLKLATGMRNVPEKQT"
        )

    def test_exact_match(self):
        self.assertEqual(0, ere.hamming_to_cluster(self.seq, "HK68"))

    def test_raises_error_with_len_mismatch(self):
        with self.assertRaises(ValueError):
            ere.hamming_to_cluster("ABCD", "HK68")


class TestHammingToAllClusters(unittest.TestCase):
    def setUp(self):
        self.seq = (
            "QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILD"
            "GINCTLIDALLGDPHCDVFQDETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEF"
            "ITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIW"
            "GVHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPG"
            "DVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKI"
            "TYGACPKYVKQNTLKLATGMRNVPEKQT"
        )

    def test_len_of_return_value(self):
        """Return value should be length 16 -- the number of clusters."""
        rv = ere.hamming_to_all_clusters(self.seq, strict_len=False)
        self.assertEqual(len(ere.influenza.clusters), len(rv))

    def test_returns_list(self):
        self.assertIsInstance(
            ere.hamming_to_all_clusters(self.seq, strict_len=False), list
        )

    def test_hk68_hd_0(self):
        rv = ere.hamming_to_all_clusters(self.seq, strict_len=False)
        self.assertEqual(0, dict(rv)["HK68"])


class TestClusterTransition(unittest.TestCase):
    def test_unknown_transition(self):
        """
        Unknown transition should generate a ValueError.
        """
        with self.assertRaisesRegex(ValueError, "unrecognised cluster transition"):
            ere.ClusterTransition("HK68", "BK79")

    def test_known_transition_str(self):
        """
        Should be able to pass clusters as strings.
        """
        ere.ClusterTransition("HK68", "EN72")

    def test_known_transition_cluster(self):
        """
        Should be able to pass clusters as Cluster instances.
        """
        ere.ClusterTransition(ere.Cluster("HK68"), ere.Cluster("EN72"))

    def test_hk68_preceding(self):
        """
        HK68 has no preceding transitions.
        """
        result = tuple(ere.ClusterTransition("HK68", "EN72").preceding_transitions)
        self.assertEqual(tuple(), result)

    def test_EN72_preceding(self):
        """
        Simple test case.
        """
        result = tuple(ere.ClusterTransition("EN72", "VI75").preceding_transitions)
        self.assertEqual((ere.ClusterTransition("HK68", "EN72"),), result)

    def test_SI87BE89_before_SI87BE92(self):
        """
        Slightly more complex test case.
        """
        si87_be89 = ere.ClusterTransition("SI87", "BE89")
        si87_be92 = ere.ClusterTransition("SI87", "BE92")
        self.assertIn(si87_be89, tuple(si87_be92.preceding_transitions))

    def test_unpackable(self):
        """
        Should be able to unpack the clusters in a transition.
        """
        c0, c1 = ere.ClusterTransition("SI87", "BE89")
        self.assertEqual(ere.Cluster("SI87"), c0)
        self.assertEqual(ere.Cluster("BE89"), c1)

    def test_from_tuple(self):
        """
        Should be able to make an instance from a tuple.
        """
        ct = ere.ClusterTransition.from_tuple(("HK68", "EN72"))
        self.assertIsInstance(ct, ere.ClusterTransition)


class TestTranslateSegment(unittest.TestCase):
    """Tests for eremitalpa.influenza.translate_segment"""

    segments = "PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"

    @classmethod
    def setUpClass(cls):
        root_dir = Path(ere.__file__).parent.parent
        test_data_dir = root_dir.joinpath("data", "flu", "a_tx37_2024")
        cls.nt_seqs = {
            segment: ere.load_fasta(test_data_dir.joinpath(f"{segment}.fasta"))[
                "EPI_ISL_19027114"
            ]
            for segment in cls.segments
        }

    def test(self):
        """Test calling translate_segment on each segment."""
        for segment in self.segments:
            with self.subTest(segment=segment):
                ere.translate_segment(self.nt_seqs[segment], segment)

    def test_single_stop_at_end(self):
        """Test that there is a single stop codon at the end of each protein."""
        for segment in self.segments:
            proteins = ere.translate_segment(self.nt_seqs[segment], segment)
            for protein, aa_seq in proteins.items():
                with self.subTest(segment=segment, protein=protein):
                    self.assertEqual(aa_seq[-1], "*")
                    if segment != "MP":
                        # Exclude MP which does have an internal stop codon at nt site 757-759
                        self.assertEqual(aa_seq.count("*"), 1)


if __name__ == "__main__":
    unittest.main()
