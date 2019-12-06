#!/usr/bin/env python3

import unittest
import eremitalpa as ere


class TestSloppyTranslate(unittest.TestCase):

    def test_atg(self):
        self.assertEqual('M', ere.sloppy_translate('ATG'))

    def test_question(self):
        self.assertEqual('X', ere.sloppy_translate('?TG'))


class TestFindMutations(unittest.TestCase):

    def test_no_mutations(self):
        self.assertEqual(0, len(ere.find_mutations('ABC', 'ABC')))

    def test_one_mutation(self):
        self.assertEqual([["D", 3, "C"]], ere.find_mutations("ABD", "ABC"))

    def test_raises_with_len_mismatch(self):
        with self.assertRaises(ValueError):
            ere.find_mutations("ABC", "AB")


if __name__ == "__main__":
    unittest.main()
