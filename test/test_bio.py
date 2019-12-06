#!/usr/bin/env python3

import eremitalpa
import unittest


class TestSloppyTranslate(unittest.TestCase):

    def test_atg(self):
        self.assertEqual('M', eremitalpa.sloppy_translate('ATG'))

    def test_question(self):
        self.assertEqual('X', eremitalpa.sloppy_translate('?TG'))


if __name__ == "__main__":
    unittest.main()
