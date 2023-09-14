#!/usr/bin/env python3

from unittest import main, TestCase

import eremitalpa as ere


class TestSplitPairs(TestCase):
    def test_raises_error_if_more_than_2_eq(self):
        """
        Should raise NotImplementedError if values are passed that contain more than 2
        values thar are equivalent.
        """
        with self.assertRaises(NotImplementedError):
            ere.split_pairs((1, 2, 3, 2, 2))

    def test_case_a(self):
        """
        Simple test case.
        """
        self.assertEqual([1, 4.5, 8, 5.5], ere.split_pairs((1, 5, 8, 5)))


if __name__ == "__main__":
    main()
