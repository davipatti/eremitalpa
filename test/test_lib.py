#!/usr/bin/env python3

from unittest import main, TestCase

import pandas as pd

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


class TestCalendarMonthsDiff(TestCase):
    def test_case_a(self):
        self.assertEqual(
            0,
            ere.cal_months_diff(pd.Timestamp("2016-01-01"), pd.Timestamp("2016-01-01")),
        )

    def test_case_b(self):
        self.assertEqual(
            0,
            ere.cal_months_diff(pd.Timestamp("2016-01-10"), pd.Timestamp("2016-01-01")),
        )

    def test_case_c(self):
        self.assertEqual(
            1,
            ere.cal_months_diff(pd.Timestamp("2016-02-01"), pd.Timestamp("2016-01-01")),
        )

    def test_case_d(self):
        self.assertEqual(
            13,
            ere.cal_months_diff(pd.Timestamp("2017-02-01"), pd.Timestamp("2016-01-01")),
        )

    def test_case_e(self):
        self.assertEqual(
            13,
            ere.cal_months_diff(pd.Timestamp("2017-02-01"), pd.Timestamp("2016-01-31")),
        )


if __name__ == "__main__":
    main()
