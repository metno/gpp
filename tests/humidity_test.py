from __future__ import print_function
import unittest
import gridpp
import numpy as np
import os


class Test(unittest.TestCase):
    def test_invalid_input(self):
        """Check that dimension missmatch results in error"""
        with self.assertRaises(Exception) as e:
            gridpp.relative_humidity([293.15], [290, 290])
        with self.assertRaises(Exception) as e:
            gridpp.dewpoint([293.15], [0.9, 0.9])
        with self.assertRaises(Exception) as e:
            gridpp.wetbulb([293.15], [101325], [0.9, 0.9])

    def test_relative_humidity(self):
        # NOTE: RH above > 100 C are 1 in the implementation
        t = [293.15, 293.15, 300, 400]
        td = [293.15, 289.783630, 300, 370]
        rh = [1, 0.817590594291687, 1, 1]
        for i in range(len(t)):
            self.assertAlmostEqual(gridpp.relative_humidity(t[i], td[i]), rh[i], 4)
        np.testing.assert_almost_equal(gridpp.relative_humidity(t, td), rh, 4)

    def test_relative_humidity_invalid(self):
        t = [np.nan, np.nan, 293.15]
        td = [293.15, np.nan, np.nan]
        for i in range(len(t)):
            self.assertTrue(np.isnan(gridpp.relative_humidity(t[i], td[i])))
        self.assertTrue(np.isnan(gridpp.relative_humidity(t, td)).all())

    def test_dewpoint(self):
        t = [293.15, 293.15, 300]
        rh = [1, 0.8, 1]
        td = [293.15, 289.783630, 300]
        for i in range(len(t)):
            self.assertAlmostEqual(gridpp.dewpoint(t[i], rh[i]), td[i], 4)
        np.testing.assert_almost_equal(gridpp.dewpoint(t, rh), td, 4)

    def test_dewpoint_invalid(self):
        t = [np.nan, np.nan, 293.15]
        rh = [293.15, np.nan, np.nan]
        for i in range(len(t)):
            self.assertTrue(np.isnan(gridpp.dewpoint(t[i], rh[i])))
        self.assertTrue(np.isnan(gridpp.dewpoint(t, rh)).all())

    def test_wetbulb(self):
        t = [270, 300, 270, 240]
        p = [100000, 101000, 100000, 50000]
        rh = [0.8, 0.7, 1, 0.9]
        ans = [269.02487,296.13763,269.92218,239.83798]
        for i in range(len(t)):
            self.assertAlmostEqual(gridpp.wetbulb(t[i], p[i], rh[i]), ans[i], 4)
        np.testing.assert_almost_equal(gridpp.wetbulb(t, p, rh), ans, 4)

    def test_wetbulb_invalid(self):
        t = [np.nan, np.nan, 293.15, 293.15, np.nan, 293.15]
        p = [101325, 101325, 101325, np.nan, np.nan, np.nan]
        rh = [0.9, np.nan, np.nan, 0.9, np.nan, 0]
        for i in range(len(t)):
            self.assertTrue(np.isnan(gridpp.wetbulb(t[i], p[i], rh[i])))
        self.assertTrue(np.isnan(gridpp.wetbulb(t, p, rh)).all())


if __name__ == '__main__':
    unittest.main()
