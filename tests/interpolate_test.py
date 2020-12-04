from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_basic(self):
        x = [0, 1, 2]
        y = [0, 2, 1]

        # Test numbers in interpolator points
        self.assertEqual(gridpp.interpolate(0, x, y), 0)
        self.assertEqual(gridpp.interpolate(2, x, y), 1)
        self.assertEqual(gridpp.interpolate(1, x, y), 2)

        # Test numbers in between interpolator points
        self.assertAlmostEqual(gridpp.interpolate(0.5, x, y), 1)
        self.assertAlmostEqual(gridpp.interpolate(0.9, x, y), 1.8)
        self.assertAlmostEqual(gridpp.interpolate(1.5, x, y), 1.5)

        # Test outside range
        self.assertEqual(gridpp.interpolate(-1, x, y), 0)
        self.assertEqual(gridpp.interpolate(3, x, y), 1)

    def test_single(self):
        x = [0]
        y = [0]

        self.assertEqual(gridpp.interpolate(0, x, y), 0)
        self.assertEqual(gridpp.interpolate(-1, x, y), 0)
        self.assertEqual(gridpp.interpolate(1, x, y), 0)

    def test_empty(self):
        x = []
        y = []

        self.assertTrue(np.isnan(gridpp.interpolate(0, x, y)))
        self.assertTrue(np.isnan(gridpp.interpolate(-1, x, y)))
        self.assertTrue(np.isnan(gridpp.interpolate(1, x, y)))

    def test_nan(self):
        x = [0]
        y = [0]

        self.assertTrue(np.isnan(gridpp.interpolate(np.nan, x, y)))

    def test_duplicates_edge(self):
        """ Check multiple identical x-values on edge of curve """
        self.assertAlmostEqual(0.9, gridpp.interpolate(1, [0, 0, 0.5, 0.5, 1, 1], [0, 0.1, 0.4, 0.6, 0.9, 1]))
        self.assertAlmostEqual(0.1, gridpp.interpolate(0, [0, 0, 0.5, 0.5, 1, 1], [0, 0.1, 0.4, 0.6, 0.9, 1]))
        self.assertAlmostEqual(0.5, gridpp.interpolate(0.5, [0, 0, 0.5, 0.5, 1, 1], [0, 0.1, 0.4, 0.6, 0.9, 1]))

    def test_duplicates_middle(self):
        """ Check multiple identical x-values in the middle """
        self.assertAlmostEqual(0.4, gridpp.interpolate(0.499999, [0, 0, 0.5, 0.5, 1, 1], [0, 0.1, 0.4, 0.6, 0.9, 1]), 5)
        self.assertAlmostEqual(0.6, gridpp.interpolate(0.500001, [0, 0, 0.5, 0.5, 1, 1], [0, 0.1, 0.4, 0.6, 0.9, 1]), 5)

    def test_duplicates_all(self):
        """ Check if all x-values are identical """
        self.assertAlmostEqual(0.5, gridpp.interpolate(0, [0, 0, 0, 0, 0, 0], [0, 0.1, 0.4, 0.6, 0.9, 1]), 5)

if __name__ == '__main__':
    unittest.main()
