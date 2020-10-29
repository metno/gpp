from __future__ import print_function
import unittest
import gridpp
import numpy as np
import psutil
import os


def memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss


class Test(unittest.TestCase):
    def setUp(self):
        self.xs = [0, -1, 1, 0, 1]
        self.ys = [0, -1, 1, 1, 0]
        self.speeds = [0, np.sqrt(2), np.sqrt(2), 1, 1]
        self.directions = [180, 45, 225, 180, 270]

    def test_speed(self):
        for i in range(len(self.xs)):
            self.assertAlmostEqual(self.speeds[i], gridpp.wind_speed(self.xs[i], self.ys[i]))
        np.testing.assert_array_almost_equal(self.speeds, gridpp.wind_speed(self.xs, self.ys))

    def test_empty_input(self):
        np.testing.assert_array_almost_equal(gridpp.wind_speed([], []), [])
        np.testing.assert_array_almost_equal(gridpp.wind_direction([], []), [])

    def test_direction(self):
        for i in range(len(self.xs)):
            self.assertAlmostEqual(self.directions[i], gridpp.wind_direction(self.xs[i], self.ys[i]))
        np.testing.assert_array_almost_equal(self.directions, gridpp.wind_direction(self.xs, self.ys))

    def test_missing(self):
        """Check that if one or more values are missing, the result is NaN"""
        for func in [gridpp.wind_speed, gridpp.wind_direction]:
            with self.subTest(func=func):
                self.assertTrue(np.isnan(func(0, np.nan)))
                self.assertTrue(np.isnan(func(np.nan, 0)))
                self.assertTrue(np.isnan(func(np.nan, np.nan)))
                np.testing.assert_array_almost_equal(func([0, np.nan, np.nan], [np.nan, 0, np.nan]), [np.nan, np.nan, np.nan])

    def test_dimension_mismatch(self):
        for func in [gridpp.wind_speed, gridpp.wind_direction]:
            with self.subTest(func=func):
                with self.assertRaises(Exception) as e:
                    func([0], [0, 1])
                with self.assertRaises(Exception) as e:
                    func([], [0, 1])
                with self.assertRaises(Exception) as e:
                    func([0], [])


if __name__ == '__main__':
    unittest.main()
