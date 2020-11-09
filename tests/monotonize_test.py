from __future__ import print_function
import unittest
import gridpp
import time
import numpy as np
import sys

class MonotonizeTest(unittest.TestCase):
    def test_empty_input(self):
        """Check for exception on empty curve"""
        with self.assertRaises(Exception) as e:
            gridpp.monotonize_curve([[], []])
        with self.assertRaises(Exception) as e:
            gridpp.monotonize_curve([[1, 2], []])
        with self.assertRaises(Exception) as e:
            gridpp.monotonize_curve([[], [1, 2]])

    def test_size_mismatch_input(self):
        """Check for exception on mismatch between sizes of x and y in curve"""
        with self.assertRaises(Exception) as e:
            gridpp.monotonize_curve([[1, 2, 3], [1, 2]])
        with self.assertRaises(Exception) as e:
            gridpp.monotonize_curve([[1, 2], [1, 2, 3]])

    def test_ok(self):
        x = [1, 2, 3]
        y = [1, 2, 3]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        np.testing.assert_array_equal(curve, curve_monotonic)

    def test_y_repeat(self):
        x = [1, 2, 3]
        y = [1, 1, 3]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        np.testing.assert_array_equal(curve, curve_monotonic)

    def test_x_repeat(self):
        x = [0, 1, 1, 3]
        y = [0, 1, 2, 3]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        np.testing.assert_array_equal([[0, 3], [0, 3]], curve_monotonic)

    def test_x_repeat_lower(self):
        """ Check that only the upper two points are kept """
        x = [0, 0, 1, 3]
        y = [0, 1, 2, 3]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        np.testing.assert_array_equal([[1, 3], [2, 3]], curve_monotonic)

    def test_x_repeat_upper(self):
        """ Check that only the lower two points are kept """
        x = [0, 1, 3, 3]
        y = [0, 1, 2, 3]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        np.testing.assert_array_equal([[0, 1], [0, 1]], curve_monotonic)

    def test_knot(self):
        """ """
        x = [0, 3, 2, 1, 5]
        y = [0, 1, 1, 2, 3]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        np.testing.assert_array_equal([[0, 5], [0, 3]], curve_monotonic)

    def test_upper_knot(self):
        """ Check that only the lower two points are kept """
        x = [0, 1, 3, 2]
        y = [0, 1, 2, 3]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        np.testing.assert_array_equal([[0, 1], [0, 1]], curve_monotonic)

    def test_x_double_repeat(self):
        x = [-1, 0, 2, 1, 3, 5, 4, 6]
        y = [i for i in range(len(x))]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        # print(curve, curve_monotonic)
        # np.testing.assert_array_equal([[0, 3, 6], [0, 3, 6]], curve_monotonic)

    def test_lower_knot(self):
        x = [-8, -9, -7, -6, -3, -1, 0, 1, 2, 3]
        y = [0, 0, 1, 2, 3, 5, 3, 6, 7, 9]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        np.testing.assert_array_equal([[-7, -6, -3, -1, 0, 1, 2, 3], [1, 2, 3, 5, 3, 6, 7, 9]], curve_monotonic)

    def test_two_knots_in_a_row(self):
        x = [0, 10, 20, 30, 25, 32, 31, 33]
        y = [0, 1, 2, 3, 4, 5, 6, 7]
        curve = [x, y]
        curve_monotonic = gridpp.monotonize_curve(curve)
        np.testing.assert_array_equal([[0, 10, 20, 33], [0, 1, 2, 7]], curve_monotonic)

    def test_with_missing(self):
        """Check curves with missing values"""
        q = np.nan
        m = -1
        x = [[0, 1, 2, 3, 4, 5, 6],
             # Missing first value
             [q, 0, 1, 2, 3, 4, 5, 6],
             [q, 0, 1, 2, 3, 4, 5, 6],
             [m, 0, 1, 2, 3, 4, 5, 6],
             # Two consecutive missing values
             [0, q, q, 1, 2, 3, 4, 5, 6],
             [0, q, q, 1, 2, 3, 4, 5, 6],
             [0, m, m, 1, 2, 3, 4, 5, 6],
             [0, m, q, 1, 2, 3, 4, 5, 6],
             # Two non-consecutive missing values
             [0, 1, q, 2, 3, q, 4, 5, 6],
             [0, 1, q, 2, 3, q, 4, 5, 6],
             [0, 1, m, 2, 3, m, 4, 5, 6],
             # Missing last value
             [0, 1, 2, 3, 4, 5, 6, q],
             [0, 1, 2, 3, 4, 5, 6, q],
             [0, 1, 2, 3, 4, 5, 6, m]]
        y = [[0, 1, 2, 3, 4, 5, 6],
             # Missing first value
             [m, 0, 1, 2, 3, 4, 5, 6],
             [q, 0, 1, 2, 3, 4, 5, 6],
             [q, 0, 1, 2, 3, 4, 5, 6],
             # Two consecutive missing values
             [0, m, m, 1, 2, 3, 4, 5, 6],
             [0, q, m, 1, 2, 3, 4, 5, 6],
             [0, q, q, 1, 2, 3, 4, 5, 6],
             [0, q, q, 1, 2, 3, 4, 5, 6],
             # Two non-consecutive missing values
             [0, 1, m, 2, 3, m, 4, 5, 6],
             [0, 1, q, 2, 3, m, 4, 5, 6],
             [0, 1, q, 2, 3, q, 4, 5, 6],
             # Missing last value
             [0, 1, 2, 3, 4, 5, 6, m],
             [0, 1, 2, 3, 4, 5, 6, q],
             [0, 1, 2, 3, 4, 5, 6, q]]
        for i in range(len(x)):
            x0 = x[i]
            y0 = y[i]
            curve = np.array([x0, y0], 'float32')
            curve_monotonic = gridpp.monotonize_curve(curve)
            np.testing.assert_array_equal([x[0], y[0]], curve_monotonic)


if __name__ == '__main__':
    unittest.main()
