from __future__ import print_function
import unittest
import gridpp
import numpy as np


lats = [60, 60, 60, 60, 60, 70]
lons = [10,10.1,10.2,10.3,10.4, 10]

"""Simple check
20 21 22 23 24
15 16 17 18 19
10 11 12 13 nan
5  6  7  nan  9
0  1  2  3  4
"""
values = np.reshape(range(25), [5, 5]).astype(float)
values[1, 3] = np.nan
values[2, 4] = np.nan
values = np.array(values)

class Test(unittest.TestCase):
    def test_invalid_arguments(self):
        field = np.ones([5, 5])
        halfwidth = -1
        quantiles = [-0.1, 1.1, np.nan, np.inf]

        for quantile in quantiles:
            with self.assertRaises(ValueError) as e:
                gridpp.neighbourhood_quantile(field, quantile, halfwidth)

    def test_empty_argument(self):
        halfwidth = 3
        for quantile in [0, 0.5, 1]:
            output = gridpp.neighbourhood_quantile([[]], quantile, halfwidth)
            self.assertEqual(len(output.shape), 2)
            self.assertEqual(output.shape[0], 0)
            self.assertEqual(output.shape[1], 0)

    def test_missing(self):
        """Checks that missing values are handled correctly"""
        empty = np.zeros([5, 5])
        empty[0:3, 0:3] = np.nan
        output = gridpp.neighbourhood_quantile(empty, 0.5, 1)
        self.assertTrue(np.isnan(np.array(output)[0:2,0:2]).all())

    def test_quantile(self):
        output = np.array(gridpp.neighbourhood_quantile(values, 0.5, 1))
        self.assertEqual(output[2][2], 12.5)
        self.assertEqual(output[2][3], 13)
        self.assertEqual(output[0][4], 4)


if __name__ == '__main__':
    unittest.main()
