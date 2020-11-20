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
        """Check that exception is thrown for invalid arguments"""
        field = np.ones([5, 5])
        halfwidth = -1
        quantiles = [-0.1, 1.1, np.nan]
        thresholds = [0, 1]

        for quantile in quantiles:
            with self.assertRaises(ValueError) as e:
                gridpp.neighbourhood_quantile_fast(field, quantile, halfwidth, thresholds)

    def test_empty(self):
        for quantile in np.arange(0.1,0.9,0.1):
            for num_thresholds in [1, 2]:
                thresholds = gridpp.get_neighbourhood_thresholds(values, num_thresholds)
                output = gridpp.neighbourhood_quantile_fast([[]], 0.9, 1, thresholds)
                self.assertEqual(len(output.shape), 2)
                self.assertEqual(output.shape[0], 0)
                self.assertEqual(output.shape[1], 0)

    def test_single_threshold(self):
        """Checks what happens when a single threshold is provided"""
        thresholds = [0]
        field = np.reshape(np.arange(9), [3, 3])
        for halfwidth in [0, 1, 2]:
            output = gridpp.neighbourhood_quantile_fast(field, 0.9, halfwidth, thresholds)
            np.testing.assert_array_equal(output, np.zeros([3, 3]))

    def test_two_thresholds(self):
        """Checks what happens when a single threshold is provided"""
        thresholds = [0, 1]
        field = np.reshape(np.arange(9), [3, 3])
        for halfwidth in [0, 1, 2]:
            output = gridpp.neighbourhood_quantile_fast(field, 0.9, 0, thresholds)
            self.assertTrue(((output >= 0) & (output <= 1)).all())

    def test_missing(self):
        empty = np.zeros([5, 5])
        empty[0:3, 0:3] = np.nan
        thresholds = [0, 1]
        output = gridpp.neighbourhood_quantile_fast(empty, 0.5, 1, thresholds)
        self.assertTrue(np.isnan(np.array(output)[0:2,0:2]).all())

    def test_quantile(self):
        thresholds = gridpp.get_neighbourhood_thresholds(values, 100)
        output = np.array(gridpp.neighbourhood_quantile_fast(values, 0.5, 1, thresholds))
        self.assertEqual(output[2][2], 12)   # Should be 12.5
        self.assertEqual(output[2][3], 12.5) # Should be 13

        output = np.array(gridpp.neighbourhood_quantile_fast(np.full([100,100], np.nan), 0.5, 1, thresholds))
        self.assertTrue(np.isnan(np.array(output)).all())

        output = np.array(gridpp.neighbourhood_quantile_fast(np.zeros([100,100]), 0.5, 1, thresholds))
        self.assertTrue((np.array(output) == 0).all())

        output = np.array(gridpp.neighbourhood_quantile(values, 0.5, 1))
        self.assertEqual(output[2][2], 12.5)
        self.assertEqual(output[2][3], 13)
        self.assertEqual(output[0][4], 4)

    def test_3d(self):
        np.random.seed(1000)
        values = np.random.rand(200, 200)
        values3 = np.zeros([200, 200, 5])
        for i in range(5):
            values3[:, :, i] = values
        halfwidths = [0, 1, 5]
        quantile = 0.5
        thresholds = [0, 0.25, 0.5, 0.75, 1]
        for halfwidth in halfwidths:
            output_2d = gridpp.neighbourhood_quantile_fast(values, quantile, halfwidth, thresholds)
            output_3d = gridpp.neighbourhood_quantile_fast(values3, quantile, halfwidth, thresholds)
            np.testing.assert_array_almost_equal(output_2d, output_3d)

    def test_varying_quantile(self):
        """ For now check that this runs """
        values = np.array([[0, 1], [2, 3], [4, 5]])
        halfwidth = 1
        quantiles = np.ones(values.shape) * 0.5
        thresholds = [0, 0.25, 0.5, 0.75, 1]
        gridpp.neighbourhood_quantile_fast(values, quantiles, halfwidth, thresholds)

        values = np.nan *np.zeros(values.shape)
        np.testing.assert_array_equal(values, gridpp.neighbourhood_quantile_fast(values, quantiles, halfwidth, thresholds))


if __name__ == '__main__':
    unittest.main()
