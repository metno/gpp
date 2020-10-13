from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_invalid_arguments(self):
        """Check that exception is thrown for invalid arguments"""
        field = np.ones([5, 5])
        nums = [-1, 0]
        for num in nums:
            with self.assertRaises(ValueError) as e:
                gridpp.get_neighbourhood_thresholds(field, num)

    def test_empty_argument(self):
        """Check that no thresholds are returned for an empty input"""
        for num in [1, 5]:
            np.testing.assert_array_equal(gridpp.get_neighbourhood_thresholds([], num), [])

    def test_high_num(self):
        """Check proper results when num > size of the input"""
        field = np.reshape(np.arange(4), [2, 2])
        nums = [4, 5, 6]
        for num in nums:
            q = gridpp.get_neighbourhood_thresholds(field, num)
            np.testing.assert_array_equal(q, [0, 1, 2, 3])

    def test_duplicates(self):
        """Check a real example with duplicates"""
        field = np.reshape([0, 0, 2, 3, 4, 5, 6, 11, 8, 9, 10, 11], [3, 4])
        nums = range(1, 5)
        for num in nums:
            q = gridpp.get_neighbourhood_thresholds(field, num)
            self.assertTrue((q >= 0).all())
            self.assertTrue((q <= 11).all())

    def test_3d(self):
        np.random.seed(1000)
        values = np.random.rand(10, 10)
        values3 = np.zeros([10, 10, 5])
        for i in range(5):
            values3[:, :, i] = values

        nums = [1, 5]
        for num in nums:
            output_2d = gridpp.get_neighbourhood_thresholds(values, num)
            output_3d = gridpp.get_neighbourhood_thresholds(values3, num)
            np.testing.assert_array_almost_equal(output_2d, output_3d)


if __name__ == '__main__':
    unittest.main()
