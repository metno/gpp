from __future__ import print_function
import unittest
import gridpp
import numpy as np


lons, lats = np.meshgrid([0, 1, 2], [0, 1, 2])
grid = gridpp.Grid(lats, lons)
points = gridpp.Points([0, 1], [0, 1])


class Test(unittest.TestCase):
    def test_invalid_radii(self):
        values = np.zeros([3, 3])
        radii = [-1, -1]
        value = 1
        for outside in [False, True]:
            with self.assertRaises(Exception) as e:
                gridpp.fill(grid, values, points, radii, value, outside)

    def test_invalid_number_of_radii(self):
        values = np.zeros([3, 3])
        for outside in [False, True]:
            with self.assertRaises(Exception) as e:
                gridpp.fill(grid, values, points, [1], value, outside)

    def test_dimension_mismatch(self):
        values = np.zeros([3, 2])
        radii = [1, 1]
        value = 1
        for outside in [False, True]:
            with self.assertRaises(Exception) as e:
                gridpp.fill(grid, values, points, radii, value, outside)

    def test_1(self):
        lons, lats = np.meshgrid([0, 1, 2, 3, 4], [0, 1, 2, 3, 4])
        grid = gridpp.Grid(lats * 1000, lons * 1000, np.zeros([5, 5]), np.zeros([5, 5]), gridpp.Cartesian)
        points = gridpp.Points([0, 0, 3000], [0, 3000, 3000], [0, 0, 0], [0, 0, 0], gridpp.Cartesian)
        values = np.zeros([5, 5])
        radii = [1010, 10, 2010]
        value = 1
        outside = False
        output = gridpp.fill(grid, values, points, radii, value, outside)
        np.testing.assert_array_almost_equal(output, [[1, 1, 0, 1, 0], [1, 0, 0, 1, 0],[0, 0, 1, 1, 1], [0, 1, 1, 1, 1], [0, 0, 1, 1, 1]])

        outside = True
        output = gridpp.fill(grid, values, points, radii, value, outside)
        np.testing.assert_array_almost_equal(output, [[0, 0, 1, 0, 1], [0, 1, 1, 0, 1],[1, 1, 0, 0, 0], [1, 0, 0, 0, 0], [1, 1, 0, 0, 0]])


if __name__ == '__main__':
    unittest.main()
