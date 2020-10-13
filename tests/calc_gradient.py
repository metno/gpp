from __future__ import print_function
import unittest
import gridpp
import numpy as np

lats, lons = np.meshgrid([0, 1, 2], [0, 1, 2])
elevs = np.zeros([3, 3])
elevs[0, 0] = -10
elevs[1, 1] = 10
points = gridpp.Points([-1, 0.9], [-1, 0.9], [-5, 5])
grid = gridpp.Grid(lats, lons, elevs)


class Test(unittest.TestCase):
    def test_dimension_mismatch(self):
        gradient = 0
        values = np.zeros([3, 2])
        with self.assertRaises(Exception) as e:
            output = gridpp.simple_gradient(grid, points, values, gradient)
        with self.assertRaises(Exception) as e:
            output = gridpp.simple_gradient(grid, grid, values, gradient)

    def test_missing_gradient(self):
        values = np.zeros([3, 3])
        gradient = np.nan
        output = gridpp.simple_gradient(grid, points, values, gradient)
        np.testing.assert_array_almost_equal(output, [np.nan, np.nan])

    def test_infinite_gradient(self):
        values = np.zeros([3, 3])
        gradient = np.inf
        output = gridpp.simple_gradient(grid, points, values, gradient)
        np.testing.assert_array_almost_equal(output, [np.inf, -np.inf])

    def test_missing_values(self):
        values = np.zeros([3, 3])
        values[1, 1] = np.nan
        gradient = 1
        output = gridpp.simple_gradient(grid, points, values, gradient)
        np.testing.assert_array_almost_equal(output, [5, np.nan])

    def test_no_grid_elev(self):
        """Check that output is missing when grid does not have elevation defined"""
        grid0 = gridpp.Grid(lats, lons)
        values = np.reshape(np.arange(9), [3, 3])
        for gradient in [0, 1]:
            output = gridpp.simple_gradient(grid0, points, values, gradient)
            np.testing.assert_array_almost_equal(output, [np.nan, np.nan])

    def test_no_point_elev(self):
        """Check that output is missing when points do not have elevation defined"""
        points0 = gridpp.Points([-1, 0.9], [-1, 0.9])
        values = np.reshape(np.arange(9), [3, 3])
        for gradient in [0, 1]:
            output = gridpp.simple_gradient(grid, points0, values, gradient)
            np.testing.assert_array_almost_equal(output, [np.nan, np.nan])

    def test_1(self):
        values = np.zeros([3, 3])
        values[0, 0] = 4
        values[1, 1] = 3
        # Point 0: Move up 5 meters
        # Point 1: Move down 5 meters
        gradient = 0
        output = gridpp.simple_gradient(grid, points, values, gradient)
        np.testing.assert_array_almost_equal(output, [4, 3])

        gradient = 1
        output = gridpp.simple_gradient(grid, points, values, gradient)
        np.testing.assert_array_almost_equal(output, [9, -2])

    def test_grid_to_grid(self):

        # Output grid has one row and column outside domain
        values = np.zeros([3, 3])
        values[0, 0] = 4
        values[1, 1] = 3
        lats0, lons0 = np.meshgrid([-0.1, 0.1, 1.1], [-0.1, 0.1, 1.1])
        elevs0 = np.reshape(np.arange(9), [3, 3])
        grid0 = gridpp.Grid(lats0, lons0, elevs0)
        gradient = 0
        output = gridpp.simple_gradient(grid, grid0, values, gradient)
        np.testing.assert_array_almost_equal(output, [[4, 4, 0], [4, 4, 0], [0, 0, 3]])

        gradient = 1
        output = gridpp.simple_gradient(grid, grid0, values, gradient)
        # Elevation changes:
        # -10 -> 0, -10 -> 1,  0 -> 2
        # -10 -> 3, -10 -> 4,  0 -> 5
        #   0 -> 6,   0 -> 7, 10 -> 8
        np.testing.assert_array_almost_equal(output, [[14, 15, 2], [17, 18, 5], [6, 7, 1]])



if __name__ == '__main__':
    unittest.main()
