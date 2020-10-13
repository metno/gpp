from __future__ import print_function
import unittest
import gridpp
import numpy as np


lats = [60, 60, 60, 60, 60, 70]
lons = [10,10.1,10.2,10.3,10.4, 10]


class Test(unittest.TestCase):
    def test_grid_to_point(self):
        """Check that grid to point interpolation works
        50  6  7  8
        40  3  4  5
        30  0  1  2
            0 10 20
        """
        for num_treads in [1, 2]:
            gridpp.set_omp_threads(num_treads)
            lons, lats = np.meshgrid([0, 10, 20], [30, 40, 50])
            grid = gridpp.Grid(lats, lons)
            values = np.reshape(range(9), lons.shape)
            points = gridpp.Points([25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55], [-1, -1, -1, 0, 0, 0, 10, 10, 10, 20, 20, 20, 21, 21, 21])
            output = gridpp.nearest(grid, points, values)
            np.testing.assert_array_equal(output, (0, 3, 6, 0, 3, 6, 1, 4, 7, 2, 5, 8, 2, 5, 8))
            output = gridpp.bilinear(grid, points, values)
            np.testing.assert_array_equal(output, (0, 3, 6, 0, 3, 6, 1, 4, 7, 2, 5, 8, 2, 5, 8))

    def test_dimension_mismatch(self):
        lons, lats = np.meshgrid([0, 10, 20], [30, 40, 50])
        grid = gridpp.Grid(lats, lons)
        points = gridpp.Points([0, 1], [0, 1])
        values = np.zeros([3, 2])
        with self.assertRaises(Exception) as e:
            gridpp.nearest(grid, grid, values)

        with self.assertRaises(Exception) as e:
            gridpp.nearest(grid, points, values)

        values3 = np.zeros([3, 2, 3])
        with self.assertRaises(Exception) as e:
            gridpp.nearest(grid, grid, values3)

        with self.assertRaises(Exception) as e:
            gridpp.nearest(grid, points, values3)

    def test_one_row(self):
        gridpp.set_omp_threads(1)
        lons, lats = np.meshgrid([0], [30, 40, 50])
        grid = gridpp.Grid(lats, lons)
        values = np.zeros(lons.shape)
        values[:] = np.reshape(range(3), lons.shape)
        points = gridpp.Points([25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55], [-1, -1, -1, 0, 0, 0, 10, 10, 10, 20, 20, 20, 21, 21, 21])
        output = gridpp.nearest(grid, points, values)
        np.testing.assert_array_equal(output, (0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2))

    def test_grid_to_grid(self):
        """Check that grid to grid interpolation works"""
        lons1, lats1 = np.meshgrid([0, 10, 20], [30, 40, 50])
        lons2, lats2 = np.meshgrid([0, 20], [30, 50])
        grid1 = gridpp.Grid(lats1, lons1)
        grid2 = gridpp.Grid(lats2, lons2)
        values = np.reshape(range(9), lons1.shape)
        output = gridpp.nearest(grid1, grid2, values)
        np.testing.assert_array_equal(output, [[0, 2], [6, 8]])

    def test_grid_to_grid_3d(self):
        """Check that grid to grid interpolation for 3D fields works"""
        lons1, lats1 = np.meshgrid([0, 10, 20], [30, 40, 50])
        lons2, lats2 = np.meshgrid([0, 20], [30, 50])
        grid1 = gridpp.Grid(lats1, lons1)
        grid2 = gridpp.Grid(lats2, lons2)
        values = np.reshape(range(18), [2, lons1.shape[0], lons1.shape[1]])
        output = gridpp.nearest(grid1, grid2, values)
        np.testing.assert_array_equal(output, [[[0, 2], [6, 8]], [[9, 11], [15, 17]]])

    def test_grid_to_points_3d(self):
        """Check that grid to point interpolation for 3D fields works"""
        lons, lats = np.meshgrid([0, 10, 20], [0, 10, 20])
        grid = gridpp.Grid(lats, lons)
        points = gridpp.Points([-1, 6], [-1, 6])
        values = np.reshape(range(18), [2, lons.shape[0], lons.shape[1]])
        output = gridpp.nearest(grid, points, values)
        np.testing.assert_array_equal(output, [[0, 4], [9, 13]])


if __name__ == '__main__':
    unittest.main()
