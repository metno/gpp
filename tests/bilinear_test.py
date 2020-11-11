from __future__ import print_function
import unittest
import gridpp
import numpy as np


lats = [60, 60, 60, 60, 60, 70]
lons = [10,10.1,10.2,10.3,10.4, 10]


class BilinearTest(unittest.TestCase):
    def test_neighbourhood(self):
        """Simple check
        50  6  7  8
        40  3  4  5
        30  0  1  2
            0 10 20
        """
        lons, lats = np.meshgrid([0, 10, 20], [30, 40, 50])
        grid = gridpp.Grid(lats, lons)
        values = np.zeros(lons.shape)
        values[:] = np.reshape(range(9), lons.shape)
        points = gridpp.Points([25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55], [-1, -1, -1, 0, 0, 0, 10, 10, 10, 20, 20, 20, 21, 21, 21])
        output = gridpp.nearest(grid, points, values)
        np.testing.assert_array_equal(output, (0, 3, 6, 0, 3, 6, 1, 4, 7, 2, 5, 8, 2, 5, 8))
        output = gridpp.bilinear(grid, points, values)
        np.testing.assert_array_equal(output, (0, 3, 6, 0, 3, 6, 1, 4, 7, 2, 5, 8, 2, 5, 8))

    def test_incompatible_sizes(self):
        """ Check for exception when input is different size than grid"""
        lats = [[0, 0, 0, 0], [1, 1, 1, 1]]
        lons = [[0, 1, 2, 3], [0, 1, 2, 3]]
        grid = gridpp.Grid(lats, lons)  # 2 x 4 grid
        points = gridpp.Points([0, 1], [0, 1])
        for N in [1, 3, 5]:
            values = np.zeros([N, N])
            with self.assertRaises(Exception) as e:
                gridpp.bilinear(grid, grid, values)
            with self.assertRaises(Exception) as e:
                gridpp.bilinear(grid, points, values)

    def test_rotation(self):
        """Check that the interpolation gives the same result when grid is rotated"""
        values = np.array([[0, 0, 0, 0], [0, 0, 1, 0], [0, 2, 3, 0], [0, 0, 0, 0]])
        rotations = range(0, 360, 10)
        lons0 = np.array([[-24, -12, 0, 12], [-24, -12, 0, 12], [-24, -12, 0, 12], [-24, -12, 0, 12]])
        lats0 = np.array([[-6, -6, -6, -6], [-3, -3, -3, -3], [3, 3, 3, 3], [6, 6, 6, 6]])
        lons_p0 = -6
        lats_p0 = 0
        lons = np.zeros([4, 4])
        lats = np.zeros([4, 4])
        for rotation in rotations:
            angle = rotation * 2 * 3.14159265 / 360
            lons_p = lons_p0 * np.cos(angle) - lats_p0 * np.sin(angle)
            lats_p = lons_p0 * np.sin(angle) + lats_p0 * np.cos(angle)
            points = gridpp.Points([lats_p], [lons_p])
            # print("Rotation: %d" % rotation)

            for i in range(lons0.shape[0]):
                for j in range(lons0.shape[1]):
                    lons[i, j] = lons0[i, j] * np.cos(angle) - lats0[i, j] * np.sin(angle)
                    lats[i, j] = lons0[i, j] * np.sin(angle) + lats0[i, j] * np.cos(angle)
            if 0:
                import matplotlib.pylab as mpl
                mpl.plot(lons, lats, '-ko')
                mpl.plot(lons_p, lats_p, '-ro')
                mpl.show()
            # mpl.plot(lons, lats)
            # mpl.show()
            # sys.exit()
            grid = gridpp.Grid(lats, lons)
            output = gridpp.bilinear(grid, points, values)[0]
            self.assertAlmostEqual(output, 1.5, 5)

    def test_parallelogram(self):
        """Skew a square"""
        values = np.array([[0, 1], [2, 3]])
        for skew in [-1,-0.1, 0, 0.1,1]:
            """Move the bottom edge from left to right"""
            lats = [[-1, -1], [1, 1]]
            lons = [[-1+skew, 1+skew], [-1-skew, 1-skew]]
            grid = gridpp.Grid(lats, lons)
            output = gridpp.bilinear(grid, gridpp.Points([0], [0]), values)[0]
            self.assertEqual(output, 1.5)

            """Move the left edge from bottom to top"""
            lats = [[-1+skew, -1-skew], [1+skew, 1-skew]]
            lons = [[-1, 1], [-1, 1]]
            grid = gridpp.Grid(lats, lons)
            output = gridpp.bilinear(grid, gridpp.Points([0], [0]), values)[0]
            self.assertEqual(output, 1.5)

    def test_non_parallelogram(self):
        """Non-parallelogram quadrilateral"""
        values = np.array([[0, 1], [2, 3]])
        lats = [[0, 0], [1, 2]]
        lons = [[0, 1], [0.1, 1]]
        grid = gridpp.Grid(lats, lons)
        output = gridpp.bilinear(grid, gridpp.Points([0.5], [0.5]), values)[0]
        # TODO: This value has not been confirmed by an independent calculation
        self.assertAlmostEqual(output, 1.1570623, 5)

    def test_with_missing_values(self):
        """Check that nearest neighbour is returned when a missing value is in the box"""
        values = np.array([[0, np.nan, 2, 3], [0, 1, 2, 3]])
        lats = [[0, 0, 0, 0], [1, 1, 1, 1]]
        lons = [[0, 1, 2, 3], [0, 1, 2, 3]]
        grid = gridpp.Grid(lats, lons)
        # The first three points have a missing value in its box, whereas the 4th does not
        output = gridpp.bilinear(grid, gridpp.Points([0, 0.3, 0.9, 0.5], [0, 0.3, 0.9, 2.5]), values)
        np.testing.assert_array_equal([0, 0, 1, 2.5], output)

    def test_grid_to_grid(self):
        """Check grid to grid itnerpolation"""
        values = np.reshape(np.arange(4), [2, 2])
        lons1, lats1 = np.meshgrid([0, 1], [0, 1])
        lons2, lats2 = np.meshgrid([0, 0.5, 1], [0, 0.5, 1])
        grid1 = gridpp.Grid(lats1, lons1)
        grid2 = gridpp.Grid(lats2, lons2)
        output = gridpp.bilinear(grid1, grid2, values)
        np.testing.assert_array_equal([[0, 0.5, 1], [1, 1.5, 2], [2, 2.5, 3]], output)

    def check(self):
        N = 4
        lons, lats = np.meshgrid(np.linspace(0, 1, N), np.linspace(0, 1, N))
        lons1, lats1 = np.meshgrid(np.linspace(0, 1, 101), np.linspace(0, 1, 101))
        grid = gridpp.Grid(lats, lons)
        grid1 = gridpp.Grid(lats1, lons1)
        values = np.random.rand(lons.shape[0], lons.shape[1])
        print(lons.shape, values.shape)
        values1 = gridpp.bilinear(grid, grid1, values)
        import matplotlib.pylab as mpl
        mpl.pcolormesh(lons1, lats1, values1)
        mpl.show()



if __name__ == '__main__':
    unittest.main()
