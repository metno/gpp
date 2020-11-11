from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_size(self):
        self.assertEqual(2, gridpp.Points([0,1], [0,1]).size())

    def test_attributes(self):
        """Test that lat, lons, etc are set"""
        points = gridpp.Points([0, 1], [2, 3], [4, 5], [0.1, 0.2])
        np.testing.assert_array_almost_equal(points.get_lats(), [0, 1])
        np.testing.assert_array_almost_equal(points.get_lons(), [2, 3])
        np.testing.assert_array_almost_equal(points.get_elevs(), [4, 5])
        np.testing.assert_array_almost_equal(points.get_lafs(), [0.1, 0.2])

    def test_get_neighbours(self):
        points = gridpp.Points([0, 1000, 2000], [0, 1000, 2000], [0,0,0], [0,0,0], gridpp.Cartesian)
        np.testing.assert_array_almost_equal(np.sort(points.get_neighbours(0, 0, 1500)), [0, 1])
        np.testing.assert_array_almost_equal(np.sort(points.get_neighbours(900, 900, 1600)), [0, 1, 2])
        np.testing.assert_array_almost_equal(np.sort(points.get_neighbours(-100, -100, 100)), [])

    def test_get_nearest_neighbour(self):
        points = gridpp.Points([0, 1000, 2000], [0, 1000, 2000], [0,0,0], [0,0,0], gridpp.Cartesian)
        self.assertEqual(points.get_nearest_neighbour(-100, -100), 0)
        self.assertEqual(points.get_nearest_neighbour(0, 0), 0)
        self.assertEqual(points.get_nearest_neighbour(900, 900), 1)
        self.assertEqual(points.get_nearest_neighbour(2100, 2100), 2)

    def test_get_neighbours_with_distance(self):
        points = gridpp.Points([0, 1000, 2000], [0, 1000, 2000], [0,0,0], [0,0,0], gridpp.Cartesian)
        indices, distances = points.get_neighbours_with_distance(0, 0, 1500)
        self.assertEqual(len(indices), 2)
        np.testing.assert_array_almost_equal(distances, [0, np.sqrt(2)*1000], 4)

    def test_empty(self):
        points = gridpp.Points()
        self.assertEqual(points.get_nearest_neighbour(0, 0), -1)
        np.testing.assert_array_equal(points.get_closest_neighbours(0, 0, 5), [])
        np.testing.assert_array_equal(points.get_neighbours(0, 0, 1000), [])
        self.assertEqual(points.get_num_neighbours(0, 0, 1000), 0)
        np.testing.assert_array_equal(points.get_lats(), [])
        np.testing.assert_array_equal(points.get_lons(), [])
        np.testing.assert_array_equal(points.get_elevs(), [])
        np.testing.assert_array_equal(points.get_lafs(), [])
        self.assertEqual(points.size(), 0)

    def test_get_in_domain_indices(self):
        lons, lats = np.meshgrid([0, 1, 2], [0, 1, 2])
        grid = gridpp.Grid(lats, lons)
        # 4 corners, 4 edges, inside, outside each edge
        lats = np.array([0, 0, 2, 2, 0, 1, 1, 2, 1, -1, 1, 1, 3])
        lons = np.array([0, 2, 0, 2, 1, 0, 2, 1, 1, 1, -1, 3, 1])
        N = len(lats)
        points = gridpp.Points(lats, lons)
        indices = points.get_in_domain_indices(grid)
        Iinside = [0, 1, 2, 3, 4, 5, 6, 7, 8]
        np.testing.assert_array_equal(indices, Iinside)

        points1 = points.get_in_domain(grid)
        np.testing.assert_array_equal(np.sort(points1.get_lats()), np.sort(lats[Iinside]))
        np.testing.assert_array_equal(np.sort(points1.get_lons()), np.sort(lons[Iinside]))

    def test_get_in_domain_indices_empty_points(self):
        lons, lats = np.meshgrid([0, 1, 2], [0, 1, 2])
        grid = gridpp.Grid(lats, lons)
        points_empty = gridpp.Points([], [])
        indices = points_empty.get_in_domain_indices(grid)
        np.testing.assert_array_equal(indices, [])

    def test_get_in_domain_indices_empty_grid(self):
        grid = gridpp.Grid()
        points_empty = gridpp.Points([0], [0])
        indices = points_empty.get_in_domain_indices(grid)
        np.testing.assert_array_equal(indices, [])

    def test_copy_constructor(self):
        lats = [0, 1, 2]
        lons = [3, 4, 5]
        elevs = [6, 7, 8]
        lafs = [0.1, 0.2, 0.3]
        points = gridpp.Points(lats, lons, elevs, lafs)
        points2 = points
        for p in [points, points2]:
            np.testing.assert_array_almost_equal(p.get_lats(), lats)
            np.testing.assert_array_almost_equal(p.get_lons(), lons)
            np.testing.assert_array_almost_equal(p.get_elevs(), elevs)
            np.testing.assert_array_almost_equal(p.get_lafs(), lafs)

    def test_get_nearest_neighbour_no_match(self):
        """Check that an exact match is removed"""
        points = gridpp.Points([0, 1, 2, 2, 4], [0]*5)
        self.assertAlmostEqual(points.get_nearest_neighbour(0, 0, False), 1)
        self.assertAlmostEqual(points.get_nearest_neighbour(0, 0, True), 0)
        # Check that multiple of identical matches are removed
        self.assertEqual(points.get_nearest_neighbour(2, 0, False), 1)
        self.assertTrue(points.get_nearest_neighbour(2, 0, True) in [2, 3])


if __name__ == '__main__':
    unittest.main()
