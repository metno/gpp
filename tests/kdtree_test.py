from __future__ import print_function
import unittest
import gridpp
import numpy as np


class KDTreeTest(unittest.TestCase):
    def test(self):
        lats = [60, 61, 62]
        lons = [10, 10, 12]
        tree = gridpp.KDTree(lats, lons)
        np.testing.assert_array_equal(tree.get_neighbours(60, 10, 1), [0])
        np.testing.assert_array_equal(tree.get_neighbours(60, 10, 112000), [0, 1])

    def test_flat(self):
        tree = gridpp.KDTree([0, 1000, 2000], [0, 1000, 2000], True)
        q, dist = tree.get_neighbours_with_distance(100, 100, 1000)
        self.assertEqual(dist, 100 * np.sqrt(2))

    """ Check that duplicate points are found """
    def test_duplicate_points(self):
        tree = gridpp.KDTree([50, 50, 51], [0, 0, 10])
        I = tree.get_neighbours(50, 0.001, 1000)
        self.assertEqual(len(I), 2)
        self.assertTrue(0 in I)
        self.assertTrue(1 in I)

    """ Check that duplicate points that are equal to lookup point are found """
    def test_duplicate_points_identical(self):
        tree = gridpp.KDTree([50, 50, 51], [0, 0, 10])
        I = tree.get_neighbours(50, 0, 1000)
        self.assertEqual(len(I), 2)
        self.assertTrue(0 in I)
        self.assertTrue(1 in I)
    def test_pole(self):
        tree = gridpp.KDTree([89, 89, 90], [0, 180, 0])
        I, dist = tree.get_neighbours_with_distance(90, 0, 1000)
        self.assertEqual(len(I), 1)
        self.assertAlmostEqual(dist[0], 0)
        I, dist = tree.get_neighbours_with_distance(90, 10, 1000)
        self.assertEqual(len(I), 1)
        self.assertAlmostEqual(dist[0], 0)
        I, dist = tree.get_neighbours_with_distance(89.999, 0, 1000)
        self.assertEqual(len(I), 1)

    def test_pole_with_duplicate_points(self):
        tree = gridpp.KDTree([89, 89, 90, 90], [0, 180, 0, 10])
        I, dist = tree.get_neighbours_with_distance(90, 0, 1000)
        self.assertEqual(len(I), 2)
        self.assertTrue(2 in I)
        self.assertTrue(3 in I)
        self.assertAlmostEqual(dist[0], 0)
        self.assertAlmostEqual(dist[1], 0)

    def test_rad2deg(self):
        self.assertAlmostEqual(gridpp.KDTree_rad2deg(1), 180 / 3.14159265, 5)
        self.assertAlmostEqual(gridpp.KDTree_rad2deg(-1), -180 / 3.14159265, 5)
        self.assertAlmostEqual(gridpp.KDTree_rad2deg(0), 0, 5)

    def test_calc_distance(self):
        self.assertAlmostEqual(0, gridpp.KDTree_calc_distance(60,10,60,10));
        self.assertAlmostEqual(20037508, gridpp.KDTree_calc_distance(90,10,-90,10));
        self.assertAlmostEqual(20037508, gridpp.KDTree_calc_distance(0,0,0,180));
        self.assertAlmostEqual(16879114, gridpp.KDTree_calc_distance(60.5,5.25,-84.75,-101.75));
        self.assertAlmostEqual(124080.79, gridpp.KDTree_calc_distance(60,10,61,11), delta=0.1);

    def test_calc_distance_fast(self):
        self.assertEqual(0, gridpp.KDTree_calc_distance_fast(60,10,60,10));
        self.assertAlmostEqual(20037508, gridpp.KDTree_calc_distance_fast(90,10,-90,10), delta=100);
        self.assertAlmostEqual(20037508, gridpp.KDTree_calc_distance_fast(0,0,0,180), delta=100);
        # self.assertAlmostEqual(16879114, gridpp.KDTree_calc_distance_fast(60.5,5.25,-84.75,-101.75), delta=100);
        self.assertAlmostEqual(124080.79, gridpp.KDTree_calc_distance_fast(60,10,61,11), delta=100);

    def test_radius_match(self):
        """Check that points right on the radius edge count as a match"""
        points = gridpp.Points([0, 1000, 2000], [0, 0, 0], [0, 0, 0], [0, 0, 0], gridpp.Cartesian)
        I = points.get_neighbours(900, 0, 501)
        np.testing.assert_array_equal(I, [1])
        I = points.get_neighbours(900, 0, 99.99)
        np.testing.assert_array_equal(I, [])
        I = points.get_neighbours(0, 0, 1000)
        np.testing.assert_array_equal(I, [0])
        I = points.get_neighbours(0, 0, 1001)
        np.testing.assert_array_equal(I, [0, 1])
        I = points.get_neighbours(0, 0, 1001, False)
        np.testing.assert_array_equal(I, [1])


if __name__ == '__main__':
    unittest.main()
