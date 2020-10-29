from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_basic(self):
        x = [0, 1000, 2000, 3000, np.nan]
        structure_corr = dict()
        barnes = gridpp.BarnesStructure(2000)
        structure_corr[barnes] = [1, 0.8824968934059143, 0.6065306663513184, 0.32465246319770813, 0]
        structure_corr[gridpp.CressmanStructure(2000)] = [1, 0.6, 0, 0, 0]
        structure_corr[gridpp.CrossValidation(barnes, 1000)] = [0, 0, 0.6065306663513184, 0.32465246319770813, 0]
        N = len(x)
        for structure, corr in structure_corr.items():
            is_cv = isinstance(structure, gridpp.CrossValidation)
            for i in range(N):
                with self.subTest(structure=type(structure), i=i):
                    p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
                    p2 = gridpp.Point(x[i], 0, 0, 0, gridpp.Cartesian)
                    funcs = [structure.corr, structure.corr_background]
                    if is_cv:
                        funcs = [structure.corr_background]
                    for func in funcs:
                        self.assertAlmostEqual(corr[i], func(p1, p2))

                        # Check that changing the order does not change the results
                        self.assertAlmostEqual(corr[i], func(p2, p1))

                        # Check identical points
                        if not is_cv and not np.isnan(x[i]):
                            self.assertAlmostEqual(1, func(p2, p2))

    def test_invalid_h(self):
        structures = list()
        structures = [gridpp.BarnesStructure, gridpp.CressmanStructure]
        for structure in structures:
            for h in [-1, np.nan]:
                with self.assertRaises(Exception) as e:
                    s = structure(h)
                with self.assertRaises(Exception) as e:
                    s = structure(h, 100)

    def test_barnes_h(self):
        hs = [-1, np.nan]
        for h in hs:
            with self.subTest(h=h):
                with self.assertRaises(Exception) as e:
                    structure = gridpp.BarnesStructure(h)

    def test_barnes_hmax(self):
        hmaxs = [0, 1000, 2000, 10000]
        p0 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        dist_ans = {0:1, 1000:0.8824968934059143, 2000:0.6065306663513184, 3000:0.32465246319770813}
        for hmax in hmaxs:
            for dist, ans in dist_ans.items():
                with self.subTest(hmax=hmax, dist=dist):
                    structure = gridpp.BarnesStructure(2000, 0, 0, hmax)
                    corr = structure.corr(p0, gridpp.Point(dist, 0, 0, 0, gridpp.Cartesian))
                    if dist > hmax:
                        self.assertEqual(0, corr)
                    else:
                        self.assertEqual(ans, corr)

    def test_invalid_elevation(self):
        """Check that point elevations are ignored if one is missing"""
        structures = [gridpp.BarnesStructure, gridpp.CressmanStructure]
        h = 2000
        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(1000, 0, 0, 0, gridpp.Cartesian)
        p3 = gridpp.Point(1000, 0, float('nan'), 0, gridpp.Cartesian)
        for structure in structures:
            with self.subTest(structure=structure):
                s1 = structure(h, 0)
                s2 = structure(h, 100)
                self.assertAlmostEqual(s1.corr(p1, p3), s1.corr(p1, p2))
                self.assertAlmostEqual(s2.corr(p1, p3), s2.corr(p1, p2))

    def test_invalid_v(self):
        structures = list()
        structures = [gridpp.BarnesStructure, gridpp.CressmanStructure]
        h = 2000
        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(1000, 0, 0, 0, gridpp.Cartesian)
        for structure in structures:
            for v in [-1, np.nan]:
                with self.assertRaises(Exception) as e:
                    s = structure(h, v)

    def test_invalid_w(self):
        structures = [gridpp.BarnesStructure, gridpp.CressmanStructure]
        h = 2000
        v = 100
        for structure in structures:
            for w in [-1, np.nan]:
                with self.assertRaises(Exception) as e:
                    s = structure(h, v, w)

    def test_invalid_hmax(self):
        with self.assertRaises(Exception) as e:
            structure = gridpp.BarnesStructure(2000, 100, 0, -1)

    def test_invalid_cv(self):
        barnes = gridpp.BarnesStructure(2000)
        for dist in [-1, np.nan]:
            with self.assertRaises(Exception) as e:
                structure = gridpp.CrossValidation(barnes, dist)


if __name__ == '__main__':
    unittest.main()
