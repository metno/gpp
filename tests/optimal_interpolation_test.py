from __future__ import print_function
import unittest
import gridpp
import numpy as np
import collections


class Test(unittest.TestCase):
    def test_invalid_arguments(self):
        """ Check that exception is thrown on invalid input values """

        # Set up struct with valid input arguments
        ok_args = collections.OrderedDict({
                'grid' : gridpp.Grid([[0,0,0]], [[0,2500,10000]], [[0,0,0]], [[0,0,0]], gridpp.Cartesian),
                'background' : np.zeros([1, 3]),
                'points' : gridpp.Points([0], [2500], [0], [0], gridpp.Cartesian),
                'pobs' : [1],
                'pratios' : [0.1],
                'pbackground' : [0],
                'structure' : gridpp.BarnesStructure(2500),
                'max_points' : 10
        })

        # Set up struct with invalid input arguments that will be substituted into ok_args one at a
        # time in order to look for an exception being raised. Use an array of different invalid
        # arguments for each key.
        x = np.zeros([3,2])
        invalid_args = {
                # Grid size mismatch, and coordinate-type mismatch
                'grid' : [gridpp.Grid(x, x, x, x, gridpp.Cartesian), gridpp.Grid([[0,0,0]], [[0,2500,10000]])],
                # Points size mismatch, and coordinate-type mismatch
                'points' : [gridpp.Points([0, 1], [0, 2500], [0, 0], [0, 0], gridpp.Cartesian), gridpp.Points([0], [2500])],
                'pratios' : [np.zeros(11)],
                'pobs' : [np.zeros([11])],
                'background' : [np.zeros([2, 11])],
                'pbackground' : [np.zeros(21)],
                'max_points' : [-1]
        }

        for key in invalid_args.keys():
            for arg in invalid_args[key]:
                args0 = ok_args.copy()
                args0[key] = arg
                q = [args0[f] for f in args0]
                with self.subTest(key=key, arg=arg):
                    with self.assertRaises(ValueError) as e:
                        output = gridpp.optimal_interpolation(*q)


    def test_simple_1d(self):
        N = 3
        y = [[0, 0, 0]]
        x = [[0, 2500, 10000]]
        grid = gridpp.Grid(y, x, y, y, gridpp.Cartesian)
        points = gridpp.Points([0], [2500], [0], [0], gridpp.Cartesian)
        pratios = [0.1]
        structure = gridpp.BarnesStructure(2500)
        pobs = [1]
        background = np.zeros([1, N])
        pbackground = [0]
        max_points = 10
        output = gridpp.optimal_interpolation(grid, background, points, pobs, pratios, pbackground, structure, max_points)
        np.testing.assert_array_almost_equal(output, np.array([[np.exp(-0.5)/1.1, 1/1.1, np.exp(-0.5*9)/1.1]]))

    def test_simple_1d_full(self):
        N = 3
        y = [[0, 0, 0]]
        x = [[0, 2500, 10000]]
        grid = gridpp.Grid(y, x, y, y, gridpp.Cartesian)
        points = gridpp.Points([0], [2500], [0], [0], gridpp.Cartesian)
        bvariance = np.ones([1, N])
        obs_variance = [0.1]
        bvariance_at_points = [1]
        structure = gridpp.BarnesStructure(2500)
        pobs = [1]
        background = np.zeros([1, N])
        background_at_points = [0]
        max_points = 10
        output, sigma = gridpp.optimal_interpolation_full(grid, background, bvariance, points, pobs,
                obs_variance, background_at_points, bvariance_at_points,
                structure, max_points)
        # np.testing.assert_array_almost_equal(output, np.array([[np.exp(-0.5)/1.1, 1/1.1, np.exp(-0.5*9)/1.1]]))
        # np.testing.assert_array_almost_equal(sigma, np.array([[0, np.sqrt(0.1/1.1), 1]]))
        self.assertAlmostEqual(sigma[0, 1], 0.1/1.1)

    def test_cross_validation(self):
        y = np.array([0, 1000, 2000, 3000])
        N = len(y)
        obs = np.array([0, 1, 2, 3])
        background = np.zeros(N)
        points = gridpp.Points(y, np.zeros(N), np.zeros(N), np.zeros(N), gridpp.Cartesian)
        ratios = np.ones(N)
        Icv = [0, 2, 3]
        points_cv = gridpp.Points(y[Icv], np.zeros(N-1), np.zeros(N-1), np.zeros(N-1), gridpp.Cartesian)
        structure = gridpp.BarnesStructure(1000, 0)
        structure_cv = gridpp.CrossValidation(structure, 750)

        analysis = gridpp.optimal_interpolation(points, background, points_cv, obs[Icv], ratios[Icv], background[Icv], structure, 100)
        analysis_cv = gridpp.optimal_interpolation(points, background, points, obs, ratios, background, structure_cv, 100)
        # print(analysis, analysis_cv)

    def test_cross_validation_grid(self):
        """ Check that the CV structure function works """
        np.random.seed(1000)
        y, x = np.meshgrid(np.arange(0, 3500, 500), np.arange(0, 3500, 500))
        Y = y.shape[0]
        X = y.shape[1]
        grid = gridpp.Grid(y, x, np.zeros(x.shape), np.zeros(x.shape), gridpp.Cartesian)
        background = np.random.rand(Y, X) * 0

        obs = np.array([10, 20, 30])
        x_o = np.array([1000, 2000, 3000])
        y_o = np.array([1000, 2000, 3000])
        N = len(obs)
        points = gridpp.Points(y_o, x_o, np.zeros(N), np.zeros(N), gridpp.Cartesian)

        background_o = gridpp.nearest(grid, points, background)

        ratios = np.ones(N)
        k = 0
        ii = np.arange(N).astype('int') != k
        points_cv = gridpp.Points(y_o[ii], x_o[ii], np.zeros(N-1), np.zeros(N-1), gridpp.Cartesian)
        structure = gridpp.BarnesStructure(1000, 0)
        structure_cv = gridpp.CrossValidation(structure, 750)

        analysis = gridpp.optimal_interpolation(grid, background, points_cv, obs[ii], ratios[ii],
                background_o[ii], structure, 100)
        analysis_cv = gridpp.optimal_interpolation(points, background_o, points, obs, ratios,
                background_o, structure_cv, 100)

        self.assertAlmostEqual(gridpp.nearest(grid, points, analysis)[k], analysis_cv[k])

    def test_missing_values(self):
        """Check that missing values are not used in OI"""
        obs = np.array([1, np.nan, 2, 3, np.nan, np.nan, 4, np.nan])
        N = len(obs)
        y = np.arange(0, N*1000, 1000)
        background = np.zeros(N)
        points = gridpp.Points(y, np.zeros(N), np.zeros(N), np.zeros(N), gridpp.Cartesian)
        ratios = np.ones(N)
        structure = gridpp.BarnesStructure(1000, 0)
        analysis = gridpp.optimal_interpolation(points, background, points, obs, ratios, background, structure, 100)

        I = np.where(np.isnan(y) == 0)[0]
        points1 = gridpp.Points(y[I], np.zeros(len(I)), np.zeros(len(I)), np.zeros(len(I)), gridpp.Cartesian)
        analysis1 = gridpp.optimal_interpolation(points, background, points1, obs[I], ratios[I], background[I], structure, 100)
        np.testing.assert_array_almost_equal(analysis, analysis1)


if __name__ == '__main__':
    unittest.main()
