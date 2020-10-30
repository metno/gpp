from __future__ import print_function
import unittest
import gridpp
import numpy as np


class QuantileMappingTest(unittest.TestCase):
    """ Test that the curve reproduces the reference values """
    def test_train(self):
        refs = [[1, 2, 3], [1, 1, 1]]
        fcst = [2, 3, 4]
        for ref in refs:
            curve = gridpp.quantile_mapping_curve(ref, fcst)
            np.testing.assert_array_equal(curve[0], fcst)
            np.testing.assert_array_equal(curve[1], ref)
            new_fcst = gridpp.apply_curve(fcst, curve, gridpp.OneToOne, gridpp.OneToOne)
            np.testing.assert_array_equal(ref, new_fcst)

    """ Test that the curve doesn't allow negative correlations """
    def test_negative(self):
        ref = [1, 0, -1]
        fcst = [2, 3, 4]
        curve = gridpp.quantile_mapping_curve(ref, fcst)
        np.testing.assert_array_equal(curve[0], fcst)
        np.testing.assert_array_equal(curve[1], np.sort(ref))
        new_fcst = gridpp.apply_curve(fcst, curve, gridpp.OneToOne, gridpp.OneToOne)
        np.testing.assert_array_equal(np.sort(ref), new_fcst)

    def test_quantiles(self):
        ref = np.arange(11)
        fcst = ref + 2
        quantiles = [0.1, 0.9]
        curve = gridpp.quantile_mapping_curve(ref, fcst, quantiles)
        np.testing.assert_array_equal([1, 9], curve[1])
        np.testing.assert_array_equal([3, 11], curve[0])

    def test_single_point(self):
        ref = [1]
        fcst = [2]
        curve = gridpp.quantile_mapping_curve(ref, fcst)
        np.testing.assert_array_almost_equal(curve, [fcst, ref])

    def test_dimension_mismatch(self):
        ref = [1, 2]
        fcst = [1, 2, 3]
        with self.assertRaises(Exception) as e:
            curve = gridpp.quantile_mapping_curve(ref, fcst)
        with self.assertRaises(Exception) as e:
            curve = gridpp.quantile_mapping_curve(ref, fcst, [0.1, 0.9])

    def test_empty_curve(self):
        ref = []
        fcst = []
        np.testing.assert_array_almost_equal([[],[]], gridpp.quantile_mapping_curve(ref, fcst))
        np.testing.assert_array_almost_equal([[],[]], gridpp.quantile_mapping_curve(ref, fcst, [0.1, 0.9]))

    def test_invalid_argument(self):
        ref = np.arange(11)
        fcst = ref + 2
        quantiles = [0.1, -1]
        with self.assertRaises(Exception) as e:
            curve = gridpp.quantile_mapping_curve(ref, fcst, quantiles)


if __name__ == '__main__':
    unittest.main()
