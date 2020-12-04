from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_get_statistic(self):
        self.assertEqual(gridpp.get_statistic("mean"), gridpp.Mean)
        self.assertEqual(gridpp.get_statistic("min"), gridpp.Min)
        self.assertEqual(gridpp.get_statistic("max"), gridpp.Max)
        self.assertEqual(gridpp.get_statistic("median"), gridpp.Median)
        self.assertEqual(gridpp.get_statistic("quantile"), gridpp.Quantile)
        self.assertEqual(gridpp.get_statistic("std"), gridpp.Std)
        self.assertEqual(gridpp.get_statistic("sum"), gridpp.Sum)

    def test_unknown_statistic(self):
        self.assertEqual(gridpp.get_statistic("mean1"), gridpp.Unknown)

    """ Check that it doesn't cause any errors """
    def test_version(self):
        gridpp.version()

    def test_clock(self):
        time = gridpp.clock()
        self.assertTrue(time > 0)

    def test_is_valid(self):
        self.assertTrue(gridpp.is_valid(1))
        self.assertTrue(gridpp.is_valid(-1))
        self.assertTrue(gridpp.is_valid(-999))  # Check that the old missing value indicator is valid now
        self.assertFalse(gridpp.is_valid(np.nan))

    def test_calc_statistic_mean(self):
        self.assertEqual(gridpp.calc_statistic([0, 1, 2], gridpp.Mean), 1)
        self.assertEqual(gridpp.calc_statistic([0, 1, np.nan], gridpp.Mean), 0.5)
        self.assertEqual(gridpp.calc_statistic([np.nan, 1, np.nan], gridpp.Mean), 1)
        self.assertTrue(np.isnan(gridpp.calc_statistic([np.nan, np.nan, np.nan], gridpp.Mean)))
        self.assertTrue(np.isnan(gridpp.calc_statistic([], gridpp.Mean)))

    def test_calc_statistic_sum(self):
        self.assertEqual(gridpp.calc_statistic([0, 1, 2], gridpp.Sum), 3)
        self.assertEqual(gridpp.calc_statistic([0, 1, np.nan], gridpp.Sum), 1)
        self.assertEqual(gridpp.calc_statistic([np.nan, 1, np.nan], gridpp.Sum), 1)
        self.assertTrue(np.isnan(gridpp.calc_statistic([np.nan, np.nan, np.nan], gridpp.Sum)))
        self.assertTrue(np.isnan(gridpp.calc_statistic([], gridpp.Sum)))

    def test_calc_quantile(self):
        self.assertEqual(gridpp.calc_quantile([0, 1, 2], 0), 0)
        self.assertEqual(gridpp.calc_quantile([0, 1, 2], 0.5), 1)
        self.assertEqual(gridpp.calc_quantile([0, 1, 2], 1), 2)
        self.assertEqual(gridpp.calc_quantile([0, np.nan, 2], 1), 2)
        self.assertEqual(gridpp.calc_quantile([0, np.nan, 2], 0), 0)
        self.assertEqual(gridpp.calc_quantile([0, np.nan, 2], 0.5), 1)
        for quantile in [0, 0.5, 1]:
            self.assertTrue(np.isnan(gridpp.calc_quantile([np.nan, np.nan, np.nan], quantile)))
            self.assertTrue(np.isnan(gridpp.calc_quantile([np.nan], quantile)))
        # BUG: This should work:
        # self.assertTrue(np.isnan(gridpp.calc_quantile([], 0.5)))

        self.assertEqual(gridpp.calc_quantile([[0, 1, 2]], 0), [0])
        self.assertEqual(gridpp.calc_quantile([[0, 1, 2]], 0.5), [1])
        self.assertEqual(gridpp.calc_quantile([[0, 1, 2]], 1), [2])
        self.assertEqual(gridpp.calc_quantile([[0, np.nan, 2]], 1), [2])
        self.assertEqual(gridpp.calc_quantile([[0, np.nan, 2]], 0), [0])
        self.assertEqual(gridpp.calc_quantile([[0, np.nan, 2]], 0.5), [1])
        quantile_of_nan_list = gridpp.calc_quantile([[np.nan, np.nan, np.nan]], 0.5)
        self.assertEqual(len(quantile_of_nan_list), 1)
        self.assertTrue(np.isnan(quantile_of_nan_list[0]))

    def test_num_missing_values(self):
        self.assertEqual(gridpp.num_missing_values([[0, np.nan, 1, np.nan]]), 2)
        self.assertEqual(gridpp.num_missing_values([[np.nan, np.nan]]), 2)
        self.assertEqual(gridpp.num_missing_values([[0, 0, 1, 1]]), 0)
        self.assertEqual(gridpp.num_missing_values([[0, np.nan], [1, np.nan]]), 2)
        self.assertEqual(gridpp.num_missing_values([[np.nan, np.nan], [np.nan, np.nan]]), 4)
        self.assertEqual(gridpp.num_missing_values([[0, 0], [1, 1]]), 0)
        self.assertEqual(gridpp.num_missing_values([[]]), 0)

    def test_calc_statistics_2d(self):
        values = np.reshape(np.arange(9), [3, 3])
        output = gridpp.calc_statistic(values, gridpp.Mean)
        np.testing.assert_array_almost_equal(output, [1, 4, 7])

    def test_warning(self):
        gridpp.warning("test")

    def test_error(self):
        with self.assertRaises(Exception) as e:
            gridpp.error("test")

    def test_not_implemented_error(self):
        with self.assertRaises(Exception) as e:
            gridpp.not_implemented_error("test")

    def test_get_index(self):
        self.assertEqual(2, gridpp.get_lower_index(1, [0, 0, 1, 1]))
        self.assertEqual(3, gridpp.get_upper_index(1, [0, 0, 1, 1]))
        self.assertEqual(0, gridpp.get_lower_index(0, [0, 0, 1, 1]))
        self.assertEqual(1, gridpp.get_upper_index(0, [0, 0, 1, 1]))


if __name__ == '__main__':
    unittest.main()
