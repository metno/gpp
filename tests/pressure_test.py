from __future__ import print_function
import unittest
import gridpp
import numpy as np
import os


class Test(unittest.TestCase):
    def test_missing_values(self):
        self.assertTrue(np.isnan(gridpp.pressure(np.nan, 0, 101325)))
        self.assertTrue(np.isnan(gridpp.pressure(0, np.nan, 101325)))
        self.assertTrue(np.isnan(gridpp.pressure(0, 100, np.nan)))

    def test_no_elev_diff(self):
        for pressure in [101325, 500]:
            self.assertAlmostEqual(gridpp.pressure(0, 0, pressure), pressure)

    def test_standard(self):
        self.assertAlmostEqual(gridpp.pressure(0, 1000, 101325, 288.15), 89996.7, 0)
        self.assertAlmostEqual(gridpp.pressure(1000, 0, 89996.7, 288.15), 101325, 0)

    def test_temperature(self):
        self.assertAlmostEqual(gridpp.pressure(0, 1000, 101325, 258.15), 88765.2, 0)
        self.assertAlmostEqual(gridpp.pressure(1000, 0, 88765.2, 258.15), 101325, 0)

    def test_no_pressure(self):
        self.assertAlmostEqual(gridpp.pressure(0, 0, 0), 0)
        self.assertAlmostEqual(gridpp.pressure(0, 1000, 0), 0)

    def test_no_temperature(self):
        self.assertTrue(np.isnan(gridpp.pressure(0, 0, 0, 0)))
        self.assertTrue(np.isnan(gridpp.pressure(0, 0, 101325, 0)))

    def test_vector(self):
        """Check that vector function gives same answer as scalar"""
        ielev = [0, 100, 200, np.nan]
        oelev = [1000, 900, 800, np.nan]
        pressure = [1e5, 1.1e5, 1.2e5, np.nan]
        temperature = [280, 290, 300, np.nan]
        truth = [gridpp.pressure(ielev[i], oelev[i], pressure[i], temperature[i]) for i in range(len(ielev))]
        np.testing.assert_array_almost_equal(truth, gridpp.pressure(ielev, oelev, pressure, temperature))

    def test_dimension_mismatch(self):
        """Check for exception when vector arguments have different sizes"""
        ielev = [0, 100, 200]
        oelev = [1000, 900, 800]
        pressure = [1e5, 1.1e5, 1.2e5]
        temperature = [280, 290, 300]
        with self.assertRaises(Exception) as e:
            gridpp.pressure(ielev + [100], oelev, pressure, temperature)
        with self.assertRaises(Exception) as e:
            gridpp.pressure(ielev, oelev + [100], pressure, temperature)
        with self.assertRaises(Exception) as e:
            gridpp.pressure(ielev, oelev, pressure + [100], temperature)
        with self.assertRaises(Exception) as e:
            gridpp.pressure(ielev, oelev, pressure, temperature + [100])


if __name__ == '__main__':
    unittest.main()
