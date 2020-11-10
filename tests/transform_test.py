from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_identity(self):
        """ Check that identity gives same answers when transformed """
        transform = gridpp.Identity()
        for ar in [1, [1, 1], [[1, 1], [1, 1], [1, 1]], [[[1,1],[1,1]], [[1,1],[1,1]]]]:
            np.testing.assert_equal(ar, transform.forward(ar))
            np.testing.assert_equal(ar, transform.backward(ar))

    def test_boxcox(self):
        threshold = 0.1
        transform = gridpp.BoxCox(threshold)
        input = [0, 1, 2, 3]
        answer = [-10, 0, 0.7177340984, 1.1612319946]
        for k in range(len(input)):
            # Check scalars
            self.assertAlmostEqual(answer[k], transform.forward(input[k]), 5)
            self.assertAlmostEqual(input[k], transform.backward(answer[k]), 5)

            # Check vector versions
            shapes = [[1], [2,2], [2, 0], [3,3,3], [3, 3, 0]]
            for shape in shapes:
                i = input[k] * np.ones(shape)
                a = answer[k] * np.ones(shape)
                output = transform.forward(i)
                np.testing.assert_almost_equal(output, a, 5)

                output = transform.backward(a)
                np.testing.assert_almost_equal(output, i, 5)

    def test_missing_values(self):
        """ Check that invalid values are transformed properly """
        transforms = [gridpp.BoxCox(0.1), gridpp.Log()]
        for transform in transforms:
            with self.subTest(structure=type(transform)):
                input = [1, np.nan, 3]
                output = transform.forward(input)
                np.testing.assert_equal(np.isnan(output), np.isnan(input))
                output = transform.backward(output)
                np.testing.assert_almost_equal(input, output, 5)

    def test_log(self):
        transform = gridpp.Log()
        input = [np.exp(-1), 1, np.exp(1)]
        answer = [-1, 0, 1]
        for k in range(len(input)):
            # Check scalars
            self.assertAlmostEqual(answer[k], transform.forward(input[k]), 5)
            self.assertAlmostEqual(input[k], transform.backward(answer[k]), 5)

            # Check vector versions
            shapes = [[1], [2,2], [2, 0], [3,3,3], [3, 3, 0]]
            for shape in shapes:
                i = input[k] * np.ones(shape)
                a = answer[k] * np.ones(shape)
                output = transform.forward(i)
                np.testing.assert_almost_equal(output, a, 5)

                output = transform.backward(a)
                np.testing.assert_almost_equal(output, i, 5)


if __name__ == '__main__':
    unittest.main()
