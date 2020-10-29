from __future__ import print_function
import unittest
import gridpp
import numpy as np
import psutil
import os


def memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss


class Test(unittest.TestCase):
    def test_int_output(self):
        """Checks if ints are passed back
        """
        output = gridpp.test_ivec_output()
        self.assertEqual(type(output[0]), np.int32)
        output = gridpp.test_ivec2_output()
        self.assertEqual(type(output[0][0]), np.int32)

    def test_float_output(self):
        output = gridpp.test_vec_output()
        self.assertEqual(type(output[0]), np.float32)
        output = gridpp.test_vec2_output()
        self.assertEqual(type(output[0][0]), np.float32)
        output = gridpp.test_vec3_output()
        self.assertEqual(type(output[0][0][0]), np.float32)

    def test_vec_input(self):
        """ Test that list, tuples and numpy arrays work"""
        ar = [1, 2, 3]

        for func in [gridpp.test_vec_input, gridpp.test_ivec_input]:
            self.assertEqual(func(ar), 6)
            self.assertEqual(func((1, 2, 3)), 6)
            self.assertEqual(func(np.array(ar)), 6)
            self.assertEqual(func(np.array(ar).astype('float32')), 6)
            self.assertEqual(func(np.array(ar).astype('float64')), 6)
            self.assertEqual(func(np.array(ar).astype('int32')), 6)

    def test_vec2_input(self):
        ar = [[1,2], [2,3], [3,4]]
        self.assertEqual(gridpp.test_vec2_input(ar), 15)
        self.assertEqual(gridpp.test_vec2_input(np.array(ar)), 15)
        self.assertEqual(gridpp.test_vec2_input(np.array(ar).astype('float32')), 15)
        self.assertEqual(gridpp.test_vec2_input(np.array(ar).astype('float64')), 15)
        self.assertEqual(gridpp.test_vec2_input(np.array(ar).astype('int32')), 15)

    def test_vec3_input(self):
        ar = [[[1,2],[2,3]], [[2,3],[3,4]], [[3,4],[4,5]]]
        self.assertEqual(gridpp.test_vec3_input(ar), 36)
        self.assertEqual(gridpp.test_vec3_input(np.array(ar)), 36)
        self.assertEqual(gridpp.test_vec3_input(np.array(ar).astype('float32')), 36)
        self.assertEqual(gridpp.test_vec3_input(np.array(ar).astype('float64')), 36)
        self.assertEqual(gridpp.test_vec3_input(np.array(ar).astype('int32')), 36)

    def test_vec_argout(self):
        n, distances = gridpp.test_vec_argout()
        self.assertEqual(len(distances), 10)

    def test_vec2_argout(self):
        n, distances = gridpp.test_vec2_argout()
        self.assertEqual(distances.shape[0], 10)
        self.assertEqual(distances.shape[1], 10)

    def test_invalid_dimension_error(self):
        """ Check that an invalid number of dimensions is detected """
        for func in [gridpp.test_vec2_input, gridpp.test_vec3_input]:
            with self.assertRaises(Exception) as e:
                func(np.zeros([5]))
        for func in [gridpp.test_vec_input, gridpp.test_vec3_input]:
            with self.assertRaises(Exception) as e:
                func(np.zeros([5, 2]))
        for func in [gridpp.test_vec_input, gridpp.test_vec2_input]:
            with self.assertRaises(Exception) as e:
                func(np.zeros([5, 2, 3]))

    def test_vec_output(self):
        ar = [-1, -1, -1]
        np.testing.assert_array_equal(gridpp.test_vec_output(), ar)
        np.testing.assert_array_equal(gridpp.test_vec2_output(), [ar, ar, ar])
        np.testing.assert_array_equal(gridpp.test_vec3_output(), [[ar, ar, ar], [ar, ar, ar], [ar, ar, ar]])

    def test_ivec_output(self):
        ar = [-1, -1, -1]
        np.testing.assert_array_equal(gridpp.test_ivec_output(), ar)
        np.testing.assert_array_equal(gridpp.test_ivec2_output(), [ar, ar, ar])
        np.testing.assert_array_equal(gridpp.test_ivec3_output(), [[ar, ar, ar], [ar, ar, ar], [ar, ar, ar]])

    def test_test_array(self):
        """Test test_array to complete code coverage. Not a needed function."""
        ar = [1, 2, 3]
        gridpp.test_array(ar)

    def test_zero_size_dimension(self):
        """Test arrays with a 0-length dimension"""
        self.assertEqual(0, gridpp.test_vec_input(np.zeros([0])))
        self.assertEqual(0, gridpp.test_vec_input([]))

        # 2D
        self.assertEqual(0, gridpp.test_vec2_input(np.zeros([3, 0])))
        self.assertEqual(0, gridpp.test_vec2_input(np.zeros([0, 3])))
        self.assertEqual(0, gridpp.test_vec2_input(np.zeros([0, 0])))

        # 3D
        self.assertEqual(0, gridpp.test_vec3_input(np.zeros([3, 3, 0])))
        self.assertEqual(0, gridpp.test_vec3_input(np.zeros([3, 0, 3])))
        self.assertEqual(0, gridpp.test_vec3_input(np.zeros([0, 3, 3])))
        self.assertEqual(0, gridpp.test_vec3_input(np.zeros([3, 0, 0])))
        self.assertEqual(0, gridpp.test_vec3_input(np.zeros([0, 3, 0])))
        self.assertEqual(0, gridpp.test_vec3_input(np.zeros([3, 0, 0])))


if __name__ == '__main__':
    unittest.main()
