from __future__ import print_function
import unittest
import gridpp
import numpy as np
import os
# import psutil


def memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss


class MemoryTest(unittest.TestCase):
    def test_memory_leak_objects(self):
        """Checks if there is a memory leak when creating gridpp objects"""
        return
        N = 1000
        last = None
        min = None
        max = None
        for i in range(5):
            lats = np.random.rand(N * N)
            lons = np.random.rand(N * N)
            lats2 = np.random.rand(N, N)
            lons2 = np.random.rand(N, N)
            points = gridpp.Points(lats, lons)
            grid = gridpp.Grid(lats2, lons2)
            curr = memory_usage()
            if min is None or max is None:
                min = curr
                max = curr
            if curr < min:
                min = curr
            if curr > max:
                max = curr
            last = curr
        self.assertTrue(max / min < 2)


    def test_memory_leak_input(self):
        """Checks if there is a memory leak when passing vectors"""
        return
        N = 1000
        last = None
        min = None
        max = None
        for i in range(5):
            vec = np.random.rand(N * N)
            vec2 = np.random.rand(N, N)
            output2 = gridpp.neighbourhood(vec2, 7, gridpp.Mean)
            output = gridpp.dewpoint(vec, vec)
            curr = memory_usage()
            if min is None or max is None:
                min = curr
                max = curr
            if curr < min:
                min = curr
            if curr > max:
                max = curr
            last = curr
        self.assertTrue(max / min < 2)


if __name__ == '__main__':
    unittest.main()
