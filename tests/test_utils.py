import unittest
from math import sqrt, radians

import numpy as np

from doublebeam.core.utils import horizontal_distance, angle


class TestHorizontalDistance(unittest.TestCase):

    def test_function(self):
        def error_msg(p1, p2, calculated, expected):
            return (f"Wrong horizontal distance between points {p1}, {p2}. Got "
                f"{calculated}, expected {expected}")

        # data consisting of two test points and their distance
        data = [((0, 0, 0), (0, 0, 0), 0),
                ((1, 0, 0), (0, 0, 0), 1),
                ((0, 0, 0), (0, -1, 2), 1),
                ((0, 0, 1), (0, 0, -1), 0),
                ((-1, -1, -42.5), (1, 1, 2), sqrt(2)*2)]
        # convert points to numpy arrays
        data = [(np.array(p1), np.array(p2), d) for p1, p2, d in data]
        for p1, p2, distance in data:
            with self.subTest(p1=p1, p2=p2, distance=distance):
                calculated = horizontal_distance(p1, p2)
                self.assertEqual(calculated, distance,
                                 msg=error_msg(p1, p2, calculated, distance))


class TestAngle(unittest.TestCase):

    def test_90_degrees(self):
        a = np.array((1, 0, 0))
        b = np.array((0, 1, 0))
        self.assertEqual(angle(a, b), radians(90))
        self.assertEqual(angle(a, b, acute=False), radians(90))

    def test_parallel(self):
        a = np.array((0.8, 0.2, 0.4))
        self.assertEqual(angle(a, a), 0)
        self.assertEqual(angle(a, a, acute=False), radians(180))

    def test_float_clipping(self):
        # these values create an acos argument >1 that needs to be clipped
        a = np.array((0.874469283050132, 0.262553720250597, 0.477968688795641))
        self.assertEqual(angle(a, a), 0)
        self.assertEqual(angle(a, a, acute=False), radians(180))

    def test_acute_obtuse(self):
        a = np.array((1, 0, 0))
        b = np.array((1, 1, 0))
        self.assertAlmostEqual(angle(a, b), radians(45))
        self.assertAlmostEqual(angle(a, b, acute=False), radians(135))
