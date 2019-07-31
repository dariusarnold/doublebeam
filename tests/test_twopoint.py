import itertools
import unittest
from math import radians, sin
from typing import Tuple

import numpy as np

from doublebeam.models import VelocityModel3D
from doublebeam.raytracing.twopoint import TwoPointRayTracing
from doublebeam.raytracing.initial_value import NumericRayTracer3D
from doublebeam.raytracing.ray import Ray3D


class TestTwoPointRayTracing(unittest.TestCase):

    def setUp(self) -> None:
        self.vm = VelocityModel3D.from_file("/home/darius/git/double-beam/fang2019model.txt")
        self.rt = TwoPointRayTracing(self.vm)
        self.sources = [np.array((469/2, 0, 500)),
                        np.array((868/2, 0, 500)),
                        np.array((2159/2, 0, 500))]
        self.receivers = [np.array((469, 0, 0)),
                          np.array((868, 0, 0)),
                          np.array((2159, 0, 0))]

    def test_correct_values(self):
        """
        Test if two point ray tracing gives the same slowness as calculated
        from fig. 9 in Fang2019.
        """
        slownesses = [self.rt.trace(source, receiver) for source, receiver in zip(self.sources, self.receivers)]
        slownesses_expected = [sin(radians(i)) / 3000 for i in (30, 50, 85)]
        for p, p_expected in zip(slownesses, slownesses_expected):
            with self.subTest(p=p, p_expected=p_expected):
                self.assertAlmostEqual(p[0], p_expected, places=5)

    def test_source_receiver_swap(self):
        """
        Swapping receiver and source should give the same result.
        """
        slownesses = [self.rt.trace(receiver, source) for source, receiver in zip(self.sources, self.receivers)]
        # negative slowness since the wave travels opposite to the x direction
        slownesses_expected = [-sin(radians(i)) / 3000 for i in (30, 50, 85)]
        for p, p_expected in zip(slownesses, slownesses_expected):
            with self.subTest(p=p, p_expected=p_expected):
                self.assertAlmostEqual(p[0], p_expected, places=5)

    def test_3D_case(self):
        """
        Test if the algorithm works with a general 3D problem where the ray
        moves outside of the xz-plane.
        """
        source = np.array((0, 0, 500))
        receiver = np.array((100, 200, 0))
        slowness = self.rt.trace(source, receiver)
        nrt = NumericRayTracer3D(self.vm)
        ray = Ray3D(source, slowness)
        nrt.trace_stack(ray, "TTTT")
        np.testing.assert_allclose(ray.last_point, receiver, atol=1e-10)

    def test_3D_case_swapped(self):
        """
        Test if the algorithm works with a general 3D problem with the source
        above the receiver
        """
        receiver = np.array((0, 0, 500))
        source = np.array((100, 200, 0))
        slowness = self.rt.trace(source, receiver)
        nrt = NumericRayTracer3D(self.vm)
        ray = Ray3D(source, slowness)
        nrt.trace_stack(ray, "TTTT")
        np.testing.assert_allclose(ray.last_point, receiver, atol=1e-9)

    def test_raise_when_out_of_depth(self):
        """
        Test if a ValueError is raised when either source, receiver or both are
        outside of the vertical boundaries of the model.
        """
        top, bottom = self.vm.vertical_boundaries()
        to_high = np.array((0, 0, top-1))
        to_low = np.array((0, 0, bottom+1))
        in_model = np.array((0, 0, top + (bottom-top)/2))
        permutations = itertools.permutations((to_high, to_low, in_model), 2)
        for point_a, point_b in permutations:
            with self.subTest(point_a=point_a, point_b=point_b):
                with self.assertRaises(ValueError):
                    self.rt.trace(point_a, point_b)


class TestMethod_v_M(unittest.TestCase):

    def setUp(self) -> None:
        layers = [(0, 100, 2600, -4), (100, 200, 2400, 0), (200, 300, 2400, 1),
                  (300, 400, 2700, 0), (400, 500, 2250, 1.5)]
        self.vm = VelocityModel3D(layers)
        self.tp = TwoPointRayTracing(self.vm)
        # format for values: (top depth, bottom depth, expected return value
        test_data = [(0, 500, 3000),
                     (50, 500, 3000),
                     (50, 450, 2925),
                     (300, 400, 2850),
                     (50, 150, 2400),
                     (25, 150, 2500),
                     (301, 302, 2700),
                     (0, 1, 2600),
                     (50, 250, 2650)]
        self.test_data = [(*self.convert_testdata(z1, z2), v) for z1, z2, v in test_data]

    @staticmethod
    def convert_testdata(z1: float, z2: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Helper method to convert depth values to full coordinates represented
        as np array of 3 values with depth as last entry
        """
        return np.array((0, 0, z1)), np.array((0, 0, z2))

    @staticmethod
    def msg_(p1: np.ndarray, p2: np.ndarray, got_velocity: float,
             expected_velocity: float) -> str:
        return (f"Error between {p1}, {p2}: expected {expected_velocity} m/s, "
                f"got {got_velocity} m/s")

    def test_method(self):
        for p1, p2, expected in self.test_data:
            with self.subTest(p1=p1, p2=p2, expected=expected):
                actual = self.tp._v_M(p1, p2)
                self.assertEqual(actual, expected, msg=self.msg_(p1, p2, actual, expected))

