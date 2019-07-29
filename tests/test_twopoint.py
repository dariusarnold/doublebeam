import unittest
from math import radians, sin

import numpy as np

from doublebeam.core.models import VelocityModel3D
from doublebeam.core.raytracing.twopoint import TwoPointRayTracing
from doublebeam.core.raytracing.initial_value import NumericRayTracer3D
from doublebeam.core.raytracing.ray import Ray3D


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
        ray = Ray3D(*source, slowness)
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
        ray = Ray3D(*source, slowness)
        nrt.trace_stack(ray, "TTTT")
        np.testing.assert_allclose(ray.last_point, receiver, atol=1e-10)