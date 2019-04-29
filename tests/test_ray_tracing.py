import unittest
from math import radians

from doublebeam.core.raytracing import ray_trace_euler, MockVelocityModel1D, Ray2D, ray_trace_scipy


class TestRayTracingEuler(unittest.TestCase):

    def test_last_point(self):
        """This just tests if the last point of the ray is the same as
        previously calculated. Changes in the structure of the function should
        not change the result."""
        last_point_expected = (9.298386897144244, -0.0029832008545231554)
        vm = MockVelocityModel1D(1)
        ray = Ray2D(0, 0, radians(20))
        s_end = 20
        traced_points = ray_trace_euler(ray, vm, s_end)[-1]
        self.assertAlmostEqual(last_point_expected[0], traced_points[0])
        self.assertAlmostEqual(last_point_expected[1], traced_points[1])


class TestRayTracingScipy(unittest.TestCase):

    def test_last_point(self):
        """This just tests if the last point of the ray is the same as
        previously calculated. Changes in the structure of the function should
        not change the result."""
        last_point_expected = (6.506034140591608, -8.331009493378616e-16)
        vm = MockVelocityModel1D(1)
        ray = Ray2D(0, 0, radians(20))
        s_end = 20
        traced_points = ray_trace_scipy(ray, vm, s_end)[-1]
        self.assertAlmostEqual(last_point_expected[0], traced_points[0])
        self.assertAlmostEqual(last_point_expected[1], traced_points[1])
