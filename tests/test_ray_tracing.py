import unittest
from math import radians

from doublebeam.core.raytracing.raytracing import VelocityModel1D, Ray2D, RayTracer2D


class TestRayTracingScipy(unittest.TestCase):

    def test_last_point(self):
        """This just tests if the last point of the ray is the same as
        previously calculated. Changes in the structure of the function should
        not change the result."""
        last_point_expected = (6.458938131869513, -8.331009493378616e-16)
        vm = VelocityModel1D.from_string("0, 1, 3, 4, 0, 0, 1, 1\n1, 101, 6, 156, 0, 0, 1, 1")
        ray = Ray2D(0, 0, radians(20))
        s_end = 20
        tracer = RayTracer2D(vm)
        ray = tracer.ray_trace(ray, s_end)
        self.assertAlmostEqual(last_point_expected[0], ray.path[0][-1])
        self.assertAlmostEqual(last_point_expected[1], ray.path[1][-1])
