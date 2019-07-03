import math
import unittest
from math import radians

from doublebeam.core.raytracing.raytracing import VelocityModel3D, Ray3D, NumericRayTracer3D


class TestRay3D(unittest.TestCase):

    def setUp(self) -> None:
        self.angles_down = (0, 5, 10, 30, 42, 44.99999)

    def test_direction_up(self):
        angles_up = [x + 90 for x in self.angles_down if x != 0]
        rays_up = [Ray3D(0, 0, 0, radians(i), 0) for i in angles_up]
        for ray in rays_up:
            with self.subTest(ray=ray):
                self.assertEqual(ray.direction, "up", msg=f"Failed test for "
                f"upgoing ray with angle {math.degrees(ray.theta)}")

    def test_direction_down(self):

        rays_down = [Ray3D(0, 0, 0, radians(i), 0) for i in self.angles_down]
        for ray in rays_down:
            with self.subTest(ray=ray):
                self.assertEqual(ray.direction, "down", msg=f"Failed test for "
                f"downgoing ray with angle {math.degrees(ray.theta)}")

    def test_direction_horizontal(self):
        ray_horizontal = Ray3D(0, 0, 0, radians(90), 0)
        self.assertEqual(ray_horizontal.direction, "horizontal")


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
