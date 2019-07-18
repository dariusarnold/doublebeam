import math
import unittest
from math import radians

import numpy as np

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
        last_point_expected = (9403.354242360037, 0, 0)
        layers = VelocityModel3D.convert_to_gradient_intercept([(0, 1000, 3000, 4000),
                                                                (1000, 101000, 6000, 156000)])
        vm = VelocityModel3D(layers)
        ray = Ray3D(0, 0, 0, radians(20), 0)
        tracer = NumericRayTracer3D(vm)
        tracer.trace_stack(ray, "TT", 10)
        self.assertAlmostEqual(last_point_expected[0], ray.last_point[0], places=4,
                               msg=f"x position wrong, got {ray.last_point[0]},"
                               f" expected {last_point_expected[0]}")
        self.assertAlmostEqual(last_point_expected[1], ray.last_point[1], places=4,
                               msg=f"z position wrong, got {ray.last_point[1]},"
                               f" expected {last_point_expected[1]}")
        self.assertAlmostEqual(last_point_expected[2], ray.last_point[2], places=4,
                               msg=f"z position wrong, got {ray.last_point[2]},"
                               f" expected {last_point_expected[2]}")
        self.assertEqual(len(ray.path), 3, msg="Wrong nuber of ray paths")
        self.assertAlmostEqual(ray.last_time, 1.864023576678209, places=4,
                               msg="Wrong travel time for ray")


class TestAnalyticalRayTracingConstantVelocity(unittest.TestCase):

    def test_correct_arrival(self):
        """
        Test if the analytic ray tracer gives the same endpoint of the ray as the numeric one
        for a medium composed of constant velocity layers.
        """
        layers = [(0, 100, 1800, 0), (100, 200, 2400, 0), (200, 300, 2400, 0),
                  (300, 400, 2700, 0), (400, 500, 2250, 0)]
        vm = VelocityModel3D(layers)
        nrt = NumericRayTracer3D(vm)
        ray = Ray3D(0, 0, 0, radians(17.4576), radians(0))
        nrt.trace_stack(ray, "TTTTRTTTT")
        expected_last_point = [4.19155952e+02,  0.00000000e+00, -3.05311332e-15]
        np.testing.assert_array_almost_equal(ray.last_point, expected_last_point)
