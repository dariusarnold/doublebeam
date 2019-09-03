import math
import unittest
from math import radians
from pathlib import Path

import numpy as np

from doublebeam.raytracing.initial_value import VelocityModel3D, KinematicRayTracer3D, DynamicRayTracer3D
from doublebeam.raytracing.ray import Ray3D, GaussBeam


class TestRay3D(unittest.TestCase):

    def setUp(self) -> None:
        self.angles_down = (0, 5, 10, 30, 42, 44.99999)

    def test_direction_up(self):
        angles_up = [x + 90 for x in self.angles_down if x != 0]
        rays_up = [Ray3D.from_angle(np.array((0, 0, 0)), radians(i), 0, 2000) for i in angles_up]
        for ray, theta in zip(rays_up, angles_up):
            with self.subTest(ray=ray, theta=theta):
                self.assertEqual(ray.direction, "up", msg=f"Failed test for "
                f"upgoing ray with angle {math.degrees(theta)}")

    def test_direction_down(self):

        rays_down = [Ray3D.from_angle(np.array((0, 0, 0)), radians(i), 0, 2000) for i in self.angles_down]
        for ray, theta in zip(rays_down, self.angles_down):
            with self.subTest(ray=ray):
                self.assertEqual(ray.direction, "down", msg=f"Failed test for "
                f"downgoing ray with angle {math.degrees(theta)}")

    def test_direction_horizontal(self):
        ray_horizontal = Ray3D.from_angle(np.array((0, 0, 0)), radians(90), 0, 2000)
        self.assertEqual(ray_horizontal.direction, "horizontal")


class TestRayTracingScipy(unittest.TestCase):

    def test_last_point(self):
        """
        This just tests if the last point of the ray is the same as
        previously calculated. Changes in the structure of the function should
        not change the result.
        """
        last_point_expected = (9403.354242360037, 0, 0)
        layers = VelocityModel3D.convert_to_gradient_intercept([(0, 1000, 3000, 4000),
                                                                (1000, 101000, 6000, 156000)])
        vm = VelocityModel3D(layers)
        ray = Ray3D.from_angle(np.array((0, 0, 0)), radians(20), 0, vm.eval_at(0))
        tracer = KinematicRayTracer3D(vm)
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


    def test_ray_path(self):
        """
        Test if the same endpoint is reached as for the rays given in fig. 9
        from the below publication.
        Slowness values are computed by the two point ray tracing algorithm
        introduced in "A fast and robust two point ray tracing method in layered
        media with constant or linearly varying layer velocity" (Fang, Chen).
        """
        slownesses = [
            [0.000166674323178,  0., 0.0005299638150872],
            [0.0002545148149717, 0., 0.0004938260668176],
            [0.0003320005004714, 0., 0.0004454409534331],
            [0.000333271179152, 0., 0.0004444910532905]
        ]
        source = np.array((0, 0, 0))
        expected_endpoints = (469, 868, 2159, 2411)
        expected_endpoints = [np.array((x, 0, 0)) for x in expected_endpoints]
        vm = VelocityModel3D.from_file("/home/darius/git/doublebeam/fang2019model.txt")
        ray_tracer = KinematicRayTracer3D(vm)
        for slowness, target in zip(slownesses, expected_endpoints):
            with self.subTest(target=target):
                ray = Ray3D(source, np.array(slowness))
                ray_tracer.trace_stack(ray, "TTTTRTTTT")
                np.testing.assert_allclose(ray.last_point, target, atol=1e-6)


class TestRayTracingSingleLayer(unittest.TestCase):

    def setUp(self) -> None:
        layer_gradient = [(0, 100, 2000, 1)]
        layer_constant = [(0, 100, 2000, 0)]
        self.vm_gradient = VelocityModel3D(layer_gradient)
        self.vm_constant = VelocityModel3D(layer_constant)
        self.rt_gradient = KinematicRayTracer3D(self.vm_gradient)
        self.rt_constant = KinematicRayTracer3D(self.vm_constant)

    def test_constant_layer_reflect_bottom(self):
        """
        Shoot a ray downward and reflect from bottom in constant velocity layer
        """
        ray = Ray3D.from_angle((0, 0, 0), radians(45), radians(45),
                               self.vm_constant.eval_at(0))
        self.rt_constant.trace_stack(ray, "R")
        np.testing.assert_allclose(ray.last_point, np.array((math.sqrt(2)*100, math.sqrt(2)*100, 0.)))

    def test_gradient_layer_reflect_bottom(self):
        """
        Shoot a ray downward and reflect from bottom in velocity gradient layer
        """
        ray = Ray3D.from_angle((0, 0, 0), radians(45), radians(0),
                               self.vm_gradient.eval_at(0))
        self.rt_gradient.trace_stack(ray, "R")
        np.testing.assert_allclose(ray.last_point, np.array((210.54093570105533, 0, 0)),
                                   atol=1E-14)

    def test_constant_layer_reflect_top(self):
        """
        Shoot a ray upward and reflect from top in constant velocity layer
        """
        ray = Ray3D.from_angle((0, 0, 100), radians(135), radians(0),
                               self.vm_constant.eval_at(100))
        self.rt_constant.trace_stack(ray, "R")
        np.testing.assert_allclose(ray.last_point, np.array((200, 0, 100)), atol=1E-14)

    def test_gradient_layer_reflect_top(self):
        """
        Shoot a ray upward and reflect from top in velocity gradient layer
        """
        ray = Ray3D.from_angle((0, 0, 100), radians(135), radians(0),
                               self.vm_gradient.eval_at(100))
        self.rt_gradient.trace_stack(ray, "R")
        np.testing.assert_allclose(ray.last_point, np.array((190.89968, 0, 100.)),
                                   atol=1E-14)


class TestExceptionWhenOutsideVerticalBoundaries(unittest.TestCase):

    def setUp(self) -> None:
        layers = [(0, 100, 1800, 0), (100, 200, 2400, 0), (200, 300, 2400, 0),
                  (300, 400, 2700, 0), (400, 500, 2250, 0)]
        self.vm = VelocityModel3D(layers)
        self.ray_tracer = KinematicRayTracer3D(self.vm)

    def test_raise_when_above_model(self):
        """
        Test if a ValueError is raised when ray starts above the model.
        """
        top, _ = self.vm.vertical_boundaries()
        too_high = np.array((0, 0, top-1))
        too_high = Ray3D(too_high, np.array((0, 0, 0)))
        with self.assertRaises(ValueError):
            self.ray_tracer.trace_stack(too_high)

    def test_raise_when_below_model(self):
        """
        Test if a ValueError is raised when ray starts below the model.
        """
        _, bottom = self.vm.vertical_boundaries()
        too_low = np.array((0, 0, bottom+1))
        too_low = Ray3D(too_low, np.array((0, 0, 0)))
        with self.assertRaises(ValueError):
            self.ray_tracer.trace_stack(too_low)


class TestAnalyticalRayTracingConstantVelocity(unittest.TestCase):

    def test_correct_arrival(self):
        """
        Test if the analytic ray tracer gives the same endpoint of the ray as the numeric one
        for a medium composed of constant velocity layers.
        """
        layers = [(0, 100, 1800, 0), (100, 200, 2400, 0), (200, 300, 2400, 0),
                  (300, 400, 2700, 0), (400, 500, 2250, 0)]
        vm = VelocityModel3D(layers)
        nrt = KinematicRayTracer3D(vm)
        ray = Ray3D.from_angle(np.array((0, 0, 0)), radians(17.4576), radians(0), vm.eval_at(0))
        nrt.trace_stack(ray, "TTTTRTTTT")
        expected_last_point = [4.19155952e+02,  0.00000000e+00, -3.05311332e-15]
        np.testing.assert_array_almost_equal(ray.last_point, expected_last_point)


class TestDynamicRayTracingSameResultAsKinematic(unittest.TestCase):
    """
    This tests if the dynamic ray tracer achieves the same result for the ray
    path as the kinematic one.
    # TODO find way to parametrize the test instead of duplicating
    """

    def test_last_point(self):
        """
        This just tests if the last point of the ray is the same as
        previously calculated. Changes in the structure of the function should
        not change the result.
        """
        last_point_expected = (9403.354242360037, 0, 0)
        layers = VelocityModel3D.convert_to_gradient_intercept([(0, 1000, 3000, 4000),
                                                                (1000, 101000, 6000, 156000)])
        vm = VelocityModel3D(layers)
        beam = GaussBeam.from_angle(np.array((0, 0, 0)), radians(20), 0, vm.eval_at(0), 10, 40)
        tracer = DynamicRayTracer3D(vm)
        tracer.trace_stack(beam, "TT", 10)
        self.assertAlmostEqual(last_point_expected[0], beam.last_point[0], places=4,
                               msg=f"x position wrong, got {beam.last_point[0]},"
                               f" expected {last_point_expected[0]}")
        self.assertAlmostEqual(last_point_expected[1], beam.last_point[1], places=4,
                               msg=f"z position wrong, got {beam.last_point[1]},"
                               f" expected {last_point_expected[1]}")
        self.assertAlmostEqual(last_point_expected[2], beam.last_point[2], places=4,
                               msg=f"z position wrong, got {beam.last_point[2]},"
                               f" expected {last_point_expected[2]}")
        self.assertEqual(len(beam.path), 3, msg="Wrong nuber of ray paths")
        self.assertAlmostEqual(beam.last_time, 1.864023576678209, places=4,
                               msg="Wrong travel time for ray")


    def test_ray_path(self):
        """
        Test if the same endpoint is reached as for the rays given in fig. 9
        from the below publication.
        Slowness values are computed by the two point ray tracing algorithm
        introduced in "A fast and robust two point ray tracing method in layered
        media with constant or linearly varying layer velocity" (Fang, Chen).
        """
        slownesses = [
            [0.000166674323178,  0., 0.0005299638150872],
            [0.0002545148149717, 0., 0.0004938260668176],
            [0.0003320005004714, 0., 0.0004454409534331],
            [0.000333271179152, 0., 0.0004444910532905]
        ]
        source = np.array((0, 0, 0))
        expected_endpoints = (469, 868, 2159, 2411)
        expected_endpoints = [np.array((x, 0, 0)) for x in expected_endpoints]
        vm = VelocityModel3D.from_file("/home/darius/git/doublebeam/fang2019model.txt")
        ray_tracer = DynamicRayTracer3D(vm)
        for slowness, target in zip(slownesses, expected_endpoints):
            with self.subTest(target=target):
                beam = GaussBeam(source, np.array(slowness), 10, 40)
                ray_tracer.trace_stack(beam, "TTTTRTTTT")
                np.testing.assert_allclose(beam.last_point, target, atol=1e-6)


class TestDynamicRayTracingOneLayer(unittest.TestCase):

    def setUp(self) -> None:

        self.vm = VelocityModel3D([(0, 10, 2000, 1)])
        self.drt = DynamicRayTracer3D(self.vm)

    def test_credibility(self):
        """
        Basic checking of output format
        """
        source = np.array((0, 0, 0))
        beam = GaussBeam.from_angle(source, radians(20), radians(0), self.vm.eval_at(source), 10, 40)
        self.drt.trace_stack(beam)
        P, Q = beam.P, beam.Q
        self.assertEqual(len(P), 1)
        self.assertEqual(len(Q), 1)
        self.assertEqual(P[0].shape, Q[0].shape)


class TestForRegressionDynamicRayTracing(unittest.TestCase):

    def setUp(self) -> None:
        self.P_desired = np.load(Path("tests/data/P_analytic.npy"))
        self.Q_desired = np.load(Path("tests/data/Q_analytic.npy"))
        vm = VelocityModel3D([(0, 10, 2000, 1)])
        self.drt = DynamicRayTracer3D(vm)
        self.beam = GaussBeam.from_angle((0, 0, 0), radians(20), radians(0), vm.eval_at(0), 10, 40)
        self.drt.trace_stack(self.beam)
        self.P_actual, self.Q_actual = self.beam.P, self.beam.Q

    def test_P_by_comparison(self):
        np.testing.assert_allclose(self.P_actual, self.P_desired)

    def test_Q_by_comparison(self):
        np.testing.assert_allclose(self.Q_actual, self.Q_desired)


class TestDynamicRayTracingMultipleLayers(unittest.TestCase):

    def setUp(self) -> None:
        self.vm = VelocityModel3D([(0, 10, 2000, 1),
                                   (10, 20, 2000, -1)])
        self.drt = DynamicRayTracer3D(self.vm)

    def test_credibility(self):
        """
        Basic checking of output format
        """
        source = np.array((0, 0, 0.))
        beam = GaussBeam.from_angle(source, radians(20), radians(0), self.vm.eval_at(source), 10, 40)
        self.drt.trace_stack(beam, "TRT")
        P, Q = beam.P, beam.Q
        self.assertEqual((4, 4), (len(P), len(Q)))
        for p, q in zip(P, Q):
            with self.subTest():
                self.assertEqual(p.shape, q.shape)

    def test_regression(self):
        """
        Load previous result from file and compare with current one.
        """
        def generate_data():
            source = np.array((0, 0, 0.))
            beam = GaussBeam.from_angle(source, radians(20), radians(0), self.vm.eval_at(source), 10, 40)
            self.drt.trace_stack(beam, "TRT")
            return beam.P, beam.Q
        P_expected = list(np.load(Path("tests/data/P_multilayer.npy"), allow_pickle=True))
        Q_expected = list(np.load(Path("tests/data/Q_multilayer.npy"), allow_pickle=True))
        P_actual, Q_actual = generate_data()

        # do segment wise comparison of the matrices
        with self.subTest("Comparing matrix P"):
            for segment_actual, segment_expected in zip(P_actual, P_expected):
                np.testing.assert_allclose(segment_actual, segment_expected)
        with self.subTest("Comparing matrix Q"):
            for segment_actual, segment_expected in zip(Q_actual, Q_expected):
                np.testing.assert_allclose(segment_actual, segment_expected)
