import unittest

import numpy as np
from utils import TempFile

from doublebeam.core.models import LinearVelocityLayer, ConstantVelocityLayer, VelocityModel1D, evaluate_at_linear


class TestLinearVelocityLayer(unittest.TestCase):

    def test_evaluation_of_layer(self):
        layer = np.array((0, 10, 3., 6., 2., 4., 2.5, 3.5), dtype=LinearVelocityLayer)
        p_vel = evaluate_at_linear(layer, 5, prop="p")
        self.assertEqual(p_vel, 4.5)


class TestLoadVelocityModel1DLinearVelocity(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # create temp file with example model data
        content = """#Depth(km), vp (km/s), vs (km/s), density (g/cm^3)
                     0.00,5.8000,3.3600,2.700
                     ,5.8000,3.4,2750
                     10.00,5.8000,3.3600,2.750
                     ,6.0,4.0,2.75
                     35.00,8.0400,4.4700,2.800"""
        cls.path, cls.file = TempFile(content)

    def test_file_load_success(self):
        pass


class TestLoadVelocityModel1D(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # create temp file with example model data
        content = """#Depth(km), vp (km/s), vs (km/s), density (g/cm^3)
                     0.00,10.,5.8000,3.3600,2.700
                     10.00.35.,5.8000,3.3600,2.750
                     35.00,50.,8.0400,4.4700,2.800"""
        cls.path, cls.f = TempFile(content)
        cls.layers = np.array([(0, 10, 5.8, 3.36, 2.7), (10., 35, 5.8, 3.36, 2.75),
                               (35., 50, 8.04, 4.47, 2.8)], dtype=ConstantVelocityLayer)

    def setUp(self):
        self.v = VelocityModel1D(self.layers)

    def test_file_load_success(self):
        """
        Test if loading from a file populates the data structure correctly
        """
        np.testing.assert_array_equal(self.v.layers, self.layers)

    def test_prop_access_success(self):
        """
        Test successful access to the models properties
        """
        self.assertEqual(self.v.eval_at(5, "r"), 2.7)
        self.assertEqual(self.v.eval_at(555, "s"), 4.47)

    def test_prop_access_out_of_range(self):
        """
        Test for correct exception when wrong depth is used
        """
        with self.assertRaises(LookupError):
            self.v.eval_at(-555, "p")
