import unittest

import numpy as np
from numpy.testing import assert_array_equal
from utils import TempFile

from doublebeam.core.models import LinearVelocityLayer, ConstantVelocityLayer, VelocityModel1D


class TestVelocityModel1DLinearLayers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # create temp file with example model data
        content = """#Depth top(km), Depth bottom (km), vp top (km/s), bp bottom (km/s), vs top (km/s), vs bottom (km/s), density top (g/cm^3), density bottom (g/cm^3)
                     0.00, 10., 5.8000, 5.8, 3.3600, 3.4, 2.700, 2.75
                     10.00, 35, 5.9000, 6.0, 3.4, 4.0, 2.750, 2.9"""
        cls.path, cls.file = TempFile(content)
        # manually parsed to compare
        cls.layers = np.array([(0., 10., 5.8, 5.8, 3.36, 3.4, 2.7, 2.75),
                               (10., 35., 5.9, 6., 3.4, 4., 2.75, 2.9)], dtype=LinearVelocityLayer)

    @classmethod
    def tearDownClass(cls):
        cls.file.close()

    def setUp(self):
        self.v = VelocityModel1D(self.layers)

    def test_init(self):
        """Test if creating a VelocityModel populates the data structure correctly"""
        assert_array_equal(self.v.layers, self.layers)

    def test_file_load(self):
        loaded_model = VelocityModel1D.from_file(self.path)
        assert_array_equal(loaded_model.layers, self.v.layers)

    def test_prop_access_success(self):
        msg = "{} evaluation failed or incorrect"
        self.assertEqual(self.v.eval_at(0, "d"), 2.7, msg=msg.format("Density"))
        self.assertEqual(self.v.eval_at(5, "s"), 3.38, msg=msg.format("S velocity"))
        self.assertEqual(self.v.eval_at(1, "p"), 5.8, msg=msg.format("P velocity"))

    def test_layer_top_inclusivity(self):
        # since upper boundary is inclusive, value of lower layer should be
        # returned here
        self.assertEqual(self.v.eval_at(10., "p"),  5.9)

    def test_raises_when_depth_out_of_model_range(self):
        with self.assertRaises(LookupError, msg="Model evaluation at negative depth not raising correct exception"):
            self.v.eval_at(-5.5, "p")
        with self.assertRaises(LookupError, msg="Model evaluation below lowest layer not raising correct exception"):
            self.v.eval_at(999.9, "p")

    def test_interface_depths_correct(self):
        assert_array_equal(self.v.interface_depths, np.array((0., 10., 35.)))


class TestVelocityModel1DConstantLayers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # create temp file with example model data
        content = """#Depth(km), vp (km/s), vs (km/s), density (g/cm^3)
                     0.00, 10., 5.8000, 3.3600, 2.700
                     10.00, 35., 5.8000, 3.3600, 2.750
                     35.00, 50., 8.0400, 4.4700, 2.800"""
        cls.path, cls.f = TempFile(content)
        # manually parsed to compare
        cls.layers = np.array([(0, 10, 5.8, 3.36, 2.7), (10., 35, 5.8, 3.36, 2.75),
                               (35., 50, 8.04, 4.47, 2.8)], dtype=ConstantVelocityLayer)

    def setUp(self):
        self.v = VelocityModel1D(self.layers)

    def test_init(self):
        """Test if creating a VelocityModel populates the data structure correctly"""
        assert_array_equal(self.v.layers, self.layers)

    def test_file_load(self):
        loaded_model = VelocityModel1D.from_file(self.path)
        assert_array_equal(loaded_model.layers, self.v.layers)

    def test_prop_access_success(self):
        """
        Test successful access to the models properties
        """
        msg = "{} evaluation failed or incorrect"
        self.assertEqual(self.v.eval_at(5, "d"), 2.7, msg=msg.format("Density"))
        self.assertEqual(self.v.eval_at(40, "s"), 4.47, msg=msg.format("S velocity"))
        self.assertEqual(self.v.eval_at(20., "p"), 5.8, msg=msg.format("P velocity"))

    def test_layer_top_inclusivity(self):
        # evaluation at layer border should yield value of lower layer
        self.assertEqual(self.v.eval_at(35., "p"), 8.04)

    def test_raises_when_depth_out_of_model_range(self):
        """
        Test for correct exception when wrong depth (outside of models depth
         range) is used
        """
        with self.assertRaises(LookupError, msg="Model evaluation at negative depth not raising correct exception"):
            self.v.eval_at(-555, "p", )
        with self.assertRaises(LookupError, msg="Model evaluation below lowest layer not raising correct exception"):
            self.v.eval_at(999, "s")

    def test_interface_depths_correct(self):
        assert_array_equal(self.v.interface_depths, np.array((0., 10., 35., 50.)))
