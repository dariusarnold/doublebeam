import unittest

import numpy as np
from numpy.testing import assert_array_equal
from utils import TempFile

from doublebeam.core.models import LinearVelocityLayer, ConstantVelocityLayer, VelocityModel1D


class TestVelocityModel1DLinearLayers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # create temp file with example model data
        cls.content = """#Depth top(km), Depth bottom (km), vp top (km/s), bp bottom (km/s), vs top (km/s), vs bottom (km/s), density top (g/cm^3), density bottom (g/cm^3)
                     0.00, 10., 5.8000, 5.8, 3.3600, 3.4, 2.700, 2.75
                     10.00, 35, 5.9000, 6.0, 3.4, 4.0, 2.750, 2.9"""
        cls.path, cls.file = TempFile(cls.content)
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

    def test_init_converting(self):
        """Test if creating a VelocityModel from a sequence other than the expected array with dtype
        LinearVelocityLayer works"""
        layers_list = list(self.layers)
        velocity_model_from_list = VelocityModel1D(layers_list)
        assert_array_equal(self.v.layers, velocity_model_from_list.layers)

    def test_file_load(self):
        loaded_model = VelocityModel1D.from_file(self.path)
        assert_array_equal(loaded_model.layers, self.v.layers)

    def test_string_load(self):
        mod = VelocityModel1D.from_string(self.content)
        assert_array_equal(self.v.layers, mod.layers)

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

    def test_one_interface_crossed_returns_true(self):
        # both orders of depth should work:
        self.assertEqual(self.v.interface_crossed(9, 11), True)
        self.assertEqual(self.v.interface_crossed(10.5, 9.5), True)
        # even small steps across the interface should be detected
        z1, z2 = np.nextafter(10, 9), np.nextafter(10, 11)
        if z1 != z2:
            self.assertEqual(self.v.interface_crossed(z1, z2), True)
            self.assertEqual(self.v.interface_crossed(z1, z2), self.v.interface_crossed(z2, z1))
        else:
            print("Assertion skipped due to float precision")
        # top and bottom of the model also count as interfaces
        self.assertEqual(self.v.interface_crossed(-0.1, 2), True)
        self.assertEqual(self.v.interface_crossed(34.99, 36), True)

    def test_multiple_interfaces_crossed_returns_true(self):
        """As long as at least one interface between the depths, the function
        should return True"""
        self.assertEqual(self.v.interface_crossed(9, 36), True)
        self.assertEqual(self.v.interface_crossed(36, 9), True)
        self.assertEqual(self.v.interface_crossed(-1, 1000), True)

    def test_no_interfaces_crossed_returns_false(self):
        self.assertEqual(self.v.interface_crossed(1, 2), False)
        z1 = np.nextafter(10, 11)
        z2 = np.nextafter(z1, 11)
        self.assertEqual(self.v.interface_crossed(z1, z2), False)


class TestVelocityModel1DConstantLayers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # create temp file with example model data
        cls.content = """#Depth(km), vp (km/s), vs (km/s), density (g/cm^3)
                     0.00, 10., 5.8000, 3.3600, 2.700
                     10.00, 35., 5.8000, 3.3600, 2.750
                     35.00, 50., 8.0400, 4.4700, 2.800"""
        cls.path, cls.f = TempFile(cls.content)
        # manually parsed to compare
        cls.layers = np.array([(0, 10, 5.8, 3.36, 2.7), (10., 35, 5.8, 3.36, 2.75),
                               (35., 50, 8.04, 4.47, 2.8)], dtype=ConstantVelocityLayer)

    def setUp(self):
        self.v = VelocityModel1D(self.layers)

    def test_init(self):
        """Test if creating a VelocityModel populates the data structure correctly"""
        assert_array_equal(self.v.layers, self.layers)

    def test_init_converting(self):
        """Test if creating a VelocityModel from a sequence other than the expected array with dtype
        ConstantVelocityLayer works"""
        layers_list = list(self.layers)
        velocity_model_from_list = VelocityModel1D(layers_list)
        assert_array_equal(self.v.layers, velocity_model_from_list.layers)

    def test_file_load(self):
        loaded_model = VelocityModel1D.from_file(self.path)
        assert_array_equal(loaded_model.layers, self.v.layers)

    def test_string_load(self):
        mod = VelocityModel1D.from_string(self.content)
        assert_array_equal(self.v.layers, mod.layers)

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

    def test_one_interface_crossed_returns_true(self):
        # both orders of depth should work:
        self.assertEqual(self.v.interface_crossed(9, 11), True)
        self.assertEqual(self.v.interface_crossed(10.5, 9.5), True)
        # even small steps across the interface should be detected
        z1, z2 = np.nextafter(35, 34), np.nextafter(35, 36)
        if z1 != z2:
            self.assertEqual(self.v.interface_crossed(z1, z2), True)
            self.assertEqual(self.v.interface_crossed(z1, z2), self.v.interface_crossed(z2, z1))
        else:
            print("Assertion skipped due to float precision")
        # top and bottom of the model also count as interfaces
        self.assertEqual(self.v.interface_crossed(-0.1, 2), True)
        self.assertEqual(self.v.interface_crossed(49.99, 51), True)

    def test_multiple_interfaces_crossed_returns_true(self):
        """As long as at least one interface between the depths, the function
        should return True"""
        self.assertEqual(self.v.interface_crossed(9, 36), True)
        self.assertEqual(self.v.interface_crossed(36, 9), True)
        self.assertEqual(self.v.interface_crossed(-1, 1000), True)

    def test_no_interfaces_crossed_returns_false(self):
        self.assertEqual(self.v.interface_crossed(1, 2), False)
        z1 = np.nextafter(10, 11)
        z2 = np.nextafter(z1, 11)
        self.assertEqual(self.v.interface_crossed(z1, z2), False)
