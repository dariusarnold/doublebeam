import unittest

import numpy as np
from numpy.testing import assert_array_equal
from utils import TempFile

from doublebeam.core.models import LinearVelocityLayer, VelocityModel3D


class TestVelocityModel3DLinearLayers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # create temp file with example model data
        cls.content = """#Depth top(m), Depth bottom (m), vp top (m/s), vp bottom (m/s)
                     0.00, 10000., 5800, 5800
                     10000.00, 11000, 5000, 6000
                     11000, 12000, 7000, 6000"""
        cls.path, cls.file = TempFile(cls.content)
        # manually parsed to compare
        cls.layers = np.array([(0., 10000., 5800, 5800),
                               (10000., 11000., 5000, 6000),
                               (11000, 12000, 7000, 6000)])
        cls.layers_correct = np.array([(0, 10000, 5800, 0),
                                       (10000, 11000, -5000, 1),
                                       (11000, 12000, 18000, -1)],
                                      dtype=LinearVelocityLayer)

    @classmethod
    def tearDownClass(cls):
        cls.file.close()

    def setUp(self):
        self.v = VelocityModel3D(self.layers_correct)

    def test_converting(self):
        """Test iconverting the values to intercept+gradient"""
        layers_converted = VelocityModel3D.convert_to_gradient_intercept(self.layers)
        assert_array_equal(layers_converted, self.layers_correct)

    def test_file_load(self):
        loaded_model = VelocityModel3D.from_file(self.path)
        assert_array_equal(loaded_model.layers, self.layers_correct)

    def test_string_load(self):
        mod = VelocityModel3D.from_string(self.content)
        assert_array_equal(self.v.layers, mod.layers)

    def test_prop_access_success(self):
        msg = "{} evaluation failed or incorrect"
        self.assertEqual(self.v.eval_at(0, 0, 0), 5800)
        self.assertEqual(self.v.eval_at(0, 0, 5000), 5800)
        self.assertEqual(self.v.eval_at(0, 0, 9999), 5800)
        self.assertEqual(self.v.eval_at(0, 0, 10001), 5001)
        self.assertEqual(self.v.eval_at(0, 0, 10500), 5500)

    def test_layer_top_inclusivity(self):
        # since upper boundary is inclusive, value of lower layer should be
        # returned here
        self.assertEqual(self.v.eval_at(0, 0, 10000),  5000)

    def test_raises_when_depth_out_of_model_range(self):
        with self.assertRaises(LookupError, msg="Model evaluation at negative depth not raising correct exception"):
            self.v.eval_at(0, 0, -5.5)
        with self.assertRaises(LookupError, msg="Model evaluation below lowest layer not raising correct exception"):
            self.v.eval_at(0, 0, 99999999999)

    def test_interface_depths_correct(self):
        assert_array_equal(self.v.interface_depths, np.insert(self.layers_correct["bot_depth"], 0, 0))

    def test_one_interface_crossed_returns_true(self):
        # both orders of depth should work:
        self.assertEqual(self.v.interface_crossed(9000, 10500), True)
        self.assertEqual(self.v.interface_crossed(10500, 9500), True)
        # even small steps across the interface should be detected
        z1, z2 = np.nextafter(10000, 9000), np.nextafter(10000, 11000)
        if z1 != z2:
            self.assertEqual(self.v.interface_crossed(z1, z2), True)
            self.assertEqual(self.v.interface_crossed(z1, z2), self.v.interface_crossed(z2, z1))
        else:
            print("Assertion skipped due to float precision")
        # top and bottom of the model also count as interfaces
        self.assertEqual(self.v.interface_crossed(-0.1, 2), True)
        self.assertEqual(self.v.interface_crossed(10999, 11001), True)

    def test_multiple_interfaces_crossed_returns_true(self):
        """As long as at least one interface between the depths, the function
        should return True"""
        self.assertEqual(self.v.interface_crossed(9000, 12000), True)
        self.assertEqual(self.v.interface_crossed(36000, 9), True)
        self.assertEqual(self.v.interface_crossed(-1, 1000000), True)

    def test_no_interfaces_crossed_returns_false(self):
        self.assertEqual(self.v.interface_crossed(1, 2), False)
        z1 = np.nextafter(10000, 11000)
        z2 = np.nextafter(z1, 11000)
        self.assertEqual(self.v.interface_crossed(z1, z2), False)

    def test_interface_velocities_evaluation(self):
        """Test if method returns velocity from above and below the closest
        interface to the depth."""
        for z in (9999, 10000, 10001):
            with self.subTest(z=z):
                v_top, v_bottom = self.v.interface_velocities(z)
                self.assertEqual(v_top, 5800)
                self.assertEqual(v_bottom, 5000)

    def test_interface_velocites_top(self):
        """For depths above the top layers mid point, the interface between the
        first and the second layer (from the top) should be returned, since the
        interface between the first layer and the air is not defined."""
        v_top, v_bottom = self.v.interface_velocities(10)
        self.assertEqual(v_top, 5800)
        self.assertEqual(v_bottom, 5000)

    def test_interface_velocities_bottom(self):
        """For depths below the bottom layers mid point, the interface between
        it and the layer above should be used."""
        v_top, v_bottom = self.v.interface_velocities(11999)
        self.assertEqual(v_top, 6000)
        self.assertEqual(v_bottom, 7000)


class TestVelocityModel3DLayerIndex(unittest.TestCase):

    def setUp(self) -> None:
        self.layers = [(0, 100, 1800, 4), (100, 200, 2400, 0),
                       (200, 300, 2400, 1), (300, 400, 2700, 0),
                       (400, 500, 2250, 1.5)]
        self.interface_depths = [0, 100, 200, 300, 400, 500]
        self.num_layers = len(self.layers)
        self.vm = VelocityModel3D(self.layers)

    @staticmethod
    def error_msg(depth: float, index_expected: int, index_got: int):
        return f"Got wrong index for depth {depth}. Expected " \
            f"{index_expected}, got {index_got}"

    def testIndexInLayers(self):
        depths = [d+50 for d in self.interface_depths[:-1]]
        for i, depth in zip(range(self.num_layers), depths):
            with self.subTest(i=i):
                index = self.vm.layer_index(depth)
                self.assertEqual(index, i, msg=self.error_msg(depth, i, index))

    def testIndexOnInterface(self):
        """Upper border of a layer is inclusive, so the index of the lower layer
        should be returned."""
        depths = self.interface_depths[:-1]
        for i, depth in zip(range(self.num_layers-1), depths):
            with self.subTest(i=i):
                index = self.vm.layer_index(depth)
                self.assertEqual(index, i, msg=self.error_msg(depth, i, index))

    def testSpecialCaseBottomOfModel(self):
        depth = self.interface_depths[-1]
        index = self.vm.layer_index(depth)
        self.assertEqual(index, self.num_layers-1,
                         msg=self.error_msg(depth, self.num_layers-1, index))


class TestVelocityModel3DNumberOfInterfaces(unittest.TestCase):

    def setUp(self) -> None:
        self.layers = [(0, 100, 1800, 4), (100, 200, 2400, 0),
                       (200, 300, 2400, 1), (300, 400, 2700, 0),
                       (400, 500, 2250, 1.5)]
        self.vm = VelocityModel3D(self.layers)
        # combination of two points and expected result
        self.points = [((0, 0, 0), (0, 0, 500), 4),
                       ((0, 0, 50), (0, 0, 450), 4),
                       ((0, 0, 50), (0, 0, 150), 1),
                       ((0, 0, 200), (0, 0, 250), 0)]
        # transform coordinates to numpy arrays
        self.points = [(np.array(a), np.array(b), num) for a, b, num in self.points]

    def generate_error_msg(self, a: np.ndarray, b: np.ndarray, expected: int,
                           actual: int) -> str:
        return f"Error using points {a}, {b}. Expected number of layers " \
            f"{expected}, got {actual}!"

    def test_correct_values(self):
        for a, b, correct in self.points:
            with self.subTest(a=a, b=b, correct=correct):
                n = self.vm.num_of_interfaces_between(a, b)
                self.assertEqual(n, correct, msg=self.generate_error_msg(a, b, correct, n))

    def test_reverse(self):
        """
        Reversing upper and lower points should not change the result
        """
        for a, b, correct in self.points:
            with self.subTest(a=a, b=b, correct=correct):
                n = self.vm.num_of_interfaces_between(b, a)
                self.assertEqual(n, correct, msg=self.generate_error_msg(a, b, correct, n))

