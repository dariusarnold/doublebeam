import unittest
from math import sqrt, radians, degrees

import numpy as np

from doublebeam.utils import horizontal_distance, angle, angle_clockwise, length, unit_vector


class TestHorizontalDistance(unittest.TestCase):

    def test_function(self):
        def error_msg(p1, p2, calculated, expected):
            return (f"Wrong horizontal distance between points {p1}, {p2}. Got "
                f"{calculated}, expected {expected}")

        # data consisting of two test points and their distance
        data = [((0, 0, 0), (0, 0, 0), 0),
                ((1, 0, 0), (0, 0, 0), 1),
                ((0, 0, 0), (0, -1, 2), 1),
                ((0, 0, 1), (0, 0, -1), 0),
                ((-1, -1, -42.5), (1, 1, 2), sqrt(2)*2)]
        # convert points to numpy arrays
        data = [(np.array(p1), np.array(p2), d) for p1, p2, d in data]
        for p1, p2, distance in data:
            with self.subTest(p1=p1, p2=p2, distance=distance):
                calculated = horizontal_distance(p1, p2)
                self.assertEqual(calculated, distance,
                                 msg=error_msg(p1, p2, calculated, distance))


class TestAngle(unittest.TestCase):

    def test_90_degrees(self):
        a = np.array((1, 0, 0))
        b = np.array((0, 1, 0))
        self.assertEqual(angle(a, b), radians(90))
        self.assertEqual(angle(a, b, acute=False), radians(90))

    def test_parallel(self):
        a = np.array((0.8, 0.2, 0.4))
        self.assertEqual(angle(a, a), 0)
        self.assertEqual(angle(a, a, acute=False), radians(180))

    def test_float_clipping(self):
        # these values create an acos argument >1 that needs to be clipped
        a = np.array((0.874469283050132, 0.262553720250597, 0.477968688795641))
        self.assertEqual(angle(a, a), 0)
        self.assertEqual(angle(a, a, acute=False), radians(180))

    def test_acute_obtuse(self):
        a = np.array((1, 0, 0))
        b = np.array((1, 1, 0))
        self.assertAlmostEqual(angle(a, b), radians(45))
        self.assertAlmostEqual(angle(a, b, acute=False), radians(135))

    def test_order_of_input_vectors_acute(self):
        a = np.array((2, 1, 0))
        b = np.array((1, 0, 0))
        angle1_acute = angle(a, b)
        angle2_acute = angle(b, a)
        self.assertEqual(angle1_acute, angle2_acute)

    def test_order_of_input_vectors_obtuse(self):
        a = np.array((2, 1, 0))
        b = np.array((1, 0, 0))
        angle1_obtuse = angle(a, b, acute=False)
        angle2_obtuse = angle(b, a, acute=False)
        self.assertEqual(angle1_obtuse, angle2_obtuse)


class TestClockwiseAngle(unittest.TestCase):

    @staticmethod
    def _msg(x: float, angle: float) -> str:
        return f"Error for {degrees(angle)}°. Got {degrees(x)}°"

    def setUp(self) -> None:
        self.x_axis = np.array((1, 0, 0))
        self.test_data = {radians(0): np.array((1, 0, 0)),
                          radians(45): np.array((1, 1, 0)),
                          radians(90): np.array((0, 1, 0)),
                          radians(135): np.array((-1, 1, 0)),
                          radians(180): np.array((-1, 0, 0)),
                          radians(225): np.array((-1, -1, 0)),
                          radians(270): np.array((0, -1, 0)),
                          radians(315): np.array((1, -1, 0))}

    def test_normal(self):
        """
        Test if function calculates angle between x axis and given vector
        correctly by checking against manually computed test data
        """
        for angle, vector in self.test_data.items():
            with self.subTest(angle=angle, vector=vector):
                x = angle_clockwise(self.x_axis, vector)
                self.assertEqual(x, angle, msg=self._msg(x, angle))

    def test_swapped(self):
        """
        If the vectors are swapped, the "other angle" has to be returned. The
        value of this other angle is 360° - first_angle.
        """
        for angle, vector in self.test_data.items():
            with self.subTest(angle=angle, vector=vector):
                x = angle_clockwise(vector, self.x_axis)
                # modulo 360° since for 0° the function will not return 360°
                expected = (radians(360) - angle) % (2*np.pi)
                self.assertEqual(x, expected, msg=self._msg(x, expected))


class TestUnitVector(unittest.TestCase):

    def test_normal_case(self):
        vectors = np.array([(1, 2, 3),
                            (4, 5, 6),
                            (-1, -1, 0),
                            (.1, 0, -.2),
                            (1, 0, 0),
                            (0, 1, 0),
                            (0, 0, -1)])
        for vector in vectors:
            with self.subTest(vector=vector):
                u = unit_vector(vector)
                self.assertAlmostEqual(length(u), 1, places=7,
                                 msg=f"Length != 1: Got {length(u)}")
                self.assertAlmostEqual(angle(vector, u), 0, places=7,
                                 msg=f"Angle not 0°: {degrees(angle(vector, u))}°")
