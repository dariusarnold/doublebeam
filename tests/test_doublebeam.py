import unittest
from math import radians

import numpy as np

from doublebeam.utils import angle, pairwise, angle_clockwise, generate_vector_arc


class TestPhiHatGeneration(unittest.TestCase):

    def setUp(self) -> None:
        self.number_of_vectors = (12, 64)
        self.directions = np.array((1, 0, 0)), np.array((1, 1, 0))
        self.vectors = [generate_vector_arc(num, dir) for num, dir in
                        zip(self.number_of_vectors, self.directions)]

    def test_lengths(self):
        """Test if all vectors are unit vectors (length 1)"""
        lengths = [np.linalg.norm(vectors, axis=1) for vectors in self.vectors]
        for length in lengths:
            with self.subTest():
                np.testing.assert_allclose(length, np.ones_like(length))

    def test_number(self):
        """Test if the specified number of vectors was returned"""
        for vector, number in zip(self.vectors, self.number_of_vectors):
            with self.subTest():
                self.assertEqual(len(vector), number)

    def test_angles(self):
        """Test if all angles are equal to the expected angle you get when you
        divide 180° degrees into n parts."""
        expected_angles = [radians(180 / (num - 1)) for num in self.number_of_vectors]
        angles_got = [np.array([angle(v1, v2) for v1, v2 in pairwise(vectors)])
                      for vectors in self.vectors]
        for angles_got_, expected_angle_ in zip(angles_got, expected_angles):
            with self.subTest():
                np.testing.assert_allclose(angles_got_, np.full_like(angles_got_, expected_angle_))

    def test_orientation(self):
        """Test if the first vector is 90° degrees to the left from direction
        and the last vector 90° to the right."""
        for direction, vectors in zip(self.directions, self.vectors):
            with self.subTest("Testing if first vector is 90° anticlockwise"):
                self.assertAlmostEqual(angle_clockwise(direction, vectors[0]), radians(270), places=14)
            with self.subTest("Testing is last vector has angle 180° to axis"):
                self.assertAlmostEqual(angle(direction, vectors[-1]), radians(90), places=14)
