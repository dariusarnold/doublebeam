import unittest
from math import radians, sin

import numpy as np

from doublebeam.core.models import VelocityModel3D
from twopoint_manual import TwoPointRayTracing


class TestTwoPointRayTracing(unittest.TestCase):

    def setUp(self) -> None:
        self.vm = VelocityModel3D.from_file("/home/darius/git/double-beam/fang2019model.txt")
        self.twopoint = TwoPointRayTracing(self.vm)

    @unittest.skip
    def test_theta_30_degrees(self):
        """Test example ray from fig. 9. Slowness can be calculated from the
        angle theta and the velocity using eq. 13"""
        p = self.twopoint.trace((469., 0., 0), (469/2, 0., 500.))
        #p = self.twopoint.trace((469 / 2, 0., 500.), (469., 0., 0))
        p_correct = sin(radians(30)) / 3000
        self.assertAlmostEqual(p, p_correct, places=5)

    @unittest.skip
    def test_theta_50_degree(self):
        """Test example ray from fig. 9. Slowness cam be calcilated from the
        angle theta and the velocity using eq. 13"""
        p = self.twopoint.trace((434, 0, 500), (868, 0, 0))
        p_correct = sin(radians(50)) / 3000
        self.assertAlmostEqual(p, p_correct, places=5)


class TestNewTwoPointRayTracing(unittest.TestCase):

    def setUp(self) -> None:
        self.vm = VelocityModel3D.from_file("/home/darius/git/double-beam/fang2019model.txt")
        self.rt = TwoPointRayTracing(self.vm)
        self.sources = [np.array((469/2, 0, 500)),
                        np.array((868/2, 0, 500)),
                        np.array((2159/2, 0, 500))]
        self.receivers = [np.array((469, 0, 0)),
                          np.array((868, 0, 0)),
                          np.array((2159, 0, 0))]

    def test_correct_values(self):
        """Test if two point ray tracing gives the same slowness as calculated
        from fig. 9 in Fang2019"""
        slownesses = [self.rt.trace(source, receiver) for source, receiver in zip(self.sources, self.receivers)]
        slownesses_expected = [sin(radians(i)) / 3000 for i in (30, 50, 85)]
        for p, p_expected in zip(slownesses, slownesses_expected):
            with self.subTest(p=p, p_expected=p_expected):
                self.assertAlmostEqual(p, p_expected, places=5)

    def test_source_receiver_swap(self):
        """Swapping receiver and source should give the same result"""
        slownesses = [self.rt.trace(receiver, source) for source, receiver in zip(self.sources, self.receivers)]
        slownesses_expected = [sin(radians(i)) / 3000 for i in (30, 50, 85)]
        for p, p_expected in zip(slownesses, slownesses_expected):
            with self.subTest(p=p, p_expected=p_expected):
                self.assertAlmostEqual(p, p_expected, places=5)
