import unittest
from math import radians, sin

from doublebeam.core.models import VelocityModel3D
from doublebeam.core.raytracing.twopoint import TwoPointRayTracing


class TestTwoPointRayTracing(unittest.TestCase):

    def setUp(self) -> None:
        self.vm = VelocityModel3D.from_file("/home/darius/git/double-beam/fang2019model.txt")
        self.twopoint = TwoPointRayTracing(self.vm)

    def test_theta_30_degrees(self):
        """Test example ray from fig. 9. Slowness can be calculated from the
        angle theta and the velocity using eq. 13"""
        p = self.twopoint.trace((469., 0., 0), (469/2, 0., 500.))
        #p = self.twopoint.trace((469 / 2, 0., 500.), (469., 0., 0))
        p_correct = sin(radians(30)) / 3000
        self.assertAlmostEqual(p, p_correct, places=5)

    def test_theta_50_degree(self):
        """Test example ray from fig. 9. Slowness cam be calcilated from the
        angle theta and the velocity using eq. 13"""
        p = self.twopoint.trace((434, 0, 500), (868, 0, 0))
        p_correct = sin(radians(50)) / 3000
        self.assertAlmostEqual(p, p_correct, places=5)
