import unittest
from typing import TextIO, Tuple
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np

from doublebeam.io.input import load_wfdata
from doublebeam.core.models import VelocityModel1D, VelocityLayer


def TempFile(content: str) -> Tuple[Path, TextIO]:
    """
    Create a temporary file with the text given in content and return it.
    The file will be kept on disk as long as variable f is alive.
    :param content: Will be written to temp file
    :return: Path to the temporary file, file object
    """
    f = NamedTemporaryFile(mode="w")
    f.write(content)
    f.flush()
    return Path(f.name), f



class TestLoadWfData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.wf_file_path, cls.file = TempFile("""1.7403104935329927E-009   1.0000000195414812E-024
        3.4806209870659854E-009  -1.9548053418115215E-023
        5.2209314805989781E-009  -1.4871093191296816E-023)""")
    @classmethod
    def tearDownClass(cls):
        cls.file.close()

    def test_file_load_success(self):
        """
        Test if loading waveform data from file works correctly
        :return:
        """
        t = np.array([1.7403104935329927E-009, 3.4806209870659854E-009, 5.2209314805989781E-009])
        x = np.array([1.0000000195414812E-024, -1.9548053418115215E-023, -1.4871093191296816E-023])
        lt, lx = load_wfdata(self.wf_file_path)
        np.testing.assert_array_equal(t, lt)
        np.testing.assert_array_equal(x, lx)


class TestLoadVelocityModel1D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        content ="""#Depth(km), vp (km/s), vs (km/s), density (g/cm³)
                    0.00,5.8000,3.3600,2.700
                    10.00,5.8000,3.3600,2.750
                    35.00,8.0400,4.4700,2.800"""
        cls.path, cls.f = TempFile(content)
        cls.layers = np.array([(0, 5.8, 3.36, 2.7), (10., 5.8, 3.36, 2.75), (35., 8.04, 4.47, 2.8)], dtype=VelocityLayer)

    @classmethod
    def tearDownClass(cls):
        cls.f.close()

    def setUp(self):
        self.v = VelocityModel1D(self.layers)

    def test_file_load_success(self):
        """
        Test if loading from a file populates the data structure correctly
        """
        np.testing.assert_array_equal(self.v.layers, self.layers)

    def test_prop_access_success(self):
        """
        Test successfull access to the models properties
        """
        self.assertEqual(self.v.eval_at(5, "r"), 2.7)
        self.assertEqual(self.v.eval_at(555, "s"), 4.47)

    def test_prop_access_out_of_range(self):
        """
        Test for correct exception when wrong depth is used
        """
        with self.assertRaises(LookupError):
            self.v.eval_at(-555, "p")