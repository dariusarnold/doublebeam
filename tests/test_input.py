import unittest

import numpy as np

from utils import TempFile
from doublebeam.io.input import load_wfdata


class TestLoadWfData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # create temp file with example waveform data
        content = """1.7403104935329927E-009   1.0000000195414812E-024
                     3.4806209870659854E-009  -1.9548053418115215E-023
                     5.2209314805989781E-009  -1.4871093191296816E-023"""
        cls.wf_file_path, cls.file = TempFile(content)

    def test_file_load_success(self):
        """
        Test if loading waveform data from file works correctly
        """
        t = np.array([1.7403104935329927E-009, 3.4806209870659854E-009, 5.2209314805989781E-009])
        x = np.array([1.0000000195414812E-024, -1.9548053418115215E-023, -1.4871093191296816E-023])
        lt, lx = load_wfdata(self.wf_file_path)
        np.testing.assert_array_equal(t, lt)
        np.testing.assert_array_equal(x, lx)


