import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np

from doublebeam.io.input import load_wfdata
from doublebeam.core.models import VelocityModel1D


def TempFile(content):
    """
    Create a temporary file with the text given in content and return it.
    The file will be kept on disk as long as variable f is alive.
    :param content: Will be written to temp file
    :return: Path to the temporary file
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
