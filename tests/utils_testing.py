import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Tuple, TextIO

import numpy as np


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


class TestCase(unittest.TestCase):
    """
    Custom TestCase that delegates array asserts to the corresponding numpy
    functions.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.addTypeEqualityFunc(np.ndarray, self._array_equal)

    @staticmethod
    def _array_equal(x, y, msg=None):
        # unittest passes None for empty message, but numpy expects empty string
        msg = "" if msg is None else msg
        np.testing.assert_array_equal(x, y, err_msg=msg)

    def assertAllClose(self, actual: np.ndarray, desired: np.ndarray, rtol: float = 1E-7,
                       atol: float = 0, msg: str = None):
        # unittest passes None for empty message, but numpy expects empty string
        msg = "" if msg is None else msg
        np.testing.assert_allclose(actual, desired, rtol, atol, err_msg=msg)