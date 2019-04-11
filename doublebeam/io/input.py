"""
Code for reading data from disk, for example waveform data
"""

from pathlib import Path

import numpy as np


def load_wfdata_binary(filename: Path) -> (np.ndarray, np.ndarray):
    return np.load(filename)


def load_wfdata(filename: Path) -> (np.ndarray, np.ndarray):
    jumbled = np.fromfile(str(filename), sep=" ")
    # reading text files with fromfile is faster than genfromtxt and loadtxt.
    # But all values are read sequentially line by line into one array, so
    # we split it after reading
    times = jumbled[::2]
    amplitudes = jumbled[1::2]
    # TODO own class for wf data that handles conversion to pressure
    return times, amplitudes


def get_file_iterator(path, extension):
    """
    Return an iterator that yields all files in the directory given by path
    which have a certain file extension
    :param path: Path to folder in which to search for files. If a relative
    path is given, it will interpreted from the current working directory.
    :type path: pathlib.Path
    :param extension: Extension to look for without leading dot
    :type extension: str
    :return:
    :rtype:
    """
    if not path.is_absolute():
        path = Path.cwd() / path
    for file in path.glob(f"*.{extension}"):
        yield file
