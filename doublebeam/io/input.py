"""
Code for reading data from disk, for example waveform data
"""
import ast
import re
from pathlib import Path
from typing import Union, Tuple

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


def read_stations(filepath: Path) -> np.ndarray:
    """
    Read file defining receiver positions.
    :param filepath:
    :return: Array containing all receiver positions with shape (N, 3).
    N is the number of receivers defined in the file (first line), 3 are the
    x, y, z coordinates of the individual receivers.
    """
    stations = np.genfromtxt(filepath, skip_header=1, usecols=(1, 2, 3),
                             dtype=np.float32)
    return stations


def read_sources(filepath: Path) -> np.ndarray:
    """
    Read sources from text file and return (N, 3) array of source coordinates.
    """
    int_regex = r"(\d+)"
    float_regex = r"(?:\d+(?:\.\d*)?|\.\d+)"
    number_of_sources_regex = fr"nsrc\s+=\s+(?P<nsrc>{int_regex})"
    x_regex = fr"xsource\s+=\s+(?P<xsource>{float_regex})"
    y_regex = fr"ysource\s+=\s+(?P<ysource>{float_regex})"
    z_regex = fr"zsource\s+=\s+(?P<zsource>{float_regex})"
    # save name of regex group as key, which can later be used to retrieve the
    # matched value
    regexes = {"nsrc": number_of_sources_regex,
               "xsource": x_regex,
               "ysource": y_regex,
               "zsource": z_regex}

    def match_line(line: str) -> Union[Tuple[str, Union[int, float]],
                                       Tuple[None, None]]:
        """
        Match all regexes against a line and return either the first matching
        one or a tuple(None, None) if no regexes matched the line.
        :return: A tuple of the key to the regex that matched and an interpreted
        value (int or float)
        """
        for key, regex in regexes.items():
            match = re.search(regex, line)
            if match:
                value = ast.literal_eval(match.group(key))
                return key, value
        return None, None

    x_values, y_values, z_values = [], [], []
    nsrc = None
    with open(filepath, "r") as f:
        for line in f:
            key, value = match_line(line)
            if key is None:
                continue
            elif key == "nsrc" and nsrc is None:
                nsrc = value
            elif key == "xsource":
                x_values.append(value)
            elif key == "ysource":
                y_values.append(value)
            elif key == "zsource":
                z_values.append(value)

        sources = np.array([(x, y, z) for x, y, z in
                            zip(x_values, y_values, z_values)], dtype=np.float32)
        if nsrc is None:
            raise ValueError(f"Number of sources (nsrc) not specified in file"
                             f"{filepath}")
        if len(sources) != nsrc:
            raise ValueError(f"Expected {nsrc} stations in file {filepath} but "
                             f"got only {len(sources)}.")
        return sources
