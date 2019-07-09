"""
Code for reading data from disk, for example waveform data
"""
import ast
import re
from pathlib import Path
from typing import Union, Tuple, List

import numpy as np


def load_wfdata_binary(filename: Path) -> (np.ndarray, np.ndarray):
    return np.load(filename)


def load_wfdata(filename: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load waveform data from file
    :param filename: Path of file containing waveform data
    :return: Tuple of timesteps and amplitudes
    """
    # TODO rename: load_timeseries?
    times, amplitudes = np.loadtxt(filename, unpack=True)
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
    # TODO rename: read_receiverfile
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


def read_header_from_file(filepath: Path) -> List[np.ndarray]:
    """
    Parse header from seismogram file and return source and receiver position.
    :param filepath: Full filename
    :return: Tuple of (source_position, receiver_position)
    """
    with open(filepath, "r") as f:
        # read first 2 lines
        header: List[str] = [next(f).rstrip("\n") for _ in range(2)]
    positions = []
    for line in header:
        # throw away source/receiver, only keep digits
        line = line.split(":")[-1]
        # remove outer square brackets
        line = re.sub("[\[\]]", "", line)
        vec = np.fromstring(line, dtype=float, sep=" ")
        positions.append(vec)
    return positions


def load_seismograms(seismograms_path: Path, seismogram_filename_template: str)\
        -> Tuple[np.ndarray, np.ndarray, List[np.ndarray], np.ndarray]:
    """
    Load seismograms from the given path.
    This loads all seismograms from the path and returns them as well as
    additional information.
    :param seismograms_path: Path to seismogram files
    :param seismogram_filename_template: String where the numeric part of the
    filename is replaced by a star. Eg. if the seismograms where saved as
    receiver_001.txt, receiver_002.txt, ..., pass "receiver_*.txt".
    :return Numpy array consisting of all loaded seismograms, numpy array
    containing the common timesteps for all seismograms, a list of receiver
    positions for these seismograms, and the source position.
    """
    seismograms = []
    receiver_positions = []
    seismogram_filenames = seismograms_path.glob(seismogram_filename_template)
    for fname in sorted(seismogram_filenames):
        source_pos, receiver_pos = read_header_from_file(fname)
        receiver_positions.append(receiver_pos)
        time, seismic_data = load_wfdata(fname)
        seismograms.append(seismic_data)
    return np.array(seismograms), time, receiver_positions, source_pos