"""
Conversion between text and binary format for input data
"""

from pathlib import Path

import numpy as np

from doublebeam.io.input import load_wfdata


def wfdata_text_to_binary(input_filename: Path, output_filename=None):
    """
    Convert data to binary file format
    :param input_filename: Filepath of input file
    :param output_filename: Filepath for output file. If not given, use
    input_filename
    :return:
    """
    t, x = load_wfdata(input_filename)
    if output_filename is None:
        output_filename = input_filename
    np.save(output_filename.as_posix(), np.array((t, x)))
