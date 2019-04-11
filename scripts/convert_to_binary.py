#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

try:
    import doublebeam.io.convert
except ImportError as e:
    # hack to import the parent directory where the doublebeam module lives
    # only required if the cwd isn't set to the double-beam folder (e.g. by an IDE)
    parent_directory = Path(sys.path[0]).parent
    sys.path.insert(0, parent_directory.as_posix())
    try:
        import doublebeam.io.convert
    except ImportError:
        raise e


def convert_to_binary(files):
    """
    Convert
    :param files:
    :return:
    """

    for file in files:
        path = Path(file)
        doublebeam.io.convert.wfdata_text_to_binary(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert waveform data to binary format to decrease load time and disk space usage. Binary files will have the same name as their text version with .npy suffix added.")
    parser.add_argument("files", nargs="*", help="List of waveform data in text format. Special symbols like * will be expanded by the shell.")

    args = parser.parse_args()

    convert_to_binary(args.files)
