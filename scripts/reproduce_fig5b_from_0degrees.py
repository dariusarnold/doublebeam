"""Try to reproduce fig. 5b from Hu2018a using the data generated"""

from pathlib import Path

import numpy as np

from doublebeam.io.input import read_sources, load_seismograms
from doublebeam.plotting import plot_seismogram_gather


def find_source() -> int:
    """
    Find the source and return its folder index
    """
    sources = read_sources("/home/darius/daten/masterarbeit/sources.txt")
    wanted_x_position = 11200
    source_index = np.argwhere(sources.T[0]==wanted_x_position)
    # +1 since numpy indexing starts with 0 while the folder indexing starts
    # with 1
    return int(source_index + 1)

def find_receivers():
    """
    Find the receivers of the bottom left line
    """
    source_index = str(find_source()).zfill(3)
    seismograms, time, receiver_positions, source_position = load_seismograms(Path(f"/home/darius/daten/masterarbeit/output/source_{source_index}/"),
                     "receiver_*.txt")
    return seismograms, time, receiver_positions, source_position

def plot_data():
    seismograms, time, receiver_positions, source_position = find_receivers()
    # A line is 40 shots long in the simulation but we also need to reverse
    # the order of the traces so the first trace is the one on the bottom
    line = seismograms[:40][::-1]
    plot_seismogram_gather(line)


if __name__ == '__main__':
    plot_data()