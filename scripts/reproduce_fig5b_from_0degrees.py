#
# Copyright (C) 2019-2020  Darius Arnold
#
# This file is part of doublebeam.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
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