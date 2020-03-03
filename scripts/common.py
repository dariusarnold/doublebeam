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
import re
from pathlib import Path
from typing import Tuple
import configparser
import numpy as np


def to_complex(s: str) -> complex:
    """
    Convert one entry in results txt to python complex.
    '(0.00181026,-0.000633689)' -> 0.00181026 -0.000633689j
    """
    return complex(*[float(x) for x in s[1:-1].split(",")])


def find_last_result(dir: Path) -> Path:
    """
    Find last result.txt in directory and return path to it.
    """
    # use a natural ordering for strings
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split(r"([0-9]+)", str(key))]
    dir_content = sorted(dir.glob("result*.txt"), key=alphanum_key)
    return dir_content[-1]


class Position:

    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z


class DataParameters:

    def __init__(self, data_path: Path, velocity_model_path: Path):
        self.data_path = data_path
        self.model_path = velocity_model_path


class SourceBeamCenters:
    def __init__(self, x0: float, x1: float, y0: float, y1: float, num_x: int, num_y: int):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.num_x = num_x
        self.num_y = num_y


class FractureParameters:

    def __init__(self, num_orientations: int, spacing_min: float, spacing_max: float, num_spacings):
        self.num_orientations = num_orientations
        self.spacing_min = spacing_min
        self.spacing_max = spacing_max
        self.num_spacings = num_spacings


class BeamParameters:

    def __init__(self, width: float, reference_frequency_hz: float, source_frequency_hz: float,
                 window_length: float, max_stacking_distance: float):
        self.width = width
        self.reference_frequency = reference_frequency_hz
        self.source_frequency = source_frequency_hz
        self.window_length = window_length
        self.max_stacking_distance = max_stacking_distance


class Options:

    def __init__(self, data_params: DataParameters, target: Position,
                 source_beam_centers: SourceBeamCenters,
                 fracture_params: FractureParameters, beam_params: BeamParameters):
        self.data = data_params
        self.target = target
        self.source_beam_centers = source_beam_centers
        self.fracture_params = fracture_params
        self.beam_params = beam_params


def extract_target(config: configparser.ConfigParser) -> Position:
    x = config["target"].getfloat("x")
    y = config["target"].getfloat("y")
    z = config["target"].getfloat("z")
    return Position(x, y, z)


def extract_data(config: configparser.ConfigParser) -> DataParameters:
    data_path = Path(config["data"]["path"])
    model_path = Path(config["data"]["model"])
    return DataParameters(data_path, model_path)


def extract_source_beam_center_params(config: configparser.ConfigParser) -> SourceBeamCenters:
    try:
        sbc_config_part = config["source_beam_centers"]
    except KeyError:
        sbc_config_part = config["source beam centers"]
    x0 = sbc_config_part.getfloat("x0")
    x1 = sbc_config_part.getfloat("x1")
    y0 = sbc_config_part.getfloat("y0")
    y1 = sbc_config_part.getfloat("y1")
    num_x = sbc_config_part.getint("num_x")
    num_y = sbc_config_part.getint("num_y")
    return SourceBeamCenters(x0, x1, y0, y1, num_x, num_y)


def extract_fracture_params(config: configparser.ConfigParser) -> FractureParameters:
    frac_conf_part = config["fractures"]
    num_orientations = frac_conf_part.getint("num_orientations")
    spacing_min = frac_conf_part.getfloat("spacing_min")
    spacing_max = frac_conf_part.getfloat("spacing_max")
    num_spacings = frac_conf_part.getint("num_spacings")
    return FractureParameters(num_orientations, spacing_min, spacing_max, num_spacings)


def extract_beam_params(config: configparser.ConfigParser) -> BeamParameters:
    width = config["beam"].getfloat("width")
    source_freq = config["beam"].getfloat("source_frequency")
    reference_freq = config["beam"].getfloat("reference_frequency")
    window_length = config["beam"].getfloat("window_length")
    max_stacking_distance = config["beam"].getfloat("max_stacking_distance")
    return BeamParameters(width, reference_freq, source_freq, window_length, max_stacking_distance)


def parse_file(filename: Path) -> Tuple[Options, np.ndarray]:
    with open(filename) as f:
        data = f.read()
    result_index = data.index("[result]")
    # header contains parameters used to create the results
    config = configparser.ConfigParser()
    config.read_string(data[0:result_index])
    options = Options(extract_data(config), extract_target(config),
                      extract_source_beam_center_params(config),
                      extract_fracture_params(config), extract_beam_params(config))
    # after [result] section the file contains the results of the doublebeam algorithm (stacking amplitude sigma)
    values = []
    for line in data[result_index + len("[result]\n"):].split("\n"):
        row = [to_complex(x) for x in line.split()]
        if len(row): values.append(row)
    return options, np.array(values, dtype=np.complex128)


def read_last_result(folder: Path) -> Tuple[Options, np.ndarray]:
    """
    Find last result in folder and read it, returning options used and the result data.
    """
    fname = find_last_result(folder)
    return parse_file(fname)
