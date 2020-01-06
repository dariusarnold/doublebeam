import re
import sys
from pathlib import Path
from typing import List, Tuple
import configparser

import numpy as np
import matplotlib.pyplot as plt


def to_complex(s: str) -> complex:
    return complex(*[float(x) for x in s[1:-1].split(",")])


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

    def __init__(self, width: float, frequency_hz: float, window_length: float,
                 max_stacking_distance: float):
        self.width = width
        self.frequency = frequency_hz
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
    freq = config["beam"].getfloat("frequency")
    window_length = config["beam"].getfloat("window_length")
    max_stacking_distance = config["beam"].getfloat("max_stacking_distance")
    return BeamParameters(width, freq, window_length, max_stacking_distance)


def parse_file(filename: Path) -> Tuple[Options, np.ndarray]:
    with open(filename) as f:
        data = f.read()
    # after [result] section the file contains the results of the doublebeam algorithm (stacking amplitude sigma)
    result_index = data.index("[result]")
    config = configparser.ConfigParser()
    config.read_string(data[0:result_index])
    options = Options(extract_data(config), extract_target(config),
                      extract_source_beam_center_params(config),
                      extract_fracture_params(config), extract_beam_params(config))

    values = []
    for line in data[result_index + len("[result]\n"):-1].split("\n"):
        row = [to_complex(x) for x in line.split()]
        values.append(row)
    return options, np.array(values, dtype=np.complex128)


def find_last_result(dir: Path) -> Path:
    """
    Find last result.txt in directory and return path to it.
    """
    # use a natural ordering for strings
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split(r"([0-9]+)", str(key))]
    dir_content = sorted(dir.glob("result*.txt"), key=alphanum_key)
    return dir_content[-1]


def plot_scattering_coefficient(data: np.ndarray, options: Options, fname):
    """
    :param data: (N, M) array that contains scattering coefficient sigma. First
    axis gives angle samples, second axis gives fracture spacing.
    :param min_spacing: min fracture spacing
    :param max_spacing: max fracture spacing
    :return:
    """
    fig, ax = plt.subplots(subplot_kw={"polar": True}, figsize=(6, 6), dpi=210)
    angles = np.radians(np.linspace(0, 180, data.shape[1] + 1))
    # +1 otherwise last column of data will be ignored
    radii = np.linspace(options.fracture_params.spacing_min, options.fracture_params.spacing_max,
                        data.shape[0] + 1)
    im = ax.pcolormesh(angles, radii, data)
    # add grid
    major_ticks_radius = np.linspace(options.fracture_params.spacing_min,
                                     options.fracture_params.spacing_max, 5)
    ax.set_rticks(major_ticks_radius)
    major_ticks_angle = np.linspace(0, np.pi, 5)
    ax.set_xticks(major_ticks_angle)
    # add label to radial axis
    rlabel_pos = ax.get_rlabel_position()
    ax.text(np.radians(rlabel_pos - 40), ax.get_rmax() / 1.45, "Fracture spacing (m)", rotation=0,
            ha="center",
            va="center")
    ax.grid(color="white")
    # limit to half circle
    ax.set_thetamax(180)
    # create inner "cutout" by setting origin and min/max for radial axis
    ax.set_rorigin(10)
    ax.set_ylim(options.fracture_params.spacing_min, options.fracture_params.spacing_max)
    cbar = fig.colorbar(im, ax=ax, shrink=.5, pad=.08, aspect=15, format="%.1E")
    cbar.set_label(r"$|\sigma|$")
    textbox_content = "\n".join((fr"$\omega = {options.beam_params.frequency}$ Hz",
                                 fr"$w = {options.beam_params.width} $ m",
                                 fr"window $= {options.beam_params.window_length}$ s",
                                 fr"max stack. dist. $= {options.beam_params.max_stacking_distance:.0f}$ m",
                                 f"{options.source_beam_centers.num_x}x{options.source_beam_centers.num_y} source beam centers:",
                                 fr"    $x = {options.source_beam_centers.x0:.0f}$ ... ${options.source_beam_centers.x1:.0f}$ m",
                                 fr"    $y = {options.source_beam_centers.y0:.0f}$ ... ${options.source_beam_centers.y1:.0f}$ m"))
    box_properties = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax.text(0., 0.15, textbox_content, transform=ax.transAxes, fontsize=10, verticalalignment="top",
            bbox=box_properties)

    ticks = list(cbar.get_ticks())
    # cbar.set_ticks([np.min(data), np.max(data)] + ticks)
    title = ax.set_title(f"Target x = {options.target.x} m, y = {options.target.y} m")
    title.set_position((.5, .85))
    # plt.show()
    # plt.savefig("{str(fname).split('.')[0]}.pdf", bbox_inches="tight")
    plt.savefig(f"{str(fname).split('.')[0]}.png", bbox_inches="tight")


def main():
    if len(sys.argv) == 1:
        # no arguments, open file picker
        import tkinter.filedialog, tkinter
        tkinter.Tk().withdraw()
        fname = tkinter.filedialog.askopenfilename()
    else:
        # pathname as argument
        fname = Path(sys.argv[1])
    if fname.is_dir():
        # if directory, find last result and plot that
        fname = find_last_result(fname)
        print(f"Plotting {fname}")
    options, data = parse_file(fname)
    plot_scattering_coefficient(np.abs(data), options, fname)


if __name__ == '__main__':
    main()
