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
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.interpolate import interp2d

from common import Options, find_last_result, parse_file


def plot_scattering_coefficient(data: np.ndarray, options: Options, fname):
    """
    :param data: (N, M) array that contains scattering coefficient sigma. First
    axis gives angle samples, second axis gives fracture spacing.
    :param min_spacing: min fracture spacing
    :param max_spacing: max fracture spacing
    :param fname: Filename for output file, pdf and png will be saved.
    :return:
    """
    fig, ax = plt.subplots(subplot_kw={"polar": True}, figsize=(6, 6), dpi=210)
    # +1 otherwise last column of data will be ignored
    angles = np.radians(np.linspace(0, 180, data.shape[1] + 1))
    radii = np.linspace(options.fracture_params.spacing_min, options.fracture_params.spacing_max,
                        data.shape[0] + 1)
    # vmin = 0 to force full color scale from 0 to max
    im = ax.pcolormesh(angles, radii, data, rasterized=True)#, cmap="jet")
    # set radial ticks to fracture spacings
    num_of_ticks = 5
    major_ticks_radius = np.linspace(options.fracture_params.spacing_min,
                                     options.fracture_params.spacing_max, num_of_ticks)
    ax.set_rticks(major_ticks_radius)
    # set angular ticks to 0°-180°
    major_ticks_angle = np.linspace(0, np.pi, num_of_ticks)
    ax.set_xticks(major_ticks_angle)
    # add fracture spacing label to radial axis
    rlabel_pos = ax.get_rlabel_position()
    ax.text(np.radians(rlabel_pos - 40), ax.get_rmax() / 1.45, "Fracture spacing (m)", rotation=0,
            ha="center",
            va="center")
    # enable grid
    ax.grid(color="white")
    # limit to half circle
    ax.set_thetamax(180)
    # create inner "cutout" by setting origin and min/max for radial axis
    ax.set_rorigin(10)
    ax.set_ylim(options.fracture_params.spacing_min, options.fracture_params.spacing_max)
    # force scientific notation with exponent at top
    cbar_format = ticker.ScalarFormatter()
    cbar_format.set_powerlimits((-1, 1))
    # add colorbar and label it
    cbar = fig.colorbar(im, ax=ax, shrink=.5, pad=.08, aspect=15, format=cbar_format)
    cbar.set_label(r"$|\sigma|$")
    # plot values used to generate data in a box on the figure
    textbox_content = "\n".join((fr"$\omega_s = {options.beam_params.source_frequency}$ Hz",
                                 fr"$\omega_r = {options.beam_params.reference_frequency}$ Hz",
                                 fr"$w = {options.beam_params.width} $ m",
                                 fr"window $= {options.beam_params.window_length}$ s",
                                 fr"max stack. dist. $= {options.beam_params.max_stacking_distance:.0f}$ m",
                                 f"{options.source_beam_centers.num_x}x{options.source_beam_centers.num_y} source beam centers:",
                                 fr"    $x = {options.source_beam_centers.x0:.0f}$ ... ${options.source_beam_centers.x1:.0f}$ m",
                                 fr"    $y = {options.source_beam_centers.y0:.0f}$ ... ${options.source_beam_centers.y1:.0f}$ m"
                                 ))
    box_properties = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax.text(0., 0.15, textbox_content, transform=ax.transAxes, fontsize=10, verticalalignment="top",
            bbox=box_properties)
    # add path to data as text above the plot
    path_textbox_content = "\n".join((f"{options.data.model_path}",
                                      f"{options.data.data_path}"))
    ax.text(-0.1, 1, path_textbox_content, transform=ax.transAxes, fontsize=7)
    # set title to target position
    title = ax.set_title(f"Target ({options.target.x} m, {options.target.y} m, {options.target.z} m)")
    title.set_position((.5, .85))
    # plt.show()
    plt.savefig(f"{str(fname).split('.')[0]}.pdf", bbox_inches="tight")
    plt.savefig(f"{str(fname).split('.')[0]}.png", bbox_inches="tight")


def main():
    if len(sys.argv) == 1:
        # no arguments, open file picker
        import tkinter.filedialog, tkinter
        tkinter.Tk().withdraw()
        fname = Path(tkinter.filedialog.askopenfilename())
    else:
        # pathname as argument
        fname = Path(sys.argv[1])
    if fname.is_dir():
        # if directory, find last result and plot that
        fname = find_last_result(fname)
    print(f"Plotting {fname}")
    options, data = parse_file(fname)

    # only use abs value since that is whats plotted by Zheng2013 too
    data = np.abs(data)

    interpolate = False
    scale_factor = 5
    if interpolate:
        y = np.arange(0, data.shape[0])
        x = np.arange(0, data.shape[1])
        interpolator = interp2d(x, y, data)
        x_dense = np.linspace(0, x[-1], int(len(x)*scale_factor))
        y_dense = np.linspace(0, y[-1], int(len(y)*scale_factor))
        data = interpolator(x_dense, y_dense)
        
    plot_scattering_coefficient(data, options, fname)


if __name__ == '__main__':
    main()
