import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

from common import Options, find_last_result, parse_file


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
    cbar_format = ticker.ScalarFormatter()
    # force scientific notation with exponent at top
    cbar_format.set_powerlimits((-1, 1))
    cbar = fig.colorbar(im, ax=ax, shrink=.5, pad=.08, aspect=15, format=cbar_format)
    cbar.set_label(r"$|\sigma|$")
    textbox_content = "\n".join((fr"$\omega_s = {options.beam_params.source_frequency}$ Hz",
                                 fr"$\omega_r = {options.beam_params.reference_frequency}$ Hz",
                                 fr"$w = {options.beam_params.width} $ m",
                                 fr"window $= {options.beam_params.window_length}$ s",
                                 fr"max stack. dist. $= {options.beam_params.max_stacking_distance:.0f}$ m",
                                 f"{options.source_beam_centers.num_x}x{options.source_beam_centers.num_y} source beam centers:",
                                 fr"    $x = {options.source_beam_centers.x0:.0f}$ ... ${options.source_beam_centers.x1:.0f}$ m",
                                 fr"    $y = {options.source_beam_centers.y0:.0f}$ ... ${options.source_beam_centers.y1:.0f}$ m"
                                 ))
    path_textbox_content = "\n".join((f"{options.data.model_path}",
                                      f"{options.data.data_path}"))
    box_properties = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax.text(0., 0.15, textbox_content, transform=ax.transAxes, fontsize=10, verticalalignment="top",
            bbox=box_properties)
    ax.text(-0.1, 1, path_textbox_content, transform=ax.transAxes, fontsize=7)
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
        fname = Path(tkinter.filedialog.askopenfilename())
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
