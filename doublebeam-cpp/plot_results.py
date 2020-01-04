import re
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import argparse


def to_complex(s: str) -> complex:
    return complex(*[float(x) for x in s[1:-1].split(",")])


def plot_scattering_coefficient(data: np.ndarray, min_spacing: float, max_spacing: float,
                                target_id: int, target_x: float, target_y: float, fname):
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
    radii = np.linspace(min_spacing, max_spacing, data.shape[0] + 1)
    im = ax.pcolormesh(angles, radii, data)
    # add grid
    major_ticks_radius = np.linspace(min_spacing, max_spacing, 5)
    ax.set_rticks(major_ticks_radius)
    major_ticks_angle = np.linspace(0, np.pi, 5)
    ax.set_xticks(major_ticks_angle)
    ax.grid(color="white")
    # limit to half circle
    ax.set_thetamax(180)
    # create inner "cutout" by setting origin and min/max for radial axis
    ax.set_rorigin(10)
    ax.set_ylim(min_spacing, max_spacing)
    # TODO add axis labels
    cbar = fig.colorbar(im, ax=ax, shrink=.5, pad=.08, aspect=15, format="%.1E")
    cbar.set_label(r"$|\sigma|$")
    ticks = list(cbar.get_ticks())
    # cbar.set_ticks([np.min(data), np.max(data)] + ticks)
    title = ax.set_title(f"Target {target_id}: x = {target_x} m, y = {target_y} m")
    title.set_position((.5, .85))
    #plt.show()
    #plt.savefig("{str(fname).split('.')[0]}.pdf", bbox_inches="tight")
    plt.savefig(f"{str(fname).split('.')[0]}.png", bbox_inches="tight")


def parse_file(filename: Path) -> np.ndarray:
    with open(filename) as f:
        lines = f.readlines()
    data = []
    for line in lines:
        if line.startswith("#"):
            continue
        if line == "\n":
            # stop parsing on empty line
            break
        row = [to_complex(x) for x in line.split()]
        data.append(row)
    return np.array(data, dtype=np.complex128)


def find_last_result(dir : Path) -> Path:
    """
    Find last result.txt in directory and return path to it.
    """
    # use a natural ordering for strings
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split(r"([0-9]+)", str(key))]
    dir_content = sorted(dir.glob("result*.txt"), key=alphanum_key)
    return dir_content[-1]


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
    data = parse_file(fname)
    plot_scattering_coefficient(np.abs(data), 100, 300, 0, 500, 500, fname)


if __name__ == '__main__':
    main()
