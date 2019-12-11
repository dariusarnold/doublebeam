import sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def to_complex(s: str) -> complex:
    return complex(*[float(x) for x in s[1:-1].split(",")])


def plot_scattering_coefficient(data: np.ndarray, min_spacing: float, max_spacing: float,
                                target_id: int, target_x: float, target_y: float):
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
    plt.show()
    #plt.savefig("test.pdf", bbox_inches="tight")
    #plt.savefig("test.png", bbox_inches="tight")


def parse_file(filename: Path) -> np.ndarray:
    with open(filename) as f:
        lines = f.readlines()
    data = []
    for line in lines:
        if line.startswith("#"):
            continue
        row = [to_complex(x) for x in line.split()]
        data.append(row)
    return np.array(data, dtype=np.complex128)


def main():
    if len(sys.argv) == 1:
        sys.exit("Give path to result file as positional argument.")
    data = parse_file(sys.argv[1])
    plot_scattering_coefficient(np.abs(data), 100, 300, 0, 0, 0)


if __name__ == '__main__':
    main()
