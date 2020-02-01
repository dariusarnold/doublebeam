from pathlib import Path
import subprocess

import numpy as np
import matplotlib.pyplot as plt

from common import read_last_result, parse_file


PROGRAM_PATH = Path("/home/darius/git/doublebeam/doublebeam-cpp/cmake-build-release/doublebeam")
OPTIONS_PATH = Path("/home/darius/git/doublebeam/doublebeam-cpp/options_thomas.txt")

# create grid from target_x0 to target_x1 in x direction and from target_y0 to target_y1 in
# y direction with num_x*num_y points.
target_x0 = 200
target_x1 = 1000
target_y0 = 200
target_y1 = 1000
num_x = 20
num_y = 20

x_values = np.linspace(target_x0, target_x1, num_x)
y_values = np.linspace(target_y0, target_y1, num_y)
max_result_x_parallel = np.zeros((num_x, num_y), np.complex)
max_result_y_parallel = np.zeros((num_x, num_y), np.complex)


def get_max_y_parallel_result(data: np.array):
    """
    Get max value from center of array.
    """
    return np.max(np.abs(data[:, 15:45]))


def get_max_x_parallel_result(data: np.array):
    """
    Get max value from left and right side of array.
    """
    return max(np.max(np.abs(data[:, 0:15])), np.max(np.abs(data[:, 45:])))


for ix, x in enumerate(x_values):
    for iy, y in enumerate(y_values):
        subprocess.call(args=(str(PROGRAM_PATH), str(OPTIONS_PATH), f"--target.x={x}", f"--target.y={y}"))
        _, result = read_last_result(PROGRAM_PATH.parent)
        _, result = parse_file(Path(f"/home/darius/git/doublebeam/doublebeam-cpp/grid/result{ix * num_x + iy}.txt"))
        # first axis of result contains num spacings (41), second axis contains orientation
        max_result_x_parallel[ix, iy] = get_max_x_parallel_result(result)
        max_result_y_parallel[ix, iy] = get_max_y_parallel_result(result)


np.save("result_grid_x", max_result_x_parallel)
np.save("result_grid_y", max_result_y_parallel)

im = plt.imshow(np.abs(max_result_x_parallel).T, extent=(target_x0, target_x1, target_y0, target_y1), origin="lower")
plt.xlabel("x axis (m)")
plt.ylabel("y axis (m)")
cb = plt.colorbar(im, format="%.2e")
cb.set_label(r"$|\sigma|$")
plt.title(r"Inverted compliance field for fractures parallel x axis ($|max(\sigma)|$)")
plt.show()

im = plt.imshow(np.abs(max_result_y_parallel).T, extent=(target_x0, target_x1, target_y0, target_y1), origin="lower")
plt.xlabel("x axis (m)")
plt.ylabel("y axis (m)")
cb = plt.colorbar(im, format="%.2e")
cb.set_label(r"$|\sigma|$")
plt.title(r"Inverted compliance field for fractures parallel y axis ($|max(\sigma)|$)")
plt.show()