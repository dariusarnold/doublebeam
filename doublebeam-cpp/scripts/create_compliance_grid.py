from pathlib import Path
import subprocess

import numpy as np
import matplotlib.pyplot as plt

from common import read_last_result


PROGRAM_PATH = Path("/home/darius/git/doublebeam/doublebeam-cpp/cmake-build-release/doublebeam")
OPTIONS_PATH = Path("/home/darius/git/doublebeam/doublebeam-cpp/options_thomas.txt")

# create grid from target_x0 to target_x1 in x direction and from target_y0 to target_y1 in
# y direction with num_x*num_y points.
target_x0 = 0
target_x1 = 1200
target_y0 = 0
target_y1 = 1200
num_x = 20
num_y = 20

x_values = np.linspace(target_x0, target_x1, num_x)
y_values = np.linspace(target_y0, target_y1, num_y)
max_result = np.zeros((num_x, num_y), np.complex)

for ix, x in enumerate(x_values):
    for iy, y in enumerate(y_values):
        subprocess.call(args=(str(PROGRAM_PATH), str(OPTIONS_PATH), f"--target.x={x}", f"--target.y={y}"))
        _, result = read_last_result(PROGRAM_PATH.parent)
        max_result[ix, iy] = np.max(result)

np.save("result_grid", max_result)
im = plt.imshow(np.abs(max_result).T, extent=(target_x0, target_x1, target_y0, target_y1), origin="lower")
plt.xlabel("x axis (m)")
plt.ylabel("y axis (m)")
cb = plt.colorbar(im, format="%.2e")
cb.set_label(r"$|\sigma|$")
plt.title(r"Inverted compliance field ($|max(\sigma)|$)")
plt.show()