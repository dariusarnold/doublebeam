/*
 * Copyright (C) 2019-2020  Darius Arnold
 *
 * This file is part of doublebeam.
 *
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
"""
This script reads the receiver files of the simulation and fixes the wrong
header format:
# source: [[1.12e+04]
#  [5.60e+03]
#  [1.00e+01]]
# receiver: [5216.6665 7708.3335    0.    ]
becomes
# source: [   0. 5200.   10.]
# receiver: [3200. 5600.    0.]

"""

import re
from pathlib import Path

import numpy as np

path = Path("/home/darius/daten/masterarbeit/output")

source_folders = sorted(list(path.iterdir()))
for source_folder in source_folders:
    print(f"Rewriting {str(source_folder)}")
    receiver_files = sorted(list(source_folder.iterdir()))
    for receiver_file in receiver_files:
        with receiver_file.open("r") as f:
            content = f.readlines()
        # failsafe: if already converted (and # not missing), skip file
        if content[0].startswith("#") and content[1].startswith("#") \
                and not content[2].startswith("#"):
            continue
        # else convert
        source = "".join(content[0:3])
        receiver = content[3]
        rest = content[4:]

        source = re.findall(r"\[+(.*?)\]+", source)
        receiver = re.findall(r"\[+(.*)\]+", receiver)[0].split()

        source = np.array([float(x) for x in source])
        receiver = np.array([float(x) for x in receiver])

        with receiver_file.open("w") as f:
            f.write(f"# source: {source}\n# receiver: {receiver}\n")
            f.writelines(rest)
