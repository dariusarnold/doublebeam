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
import sys
from pathlib import Path

"""
Format source file to required format.
Goes from one source per line with 

source_index x y z

to

nsrc=number_of_sources_in_file

source = source_index
xsource = x
ysource = y
zsource = z

"""

def replace_filename(path: Path, new_name: str) -> Path:
    """
    Replace filename from path with new_name.
    """
    new_path = path.with_name(new_name).with_suffix(path.suffix)
    return new_path
    

def main(source_file_path: Path):
    with open(source_file_path) as source_file:
        lines = source_file.readlines()
    number_of_sources = len(lines)
    new_source_file_path = replace_filename(source_file_path, f"{source_file_path.stem}_converted")
    with open(new_source_file_path, "w") as source_file:
        source_file.write(f"nsrc = {number_of_sources}\n")
        for line in lines:
            source_index, x, y, z = line.split()
            source_file.write(f"source = {source_index}\n")
            source_file.write(f"xsource = {x}\n")
            source_file.write(f"ysource = {y}\n")
            source_file.write(f"zsource = {z}\n\n")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Specify source file as positional argument.")
    try:
        source_file_path = Path(sys.argv[1])
    except Exception:
        sys.exit(f"Failed to find path {sys.argv[1]}")
    main(source_file_path)
