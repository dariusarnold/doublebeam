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
Add missing source index to sources.txt
"""

def main(source_file_path: Path):
    with open(source_file_path) as source_file:
        lines = source_file.readlines()
    new_lines = []
    source_index = 1
    for line in lines:
        if line != "\n" and line.split()[0].startswith("xsource"):
            new_lines.append(f"source = {source_index}\n")
            source_index += 1
        new_lines.append(line)
    with open(source_file_path, "w") as source_file:
        source_file.writelines(new_lines)



if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Specify source file as positional argument.")
    try:
        source_file_path = Path(sys.argv[1])
    except Exception:
        sys.exit(f"Failed to find path {sys.argv[1]}")
    main(source_file_path)