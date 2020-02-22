# /usr/bin/env python3

from pathlib import Path
import sys


def main(shotdata_folder: Path, amplitude_factor: float) -> None:
    """
    Read all seismograms from source folders in shotdata folder
    and write them back with their amplitude multiplied by a constant
    factor.
    """
    # iterate over all source folders
    for entry in shotdata_folder.iterdir():
        if not entry.is_dir():
            continue
        print(f"Converting {entry}")
        # iterate over all receiver files
        for filepath in entry.iterdir():
            if not filepath.is_file():
                continue
            if filepath.suffix == ".bin":
                # skip binary files
                continue
            # read file content
            with open(filepath) as file:
                data = file.readlines()
                newdata = []
            for line in data:
                if line.startswith("#"):
                    newdata.append(line)
                    continue
                # flip sign of amplitude
                t, x = [float(s) for s in line.split()]
                newdata.append(f"{t:.3f} {x*amplitude_factor}\n")
            with open(filepath, "w") as file:
                file.writelines(newdata)


help_str = """
Multiply all seismogram amplitudes by -1 to fix a bug in Born forward modeling.

Argument:
path: Path to folder containing the shotdata directory.
    Positional argument 1.
    
amplitude factor: float value. Old amplitude will be multiplied with this factor 
    before seismogram will be saved.
    Positional argument 2.  
"""

if __name__ == '__main__':
    if len(sys.argv) != 3:
        exit(help_str)
    try:
        input_path: Path = Path(sys.argv[1])
    except Exception as e:
        exit(f"Invalid path:\n{str(e)}")
    try:
        amplitude_factor = float(sys.argv[2])
    except ValueError:
        exit(f"Argument amplitude factor is not valid: {sys.argv[2]}")
    input_path = input_path / "shotdata"
    if not input_path.exists():
        exit(f"Path does not exit: {input_path}")
    if not input_path.is_dir():
        exit(f"Did not find shotdata folder at path {input_path}")
    main(input_path, amplitude_factor)
