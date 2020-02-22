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
