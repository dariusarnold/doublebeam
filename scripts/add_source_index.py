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