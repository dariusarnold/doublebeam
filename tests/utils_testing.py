from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Tuple, TextIO


def TempFile(content: str) -> Tuple[Path, TextIO]:
    """
    Create a temporary file with the text given in content and return it.
    The file will be kept on disk as long as variable f is alive.
    :param content: Will be written to temp file
    :return: Path to the temporary file, file object
    """
    f = NamedTemporaryFile(mode="w")
    f.write(content)
    f.flush()
    return Path(f.name), f
