"""Small I/O helpers shared by the pipeline modules."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TextIO


def open_text(path: Path) -> TextIO:
    """Open a text file for reading, transparently handling ``.gz`` files."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()
