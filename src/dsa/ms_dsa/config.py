"""Paths to bundled data and the project tree."""

from __future__ import annotations

from importlib.resources import as_file, files
from pathlib import Path


from dsa.config import PROJECT_ROOT, _SRC, _PKG, _root



MS_DATA_DIR = _SRC / "dsa" / "data"/ "ms_data"

ms_data_relativePath = Path("report.parquet")
