"""Parquet I/O helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

_PACKAGE_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _PACKAGE_DIR.parent.parent


def default_parquet_path() -> Path:
    """Path to ``report.parquet`` in the project root."""
    return _PROJECT_ROOT / "report.parquet"


def read_report_parquet(path: Optional[str | Path] = None) -> pd.DataFrame:
    """Load the report table from Parquet.

    Parameters
    ----------
    path
        File to read. Defaults to :func:`default_parquet_path`.
    """
    p = Path(path) if path is not None else default_parquet_path()
    return pd.read_parquet(p)
