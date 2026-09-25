from importlib.resources import as_file, files
from pathlib import Path


# Anchor on the package via importlib.resources (string name avoids import cycles).
_root = files("dsa")
try:
    _PKG = Path(_root).resolve()
except TypeError:
    with as_file(_root) as p:
        _PKG = Path(p).resolve()

_SRC = _PKG.parent
PROJECT_ROOT = _SRC.parent

DATA_DIR = PROJECT_ROOT / "src" / "dsa" / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
