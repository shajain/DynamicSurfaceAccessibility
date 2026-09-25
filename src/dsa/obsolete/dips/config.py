"""Runtime configuration: data paths and environment overrides."""

from __future__ import annotations

import os
from pathlib import Path

from dsa.config import PROJECT_ROOT




#DIPS utput files
INTERFACE_DATA_DIR = {}
sequence_mapping_file = {}
INTERFACE_DATA_DIR[8.0] = PROJECT_ROOT / "src" / "dsa" / "data" / "dips"/ "dips_output_calpha_8A" / "stage1"
INTERFACE_DATA_DIR[4.5] = PROJECT_ROOT / "src" / "dsa" / "data" / "dips"/ "dips_output_heavy_45A" / "stage1"
sequence_mapping_file[8.0] = PROJECT_ROOT / "src" / "dsa" / "data" / "dips"/ "sequence_mapping_8A.csv"
sequence_mapping_file[4.5] = PROJECT_ROOT / "src" / "dsa" / "data" / "dips"/ "sequence_mapping_45A.csv"