from dsa.config import DATA_DIR

EVALUATION1_DIR = DATA_DIR / "ms_dsa" / "evaluation1"
# lookup table with one row per uniprot lysine, written by LysinePropertyTable
LYSINE_PROPERTIES_FILE = EVALUATION1_DIR / "lysine_properties.parquet"
# one row per (peptide, lysine) of the ms data, written by MSDSAAnnotator
PEPTIDE_LYSINES_FILE = EVALUATION1_DIR / "peptide_lysines.parquet"
# one annotated ms_dsa table per run, written by MSDSAAnnotator
ANNOTATED_RUNS_DIR = EVALUATION1_DIR / "runs"
