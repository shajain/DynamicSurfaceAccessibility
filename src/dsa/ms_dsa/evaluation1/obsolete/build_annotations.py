"""
This script builds (or extends) the lysine property table for the lysines of the MS data,
and writes the peptide x lysine table and the annotated MS_DSA table of every run.
"""
from dsa.ms_dsa.config import ms_data_relativePath
from dsa.ms_dsa.evaluation1.config import LYSINE_PROPERTIES_FILE
from dsa.ms_dsa.evaluation1.lysine_properties import LysinePropertyTable
from dsa.ms_dsa.evaluation1.ms_dsa_annotator import MSDSAAnnotator


def main(ms_data_relative_path=ms_data_relativePath, rebuild: bool = False, n_workers: int = 24):
    annotator = MSDSAAnnotator.from_ms_file(ms_data_relative_path)
    lysines = annotator.lysines_of_interest()
    print(f"{sum(len(r) for r in lysines.values())} lysines of interest in {len(lysines)} uniprot ids")

    if LYSINE_PROPERTIES_FILE.exists() and not rebuild:
        lysine_table = LysinePropertyTable.load().extend(lysines, n_workers=n_workers)
    else:
        lysine_table = LysinePropertyTable.build(lysines, n_workers=n_workers)
    lysine_table.save()
    print(f"Saved {len(lysine_table)} lysines to {LYSINE_PROPERTIES_FILE}")

    annotator.save_peptide_lysines(lysine_table)
    runs = annotator.annotate(lysine_table)
    annotator.save_runs(runs)
    print(f"Saved annotated runs {list(runs)}")


if __name__ == "__main__":
    main()
