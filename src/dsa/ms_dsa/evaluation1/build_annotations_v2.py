"""
This script annotates the v2 MS_DSA table (dsa.ms_dsa.ms_dsa_v2: MS_DSA, MS_DSA_shared_charges and the
channel flags) with the lysine properties, and writes the annotated table of every run.

The lysine property table is loaded, and only extended if the v2 table has lysines it does not contain.
"""
from dsa.ms_dsa.config import ms_data_relativePath
from dsa.ms_dsa.ms_dsa_v2 import MS_DSA_Builder
from dsa.ms_dsa.evaluation1.lysine_properties import LysinePropertyTable
from dsa.ms_dsa.evaluation1.ms_dsa_annotator import MSDSAAnnotator


def main(ms_data_relative_path=ms_data_relativePath, n_workers: int = 24, rebuild_lysine_table: bool = False):
    # rebuild_lysine_table: recompute every lysine (needed after adding property columns), instead of only
    # the lysines missing from the saved table
    builder = MS_DSA_Builder(ms_data_relative_path)
    print(f"Annotating {builder.ms_dsa_path}")
    annotator = MSDSAAnnotator(builder.ms_dsa)

    if rebuild_lysine_table:
        lysine_table = LysinePropertyTable.build(annotator.lysines_of_interest(), n_workers=n_workers)
        lysine_table.save()
        annotator.save_peptide_lysines(lysine_table)
        print(f"Rebuilt the lysine property table: {len(lysine_table)} lysines")
    else:
        lysine_table = LysinePropertyTable.load()
        n_before = len(lysine_table)
        lysine_table = lysine_table.extend(annotator.lysines_of_interest(), n_workers=n_workers)
        if len(lysine_table) > n_before:
            lysine_table.save()
            annotator.save_peptide_lysines(lysine_table)
            print(f"Added {len(lysine_table) - n_before} lysines to the lysine property table")

    runs = annotator.annotate(lysine_table)
    annotator.save_runs(runs)
    print(f"Saved annotated runs {list(runs)}")


if __name__ == "__main__":
    import sys
    main(rebuild_lysine_table="--rebuild" in sys.argv)
