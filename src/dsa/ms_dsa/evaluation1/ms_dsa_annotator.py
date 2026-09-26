"""
Annotates the MS_DSA peptides of every run with the properties of their lysines,
looked up from a LysinePropertyTable.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from dsa.ms_dsa.ms_dsa import MS_DSA_Builder
from dsa.ms_dsa.config import ms_data_relativePath
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
from dsa.ms_dsa.evaluation1.lysine_properties import LysinePropertyTable, CONTACT_CLASSES, CONTACT_TYPES
from dsa.ms_dsa.evaluation1.config import PEPTIDE_LYSINES_FILE, ANNOTATED_RUNS_DIR

PEPTIDE_KEY = ["Uniprot", "peptide_seq"]
STATISTICS = {"avg": "mean", "min": "min", "max": "max"}
# statistics of the numeric lysine properties over the lysines of a peptide, written as {column}_{statistic};
# numeric properties not listed here are only in the peptide x lysine table
NUMERIC_AGGREGATIONS = {
    "asa_nz":               ["avg", "min", "max"],   # AlphaFold monomer
    "rel_asa":              ["avg", "min", "max"],
    "asa_nz_multimer_mean": ["avg", "min", "max"],   # mean over the chain copies of each lysine
    "asa_nz_multimer_min":  ["min"],                 # most buried chain copy of the most buried lysine
    "asa_nz_multimer_max":  ["max"],                 # most exposed chain copy of the most exposed lysine
    "nz_bfactor_z_mean":    ["avg", "min", "max"],
    "nz_bfactor_z_min":     ["min"],                 # lowest NZ B-factor over chains and lysines
    "nz_bfactor_z_max":     ["max"],                 # highest NZ B-factor over chains and lysines
    "avg_bfactor_z_mean":   ["avg", "min", "max"],
}
# contact-class asa_nz: aggregated over only the lysines of the peptide that have a value (the lysines in that
# class), with max of the per-lysine max, min of the per-lysine min, and average of the per-lysine mean
CLASS_ASA_PREFIXES = [f"asa_nz_pdb_{contact_type}_{contact_class}"
                      for contact_type in CONTACT_TYPES for contact_class in CONTACT_CLASSES] + ["asa_nz_pdb_no_contact"]
AVAILABLE_LYSINE_AGGREGATIONS = {f"{prefix}_{statistic}": [name]
                                 for prefix in CLASS_ASA_PREFIXES
                                 for statistic, name in (("max", "max"), ("min", "min"), ("mean", "avg"))}
# boolean properties that hold for a peptide only if they hold for all its lysines
ALL_LYSINE_FLAGS = ["no_contact_pdb", "no_contact_af"]


class MSDSAAnnotator:
    """
    Maps every (Uniprot, peptide) of the MS_DSA table to the uniprot positions of its lysines,
    and aggregates the lysine properties over the lysines of each peptide:
        boolean properties -> True if any lysine is True (NA if all are NA); ALL_LYSINE_FLAGS: True if all lysines
                              are True (NA if some lysine is NA and none is False)
        numeric properties -> {column}_{statistic}, as listed in NUMERIC_AGGREGATIONS and
                              AVAILABLE_LYSINE_AGGREGATIONS
        purely_intra_{pdb,af} -> an intra-chain class for the peptide and no lysine with pdb_inter_any_bond

    Rows with several uniprot ids in Protein.Ids are split into one row per uniprot id.
    A peptide is located by exact match of its first occurrence in the uniprot sequence.

    Usage:
        annotator = MSDSAAnnotator.from_ms_file()
        table = LysinePropertyTable.build(annotator.lysines_of_interest())
        runs = annotator.annotate(table)       # {run_index: DataFrame}
    """

    def __init__(self, ms_dsa: pd.DataFrame, require_all_lysines: bool = True):
        # require_all_lysines: numeric aggregates are NA unless every lysine of the peptide has a value
        self.require_all_lysines = require_all_lysines
        self.ms_dsa = split_multi_uniprot_id_rows(ms_dsa)
        self.peptide_lysines = self.map_peptides_to_lysines()

    @classmethod
    def from_ms_file(cls, ms_data_relative_path: Path = ms_data_relativePath, **kwargs) -> "MSDSAAnnotator":
        return cls(MS_DSA_Builder(ms_data_relative_path).ms_dsa, **kwargs)

    def map_peptides_to_lysines(self) -> pd.DataFrame:
        """One row per (Uniprot, peptide_seq, lysine); peptides not found in the sequence have no rows."""
        rows = []
        for uniprot_id, peptide_seq in self.ms_dsa[PEPTIDE_KEY].drop_duplicates().itertuples(index=False):
            start = peptide_start_in_uniprot(peptide_seq, uniprot_id)
            if start is None:
                continue
            rows.extend({"Uniprot": uniprot_id, "peptide_seq": peptide_seq, "peptide_start_index": start,
                         "position_in_peptide": i + 1, "unp_resnum": start + i}
                        for i, aa in enumerate(peptide_seq) if aa == "K")
        peptide_lysines = pd.DataFrame(rows, columns=PEPTIDE_KEY + ["peptide_start_index", "position_in_peptide", "unp_resnum"])
        peptide_lysines["unp_resnum"] = peptide_lysines["unp_resnum"].astype("Int64")
        return peptide_lysines

    def lysines_of_interest(self) -> dict[str, list[int]]:
        return {uniprot_id: sorted(resnums.unique().tolist())
                for uniprot_id, resnums in self.peptide_lysines.groupby("Uniprot")["unp_resnum"]}

    def peptide_lysine_properties(self, lysine_table: LysinePropertyTable) -> pd.DataFrame:
        """The peptide x lysine table with the properties of each lysine."""
        return self.peptide_lysines.merge(lysine_table.table, how="left",
                                          left_on=["Uniprot", "unp_resnum"], right_on=LysinePropertyTable.KEY
                                          ).drop(columns="uniprot_id")

    def peptide_properties(self, lysine_table: LysinePropertyTable) -> pd.DataFrame:
        """One row per (Uniprot, peptide_seq) with the lysine properties aggregated over its lysines."""
        lysines = self.peptide_lysine_properties(lysine_table)
        grouped = lysines.groupby(PEPTIDE_KEY, sort=False)
        columns = {"peptide_start_index": grouped["peptide_start_index"].first(),
                   "unp_lysines": grouped["unp_resnum"].agg(lambda s: s.astype(int).tolist()),
                   "n_lysines": grouped.size()}
        keys = [lysines[k] for k in PEPTIDE_KEY]
        for column in lysine_table.property_columns:
            if lysines[column].dtype == "boolean":
                as_float = lysines[column].astype("Float64").groupby(keys, sort=False)
                if column in ALL_LYSINE_FLAGS:
                    # False if any lysine is False, True if all lysines are True, NA otherwise
                    lowest = as_float.min()
                    has_missing = lysines[column].isna().groupby(keys, sort=False).any()
                    columns[column] = (lowest > 0).astype("boolean").mask(lowest.isna() | ((lowest > 0) & has_missing))
                else:
                    # True if any lysine is True, False if all known lysines are False, NA if no lysine is known
                    columns[column] = (as_float.max() > 0).astype("boolean")
        for column, statistics in NUMERIC_AGGREGATIONS.items():
            has_missing = lysines[column].isna().groupby(keys, sort=False).any()
            for name in statistics:
                values = grouped[column].agg(STATISTICS[name])
                columns[f"{column}_{name}"] = values.mask(has_missing) if self.require_all_lysines else values
        for column, statistics in AVAILABLE_LYSINE_AGGREGATIONS.items():
            if column in lysines:
                for name in statistics:
                    columns[f"{column}_{name}"] = grouped[column].agg(STATISTICS[name])
        peptides = pd.DataFrame(columns)
        self.add_purely_intra(peptides)
        return peptides.reset_index()

    @staticmethod
    def add_purely_intra(peptides: pd.DataFrame):
        # from the peptide flags: a lysine with an intra-chain class and no lysine with an inter-chain bond
        if "pdb_inter_any_bond" not in peptides:
            return
        inter = peptides["pdb_inter_any_bond"]
        for source in ("pdb", "af"):
            intra = peptides[[f"{source}_intra_{contact_class}" for contact_class in CONTACT_CLASSES]]
            known = intra.notna().all(axis=1) & inter.notna()
            peptides[f"purely_intra_{source}"] = (intra.fillna(False).any(axis=1) & ~inter.fillna(False)).astype("boolean").mask(~known)

    def annotate(self, lysine_table: LysinePropertyTable) -> dict[int, pd.DataFrame]:
        """One annotated MS_DSA table per run; peptide_mapped is False for peptides not found in the sequence."""
        annotated = self.ms_dsa.merge(self.peptide_properties(lysine_table), how="left", on=PEPTIDE_KEY)
        annotated["peptide_mapped"] = annotated["n_lysines"].notna()
        annotated["n_lysines"] = annotated["n_lysines"].fillna(0).astype(int)
        return {run_index: run_df.reset_index(drop=True) for run_index, run_df in annotated.groupby("Run.Index")}

    def save_peptide_lysines(self, lysine_table: LysinePropertyTable, path: Path = PEPTIDE_LYSINES_FILE):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.peptide_lysine_properties(lysine_table).to_parquet(path, index=False)

    @staticmethod
    def save_runs(runs: dict[int, pd.DataFrame], directory: Path = ANNOTATED_RUNS_DIR):
        directory.mkdir(parents=True, exist_ok=True)
        for run_index, run_df in runs.items():
            run_df.to_parquet(directory / f"run_{run_index}.parquet", index=False)

    @staticmethod
    def load_run(run_index: int, directory: Path = ANNOTATED_RUNS_DIR) -> pd.DataFrame:
        return pd.read_parquet(directory / f"run_{run_index}.parquet")


def split_multi_uniprot_id_rows(ms_dsa: pd.DataFrame) -> pd.DataFrame:
    ms_dsa = ms_dsa.copy()
    ms_dsa["Uniprot"] = ms_dsa["Protein.Ids"].str.split(";")
    ms_dsa = ms_dsa.explode("Uniprot", ignore_index=True)
    ms_dsa["Uniprot"] = ms_dsa["Uniprot"].str.strip().str.upper()
    return ms_dsa


def peptide_start_in_uniprot(peptide_seq: str, uniprot_id: str) -> int|None:
    # 1-based position of the first exact match of the peptide in the uniprot sequence
    sequence = UniprotToSequence.get_sequence(uniprot_id)
    if not isinstance(peptide_seq, str) or not sequence:
        return None
    index = sequence.find(peptide_seq)
    return index + 1 if index != -1 else None

