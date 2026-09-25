import importlib
import os
from typing import Literal

import numpy as np
import pandas as pd

from dsa.dips.config import INTERFACE_DATA_DIR, PROJECT_ROOT
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from pypdb import get_all_info
import requests
from dsa.dips.mapping.dips_to_uniprot import map_dips_to_uniprot





class BindingInterface:
    def __init__(self, angstrom_threshold: float=8.0):
        interface_data_dir = INTERFACE_DATA_DIR[angstrom_threshold]
        df_train = pd.read_csv(os.path.join(interface_data_dir, "dips_train.csv"))
        df_test = pd.read_csv(os.path.join(interface_data_dir, "dips_val.csv"))
        self.df_ppi = pd.concat([df_train, df_test], ignore_index=True)
        self.df_ppi["pdb_id"] = self.df_ppi["pair_key"].apply(self.extract_pdb_id)
        #self.df_sequence_mapping = self._build_sequence_mapping_df()
        self.uniprot_mapping = map_dips_to_uniprot(self.df_ppi)
        # self.aligner = PairwiseAligner()
        # self.aligner.mode = "global"
        # self.aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        # self.aligner.open_gap_score = -10.0
        # self.aligner.extend_gap_score = -0.5

    # def _build_sequence_mapping_df(self):
    #     sequence_to_indices = {}

    #     for row_idx, row in enumerate(self.df_ppi.itertuples(index=False)):
    #         seq_a = row.seq_a
    #         pdb_id = self.extract_pdb_id(row.pair_key)
    #         if isinstance(seq_a, str) and seq_a:
    #             if seq_a not in sequence_to_indices:
    #                 sequence_to_indices[seq_a] = {"seq_a": [], "seq_b": [], 'pdb_id': [pdb_id]}
    #             sequence_to_indices[seq_a]["seq_a"].append(row_idx)
    #             sequence_to_indices[seq_a]["pdb_id"].append(pdb_id)

    #         seq_b = row.seq_b
    #         if isinstance(seq_b, str) and seq_b:
    #             if seq_b not in sequence_to_indices:
    #                 sequence_to_indices[seq_b] = {"seq_a": [], "seq_b": [], 'pdb_id': [pdb_id]}
    #             sequence_to_indices[seq_b]["seq_b"].append(row_idx)
    #             sequence_to_indices[seq_b]["pdb_id"].append(pdb_id)    
    #     rows = [
    #         {
    #             "sequence": sequence,
    #             "seq_a": chain_indices["seq_a"],
    #             "seq_b": chain_indices["seq_b"],
    #             "pdb_id": chain_indices["pdb_id"],
    #         }
    #         for sequence, chain_indices in sequence_to_indices.items()
    #     ]
    #     df_sequence_mapping = pd.DataFrame(rows)
    #     df_pdb_sequences = pd.DataFrame(
    #         {
    #             "pdb_id": [
    #                 self.extract_pdb_id(row.pair_key)
    #                 for row in self.df_ppi.itertuples(index=False)
    #             ],
    #             "seq_a": self.df_ppi["seq_a"],
    #             "seq_b": self.df_ppi["seq_b"],
    #         }
    #     )
    #     self.uniprot_mapping = map_dips_to_uniprot(df_pdb_sequences)
    #     return df_sequence_mapping

    def extract_pdb_id(self, pdb_file_name: str):
        return pdb_file_name.split('_', 1)[1].split('.', 1)[0]


  
    def get_interface_residues(self, pair_index: int):
        row = self.df_ppi.iloc[pair_index]
        seq_a = row["seq_a"]
        seq_b = row["seq_b"]
        cmap = np.load(PROJECT_ROOT / row["contact_map"])


if __name__ == "__main__":
    BI = BindingInterface()
    df1 = BI.df_ppi[BI.df_ppi["label"]==1]
    df0 = BI.df_ppi[BI.df_ppi["label"]==0]
    print(df1.columns)
    print(df1.head()["n_contacts"])
    print(len(df1["seq_a"].iloc[0]))
    print(len(df1["seq_b"].iloc[0]))
    print(df1["seq_a"].nunique())
    print(df1["seq_b"].nunique())
    print(df1.shape)
    print(df1.len_a.max())
    print(df1.len_b.max())
    print(df1.len_a.mean())
    print(df1.len_b.mean())
   


  # def _sequence_similarity(self, query, sequence):
    #     if not query or not sequence:
    #         return 0.0
    #     alignments = self.aligner.align(query, sequence)
    #     if len(alignments) == 0:
    #         return 0.0
    #     best_alignment = alignments[0]

    #     matches = 0
    #     alignment_length = 0
    #     coordinates = best_alignment.coordinates
    #     for i in range(coordinates.shape[1] - 1):
    #         q_start, q_end = coordinates[0, i], coordinates[0, i + 1]
    #         s_start, s_end = coordinates[1, i], coordinates[1, i + 1]
    #         q_span = q_end - q_start
    #         s_span = s_end - s_start

    #         if q_span > 0 and s_span > 0:
    #             span = min(q_span, s_span)
    #             q_block = query[q_start : q_start + span]
    #             s_block = sequence[s_start : s_start + span]
    #             matches += sum(aa == bb for aa, bb in zip(q_block, s_block))
    #             alignment_length += max(q_span, s_span)
    #         else:
    #             alignment_length += max(q_span, s_span)

    #     if alignment_length == 0:
    #         return 0.0
    #     return matches / alignment_length

  

    # def find_closest_sequence(
    #     self,
    #     query: str,
    #     *,
    #     chain: Literal["a", "b", "both"] = "both",
    # ) -> tuple[pd.Series, float, Literal["seq_a", "seq_b"]]:
    #     """Return the dataframe row whose sequence best matches ``query``.

    #     Similarity is a normalized global alignment score in [0, 1] (higher is closer),
    #     using a gap-aware protein sequence matcher.
    #     Searches ``seq_a`` and/or ``seq_b`` depending on ``chain``.

    #     Returns
    #     -------
    #     row : pd.Series
    #         Full PPI row for the best match.
    #     score : float
    #         SequenceMatcher ratio for the winning (sequence, query) pair.
    #     which : {"seq_a", "seq_b"}
    #         Which column matched best (when ``chain=="both"``, the higher score wins).
    #     """
    #     q = query.strip().upper()
    #     if not q:
    #         raise ValueError("query sequence is empty")

    #     best_score = -1.0
    #     best_idx = None
    #     best_col = None

    #     cols = []
    #     if chain in ("a", "both"):
    #         cols.append("seq_a")
    #     if chain in ("b", "both"):
    #         cols.append("seq_b")
    #     if not cols:
    #         raise ValueError("chain must be 'a', 'b', or 'both'")

    #     for row in self.df_sequence_index.itertuples(index=False):
    #         sequence = row.sequence
    #         if not isinstance(sequence, str) or not sequence:
    #             continue

    #         score = self._sequence_similarity(q, sequence.upper())
    #         for col in cols:
    #             indices = getattr(row, col)
    #             if not isinstance(indices, list) or not indices:
    #                 continue
    #             if score > best_score:
    #                 best_score = score
    #                 best_idx = indices[0]
    #                 best_col = col

    #     if best_idx is None or best_col is None:
    #         raise ValueError("no comparable sequences in dataframe")

    #     return self.df_ppi.iloc[best_idx], best_score, best_col
    # def get_gene_name(uniprot_id):
    #     url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    #     data = requests.get(url).json()
    #     genes = data.get("genes", [])
    #     return genes[0].get("geneName", {}).get("value", "") if genes else ""
    #     # returns protein name, organism, UniProt IDs, etc

