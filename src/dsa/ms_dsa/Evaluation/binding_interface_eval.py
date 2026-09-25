from Bio.Seq import Seq
from dsa.binding_interfaces.contacts.protein_contacts import ProteinContactsStore
from dsa.ms_dsa.ms_dsa import MS_DSA_Builder
from dsa.ms_dsa.config import ms_data_relativePath
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
import numpy as np
from matplotlib import pyplot as plt
import warnings
from dsa.config import FIGURES_DIR, DATA_DIR
import pandas as pd
from dsa.sasa.compute_sasa import SASAMonomerStore
from matplotlib_venn import venn2
from dsa.binding_interfaces.contacts.bond_characterization.library import BondLibrary
from collections import defaultdict

seq_dsa_dict = {"peptide_seq": [], "MS_DSA": [], "Protein.Names": [], "Genes": [], "Protein.Ids": [], "Protein.Group": []}
class BIEvaluation:
    BOND_TYPES = BondLibrary.bond_types_implemented()
    BOND_METHOD = BondLibrary.default_method
    SASA_MONOMER_STORE = SASAMonomerStore()
    SASA_TYPE = "asa_nz"
    def __init__(self, ms_data_relativePath: str=ms_data_relativePath, distance_cutoff: float=6.5, contact_filter_key: str = "lysine_only", max_uniprot_ids: int=None):
        # self.cutoff = cutoff
        # self.lysine_only = lysine_only
        self.ms_data_relativePath = ms_data_relativePath
        self.ms_dsa = MS_DSA_Builder(self.ms_data_relativePath).ms_dsa
        self.ms_dsa = self.split_multi_uniprot_id_rows(self.ms_dsa)
        self.run_indices = self.ms_dsa["Run.Index"].unique()
        self.ms_dsa_by_run = {run_index: self.ms_dsa[self.ms_dsa["Run.Index"] == run_index] for run_index in self.run_indices}
        self.contacts_pdb = ProteinContactsStore("pdb", distance_cutoff, contact_filter_key)
        self.contacts_afold = ProteinContactsStore("alphafold", distance_cutoff, contact_filter_key)
        self.max_uniprot_ids = max_uniprot_ids
        self.uniprot_ids_processed = set()
        #self.ms_dsa_avg_over_runs = self.build_ms_dsa_avg_over_runs()
        for run_index in self.run_indices:
            self.ms_dsa_by_run[run_index] = self.add_binding_interface_and_sasa_scores_info(self.ms_dsa_by_run[run_index])
        #self.sasa_df = self.build_sasa_df()
        #self.combined_df = self.combine_ms_dsa_and_sasa_df()
    # def combine_ms_dsa_and_sasa_df(self):
    #     return pd.merge(self.sasa_df, self.ms_dsa, left_on=["Sequence"], right_on=["peptide_seq"], how="inner")

    # def build_sasa_df(self):
    #     sasa_df_raw = pd.read_csv(DATA_DIR / "sasa" / "afold.csv")
    #     sasa_df = sasa_df_raw[["Sequence", "Keys", "RelASA", "asa"]].groupby(["Sequence", "Keys"]).agg({"RelASA": "mean", "asa": "mean"}).reset_index()
    #     return sasa_df
    
    # def build_ms_dsa_avg_over_runs(self):

    #     for (seq, run_index), seq_df in self.ms_dsa.groupby(["peptide_seq", "Run.Index"]):


    def add_binding_interface_and_sasa_scores_info(self, ms_dsa: pd.DataFrame):
        ms_dsa["unp_lysines"] = pd.Series([[]] * len(ms_dsa), dtype=object, index=ms_dsa.index)
        # ms_dsa["peptide_start_index"] = None # This column will be added to the ms_dsa dataframe
        # ms_dsa["is_inter_chain_contacts"] = False # This column will be added to the ms_dsa dataframe
        # ms_dsa["is_intra_chain_contacts"] = False # This column will be added to the ms_dsa dataframe
        # ms_dsa["asa_avg"] = np.nan # This column will be added to the ms_dsa dataframe
        # ms_dsa["asa_min"] = np.nan # This column will be added to the ms_dsa dataframe
        # ms_dsa["asa_max"] = np.nan # This column will be added to the ms_dsa dataframe
        # ms_dsa["rel_asa_avg"] = np.nan # This column will be added to the ms_dsa dataframe
        # ms_dsa["rel_asa_min"] = np.nan # This column will be added to the ms_dsa dataframe
        # ms_dsa["rel_asa_max"] = np.nan # This column will be added to the ms_dsa dataframe
        # ms_dsa["asa_max"] = np.nan # This column will be added to the ms_dsa dataframe
        ms_dsa["skipped"] = True # This column will be added to the ms_dsa dataframe
        iter_count = 0
        sasa_mismatch_dict = defaultdict(bool)
        for uniprot_id, unp_df in ms_dsa.groupby("Uniprot"):
            contacts_pdb = self.contacts_pdb.load(uniprot_id)
            contacts_afold = self.contacts_afold.load(uniprot_id)
            if contacts_pdb is None or contacts_afold is None:
                warnings.warn(f"Alphafold or PDB Contacts not found for uniprot_id {uniprot_id}")
                continue
            for index, row in unp_df.iterrows():
                peptide_seq = row["peptide_seq"]
                #uniprot_ids = [id.strip().upper() for id in uniprot_ids.split(";")]
                unp_lysines, start_index = self.peptide_to_lysine_residues_in_uniprot(peptide_seq, uniprot_id)
                ms_dsa.at[index, "peptide_start_index"] = start_index
                if unp_lysines:
                    ms_dsa.loc[[index], "unp_lysines"] = [unp_lysines]
                ms_dsa.at[index, "intra_chain_pdb"] = any(contacts_pdb.does_residue_satisfy(unp_lysines, contact_type="intra_chain"))
                ms_dsa.at[index, "inter_chain_pdb"] = any(contacts_pdb.does_residue_satisfy(unp_lysines, contact_type="inter_chain"))
                ms_dsa.at[index, "intra_chain_alphafold"] = any(contacts_afold.does_residue_satisfy(unp_lysines, contact_type="intra_chain"))
                # ms_dsa.at[index, "intra_chain"] = ms_dsa.at[index, "intra_chain_pdb"] or ms_dsa.at[index, "intra_chain_alphafold"]
                # ms_dsa.at[index, "inter_chain"] = ms_dsa.at[index, "inter_chain_pdb"] 
                # ms_dsa.at[index, "intra_chain_strict
                for bond_type in self.BOND_TYPES:
                    ms_dsa.at[index, f"{bond_type}_pdb"] = any(contacts_pdb.does_residue_satisfy(unp_lysines, bond_type=bond_type, method=self.BOND_METHOD))
                    ms_dsa.at[index, f"{bond_type}_alphafold"] = any(contacts_afold.does_residue_satisfy(unp_lysines, bond_type=bond_type, method=self.BOND_METHOD))
                
                # sasa_monomer = sasa_monomer_store.load(uniprot_id)
                # if sasa_monomer is None:
                #     warnings.warn(f"SASA not found for uniprot_id {uniprot_id}")
                #     continue
                # sasa_sequence_consistency = np.all(np.array(sasa_monomer.get_residue_names(lysine_indices)) == 'K')
                # if not sasa_sequence_consistency:
                #     warnings.warn(f"SASA sequence consistency not found for uniprot_id {uniprot_id}")
                #     continue
                # asa = [sasa_monomer.get_sasa_for_residue(li, type="asa") for li in lysine_indices]
                # rel_asa = [sasa_monomer.get_sasa_for_residue(li, type="rel_asa") for li in lysine_indices]
                #asa_nz = [sasa_monomer.get_sasa_for_residue(li, type="asa_nz") for li in lysine_indices]
                asa_nz = self.get_solvent_accessiblity(unp_lysines, sasa_type=self.SASA_TYPE, uniprot_id=uniprot_id)
                sasa_mismatch = False
                if asa_nz and any(np.isnan(asa_nz)):
                    sasa_mismatch = True
                    sasa_mismatch_dict[uniprot_id] = True
                try:    
                    ms_dsa.at[index, "asa_nz_avg"] = np.nan if sasa_mismatch else np.mean(asa_nz) 
                    ms_dsa.at[index, "asa_nz_min"] = np.nan if sasa_mismatch else np.min(asa_nz) 
                    ms_dsa.at[index, "asa_nz_max"] = np.nan if sasa_mismatch else np.max(asa_nz) 
                except:
                    warnings.warn(f"SASA mismatch for index {index} and uniprot_id {uniprot_id}")
                # is_intra_chain_contact = is_intra_chain_contact_pdb or is_intra_chain_contact_AF
                # ms_dsa.at[index, "is_intra_chain_contact"] = is_intra_chain_contact
                # ms_dsa.at[index, "is_strictly_inter_chain_contact_pdb"] = is_inter_chain_contact_pdb and not is_intra_chain_contact
                # ms_dsa.at[index, "is_strictly_intra_chain_contact_pdb"] = is_intra_chain_contact and not is_inter_chain_contact_pdb
                # ms_dsa.at[index, "peptide_start_index"] = start_indices
                ms_dsa.at[index, "skipped"] = False
            if iter_count % 100 == 0:
                print(f"Processed {iter_count} uniprot_ids")
                keys = list(sasa_mismatch_dict.keys())
                print(f"Number of uniprot_ids withSASA mismatch: {len(sasa_mismatch_dict)}")
                print(f"example mismatch: {keys[0] if len(keys) > 0 else 'None'}")
            iter_count += 1
            self.uniprot_ids_processed.add(uniprot_id)
            if self.max_uniprot_ids is not None and iter_count >= self.max_uniprot_ids:
                break
        return ms_dsa
            
    def get_solvent_accessiblity(self, unp_lysines: list[tuple[int, str]], sasa_type: str, uniprot_id: str):
        sasa_monomer = self.SASA_MONOMER_STORE.load(uniprot_id)
        if sasa_monomer is None:
            warnings.warn(f"SASA not found for uniprot_id {uniprot_id}")
            return None
        return sasa_monomer.get_sasa_for_residue(unp_lysines, type=sasa_type)

    def save_processed_ms_dsa(self):
        [self.ms_dsa_by_run[run_index].to_csv(f"ms_dsa_with_binding_interface_and_sasa_scores_{run_index}.csv", index=False) for run_index in self.run_indices]
        return

    @classmethod
    def load_processed_ms_dsa(cls, run_index: int):
        return pd.read_csv(f"ms_dsa_with_binding_interface_and_sasa_scores_{run_index}.csv")

    def plot_boxplot(self):
        ms_dsa = self.ms_dsa[self.ms_dsa["skipped"] == False]
        print(ms_dsa["is_inter_chain_contact_pdb"].dtype)
        print(ms_dsa["is_inter_chain_contact_pdb"].value_counts(dropna=False))
        inter = ms_dsa["MS_DSA"][ms_dsa["is_inter_chain_contact_pdb"]]
        not_inter = ms_dsa["MS_DSA"][~ms_dsa["is_inter_chain_contact_pdb"]]
        intra = ms_dsa["MS_DSA"][ms_dsa["is_intra_chain_contact_pdb"]]
        not_intra = ms_dsa["MS_DSA"][~ms_dsa["is_intra_chain_contact_pdb"]]
        intra_AF = ms_dsa["MS_DSA"][ms_dsa["is_intra_chain_contact_AF"]]
        not_intra_AF = ms_dsa["MS_DSA"][~ms_dsa["is_intra_chain_contact_AF"]]
        print(inter.count(), not_inter.count(), intra.count(), not_intra.count(), intra_AF.count(), not_intra_AF.count())
        data = [inter.dropna().values, not_inter.dropna().values, intra.dropna().values, not_intra.dropna().values, intra_AF.dropna().values, not_intra_AF.dropna().values]
        labels = ["Inter", "Not Inter", "Intra", "Not Intra", "Intra AF", "Not Intra AF"]
        plt.boxplot(data, labels=labels)
        plt.xlabel("Binding interface type")
        plt.ylabel("MS_DSA")
        plt.title("MS_DSA vs Binding interface type")
        plt.legend()
        plt.show()

    def plot_dsa_vs_sasa_with_binding_annotation(self, sasa_type: str="asa_nz", statistic: str="avg", run_index: int=None):
        ms_dsa = self.ms_dsa_by_run[run_index]
        ms_dsa = ms_dsa[ms_dsa["skipped"] == False]
        sasa_column = f"{sasa_type}_{statistic}"
        inter = ms_dsa[ms_dsa["is_strictly_inter_chain_contact_pdb"]]
        print(f"Inter: {inter.shape}")
        print(f"% of Inter less than 0.9: {len(inter[inter['MS_DSA'] < 0.9]) / len(inter)}")
        #not_inter = self.combined_df[~self.combined_df["is_inter_chain_contacts"]]
        intra_pdb = ms_dsa[ms_dsa["is_intra_chain_contact_pdb"]]
        print(f"Intra: {intra_pdb.shape}")
        print(f"% of Intra less than 0.9: {len(intra_pdb[intra_pdb['MS_DSA'] < 0.9]) / len(intra_pdb)}")
        intra_AF = ms_dsa[ms_dsa["is_intra_chain_contact_AF"]]
        print(f"Intra AF: {intra_AF.shape}")
        print(f"% of Intra AF less than 0.9: {len(intra_AF[intra_AF['MS_DSA'] < 0.9]) / len(intra_AF)}")
        intra = ms_dsa[ms_dsa["is_intra_chain_contact"]]
        print(f"Intra: {intra.shape}")
        print(f"% of Intra less than 0.9: {len(intra[intra['MS_DSA'] < 0.9]) / len(intra)}")
        not_contact = ms_dsa[~ms_dsa["is_inter_chain_contact_pdb"] & ~ms_dsa["is_intra_chain_contact"]]
        print(f"Not Contact: {not_contact.shape}")
        print(f"% of Not Contact less than 0.9: {len(not_contact[not_contact['MS_DSA'] < 0.9]) / len(not_contact)}")
        plt.scatter(not_contact["MS_DSA"].values, not_contact[sasa_column].values, color="black", label="Not Contact", s=10, alpha=0.5)
        # plt.scatter(not_inter["MS_DSA"].dropna().values, not_inter["RelASA"].dropna().values, color="black", label="Not Inter")
        plt.scatter(intra["MS_DSA"].values, intra[sasa_column].values, color="green", label="Intra", s=10, alpha=0.5)
        # plt.scatter(not_intra["MS_DSA"].dropna().values, not_intra["RelASA"].dropna().values, color="black", label="Not Intra")
        plt.scatter(inter["MS_DSA"].values, inter[sasa_column].values, color="red", label="Inter", s=10, alpha=0.5)
        plt.xlabel("MS_DSA")
        plt.ylabel(sasa_column)
        plt.title(f"MS_DSA vs {sasa_column}")
        plt.legend()
        plt.show()

    def histograms(self, run_index: int=None):
        ms_dsa = self.ms_dsa_by_run[run_index]
        ms_dsa = ms_dsa[ms_dsa["skipped"] == False]
        inter = ms_dsa[ms_dsa["is_strictly_inter_chain_contact_pdb"]]
        # intra_pdb = ms_dsa[ms_dsa["is_intra_chain_contact_pdb"]]
        # intra_AF = ms_dsa[ms_dsa["is_intra_chain_contact_AF"]]
        intra = ms_dsa[ms_dsa["is_intra_chain_contact"]]
        not_contact = ms_dsa[~ms_dsa["is_inter_chain_contact_pdb"] & ~ms_dsa["is_intra_chain_contact"]]
        plt.hist(inter["MS_DSA"].values, bins=100, label="Inter", alpha=0.5)
        plt.hist(intra["MS_DSA"].values, bins=100, label="Intra", alpha=0.5)
        plt.hist(not_contact["MS_DSA"].values, bins=100, label="Not Contact", alpha=0.5)
        plt.xlabel("MS_DSA")
        plt.ylabel("Frequency")
        plt.title("MS_DSA Distribution")
        plt.legend()
        plt.show()
      
    def intra_chain_consistency(self, run_index: int=None):
        ms_dsa = self.ms_dsa_by_run[run_index]
        ms_dsa = ms_dsa[ms_dsa["skipped"] == False]
        A = set(ms_dsa.index[ms_dsa["is_intra_chain_contact_pdb"]])
        B = set(ms_dsa.index[ms_dsa["is_intra_chain_contact_AF"]])
        venn2(subsets=(A, B), set_labels=("Intra (PDB)", "Intra (AF)"))
        plt.title("Intra chain consistency")
        plt.legend()
        plt.show()


    # @classmethod
    # def peptide_to_lysine_residues_in_uniprot(cls, peptide_seq: str, uniprot_id: str):
    #     unp_residues = [None]
    #     start_index = None
    #     sequence = UniprotToSequence.get_sequence(uniprot_id)
    #     if sequence is not None and len(sequence) > 0:
    #         idx = sequence.find(peptide_seq) + 1
    #         if start_index is not None:
    #             unp_residues = [(i + start_index, "K") for i in range(len(peptide_seq)) if peptide_seq[i] == "K"]
    #     return unp_residues, start_index

    @classmethod
    def peptide_to_lysine_residues_in_uniprot(cls, peptide_seq: str, uniprot_id: str):
        sequence = UniprotToSequence.get_sequence(uniprot_id)
        if not isinstance(peptide_seq, str) or not sequence:
            return [], None
        idx = sequence.find(peptide_seq)
        if idx == -1:
            return [], None
        start_index = idx + 1
        unp_residues = [(i + start_index, "K") for i, aa in enumerate(peptide_seq) if aa == "K"]
        return unp_residues, start_index

    #@classmethod
    # def split_multi_uniprot_id_rows(cls, ms_dsa: pd.DataFrame):
    #     #create empty dataframe with the same columns as ms_dsa
    #     #new_ms_dsa = pd.DataFrame(columns=ms_dsa.columns, dtype=ms_dsa.dtypes)
    #     new_ms_dsa = ms_dsa.iloc[0:0].copy()
    #     new_ms_dsa["Uniprot"] = pd.Series(dtype="string")
    #     for _, row in ms_dsa.iterrows():
    #         uniprot_ids = [id.strip().upper() for id in row["Protein.Ids"].split(";")]
    #         for uniprot_id in uniprot_ids:
    #             new_row = row.copy()
    #             new_row["Uniprot"] = uniprot_id
    #             new_ms_dsa = new_ms_dsa.append(new_row, ignore_index=True)
    #     new_ms_dsa = new_ms_dsa.reset_index(drop=True)
    #     return new_ms_dsa
    #@classmethod
    def split_multi_uniprot_id_rows(self, ms_dsa: pd.DataFrame):
        new_ms_dsa = ms_dsa.copy()
        new_ms_dsa["Uniprot"] = new_ms_dsa["Protein.Ids"].str.split(";")
        new_ms_dsa = new_ms_dsa.explode("Uniprot", ignore_index=True)
        new_ms_dsa["Uniprot"] = new_ms_dsa["Uniprot"].str.strip().str.upper()
        return new_ms_dsa
    
    # @classmethod
    # def align_peptides(cls, peptide: str, sequence: str):
    #     return sequence.find(peptide) + 1
    # @classmethod
    # def lysine_index(cls, peptide_seq: str, start_index: int):
    #     seq = peptide_seq.strip().upper()
    #     if not seq:
    #         raise ValueError("query sequence is empty")
    #     seq_obj = Seq(seq)
    #     #return all indices of lysine
    #     indices = [i + start_index for i, res in enumerate(str(seq_obj)) if res == "K"]
    #     return indices

    
        
if __name__ == "__main__":
    binding_interface_eval = BIEvaluation()
    binding_interface_eval.plot_boxplot()
    binding_interface_eval.plot_dsa_vs_sasa_with_binding_annotation(sasa_column="asa")
    binding_interface_eval.plot_dsa_vs_sasa_with_binding_annotation(sasa_column="RelASA")