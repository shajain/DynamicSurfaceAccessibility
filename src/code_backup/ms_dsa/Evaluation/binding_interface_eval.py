from Bio.Seq import Seq
from dsa.binding_interfaces.binding_interfaces import BindingInterfaces
from dsa.ms_dsa.ms_dsa import MS_DSA_Builder
from dsa.ms_dsa.config import ms_data_relativePath
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
import numpy as np
from matplotlib import pyplot as plt
import warnings
from dsa.config import FIGURES_DIR, DATA_DIR
import pandas as pd
seq_dsa_dict = {"peptide_seq": [], "MS_DSA": [], "Protein.Names": [], "Genes": [], "Protein.Ids": [], "Protein.Group": []}
class BIEvaluation:
    def __init__(self, ms_data_relativePath: str=ms_data_relativePath):
        self.ms_data_relativePath = ms_data_relativePath
        self.ms_dsa = MS_DSA_Builder(self.ms_data_relativePath).ms_dsa
        self.add_binding_interface_and_sasa_scores_info()
        self.sasa_df = self.build_sasa_df()
        self.combined_df = self.combine_ms_dsa_and_sasa_df()
    def combine_ms_dsa_and_sasa_df(self):
        return pd.merge(self.sasa_df, self.ms_dsa, left_on=["Sequence"], right_on=["peptide_seq"], how="inner")

    def build_sasa_df(self):
        sasa_df_raw = pd.read_csv(DATA_DIR / "sasa" / "afold.csv")
        sasa_df = sasa_df_raw[["Sequence", "Keys", "RelASA", "asa"]].groupby(["Sequence", "Keys"]).agg({"RelASA": "mean", "asa": "mean"}).reset_index()
        return sasa_df



    def add_binding_interface_and_sasa_scores_info(self):
        self.ms_dsa["peptide_start_index"] = None # This column will be added to the ms_dsa dataframe
        self.ms_dsa["is_inter_chain_contacts"] = False # This column will be added to the ms_dsa dataframe
        self.ms_dsa["is_intra_chain_contacts"] = False # This column will be added to the ms_dsa dataframe
        for index, row in self.ms_dsa.iterrows():
            peptide_seq = row["peptide_seq"]
            uniprot_ids = row["Protein.Ids"]
            uniprot_ids = [id.strip().upper() for id in uniprot_ids.split(";")]
            start_indices = []
            is_inter_chain_contacts = False
            is_intra_chain_contacts = False
            for uniprot_id in uniprot_ids:
                sequence = UniprotToSequence.get_sequence(uniprot_id)
                if sequence is None or len(sequence) == 0:
                    warnings.warn(f"Sequence not found for uniprot_id {uniprot_id}")
                    continue
                start_index = self.align_peptides(peptide_seq, sequence)
                if start_index is None:
                    warnings.warn(f"Peptide not found in sequence {sequence}")
                    continue
                lysine_indices = self.lysine_index(peptide_seq, start_index) 
                if np.max(lysine_indices) > len(sequence)-1 or np.any([sequence[li] != "K" for li in lysine_indices]):
                    warnings.warn(f"Lysine not correctly mapped to sequence {sequence[start_index:start_index+len(peptide_seq)]}")
                    continue
                start_indices.append(start_index)
                binding_interfaces_inter = BindingInterfaces.binding_interface(uniprot_id, "pdb", 8.0, "inter-chain")
                binding_interfaces_intra = BindingInterfaces.binding_interface(uniprot_id, "pdb", 8.0, "intra-chain")
                if not is_inter_chain_contacts:
                    is_inter_chain_contacts = np.any([li in [bi[0] for bi in binding_interfaces_inter] for li in lysine_indices]) 
                if not is_intra_chain_contacts:
                    is_intra_chain_contacts = np.any([li in [bi[0] for bi in binding_interfaces_intra] for li in lysine_indices])
            self.ms_dsa.at[index, "is_inter_chain_contacts"] = is_inter_chain_contacts
            self.ms_dsa.at[index, "is_intra_chain_contacts"] = is_intra_chain_contacts
            self.ms_dsa.at[index, "peptide_start_index"] = start_indices

    def plot_boxplot(self):
        print(self.ms_dsa["is_inter_chain_contacts"].dtype)
        print(self.ms_dsa["is_inter_chain_contacts"].value_counts(dropna=False))

        inter = self.ms_dsa["MS_DSA"][self.ms_dsa["is_inter_chain_contacts"]]
        not_inter = self.ms_dsa["MS_DSA"][~self.ms_dsa["is_inter_chain_contacts"]]
        intra = self.ms_dsa["MS_DSA"][self.ms_dsa["is_intra_chain_contacts"]]
        not_intra = self.ms_dsa["MS_DSA"][~self.ms_dsa["is_intra_chain_contacts"]]

        print(inter.count(), not_inter.count(), intra.count(), not_intra.count())
        data = [inter.dropna().values, not_inter.dropna().values, intra.dropna().values, not_intra.dropna().values]
        labels = ["Inter", "Not Inter", "Intra", "Not Intra"]
        plt.boxplot(data, labels=labels)
        plt.xlabel("Binding interface type")
        plt.ylabel("MS_DSA")
        plt.title("MS_DSA vs Binding interface type")
        plt.show()

    def plot_dsa_vs_sasa_with_binding_annotation(self, sasa_column: str="RelASA"):
        inter = self.combined_df[self.combined_df["is_inter_chain_contacts"]]
        print(f"Inter: {inter.shape}")
        print(f"% of Inter less than 0.9: {len(inter[inter['MS_DSA'] < 0.9]) / len(inter)}")
        #not_inter = self.combined_df[~self.combined_df["is_inter_chain_contacts"]]
        intra = self.combined_df[self.combined_df["is_intra_chain_contacts"]]
        print(f"Intra: {intra.shape}")
        print(f"% of Intra less than 0.9: {len(intra[intra['MS_DSA'] < 0.9]) / len(intra)}")
        not_contact = self.combined_df[~self.combined_df["is_inter_chain_contacts"] & ~self.combined_df["is_intra_chain_contacts"]]
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

      
    

    
    @staticmethod
    def align_peptides(peptide: str, sequence: str):
        return sequence.find(peptide)
    @staticmethod
    def lysine_index(peptide_seq: str, start_index: int):
        seq = peptide_seq.strip().upper()
        if not seq:
            raise ValueError("query sequence is empty")
        seq_obj = Seq(seq)
        #return all indices of lysine
        indices = [i + start_index for i, res in enumerate(str(seq_obj)) if res == "K"]
        return indices
        
if __name__ == "__main__":
    binding_interface_eval = BIEvaluation()
    binding_interface_eval.plot_boxplot()
    binding_interface_eval.plot_dsa_vs_sasa_with_binding_annotation(sasa_column="asa")
    binding_interface_eval.plot_dsa_vs_sasa_with_binding_annotation(sasa_column="RelASA")