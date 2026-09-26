from Bio.Seq import Seq
from dsa.binding_interfaces.binding_interfaces import BindingInterfaces, BindingInterfacesStore
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

seq_dsa_dict = {"peptide_seq": [], "MS_DSA": [], "Protein.Names": [], "Genes": [], "Protein.Ids": [], "Protein.Group": []}
class BIEvaluation:
    def __init__(self, ms_data_relativePath: str=ms_data_relativePath, cutoff: float=6.5, lysine_only: bool=True):
        self.cutoff = cutoff
        self.lysine_only = lysine_only
        self.ms_data_relativePath = ms_data_relativePath
        self.ms_dsa = MS_DSA_Builder(self.ms_data_relativePath).ms_dsa
        self.run_indices = self.ms_dsa["Run.Index"].unique()
        self.ms_dsa_by_run = {run_index: self.ms_dsa[self.ms_dsa["Run.Index"] == run_index] for run_index in self.run_indices}
        #self.ms_dsa_avg_over_runs = self.build_ms_dsa_avg_over_runs()
        for run_index in self.run_indices:
            self.ms_dsa_by_run[run_index] = self.add_binding_interface_and_sasa_scores_info(self.ms_dsa_by_run[run_index])
        #self.sasa_df = self.build_sasa_df()
        #self.combined_df = self.combine_ms_dsa_and_sasa_df()
    def combine_ms_dsa_and_sasa_df(self):
        return pd.merge(self.sasa_df, self.ms_dsa, left_on=["Sequence"], right_on=["peptide_seq"], how="inner")

    def build_sasa_df(self):
        sasa_df_raw = pd.read_csv(DATA_DIR / "sasa" / "afold.csv")
        sasa_df = sasa_df_raw[["Sequence", "Keys", "RelASA", "asa"]].groupby(["Sequence", "Keys"]).agg({"RelASA": "mean", "asa": "mean"}).reset_index()
        return sasa_df
    
    # def build_ms_dsa_avg_over_runs(self):

    #     for (seq, run_index), seq_df in self.ms_dsa.groupby(["peptide_seq", "Run.Index"]):


    def add_binding_interface_and_sasa_scores_info(self, ms_dsa: pd.DataFrame):
        sasa_monomer_store = SASAMonomerStore()
        ms_dsa["peptide_start_index"] = None # This column will be added to the ms_dsa dataframe
        ms_dsa["is_inter_chain_contacts"] = False # This column will be added to the ms_dsa dataframe
        ms_dsa["is_intra_chain_contacts"] = False # This column will be added to the ms_dsa dataframe
        ms_dsa["asa_avg"] = np.nan # This column will be added to the ms_dsa dataframe
        ms_dsa["asa_min"] = np.nan # This column will be added to the ms_dsa dataframe
        ms_dsa["asa_max"] = np.nan # This column will be added to the ms_dsa dataframe
        ms_dsa["rel_asa_avg"] = np.nan # This column will be added to the ms_dsa dataframe
        ms_dsa["rel_asa_min"] = np.nan # This column will be added to the ms_dsa dataframe
        ms_dsa["rel_asa_max"] = np.nan # This column will be added to the ms_dsa dataframe
        ms_dsa["asa_max"] = np.nan # This column will be added to the ms_dsa dataframe
        ms_dsa["skipped"] = True # This column will be added to the ms_dsa dataframe
        for index, row in ms_dsa.iterrows():
            peptide_seq = row["peptide_seq"]
            uniprot_ids = row["Protein.Ids"]
            uniprot_ids = [id.strip().upper() for id in uniprot_ids.split(";")]
            start_indices = []
            is_inter_chain_contact_pdb = False
            is_intra_chain_contact_pdb   = False
            is_intra_chain_contact_AF = False
            
            if len(uniprot_ids) > 1:
                warnings.warn(f"Multiple uniprot_ids found for {peptide_seq}")
                continue
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
                lysine_tuples = [(li, sequence[li-1]) for li in lysine_indices]
                if np.max(lysine_indices) > len(sequence) or np.any([AA != "K" for li, AA in lysine_tuples]):
                    warnings.warn(f"Lysine not correctly mapped to sequence {sequence[start_index:start_index+len(peptide_seq)]}")
                    continue
                start_indices.append(start_index)
                binding_interfaces = BindingInterfacesStore(structure_db="pdb", cutoff=self.cutoff, lysine_only=self.lysine_only).load(uniprot_id)
                binding_interfaces_AF = BindingInterfacesStore(structure_db="alphafold", cutoff=self.cutoff, lysine_only=self.lysine_only).load(uniprot_id)
                if binding_interfaces is None or binding_interfaces.is_empty():
                    warnings.warn(f"Binding interfaces not found for uniprot_id {uniprot_id} in pdb")
                    continue
                if binding_interfaces_AF is None or binding_interfaces_AF.is_empty():
                    warnings.warn(f"Binding interfaces not found for uniprot_id {uniprot_id} in alphafold")
                    continue
                inter_chain = binding_interfaces.get_interfaces_by_type("inter-chain")
                intra_chain = binding_interfaces.get_interfaces_by_type("intra-chain")
                intra_chain_AF = binding_interfaces_AF.get_interfaces_by_type("intra-chain")
                inter_chain_lysines = [lys_pair for lys_pair in lysine_tuples if lys_pair in inter_chain]
                intra_chain_lysines = [lys_pair for lys_pair in lysine_tuples if lys_pair in intra_chain]
                intra_chain_lysines_AF = [lys_pair for lys_pair in lysine_tuples if lys_pair in intra_chain_AF]
                is_inter_chain_contact_pdb = len(inter_chain_lysines) > 0 or is_inter_chain_contact_pdb 
                is_intra_chain_contact_pdb = len(intra_chain_lysines) > 0 or is_intra_chain_contact_pdb
                is_intra_chain_contact_AF = len(intra_chain_lysines_AF) > 0 or is_intra_chain_contact_AF
                sasa_monomer = sasa_monomer_store.load(uniprot_id)
                if sasa_monomer is None:
                    warnings.warn(f"SASA not found for uniprot_id {uniprot_id}")
                    continue
                sasa_sequence_consistency = np.all(np.array(sasa_monomer.get_residue_names(lysine_indices)) == 'K')
                if not sasa_sequence_consistency:
                    warnings.warn(f"SASA sequence consistency not found for uniprot_id {uniprot_id}")
                    continue
                asa = [sasa_monomer.get_sasa_for_residue(li, type="asa") for li in lysine_indices]
                rel_asa = [sasa_monomer.get_sasa_for_residue(li, type="rel_asa") for li in lysine_indices]
                asa_nz = [sasa_monomer.get_sasa_for_residue(li, type="asa_nz") for li in lysine_indices]
                ms_dsa.at[index, "asa_avg"] = np.mean(asa)
                ms_dsa.at[index, "asa_min"] = np.min(asa)
                ms_dsa.at[index, "asa_max"] = np.max(asa)
                ms_dsa.at[index, "rel_asa_avg"] = np.mean(rel_asa)
                ms_dsa.at[index, "rel_asa_min"] = np.min(rel_asa)
                ms_dsa.at[index, "rel_asa_max"] = np.max(rel_asa)
                ms_dsa.at[index, "asa_nz_avg"] = np.mean(asa_nz)
                ms_dsa.at[index, "asa_nz_min"] = np.min(asa_nz)
                ms_dsa.at[index, "asa_nz_max"] = np.max(asa_nz)
                ms_dsa.at[index, "is_inter_chain_contact_pdb"] = is_inter_chain_contact_pdb
                ms_dsa.at[index, "is_intra_chain_contact_pdb"] = is_intra_chain_contact_pdb
                ms_dsa.at[index, "is_intra_chain_contact_AF"] = is_intra_chain_contact_AF
                is_intra_chain_contact = is_intra_chain_contact_pdb or is_intra_chain_contact_AF
                ms_dsa.at[index, "is_intra_chain_contact"] = is_intra_chain_contact
                ms_dsa.at[index, "is_strictly_inter_chain_contact_pdb"] = is_inter_chain_contact_pdb and not is_intra_chain_contact
                ms_dsa.at[index, "is_strictly_intra_chain_contact_pdb"] = is_intra_chain_contact and not is_inter_chain_contact_pdb
                ms_dsa.at[index, "peptide_start_index"] = start_indices
                ms_dsa.at[index, "skipped"] = False
        return ms_dsa
            

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
        indices = [i + start_index + 1 for i, res in enumerate(str(seq_obj)) if res == "K"]
        return indices
        
if __name__ == "__main__":
    binding_interface_eval = BIEvaluation()
    binding_interface_eval.plot_boxplot()
    binding_interface_eval.plot_dsa_vs_sasa_with_binding_annotation(sasa_column="asa")
    binding_interface_eval.plot_dsa_vs_sasa_with_binding_annotation(sasa_column="RelASA")