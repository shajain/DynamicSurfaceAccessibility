import pandas as pd
from Bio.Seq import Seq
import numpy as np
from dsa.ms_dsa.config import MS_DATA_DIR, ms_data_relativePath


class MS_DSA_Builder:
    def __init__(self, ms_data_relativePath: str=ms_data_relativePath):
        self.ms_data_relativePath = ms_data_relativePath
        self.ms_path = MS_DATA_DIR / ms_data_relativePath
        assert self.ms_path.exists(), f"MS data file not found at {self.ms_path}"
        self.ms_dsa_path = MS_DATA_DIR / "ms_dsa.csv"
        self.ms = pd.read_parquet(self.ms_path)
        if self.ms_dsa_path.exists():
            self.ms_dsa = self.load_ms_dsa()
        else:
            self.ms_dsa = self.build_ms_dsa()
            self.save_ms_dsa_df()

    def load_ms(self):
        return pd.read_parquet(self.ms_path)
    
    def load_ms_dsa(self):
        return pd.read_csv(self.ms_dsa_path)

    def save_ms_dsa(self):
        self.ms_dsa.to_csv(self.ms_dsa_path, index=False)
    

    def build_ms_dsa_df(self):
        df_ms = self.keep_sequences_with_lysine(self.df_ms)
        seq_dsa_dict = {"peptide_seq": [], "MS_DSA": [], "Protein.Names": [], "Genes": [], "Protein.Ids": [], "Protein.Group": []}
        for seq, seq_df in df_ms.groupby("Stripped.Sequence"):
            seq_dsa_dict["peptide_seq"].append(seq)
            assert seq_df["Protein.Names"].nunique() == 1
            assert seq_df["Genes"].nunique() == 1
            assert seq_df["Protein.Ids"].nunique() == 1
            assert seq_df["Protein.Group"].nunique() == 1
            protein_names = seq_df["Protein.Names"].iloc[0]
            genes = seq_df["Genes"].iloc[0]
            protein_ids = seq_df["Protein.Ids"].iloc[0]
            protein_group = seq_df["Protein.Group"].iloc[0]
            seq_dsa_dict["Protein.Names"].append(protein_names)
            seq_dsa_dict["Genes"].append(genes)
            seq_dsa_dict["Protein.Ids"].append(protein_ids)
            seq_dsa_dict["Protein.Group"].append(protein_group)
            dsa_list = []
            for _, run_df in seq_df.groupby("Run"):
                ms1_area_light = run_df[run_df["Channel"] == '0']['Ms1.Area'].sum()
                ms1_area_heavy = run_df[run_df["Channel"] == '8']['Ms1.Area'].sum()
                dsa = ms1_area_light / (ms1_area_light + ms1_area_heavy)
                dsa_list.append(dsa)
            dsa_mean = np.mean(dsa_list)
            seq_dsa_dict["MS_DSA"].append(dsa_mean)
        #create dataframe frm dict the key is one column and the value is another column
        df_ms_dsa = pd.DataFrame(seq_dsa_dict)
        return df_ms_dsa

    def get_ms_dsa_df(self):
        return self.ms_dsa

    def get_uniprot_ids(self):
        ids = self.ms_dsa["Protein.Ids"].tolist()
        ids = [id.split(";") for id in ids]
        ids = [id for sublist in ids for id in sublist]
        return list(set(ids))

    @staticmethod
    def keep_sequences_with_lysine(ms: pd.DataFrame):
        return ms[ms["Stripped.Sequence"].apply(MS_DSA_Builder.query_contains_lysine)]

    @staticmethod
    def query_contains_lysine(query: str):
        q = query.strip().upper()
        if not q:
            raise ValueError("query sequence is empty")
        seq_obj = Seq(q)
        return "K" in str(seq_obj)

    @staticmethod
    def columns_to_keep():
        cols = [ 'Run.Index', 'Run', 'Precursor.Id', 'Stripped.Sequence', 'Modified.Sequence', 'Channel', 'Precursor.Charge', 'Protein.Ids', 'Protein.Group', 'Protein.Names', 'Genes', 'Ms1.Area', 'Ms1.Normalised']
        return cols
    
    @staticmethod
    def get_all_uniprot_ids(ms_data_dir = MS_DATA_DIR):
        uniprot_ids = set()
        for file in MS_DATA_DIR.glob("**/*parquet"):
            #Get relative path from ms_data_dir
            relative_path = file.relative_to(ms_data_dir)
            ids = MS_DSA_Builder(relative_path).get_uniprot_ids()
            uniprot_ids.update(ids)
        return list(uniprot_ids)

if __name__ == "__main__":
    ms_dsa = MS_DSA_Builder(ms_data_relativePath=ms_data_relativePath)
    print(f"Uniprot IDs in ms_dsa: {ms_dsa.ms_dsa_path}")
    print(ms_dsa.get_uniprot_ids())
    print(f"Uniprot IDs in all ms_dsa files under: {MS_DATA_DIR}")
    print(MS_DSA_Builder.get_all_uniprot_ids())