import pandas as pd
from Bio.Seq import Seq
import numpy as np
from matplotlib import pyplot as plt
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
            self.save_ms_dsa()

    def load_ms(self):
        return pd.read_parquet(self.ms_path)
    
    def load_ms_dsa(self):
        return pd.read_csv(self.ms_dsa_path)

    def save_ms_dsa(self):
        self.ms_dsa.to_csv(self.ms_dsa_path, index=False)
    

    def build_ms_dsa(self):
        self.light_areas = []
        self.heavy_areas = []
        df_ms = self.keep_sequences_with_lysine(self.ms)
        seq_dsa_dict = {"peptide_seq": [], "Run":[], "Run.Index":[], "MS_DSA": [], "Protein.Names": [], "Genes": [], "Protein.Ids": [], "Protein.Group": [], "multiple_entries":[]}
        for (seq, run_index), seq_df in df_ms.groupby(["Stripped.Sequence", "Run.Index"]):
            assert seq_df["Protein.Names"].nunique() == 1
            assert seq_df["Genes"].nunique() == 1
            assert seq_df["Protein.Ids"].nunique() == 1
            assert seq_df["Protein.Group"].nunique() == 1
            assert seq_df["Run"].nunique() == 1
            run = seq_df["Run"].iloc[0]
            seq_dsa_dict["peptide_seq"].append(seq)
            seq_dsa_dict["Run.Index"].append(run_index)
            seq_dsa_dict["Run"].append(run)
            protein_names = seq_df["Protein.Names"].iloc[0]
            genes = seq_df["Genes"].iloc[0]
            protein_ids = seq_df["Protein.Ids"].iloc[0]
            protein_group = seq_df["Protein.Group"].iloc[0]
            seq_dsa_dict["Protein.Names"].append(protein_names)
            seq_dsa_dict["Genes"].append(genes)
            seq_dsa_dict["Protein.Ids"].append(protein_ids)
            seq_dsa_dict["Protein.Group"].append(protein_group)
            #dsa_list = []
            #for _, run_df in seq_df.groupby("Run"):
            # ms1_areas_light = run_df[run_df["Channel"] == '0']['Ms1.Area'].values
            # ms1_areas_heavy = run_df[run_df["Channel"] == '8']['Ms1.Area'].values
            ms1_areas_light = seq_df[seq_df["Channel"] == '0']['Ms1.Area'].values
            ms1_areas_heavy = seq_df[seq_df["Channel"] == '8']['Ms1.Area'].values
            if len(ms1_areas_light) > 1 or len(ms1_areas_heavy) > 1:
                print(f"Warning: more than one light or heavy area found for sequence {seq}")
                print(f"Light areas: {ms1_areas_light}")
                print(f"Heavy areas: {ms1_areas_heavy}")
                seq_dsa_dict["multiple_entries"].append(True)
            else:
                seq_dsa_dict["multiple_entries"].append(False)
            # self.light_areas.extend(ms1_areas_light)
            # self.heavy_areas.extend(ms1_areas_heavy)
            # # print(f"number of light ms1 areas = 0: {(ms1_areas_light == 0.0).sum()}")
            # print(f"number of heavy ms1 areas = 0: {(ms1_areas_heavy == 0.0).sum()}")
            # light_sum = np.sum(ms1_areas_light[ms1_areas_light > 0.0])
            # heavy_sum = np.sum(ms1_areas_heavy[ms1_areas_heavy > 0.0])
            light_sum = np.sum(ms1_areas_light)
            heavy_sum = np.sum(ms1_areas_heavy)
            dsa = light_sum / (light_sum + heavy_sum)
            # dsa_list.append(dsa)
            # dsa_mean = np.mean(dsa_list(dsa_list > 0.0))
            #seq_dsa_dict["MS_DSA"].append(dsa_mean)
            seq_dsa_dict["MS_DSA"].append(dsa)
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

def plot_light_heavy_areas(light_areas: list, heavy_areas: list):
    bins = np.linspace(0, 1.001, 1000)
    plt.hist(light_areas, bins=bins, label="light")
    plt.hist(heavy_areas, bins=bins, label="heavy")
    plt.legend()
    plt.xlim(0, 1)
    plt.show()
    #plt.savefig("light_heavy_areas.png")
    #plt.close()

if __name__ == "__main__":
    ms_dsa = MS_DSA_Builder(ms_data_relativePath=ms_data_relativePath)
    plot_light_heavy_areas(ms_dsa.light_areas, ms_dsa.heavy_areas)
    print(f"Uniprot IDs in ms_dsa: {ms_dsa.ms_dsa_path}")
    print(ms_dsa.get_uniprot_ids())
    print(f"Uniprot IDs in all ms_dsa files under: {MS_DATA_DIR}")
    print(MS_DSA_Builder.get_all_uniprot_ids())