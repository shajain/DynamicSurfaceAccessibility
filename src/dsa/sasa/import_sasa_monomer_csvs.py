"""
This script is used to import the monomer SASA csv files written by compute_sasa.py into the
SASAMonomerStore, without re-computing SASA.
"""
import pandas as pd
from dsa.sasa.config import SASA_MONOMER_CSV_DIR
from dsa.sasa.sasa_monomer import SASAMonomer, SASAMonomerStore
from dsa.misc.compute import ComputeManager


def load_sasa_monomer_from_csv(uniprot_id: str) -> SASAMonomer|None:
    # ComputeManager treats None as an unsuccessful key
    try:
        sasa_df = pd.read_csv(SASA_MONOMER_CSV_DIR / f"{uniprot_id}.csv")
        return SASAMonomer.from_sasa_df(uniprot_id, sasa_df)
    except Exception as e:
        print(f"Error importing SASA csv for {uniprot_id}: {e}")
        return None


def main(update_interval = 100, show_progress = True):
    uniprot_ids = [file.stem for file in SASA_MONOMER_CSV_DIR.glob("*.csv")]

    store = SASAMonomerStore()

    compute_manager = ComputeManager(uniprot_ids,
                                     compute_function=load_sasa_monomer_from_csv,
                                     store=store,
                                     update_interval=update_interval,
                                     show_progress=show_progress)
    compute_manager.compute()


if __name__ == "__main__":
    main(update_interval=100, show_progress=True)
