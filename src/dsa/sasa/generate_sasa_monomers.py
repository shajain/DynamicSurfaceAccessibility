"""
This script is used to compute the SASA of the AlphaFold monomers of the uniprot ids in the MS data.
"""
from dsa.uniprot.config import ALL_UNIPROT_IDS
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure.structure_store import StructureFileStore
from dsa.sasa.sasa_monomer import SASAMonomer, SASAMonomerStore
from dsa.misc.compute import ComputeManager


def compute_sasa_monomer(uniprot_id: str) -> SASAMonomer|None:
    # ComputeManager treats None as an unsuccessful key
    try:
        return SASAMonomer(uniprot_id)
    except Exception as e:
        print(f"Error computing SASA for {uniprot_id}: {e}")
        return None


def main(update_interval = 25, show_progress = True):
    structure_file_store = ClassFactory.create(StructureFileStore, SASAMonomer.DATABASE)
    uniprot_ids = structure_file_store.get_saved_structures(from_structures=ALL_UNIPROT_IDS)

    store = SASAMonomerStore()

    compute_manager = ComputeManager(uniprot_ids,
                                     compute_function=compute_sasa_monomer,
                                     store=store,
                                     update_interval=update_interval,
                                     show_progress=show_progress)
    compute_manager.compute()


if __name__ == "__main__":
    main(update_interval=25, show_progress=True)
