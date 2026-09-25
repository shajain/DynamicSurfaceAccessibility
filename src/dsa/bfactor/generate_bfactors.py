"""
This script is used to compute the lysine B-factors of the PDB structures,
the same set of structures for which contacts and multimer SASA are generated.
"""
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappingsStore
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure.structure_store import StructureFileStore
from dsa.bfactor.lysine_bfactor import LysineBFactor, LysineBFactorStore
from dsa.misc.compute import ComputeManager
PDB_STRUCTURES_TO_PROCESS = UniprotStructureMappingsStore("pdb").get_mapper().get_mapped_structures()


def compute_lysine_bfactor(structure_id: str) -> LysineBFactor|None:
    # ComputeManager treats None as an unsuccessful key
    try:
        return LysineBFactor(structure_id)
    except Exception as e:
        print(f"Error computing B-factors for {structure_id}: {e}")
        return None


def main(update_interval = 25, show_progress = True):
    structure_file_store = ClassFactory.create(StructureFileStore, LysineBFactor.DATABASE)
    structure_ids = structure_file_store.get_saved_structures(from_structures=PDB_STRUCTURES_TO_PROCESS)

    store = LysineBFactorStore()

    compute_manager = ComputeManager(structure_ids,
                                     compute_function=compute_lysine_bfactor,
                                     store=store,
                                     update_interval=update_interval,
                                     show_progress=show_progress)
    compute_manager.compute()


if __name__ == "__main__":
    main(update_interval=25, show_progress=True)
