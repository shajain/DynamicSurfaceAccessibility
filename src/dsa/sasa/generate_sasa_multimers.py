"""
This script is used to compute the SASA of the biological assemblies of the PDB structures,
the same set of structures for which contacts are generated.
"""
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappingsStore
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure.structure_store import StructureFileStore
from dsa.sasa.sasa_multimer import SASAMultimer, SASAMultimerStore
from dsa.misc.compute import ComputeManager
PDB_STRUCTURES_TO_PROCESS = UniprotStructureMappingsStore("pdb").get_mapper().get_mapped_structures()


def compute_sasa_multimer(structure_id: str) -> SASAMultimer|None:
    # ComputeManager treats None as an unsuccessful key
    try:
        return SASAMultimer(structure_id)
    except Exception as e:
        print(f"Error computing SASA for {structure_id}: {e}")
        return None


def main(update_interval = 25, show_progress = True):
    structure_file_store = ClassFactory.create(StructureFileStore, SASAMultimer.DATABASE)
    structure_ids = structure_file_store.get_saved_structures(from_structures=PDB_STRUCTURES_TO_PROCESS)

    store = SASAMultimerStore()

    compute_manager = ComputeManager(structure_ids,
                                     compute_function=compute_sasa_multimer,
                                     store=store,
                                     update_interval=update_interval,
                                     show_progress=show_progress)
    compute_manager.compute()


if __name__ == "__main__":
    main(update_interval=25, show_progress=True)
