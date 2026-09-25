"""
This script is used to generate the contacts for the structures in the database.
"""
from dsa.binding_interfaces.config import ALPHAFOLD_STRUCTURES_TO_PROCESS
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappingsStore
from dsa.binding_interfaces.structure.structure_store import StructureFileStore
from dsa.binding_interfaces.factory import ClassFactory
PDB_STRUCTURES_TO_PROCESS = UniprotStructureMappingsStore("pdb").get_mapper().get_mapped_structures()
from dsa.binding_interfaces.contacts.structure_contacts import StructureContactsStore, StructureContacts
#from dsa.binding_interfaces.contacts.filters import lysine_only_filter
from dsa.misc.compute import ComputeManager

def main(database='pdb', distance_cutoff = 6.5, contact_filter_key = "lysine_only", update_interval = 25, show_progress = True):
    # database = "pdb"
    # #database = "alphafold"
    # distance_cutoff = 6.5
    # contact_filter = lysine_only_filter
    # update_interval = 25
    # show_progress = True
    if database == "pdb":
        structure_ids = PDB_STRUCTURES_TO_PROCESS
    elif database == "alphafold":
        structure_ids = ALPHAFOLD_STRUCTURES_TO_PROCESS
    else:
        raise ValueError(f"Invalid database: {database}")

    structure_file_store = ClassFactory.create(StructureFileStore, database)
    structure_ids = structure_file_store.get_saved_structures(from_structures=structure_ids)
    
    store = StructureContactsStore(database, distance_cutoff, contact_filter_key)

    contact_object_generator = lambda structure_id: StructureContacts(structure_id, database, distance_cutoff, contact_filter_key)

    compute_manager = ComputeManager(structure_ids, 
                                     compute_function=contact_object_generator, 
                                     store=store,
                                     update_interval=update_interval,
                                     show_progress=show_progress)
    compute_manager.compute()

        

if __name__ == "__main__":
    db = "alphafold"
    #db = "pdb"
    main(database=db, distance_cutoff=6.5, contact_filter_key="lysine_only", update_interval=25, show_progress=True)
