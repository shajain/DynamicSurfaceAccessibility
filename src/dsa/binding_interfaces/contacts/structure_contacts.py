from copyreg import pickle
from dsa.binding_interfaces.contacts.contact import Contact
from dsa.binding_interfaces.contacts.contact_extractor import ContactExtractor
from dsa.binding_interfaces.contacts.filters import ContactFilter
from dsa.binding_interfaces.structure.structure import Structure
#from dsa.binding_interfaces.structure.structure_store import StructureStore
from dsa.binding_interfaces.config import STRUCTURE_CONTACTS_DIR
from dsa.binding_interfaces.config import ALPHAFOLD_STRUCTURES_TO_PROCESS
from dsa.misc.file_strings import FileStrings
import pandas as pd
from dsa.misc.store import CompressedPickleStore
from dsa.misc.compute import ComputeManager
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappingsStore
from dsa.binding_interfaces.factory import ClassFactory
import numpy as np
# from collections import defaultdict
# from dsa.misc.utilities import defaultdict_to_dict

PDB_STRUCTURES_TO_PROCESS = UniprotStructureMappingsStore("pdb").get_mapper().get_mapped_structures()


class StructureContacts:
    BOND_TYPES = ["salt_bridge", "hydrogen_bond", "covalent_isopeptide", "cation_pi"]
    def __init__(self, structure_id: str, 
                       database:str, 
                       distance_cutoff: float = 6.5,
                       contact_filter_key: str = "lysine_only"):
        # self.structure_id = structure_id
        # self.database = database
        # self.contact_filter = contact_filter
        self.structure = ClassFactory.create(Structure, key=database, structure_id=structure_id)
        self.contact_extractor = ContactExtractor(distance_cutoff, contact_filter_key)
        self.contacts = self.extract() 

    def extract(self):
        return self.contact_extractor.extract(self.structure)

    def __str__(self):
        return f"{self.structure}-{self.contact_extractor}"

    # def __getstate__(self):
    #     state = self.__dict__.copy()
    #     # Convert contacts to a pandas DataFrame
    #     state["contacts"] = pd.DataFrame(self.contacts)
    #     return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        # self.contacts = state["contacts"].to_dict(orient="records")
        # self.contacts = [Contact(**c) for c in self.contacts]
        [c.__setstate_structure_context__(self.structure) for c in self.contacts]

    # def is_inter_chain(self):
    #     return [c.is_inter_chain() for c in self.contacts]

    # def is_residue_bond_type(self, bond_type: str):
    #     return [(c.is_residue1_bond_type(bond_type), c.is_residue2_bond_type(bond_type)) for c in self.contacts]

    def filter_contacts(self, contact_type: str=None, bond_type: str=None, method: str=None, distance_threshold: float=None) -> set[tuple[int, str]]:
        filtered_contacts = []
        for c in self.contacts:
            if not contact_type or c.satisfies_contact_type(contact_type):
                if not bond_type or c.satisfies_bond_type(bond_type, method):
                    if not distance_threshold or c.is_distance_within(distance_threshold):
                        filtered_contacts.append(c)
        return filtered_contacts

    def summary(self) -> str:
        return f"Number of contacts: {len(self.contacts)}"


class StructureContactsStore(CompressedPickleStore):
    def __init__(self, database: str, distance_cutoff: float = 6.5, contact_filter_key: str = "lysine_only"):
        self.distance_cutoff = distance_cutoff
        self.contact_filter = ClassFactory.create(interface=ContactFilter, key=contact_filter_key)
        self.contact_extractor = ContactExtractor(distance_cutoff, contact_filter_key)
        self.database = database
        super().__init__(self.get_directory())

    def get_directory(self):
        contact_extractor_string = f"{self.contact_extractor}"
        return STRUCTURE_CONTACTS_DIR / f"{self.database}_{contact_extractor_string}"


if __name__ == "__main__":
    database = "pdb"
    #database = "alphafold"
    distance_cutoff = 6.5
    contact_filter_key = "lysine_only"
    update_interval = 25
    show_progress = True

    if database == "pdb":
        structure_ids = PDB_STRUCTURES_TO_PROCESS
    elif database == "alphafold":
        structure_ids = ALPHAFOLD_STRUCTURES_TO_PROCESS
    else:
        raise ValueError(f"Invalid database: {database}")

    store = StructureContactsStore(database, distance_cutoff, contact_filter_key)

    contact_object_generator = lambda structure_id: StructureContacts(structure_id, database, distance_cutoff, contact_filter_key)

    compute_manager = ComputeManager(structure_ids, 
                                     compute_function=contact_object_generator, 
                                     store=store,
                                     update_interval=35,
                                     show_progress=show_progress)
    compute_manager.compute()

        