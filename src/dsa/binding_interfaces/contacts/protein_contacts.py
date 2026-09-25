from pickle import FALSE
from dsa.binding_interfaces.config import ALL_UNIPROT_IDS
from dsa.binding_interfaces.contacts import structure_contacts
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappingsStore
#from dsa.binding_interfaces.structure.structure import Structure
from dsa.binding_interfaces.contacts.structure_contacts import StructureContactsStore, StructureContacts
from dsa.binding_interfaces.contacts.uniprot_aligned_structure_contacts import UniprotAlignedStructureContacts
from dsa.binding_interfaces.structure_uniprot_alignment.aligner import StructureSequenceAligner   
from typing import Literal
from dsa.binding_interfaces.config import PROTEIN_CONTACTS_DIR
from Bio.SeqUtils import seq1
from pathlib import Path
import json
import warnings
from collections import defaultdict
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
from dsa.misc.store import CompressedPickleStore
from dsa.binding_interfaces.contacts.filters import ContactFilter
from dsa.binding_interfaces.contacts.contact_extractor import ContactExtractor
from itertools import product
from dsa.misc.utilities import defaultdict_to_dict
from dsa.misc.utilities import nested_dict_to_table
import numpy as np
from dsa.binding_interfaces.contacts.bond_characterization.library import BondLibrary



# The main purpose of this code is to get all binding interface for a uniprotid.
# To this end it should use use UniprotStructureMappings to get structure IDs
# corresponding to the uniprotid and then use Structure class to extract the saved binding interfaces
# for each structure. Next it should aggregate the binding interfaces for all the structures.
# The interfaces should be stored in a dictionary with  keys coresponding to the
# interface types. These keys are also present in the binding interface from Structure class.
# The value should be a list of dictionaries. Each dictionary in the list should have the following keys:
# res_num: int, res_name: str, chain: str, pdb_res_num: str, pdb_res_name: str, structure_id: str. Note
# this is different from the current impentation in this file.
# One class should be responsible for extracting the binding interfaces. Another class
# should be saving the results under BI_PROTEIN_DIR. One file for each uniprotid should be created.
# The file name should be the {uniprotid}.json. The class should also have the logic to load the
# saved binding interfaces from the file for a given uniprotid.
#Modify binding interfaces to only process one structure at a time. 
# Have a if __name__ = __main__ block to process uniprot_ids from ALL_UNIPROT_IDS minus any 
# uniporotid whoose interfaces are already saved. call this list uniport_ids_to_process  
# Pick a random UniProt_id in each iteration before attempting to extract interfaces check 
# again if the interface is alreayd saved. This is because I will run this parallely with 
# multiple threads. After every. 50 saved uniprot ids update the uniport_ids_to_process  
# by removig those that were saved by another thread. 


# Interface types as produced by the Structure / BindingInterfaceExtractor binding interfaces.
INTERFACE_TYPES = ["all", "residue-residue", "residue-other", "intra-chain", "inter-chain"]


class ProteinContactsStore(CompressedPickleStore):
    def __init__(self, structure_database: str = "pdb", distance_cutoff: float = 6.5, contact_filter_key: str = "lysine_only"):
        self.structure_database = structure_database
        self.contact_filter = ClassFactory.create(interface=ContactFilter, key=contact_filter_key)
        self.contact_extractor = ContactExtractor(distance_cutoff, contact_filter_key)
        self.contact_filter_key = contact_filter_key
        super().__init__(directory=self.get_directory(), object_type=ProteinContacts)
    

    def get_directory(self) -> Path:
        contact_extractor_str = f"{self.contact_extractor}"
        return PROTEIN_CONTACTS_DIR / "_".join([self.structure_database, contact_extractor_str])

    # def file_path(self, uniprot_id: str) -> Path:
    #     return self.directory / f"{uniprot_id}.json"

    # def is_saved(self, uniprot_id: str) -> bool:
    #     return self.file_path(uniprot_id).exists()

    # def save(self, uniprot_id: str, binding_interfaces: BindingInterfaces) -> Path:
    #     #interfaces = binding_interfaces.get_interfaces()
    #     file_path = self.file_path(uniprot_id)
    #     file_path.parent.mkdir(parents=True, exist_ok=True)
    #     with open(file_path, "w") as f:
    #         json.dump(binding_interfaces.__dict__, f, indent=4)
    #     return file_path

    # def load(self, uniprot_id: str) -> BindingInterfaces | None:
    #     file_path = self.file_path(uniprot_id)
    #     if not file_path.exists():
    #         return None
    #     with open(file_path, "r") as f:
    #         try:
    #             object_dict = json.load(f)
    #         except json.JSONDecodeError:
    #             warnings.warn(f"JSONDecodeError for {file_path}")
    #             return None
    #         binding_interfaces = BindingInterfaces.factory_method(**object_dict)
    #         return binding_interfaces

    # def saved_uniprot_ids(self) -> list[str]:
    #     """Look at all the files in the directory and return the uniprot ids of the files that are saved."""
    #     directory = self.directory
    #     uniprot_ids = [file.stem for file in directory.glob("*.json")]
    #     return uniprot_ids

    # def unsaved_uniprot_ids(self, uniprot_ids: list[str]) -> list[str]:
    #     """Select the UniProt ids from the list whose interfaces are not saved yet."""
    #     return [uid for uid in uniprot_ids if not self.is_saved(uid)]

    # def corrupted_uniprot_ids(self) -> list[str]:
    #     saved_uniprot_ids = self.saved_uniprot_ids()
    #     corrupted_uniprot_ids = [id for id in saved_uniprot_ids if self.load(id) is None]
    #     return corrupted_uniprot_ids

    # def empty_uniprot_ids(self) -> list[str]:
    #     """Look at all the files in the directory and return the uniprot ids with no interfaces."""
    #     uniprot_ids = self.saved_uniprot_ids()
    #     empty_uniprot_ids = []
    #     for unp in uniprot_ids:
    #         binding_interfaces = self.load(unp)
    #         if binding_interfaces is not None and binding_interfaces.is_empty():
    #                 empty_uniprot_ids.append(unp)
    #     #empty_uniprot_ids = [unp for unp in uniprot_ids if self.load(unp).is_empty()]
    #     return empty_uniprot_ids
    
    # def remove_corrupted_uniprot_ids(self) -> list[str]:
    #     corrupted_uniprot_ids = self.corrupted_uniprot_ids()
    #     for unp in corrupted_uniprot_ids:
    #         self.delete(unp)
    #     return corrupted_uniprot_ids

    # def delete(self, uniprot_id: str) -> None:
    #     file_path = self.file_path(uniprot_id)
    #     if file_path.exists():
    #         file_path.unlink()
    #     else:
    #         warnings.warn(f"File {file_path} does not exist")

    # def extract_and_save(self, uniprot_id: str) -> dict[str, list[dict]] | None:
    #     """Extract and save the binding interfaces for a UniProt id unless they are already saved.

    #     Returns the extracted interfaces, or None if they were already saved (and not overwritten).
    #     A `BindingInterfaces` extractor can be supplied to be reused across many calls.
    #     """
    #     binding_interfaces = BindingInterfaces.factory_method(uniprot_id, self.structure_db, self.cutoff, self.lysine_only)
    #     if binding_interfaces is not None:
    #         self.save(uniprot_id, binding_interfaces)
    #     else:
    #         warnings.warn(f"Could not extract binding interfaces for {uniprot_id}")
    #     return binding_interfaces


class ProteinContacts:
    """Extracts and aggregates binding interfaces for UniProt ids across all their structures.

    For each UniProt id the structure ids are looked up through `UniprotStructureMappings`,
    the saved per-structure contacts are loaded through `StructureContactsStore`, and the
    contacting residues are aligned to the UniProt sequence using via 'Structure' class. 
    """
    CONTACT_TYPES = ["intra_chain", "inter_chain"]
    BOND_TYPES = ["salt_bridge", "hydrogen_bond", "covalent_isopeptide", "cation_pi"]
    

    def __init__(self, uniprot_id: str,
                 structure_database: str = "pdb",
                 distance_cutoff: float = 6.5,
                 contact_filter_key: str = "lysine_only"):
        self.structure_database = structure_database
        self.distance_cutoff = distance_cutoff
        self.contact_filter_key = contact_filter_key
        self.uniprot_id = uniprot_id
        self.structure_ids = self.get_structure_ids()
        self.structure_contacts = self.load_structure_contacts()
        self.distance_thresholds = [float(dt) for dt in np.arange(3.0, self.distance_cutoff+0.01, 0.5)]
        self.grouped_residues = self.group_residues()

    def get_structure_ids(self) -> list[str]:
        structure_mapper = UniprotStructureMappingsStore.get_mapper(structure_database=self.structure_database)
        structure_ids = structure_mapper.structures_for_uniprot_id(self.uniprot_id)
        return structure_ids


    def load_structure_contacts(self) -> dict[str, dict[str, list]] | None:
        """Load the saved per-structure binding interfaces (chain -> ctype -> [(pdb_res_num, pdb_res_name)])."""
        contacts_from_structures = {}
        for structure_id in self.structure_ids:
            try:
                structure_contacts = UniprotAlignedStructureContacts(uniprot_id=self.uniprot_id, 
                                                                    structure_id=structure_id, 
                                                                    structure_database=self.structure_database, 
                                                                    distance_cutoff=self.distance_cutoff, 
                                                                    contact_filter_key=self.contact_filter_key)
            except Exception as e:
                warnings.warn(f"Could not load saved contacts for {self.structure_database} structure {structure_id} structure: {e}")
                continue
            contacts_from_structures[structure_id] = structure_contacts
        return contacts_from_structures

    def filter_contacts(self, contact_type: str=None, bond_type: str=None, method: str=None, distance_threshold: float=None) -> set[tuple[int, str]]:
        filtered_contacts = defaultdict(list)
        for structure_id, contacts in self.structure_contacts.items():
            filtered_contacts[structure_id] = contacts.filter_contacts(contact_type=contact_type, bond_type=bond_type, method=method, distance_threshold=distance_threshold)
        return filtered_contacts

    # def restrict_to_uniprot_id(self, contacts_from_structures: dict[str, StructureContacts]):
    #     contacts_from_structures = {structure_id: structure_contacts.restrict_to_uniprot_id(self.uniprot_id) for structure_id, structure_contacts in contacts_from_structures.items()}
    #     return contacts_from_structures


    # def initialize_empty_residue_groups(self):
    #     residue_groups = {}
    #     for ctype in self.CONTACT_TYPES:
    #         residue_groups[ctype] = {}
    #         for btype in self.BOND_TYPES:
    #             residue_groups[ctype][btype] = set()
    #         for dist_threshold in self.distance_thresholds:
    #             residue_groups[ctype][dist_threshold] = set()
    #     return residue_groups

    def group_residues(self):
        #grouped_residues = defaultdict(lambda: defaultdict(set))
        grouped_residues = {}
        for ctype in self.CONTACT_TYPES:
            grouped_residues[ctype] = {}
            for bond_type in BondLibrary.bond_types_implemented():
                grouped_residues[ctype][bond_type] = {}
                for method in BondLibrary.methods_implemented(bond_type):
                    grouped_residues[ctype][bond_type][method] = set()
                    for structure_contacts in self.structure_contacts.values():
                        residues = structure_contacts.get_residues_satisfying(contact_type=ctype, bond_type=bond_type, method=method)
                        grouped_residues[ctype][bond_type][method].update(residues)
        for ctype, dist_threshold in product(self.CONTACT_TYPES, self.distance_thresholds):
            grouped_residues[ctype][dist_threshold] = set()
            for structure_contacts in self.structure_contacts.values():
                residues = structure_contacts.get_residues_satisfying(contact_type=ctype, distance_threshold=dist_threshold)
                grouped_residues[ctype][dist_threshold].update(residues)
        return grouped_residues

    def does_residue_satisfy(self, unp_residues: list[tuple[int, str]], 
                                    contact_type: str=None, 
                                    bond_type: str=None, 
                                    method: str=None) -> bool:
        satisfied = [False] * len(unp_residues)
        contacts_types = [contact_type] if contact_type is not None else self.CONTACT_TYPES
        bond_types = [bond_type] if bond_type is not None else BondLibrary.bond_types_implemented()
        for ct in contacts_types:
            for bt in bond_types:
                methods = [method] if method is not None else BondLibrary.methods_implemented(bt)
                for m in methods:
                    contacts = self.grouped_residues[ct][bt][m]
                    for i, unp_residue in enumerate(unp_residues):
                        if unp_residue in contacts:
                            satisfied[i] = True
        return satisfied

    def summary_table(self):
        table = nested_dict_to_table(self.grouped_residues)
        return table
        


    
        

    # def group_residues_by_bond_types(self, contacts_from_structures: dict[str, StructureContacts]):
    #     for i, structure_contacts in enumerate(contacts_from_structures.values()):
    #         res_bt = structure_contacts.group_uniprot_residues_by_distances()
    #         if i == 0:
    #             residues_by_bond_types = res_bt
    #         else:
    #             for inter_or_intra in res_bt.keys():
    #                 for bond_type in res_bt[inter_or_intra].keys():
    #                     residues_by_bond_types[inter_or_intra][bond_type].union(res_bt[inter_or_intra][bond_type])
    #     return residues_by_bond_types


    # def group_residues_by_distances(self, contacts_from_structures: dict[str, StructureContacts]):
    #     for i, structure_contacts in enumerate(contacts_from_structures.values()):
    #         res_dist = structure_contacts.group_uniprot_residues_by_bond_types()
    #         if i == 0:
    #             residues_by_distances = res_dist
    #         else:
    #             for inter_or_intra in res_dist.keys():
    #                 for distance in res_dist[inter_or_intra].keys():
    #                     residues_by_distances[inter_or_intra][distance].union(res_dist[inter_or_intra][distance])
    #     return residues_by_distances

    # def filter_by_uniprot_id_and_map(self):
    #     uniprot_contacts = {}
    #     for structure_id, structure_contacts in self.contacts_from_structures.items():
    #         uniprot_contacts[structure_id] = structure_contacts.filter_by_uniprot_id_and_map(self.uniprot_id)
    #     return uniprot_contacts

    # def align_contacts_to_uniprot(self, contacts_from_structures: dict[str, StructureContacts]) -> dict[str, list[tuple[tuple[int, str], tuple[int, str]]]]:
    #     aligned_contacts = {}
    #     for structure_id, structure_contacts in contacts_from_structures.items():
    #         aligner = ClassFactory.create(interface=StructureSequenceAligner, 
    #                                         key=self.structure_database,
    #                                         structure = structure_contacts.structure)
            
    #         aligned_contacts[structure_id] = self._align_contact_to_uniprot(structure_contacts.contacts, aligner)
    #     return aligned_contacts

    # def _align_contact_to_uniprot(self, contacts, aligner: StructureSequenceAligner)->list[tuple[tuple[int, str], tuple[int, str]]]:
    #     unp_contacts = []
    #     for contact in contacts:
    #         chain1, residue1 = contact.author_chain1, contact.residue1
    #         chain2, residue2 = contact.author_chain2, contact.residue2
    #         unp_residue1 = aligner.map(chain1, self.uniprot_id, residue1)
    #         unp_residue2 = aligner.map(chain2, self.uniprot_id, residue2)
    #         unp_contacts.append((unp_residue1, unp_residue2))
    #     return unp_contacts

    # def group_uniprot_residues_by_bond_types(self):
    #     positions_grouped = {}
    #     positions_grouped['inter-chain']={bp:set() for bp in self.BOND_TYPES}
    #     positions_grouped['intra-chain']={bp:set() for bp in self.BOND_TYPES}
    #     for structure_id in self.structure_ids:
    #         structure_contacts = self.contacts_from_structures[structure_id]
    #         unp_contacts = self.uniprot_aligned_contacts[structure_id]
    #         for i in range(len(unp_contacts)):
    #             unp_contact = unp_contacts[i]
    #             structure_contact = structure_contacts[i]
    #             inter_or_intra = "intra-chain" if structure_contact.is_intra_chain() else "inter-chain"
    #             if unp_contact[0] is not None:
    #                 for bp in structure_contact.get_residue1_bond_types():
    #                     positions_grouped[inter_or_intra][bp].add(unp_contact[0])
    #             if unp_contact[1] is not None:
    #                 for bp in structure_contact.get_residue2_bond_types():
    #                     positions_grouped[inter_or_intra][bp].add(unp_contact[1])
    #     return positions_grouped

    # def group_uniprot_residues_by_distances(self):
    #     distance_thresholds = np.arange(3.0, self.distance_cutoff, 0.5)
    #     positions_grouped = {}
    #     positions_grouped['intra-chain'] = {dist:set() for dist in distance_thresholds}
    #     positions_grouped['inter-chain'] = {dist:set() for dist in distance_thresholds}
    #     for structure_id in self.structure_ids:
    #         structure_contacts = self.contacts_from_structures[structure_id]
    #         unp_contacts = self.uniprot_aligned_contacts[structure_id]
    #         for i in range(len(unp_contacts)):
    #             unp_contact = unp_contacts[i]
    #             structure_contact = structure_contacts[i]
    #             inter_or_intra = "intra-chain" if structure_contact.is_intra_chain() else "inter-chain"
    #             distance = structure_contact.distance
    #             for dist_threshold in distance_thresholds:
    #                 if distance <= dist_threshold:
    #                     positions_grouped[inter_or_intra][dist_threshold].add(unp_contact)
    #     return positions_grouped


    # def is_empty(self) -> bool:
    #     return all(len(interfaces) == 0 for interfaces in self.contacts.values())

    # def extract_contacts(self) -> dict[str, list[dict]]:
    #     structure_mapper = UniprotStructureMappingsStore.get_mapper(structure_database=self.structure_database)
    #     structure_ids = structure_mapper.structures_for_uniprot_id(self.uniprot_id)
    #     # Process one structure at a time; the mapper is built and discarded per structure
    #     # so that memory stays bounded when running many extractions in parallel.
    #     for structure_id in structure_ids:
    #         self.extract_from_saved_structure_contacts(structure_id)

    # def get_contacts(self) -> dict[str, list[dict]]:
    #     return self.contacts

    # def set_contacts(self, contacts: dict[str, list[dict]]):
    #     self.contacts = contacts

    # def extract_from_saved_structure_contacts(self, structure_id: str):
    #     """Map and aggregate the contacts from a single structure for a given UniProt id."""
    #     structure_contacts = self.structure_contacts_store[self.structure_database].load(structure_id)
    #     if structure_contacts is None:
    #         warnings.warn(f"Saved contacts not found for {self.structure_database} structure {structure_id} structure or error in loading")
    #         return
    #     all_contacts = [len(contacts[c]["all"]) for c in contacts.keys()]
    #     if np.max(all_contacts) > 0:
    #         print(f"Structure has contacts for uniprot id {self.uniprot_id}")
    #     aligner = self.build_aligner(structure_id)
    #     if aligner is None:
    #         warnings.warn(f"Could not build PDB_Uniprot_Mapper for structure {structure_id}")
    #         return
    #     if contacts is None:
    #         warnings.warn(f"Could not load contacts for structure {structure_id}")
    #         return
    #     chains = aligner.chains_containing_uniprot(self.uniprot_id)
    #     for chain in chains:   
    #         self.build_from_chain(contacts, aligner, chain, structure_id)

    # def build_from_chain(self, contacts: dict[str, dict[str, list]], aligner: Structure_Sequence_Aligner, chain: str, structure_id: str):
    #     if not aligner.is_alignment_successful(chain, self.uniprot_id):
    #         warnings.warn(f"Alignment not successful for chain {chain} in structure {structure_id} with uniprot id {self.uniprot_id}")
    #         return
    #     contacts = contacts.get(chain, {})
    #     if not contacts:
    #         print(f"No contacts found for chain {chain} in structure {structure_id}")
    #         return
    #     for ctype, residues in contacts.items():
    #         for pdb_residue in residues:
    #             unp_residue = aligner.map(chain, self.uniprot_id, pdb_residue, remove_mismatch=True)       
    #             if unp_residue is not None:
    #                 contact_info = {
    #                     "residue": unp_residue,
    #                     "pdb_residue": pdb_residue,
    #                     "chain": chain,   
    #                     "structure_id": structure_id,
    #                 }
    #                 self.contacts[ctype].append(contact_info)

    # def build_aligner(self, structure_id: str) -> Structure_Sequence_Aligner | None:
    #     """Build the PDB -> UniProt residue mapper for a single structure."""
    #     try:
    #         return Structure_Sequence_Aligner(self.structure_db, structure_id)
    #     except Exception as e:
    #         warnings.warn(f"Could not build aligner for {structure_id}: {e}")
    #         return None
    
    # def get_interfaces_by_type(self, type: Literal[*INTERFACE_TYPES] = "inter-chain", count: bool = False, full_contact_info: bool = False ) -> list[dict]:
    #     if full_contact_info:
    #         contacts = self.contacts.get(type, [])
    #     else:
    #         contacts = [tuple(contact["residue"]) for contact in self.contacts.get(type, [])]
    #         contacts = list(set(contacts))
    #     if count:
    #         return len(contacts)
    #     return contacts
    

    # def lysine_nz_proportion_in_structures(self) -> list[str]:
    #     """Return the proportion of lysine having zeta-nitrogen in the uniprot_id associated chains per structure."""
    #     lysine_nz_proportions = defaultdict(float)
    #     structure_ids = self.structure_ids
    #     if len(structure_ids) == 0:
    #         return np.nan, np.nan, {}
    #     for structure_id in self.structure_ids:
    #         aligner = self.build_aligner(structure_id)
    #         if aligner is None:
    #             warnings.warn(f"Could not build aligner for structure {structure_id}")
    #             continue
    #         chains = aligner.chains_containing_uniprot(self.uniprot_id)
    #         nz_prop_per_chain = Structure(structure_id, database=self.structure_db).lysine_nz_proportion(chains=chains)
    #         if not nz_prop_per_chain:
    #             warnings.warn(f"input chains not present in structure {structure_id} for uniprot id {self.uniprot_id}")
    #             continue
    #         #assert set(nz_prop_per_chain.keys()) == set(chains), f"chains returned are not same as input chains"
    #         if set(nz_prop_per_chain.keys()) != set(chains):
    #             warnings.warn(f"chains returned are not same as input chains for structure {structure_id} for uniprot id {self.uniprot_id}")
    #         nz_prop = sum(nz_prop_per_chain.values()) / len(chains)  # average nz proportion per chain
    #         lysine_nz_proportions[structure_id] = nz_prop
    #     if len(lysine_nz_proportions) == 0:
    #         warnings.warn(f"No lysine nz proportions found for uniprot id {self.uniprot_id}")
    #         return np.nan, np.nan, {}
    #     average_nz_prop = sum(lysine_nz_proportions.values()) / len(lysine_nz_proportions)
    #     max_nz_prop = max(lysine_nz_proportions.values())
    #     return average_nz_prop, max_nz_prop, lysine_nz_proportions

    # @classmethod
    # def factory_method(cls, uniprot_id: str, 
    #                    structure_db: str, 
    #                    cutoff: float = 8.0, 
    #                    lysine_only: bool = False,
    #                    interfaces: dict[str, list[dict]] = None) -> BindingInterfaces:
    #     binding_interfaces = BindingInterfaces(uniprot_id, structure_db, cutoff, lysine_only)
    #     if interfaces is not None:
    #         binding_interfaces.set_interfaces(interfaces)
    #     else:
    #         binding_interfaces.build_interfaces()
    #     return binding_interfaces


# if __name__ == "__main__":
#     import random

#     structure_db = "pdb"
#     #structure_db = "alphafold"
#     cutoff = 6.5
#     lysine_only = True
#     refresh_interval = 10
#     save = True
#     # Whether to process uniprot ids with no interfaces.
#     process_empty_uniprot_ids = False   

#     binding_interfaces_store = BindingInterfacesStore(structure_db, cutoff, lysine_only)

#     # UniProt ids that still need processing: all ids minus those already saved.
#     uniprot_ids_to_process = binding_interfaces_store.unsaved_uniprot_ids(ALL_UNIPROT_IDS)
#     empty_uniprot_ids = binding_interfaces_store.empty_uniprot_ids()
#     if process_empty_uniprot_ids:
#         uniprot_ids_to_process = uniprot_ids_to_process + empty_uniprot_ids
#     print(f"{len(uniprot_ids_to_process)} uniprot ids to process")
#     saved_count = 0
#     # A single lightweight extractor reused across iterations.
#     while uniprot_ids_to_process:
#         # Pick a random uniprot id so that parallel threads are unlikely to collide.
#         uniprot_id = uniprot_ids_to_process.pop(random.randrange(len(uniprot_ids_to_process)))
#         # The store re-checks if it is already saved (possibly by another thread) and only
#         # extracts when needed; returns None if it was already saved.
#         if save:
#             interfaces = binding_interfaces_store.extract_and_save(uniprot_id)
#         else:
#             interfaces = BindingInterfaces.factory_method(uniprot_id, structure_db, cutoff, lysine_only)
#         if interfaces is None:
#             continue
#         saved_count += 1
#         if not save:
#             seq = UniprotToSequence.get_sequence(uniprot_id)
#             structure_mapper = UniprotStructureMappingsStore.get_mapper(structure_db=structure_db)
#             structure_ids = structure_mapper.structures_for_uniprot_id(uniprot_id)
#             print(f"{uniprot_id}:")
#             print(f"Sequence length: {len(seq)}")
#             print(f"{len(structure_ids)} structures")
#             print(f"{interfaces.get_interfaces_by_type("all", count=True)} all interfaces")
#             print(f"{interfaces.get_interfaces_by_type("inter-chain", count=True)} inter-chain interfaces")
#             print(f"{interfaces.get_interfaces_by_type("intra-chain", count=True)} intra-chain interfaces")
#         # Periodically drop ids that were saved by other threads.
#         if saved_count % refresh_interval == 0:
#             uniprot_ids_to_process = binding_interfaces_store.unsaved_uniprot_ids(uniprot_ids_to_process)
#             if process_empty_uniprot_ids:
#                 uniprot_ids_to_process = uniprot_ids_to_process + empty_uniprot_ids
#             print(f"Refreshed queue: {len(uniprot_ids_to_process)} uniprot ids remaining")
        
