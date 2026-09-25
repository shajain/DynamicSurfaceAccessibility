from pickle import FALSE
from dsa.binding_interfaces.config import ALL_UNIPROT_IDS
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappingsStore
from dsa.binding_interfaces.process_structure.Structure import Structure
from dsa.binding_interfaces.structure_alignment.interface import Structure_Sequence_Aligner
from typing import Literal
from dsa.binding_interfaces.config import BI_PROTEIN_DIR
from Bio.SeqUtils import seq1
from pathlib import Path
import json
import warnings
from collections import defaultdict
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
import numpy as np

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


class BindingInterfacesStore:
    """Persists aggregated, UniProt wise binding interfaces.

    The store is associated with a particular `structure_db` and `cutoff`. One JSON file is
    written per UniProt id under BI_PROTEIN_DIR. The file name is `{uniprot_id}.json`. Files
    are nested in a `{structure_db}_{cutoff}` sub directory so that interfaces computed with
    different databases/cutoffs do not collide while keeping the requested `{uniprot_id}.json`
    file name.
    """

    def __init__(self, structure_db: str = "pdb", cutoff: float = 8.0, lysine_only: bool = False):
        self.structure_db = structure_db
        self.cutoff = cutoff
        self.lysine_only = lysine_only
    
    @property
    def directory(self) -> Path:
        lysine_str = "lysine_only" if self.lysine_only else ""
        return BI_PROTEIN_DIR / "_".join([self.structure_db, str(self.cutoff), lysine_str])

    def file_path(self, uniprot_id: str) -> Path:
        return self.directory / f"{uniprot_id}.json"

    def is_saved(self, uniprot_id: str) -> bool:
        return self.file_path(uniprot_id).exists()

    def save(self, uniprot_id: str, binding_interfaces: BindingInterfaces) -> Path:
        #interfaces = binding_interfaces.get_interfaces()
        file_path = self.file_path(uniprot_id)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "w") as f:
            json.dump(binding_interfaces.__dict__, f, indent=4)
        return file_path

    def load(self, uniprot_id: str) -> BindingInterfaces | None:
        file_path = self.file_path(uniprot_id)
        if not file_path.exists():
            return None
        with open(file_path, "r") as f:
            try:
                object_dict = json.load(f)
            except json.JSONDecodeError:
                warnings.warn(f"JSONDecodeError for {file_path}")
                return None
            binding_interfaces = BindingInterfaces.factory_method(**object_dict)
            return binding_interfaces

    def saved_uniprot_ids(self) -> list[str]:
        """Look at all the files in the directory and return the uniprot ids of the files that are saved."""
        directory = self.directory
        uniprot_ids = [file.stem for file in directory.glob("*.json")]
        return uniprot_ids

    def unsaved_uniprot_ids(self, uniprot_ids: list[str]) -> list[str]:
        """Select the UniProt ids from the list whose interfaces are not saved yet."""
        return [uid for uid in uniprot_ids if not self.is_saved(uid)]

    def corrupted_uniprot_ids(self) -> list[str]:
        saved_uniprot_ids = self.saved_uniprot_ids()
        corrupted_uniprot_ids = [id for id in saved_uniprot_ids if self.load(id) is None]
        return corrupted_uniprot_ids

    def empty_uniprot_ids(self) -> list[str]:
        """Look at all the files in the directory and return the uniprot ids with no interfaces."""
        uniprot_ids = self.saved_uniprot_ids()
        empty_uniprot_ids = []
        for unp in uniprot_ids:
            binding_interfaces = self.load(unp)
            if binding_interfaces is not None and binding_interfaces.is_empty():
                    empty_uniprot_ids.append(unp)
        #empty_uniprot_ids = [unp for unp in uniprot_ids if self.load(unp).is_empty()]
        return empty_uniprot_ids
    
    def remove_corrupted_uniprot_ids(self) -> list[str]:
        corrupted_uniprot_ids = self.corrupted_uniprot_ids()
        for unp in corrupted_uniprot_ids:
            self.delete(unp)
        return corrupted_uniprot_ids

    def delete(self, uniprot_id: str) -> None:
        file_path = self.file_path(uniprot_id)
        if file_path.exists():
            file_path.unlink()
        else:
            warnings.warn(f"File {file_path} does not exist")

    def extract_and_save(self, uniprot_id: str) -> dict[str, list[dict]] | None:
        """Extract and save the binding interfaces for a UniProt id unless they are already saved.

        Returns the extracted interfaces, or None if they were already saved (and not overwritten).
        A `BindingInterfaces` extractor can be supplied to be reused across many calls.
        """
        binding_interfaces = BindingInterfaces.factory_method(uniprot_id, self.structure_db, self.cutoff, self.lysine_only)
        if binding_interfaces is not None:
            self.save(uniprot_id, binding_interfaces)
        else:
            warnings.warn(f"Could not extract binding interfaces for {uniprot_id}")
        return binding_interfaces


class BindingInterfaces:
    """Extracts and aggregates binding interfaces for UniProt ids across all their structures.

    For each UniProt id the structure ids are looked up through `UniprotStructureMappings`,
    the saved per-structure binding interfaces are loaded through `Structure`, and the
    contacting residues are mapped to the UniProt sequence using `PDB_Uniprot_Mapper` from
    `structure_alignment`. Residues are aggregated per interface type as a list of
    dictionaries with the keys in `RESIDUE_KEYS`.
    """

    INTERFACE_TYPES = INTERFACE_TYPES

    def __init__(self, uniprot_id: str,
                 structure_db: str = "pdb",
                 cutoff: float = 8.0,
                 lysine_only: bool = False):
        self.structure_db = structure_db
        self.cutoff = cutoff
        self.lysine_only = lysine_only
        self.uniprot_id = uniprot_id
        self.interfaces: dict[str, list[dict]] = {ctype: [] for ctype in INTERFACE_TYPES}

    @property
    def structure_ids(self) -> list[str]:
        structure_mapper = UniprotStructureMappingsStore.get_mapper(structure_db=self.structure_db)
        structure_ids = structure_mapper.structures_for_uniprot_id(self.uniprot_id)
        return structure_ids

    def is_empty(self) -> bool:
        return all(len(interfaces) == 0 for interfaces in self.interfaces.values())

    def build_interfaces(self) -> dict[str, list[dict]]:
        structure_mapper = UniprotStructureMappingsStore.get_mapper(structure_db=self.structure_db)
        structure_ids = structure_mapper.structures_for_uniprot_id(self.uniprot_id)
        # Process one structure at a time; the mapper is built and discarded per structure
        # so that memory stays bounded when running many extractions in parallel.
        for structure_id in structure_ids:
            self.build_from_structure(structure_id)

    def get_interfaces(self) -> dict[str, list[dict]]:
        return self.interfaces

    def set_interfaces(self, interfaces: dict[str, list[dict]]):
        self.interfaces = interfaces

    def build_from_structure(self, structure_id: str):
        """Map and aggregate the binding-interface residues of a single structure for a UniProt id."""
        contacts = self.structure_contacts(structure_id)
        if contacts is None:
            warnings.warn(f"Could not load contacts for structure {structure_id}")
            return
        all_contacts = [len(contacts[c]["all"]) for c in contacts.keys()]
        if np.max(all_contacts) > 0:
            print(f"Structure has contacts for uniprot id {self.uniprot_id}")
        aligner = self.build_aligner(structure_id)
        if aligner is None:
            warnings.warn(f"Could not build PDB_Uniprot_Mapper for structure {structure_id}")
            return
        if contacts is None:
            warnings.warn(f"Could not load contacts for structure {structure_id}")
            return
        chains = aligner.chains_containing_uniprot(self.uniprot_id)
        for chain in chains:   
            self.build_from_chain(contacts, aligner, chain, structure_id)

    def build_from_chain(self, contacts: dict[str, dict[str, list]], aligner: Structure_Sequence_Aligner, chain: str, structure_id: str):
        if not aligner.is_alignment_successful(chain, self.uniprot_id):
            warnings.warn(f"Alignment not successful for chain {chain} in structure {structure_id} with uniprot id {self.uniprot_id}")
            return
        contacts = contacts.get(chain, {})
        if not contacts:
            print(f"No contacts found for chain {chain} in structure {structure_id}")
            return
        for ctype, residues in contacts.items():
            for pdb_residue in residues:
                unp_residue = aligner.map(chain, self.uniprot_id, pdb_residue, remove_mismatch=True)       
                if unp_residue is not None:
                    contact_info = {
                        "residue": unp_residue,
                        "pdb_residue": pdb_residue,
                        "chain": chain,   
                        "structure_id": structure_id,
                    }
                    self.interfaces[ctype].append(contact_info)

    def structure_contacts(self, structure_id: str) -> dict[str, dict[str, list]] | None:
        """Load the saved per-structure binding interfaces (chain -> ctype -> [(pdb_res_num, pdb_res_name)])."""
        try:
            return Structure.extract_binding_interfaces(structure_id, database=self.structure_db,
                                                        cutoff=self.cutoff, lysine_only=self.lysine_only, only_from_saved_file=True)
        except Exception as e:
            warnings.warn(f"Could not load saved binding interfaces for {structure_id}: {e}")
            return None

    def build_aligner(self, structure_id: str) -> Structure_Sequence_Aligner | None:
        """Build the PDB -> UniProt residue mapper for a single structure."""
        try:
            return Structure_Sequence_Aligner(self.structure_db, structure_id)
        except Exception as e:
            warnings.warn(f"Could not build aligner for {structure_id}: {e}")
            return None
    
    def get_interfaces_by_type(self, type: Literal[*INTERFACE_TYPES] = "inter-chain", count: bool = False, full_contact_info: bool = False ) -> list[dict]:
        if full_contact_info:
            contacts = self.interfaces.get(type, [])
        else:
            contacts = [tuple(contact["residue"]) for contact in self.interfaces.get(type, [])]
            contacts = list(set(contacts))
        if count:
            return len(contacts)
        return contacts

    def lysine_nz_proportion_in_structures(self) -> list[str]:
        """Return the proportion of lysine having zeta-nitrogen in the uniprot_id associated chains per structure."""
        lysine_nz_proportions = defaultdict(float)
        structure_ids = self.structure_ids
        if len(structure_ids) == 0:
            return np.nan, np.nan, {}
        for structure_id in self.structure_ids:
            aligner = self.build_aligner(structure_id)
            if aligner is None:
                warnings.warn(f"Could not build aligner for structure {structure_id}")
                continue
            chains = aligner.chains_containing_uniprot(self.uniprot_id)
            nz_prop_per_chain = Structure(structure_id, database=self.structure_db).lysine_nz_proportion(chains=chains)
            if not nz_prop_per_chain:
                warnings.warn(f"input chains not present in structure {structure_id} for uniprot id {self.uniprot_id}")
                continue
            #assert set(nz_prop_per_chain.keys()) == set(chains), f"chains returned are not same as input chains"
            if set(nz_prop_per_chain.keys()) != set(chains):
                warnings.warn(f"chains returned are not same as input chains for structure {structure_id} for uniprot id {self.uniprot_id}")
            nz_prop = sum(nz_prop_per_chain.values()) / len(chains)  # average nz proportion per chain
            lysine_nz_proportions[structure_id] = nz_prop
        if len(lysine_nz_proportions) == 0:
            warnings.warn(f"No lysine nz proportions found for uniprot id {self.uniprot_id}")
            return np.nan, np.nan, {}
        average_nz_prop = sum(lysine_nz_proportions.values()) / len(lysine_nz_proportions)
        max_nz_prop = max(lysine_nz_proportions.values())
        return average_nz_prop, max_nz_prop, lysine_nz_proportions

    @classmethod
    def factory_method(cls, uniprot_id: str, 
                       structure_db: str, 
                       cutoff: float = 8.0, 
                       lysine_only: bool = False,
                       interfaces: dict[str, list[dict]] = None) -> BindingInterfaces:
        binding_interfaces = BindingInterfaces(uniprot_id, structure_db, cutoff, lysine_only)
        if interfaces is not None:
            binding_interfaces.set_interfaces(interfaces)
        else:
            binding_interfaces.build_interfaces()
        return binding_interfaces


if __name__ == "__main__":
    import random

    structure_db = "pdb"
    #structure_db = "alphafold"
    cutoff = 6.5
    lysine_only = True
    refresh_interval = 10
    save = True
    # Whether to process uniprot ids with no interfaces.
    process_empty_uniprot_ids = False   

    binding_interfaces_store = BindingInterfacesStore(structure_db, cutoff, lysine_only)

    # UniProt ids that still need processing: all ids minus those already saved.
    uniprot_ids_to_process = binding_interfaces_store.unsaved_uniprot_ids(ALL_UNIPROT_IDS)
    empty_uniprot_ids = binding_interfaces_store.empty_uniprot_ids()
    if process_empty_uniprot_ids:
        uniprot_ids_to_process = uniprot_ids_to_process + empty_uniprot_ids
    print(f"{len(uniprot_ids_to_process)} uniprot ids to process")
    saved_count = 0
    # A single lightweight extractor reused across iterations.
    while uniprot_ids_to_process:
        # Pick a random uniprot id so that parallel threads are unlikely to collide.
        uniprot_id = uniprot_ids_to_process.pop(random.randrange(len(uniprot_ids_to_process)))
        # The store re-checks if it is already saved (possibly by another thread) and only
        # extracts when needed; returns None if it was already saved.
        if save:
            interfaces = binding_interfaces_store.extract_and_save(uniprot_id)
        else:
            interfaces = BindingInterfaces.factory_method(uniprot_id, structure_db, cutoff, lysine_only)
        if interfaces is None:
            continue
        saved_count += 1
        if not save:
            seq = UniprotToSequence.get_sequence(uniprot_id)
            structure_mapper = UniprotStructureMappingsStore.get_mapper(structure_db=structure_db)
            structure_ids = structure_mapper.structures_for_uniprot_id(uniprot_id)
            print(f"{uniprot_id}:")
            print(f"Sequence length: {len(seq)}")
            print(f"{len(structure_ids)} structures")
            print(f"{interfaces.get_interfaces_by_type("all", count=True)} all interfaces")
            print(f"{interfaces.get_interfaces_by_type("inter-chain", count=True)} inter-chain interfaces")
            print(f"{interfaces.get_interfaces_by_type("intra-chain", count=True)} intra-chain interfaces")
        # Periodically drop ids that were saved by other threads.
        if saved_count % refresh_interval == 0:
            uniprot_ids_to_process = binding_interfaces_store.unsaved_uniprot_ids(uniprot_ids_to_process)
            if process_empty_uniprot_ids:
                uniprot_ids_to_process = uniprot_ids_to_process + empty_uniprot_ids
            print(f"Refreshed queue: {len(uniprot_ids_to_process)} uniprot ids remaining")
        
