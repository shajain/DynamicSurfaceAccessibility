from pathlib import Path, PosixPath
from typing_extensions import List
import gemmi
from collections import defaultdict
#from dsa.binding_interfaces.config import PDB_DIR, test_structure_id, BI_STRUCTURE_DIR, STRUCTURE_DIR
#from dsa.binding_interfaces.structure_uniprot_alignment.interval_dict import IntervalDict
#from dsa.binding_interfaces.process_structure.interface_extraction import BindingInterfaceExtractor
#from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappingsStore
#from dsa.binding_interfaces.mappings.uniprot_to_alphaFold import get_alphafold_structure, uniprot_ids_with_saved_alphafold_structures
#from dsa.uniprot.config import ALL_UNIPROT_IDS
import pdb
import json
import random
import warnings
#from dsa.binding_interfaces.process_structure.sequence_aligner import SequenceAligner
from Bio.SeqUtils import seq1
import numpy as np
import re
# from dsa.misc.store import Store
# from dsa.misc.file_names import FileNames
from dsa.binding_interfaces.structure.structure_store import StructureFileStore
from abc import ABC, abstractmethod
from dsa.binding_interfaces.structure_uniprot_alignment.aligner import StructureSequenceAligner
from dsa.binding_interfaces.structure.structure_model import StructureModel
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure_uniprot_alignment.pdb_aligner import PDBSequenceAligner
from dsa.binding_interfaces.structure_uniprot_alignment.alphafold_aligner import AlphaFoldSequenceAligner
#from dsa.binding_interfaces.structure_uniprot_alignment.aligner import StructureSequenceAligner


class Structure(ABC):
    def __init__(self, structure_id: str, database: str):
        self.structure_id = structure_id
        self.database = database
        self.store = ClassFactory.create(interface=StructureFileStore, key=database)
        self.structure = self.store.load(structure_id)
        self.structure.setup_entities()
        self.structure.add_entity_ids()
        self.uniprot_to_entities = self.map_uniprot_to_entities()
        self.model = StructureModel(self.structure[0])
        # self.entity_to_subchains = self.entity_to_subchains()
        self.bio_model = StructureModel(self.create_biological_model())
        self.uniprot_aligner = None
        # self.build_sequence_aligners()

    # @property
    # def entities(self) -> list[str]:
    #     return list(self.entity_to_subchains.keys())

    @property
    def entities(self) -> list[str]:
        return [e.name for e in self.structure.entities]

    @property
    def chains(self) -> list[str]:
        return self.model.chains

    def setup_aligner_if_none(self):
        if self.uniprot_aligner is None:
            self.uniprot_aligner = ClassFactory.create(interface=StructureSequenceAligner, 
                                                        key=self.database, 
                                                        structure_id = self.structure_id,
                                                        structure_model = self.model, 
                                                        uniprot_to_entity = self.uniprot_to_entities)
        
        

    def is_chain_aligned_with_uniprot(self, chain: str, uniprot_id: str) -> bool:
        self.setup_aligner_if_none()
        return self.uniprot_aligner.is_chain_aligned(chain, uniprot_id)
    
    def map_residue_to_uniprot(self, chain: str, residue: str, uniprot_id: str=None) -> tuple[int, str]:
        self.setup_aligner_if_none()
        return self.uniprot_aligner.map(chain, residue, uniprot_id)

    def uniprot_to_chains(self, uniprot_id: str) -> list[str]:
        self.setup_aligner_if_none()
        return self.uniprot_aligner.uniprot_to_chains(uniprot_id)

    def biochain_to_chain(self, bio_chain: str) -> str:
        bio_subchains = self.bio_model.get_subchains_in_chain(bio_chain)
        chains = []
        for bio_subchain in bio_subchains:
            subchain = re.sub(r'\d+', '', bio_subchain)
            chain = self.model.get_chain_containing_subchain(subchain)
            chains.append(chain)
        assert len(set(chains)) == 1, f"Multiple chains found for bio_chain {bio_chain}"
        return chains[0]

    @property
    def gemmi_structure(self) -> gemmi.Structure:
        return self.structure

    @property
    def gemmi_bio_model(self) -> gemmi.Model:
        return self.bio_model.model

    # @property
    # def subchains(self) -> list[str]:
    #     return [chain for chains in self.entity_to_subchains.values() for chain in chains]

    # @property
    # def protein_subchains(self) -> list[str]:
    #     chains = [c for entities in self.uniprot_to_entity.values() 
    #                 for entity in entities 
    #                 for c in self.entity_to_subchains[entity]]
    #     return list(set(chains))

    # @property
    # v1
    # def asym_to_author_chain(self) -> dict[str, str]:
    #     mapping = {}
    #     for chain in self.structure[0]:
    #         # subchains() queries the underlying C++ structure directly
    #         for sub in chain.subchains():
    #             if sub and sub[0].subchain:
    #                 mapping[sub[0].subchain] = chain.name
    #     return mapping

    # def asym_to_author_chain(self) -> dict[str, str]:
    # v1
    #     mapping = {}
    #     for chain in self.structure[0]:
    #         for residue in chain:
    #             if residue.subchain:
    #                 mapping[residue.subchain] = chain.name
    #                 break
    #     return mapping

    # @property
    # v1
    # def bio_to_author_chain(self) -> dict[str, str]:
    #     mapping = {}
    #     asym_to_author = self.asym_to_author_chain
    #     for chain in self.bio_structure:
    #         for residue in chain:
    #             if residue.subchain and residue.entity_type == gemmi.EntityType.Polymer:
    #                 asym_id = re.sub(r'\d+', '', residue.subchain)
    #                 mapping[chain.name] = asym_to_author[asym_id]
    #                 break
    #     return mapping

    # @property
    # v2
    # def bio_to_author_chain(self) -> dict[str, str]:
    #     mapping = {}
    #     asym_to_author = self.asym_to_author_chain
    #     for chain in self.bio_structure:
    #         for sub in chain.subchains():
    #             if sub and sub[0].subchain:
    #                 asym_id = re.sub(r'\d+', '', sub[0].subchain)
    #                 mapping[chain.name] = asym_to_author[asym_id]
    #                 break
    #     return mapping

    # @property 
    # def entity_to_uniprot(self) -> dict[str, str]:
    #     e2u = defaultdict(str)
    #     for unp, entities in self.uniprot_to_entity.items():
    #         for entity in entities:
    #             e2u[entity] = unp
    #     return e2u
        
    @property
    def uniprot_ids(self) -> list[str]:
        return list(self.uniprot_to_entity.keys())

    # @property
    # def gemmi_structure(self) -> gemmi.Structure:
    #     return self.structure

    # @property
    # def biological_structure(self) -> gemmi.Structure:
    #     return self.bio_structure

    def create_biological_model(self, assembly_index: int = 0):
        if not self.structure.assemblies:
            #print("No biological assembly annotations found, using structure as-is.")
            return self.structure[0]
        assembly = self.structure.assemblies[assembly_index]
        #print(f"Using biological assembly: '{assembly.name}'")
        bio_model = gemmi.make_assembly(assembly, self.structure[0], gemmi.HowToNameCopiedChain.AddNumber)
        return bio_model

    def get_bio_chain_residue_atom_objects(self, chain: str, 
                                            residue: str, 
                                            atom: str) -> tuple[gemmi.Chain, gemmi.Residue, gemmi.Atom]:
        return self.bio_model.get_chain_reisdue_atom_objects(chain, residue, atom)

    # @property
    # def entity_to_subchains(self) -> dict[str, list[str]]:
    #     mapping = defaultdict(list)
    #     for ent in self.structure.entities:
    #         mapping[ent.name].extend(list(ent.subchains))
    #     return mapping


    # def map_residue_to_uniprot(self, chain: str, residue: gemmi.Residue) -> tuple[int, str]:
    #     entity = self.residue_entity(residue)
    #     return self.entity_to_uniprot[entity.name]

    # @abstractmethod
    # def map_uniprot_to_internal_ids(self, id_type: str = 'entity_id') -> dict[str, list[str]]:
    #     pass

    @abstractmethod
    def map_uniprot_to_entities(self) -> dict[str, list[str]]:
        pass

    # def chain_sequence_grouped_by_uniprot(self) -> dict[str, dict[str, dict[str, str]]]:
    #     chain_sequences: dict[str, dict[str, dict[str, str]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    #     for chain in self.structure[0]:
    #         for res in chain:
    #             if res.entity_type != gemmi.EntityType.Polymer:
    #                 continue
    #             uniprot_id = self.entity_to_uniprot[res.entity_id]
    #             res_seq_id, res_name = str(res.seqid), res.name
    #             chain_sequences[uniprot_id][chain.name][res_seq_id] = res_name
    #     return chain_sequences

    # def lysine_nz_proportion(self, chains: list[str]) -> float:
    #     nz_lysines = defaultdict(float)
    #     total_lysines = defaultdict(int)
    #     nz_proportions = defaultdict(float)
    #     for chain in self.structure[0]:
    #         if chain.name not in chains:
    #             continue
    #         for res in chain:
    #             if res.name in LYSINE_NAMES:
    #                 total_lysines[chain.name] += 1
    #                 for atom in res:
    #                     if atom.name == 'NZ':
    #                         nz_lysines[chain.name] += 1
    #                         break
    #         if total_lysines[chain.name] == 0:
    #             warnings.warn(f"No lysines found in chain {chain.name}")
    #             nz_proportions[chain.name] = 0.0
    #         else:
    #             nz_proportions[chain.name] = nz_lysines[chain.name] / total_lysines[chain.name]
    #     if set(nz_proportions.keys()) != set(chains):
    #         warnings.warn(f"Some chains are not present in the structure")
    #     return nz_proportions

    # def residue_entity(self, residue: gemmi.Residue) -> gemmi.Entity:
    #     entity = next((e for e in self.structure.entities if e.name == residue.entity_id), None)
    #     return entity

    # def residue_polymer_type(self, res: gemmi.Residue) -> gemmi.PolymerType:
    #     entity = self.residue_entity(res)
    #     return entity.polymer_type if entity is not None else None

    def __getstate__(self):
        state = {"structure_id": self.structure_id}
        return state

    def __setstate__(self, state):
        self.__init__(state["structure_id"])

    # @classmethod
    # def factory(cls, structure_id: str, database: str):
    #     if database == "pdb":
    #         return PDBStructure(structure_id)
    #     elif database == "alphafold":
    #         return AlphaFoldStructure(structure_id)
    #     else:
    #         raise ValueError(f"Invalid database: {database}")


    # @classmethod
    # def extract_position_and_name(cls, residue: gemmi.Residue) -> tuple[str, str]:
    #     if cls.residue_representation['position'] == 'int_as_str':
    #         res_num = str(residue.seqid.num)
    #     elif cls.residue_representation['position'] == 'alphanumeric':
    #         res_num = str(residue.seqid)
    #     else:
    #         raise ValueError(f"Invalid residue representation: {cls.residue_representation['position']}")
    #     if cls.residue_representation['name'] == 'three_letter':
    #         res_name = residue.name
    #     elif cls.residue_representation['name'] == 'single_letter':
    #         info = gemmi.find_tabulated_residue(residue.name)
    #         res_name = info.one_letter_code if info is not None else 'X'
    #     else:
    #         raise ValueError(f"Invalid residue representation: {cls.residue_representation['name']}")
    #     return res_num, res_name


@ClassFactory.register(key="pdb", interface=Structure)
class PDBStructure(Structure):
    def __init__(self, structure_id: str):
        super().__init__(structure_id, "pdb")

    # def map_uniprot_to_internal_ids(self, id_type: str = 'entity_id') -> dict[str, list[str]]:
    #     #id_type can be 'entity_id' or 'id'. 'id' is the ref_id in the _struct_ref_seq table.
    #     uniprot_to_entity = defaultdict(list)
    #     block = self.store.load_block(self.structure_id)
    #     ID_table = block.find('_struct_ref.', [id_type, 'db_name', 'pdbx_db_accession'])
    #     for row in ID_table:
    #         if row[1].upper() == 'UNP':
    #             uniprot_to_entity[row[2]].append(row[0])            
    #     return uniprot_to_entity

    def map_uniprot_to_entities(self) -> dict[str, list[str]]:
        uniprot_to_entity = defaultdict(list)
        block = self.store.load_block(self.structure_id)
        ID_table = block.find('_struct_ref.', ['entity_id', 'db_name', 'pdbx_db_accession'])
        for row in ID_table:
            if row[1].upper() == 'UNP':
                uniprot_to_entity[row[2]].append(row[0])            
        return uniprot_to_entity

    def __getstate__(self):
        state = super().__getstate__()
        return state
    

@ClassFactory.register(key="alphafold", interface=Structure)
class AlphaFoldStructure(Structure):
    def __init__(self, structure_id: str):
        super().__init__(structure_id, "alphafold")

    def map_uniprot_to_entities(self) -> dict[str, list[str]]:
        uniprot_to_entities = {self.structure_id: self.entities}
        return uniprot_to_entities

    def __getstate__(self):
        state = super().__getstate__()
        return state

    
    


    # @staticmethod
    # def select_unprocessed_structures(structure_ids: list[str], structure_db: str="pdb", cutoff: float=8.0, lysine_only: bool = False) -> list[str]:
    #     if structure_db == 'pdb':   
    #         structure_ids = [st.lower() for st in structure_ids] 
    #     saved_structures = Structure.get_saved_structures(structure_db, from_structures=structure_ids)
    #     eligible_structures = list(set(structure_ids) & set(saved_structures))
    #     processed_structures = Structure.processed_structures(structure_db, cutoff, lysine_only)
    #     return list(set(eligible_structures) - set(processed_structures))

    # @staticmethod
    # def save_binding_interfaces(contact_residues: dict[str, IntervalDict], structure_id: str, structure_db: str, cutoff: float, lysine_only: bool = False):
    #     file_path = Structure.BI_path(structure_id, structure_db, cutoff, lysine_only)
    #     Path(file_path).parent.mkdir(parents=True, exist_ok=True)
    #     with open(file_path, "w") as f:
    #         json.dump(contact_residues, f)

    # @staticmethod
    # def load_binding_interfaces(structure_id: str, structure_db: str, cutoff: float, lysine_only: bool = False):
    #     is_saved = Structure.is_saved_binding_interfaces(structure_id, structure_db, cutoff, lysine_only)
    #     if is_saved:
    #         file_path = Structure.BI_path(structure_id, structure_db, cutoff, lysine_only)
    #         with open(file_path, "r") as f:
    #             return json.load(f)
    #     return None

    # @staticmethod
    # def get_structure_path(structure_id:str, structure_db:str) -> Path:
    #     if structure_db == 'pdb':
    #         structure_id = structure_id.lower()
    #         if not structure_id.endswith('.cif'):
    #             structure_id = structure_id + '.cif'
    #         return PDB_DIR / structure_id
    #     elif structure_db == 'alphafold':
    #         return get_alphafold_structure(uniprot_id=structure_id, format="pdb")
    #     ValueError(f"Invalid database: {structure_db}")
    #     return 

    # @staticmethod
    # def is_saved_structure(structure_id:str, structure_db:str) -> bool:
    #     file_path = Structure.get_structure_path(structure_id, structure_db)
    #     return file_path.exists()

    # @staticmethod
    # def is_saved_binding_interfaces(structure_id:str, structure_db:str, cutoff:float, lysine_only: bool = False) -> bool:
    #     file_path = Structure.BI_path(structure_id, structure_db, cutoff, lysine_only)
    #     return file_path.exists()

    # @staticmethod
    # def processed_structures(structure_db: str, cutoff: float, lysine_only: bool = False) -> list[str]:
    #     lysine_str = "lysine_only" if lysine_only else ""
    #     path = BI_STRUCTURE_DIR / "_".join([structure_db, str(cutoff),lysine_str])
    #     return [file.name.split("_")[0] for file in path.glob("*.json")]
            
    # @staticmethod
    # def BI_path(structure_id: str, structure_db: str, cutoff: float, lysine_only: bool = False) -> Path:
    #     db = structure_db
    #     id = structure_id
    #     lysine_str = "lysine_only" if lysine_only else ""
    #     path = BI_STRUCTURE_DIR / "_".join([db, str(cutoff),lysine_str])/ ("_".join([id, db, str(cutoff), lysine_str]) + ".json")
    #     return path

    # def __str__(self):
    #     string = f"{self.database}-{self.id}"
    #     return string

    
    # @staticmethod
    # def canonicalize_structure_id(structure_id: str, structure_db: str) -> str:
    #     if structure_db == "pdb":
    #         return structure_id.lower()
    #     return structure_id



