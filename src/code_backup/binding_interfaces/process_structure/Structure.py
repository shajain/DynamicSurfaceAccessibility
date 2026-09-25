from pathlib import Path, PosixPath
from typing_extensions import List
import gemmi
from collections import defaultdict
from dsa.binding_interfaces.config import PDB_DIR, test_structure_id, BI_STRUCTURE_DIR
from dsa.binding_interfaces.process_structure.IntervalDict import IntervalDict
from dsa.binding_interfaces.process_structure.interface_extraction import BindingInterfaceExtractor
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappings
import pdb
import json
import random
import warnings
#from dsa.binding_interfaces.process_structure.sequence_aligner import SequenceAligner
from Bio.SeqUtils import seq1
import numpy as np


class Structure:
    SavedStructures = defaultdict(list)
    residue_representation = {'position': 'int_as_str', 'name': 'three_letter'}
    # residue_representation = {'position': 'int_as_str', 'name': 'single_letter'}
    # residue_representation = {'position': 'alphanumeric', 'name': 'three_letter'}
    # residue_representation = {'position': 'alphanumeric', 'name': 'single_letter'}
    def __init__(self, structure_id: str, database: str = "pdb"):
        if database == "pdb":
            self.structure_dir = PDB_DIR
        else:   
            raise ValueError(f"Invalid database: {database}")
        self.structure_id = structure_id
        self.structure_path = Structure.get_structure_path(structure_id, database)
        if not self.structure_path.exists():
            raise FileNotFoundError(f"Structure file not found: {self.structure_path}")
        self.block = gemmi.cif.read(self.structure_path.as_posix()).sole_block()
        self.structure = gemmi.read_structure(self.structure_path.as_posix())  # or .cif.gz
        self.structure.setup_entities()
        self.uniprot_to_entity = self.extract_uniprot_to_internalIDs_mappings(id_type='entity_id')
        self.entity_to_chains = self.entity_to_chains()
        # self.sequence_aligners = defaultdict(SequenceAligner)
        # self.build_sequence_aligners()

    @property
    def entities(self) -> list[str]:
        return list(self.entity_to_chains.keys())

    @property
    def chains(self) -> list[str]:
        return [chain for chains in self.entity_to_chains.values() for chain in chains]

    @property
    def protein_chains(self) -> list[str]:
        chains = [c for entities in self.uniprot_to_entity.values() 
                    for entity in entities 
                    for c in self.entity_to_chains[entity]]
        return list(set(chains))
    

    @property 
    def entity_to_uniprot(self) -> dict[str, str]:
        e2u = defaultdict(str)
        for unp, entities in self.uniprot_to_entity.items():
            for entity in entities:
                e2u[entity] = unp
        return e2u
        
    @property
    def uniprotIds(self) -> list[str]:
        return list(self.uniprot_to_entity.keys())

    def initialize_binding_interfaces(self, cutoff: float = 8.0) -> BindingInterfaceExtractor:
        self.BI_extractor = BindingInterfaceExtractor(self.structure, self.block, cutoff=cutoff)
        return self.BI_extractor



    def extract_uniprot_to_internalIDs_mappings(self, id_type: str = 'entity_id') -> dict[str, list[str]]:
        #id_type can be 'entity_id' or 'id'. 'id' is the ref_id in the _struct_ref_seq table.
        uniprot_to_entity = defaultdict(list)
        ID_table = self.block.find('_struct_ref.', [id_type, 'db_name', 'pdbx_db_accession'])
        for row in ID_table:
            if row[1].upper() == 'UNP':
                uniprot_to_entity[row[2]].append(row[0])
        return uniprot_to_entity

    # def unaligned_chains(self) -> int:
    #     unaligned_chains = [chain for chain, aligner in self.sequence_aligners.items() if not aligner.successful_alignment]
    #     return unaligned_chains
    
    def entity_to_chains(self) -> dict[str, list[str]]:
        entity_to_chains = defaultdict(list)
        for ent in self.structure.entities:
            entity_to_chains[ent.name].extend(list(ent.subchains))
        return entity_to_chains

    
    # def build_sequence_aligners(self) -> dict[str, SequenceAligner]:
    #     sequence_aligner = defaultdict(IntervalDict)
    #     tab = self.block.find('_struct_ref_seq.',
    #                            ['pdbx_strand_id', 'pdbx_auth_seq_align_beg', 'pdbx_seq_align_beg_ins_code',
    #                             'pdbx_auth_seq_align_end', 'pdbx_seq_align_end_ins_code', 
    #                             #'seq_align_beg', 'seq_align_end',
    #                             'db_align_beg', 'pdbx_db_align_beg_ins_code', 'db_align_end', 'pdbx_db_align_end_ins_code'])
    #     for row in tab:
    #         if row[1].isdigit() and row[3].isdigit() and row[5].isdigit() and row[7].isdigit():
    #             beg, end, unp_beg, unp_end = int(row[1]), int(row[3]), int(row[5]), int(row[7])
    #         else:
    #             continue
    #         if row[2]!=row[6] or row[4]!=row[8]:
    #             print(f"insertion code is not identical for pdb and uniprot chains {row[0]}")
    #         offset = unp_beg - beg
    #         if end-beg != unp_end - unp_beg:
    #             print(f"length of PDB and UniProt segment is not equal for chain {row[0]}")
    #             offset = None
    #         sequence_aligner[row[0]].add(beg, end, offset) # PDB residue number -> UniProt residue number
    #     return sequence_aligner
    

    def chain_sequence_grouped_by_uniprot(self) -> dict[str, dict[str, dict[str, str]]]:
        chain_sequences: dict[str, dict[str, dict[str, str]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
        for chain in self.structure[0]:
            for res in chain:
                if res.entity_type != gemmi.EntityType.Polymer:
                    continue
                uniprot_id = self.entity_to_uniprot[res.entity_id]
                res_num, res_name = self.extract_position_and_name(res)
                chain_sequences[uniprot_id][chain.name][res_num] = res_name
        return chain_sequences

    # def build_sequence_aligners(self) -> dict[str, SequenceAligner]:
    #     for uniprot in self.uniprotIds:
    #         for entity in self.uniprot_to_entity[uniprot]:
    #             for chain in self.entity_to_chains[entity]:
    #                 if chain in self.sequence_aligners:
    #                     continue
    #                 chain_sequence = self.chain_sequence(chain)
    #                 self.sequence_aligners[chain] = SequenceAligner(chain_sequence, uniprot)
            
    

    def entitywise_contacts(self, entity: str) -> list[gemmi.Residue]:
        contacts = list[gemmi.Residue]
        for chain in self.entity_to_chains[entity]:
            residues = self.BI_extractor.chainwise_contacts(chain)
            contacts.extend(residues)
            contacts = list(set(contacts))
        return contacts

    def chainwise_contacts(self, res_id_type: str = 'full_seqid', res_name_type: str = None) -> dict[str, dict[str, list[tuple[str, str]]]]:
        contact_residues: dict[str, dict[str, list[tuple[str, str]]]] = defaultdict(lambda: defaultdict(list))
        for chain in self.protein_chains:
            cnt_res = self.BI_extractor.group_residues_by_contact_type(chain)
            for ctype in cnt_res.keys():
                for r in cnt_res[ctype]:
                    res_num, res_name = self.extract_position_and_name(r)
                    contact_residues[chain][ctype].append((res_num, res_name))
                contact_residues[chain][ctype] = list(set(contact_residues[chain][ctype]))
        return contact_residues


    # def uniprotwise_contact_residues(self, uniprots: list[str]|str = None) -> dict[str, dict[str, list[int]]]:
    #     contact_residues: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))
    #     uniprots = self.uniprotIds if uniprots is None else [uniprots] if isinstance(uniprots, str) else uniprots
    #     for uniprot in uniprots:
    #         #Pooling residues across all chains of a uniprot
    #         for entity in self.uniprot_to_entity[uniprot]:
    #             for chain in self.entity_to_chains[entity]:
    #                 #Grouping residues by contact type
    #                 cnt_res = self.BI_extractor.group_residues_by_contact_type(chain)
    #                 #Mapping function to convert PDB residue numbers to UniProt sequence numbers
    #                 seq_aligner = self.sequence_aligners[chain]
    #                 if seq_aligner.successful_alignment:
    #                     for ctype in cnt_res.keys():
    #                         cr_tuples = []
    #                         for r in cnt_res[ctype]:
    #                             seq_num, ref_aa, _ = seq_aligner.map(r.seqid.num)
    #                             cr_tuples.append((seq_num, ref_aa, r.name, chain))
    #                         contact_residues[uniprot][ctype].extend(cr_tuples)
    #         #Removing duplicates    
    #         contact_residues[uniprot] = {ctype: list(set(cr_tuples)) for ctype, cr_tuples in contact_residues[uniprot].items()}
    #     return contact_residues

    # def sequence_differences(self) -> dict[str, list[tuple[int, str, str, str, str]]]:
    #     # check struct_ref_seq_dif for annotated differences
    #     chain_to_differences = defaultdict(list)
    #     dif = self.block.find('_struct_ref_seq_dif.', 
    #         ['pdbx_pdb_strand_id', 'seq_num', 'mon_id', 
    #         'db_mon_id', 'details'])
    #     for row in dif:
    #         chain_to_differences[row[0]].append((row[1], row[2], row[3], row[4]))
    #     return chain_to_differences
    @classmethod
    def extract_position_and_name(cls, residue: gemmi.Residue) -> tuple[str, str]:
        if cls.residue_representation['position'] == 'int_as_str':
            res_num = str(residue.seqid.num)
        elif cls.residue_representation['position'] == 'alphanumeric':
            res_num = str(residue.seqid)
        else:
            raise ValueError(f"Invalid residue representation: {cls.residue_representation['position']}")
        if cls.residue_representation['name'] == 'three_letter':
            res_name = residue.name
        elif cls.residue_representation['name'] == 'single_letter':
            res_name = seq1(residue.name)
        else:
            raise ValueError(f"Invalid residue representation: {cls.residue_representation['name']}")
        return res_num, res_name

    @staticmethod
    def get_saved_structures(structure_db: str) -> list[str]:
        if not Structure.SavedStructures[structure_db]:
            files = PDB_DIR.glob("*.cif")
            structure_ids = [file.stem.lower() for file in files]
            Structure.SavedStructures[structure_db].extend(structure_ids)
        return Structure.SavedStructures[structure_db]


    @staticmethod
    def select_unprocessed_structures(structure_ids: list[str], structure_db: str="pdb", cutoff: float=8.0) -> list[str]:
        structure_ids = [st.lower() for st in structure_ids]
        saved_structures = Structure.get_saved_structures(structure_db)
        eligible_structures = list(set(structure_ids) & set(saved_structures))
        processed_structures = Structure.processed_structures(structure_db, cutoff)
        return list(set(eligible_structures) - set(processed_structures))

    @staticmethod
    def save_binding_interfaces(contact_residues: dict[str, IntervalDict], structure_id: str, structure_db: str, cutoff: float):
        file_path = Structure.BI_path(structure_id, structure_db, cutoff)
        Path(file_path).parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "w") as f:
            json.dump(contact_residues, f)

    @staticmethod
    def load_binding_interfaces(structure_id: str, structure_db: str, cutoff: float):
        is_saved = Structure.is_saved_binding_interfaces(structure_id, structure_db, cutoff)
        if is_saved:
            file_path = Structure.BI_path(structure_id, structure_db, cutoff)
            with open(file_path, "r") as f:
                return json.load(f)
        return None

    @staticmethod
    def get_structure_path(structure_id:str, structure_db:str) -> Path:
        structure_id = structure_id.lower()
        if not structure_id.endswith('.cif'):
            structure_id = structure_id + '.cif'
        return PDB_DIR / structure_id

    @staticmethod
    def is_saved_structure(structure_id:str, structure_db:str) -> bool:
        file_path = Structure.get_structure_path(structure_id, structure_db)
        return file_path.exists()

    @staticmethod
    def is_saved_binding_interfaces(structure_id:str, structure_db:str, cutoff:float) -> bool:
        file_path = Structure.BI_path(structure_id, structure_db, cutoff)
        return file_path.exists()

    @staticmethod
    def processed_structures(structure_db: str, cutoff: float) -> list[str]:
        path = BI_STRUCTURE_DIR / "_".join([structure_db, str(cutoff)])
        return [file.name.split("_")[0] for file in path.glob("*.json")]
            
    @staticmethod
    def BI_path(structure_id: str, structure_db: str, cutoff: float) -> Path:
        db = structure_db
        id = structure_id
        path = BI_STRUCTURE_DIR / "_".join([db, str(cutoff)])/ ("_".join([id, db, str(cutoff)]) + ".json")
        return path
    

    # @staticmethod
    # def not_really_polypeptide(residues: list[gemmi.Residue]) -> bool:
    #     return False
    

    @classmethod
    def extract_binding_interfaces(cls, structure_id: str, 
                                        database: str = "pdb", 
                                        cutoff: float = 8.0, 
                                        try_saved_file: bool = True, 
                                        only_from_saved_file: bool = False) -> dict[str, dict[str, list[tuple[int, str]]]]|None:
        is_BI_saved = cls.is_saved_binding_interfaces(structure_id, database, cutoff)
        is_structure_saved = Structure.is_saved_structure(structure_id, database)
        if not is_structure_saved:
            print(f"Structure {structure_id} not available in local {database} structure database")
            return None
        if is_BI_saved and (try_saved_file or only_from_saved_file):
            contact_residues = cls.load_binding_interfaces(structure_id, database, cutoff)
        if not is_BI_saved and (not only_from_saved_file):
            try:
                structure = cls(structure_id, database=database)
                structure.initialize_binding_interfaces(cutoff=cutoff)
                contact_residues = structure.chainwise_contacts()
                cls.save_binding_interfaces(contact_residues, structure_id, database, cutoff)
            except Exception as e:
                print(f"Error extracting binding interfaces for {structure_id}")
                print(f"Error: {e}")
                return None
        return contact_residues
    

if __name__ == "__main__":
    db = "pdb"
    cutoff = 8.0
    log_interval = 50
    #uniprot_to_structure_ids = UniprotStructureMappings.representative_ids(structure_db="pdb", save=False, from_saved_file=True)
    #structure_ids = [st for sts in list(uniprot_to_structure_ids.values()) for st in sts]
    structures_ids = UniprotStructureMappings.get_mapped_structures(db)
    structures_to_process = Structure.select_unprocessed_structures(structures_ids, db, cutoff)
    print(f"Generating contacts for {len(structures_to_process)} structures")
    #random.shuffle(structure_ids)
    count_processed = 0
    while len(structures_to_process) > 0:
        structure_id = structures_to_process.pop(random.randrange(len(structures_to_process)))
        is_BI_saved = Structure.is_saved_binding_interfaces(structure_id, db, cutoff)
        if is_BI_saved:
            file_path = Structure.BI_path(structure_id, db, cutoff)
            print(f"Contacts for {structure_id} already saved at {file_path}")
            continue
        contact_residues = Structure.extract_binding_interfaces(structure_id, database=db, cutoff=cutoff)
        if contact_residues is None:
            print(f"Failed to extract contacts for {structure_id}")
            continue
        count_processed += 1
        if count_processed % log_interval == 0:
            structures_to_process = Structure.select_unprocessed_structures(structures_to_process, db, cutoff)
            #print(f"Structures with saved contacts: {Structure.processed_structures(db, cutoff)}")
            print(f"Generated contacts for {count_processed} structures")
            print(f"Remaining structures to process: {len(structures_to_process)}")
            #print(f"Unaligned chains: {count_chains_unaligned} out of {count_chains_processed} chains processed")
            
        
    
    # def map_uniprot_to_chains(self) -> dict[str, str]:
    #     block = self.cif_block
    #     uniprot_to_chains = {}
    #     for row in block.find('_struct_asym.', ['id', 'entity_id']):
    #         chain_id = row[0]
    #         entity_id = row[1]
    #         uniprot_id = self.entity_to_accession.get(entity_id, None)
    #         if uniprot_id is not None and uniprot_id.startswith('UNP'):
    #             uniprot_to_chains[uniprot_id.split(':')[1]].add(chain_id)
    #     return uniprot_to_chains


    # def map_uniprot_to_chains(self) -> dict[str, list[str]]:
    #     uniprot_to_chains = defaultdict(list)
    #     for ent in self.structure.entities:
    #         uniprot_id = self.entity_to_accession.get(ent.name, None)
    #         if uniprot_id is not None and uniprot_id.startswith('UNP'):
    #             uniprot_to_chains[uniprot_id.split(':')[1]].extend(list(ent.subchains))
    #     return uniprot_to_chains

   # def build_chains_to_uniprot_sequence_alignment(self) -> dict[str, dict[gemmi.SeqId, int]]:
    #     def _ins(cell: str) -> str:
    #         if gemmi.cif.is_null(cell):
    #             return ""
    #         else:
    #             s = gemmi.cif.as_string(cell).strip()
    #             return "" if s in {".", "?"} else s
        
    #     def _seqid(num_raw: str, ins_raw: str) -> gemmi.SeqId:
    #         n = int(gemmi.cif.as_string(num_raw).strip())
    #         ins = _ins(ins_raw)
    #         return gemmi.SeqId(f"{n}{ins}" if ins else str(n))

    #     tab = self.block.find(
    #         "_pdbx_poly_seq_scheme.",
    #         ["pdb_strand_id", "auth_seq_num", "pdb_ins_code",
    #         "pdbx_seq_db_accession_code", "pdbx_seq_db_seq_num"],
    #     )
    #     sequence_alignment: dict[str, dict[gemmi.SeqId, int]] = {}
    #     if not tab:
    #         return sequence_alignment
    #     for row in tab:
    #         if gemmi.cif.is_null(row[3]) or gemmi.cif.is_null(row[4]):
    #             continue
    #         chain = gemmi.cif.as_string(row[0]).strip()
    #         sid = _seqid(row[1], row[2])
    #         acc = gemmi.cif.as_string(row[3]).strip()
    #         resnum = int(gemmi.cif.as_string(row[4]).strip())
    #         sequence_alignment[acc+':'+chain][sid] = resnum
    #     return sequence_alignment

     # def uniprotChainKeys(self, uniprot) -> list[str]:
    #     unpChainKeys = [self.uniprotChainKey(uniprot, chain) for entities in self.uniprot_to_entity[uniprot]
    #                                   for entity in entities 
    #                                   for chain in self.entity_to_chains[entity]]
    #     return unpChainKeys

    # def uniprotChainKey(self, uniprot: str, chain: str) -> str:
    #     return uniprot+':'+chain


    #self.uniprot_to_refID = self.extract_uniprot_to_internalIDs_mappings(id_type='id')
    #self.entity_to_accession = self.map_entity_to_accession()
    #self.uniprot_to_chains = self.map_uniprot_to_chains()
    #self.uniprotChainKeys = {uniprot: self.uniprotChainKeys(uniprot) for uniprot in self.uniprotIds}
    #self.BI_extractor = BindingInterfacesExtractor(self.structure, self.block, cutoff=cutoff)
    #self.entitywise_contacts = self.entitywise_contacts()
