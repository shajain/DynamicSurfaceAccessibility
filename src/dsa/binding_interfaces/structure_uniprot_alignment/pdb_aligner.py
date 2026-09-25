#from curses import A_BOLD
#from dsa.binding_interfaces.external_validation_delete_later.withPDBePISA import PDBE_API
from dsa.binding_interfaces.structure_uniprot_alignment.aligner import StructureSequenceAligner
from dsa.binding_interfaces.structure_uniprot_alignment.interval_dict import IntervalDict
from Bio import Align
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
import pandas as pd
#from dsa.binding_interfaces.config import PDB_DIR, PDB_UNIPROT_MAP_DIR, SIFT_PDB_to_uniprot_file  
from dsa.binding_interfaces.config import PDB_DIR  
import gemmi
from dsa.binding_interfaces.structure.structure_model import StructureModel
#from dsa.binding_interfaces.structure.structure import Structure
from Bio.SeqUtils import seq1
import requests
import gzip
import xml.etree.ElementTree as ET
from random import shuffle
from collections import defaultdict
import numpy as np
import warnings
import re
from typing import Literal
from dsa.binding_interfaces.structure_uniprot_alignment.utilities import extract_num_str, int_if_str
from dsa.binding_interfaces.common_functions import single_letter_residue_code, is_actually_a_number
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure_uniprot_alignment.alignment_store import PDBeUniprotAlignmentStore, SIFTMappingsStore
import numbers


@ClassFactory.register(key="pdb", interface=StructureSequenceAligner)
class PDBSequenceAligner(StructureSequenceAligner):
    pdb_residue_representation = {'postion': 'int_as_str', 'name': 'three_letter'}
    sift_alignment_store = SIFTMappingsStore()
    pdbe_alignment_store = PDBeUniprotAlignmentStore()
    def __init__(self, structure_id: str, structure_model: StructureModel, uniprot_to_entity: dict[str, list[str]]):
    #def __init__(self, structure: Structure):
        #self.structure_id = structure_id
        self.sift_mappings = self.sift_alignment_store.load_mappings(structure_id)
        #self.pdbe_mappings = self.load_pdbe_mappings(structure_id) 
        self.pdbe_mappings = self.pdbe_alignment_store.load_mappings(structure_id)
        # self.uniprot_sequences = self.initialize_uniprot_sequences()
        # self.chain_sequences = {chain: self.structure_model.get_chain_sequence_as_dict(chain) for chain in self.structure_model.chains}
        #self.mappings = defaultdict(lambda: defaultdict(dict))
        self.pdbe_found = True if self.pdbe_mappings is not None else False
        self.sift_found = True if not  self.sift_mappings.empty else False
        
        super().__init__(structure_id, structure_model, uniprot_to_entity)
        # if self.pdbe_found:
        #     self.align_using_pdbe_mappings()
        # elif self.sift_found:
        #     self.align_using_sift_mappings()
        #     self.sift_used = True
        
    @property
    def sift_used(self) -> bool:
        return self.sift_found and not self.pdbe_found
    
    @property
    def pdbe_used(self) -> bool:
        return self.pdbe_found




    # def compute_mismatch_metric(self, chain: str, uniprot_id: str, type: Literal['fraction', 'count'] = 'fraction'):
    #     mapping = self.mappings[chain][uniprot_id]
    #     mismatch_count = 0
    #     mismatch_metric = 0
    #     for pdb_res, unp_res in mapping.items():
    #         if single_letter_residue_code(pdb_res[1])!= unp_res[1]:
    #             mismatch_count += 1
    #         if type == 'fraction':
    #             mismatch_metric = mismatch_count / len(mapping)
    #         elif type == 'count':
    #                 mismatch_metric = mismatch_count
    #         else:   
    #             raise ValueError(f"Invalid type in compute_mismatch_metric: {type}")
    #     return mismatch_metric

    # def initialize_pdb_sequences(self):
    #     if self.pdb_sequences:
    #         return
    #     if self.structure is None:
    #         self.structure = Structure(self.structure_id)
    #     for chain in self.structure[0]:
    #         for res in chain:
    #             if res.entity_type != gemmi.EntityType.Polymer:
    #                 continue
    #             self.pdb_sequences[self.structure.entity_to_uniprot[res.entity_id]][chain.name][str(res.seqid)]= res.name


    # @property 
    # def accession_numbers(self) -> list[str]:
    #     ids = [id for dicts in self.mappings.values() for dct in dicts.values() for id in dct.keys()]
    #     return list(set(ids))
    

    def _generate_mapping(self, chain: str, uniprot_id: str):
        if self.pdbe_used:
            return self._generate_mapping_pdbe(chain, uniprot_id)
        elif self.sift_used:
            return self._generate_mapping_sift(chain, uniprot_id)
        else:
            warnings.warn("No mapping resource available")
        return {}



    def _generate_mapping_sift(self, chain: str, uniprot_id: str):
        mapping = {}
        mapping_df = self.sift_mappings[self.sift_mappings["SP_PRIMARY"] == uniprot_id]
        mapping_df = mapping_df[mapping_df["CHAIN"] == chain]
        structure_residues = self.chain_sequences[chain]
        if mapping_df.empty:
            return {}
        for row in mapping_df.itertuples(index=False):
            pdb_beg = self.extract_position(row.PDB_BEG)
            pdb_end = self.extract_position(row.PDB_END)
            unp_beg = self.extract_position(row.SP_BEG)
            unp_end = self.extract_position(row.SP_END)
            if pdb_beg is None or pdb_end is None or unp_beg is None or unp_end is None:
                continue
            if pdb_beg - pdb_end != unp_beg - unp_end:
            #if False:
                continue
            for pdb_resnum in range(pdb_beg, pdb_end + 1):
                pdb_resname = structure_residues.get(str(pdb_resnum), None)
                pdb_residue_key = self.structure_residue_key((pdb_resnum, pdb_resname))
                # if pdb_residue_key[1] is None or pdb_residue_key[1] == '' or pdb_residue_key[1] == 'X':
                #     continue
                #unp_resnum = unp_beg + (pdb_resnum - pdb_beg)
                unp_resnum = self.offset_adjusted_uniprot_position((pdb_resnum - pdb_beg), offset=unp_beg)
                if unp_resnum is None:
                    continue
                unp_resname = self.uniprot_sequences[uniprot_id][unp_resnum-1]
                mapping[pdb_residue_key] = (unp_resnum, unp_resname)
        return mapping

    # def align_using_sift_mappings(self):
    #     for row in self.sift_mapping.itertuples(index=False):
    #         pdb_beg = int_if_str(row.PDB_BEG)
    #         pdb_end = int_if_str(row.PDB_END)
    #         unp_beg = int_if_str(row.SP_BEG)
    #         unp_end = int_if_str(row.SP_END)
    #         if pdb_beg is None or pdb_end is None or unp_beg is None or unp_end is None:
    #             continue
    #         if pdb_beg - pdb_end != unp_beg - unp_end:
    #             continue
    #         chain = row.CHAIN
    #         uniprot_id = row.SP_PRIMARY
    #         unp_seq = self.uniprot_sequences[uniprot_id]
    #         if unp_seq is None:
    #             continue
    #         for pdb_resnum in range(pdb_beg, pdb_end + 1):
    #             pdb_resname = self.chain_sequences[chain].get(str(pdb_resnum), None)
    #             pdb_resname = single_letter_residue_code(pdb_resname)
    #             if pdb_resname is None or pdb_resname == '' or pdb_resname == 'X':
    #                 continue
    #             unp_resnum = unp_beg + (pdb_resnum - pdb_beg)
    #             unp_resname = unp_seq[unp_resnum-1]
    #             # if single_letter_resname(pdb_resname) != unp_resname:
    #             #     continue
    #             self.mappings[chain][uniprot_id][(str(pdb_resnum), pdb_resname)] = (unp_resnum, unp_resname)


    def is_pdbe_consistent_with_saved_uniprot_sequence(self, uniprot_id: str, pdbe_uniprot_sequence: list[tuple[int, str]], offset: int = None):
        unp_seq = self.uniprot_sequences[uniprot_id]
        best_offset, mismatch_count = self.best_alignment_offset(unp_seq, pdbe_uniprot_sequence, offset) 
        is_consistent = self.is_acceptable_mismatch(mismatch_count, len(pdbe_uniprot_sequence))
        if not is_consistent:
            return False, None
        return True, best_offset

    # def is_consistent_with_pdb_chain(self, chain: str, position_residue_pairs: list[tuple[int, str]]):
    #     for pdb_resnum, pdb_resname in position_residue_pairs:
    #         if pdb_resname != self.pdb_sequences[chain][pdb_resnum]:
    #             return False
    #     return True

    # def is_pdbe_conistent_with_saved_uniprot_sequence(self, uniprot_id: str, position_residue_pairs: list[tuple[int, str]], offset: int = None):
    #     # Checks if the uniprot position sequence pairs in pdbe mappings are consistent with the saved UniProt sequence
    #     unp_seq = self.uniprot_sequences[uniprot_id]
    #     unp_seq_length = len(unp_seq)
    #     offsets = [0, 1, -1, 2, -2, 3, -3]
    #     if offset is not None:
    #         offsets = [offset]
    #     best_offset = offsets[0]
    #     best_mismatch_count = len(position_residue_pairs)
    #     is_consistent = False
    #     for offset in offsets:
    #         mismatch_count = 0
    #         for p, r in position_residue_pairs:
    #             p = int(p) if isinstance(p, str) else p
    #             if 0 <= p-1+offset < unp_seq_length and r != unp_seq[p-1+offset]:
    #                 mismatch_count += 1
    #         if mismatch_count < best_mismatch_count:
    #             best_mismatch_count = mismatch_count
    #             fraction_mismatches = mismatch_count/len(position_residue_pairs)
    #             best_offset = offset
    #             is_consistent = mismatch_count <= 5
    #             if is_consistent:
    #                 break
    #     if not is_consistent and len(offsets) > 1:
    #         is_consistent, best_offset, mismatch_count, fraction_mismatches = self.is_pdbe_conistent_with_saved_uniprot_sequence(uniprot_id, position_residue_pairs, offset=best_offset)
    #         print(f"Inconsistent with saved uniprot sequence for {uniprot_id} in {self.structure_id}")
    #         print(f"Mismatch count: {mismatch_count}")
    #         print(f"Fraction mismatches: {fraction_mismatches}")
    #         print(f"Best offset: {best_offset}")
    #     if is_consistent and best_offset != 0:
    #         warnings.warn(f"Offset for {uniprot_id} is {best_offset}, not 0")
    #         print(f"Uniprot ID: {uniprot_id}")
    #         print(f"Structure ID: {self.structure_id}")
    #     return is_consistent, best_offset, mismatch_count, fraction_mismatches


    # def align_using_pdbe_mappings(self):
    #     chains = list(set(self.pdbe_mappings["chain_id"].tolist()))
    #     # int2str = lambda x: x if isinstance(x, str) else str(x)
    #     # str2int = lambda x: x if isinstance(x, int) else int(x)
    #     for chain in chains:
    #         chain_mapping = self.pdbe_mappings[self.pdbe_mappings["chain_id"] == chain]
    #         for uniprot_id in list(set(chain_mapping["uniprot_id"].tolist())):
    #             chain_uniprot_mapping = chain_mapping[chain_mapping["uniprot_id"] == uniprot_id]
    #             unp_pairs = chain_uniprot_mapping[["unp_resnum", "unp_resname"]].values.tolist()
    #             is_consistent, offset, mismatch_count, fraction_mismatches = self.is_conistent_with_saved_uniprot_sequence(uniprot_id, unp_pairs)
    #             if not is_consistent:
    #                 warnings.warn(f"Inconsistent with saved uniprot sequence for {uniprot_id} on chain {chain}")
    #                 continue
    #             # if offset != 0:
    #             #     warnings.warn(f"Offset for {uniprot_id} on chain {chain} is {offset}, not 0")
    #             #     print(f"Offset for {uniprot_id} on chain {chain} is {offset}")
    #             #     print(f"Mismatch count: {mismatch_count}")
    #             #     print(f"Fraction mismatches: {fraction_mismatches}")
    #             #     print(f"Best offset: {offset}")
    #             for row in chain_uniprot_mapping.itertuples(index=False):
    #                 #pdb_resnum = extract_num_str(int2str(row.pdb_resnum))
    #                 # if pd.isna(row.pdb_resnum):
    #                 #     warnings.warn(f"Invalid PDB residue number {row.pdb_resnum} for {uniprot_id} on chain {chain}")
    #                 #     continue
    #                 # pdb_resnum = int2str(row.pdb_resnum)
    #                 pdb_resnum = row.pdb_resnum
    #                 pdb_resname = row.pdb_resname
    #                 pdb_residue = (pdb_resnum, pdb_resname)
    #                 pdb_residue_key = self.structure_residue_key(pdb_residue)
    #                 pdb_resname_sl = single_letter_residue_code(pdb_residue_key[1]).upper()
    #                 unp_resnum = str2int(row.unp_resnum)
    #                 unp_resname = single_letter_residue_code(row.unp_resname).upper()
    #                 if pdb_resname_sl == 'X' or pdb_resname_sl == '' or pdb_residue_key[0] is None:
    #                     warnings.warn(f"invalid PDB residue name {pdb_residue_key[1]} at position {pdb_residue_key[0]} for {uniprot_id} on chain {chain}")
    #                     continue
    #                 if pdb_residue_key in self.mappings[chain][uniprot_id]:
    #                     mapped_unp_resname = self.mappings[chain][uniprot_id][pdb_residue_key][1]
    #                     if mapped_unp_resname == pdb_resname_sl:
    #                             continue
    #                 self.mappings[chain][uniprot_id][pdb_residue_key] = (unp_resnum + offset, unp_resname)

    def _generate_mapping_pdbe(self, chain: str, uniprot_id: str):
        mapping = defaultdict(tuple[int, str])
        mapping_pdbe = self.pdbe_mappings[self.pdbe_mappings["pdb_chain_id"] == chain]
        mapping_pdbe = mapping_pdbe[mapping_pdbe["uniprot_id"] == uniprot_id]
        if mapping_pdbe.empty:
            return {}
        uniprot_sequence_in_pdbe = mapping_pdbe[["unp_resnum", "unp_resname"]].values.tolist()
        is_consistent, offset = self.is_pdbe_consistent_with_saved_uniprot_sequence(uniprot_id, uniprot_sequence_in_pdbe)
        if not is_consistent:
            warnings.warn(f"PDBe mapping is inconsistent with saved uniprot sequence for {uniprot_id} on chain {chain}")
            return {}
        #offset = []
        #mapping_pdbe = self.resolve_missing_seqids(mapping_pdbe, self.chain_sequences[chain])
        for row in mapping_pdbe.itertuples(index=False):
            pdb_residue_key = self.structure_residue_key((row.pdb_resnum, row.pdb_resname))
            unp_residue = self.uniprot_residue_as_value(row.unp_resnum, row.unp_resname, offset)
            if pd.isna(row.unp_resnum) or pdb_residue_key[0] is None or unp_residue[0] is None:
                continue
            if pdb_residue_key in mapping and not self.is_mismatch(pdb_residue_key, mapping[pdb_residue_key]):
                continue
            mapping[pdb_residue_key] = unp_residue
        return mapping



    def uniprot_residue_as_value(self, unp_resnum: int, unp_resname: str, offset: int = 0) -> tuple[int, str]:
        unp_resnum = self.extract_position(unp_resnum)
        unp_resname = self.single_letter_residue_code(unp_resname)
        unp_resnum = self.offset_adjusted_uniprot_position(unp_resnum, offset)
        return (unp_resnum, unp_resname)

    # @classmethod
    # def resolve_missing_seqids(cls, mapping_pdbe: pd.DataFrame, chain_sequences: dict[str, str]) -> pd.DataFrame:
    #     chain_seqnum_to_seqid = defaultdict(list)
    #     for seqid, resname in chain_sequences.items():
    #         p = cls.extract_position(seqid)
    #         if p is not None:
    #             chain_seqnum_to_seqid[p].append(seqid)
    #     seqnums = list(chain_seqnum_to_seqid.keys())
    #     seqnums.sort()
    #     missing_seqnum_streak = 0
    #     previous_seqid_matched = None
    #     for row in mapping_pdbe.itertuples(index=False):
    #         pdbe_residue_as_key = cls.structure_residue_key((row.pdb_resnum, row.pdb_resname))
    #         if pdbe_residue_as_key[0] is not None:
    #             continue
    #         pdbe_resname_single_letter = cls.single_letter_residue_code(pdbe_residue_as_key)
    #         if previous_seqid_matched is not None:
    #             previous_seqnum_matched = cls.extract_position(previous_seqid_matched)
    #             candidate_seqids = chain_seqnum_to_seqid.get(previous_seqnum_matched, []) + chain_seqnum_to_seqid.get(previous_seqnum_matched+1, [])
    #             candidate_seqids.remove(previous_seqid_matched)
    #             for seqid in candidate_seqids:
    #                 if pdbe_resname_single_letter == cls.single_letter_residue_code(chain_sequences[seqid]):
    #             and cls.single_letter_residue_code(residue_as_key[1]) != cls.single_letter_residue_code(chain_sequences[previous_seqnum_matched]):

    #         if row.pdb_resnum 
    #         pdb_resnum = row.pdb_resnum
    #         pdb_resname = row.pdb_resname
    #         pdb_residue_key = cls.structure_residue_key((pdb_resnum, pdb_resname))
    #         if pdb_residue_key[0] is None:
    #             continue
    #         chain_sequence = chain_sequences[row.chain_id]
    #         if pdb_resnum not in chain_sequence:

    

    # @staticmethod
    # def download_mappings(structure_id: str):
    #     url = f"https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{structure_id.lower()}.xml.gz"
    #     response = requests.get(url)
    #     # decompress and parse XML
    #     try:
    #         xml_content = gzip.decompress(response.content)
    #     except Exception as e:
    #         print(f"Error decompressing {structure_id}: {e}")
    #         return None
    #     root = ET.fromstring(xml_content)
    #     # namespace used in SIFTS XML
    #     ns = {'sifts': 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'}
    #     residue_mappings = []
    #     for entity in root.findall('.//sifts:entity', ns):
    #         chain_id = entity.get('entityId')  # this is the chain ID    
    #         for residue in entity.findall('.//sifts:residue', ns):
    #             pdb_res = residue.find('sifts:crossRefDb[@dbSource="PDB"]', ns)
    #             unp_res = residue.find('sifts:crossRefDb[@dbSource="UniProt"]', ns)    
    #             if pdb_res is not None and unp_res is not None:
    #                 residue_mappings.append({
    #                     'chain_id':    chain_id,
    #                     'pdb_resnum':  pdb_res.get('dbResNum'),
    #                     'pdb_resname': pdb_res.get('dbResName'),
    #                     'uniprot_id':  unp_res.get('dbAccessionId'),
    #                     'unp_resnum':  unp_res.get('dbResNum'),
    #                     'unp_resname': unp_res.get('dbResName'),
    #                     #'db_name':     unp_res.get('dbSource'), # add later since it is not always UniProt
    #                 })
    #     if residue_mappings is not None:
    #         pd.DataFrame(residue_mappings).to_csv(PDB_UNIPROT_MAP_DIR / f"{structure_id}.tsv", sep="\t", index=False)
    #     return residue_mappings

    # @staticmethod
    # def get_structure_residue_key(structure_residue: tuple[str, str]|gemmi.Residue) -> tuple[str, str]:
    #     if isinstance(structure_residue, gemmi.Residue):
    #         structure_residue = (str(structure_residue.seqid), structure_residue.name)
    #     return structure_residue

    # @staticmethod
    # def missing_mappings(structure_ids: str|list[str]|None=None):
    #     if structure_ids  is None:
    #         structure_ids = [pdb_id.name.split(".")[0] for pdb_id in PDB_DIR.glob("*.cif")]
    #     else:
    #         if isinstance(structure_ids, str):
    #             structure_ids = [structure_ids]
    #         else:
    #             structure_ids = structure_ids
    #     mapped_structure_ids = [pdb_id.name.split(".")[0] for pdb_id in PDB_UNIPROT_MAP_DIR.glob("*.tsv")]
    #     unmapped_structure_ids = set(structure_ids) - set(mapped_structure_ids)
    #     return list(unmapped_structure_ids)

    # @staticmethod
    # def download_missing_mappings(structure_ids: str|list[str]|None=None, random_order: bool = True):
    #     unmapped_structure_ids = PDBSequenceAligner.missing_mappings(structure_ids)
    #     if random_order:
    #         shuffle(unmapped_structure_ids)
    #     count = 0
    #     for structure_id in unmapped_structure_ids:
    #         mappings = PDBSequenceAligner.download_mappings(structure_id)
    #         if mappings is None:
    #             continue
    #         count += 1
    #         if count % 100 == 0:
    #             PDBSequenceAligner.download_missing_mappings()
    #             break

    # @staticmethod
    # def load_pdbe_mappings(structure_id: str):
    #     success = False
    #     tried_downloading = False
    #     while not success and not tried_downloading:
    #         try:
    #             dtype = {'pdb_resnum': str, 'unp_resnum': int}
    #             mappings =  pd.read_csv(PDB_UNIPROT_MAP_DIR / f"{structure_id}.tsv", sep="\t", dtype=dtype)
    #             success = True
    #         except FileNotFoundError:
    #             try:
    #                 mappings = PDBSequenceAligner.download_mappings(structure_id)
    #                 tried_downloading = True
    #             except Exception as e:
    #                 mappings = None
    #                 tried_downloading = True
    #     return mappings

    # @classmethod
    # def int_or_None(cls, x: int|str|None) -> int|None:
    #     try:
    #         return int(x)
    #     except (ValueError, TypeError):
    #         return None
    

   

# how to handle numpy ints or numeric   
# def int2str(x):
#     if x is None or pd.isna(x):                 # NaN, None, pd.NA  -> None
#         return None
#     if isinstance(x, str):         # already a string (e.g. "55", "55A") -> keep
#         return x
#     if isinstance(x, float):
#         return str(int(x))
#     return str(x)    # int / numpy int  -> "55"

# def str2int(x):
#     if x is None or pd.isna(x):                 
#         return None
#     if isinstance(x, int):         
#         return x
#     return int(x)



# if __name__ == "__main__":
#     # PDB_Uniprot_Mapper.download_missing_mappings(random_order=False)
#     # print(f"Number of missing mappings: {len(PDB_Uniprot_Mapper.missing_mappings())}")
#     mapper = PDBSequenceAligner("2v89")
#     print(mapper.mappings['C']['P68431'])
#     structure = Structure("2v89")
#     print(structure.structure[0][2].name)
#     print(len([r for r in structure.structure[0][2] if r.entity_type == gemmi.EntityType.Polymer]))
#     for c, u_dict in mapper.mappings.items():
#         for u in u_dict.keys():
#             print(c, u, mapper.compute_mismatch_metric(chain=c, uniprot_id=u, type="count"))
    
    





#@staticmethod
    # def build_IntervalDict(pdb_resnums: list[int], uniprot_resnums: list[int]) -> IntervalDict:
    #     interval_dict = IntervalDict()
    #     pdb_resnums = np.array(pdb_resnums)
    #     uniprot_resnums = np.array(uniprot_resnums)
    #     ix = np.argsort(pdb_resnums)
    #     pdb_resnums = pdb_resnums[ix]
    #     uniprot_resnums = uniprot_resnums[ix]
    #     pdb_diff = np.diff(pdb_resnums)
    #     uniprot_diff = np.diff(uniprot_resnums)
    #     pdb_interval_boundaries = np.where(pdb_diff >1)[0]
    #     uniprot_interval_boundaries = np.where(uniprot_diff >1)[0]
    #     pdb_repeats = np.where(pdb_diff == 0)[0]
    #     has_repeats = len(pdb_repeats) > 0
    #     left_boundary = 0
    #     right_boundary = None
    #     for ix in uniprot_interval_boundaries + pdb_interval_boundaries:
    #         if right_boundary is not None:
    #             right_boundary = pdb_resnums[ix]
    #             offset = uniprot_resnums[ix] - pdb_resnums[ix]
    #             interval_dict.add(left_boundary, right_boundary, offset)
    #             left_boundary = None
    #         if left_boundary is None:
    #         if ix in pdb_interval_boundaries and ix in uniprot_interval_boundaries:         
    #     return interval_dict