from dsa.binding_interfaces.process_structure.IntervalDict import IntervalDict
from Bio import Align
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
import pandas as pd
from dsa.binding_interfaces.config import PDB_DIR, PDB_UNIPROT_MAP_DIR, SIFT_PDB_to_uniprot_file  
import gemmi
from dsa.binding_interfaces.process_structure.Structure import Structure
from Bio.SeqUtils import seq1
import requests
import gzip
import xml.etree.ElementTree as ET
from random import shuffle
from collections import defaultdict
import numpy as np
import warnings
import re



class PDB_Uniprot_Mapper:
    residue_representation = {'postion': 'int_as_str', 'name': 'three_letter'}
    SIFT_Mapping_Table = pd.read_csv(SIFT_PDB_to_uniprot_file, sep="\t", comment="#")
    def __init__(self, structure_id: str):
        self.structure_id = structure_id.lower()
        self.structure: Structure = Structure(self.structure_id, database="pdb")
        self.sift_mapping = self.SIFT_Mapping_Table[self.SIFT_Mapping_Table["PDB"] == self.structure_id]
        #self.uniprot_ids = list(set(self.sift_mapping["SP_PRIMARY"].tolist())) 
        self.pdbe_mappings = self.load_pdbe_mappings(self.structure_id) 
        self.uniprot_sequences = self.initialize_uniprot_sequences()
        #self.pdb_sequences = defaultdict(lambda: defaultdict(lambda:dict))
        self.pdb_sequences = self.structure.chain_sequence_grouped_by_uniprot()
        self.mappings = defaultdict(lambda: defaultdict(dict))
        self.pdbe_found = True if self.pdbe_mappings is not None else False
        self.sift_found = True if not  self.sift_mapping.empty else False
        self.sift_used = False
        if self.pdbe_found:
            self.align_using_pdbe_mappings()
        elif self.sift_found:
            self.align_using_sift_mappings()
            self.sift_used = True
        self.mismatches = self.compute_mismatches()

    def initialize_uniprot_sequences(self):
        uniprot_sequences = defaultdict(str)
        uniprot_ids = self.structure.uniprotIds
        for uniprot_id in uniprot_ids: 
            seq = UniprotToSequence.get_sequence(uniprot_id)
            if seq is not None:
                uniprot_sequences[uniprot_id] = seq
        return uniprot_sequences

    def compute_mismatches(self):
        mismatches = defaultdict(lambda:defaultdict(int))
        for chain in self.mappings:
            for uniprot_id in self.mappings[chain]:
                mismatches[chain][uniprot_id] = 0
                for pdb_res, unp_res in self.mappings[chain][uniprot_id].items():
                    if pdb_res[1]!= unp_res[1]:
                        mismatches[chain][uniprot_id] += 1
        return mismatches

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


    @property 
    def accession_numbers(self) -> list[str]:
        ids = [id for dicts in self.mappings.values() for dct in dicts.values() for id in dct.keys()]
        return list(set(ids))
    
    def map(self, chain: str, residue: tuple[int, str]) -> tuple[int, str]:
        if chain not in self.mappings:
            return None
        if residue not in self.mappings[chain]:
            return None
        return self.mappings[chain][residue]

    def align_using_sift_mappings(self):
        for row in self.sift_mapping.itertuples(index=False):
            if row.PDB_BEG - row.PDB_END != row.SP_BEG - row.SP_END:
                continue
            chain = row.CHAIN
            uniprot_id = row.SP_PRIMARY
            unp_seq = self.uniprot_sequences[uniprot_id]
            if unp_seq is None:
                continue
            for pdb_resnum in range(row.PDB_BEG, row.PDB_END + 1):
                pdb_resname = self.pdb_sequences[uniprot_id][chain].get(str(pdb_resnum), None)
                if pdb_resname is None:
                    continue
                unp_resnum = row.SP_BEG + (pdb_resnum - row.PDB_BEG)
                unp_resname = unp_seq[unp_resnum-1]
                if seq1(pdb_resname) != unp_resname:
                    continue
                self.mappings[chain][uniprot_id][(str(pdb_resnum), pdb_resname)] = (unp_resnum, unp_resname)

    def is_consistent_with_pdb_chain(self, chain: str, position_residue_pairs: list[tuple[int, str]]):
        for pdb_resnum, pdb_resname in position_residue_pairs:
            if pdb_resname != self.pdb_sequences[chain][pdb_resnum]:
                return False
        return True

    def is_conistent_with_saved_uniprot_sequence(self, uniprot_id: str, position_residue_pairs: list[tuple[int, str]], offset: int = None):
        # Checks if the uniprot position sequence pairs in pdbe mappings are consistent with the saved UniProt sequence
        unp_seq = self.uniprot_sequences[uniprot_id]
        unp_seq_length = len(unp_seq)
        offsets = [0, 1, -1, 2, -2, 3, -3]
        if offset is not None:
            offsets = [offset]
        best_offset = offsets[0]
        best_mismatch_count = len(position_residue_pairs)
        is_consistent = False
        for offset in offsets:
            mismatch_count = 0
            for p, r in position_residue_pairs:
                p = int(p) if isinstance(p, str) else p
                if 0 <= p-1+offset < unp_seq_length and r != unp_seq[p-1+offset]:
                    mismatch_count += 1
            if mismatch_count < best_mismatch_count:
                best_mismatch_count = mismatch_count
                fraction_mismatches = mismatch_count/len(position_residue_pairs)
                best_offset = offset
                is_consistent = mismatch_count <= 5
                if is_consistent:
                    break
        if not is_consistent and len(offsets) > 1:
            is_consistent, best_offset, mismatch_count, fraction_mismatches = self.is_conistent_with_saved_uniprot_sequence(uniprot_id, position_residue_pairs, offset=best_offset)
            print(f"Inconsistent with saved uniprot sequence for {uniprot_id} in {self.structure_id}")
            print(f"Mismatch count: {mismatch_count}")
            print(f"Fraction mismatches: {fraction_mismatches}")
            print(f"Best offset: {best_offset}")
        if is_consistent and best_offset != 0:
            warnings.warn(f"Offset for {uniprot_id} is {best_offset}, not 0")
            print(f"Uniprot ID: {uniprot_id}")
            print(f"Structure ID: {self.structure_id}")
        return is_consistent, best_offset, mismatch_count, fraction_mismatches

    def number_of_mismatches(self):
        mismatches = defaultdict(lambda:defaultdict(int))
        for chain in self.mappings:
            for uniprot_id in self.mappings[chain]:
                mismatches[chain][uniprot_id] = 0
                for residue in self.mappings[chain][uniprot_id]:
                    unp_resname = self.mappings[chain][uniprot_id][residue][1]
                    pdb_resname = residue[1]
                    if unp_resname!= seq1(pdb_resname):
                        mismatches[chain][uniprot_id] += 1
        return mismatches

    def align_using_pdbe_mappings(self):
        chains = list(set(self.pdbe_mappings["chain_id"].tolist()))
        int2str = lambda x: x if isinstance(x, str) else str(x)
        str2int = lambda x: x if isinstance(x, int) else int(x)
        for chain in chains:
            chain_mapping = self.pdbe_mappings[self.pdbe_mappings["chain_id"] == chain]
            for uniprot_id in list(set(chain_mapping["uniprot_id"].tolist())):
                chain_uniprot_mapping = chain_mapping[chain_mapping["uniprot_id"] == uniprot_id]
                unp_pairs = chain_uniprot_mapping[["unp_resnum", "unp_resname"]].values.tolist()
                is_consistent, offset, mismatch_count, fraction_mismatches = self.is_conistent_with_saved_uniprot_sequence(uniprot_id, unp_pairs)
                if not is_consistent:
                    warnings.warn(f"Inconsistent with saved uniprot sequence for {uniprot_id} on chain {chain}")
                    continue
                if offset != 0:
                    warnings.warn(f"Offset for {uniprot_id} on chain {chain} is {offset}, not 0")
                    print(f"Offset for {uniprot_id} on chain {chain} is {offset}")
                    print(f"Mismatch count: {mismatch_count}")
                    print(f"Fraction mismatches: {fraction_mismatches}")
                    print(f"Best offset: {offset}")
                for row in chain_uniprot_mapping.itertuples(index=False):
                    pdb_resnum = extract_num_str(int2str(row.pdb_resnum))
                    pdb_resname = seq1(row.pdb_resname)
                    unp_resnum = str2int(row.unp_resnum)
                    unp_resname = row.unp_resname
                    if pdb_resname == 'X' or pdb_resname == '':
                        warnings.warn(f"PDB residue name {row.pdb_resname} converted to invalid residue name {pdb_resname} at position {pdb_resnum} for {uniprot_id} on chain {chain}")
                        continue
                    if (pdb_resnum, pdb_resname) in self.mappings[chain][uniprot_id]:
                        mapped_unp_resname = self.mappings[chain][uniprot_id][(pdb_resnum, pdb_resname)][1]
                        if mapped_unp_resname == pdb_resname:
                                continue
                    self.mappings[chain][uniprot_id][(pdb_resnum, pdb_resname)] = (unp_resnum + offset, unp_resname)

    @staticmethod
    def download_mappings(structure_id: str):
        url = f"https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{structure_id.lower()}.xml.gz"
        response = requests.get(url)
        # decompress and parse XML
        try:
            xml_content = gzip.decompress(response.content)
        except Exception as e:
            print(f"Error decompressing {structure_id}: {e}")
            return None
        root = ET.fromstring(xml_content)
        # namespace used in SIFTS XML
        ns = {'sifts': 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'}
        residue_mappings = []
        for entity in root.findall('.//sifts:entity', ns):
            chain_id = entity.get('entityId')  # this is the chain ID    
            for residue in entity.findall('.//sifts:residue', ns):
                pdb_res = residue.find('sifts:crossRefDb[@dbSource="PDB"]', ns)
                unp_res = residue.find('sifts:crossRefDb[@dbSource="UniProt"]', ns)    
                if pdb_res is not None and unp_res is not None:
                    residue_mappings.append({
                        'chain_id':    chain_id,
                        'pdb_resnum':  pdb_res.get('dbResNum'),
                        'pdb_resname': pdb_res.get('dbResName'),
                        'uniprot_id':  unp_res.get('dbAccessionId'),
                        'unp_resnum':  unp_res.get('dbResNum'),
                        'unp_resname': unp_res.get('dbResName'),
                        #'db_name':     unp_res.get('dbSource'), # add later since it is not always UniProt
                    })
        if residue_mappings is not None:
            pd.DataFrame(residue_mappings).to_csv(PDB_UNIPROT_MAP_DIR / f"{structure_id}.tsv", sep="\t", index=False)
        return residue_mappings

    

    @staticmethod
    def missing_mappings(structure_ids: str|list[str]|None=None):
        if structure_ids  is None:
            structure_ids = [pdb_id.name.split(".")[0] for pdb_id in PDB_DIR.glob("*.cif")]
        else:
            if isinstance(structure_ids, str):
                structure_ids = [structure_ids]
            else:
                structure_ids = structure_ids
        mapped_structure_ids = [pdb_id.name.split(".")[0] for pdb_id in PDB_UNIPROT_MAP_DIR.glob("*.tsv")]
        unmapped_structure_ids = set(structure_ids) - set(mapped_structure_ids)
        return list(unmapped_structure_ids)

    @staticmethod
    def download_missing_mappings(structure_ids: str|list[str]|None=None, random_order: bool = True):
        unmapped_structure_ids = PDB_Uniprot_Mapper.missing_mappings(structure_ids)
        if random_order:
            shuffle(unmapped_structure_ids)
        count = 0
        for structure_id in unmapped_structure_ids:
            mappings = PDB_Uniprot_Mapper.download_mappings(structure_id)
            if mappings is None:
                continue
            count += 1
            if count % 100 == 0:
                PDB_Uniprot_Mapper.download_missing_mappings()
                break

    @staticmethod
    def load_pdbe_mappings(structure_id: str):
        success = False
        tried_downloading = False
        while not success and not tried_downloading:
            try:
                mappings =  pd.read_csv(PDB_UNIPROT_MAP_DIR / f"{structure_id}.tsv", sep="\t")
                success = True
            except FileNotFoundError:
                try:
                    mappings = PDB_Uniprot_Mapper.download_mappings(structure_id)
                    tried_downloading = True
                except Exception as e:
                    mappings = None
                    tried_downloading = True
        return mappings


def extract_num_str(s) -> str | None:
    m = re.search(r'-?\d+', str(s))
    return m.group() if m else None

if __name__ == "__main__":
    # PDB_Uniprot_Mapper.download_missing_mappings(random_order=False)
    # print(f"Number of missing mappings: {len(PDB_Uniprot_Mapper.missing_mappings())}")
    mapper = PDB_Uniprot_Mapper("9bq2")
    print(mapper.mappings)
    print(mapper.number_of_mismatches())




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