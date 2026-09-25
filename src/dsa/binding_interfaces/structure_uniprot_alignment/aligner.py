# from dsa.binding_interfaces.structure_uniprot_alignment.pdb_aligner import PDBSequenceAligner
# from dsa.binding_interfaces.structure_uniprot_alignment.alphafold_aligner import AlphaFoldSequenceAligner
#from dsa.binding_interfaces.structure.structure import Structure
# from dsa.binding_interfaces.structure_uniprot_alignment.pdb_aligner import PDBSequenceAligner
# from dsa.binding_interfaces.structure_uniprot_alignment.alphafold_aligner import AlphaFoldSequenceAligner
import gemmi
from dsa.binding_interfaces.structure.structure_model import StructureModel
#from dsa.binding_interfaces.structure.structure import Structure
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
from collections import defaultdict
import warnings
from abc import abstractmethod
from dsa.binding_interfaces.common_functions import single_letter_residue_code, is_actually_a_number



class StructureSequenceAligner:
    #def __init__(self, structure_id: str, structure_model: StructureModel, uniprot_to_entity: dict[str, list[str]]):
    def __init__(self, structure_id: str, structure_model: StructureModel, uniprot_to_entity: dict[str, list[str]]):
        self.structure_id = structure_id
        self.structure_model = structure_model
        self.uniprot_to_entity = uniprot_to_entity
        self.uniprot_sequences = self.initialize_uniprot_sequences()
        self.chain_sequences = {chain_id: self.structure_model.get_chain_sequence_as_dict(chain_id) for chain_id in self.structure_model.chain_ids}
        # self.mappings = defaultdict(lambda: defaultdict(dict))
        self.mappings = self.generate_mappings()
        self.alignment_success = self.determine_alignment_success()
        
            #self.aligner = self.__class__.factory(structure)
        # self.structure_db = structure_db
        # self.structure_id = structure_id
        # if self.structure_db == "pdb":
        #     self.aligner = PDBSequenceAligner(self.structure_id, structure_db=self.structure_db)
        # elif self.structure_db == "alphafold":
        #     self.aligner = AlphaFoldSequenceAligner(self.structure_id, structure_db=self.structure_db)
        # else:
        #     raise ValueError(f"Invalid structure database: {self.structure_db}")

    

    @property
    def uniprot_ids(self) -> list[str]:
        return list(self.uniprot_to_entity.keys())

    def _is_alignment_successful(self, chain: str, uniprot_id: str) -> bool:
        if uniprot_id != self.uniprot_id or chain != self.chain or len(self.structure_sequence) == 0:
            return False
        #structure_sequence = [(res.seqid.num, res.name) for res in self.structure.structure[0][chain] if res.entity_type == gemmi.EntityType.Polymer]
        count = 0
        structure_sequence_length = len(self.structure_sequence)
        for residue_in_structure in self.structure_sequence:
            if not self.is_match(residue_in_structure):
                warnings.warn(f"Mismatch found between UniProt and AlphaFold for chain {chain} and UniProt ID {uniprot_id}")
                count += 1
        fraction_mismatched = count / structure_sequence_length
        return fraction_mismatched <= 0.1

    
    def initialize_uniprot_sequences(self):
        uniprot_sequences = defaultdict(str)
        for uniprot_id in self.uniprot_ids: 
            seq = UniprotToSequence.get_sequence(uniprot_id)
            if seq is not None:
                uniprot_sequences[uniprot_id] = seq
        return uniprot_sequences


    def chains_containing_uniprot_id(self, uniprot_id) -> list[str]:
        # chain_to_uniprot = defaultdict(list)
        entities = self.uniprot_to_entity[uniprot_id]
        chains = [c for e in entities for c in self.structure_model.get_chains_containing_entity(e)]
        return list(set(chains))


    def generate_mappings(self):
        mappings = defaultdict(lambda: defaultdict(dict))
        for uniprot_id in self.uniprot_ids:
            for chain in self.chains_containing_uniprot_id(uniprot_id):
                mappings[chain][uniprot_id] = self._generate_mapping(chain, uniprot_id)
        return mappings

    def _generate_mapping(self, chain: str, uniprot_id: str) -> dict[tuple[str, str], tuple[int, str]]:
        uniprot_sequence = self.uniprot_sequences[uniprot_id]
        structure_residues = list(self.chain_sequences[chain].items())
        mapping = defaultdict(dict)
        best_offset, _ = self.best_alignment_offset(uniprot_sequence, structure_residues)
        for structure_residue in structure_residues:
            structure_residue = self.structure_residue_key(structure_residue)
            structure_residue_pos = structure_residue[0]
            structure_residue_pos = self.extract_position(structure_residue_pos)
            if not isinstance(structure_residue_pos, int):
                continue
            unp_pos = self.offset_adjusted_uniprot_position(structure_residue_pos, best_offset)
            if 1 <= unp_pos <= len(uniprot_sequence):
                unp_resname = uniprot_sequence[unp_pos-1]
                mapping[structure_residue] = (unp_pos, unp_resname)
        return mapping

    def _compute_mismatch_count(self, chain: str, uniprot_id: str):
        mapping = self.mappings[chain][uniprot_id]
        mismatch_count = 0
        for structure_residue, unp_res in mapping.items():
            if self.is_mismatch(structure_residue, unp_res):
                mismatch_count += 1
        return mismatch_count

    def determine_alignment_success(self) -> bool:
        alignment_success = defaultdict(lambda: defaultdict(bool))
        for chain, unp_dict in self.mappings.items():
            for uniprot_id in unp_dict.keys():
                if self._is_alignment_successful(chain, uniprot_id):
                    alignment_success[chain][uniprot_id] = True
                else:
                    alignment_success[chain][uniprot_id] = False
        return alignment_success

    def is_alignment_successful(self, chain: str, uniprot_id: str) -> bool:
        if uniprot_id not in self.mappings[chain]:
            return False
        return self.alignment_success[chain][uniprot_id]

    def _is_alignment_successful(self, chain: str, uniprot_id: str) -> bool:
        mapping = self.mappings.get(chain, {}).get(uniprot_id, {})
        num_residues_to_align = len(mapping)
        if num_residues_to_align == 0:
            return False
        # if chain not in self.mappings:
        #     return False
        # chain_mapping = self.mappings[chain]
        # if uniprot_id not in chain_mapping:
        #     return False
        # mapping = chain_mapping[uniprot_id]
        # if len(mapping) == 0:
        #     return False
        mismatch_count = self._compute_mismatch_count(chain, uniprot_id)
        return self.is_acceptable_mismatch(mismatch_count, num_residues_to_align)

    def map(self, chain: str, structure_residue: tuple[str, str]|gemmi.Residue, uniprot_id: str, remove_mismatch: bool = True) -> tuple[int, str]:
        if chain not in self.chains_containing_uniprot_id(uniprot_id) or not self.is_alignment_successful(chain, uniprot_id):
            return None
        structure_residue = self.structure_residue_key(structure_residue)
        unp_residue = self.mappings.get(chain, {}).get(uniprot_id, {}).get(structure_residue, None)
        if remove_mismatch and unp_residue is not None:
            if self.is_mismatch(structure_residue, unp_residue):
                return None
        return unp_residue


    # @classmethod
    # def is_conistent_with_uniprot_sequence(self, uniprot_sequence: str, residue_tuple: list[tuple]|dict):
    #     best_offset, best_mismatch_count = self.best_alignment_offset(uniprot_sequence, residue_tuple)
    #     if best_mismatch_count <= 5 and len(residue_tuple) > 15:
    #         return True
    #         return False, best_offset, best_mismatch_count, best_mismatch_count/len(residue_tuple)
        
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


    @classmethod
    def best_alignment_offset(cls, uniprot_sequence: str, residues: list[tuple|gemmi.Residue]|dict, offset: int = None):
        num_residues_to_align = len(residues)
        offsets = [0, 1, -1, 2, -2, 3, -3]
        if offset is not None:
            offsets = [offset]
        best_offset = offsets[0]
        best_mismatch_count = num_residues_to_align
        if isinstance(residues, dict):
            residues = list(residues.items())
        elif isinstance(residues, list) and all(isinstance(residue, gemmi.Residue) for residue in residues):
            residues = [(str(residue.seqid), residue.name) for residue in residues]
        for offset in offsets:
            mismatch_count = 0
            for pos, resname in residues:
                resname = cls.single_letter_residue_code(resname)
                pos = cls.extract_position(pos)
                if not isinstance(pos, int):
                    continue
                unp_pos = cls.offset_adjusted_uniprot_position(pos, offset)
                if 1 <= unp_pos <= len(uniprot_sequence) and resname != uniprot_sequence[unp_pos-1]:
                    mismatch_count += 1
            if mismatch_count < best_mismatch_count:
                best_mismatch_count = mismatch_count
                best_offset = offset
            if cls.is_acceptable_mismatch(mismatch_count, num_residues_to_align):
                break
        return best_offset, best_mismatch_count

    @classmethod
    def single_letter_residue_code(cls, residue: str|tuple[str, str]|gemmi.Residue) -> str:
        if isinstance(residue, gemmi.Residue):
            return single_letter_residue_code(residue.name)
        if isinstance(residue, tuple) and len(residue) >=2:
            residue = residue[1]
        if isinstance(residue, str):
            if len(residue) == 3:
                return single_letter_residue_code(residue)
            elif len(residue) == 1:
                return residue.upper()
        warnings.warn(f"Invalid residue: {residue}")
        return None


    
    @classmethod
    def is_mismatch(cls, structure_residue: tuple[str, str], unp_residue: tuple[int, str]) -> bool:
        return cls.single_letter_residue_code(structure_residue) != unp_residue[1]


    @classmethod
    def extract_position(cls, seqid: str|gemmi.SeqId|int|None) -> int|None:
        if seqid is is_actually_a_number(seqid):
            return int(seqid)
        #seqid = str(int(seqid)) if is_actually_a_number(seqid) else seqid
        if not isinstance(seqid, gemmi.SeqId):
            try:
                seqid = gemmi.SeqId(str(seqid))
            except (ValueError, TypeError):
                #warnings.warn(f"Invalid seqid: {seqid}")
                #seqid = None
                pass
        return seqid.num if isinstance(seqid, gemmi.SeqId) else None

    @classmethod
    def offset_adjusted_uniprot_position(cls, pos: int|str, offset: int = 0) -> int:
        if isinstance(pos, str):
            pos = int(pos)
        if not isinstance(pos, int):
            return None
        return pos+offset

    @classmethod
    def is_acceptable_mismatch(cls, mismatch_count: int, num_residues_to_align: int) -> bool:
        mismatch_fraction = mismatch_count/num_residues_to_align
        return mismatch_fraction <= 0.1 and num_residues_to_align > 10 and mismatch_count < 20


    @classmethod
    def valid_seqid_str_or_None(cls, x: str|gemmi.SeqId|int|None) -> str|None:
        try:
            return str(gemmi.SeqId(str(x)))
        except (ValueError, TypeError):
            return None
        

    @classmethod
    def structure_residue_key(cls, residue_in_structure: tuple|gemmi.Residue) -> tuple[str, str]:
        res_name = cls.single_letter_residue_code(residue_in_structure)
        id = residue_in_structure.seqid if isinstance(residue_in_structure, gemmi.Residue) else None
        if isinstance(residue_in_structure, tuple) and len(residue_in_structure) >=2:
            id = residue_in_structure[0]
        id = int(id) if is_actually_a_number(id) else id
        id = cls.valid_seqid_str_or_None(id)
            # if isinstance(id, (int):
            #     return (str(id), res_name)
            # elif (isinstance(id, str) and id != ''):
            #     return (id, res_name)
        if id is None:
            warnings.warn(f"Invalid residue id: {residue_in_structure}")
        return (id, res_name)




    # @abstractmethod
    # def is_alignment_successful(self, chain: str, uniprot_id: str) -> bool:
    #     pass

    # def map(self, chain: str, 
    #               uniprot_id: str, 
    #               structure_residue: tuple[str, str]|gemmi.Residue, 
    #               remove_mismatch: bool = True) -> tuple[int, str]:
    #     return self.aligner.map(chain, uniprot_id, structure_residue, remove_mismatch)
    # @abstractmethod
    # def map(self, chain: str, 
    #               structure_residue: tuple[str, str]|gemmi.Residue, 
    #               remove_mismatch: bool = True) -> tuple[int, str]:
    #     pass

    # @abstractmethod
    # def chains_containing_uniprot(self, uniprot_id: str) -> list[str]:
    #     pass
