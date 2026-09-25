from Bio.SeqUtils import seq1
import gemmi
from dsa.binding_interfaces.factory import ClassFactory
#from dsa.binding_interfaces.structure.structure import Structure
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
from dsa.binding_interfaces.structure_uniprot_alignment.utilities import extract_num_str, single_letter_resname, int_if_str, str2int
import warnings
from dsa.binding_interfaces.structure_uniprot_alignment.aligner import StructureSequenceAligner
from dsa.binding_interfaces.structure.structure_model import StructureModel




@ClassFactory.register(key="alphafold", interface=StructureSequenceAligner)
class AlphaFoldSequenceAligner(StructureSequenceAligner):
    def __init__(self, structure_id: str, 
                        structure_model: StructureModel, 
                        uniprot_to_entity: dict[str, list[str]]):
        super().__init__(structure_id, structure_model, uniprot_to_entity)
    # def __init__(self, structure: Structure):
    #     super().__init__(structure)

# @ClassFactory.register(key="alphafold", interface=StructureSequenceAligner)
# class AlphaFoldSequenceAligner(StructureSequenceAligner):
#     def __init__(self, structure_id: str, 
#                         structure_model: StructureModel, 
#                         uniprot_to_entity: dict[str, list[str]]):
#         self.structure_id = structure_id
#         self.structure_model = structure_model
#         self.uniprot_to_entity = uniprot_to_entity
#         self.structure = structure
#         self.chain = self.structure.chains[0]
#         self.uniprot_id = structure.uniprot_ids[0]
#         self.structure_sequence = self.get_structure_sequence()
#         self.uniprot_sequence = UniprotToSequence.get_sequence(self.uniprot_id)
#         self.alignment_success = self._is_alignment_successful(self.chain, self.uniprot_id)

#     def get_structure_sequence(self, chain: str|None = None) -> str:
#         chain = chain if chain is not None else self.chain
#         alphafold_sequences = self.structure.chain_sequence_grouped_by_uniprot()
#         structure_sequence = alphafold_sequences.get(self.uniprot_id, {}).get(chain, [])
#         if len(structure_sequence) == 0:
#             warnings.warn(f"No structure sequence found for chain {chain} and UniProt ID {self.uniprot_id}")
#         return structure_sequence

#     def is_alignment_successful(self) -> bool:
#         return self.alignment_success

#     def _is_alignment_successful(self, chain: str, uniprot_id: str) -> bool:
#         if uniprot_id != self.uniprot_id or chain != self.chain or len(self.structure_sequence) == 0:
#             return False
#         #structure_sequence = [(res.seqid.num, res.name) for res in self.structure.structure[0][chain] if res.entity_type == gemmi.EntityType.Polymer]
#         count = 0
#         structure_sequence_length = len(self.structure_sequence)
#         for residue_in_structure in self.structure_sequence:
#             if not self.is_match(residue_in_structure):
#                 warnings.warn(f"Mismatch found between UniProt and AlphaFold for chain {chain} and UniProt ID {uniprot_id}")
#                 count += 1
#         fraction_mismatched = count / structure_sequence_length
#         return fraction_mismatched <= 0.1

#     def map(self, chain: str, uniprot_id: str, structure_residue: tuple[str, str]|gemmi.Residue, remove_mismatch: bool = True) -> tuple[int, str]:
#         if chain != self.chain or uniprot_id != self.uniprot_id:
#             return None
#         if remove_mismatch and not self.is_match(structure_residue):
#             return None
#         pos, name = self.convert_residue_to_uniprot_format(structure_residue)
#         return (pos, name)

#     def convert_residue_to_uniprot_format(self, structure_residue: tuple[str, str]|gemmi.Residue) -> tuple[int, str]:
#         if isinstance(structure_residue, gemmi.Residue):
#             structure_residue = (structure_residue.seqid.num, structure_residue.name)
#         name = single_letter_resname(structure_residue[1])
#         pos = int_if_str(extract_num_str(structure_residue[0]))
#         return (pos, name)

#     def is_match(self, structure_residue: tuple[str, str]|gemmi.Residue) -> bool:
#         pos, name = self.convert_residue_to_uniprot_format(structure_residue)
#         return self.uniprot_sequence[pos-1] == name

#     def chains_containing_uniprot(self, uniprot_id: str) -> list[str]:
#         chains = [c.name for c in self.structure.structure[0]]
#         return chains


