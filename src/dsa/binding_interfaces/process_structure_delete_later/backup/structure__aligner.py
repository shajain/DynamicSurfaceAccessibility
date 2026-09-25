from dsa.binding_interfaces.structure_alignment.pdb_aligner import PDB_Uniprot_Mapper
from dsa.binding_interfaces.structure_alignm.alphafold_aligner import AlphaFold_Uniprot_Mapper



class StructureSequenceAligner:
    def __init__(self, ):
        self.structure_db = structure_db
        self.structure_id = structure_id
        if self.structure_db == "pdb":
            self.aligner = PDBSequenceAligner(self.structure_id, structure_db=self.structure_db)
        elif self.structure_db == "alphafold":
            self.aligner = AlphaFoldSequenceAligner(self.structure_id, structure_db=self.structure_db)
        else:
            raise ValueError(f"Invalid structure database: {self.structure_db}")

    def is_alignment_successful(self, chain: str, uniprot_id: str) -> bool:
        return self.aligner.is_alignment_successful(chain, uniprot_id, self.structure_id)

    def map(self, chain: str, uniprot_id: str, residue_in_structure: tuple[str, str], remove_mismatch: bool = True) -> tuple[int, str]:
        return self.aligner.map(chain, uniprot_id, residue_in_structure, remove_mismatch, self.structure_id)

    def chains_containing_uniprot(self, uniprot_id: str) -> list[str]:
        return self.aligner.chains_containing_uniprot(uniprot_id, self.structure_id)