from dsa.binding_interfaces.factory import ClassFactory
from dsa.misc.store import Store
import gemmi
from pathlib import Path
from dsa.binding_interfaces.config import STRUCTURE_DIR

class StructureFileStore(Store):
    def __init__(self, database: str, extension: str, object_type: type = gemmi.Structure):
        self.database = database
        if "cif" not in extension and "pdb" not in extension:
            raise ValueError(f"Invalid extension: {extension}")
        super().__init__(self.get_directory(), extension, object_type=object_type)
    
    def get_directory(self):
        return STRUCTURE_DIR / self.database

    def dump_file(self, file_path: Path, data: any):
        NotImplementedError("This method is not implemented")

    def load_file(self, file_path: Path) -> any:
        return gemmi.read_structure(file_path.as_posix())

    def load_structure(self, structure_id: str) -> gemmi.Structure:
        self.load(structure_id)

    def load_block(self, structure_id: str) -> gemmi.cif.Block:
        if "cif" not in self.extension:
            raise ValueError(f"Block can only be loaded from cif files, not {self.extension} files")
        file_path = self.get_file_path(structure_id)
        if file_path.exists():
            return gemmi.cif.read(file_path.as_posix()).sole_block()
        else:
            raise FileNotFoundError(f"File {file_path} does not exist")

    def get_saved_structures(self, from_structures: list[str] = None) -> list[str]:
        return super().get_saved_keys(from_keys=from_structures)
    
    # @classmethod
    # def factory(cls, database: str):
    #     if database == "pdb":
    #         return GemmiPDBStore()
    #     elif database == "alphafold":
    #         return GemmiAlphaFoldStore()
    #     else:
    #         raise ValueError(f"Invalid database: {database}")

@ClassFactory.register(key="pdb", interface=StructureFileStore)
class PDBStructureFileStore(StructureFileStore):
    def __init__(self):
        database = "pdb"
        super().__init__(database, extension="cif", object_type=gemmi.Structure)

@ClassFactory.register(key="alphafold", interface=StructureFileStore)
class AlphaFoldStructureFileStore(StructureFileStore):
    def __init__(self):
        database = "alphafold"
        super().__init__(database, extension="pdb.gz", object_type=gemmi.Structure)

    def _key_to_file_names(self, structure_id: str) -> list[str]:
        file_name_v4 = f"AF-{structure_id}-F1-model_v4.{self.extension}"
        file_name_v6 = f"AF-{structure_id}-F1-model_v6.{self.extension}"
        return [file_name_v4, file_name_v6]

    def key_to_file_name(self, structure_id: str) -> str:
        for fn in self._key_to_file_names(structure_id):
            file_path = self.directory / fn
            if file_path.exists():
                return fn
        return fn
    
    def file_name_to_key(self, file_name: str) -> str:
        return file_name.split("-")[1]






# def get_alphafold_structure(uniprot_id, format="pdb"):
#     if format != "pdb" and format != "cif":
#         raise ValueError("format must be 'pdb' or 'cif'")
#     file_name_v4 = f"AF-{uniprot_id}-F1-model_v4.{format}.gz"
#     file_name_v6 = f"AF-{uniprot_id}-F1-model_v6.{format}.gz"
#     file_path_v4 = STRUCTURE_DIR / "alphafold" / file_name_v4
#     file_path_v6 = STRUCTURE_DIR / "alphafold" / file_name_v6
#     if file_path_v4.exists():
#         return file_path_v4
#     elif file_path_v6.exists():
#         return file_path_v6
#     else:
#         return None

# def uniprot_ids_with_saved_alphafold_structures(format="pdb"):
#     if format != "pdb" and format != "cif":
#         raise ValueError("format must be 'pdb' or 'cif'")
#     folder = STRUCTURE_DIR / "alphafold"
#     files = folder.glob(f"AF-*F1-model_v*.{format}.gz")
#     uniprot_ids = [file.stem.split("-")[1] for file in files]
#     uniprot_ids = list(set(uniprot_ids))
#     return uniprot_ids