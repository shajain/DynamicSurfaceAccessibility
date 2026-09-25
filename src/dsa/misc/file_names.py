
class FileNames:

    @staticmethod
    def get_alphafold_file_name(uniprot_id: str, format: str="pdb") -> str:
        if format != "pdb" and format != "cif":
            raise ValueError("format must be 'pdb' or 'cif'")
        file_name_v4 = f"AF-{uniprot_id.upper()}-F1-model_v4.{format}.gz"
        file_name_v6 = f"AF-{uniprot_id.upper()}-F1-model_v6.{format}.gz"
        return [file_name_v4, file_name_v6] 

    @staticmethod
    def get_pdb_file_name(structure_id: str, format: str="cif") -> str:
        if format != "pdb" and format != "cif":
            raise ValueError("format must be 'pdb' or 'cif'")
        return f"{structure_id.lower()}.{format}"