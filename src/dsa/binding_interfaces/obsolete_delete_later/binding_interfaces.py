from dsa.ppi_binding_interfaces.config import SOURCES, uniprot_to_strutureIDs_file, STRUCTURE_DIR, BINDING_INTERFACES_DIR, CUTOFF_ANGSTROM
import json
from pathlib import Path

class BindingInterface:
    def __init__(self, sources: str = SOURCES, structure_dir: str = STRUCTURE_DIR, binding_interfaces_dir: str = BINDING_INTERFACES_DIR, uniprot_mappings_file: str = uniprot_to_strutureIDs_file):
        self.sources = sources
        self.structure_dir = structure_dir
        self.binding_interfaces_dir = binding_interfaces_dir
        self.uniprot_mappings_file = uniprot_mappings_file
        self.uniprot_to_strutureIDs = self._load_uniprot_to_strutureIDs()

    def _load_uniprot_to_strutureIDs(self):
        if not self.uniprot_mappings_file.exists():
            print(f"Uniprot mappings file not found at {self.uniprot_mappings_file}")
            print(f"Creating an empty uniprot mappings file at {self.uniprot_mappings_file}")
            with open(self.uniprot_mappings_file, "w") as f:
                json.dump({}, f)
            return {}
        with open(self.uniprot_mappings_file, "r") as f:
            return json.load(f)
    

    def extract_binding_interfaces(self, cutoff_angstrom: int = CUTOFF_ANGSTROM):
        for source in self.sources:
            SOURCE_DIR = self.structure_dir / source
            #read all cif files in SOURCE_DIR
            for cif_file in SOURCE_DIR.glob("*.cif"):
                #extract the binding interfaces
                binding_interfaces = self._extract_binding_interfaces(cif_file, cutoff_angstrom)
                #save the binding interfaces
                self._save_binding_interfaces(binding_interfaces, source)

    def _extract_binding_interfaces(self, cif_file: Path, cutoff_angstrom: int):
        #read the cif file
        with open(cif_file, "r") as f:
            cif_data = f.read()
        #extract the binding interfaces
        binding_interfaces = self._extract_binding_interfaces(cif_data, cutoff_angstrom)