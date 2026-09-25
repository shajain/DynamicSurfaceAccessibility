from dsa.config import PROJECT_ROOT, DATA_DIR
from pathlib import Path
from dsa.uniprot.config import ALL_UNIPROT_IDS

STRUCTURE_DIR = DATA_DIR / "structures" 

BINDING_INTERFACES_DIR = DATA_DIR / "binding_interfaces"


PDB_DIR = STRUCTURE_DIR / "pdb"
test_structure_id = "3cc0"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
PDBE_UNIPROT_MAP_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{uniprot_id}"
DEFAULT_HEADERS = {
    "User-Agent": "DynamicSurfaceAccessibility/1.0 (https://github.com/)",
    "Accept": "application/json",
}



#SOURCES = ["pdb", "AlphaFold2", "AlphaFold3", "ESM3"]
SOURCES = ["pdb"]

#Specified as angstroms * 10
#CUTOFF_ANGSTROM = [45, 60, 80]
CUTOFF_ANGSTROM = [45]

UNIPROT_STRUCTURE_MAPPINGS_DIR = DATA_DIR / "uniprot_structure_mappings"
PROTEIN_CONTACTS_DIR = BINDING_INTERFACES_DIR / "protein_contacts"
STRUCTURE_CONTACTS_DIR = BINDING_INTERFACES_DIR / "structure_contacts"
STRUCTURE_INTERFACES_DIR = BINDING_INTERFACES_DIR / "structure_interfaces"


# uniprot_to_PDB_ids_file = UNIPROT_STRUCTURE_MAPPINGS_DIR / "uniprot_to_PDB_ids.json"
# uniprot_to_alphafold_ids_file = UNIPROT_STRUCTURE_MAPPINGS_DIR / "uniprot_to_alphafold_ids.json"
small_uniprot_to_structureIDs_file = {db: BINDING_INTERFACES_DIR / f"small_uniprot_to_{db}_structureIDs.json" for db in SOURCES}
ALIGNMENT_DIR = STRUCTURE_DIR / "alignment"
PDBE_ALIGNMENT_DIR = ALIGNMENT_DIR / "pdbe"
SIFT_alignment_file =  "pdb_chain_uniprot.tsv"

ALPHAFOLD_STRUCTURES_TO_PROCESS = ALL_UNIPROT_IDS

BOND_PROPERTIES = ["salt_bridge", "hydrogen_bond", "cation_pi", "covalent_isopeptide"]