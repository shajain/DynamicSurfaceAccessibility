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

uniprot_to_strutureIDs_file = BINDING_INTERFACES_DIR / "uniprot_to_structureIDs.json"
BI_PROTEIN_DIR = BINDING_INTERFACES_DIR / "protein_wise_contacts"
BI_STRUCTURE_DIR = BINDING_INTERFACES_DIR / "structure_wise_contacts"

small_uniprot_to_structureIDs_file = {db: BINDING_INTERFACES_DIR / f"small_uniprot_to_{db}_structureIDs.json" for db in SOURCES}
PDB_UNIPROT_MAP_DIR = PDB_DIR / "pdb_uniprot_mappings"
SIFT_PDB_to_uniprot_file = PDB_UNIPROT_MAP_DIR / "pdb_chain_uniprot.tsv.gz"