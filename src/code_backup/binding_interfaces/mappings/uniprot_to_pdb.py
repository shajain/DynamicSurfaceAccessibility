from Bio.PDB import PDBList
import requests
from dsa.binding_interfaces.config import PDB_DIR, RCSB_SEARCH_URL, DEFAULT_HEADERS
from pathlib import Path
import pandas as pd
import time
from dsa.binding_interfaces.mappings.uniprot_mappings import save_uniprot_to_strutureIDs, get_mapped_uniprot_ids
from dsa.ms_dsa.ms_dsa import ALL_UNIPROT_IDS
import json



def get_pdb_ids_from_rcsb(
    uniprot_id: str,
    *,
    include_computational_models: bool = True,
) -> list[str]:
    """Return PDB entry IDs for ``uniprot_id`` using RCSB Search API v2.

    Set ``include_computational_models=False`` to list only experimental entries (works with
    :class:`Bio.PDB.PDBList`); ``True`` also returns e.g. AlphaFold DB identifiers.
    """
    sess = requests.Session()
    
    content_types = (
        ["experimental", "computational"]
        if include_computational_models
        else ["experimental"]
    )
    body = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                        "operator": "exact_match",
                        "value": "UniProt",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": content_types,
            "return_all_hits": True,
        },
    }
    got_pdbids = False
    attempts = 0
    pdbids = []
    while not got_pdbids and attempts < 2:
        r = sess.post(RCSB_SEARCH_URL, json=body, headers=DEFAULT_HEADERS, timeout=60)
        r.raise_for_status()
        try:
            payload = r.json()
            result_set = payload.get("result_set") or []
            pdbids = [result.get("identifier") for result in result_set]
            print(f"Found {len(pdbids)} PDB IDs for {uniprot_id}")
            got_pdbids = True
        except Exception as e:
            print(f"Error parsing JSON: {e} for attempt {attempts} for {uniprot_id}")
            #time.sleep(3)
            attempts += 1
            continue
    return pdbids


# def _pdb_ids_from_pdbe(uniprot_id: str, session: requests.Session | None = None) -> list[str] | None:
#     """Return PDB IDs from PDBe UniProt mapping, or None if unavailable or empty."""
#     sess = session or requests.Session()
#     url = PD_BE_UNIPROT_MAP_URL.format(uniprot_id=uniprot_id)
#     response = sess.get(url, headers=DEFAULT_HEADERS, timeout=60)
#     if response.status_code != 200:
#         return None
#     data = response.json()
#     if uniprot_id not in data or not isinstance(data[uniprot_id], dict):
#         return None
#     return list(data[uniprot_id].keys())


def download_structures_for_uniprot_id(uniprot_id, download_dir=PDB_DIR, redownload=False):
    pdbids = get_pdb_ids_from_rcsb(uniprot_id, include_computational_models=False)
    if not redownload:
        pdbids = remove_downloaded_ids(pdbids, download_dir)
    failed_pdbids = []
    for pdb_id in pdbids:
        pid = download_structure_from_pdb(pdb_id, download_dir)
        if pid is None:
            failed_pdbids.append(pdb_id)
            continue
        # try:
        #     pdbl.retrieve_pdb_file(pdb_id, pdir=download_dir, file_format="mmCif")
        # except Exception as e:
        #     print(f"Error downloading {pdb_id}: {e}")
        #     continue
    pdbids = [pdb_id for pdb_id in pdbids if pdb_id not in failed_pdbids]
    return pdbids

def download_structure_from_pdb(pdb_id, download_dir=PDB_DIR):
    pdbl = PDBList()
    try:
        pdbl.retrieve_pdb_file(pdb_id, pdir=download_dir, file_format="mmCif")
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        return None
    return pdb_id

def save_mapping_to_file(uniprot_ids):
    uniprot_mapping = {}
    mapped_uniprot_ids = get_mapped_uniprot_ids("pdb")
    uniprot_ids_remaining = [uniprot_id for uniprot_id in uniprot_ids if uniprot_id not in mapped_uniprot_ids]
    print(f"mapping {len(uniprot_ids_remaining)} uniprot ids to pdb ids")
    time.sleep(3)
    for uniprot_id in uniprot_ids_remaining:
        pdbids = get_pdb_ids_from_rcsb(uniprot_id, include_computational_models=False)
        uniprot_mapping[uniprot_id] = pdbids
    save_uniprot_to_strutureIDs(uniprot_mapping, "pdb")

def remove_downloaded_ids(pdb_ids, download_dir=PDB_DIR):
    loaded_pdbids = get_loaded_pdbStructureIDs(download_dir)
    return [pdb_id for pdb_id in pdb_ids if pdb_id.lower() not in loaded_pdbids]

def get_loaded_pdbStructureIDs(download_dir=PDB_DIR):
    pdb_ids = []
    for pdb_file in download_dir.glob("*.cif"):
        pdb_id = pdb_file.stem
        pdb_ids.append(pdb_id)
    return pdb_ids




if __name__ == "__main__":
    # Get all the UniProt IDs from all the MS-DSA files
    save_mapping_to_file(uniprot_ids=ALL_UNIPROT_IDS)
    [download_structures_for_uniprot_id(uniprot_id) for uniprot_id in ALL_UNIPROT_IDS]
