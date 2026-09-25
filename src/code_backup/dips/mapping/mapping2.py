import requests
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq

def get_sifts_data(pdb_id):
    """Get UniProt mappings and sequences for all chains in a PDB entry."""
    # Get UniProt mapping
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    mapping_resp = requests.get(url).json()
    
    # Get SEQRES sequences per chain
    seqres_url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id.lower()}"
    seqres_resp = requests.get(seqres_url).json()
    
    return mapping_resp, seqres_resp

def parse_chain_sequences(seqres_resp, pdb_id):
    """Extract {chain_id: sequence} from PDBe molecules endpoint."""
    chain_seqs = {}
    for entity in seqres_resp.get(pdb_id.lower(), []):
        seq = entity.get('sequence', '')
        for chain in entity.get('in_chains', []):
            chain_seqs[chain] = seq
    return chain_seqs  # {'A': 'MKTAY...', 'B': 'MKTAY...'}

def parse_uniprot_mapping(mapping_resp, pdb_id):
    """Extract {chain_id: {uniprot_id, residue_mapping}} from SIFTS."""
    result = {}
    uniprot_entries = mapping_resp.get(pdb_id.lower(), {}).get('UniProt', {})
    
    for uniprot_id, data in uniprot_entries.items():
        for mapping in data.get('mappings', []):
            chain = mapping['chain_id']
            if chain not in result:
                result[chain] = []
            result[chain].append({
                'uniprot_id': uniprot_id,
                'uniprot_start': mapping['unp_start'],
                'uniprot_end': mapping['unp_end'],
                'pdb_start': mapping['start']['residue_number'],
                'pdb_end': mapping['end']['residue_number'],
            })
    return result  # {'A': [{uniprot_id, uniprot_start, ...}], ...}