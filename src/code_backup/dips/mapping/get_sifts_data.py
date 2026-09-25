import requests
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
from dsa.dips.mapping.search_dips_in_pdbchains import search_dips_in_pdbchains

def get_sifts_data(pdb_id):
    """Get UniProt mappings and sequences for all chains in a PDB entry."""
    # Get UniProt mapping
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    mapping_resp = requests.get(url).json()
    
    # Get SEQRES sequences per chain
    seqres_url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id.lower()}"
    seqres_resp = requests.get(seqres_url).json()
    
    return mapping_resp, seqres_resp

def pdb_sequence_to_uniprot(pdb_id, pdb_sequence):
    mapping_resp, seqres_resp = get_sifts_data(pdb_id)
    # search the pdb_sequence in the entities in seqres_resp
    matched_entity = None
    seqs = {e.get('entity_id',''): (e.get('pdb_sequence','')) for e in seqres_resp.get(pdb_id.lower(), [])}
    match, best_chain, best_chain_seq, best_score = search_dips_in_pdbchains(pdb_sequence, seqs)
    for entity in seqres_resp.get(pdb_id.lower(), []):
        score = search_dips_in_pdbchains(pdb_sequence, entity.get('pdb_sequence'))
        if entity.get('pdb_sequence') == pdb_sequence:
            matched_entity = entity
    assert matched_entity is not None
    # find the uniprot_id in the mapping_resp corresponding to the matched_entity
    uniprot_id = None
    for uniprot_id, data in mapping_resp.get(pdb_id.lower(), {}).get('UniProt', {}).items():
        for mapping in data.get('mappings', []):
            if mapping.get('chain_id') == matched_entity.get('entity_id'):
                uniprot_id = mapping.get('uniprot_id')
                break
    assert uniprot_id is not None
    # return the uniprot_id
    return None

def parse_entities(seqres_resp, pdb_id):
    """Extract {chain_id: sequence} from PDBe molecules endpoint."""
    entities = {}
    #for entity in seqres_resp.get(pdb_id.lower(), []):
        # seq = entity.get('sequence', '')
        # pdb_seq = entity.get('pdb_sequence', '')
        # entity_id = entity.get('entity_id', '')
        # entity_name = entity.get('entity_name', '')
        # source = entity.get('source', '')
        # in_chains = entity.get('in_chains', [])
        # molecule_type = entity.get('molecule_type', '')
        # gene_name = entity.get('gene_name', '')
        # ca
        # entities[entity_id] = {
        #     'pdb_sequence': pdb_seq,
        #     'entity_name': entity_name,
        #     'source': source,
        #     'in_chains': in_chains,
        # }
    return entities  # {'1': {'pdb_sequence': 'MKTAY...', 'entity_name': 'A', 'source': 'PDB', 'in_chains': ['A']}, '2': {'pdb_sequence': 'MKTAY...', 'entity_name': 'B', 'source': 'PDB', 'in_chains': ['B']}}

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


if __name__ == "__main__":
    pdb_id = "1234"
    mapping_resp, seqres_resp = get_sifts_data(pdb_id)
    print(mapping_resp)
    print(seqres_resp)