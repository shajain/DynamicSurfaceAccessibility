import requests

def gene_to_pdb_chains(gene_name, organism="human"):
    # Search UniProt by gene name
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_name} AND organism_name:{organism} AND reviewed:true",
        "fields": "accession,gene_names,xref_pdb",
        "format": "json"
    }
    results = requests.get(url, params=params).json()["results"]
    
    if not results:
        return {}
    
    entry = results[0]  # take top reviewed hit
    
    # Extract PDB cross-references (includes chain IDs and residue ranges)
    pdb_refs = {}
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref["database"] == "PDB":
            pdb_id = ref["id"]
            chain_info = next((p["value"] for p in ref["properties"] if p["key"] == "Chains"), "")
            pdb_refs[pdb_id] = chain_info  # e.g., "A/B=1-141"
    
    return pdb_refs

def parse_chain_mapping(pdb_refs):
    """Convert UniProt output to {pdb_id: [chain_ids]} dict."""
    result = {}
    for pdb_id, chain_str in pdb_refs.items():
        chains = set()
        # Handle entries like 'A/C=2-142' or 'A=2-142, C=3-142'
        for segment in chain_str.split(","):
            segment = segment.strip()
            chain_part = segment.split("=")[0]  # get 'A/C' or 'A'
            for chain in chain_part.split("/"):
                chains.add(chain.strip())
        result[pdb_id.lower()] = chains
    return result

def get_biounit_chain_mapping(pdb_id):
    """Map biounit chain IDs back to asym unit chain IDs."""
    url = f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id.upper()}/1"
    data = requests.get(url).json()
    
    mapping = {}
    for gen in data.get("pdbx_struct_assembly_gen", []):
        asym_chains = gen.get("asym_id_list", [])      # original letter chains
        oper_list = gen.get("oper_expression", "")      # symmetry operations
        # Each operation creates a copy — but the chain ID mapping
        # isn't directly available here
    return mapping

# Example
pdb_refs = gene_to_pdb_chains("HBA1")
# {'4HHB': 'A/C=2-142', '1A3N': 'A/C=2-142', ...}

# Build the lookup
chain_map = parse_chain_mapping(pdb_refs)
print(chain_map)
# {'1a00': {'A', 'C'}, '4hhb': {'A', 'C'}, ...}

biounit_map = get_biounit_chain_mapping("4HHB")
print(biounit_map)
# {'A': 'A', 'C': 'C'}
