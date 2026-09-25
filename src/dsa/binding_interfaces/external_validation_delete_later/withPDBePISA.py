"""
Compare your computed interface contacts (mapped to UniProt positions)
against PDBe-KB precomputed interface residues (also at UniProt level).

PDBe-KB uses PISA under the hood and aggregates interface data across
ALL PDB structures for a given UniProt accession.

Usage:
    python compare_interface_residues.py

Requirements:
    pip install requests pandas
"""

import requests
import pandas as pd
from collections import defaultdict
from dsa.uniprot.config import ALL_UNIPROT_IDS
from dsa.binding_interfaces.binding_interfaces import BindingInterfacesStore
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
from Bio.SeqUtils import seq1

PDBE_API = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues/{accession}"



def fetch_pdbekb_interface_residues(accession: str) -> pd.DataFrame:
    """
    Returns a DataFrame of interface residues at UniProt positions for a
    given UniProt accession, sourced from PDBe-KB (PISA-based).

    Columns:
        uniprot_pos   : UniProt residue position
        pdb_id        : PDB entry the interface was observed in
        chain_id      : Chain in that PDB entry
        partner_type  : Type of interaction partner (protein, ligand, etc.)
        partner_id    : Identity of the partner (chain or ligand name)
    """
    url = PDBE_API.format(accession=accession)
    response = requests.get(url, timeout=30)

    if response.status_code != 200:
        print(f"  [!] No data returned for {accession} (status {response.status_code})")
        return pd.DataFrame()

    data = response.json()

    if accession not in data:
        print(f"  [!] Accession {accession} not found in response.")
        return pd.DataFrame()

    # rows = []
    # for entry in data[accession].get("data", []):
    #     partner_type = entry.get("accession", "unknown")   # e.g. protein, ligand
    #     pdb_id       = entry.get("pdb_id", "")
    #     chain_id     = entry.get("chain_id", "")
    #     partner_id   = entry.get("interaction_partner_id", "")

    #     for residue in entry.get("residues", []):
    #         rows.append({
    #             "uniprot_pos":  residue.get("startIndex"),
    #             "pdb_id":       pdb_id,
    #             "chain_id":     chain_id,
    #             "partner_type": partner_type,
    #             "partner_id":   partner_id,
    #         })

    # df = pd.DataFrame(rows)
    return data


class ContactsValidator:
    def __init__(self, uniprot_id: str):
        self.uniprot_id = uniprot_id
        data = fetch_pdbekb_interface_residues(uniprot_id)
        data = data[uniprot_id]
        self.data = data.get("data", [])
        self.sequence = data.get("sequence", "")
        self.length = len(self.sequence)

    def get_ppi_residues(self) -> list[tuple[dict]]:
        contacts = []
        for entry in self.data[self.uniprot_id].get("data", []):
            unp = entry['accession']
            for res in entry['residues']:
                if res['indexType'].lower() == 'uniprot':
                    start = res['startIndex']
                    end = res['endIndex']
                    assert seq1(res['startCode']) == self.sequence[start-1]
                    assert seq1(res['endCode']) == self.sequence[end-1]
                    interacting_pdb = res['interactingPDBEntries']
                    for i in range(start, end):
                        for interacting_pdb in interacting_pdb:
                            for chain in interacting_pdb['chainIds']:
                                contact_info = {
                                    'residue': (i, self.sequence[i-1]),
                                    'uniprot': unp,
                                    'structure_id': interacting_pdb['pdbId'],
                                    'chain': chain,   
                                }
                                contacts.append(contact_info)
        return contacts

    def get_residue

    @classmethod
    def validate_contacts(cls, uniprot_id: str, binding_interfaces: BindingInterfacesStore) -> list[dict]:
        contacts1 = cls(uniprot_id).get_ppi_residues()

        contacts in self.get_ppi_residues():
            if contact in binding_interfaces:
                return True
        return False






# ─────────────────────────────────────────────
# 4.  OPTIONAL: SAVE FULL PDBE-KB DATA TO CSV
# ─────────────────────────────────────────────

def save_pdbekb_data(accession: str, df: pd.DataFrame):
    if df.empty:
        return
    fname = f"{accession}_pdbekb_interfaces.csv"
    df.to_csv(fname, index=False)
    print(f"\n  Saved full PDBe-KB interface data to: {fname}")


# ─────────────────────────────────────────────
# 5.  MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    results = {}

    for accession in ALL_UNIPROT_IDS[:10]:
        print(f"\nFetching PDBe-KB interface residues for {accession} ...")
        df_pdbe = fetch_pdbekb_interface_residues(accession)
        print(df_pdbe)
        binding_interfaces = BindingInterfacesStore(structure_db="pdb", cutoff=8.0).load(accession)
        seq1 = UniprotToSequence.get_sequence(accession)
        seq2 = 
        # your_contacts = YOUR_CONTACTS.get(accession, set())
        # result = compare_contacts(accession, your_contacts, df_pdbe)

        # if result:
        #     results[accession] = result
        #     save_pdbekb_data(accession, df_pdbe)
    
    print("\nDone.")