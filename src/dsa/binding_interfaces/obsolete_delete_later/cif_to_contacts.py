import gemmi
import sys
from collections import defaultdict

def extract_uniprot_ids(doc):
    """Extract chain -> UniProt ID mapping from mmCIF metadata."""
    block = doc.sole_block()

    # entity_id -> uniprot accession
    entity_to_uniprot = {}
    try:
        struct_ref = block.find(
            '_struct_ref.',
            ['entity_id', 'db_name', 'db_code', 'pdbx_db_accession']
        )
        for row in struct_ref:
            if row[1].upper() == 'UNP':  # UniProt
                entity_to_uniprot[row[0]] = row[3] if row[3] != '?' else row[2]
    except Exception as e:
        print(f"Warning: could not read _struct_ref: {e}")

    # chain_id -> entity_id
    chain_to_uniprot = {}
    try:
        struct_asym = block.find('_struct_asym.', ['id', 'entity_id'])
        for row in struct_asym:
            chain_id = row[0]
            entity_id = row[1]
            uniprot = entity_to_uniprot.get(entity_id, None)
            chain_to_uniprot[chain_id] = uniprot
    except Exception as e:
        print(f"Warning: could not read _struct_asym: {e}")

    return chain_to_uniprot


def find_chain_contacts(structure, cutoff_angstrom=5.0):
    """Find residue-level contacts between different chains."""
    model = structure[0]  # First model

    # Build neighbor search over all atoms
    ns = gemmi.NeighborSearch(model, structure.cell, cutoff_angstrom).populate()

    contacts = defaultdict(set)  # (chain_A, chain_B) -> set of (resA, resB) pairs

    for chain in model:
        for residue in chain:
            for atom in residue:
                marks = ns.find_neighbors(atom, min_dist=0.1, max_dist=cutoff_angstrom)
                for mark in marks:
                    cra = mark.to_cra(model)
                    other_chain = cra.chain.name

                    if other_chain == chain.name:
                        continue  # skip same-chain contacts

                    # Canonical key: always sort chain pair
                    pair_key = tuple(sorted([chain.name, other_chain]))
                    res_a = f"{residue.name}{residue.seqid}"
                    res_b = f"{cra.residue.name}{cra.residue.seqid}"

                    if chain.name < other_chain:
                        contacts[pair_key].add((str(res_a), str(res_b)))
                    else:
                        contacts[pair_key].add((str(res_b), str(res_a)))

    return contacts

import gemmi

def extract_sequences(doc):
    block = doc.sole_block()

    # 1. Deposited sequences per entity
    print("=== Deposited sequences (per entity) ===")
    try:
        poly = block.find('_entity_poly.', 
                          ['entity_id', 'pdbx_strand_id', 
                           'pdbx_seq_one_letter_code_can'])
        for row in poly:
            entity_id, chains, seq = row[0], row[1], row[2]
            seq_clean = seq.replace('\n', '').replace(' ', '')
            print(f"Entity {entity_id} (chains: {chains})")
            print(f"  Length: {len(seq_clean)}")
            print(f"  Seq: {seq_clean[:60]}{'...' if len(seq_clean) > 60 else ''}")
    except Exception as e:
        print(f"Warning: {e}")

    # 2. Alignment to UniProt (what range of UniProt is covered)
    print("\n=== UniProt alignment ranges ===")
    try:
        struct_ref = block.find('_struct_ref.', 
                                ['id', 'entity_id', 'db_name', 
                                 'pdbx_db_accession'])
        ref_map = {row[0]: (row[1], row[2], row[3]) for row in struct_ref}

        seq_align = block.find('_struct_ref_seq.',
                               ['ref_id', 'pdbx_auth_seq_align_beg',
                                'pdbx_auth_seq_align_end',
                                'db_align_beg', 'db_align_end'])
        for row in seq_align:
            ref_id = row[0]
            entity_id, db_name, accession = ref_map.get(ref_id, ('?','?','?'))
            print(f"  Entity {entity_id} | {db_name} {accession}")
            print(f"    PDB residues {row[1]}–{row[2]}  "
                  f"→  UniProt positions {row[3]}–{row[4]}")
    except Exception as e:
        print(f"Warning: {e}")

    # 3. Mutations/conflicts vs UniProt
    print("\n=== Sequence differences vs UniProt ===")
    try:
        dif = block.find('_struct_ref_seq_dif.',
                         ['pdbx_pdb_strand_id', 'pdbx_auth_seq_num',
                          'db_mon_id', 'mon_id', 'details'])
        for row in dif:
            chain, resnum, uniprot_res, pdb_res, detail = row
            print(f"  Chain {chain} pos {resnum}: "
                  f"UniProt={uniprot_res} → PDB={pdb_res} ({detail})")
    except Exception as e:
        print(f"Warning: {e}")

def main(mmcif_path, cutoff=5.0):
    print(f"Loading {mmcif_path}...")
    doc = gemmi.cif.read(mmcif_path)
    structure = gemmi.make_structure_from_block(doc.sole_block())
    structure.setup_entities()

    # --- UniProt IDs ---
    print("\n=== Chain → UniProt IDs ===")
    chain_to_uniprot = extract_uniprot_ids(doc)
    for chain, uniprot in sorted(chain_to_uniprot.items()):
        label = uniprot if uniprot else "No UniProt (ligand/water/unknown)"
        print(f"  Chain {chain}: {label}")

    # --- Inter-chain contacts ---
    print(f"\n=== Inter-chain contacts (cutoff: {cutoff} Å) ===")
    contacts = find_chain_contacts(structure, cutoff_angstrom=cutoff)

    for (chain_a, chain_b), residue_pairs in sorted(contacts.items()):
        uniprot_a = chain_to_uniprot.get(chain_a, "?")
        uniprot_b = chain_to_uniprot.get(chain_b, "?")
        print(f"\n  {chain_a} ({uniprot_a})  ↔  {chain_b} ({uniprot_b})")
        print(f"  {len(residue_pairs)} contacting residue pairs")
        # Print first 5 as a sample
        for res_a, res_b in sorted(residue_pairs)[:5]:
            print(f"    {res_a}  —  {res_b}")
        if len(residue_pairs) > 5:
            print(f"    ... and {len(residue_pairs) - 5} more")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "structure.cif"
    cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
    main(path, cutoff)