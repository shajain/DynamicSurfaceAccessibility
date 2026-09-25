import numpy as np
from dsa.binding_interfaces.contacts.bond_characterization.common import get_lys_nz, get_residue_atoms, atom_pos
from dsa.binding_interfaces.contacts.helpers import is_lysine
from dsa.binding_interfaces.contacts.bond_characterization.library import BondLibrary
# HBOND_ACCEPTOR_ATOMS = {
#     'O', 'OG', 'OG1', 'OH', 'OD1', 'OD2', 'OE1', 'OE2',
#     'ND1', 'NE2', 'SG', 'OXT',
# }


PROTEIN_ACCEPTORS = {
    # Backbone
    "O",
    "OXT",  # C-terminal oxygens
    # Side chains
    "OD1",
    "OD2",  # Asp
    "OE1",
    "OE2",  # Glu, Gln, Asn
    "OG",
    "OG1",  # Ser, Thr
    "OH",  # Tyr
    "ND1",
    "NE2",  # His
    "SG",  # Cys (weak, but acts as acceptor)
}

NUCLEIC_ACCEPTORS = {
    # Phosphate & Sugar Backbone
    "OP1",
    "OP2",
    "OP3",
    "O3'",
    "O4'",
    "O5'",
    "O3*",
    "O4*",
    "O5*",
    # Nitrogenous Bases
    "N1",
    "N3",
    "N7",
    "O2",
    "O4",
    "O6",
}

ALL_ACCEPTORS = PROTEIN_ACCEPTORS | NUCLEIC_ACCEPTORS


@BondLibrary.register(bond_type="hydrogen_bond", method="assume_lysine")
def check_hydrogen_bonds(lys_residue, partner_residue, cutoff=3.5):
    """Distance-based proxy for hydrogen bonding between Lysine NZ and an acceptor residue.

    NZ (-NH3+) can donate up to 3 hydrogen bonds simultaneously.
    """

    
    
    bond_properties = {'exists': False, 'info': []}

    if not is_lysine(lys_residue):
        lys_residue, partner_residue = partner_residue, lys_residue
    
    if not is_lysine(lys_residue):
        return bond_properties

    hbond_acceptors = ALL_ACCEPTORS

    nz = get_lys_nz(lys_residue)
    if nz is None or nz.get("coord") is None:
        return bond_properties

    nz_pos = np.array(nz["coord"], dtype=float)
    partner_atoms = get_residue_atoms(partner_residue)
    hits = []

    for name, atom in partner_atoms.items():
        # Clean atom name if altloc keys exist (e.g., 'O_A' -> 'O')
        clean_name = name.split("_")[0] if "_" in name else name

        if clean_name in hbond_acceptors:
            pos = np.array(atom_pos(atom), dtype=float)
            d = float(np.linalg.norm(nz_pos - pos))

            if d <= cutoff:
                bond_properties['exists'] = True
                bond_properties['info'].append({
                    'distance': round(d, 2),
                    'partner_atom': name,
                    'nz_status': nz['status'],
                })
    return bond_properties


# def check_hydrogen_bonds(lys_residue, partner_residue, cutoff=3.5):
#     """
#     Distance-only proxy (no explicit H). Returns list of dicts (NZ can donate
#     to multiple acceptors simultaneously via its 3 H's).
#     """
#     nz = get_lys_nz(lys_residue)
#     if nz['coord'] is None:
#         return []

#     partner_atoms = get_residue_atoms(partner_residue)
#     hits = []
#     for name, atom in partner_atoms.items():
#         if name in HBOND_ACCEPTOR_ATOMS:
#             d = np.linalg.norm(nz['coord'] - atom_pos(atom))
#             if d <= cutoff:
#                 hits.append({
#                     'exists': True,
#                     'distance': round(float(d), 2),
#                     'partner_resname': partner_residue.name,
#                     'partner_atom': name,
#                     'nz_status': nz['status'],
#                 })
#     return hits

