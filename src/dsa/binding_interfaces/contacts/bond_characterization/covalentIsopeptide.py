from dsa.binding_interfaces.contacts.bond_characterization.common import get_lys_nz, get_residue_atoms, atom_pos
import numpy as np
from dsa.binding_interfaces.contacts.helpers import is_lysine
from dsa.binding_interfaces.contacts.bond_characterization.library import BondLibrary

# def check_covalent_isopeptide(lys_residue, partner_residue, cutoff=1.6):
#     """
#     Checks NZ against a carbonyl carbon (backbone C or Gln side-chain CD)
#     from the partner residue — candidate isopeptide/amide linkage.
#     Uses only the OBSERVED NZ position (reconstructed NZ is meaningless here,
#     since a true covalent bond means NZ position is fixed by that bond, not
#     by rotamer geometry).
#     """
#     nz = get_lys_nz(lys_residue)
#     if nz['coord'] is None or nz['status'] != 'observed':
#         return None

#     partner_atoms = get_residue_atoms(partner_residue)
#     hits = []
#     for name in ('C', 'CD'):
#         if name in partner_atoms:
#             d = np.linalg.norm(nz['coord'] - atom_pos(partner_atoms[name]))
#             if d <= cutoff:
#                 hits.append({
#                     'exists': True,
#                     'distance': round(float(d), 2),
#                     'partner_resname': partner_residue.name,
#                     'partner_atom': name,
#                     'nz_status': nz['status'],
#                 })
#     return hits or None


@BondLibrary.register(bond_type="covalent_isopeptide", method="assume_lysine")
def check_covalent_isopeptide(lys_residue, partner_residue, cutoff=1.6):


    bond_properties = {'exists':False, 'info': []}

    if not is_lysine(lys_residue):
        lys_residue, partner_residue = partner_residue, lys_residue
    
    if not is_lysine(lys_residue):
        return bond_properties

    nz = get_lys_nz(lys_residue)
    if nz["coord"] is None or nz["status"] != "observed":
        return bond_properties

    partner_atoms = get_residue_atoms(partner_residue)
    best_hit = None
    min_dist = cutoff

    # Check potential carbonyl carbons (Backbone C, Glutamate CG/CD, Aspartate CG, Gln CD)
    for name in ("C", "CD", "CG"):
        if name in partner_atoms:
            d = float(np.linalg.norm(nz["coord"] - atom_pos(partner_atoms[name])))
            if d <= min_dist:
                bond_properties['exists'] = True
                min_dist = d
                bond_properties['info'].append({
                    'distance': round(d, 2),
                    'partner_atom': name,
                    'nz_status': nz['status'],
                })

    return bond_properties
