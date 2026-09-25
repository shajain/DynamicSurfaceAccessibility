from dsa.binding_interfaces.contacts.bond_characterization.common import get_lys_nz, get_residue_atoms, atom_pos
import numpy as np
from dsa.binding_interfaces.contacts.helpers import is_lysine
from dsa.binding_interfaces.contacts.bond_characterization.library import BondLibrary

@BondLibrary.register(bond_type="salt_bridge", method="assume_lysine")
def check_salt_bridge(lys_residue, partner_residue, cutoff=4.0):
    """
    Checks NZ against carboxylate centroid of Asp/Glu, or C-terminal OXT.
    Returns dict or None.
    """


    bond_properties = {'exists': False, 'info': []}

    if not is_lysine(lys_residue):
        lys_residue, partner_residue = partner_residue, lys_residue
    
    if not is_lysine(lys_residue):
        return bond_properties

    nz = get_lys_nz(lys_residue)
    if nz['coord'] is None:
        return bond_properties

    partner_atoms = get_residue_atoms(partner_residue)
    o_names = [n for n in ('OD1', 'OD2', 'OE1', 'OE2', 'OXT') if n in partner_atoms]
    if not o_names:
        return bond_properties

    o_coords = np.array([atom_pos(partner_atoms[n]) for n in o_names])
    centroid = o_coords.mean(axis=0)
    d = np.linalg.norm(nz['coord'] - centroid)

    if d <= cutoff:
        bond_properties['exists'] = True
        bond_properties['info'].append({
            'distance': round(float(d), 2),
            'partner_atoms': o_names,
            'nz_status': nz['status'],
        })
    return bond_properties