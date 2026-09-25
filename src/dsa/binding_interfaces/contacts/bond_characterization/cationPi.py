from dsa.binding_interfaces.contacts.bond_characterization.common import get_lys_nz, get_residue_atoms, atom_pos
import numpy as np
from dsa.binding_interfaces.contacts.helpers import is_lysine
from dsa.binding_interfaces.contacts.bond_characterization.library import BondLibrary

RING_ATOMS = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
}

@BondLibrary.register(bond_type="cation_pi", method="assume_lysine")
def check_cation_pi(lys_residue, partner_residue, dist_cutoff=6.0, angle_max=40.0):
    """
    Checks NZ position relative to aromatic ring centroid + normal.
    Returns dict or None.
    """     
    bond_properties = {'exists': False, 'info': []}

    if not is_lysine(lys_residue):
        lys_residue, partner_residue = partner_residue, lys_residue
    
    if not is_lysine(lys_residue):
        return bond_properties

    if partner_residue.name not in RING_ATOMS:
        return bond_properties

    nz = get_lys_nz(lys_residue)
    if nz['coord'] is None:
        return bond_properties

    partner_atoms = get_residue_atoms(partner_residue)
    names = RING_ATOMS[partner_residue.name]
    if not all(n in partner_atoms for n in names):
        return bond_properties  # ring incomplete, skip rather than guess

    coords = np.array([atom_pos(partner_atoms[n]) for n in names])
    centroid = coords.mean(axis=0)
    v1, v2 = coords[1] - coords[0], coords[2] - coords[0]
    normal = np.cross(v1, v2)
    normal /= np.linalg.norm(normal)

    d = np.linalg.norm(nz['coord'] - centroid)
    if d > dist_cutoff:
        return bond_properties

    v = (nz['coord'] - centroid)
    v_norm = v / np.linalg.norm(v)
    angle_from_normal = np.degrees(np.arccos(np.clip(abs(np.dot(v_norm, normal)), -1, 1)))
    offset_angle = 90 - angle_from_normal  # 0 = directly above ring face

    if offset_angle <= angle_max:
        bond_properties['exists'] = True
        bond_properties['info'].append({
            'distance': round(float(d), 2),
            'offset_angle': round(float(offset_angle), 1),
            'nz_status': nz['status'],
        })
    return bond_properties