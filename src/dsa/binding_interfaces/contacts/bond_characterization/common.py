import gemmi
import numpy as np

def get_atom(residue, atom_name):
    """Return highest-occupancy atom matching atom_name, or None."""
    candidates = [a for a in residue if a.name == atom_name]
    if not candidates:
        return None
    return max(candidates, key=lambda a: a.occ)

def atom_pos(atom):
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z])

def get_residue_atoms(residue):
    """Dict of atom_name -> gemmi.Atom, deduplicated by occupancy."""
    out = {}
    for a in residue:
        if a.name not in out or a.occ > out[a.name].occ:
            out[a.name] = a
    return out

def nerf_place(a, b, c, bond_length, bond_angle_deg, dihedral_deg):
    bond_angle = np.radians(bond_angle_deg)
    dihedral = np.radians(dihedral_deg)
    bc = c - b
    bc_norm = bc / np.linalg.norm(bc)
    ab = b - a
    n = np.cross(ab, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc_norm)
    d2 = -bond_length * np.cos(bond_angle)
    d1 = bond_length * np.sin(bond_angle) * np.cos(dihedral)
    d0 = bond_length * np.sin(bond_angle) * np.sin(dihedral)
    return c + d2 * bc_norm + d1 * m + d0 * n


def get_lys_nz(lys_residue, default_chi4=180.0):
    """
    Returns dict: {'coord': np.array or None, 'status': str, 'atom': gemmi.Atom or None}
    status in: 'observed', 'reconstructed', 'unresolvable'
    """
    atoms = get_residue_atoms(lys_residue)

    if 'NZ' in atoms:
        return {'coord': atom_pos(atoms['NZ']), 'status': 'observed', 'atom': atoms['NZ']}

    if all(n in atoms for n in ('CG', 'CD', 'CE')):
        nz = nerf_place(
            atom_pos(atoms['CG']), atom_pos(atoms['CD']), atom_pos(atoms['CE']),
            bond_length=1.47, bond_angle_deg=111.0, dihedral_deg=default_chi4
        )
        return {'coord': nz, 'status': 'reconstructed', 'atom': None}

    return {'coord': None, 'status': 'unresolvable', 'atom': None}