from dsa.binding_interfaces.contacts.bond_characterization.saltBridge import check_salt_bridge
from dsa.binding_interfaces.contacts.bond_characterization.hydrogenBond import check_hydrogen_bonds
from dsa.binding_interfaces.contacts.bond_characterization.cationPi import check_cation_pi
from dsa.binding_interfaces.contacts.bond_characterization.covalentIsopeptide import check_covalent_isopeptide
from dsa.binding_interfaces.contacts.helpers import characterize_residue
from collections import defaultdict
import gemmi
from dsa.binding_interfaces.contacts.bond_characterization.library import BondLibrary
from dsa.misc.utilities import defaultdict_to_dict
# Current implelementation assumes one of the residues is a lysine
class BondCharacteristics:
    def __init__(self, residue1: gemmi.Residue, residue2: gemmi.Residue):
        self.bond_characteristics = defaultdict(lambda: defaultdict(dict))
        for bond_type in BondLibrary.bond_types_implemented():
            for method in BondLibrary.methods_implemented(bond_type):
                characteristics = BondLibrary.check_bond(residue1, residue2, bond_type, method)
                self.bond_characteristics[bond_type][method]['exists'] = characteristics['exists']
                self.bond_characteristics[bond_type][method]['info'] = characteristics['info']
        self.bond_characteristics = defaultdict_to_dict(self.bond_characteristics)
        
    def satisfies_bond_type(self, bond_type: str, method: str|None=None) -> bool:
        if bond_type not in self.bond_characteristics:
            raise ValueError(f"Bond type {bond_type} not implemented")
        if method is None:  
            return any(char['exists'] for char in self.bond_characteristics[bond_type].values())
        elif method in self.bond_characteristics[bond_type]:
            return self.bond_characteristics[bond_type][method]['exists']
        else:
            raise ValueError(f"Method {method} not implemented for bond type {bond_type}")
        

# def characterize_lysine_bond(residue1, residue2):
#     """Runs all checks, returns flat list of interaction dicts (possibly empty)."""
#     properties1 = characterize_residue(residue1)
#     properties2 = characterize_residue(residue2)
#     bond_properties2 = {}
#     bond_properties1 = {}
#     if properties1['is_lysine']:
#         bond_properties1 = _characterize_lysine_bond(residue1, residue2)
#     elif properties2['is_lysine']:
#         bond_properties2 = _characterize_lysine_bond(residue2, residue1)

#     bond_properties = {}
#     for bp_key in bond_properties1.keys() | bond_properties2.keys():
#         bond_properties[bp_key] = {
#             'exists': False,
#             'residue1_as_lysine': {'exists': False},
#             'residue2_as_lysine': {'exists': False}
#         }
#         if bond_properties1 and bond_properties2:
#             bond_properties[bp_key]['exists'] = bond_properties1[bp_key]['exists'] | bond_properties2[bp_key]['exists']
#         elif bond_properties1:
#             bond_properties[bp_key]['exists'] = bond_properties1[bp_key]['exists']
#         elif bond_properties2:
#             bond_properties[bp_key]['exists'] = bond_properties2[bp_key]['exists']
#         if bond_properties1:
#             bond_properties[bp_key]['residue1_as_lysine'] = bond_properties1[bp_key]
#         if bond_properties2:
#             bond_properties[bp_key]['residue2_as_lysine'] = bond_properties2[bp_key]
#     return bond_properties


# def _characterize_lysine_bond(lys_residue, partner_residue):
#     bond_properties = {}
#     bond_properties['salt_bridge'] = check_salt_bridge(lys_residue, partner_residue)
#     bond_properties['hydrogen_bond'] = check_hydrogen_bonds(lys_residue, partner_residue)
#     bond_properties['cation_pi'] = check_cation_pi(lys_residue, partner_residue)
#     bond_properties['covalent_isopeptide'] = check_covalent_isopeptide(lys_residue, partner_residue)
#     return bond_properties