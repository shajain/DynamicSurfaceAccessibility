from dataclasses import dataclass
import traceback
import gemmi
from dsa.binding_interfaces.structure.structure import Structure
from dsa.binding_interfaces.contacts.helpers import ALL_LYSINE_NAMES, is_lysine, is_nz_atom
from dsa.binding_interfaces.contacts.helpers import is_atom_representative, is_peptide, is_dna_or_rna
from dsa.binding_interfaces.contacts.helpers import are_residues_adjacent
#from dsa.binding_interfaces.structure_uniprot_alignment.aligner import StructureSequenceAligner
from dsa.binding_interfaces.contacts.helpers import characterize_residue
#from dsa.binding_interfaces.contacts.bond_characterization.characterize import characterize_lysine_bond
from dsa.binding_interfaces.contacts.bond_characterization.bond_characteristics import BondCharacteristics
from dsa.binding_interfaces.common_functions import single_letter_residue_code

class Contact:
    def __init__(self, contact: gemmi.ContactSearch.Result, 
                       structure: Structure,
                       #structure_uniprot_aligner: StructureSequenceAligner
                       ):
        partner1, partner2 = contact.partner1, contact.partner2
        self.structure = structure
        self.chain1 = partner1.chain
        self.residue1 = partner1.residue
        self.atom1 = partner1.atom
        self.chain2 = partner2.chain
        self.residue2 = partner2.residue
        self.atom2 = partner2.atom
        self.dist = contact.dist
        self.properties1 = characterize_residue(self.residue1)
        self.properties2 = characterize_residue(self.residue2)
        # self.bond_properties = {}
        # self.is_bond = False
        # self.bond_type = []
        #if self.has_lysine():
            # self.bond_properties = characterize_lysine_bond(self.residue1, self.residue2)
            # self.is_bond = any([bp_dict['exists'] for bp_dict in self.bond_properties.values()])
            # self.bond_type = [bp_key for (bp_key, bp_dict) in self.bond_properties.items() if bp_dict['exists']]
        self.bond_properties = BondCharacteristics(self.residue1, self.residue2)
        try:
            self.author_chain1 = structure.biochain_to_chain(self.chain1.name)
            self.author_chain2 = structure.biochain_to_chain(self.chain2.name)
        except KeyError:
            raise ValueError(f"Either Chain {self.chain1.name} or {self.chain2.name} in structure {structure.structure_id} not mapped to author chain")
        # self.entity1 = self.residue1.entity_id
        # self.entity2 = self.residue2.entity_id
        # self.uniprot1 = structure.entity_to_uniprot[self.entity1]
        # self.uniprot2 = structure.entity_to_uniprot[self.entity2]
        # self.unp_residue1 = structure_uniprot_aligner.map(self.author_chain1, self.uniprot1, self.residue1) if self.uniprot1 else None
        # self.unp_residue2 = structure_uniprot_aligner.map(self.author_chain2, self.uniprot2, self.residue2) if self.uniprot2 else None

    
    def has_amino_acid(self):
        return self.properties1['is_amino_acid'] or self.properties2['is_amino_acid']
    
    def both_are_amino_acids(self) -> bool:
        return self.properties1['is_amino_acid'] and self.properties2['is_amino_acid']

    def are_both_polymers(self) -> bool:
        return self.properties1['is_polymer'] and self.properties2['is_polymer']

    def has_lysine_and_nz(self) -> bool:
        return self.is_residue1_lysine_and_has_nz() or self.is_residue2_lysine_and_has_nz()


    def has_lysine(self) -> bool:
        return self.is_residue1_lysine() or self.is_residue2_lysine()
    
    def is_residue1_lysine(self) -> bool:
        return self.properties1['is_lysine']
    
    def is_residue2_lysine(self) -> bool:
        return self.properties2['is_lysine']

    def is_residue1_lysine_and_has_nz(self) -> bool:
        return self.properties1['is_lysine'] and self.atom1.name == 'NZ'
    
    def is_residue2_lysine_and_has_nz(self) -> bool:
        return self.properties2['is_lysine'] and self.atom2.name == 'NZ'

    # def has_lysine_and_nz(self) -> bool:
    #     return self.is_residue1_lysine_and_has_nz() or self.is_residue2_lysine_and_has_nz()
    
    def are_both_atoms_representative(self) -> bool:
        is_atom1 = is_atom_representative(self.residue1, self.atom1, self.residue_type1)
        is_atom2 = is_atom_representative(self.residue2, self.atom2, self.residue_type2)
        return is_atom1 and is_atom2

    def residue1_in_uniprot_id(self, uniprot_id:str) -> tuple[int, str] | None:
        return self.residue_in_uniprot_id(self.author_chain1, self.residue1, uniprot_id)
    
    def residue2_in_uniprot_id(self, uniprot_id:str) -> tuple[int, str] | None:
        return self.residue_in_uniprot_id(self.author_chain2, self.residue2, uniprot_id)

    def residue_in_uniprot_id(self, chain:str, residue:gemmi.Residue|tuple[str, str], uniprot_id:str) -> tuple[int, str] | None:
        res = self.structure.map_residue_to_uniprot(chain, residue, uniprot_id)
        return res


    def is_distance_within(self, distance_threshold: float) -> bool:
        return self.dist <= distance_threshold

    def satisfies_contact_type(self, contact_type: str) -> bool:
        if contact_type == "intra_chain":
            return self.chain1 == self.chain2
        elif contact_type == "inter_chain":
            if self.chain1 != self.chain2 and self.author_chain1 == self.author_chain2:
                print(f"found a case where biochains are different but author chains are the same: this is expected behavior.")
            return self.chain1 != self.chain2
        else:
            raise ValueError(f"Invalid contact type: {contact_type}")

    def satisfies_bond_type(self, bond_type: str, method: str|None=None) -> bool:
        return self.bond_properties.satisfies_bond_type(bond_type, method)

    # def satisfies_bond_type1(self, bond_type: str) -> bool:
    #     if bond_type not in self.bond_properties:
    #         ValueError(f"Invalid bond type: {bond_type}")
    #     bond_info = self.bond_properties[bond_type]['residue1_as_lysine']
    #     return 'exists' in bond_info and bond_info['exists']

    # def satisfies_bond_type2(self, bond_type: str) -> bool:
    #     if bond_type not in self.bond_properties:
    #         ValueError(f"Invalid bond type: {bond_type}")
    #     bond_info = self.bond_properties[bond_type]['residue2_as_lysine']
    #     return 'exists' in bond_info and bond_info['exists']

    def is_inter_chain(self) -> bool:
        return self.chain1 != self.chain2

    # def bond_types(self) -> list[str]:
    #     return self.bond_properties.bond_types_satisfied()
    
    # def residue1_bond_types(self) -> list[str]:
    #     return [bp_key for bp_key in self.bond_properties.keys() if self.bond_properties[bp_key]['residue1_as_lysine']['exists']]
    
    # def residue2_bond_types(self) -> list[str]:
    #     [bp_key for bp_key in self.bond_properties.keys() if self.bond_properties[bp_key]['residue2_as_lysine']['exists']]

    # def is_atom1_nz(self) -> bool:
    #     return is_nz_atom(self.atom1)
    
    # def is_atom2_nz(self) -> bool:
    #     return is_nz_atom(self.atom2)

    # def is_residue1_lysine(self) -> bool:
    #     return is_lysine(self.residue1)
    
    # def is_residue2_lysine(self) -> bool:
    #     return is_lysine(self.residue2)

    # def is_lysine_and_nz1(self) -> bool:
    #     return is_lysine_and_nz(self.residue1, self.atom1)
    
    # def is_lysine_and_nz2(self) -> bool:
    #     return is_lysine_and_nz(self.residue2, self.atom2)
    
    # def has_lysine_and_nz(self) -> bool:
    #     return self.is_lysine_and_nz1() or self.is_lysine_and_nz2()

    def are_adjacent(self) -> bool:
        if self.chain1 != self.chain2:
            return False
        are_adjacent = are_residues_adjacent(self.residue1, self.residue2)
    
    
    
    # def are_adjacent(self) -> bool:
    #     if self.chain1 != self.chain2:
    #         return False
    #     if self.unp_residue1 is None or self.unp_residue2 is None:
    #         return False
    #     uniprot_pos1 = self.unp_residue1[0]
    #     uniprot_pos2 = self.unp_residue2[0]
    #     if uniprot_pos1 is None or uniprot_pos2 is None:
    #         return False
    #     delta_seqnum = abs(uniprot_pos1 - uniprot_pos2)
    #     if delta_seqnum <= 1:
    #         return True
    #     return False

    def __getstate__(self):
        # print("--- __getstate__ was called! ---")
        # traceback.print_stack()  # This will print the stack trace
        state = self.__dict__.copy()
        state["structure_id"] = state["structure"].structure_id
        state["database"] = state["structure"].database
        del state["structure"]
        state["chain1"] = state["chain1"].name
        state["chain2"] = state["chain2"].name
        state["residue1"] = (str(state["residue1"].seqid), state["residue1"].name, single_letter_residue_code(state["residue1"]))
        state["residue2"] = (str(state["residue2"].seqid), state["residue2"].name, single_letter_residue_code(state["residue2"]))
        state["atom1"] = (state["atom1"].name, state["atom1"].altloc)
        state["atom2"] = (state["atom2"].name, state["atom2"].altloc)
        return state

    # def __setstate__(self, state):
    #     self.__dict__.update(state)
    #     self.chain1, self.residue1, self.atom1 = self.structure.get_bio_chain_residue_atom_objects(state["chain1"], state["residue1"], state["atom1"])
    #     self.chain2, self.residue2, self.atom2 = self.structure.get_bio_chain_residue_atom_objects(state["chain2"], state["residue2"], state["atom2"])

    def __setstate_structure_context__(self, structure: Structure):
        assert self.structure_id == structure.structure_id and self.database == structure.database
        del self.__dict__["structure_id"]
        del self.__dict__["database"]
        self.structure = structure
        self.chain1, self.residue1, self.atom1 = structure.get_bio_chain_residue_atom_objects(self.chain1, self.residue1, self.atom1)
        self.chain2, self.residue2, self.atom2 = structure.get_bio_chain_residue_atom_objects(self.chain2, self.residue2, self.atom2)
        if self.atom1 is None or self.atom2 is None or self.residue1 is None or self.residue2 is None or self.chain1 is None or self.chain2 is None:
            raise ValueError(f"Unable to get CRA objects for {self.chain1}, {self.residue1}, {self.atom1} or {self.chain2}, {self.residue2}, {self.atom2} from {structure.structure_id}")


    # def asdict(self):
    #     return {
    #         "chain1": self.chain1.name,
    #         "residue1": (str(self.res1.seqid), self.res1.name),
    #         "atom1": self.atom1.name,
    #         "chain2": self.chain2.name,
    #         "residue2": (str(self.res2.seqid), self.res2.name),
    #         "atom2": self.atom2.name,
    #         "dist": self.dist,
    #         "uniprot1": self.uniprot1,
    #         "author_chain1": self.author_chain1,
    #         "uniprot_residue1": self.unp_residue1,
    #         "uniprot2": self.uniprot2,
    #         "author_chain2": self.author_chain2,
    #         "uniprot_residue2": self.unp_residue2
    #     }

    # @classmethod
    # def from_gemmi_contact(cls, contact: gemmi.ContactSearch.Result, structure: Structure):
    #     partner1, partner2 = contact.partner1, contact.partner2
    #     contact = cls(partner1.chain, partner1.residue, partner1.atom, 
    #                              partner2.chain, partner2.residue, partner2.atom, 
    #                              contact.dist, structure)
    #     return contact
    
    # @classmethod
    # def from_dict(cls, contact_dict: dict, structure: Structure):
    #     # get chain object from gemmi structure
    #     chain1 = structure.get_bio_chain_object(contact_dict["chain1"])
    #     chain2 = structure.get_bio_chain_object(contact_dict["chain2"])
    #     res1 = structure.get_bio_residue_object(contact_dict["residue1"])
    #     res2 = structure.get_bio_residue_object(contact_dict["residue2"])
    #     atom1 = structure.get_bio_atom_object(res1,contact_dict["atom1"])
    #     atom2 = structure.get_bio_atom_object(res2,contact_dict["atom2"])
    #     dist = contact_dict["dist"]
    #     contact = cls(chain1, res1, atom1, chain2, res2, atom2, dist, structure)
    #     return contact


