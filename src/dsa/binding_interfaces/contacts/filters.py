
from dsa.binding_interfaces.structure.structure import Structure
from abc import ABC, abstractmethod
import gemmi
from dsa.binding_interfaces.contacts.contact import Contact
from dsa.binding_interfaces.factory import ClassFactory

class ContactFilter(ABC):

    @abstractmethod
    def __call__(self, structure: Structure, contact: gemmi.ContactSearch.Result) -> bool:
        pass
    
    def __str__(self):
        return f"{self.__class__.__name__}"

@ClassFactory.register(interface=ContactFilter, key="lysine_only")
class LysineOnlyFilter(ContactFilter):
    def __init__(self):
        pass

    def __call__(self, contact: Contact) -> bool:
        if not contact.has_lysine_and_nz():
            return False
        if not contact.are_both_polymers():
            return False
        if contact.are_adjacent():
            return False
        return True

#lysine_only_filter = LysineOnlyFilter()
# @classmethod
#     def is_representative(cls, res: gemmi.Residue, atom: gemmi.Atom, polymer_type: gemmi.PolymerType) -> bool:
#         if res.entity_type != gemmi.EntityType.Polymer:
#             return False
#         #entity = next((e for e in self.structure.entities if e.name == res.entity_id), None)
#         #polymer_type = entity.polymer_type if entity is not None else None
#         if polymer_type == gemmi.PolymerType.PeptideL or polymer_type == gemmi.PolymerType.PeptideD:
#             if res.name == 'GLY':
#                 return atom.name == 'CA'
#             else:
#                 return atom.name == 'CB'
#         if polymer_type == gemmi.PolymerType.Dna or polymer_type == gemmi.PolymerType.Rna or polymer_type == gemmi.PolymerType.DnaRnaHybrid:
#             return atom.name == "C1'"
#         return False


#     @classmethod
#     def is_adjacent(cls, chain1: str, chain2: str, res1: gemmi.Residue, res2: gemmi.Residue) -> bool:
#         if chain1 != chain2:
#             return False
#         delta_seqnum = abs(res1.seqid.num - res2.seqid.num)
#         if delta_seqnum <= 1:
#             return True
#         return False



#     def is_contact_pair(self, contact: gemmi.Contact) -> bool:
#         count = 0
#         iface = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
#         for r in self.pairwise_contacts:
#             cra1, cra2 = r.partner1, r.partner2
#             chain1, chain2 = cra1.chain.name, cra2.chain.name
#             res1, res2 = cra1.residue, cra2.residue
#             atom1, atom2 = cra1.atom, cra2.atom
#             author_chain1 = self.bioChain_to_author[cra1.chain.name]
#             author_chain2 = self.bioChain_to_author[cra2.chain.name]
#             polymer_type1 = get_polymer_type(res1, structure)
#             polymer_type2 = get_polymer_type(res2, structure)
#             if res1.entity_type != gemmi.EntityType.Polymer or res2.entity_type != gemmi.EntityType.Polymer:
#                 continue
#             if res1.entity_type == gemmi.EntityType.Water or res2.entity_type == gemmi.EntityType.Water:
#                 continue
#             if self.is_adjacent(chain1, chain2, res1, res2):
#                 continue
#             if self.lysine_only:
#                 is_zeta_nitrogen1 = atom1.name == 'NZ' and res1.name == 'LYS'
#                 is_zeta_nitrogen2 = atom2.name == 'NZ' and res2.name == 'LYS'
#                 if is_zeta_nitrogen1:
#                     iface[author_chain1][res1][author_chain2].append(res2)
#                 if is_zeta_nitrogen2:
#                     iface[author_chain2][res2][author_chain1].append(res1)
#                 if is_zeta_nitrogen1 or is_zeta_nitrogen2:
#                     count += 1
#             elif self.is_representative(res1, atom1) and self.is_representative(res2, atom2):
#                 iface[author_chain1][res1][author_chain2].append(res2)
#                 iface[author_chain2][res2][author_chain1].append(res1)
#                 count += 1
#         print(f"Number of biologically relevant contact pairs: {count}")
#         return iface

#     @classmethod
#     def ContactFactory(cls, contact: gemmi.Contact) -> Contact: