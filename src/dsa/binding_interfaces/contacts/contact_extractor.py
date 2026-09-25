import gemmi
from collections import defaultdict
from dsa.binding_interfaces.structure.structure import Structure
# from dsa.binding_interfaces.structure_uniprot_alignment.aligner import StructureSequenceAligner
from dsa.binding_interfaces.contacts.filters import ContactFilter
from dsa.binding_interfaces.contacts.contact import Contact
from dsa.binding_interfaces.factory import ClassFactory
import re

from dsa.misc.file_strings import FileStrings

class ContactExtractor:

    def __init__(self, distance_cutoff: float, contact_filter_key: str):
        self.distance_cutoff = distance_cutoff
        self.filter = ClassFactory.create(interface=ContactFilter, key=contact_filter_key)
    
    

    # def build_bio_structure(self):
    #     # ------------------------------------------------------------------
    #     # 1. Check biological assemblies are present
    #     # ------------------------------------------------------------------
    #     # if not self.structure.assemblies:
    #     #     raise ValueError(
    #     #         "No biological assembly annotations found in this file. "
    #     #         "Cannot distinguish biological contacts from crystal packing."
    #     #     )

    #     if not self.structure.assemblies:
    #         print("No biological assembly annotations found, using structure as-is.")
    #         return self.structure[0]

    #     assembly = self.structure.assemblies[self.assembly_index]
    #     print(f"Using biological assembly: '{assembly.name}'")
    #     # ------------------------------------------------------------------
    #     # 2. Build the biological assembly explicitly
    #     #    This applies all rotation/translation operators stored in the file
    #     #    to reconstruct the full functional complex.
    #     #    AddNumber renames copied chains: A -> A1, A2 etc.
    #     # ------------------------------------------------------------------
    #     bio_st = gemmi.make_assembly(assembly, self.structure[0], gemmi.HowToNameCopiedChain.AddNumber)
    #     return bio_st


    def extract(self, structure: Structure):
        # structure_uniprot_aligner = ClassFactory.create(StructureSequenceAligner, source_type=structure.sour, structure_id=structure.structure_id)
        contacts = [Contact(c, structure) for c in self.distance_search(structure)]
        filtered_contacts = [c for c in contacts if self.filter(c)]
        return filtered_contacts

    def distance_search(self, structure: Structure):
        # ------------------------------------------------------------------
        # 4. Run NeighborSearch on the biological assembly
        #    We use bio_st.cell but if the biological assembly has no unit cell
        #    (common after make_assembly), gemmi uses the bounding box instead.
        # ------------------------------------------------------------------
        bio_model = structure.gemmi_bio_model
        gemmi_structure = structure.gemmi_structure
        ns = gemmi.NeighborSearch(bio_model, gemmi_structure.cell, self.distance_cutoff).populate(include_h=False)
        # ------------------------------------------------------------------
        # 5. Run ContactSearch — ignore same-chain contacts
        # ------------------------------------------------------------------
        cs = gemmi.ContactSearch(self.distance_cutoff)
        cs.ignore = gemmi.ContactSearch.Ignore.AdjacentResidues
        #cs.ignore = gemmi.ContactSearch.Ignore.SameResidue
        contacts = cs.find_contacts(ns)
        #print(f"Total inter-chain contacts found: {len(results)}")
        # ------------------------------------------------------------------
        # 6. Filter out symmetry copies (image_idx > 0)
        #    After make_assembly, all biologically relevant chains are explicit.
        #    Any image_idx > 0 hit is a crystal packing artifact on top of that.
        # ------------------------------------------------------------------
        contacts = [c for c in contacts if c.image_idx == 0]
        # print(f"Number of biologically relevant contacts: {len(results)}")
        return contacts


    def __str__(self):
        cutoff_string = FileStrings.stringify_cutoff(self.distance_cutoff)
        string = f"{cutoff_string}-{self.filter}"
        return string


    


    