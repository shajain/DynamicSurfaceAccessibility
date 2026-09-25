"""
Find biologically relevant inter-chain contacts and extract
interface residue positions for each chain pair.

Output: dict mapping (chain_A, chain_B) -> {
    'chain_a_residues': list of (residue_name, residue_seqnum, CA_position),
    'chain_b_residues': list of (residue_name, residue_seqnum, CA_position),
    'contacts':         list of (res_a, res_b, distance)
}
"""

import gemmi
from collections import defaultdict
import re
import numpy as np

LYSINE_NAMES = {'LYS', 'LYZ', 'KCX', 'MLY', 'MLZ', 'M3L', 'ALY'}

class BindingInterfaceExtractor:
    def __init__(self, structure: gemmi.Structure, cutoff: float = 8.0, lysine_only: bool = False, assembly_index: int = 0):
        self.structure = structure
        self.cutoff = cutoff
        self.lysine_only = lysine_only
        self.assembly_index = assembly_index
        self.asym_to_author = self.asym_to_author_chain_mapping()
        self.bio_st = self.build_bio_structure()
        self.bioChain_to_author = self.build_bioChain_to_authorChain_mapping()
        self.pairwise_contacts = self.search_pairwise_contacts()
        self.interfaces = self.process_pairwise_contacts()


    # def get_ca_position(residue):
    #     """
    #     Return the CA (alpha carbon) position of a residue.
    #     Falls back to the first atom if CA is not found (e.g. for non-standard residues).
    #     """
    #     for atom in residue:
    #         if atom.name == 'CA':
    #             return atom.pos
    #     # fallback to first atom
    #     if len(residue) > 0:
    #         return residue[0].pos
    #     return None
    
    def build_bioChain_to_authorChain_mapping(self) -> str:
        bioChain_to_author = {}
        for chain in self.bio_st:
            for residue in chain:
                if residue.subchain and residue.entity_type == gemmi.EntityType.Polymer:
                    asym_id = re.sub(r'\d+', '', residue.subchain)
                    bioChain_to_author[chain.name] = self.asym_to_author[asym_id]
                    break
        return bioChain_to_author

    # def asym_to_author_chain_mapping(self):
    #     scheme = self.block.find('_pdbx_poly_seq_scheme.',
    #                 ['asym_id', 'pdb_strand_id'])
    #     asym_to_author = {row[0]: row[1] for row in scheme}
    #     return asym_to_author
    
    def asym_to_author_chain_mapping(self):
        asym_to_author = {}
        for chain in self.structure[0]:
            for residue in chain:
                if residue.subchain:
                    asym_to_author[residue.subchain] = chain.name
                    break
        return asym_to_author

    

    def build_bio_structure(self):
        # ------------------------------------------------------------------
        # 1. Check biological assemblies are present
        # ------------------------------------------------------------------
        # if not self.structure.assemblies:
        #     raise ValueError(
        #         "No biological assembly annotations found in this file. "
        #         "Cannot distinguish biological contacts from crystal packing."
        #     )

        if not self.structure.assemblies:
            print("No biological assembly annotations found, using structure as-is.")
            return self.structure[0]

        assembly = self.structure.assemblies[self.assembly_index]
        print(f"Using biological assembly: '{assembly.name}'")
        # ------------------------------------------------------------------
        # 2. Build the biological assembly explicitly
        #    This applies all rotation/translation operators stored in the file
        #    to reconstruct the full functional complex.
        #    AddNumber renames copied chains: A -> A1, A2 etc.
        # ------------------------------------------------------------------
        bio_st = gemmi.make_assembly(assembly, self.structure[0], gemmi.HowToNameCopiedChain.AddNumber)
        return bio_st

    def search_pairwise_contacts(self):
        # ------------------------------------------------------------------
        # 4. Run NeighborSearch on the biological assembly
        #    We use bio_st.cell but if the biological assembly has no unit cell
        #    (common after make_assembly), gemmi uses the bounding box instead.
        # ------------------------------------------------------------------
        ns = gemmi.NeighborSearch(self.bio_st, self.structure.cell, self.cutoff).populate(include_h=False)
        # ------------------------------------------------------------------
        # 5. Run ContactSearch — ignore same-chain contacts
        # ------------------------------------------------------------------
        cs = gemmi.ContactSearch(self.cutoff)
        #cs.ignore = gemmi.ContactSearch.Ignore.AdjacentResidues
        cs.ignore = gemmi.ContactSearch.Ignore.SameResidue
        results = cs.find_contacts(ns)
        #print(f"Total inter-chain contacts found: {len(results)}")
        # ------------------------------------------------------------------
        # 6. Filter out symmetry copies (image_idx > 0)
        #    After make_assembly, all biologically relevant chains are explicit.
        #    Any image_idx > 0 hit is a crystal packing artifact on top of that.
        # ------------------------------------------------------------------
        results = [r for r in results if r.image_idx == 0]
        # print(f"Number of biologically relevant contacts: {len(results)}")
        return results


    def is_representative(self, res: gemmi.Residue, atom: gemmi.Atom) -> bool:
        if res.entity_type != gemmi.EntityType.Polymer:
            return False
        entity = next((e for e in self.structure.entities if e.name == res.entity_id), None)
        polymer_type = entity.polymer_type if entity is not None else None
        if polymer_type == gemmi.PolymerType.PeptideL or polymer_type == gemmi.PolymerType.PeptideD:
            if res.name == 'GLY':
                return atom.name == 'CA'
            else:
                return atom.name == 'CB'
        if polymer_type == gemmi.PolymerType.Dna or polymer_type == gemmi.PolymerType.Rna or polymer_type == gemmi.PolymerType.DnaRnaHybrid:
            return atom.name == "C1'"
        return False
    def is_adjacent(self, chain1: str, chain2: str, res1: gemmi.Residue, res2: gemmi.Residue) -> bool:
        if chain1 != chain2:
            return False
        delta_seqnum = abs(res1.seqid.num - res2.seqid.num)
        if delta_seqnum <= 1:
            return True
        return False

    def process_pairwise_contacts(self) -> dict[str, dict[str, dict[str, list[gemmi.Residue]]]]:
        count = 0
        iface = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        for r in self.pairwise_contacts:
            cra1, cra2 = r.partner1, r.partner2
            chain1, chain2 = cra1.chain.name, cra2.chain.name
            res1, res2 = cra1.residue, cra2.residue
            atom1, atom2 = cra1.atom, cra2.atom
            author_chain1 = self.bioChain_to_author[cra1.chain.name]
            author_chain2 = self.bioChain_to_author[cra2.chain.name]
            if res1.entity_type != gemmi.EntityType.Polymer or res2.entity_type != gemmi.EntityType.Polymer:
                continue
            if res1.entity_type == gemmi.EntityType.Water or res2.entity_type == gemmi.EntityType.Water:
                continue
            if self.is_adjacent(chain1, chain2, res1, res2):
                continue
            if self.lysine_only:
                is_zeta_nitrogen1 = atom1.name == 'NZ' and res1.name == 'LYS'
                is_zeta_nitrogen2 = atom2.name == 'NZ' and res2.name == 'LYS'
                if is_zeta_nitrogen1:
                    iface[author_chain1][res1][author_chain2].append(res2)
                if is_zeta_nitrogen2:
                    iface[author_chain2][res2][author_chain1].append(res1)
                if is_zeta_nitrogen1 or is_zeta_nitrogen2:
                    count += 1
            elif self.is_representative(res1, atom1) and self.is_representative(res2, atom2):
                iface[author_chain1][res1][author_chain2].append(res2)
                iface[author_chain2][res2][author_chain1].append(res1)
                count += 1
        print(f"Number of biologically relevant contact pairs: {count}")
        return iface

    def group_residues_by_contact_type(self, chain: str) -> dict[str, list[gemmi.Residue]]:
        keys = ["all", "residue-residue", "residue-other", "intra-chain", "inter-chain"]
        contact_residues = defaultdict(list, {k: [] for k in keys})
        for r1 in self.interfaces[chain].keys():
            if r1.entity_type == gemmi.EntityType.Polymer:
                contact_residues["all"].append(r1)
                is_res_res, is_res_other, is_intra_chain, is_inter_chain = False, False, False, False
                for c2 in self.interfaces[chain][r1]:
                    for r2 in self.interfaces[chain][r1][c2]:
                        if r2.entity_type == gemmi.EntityType.Polymer:
                            if not is_res_res:
                                contact_residues["residue-residue"].append(r1)
                                is_res_res = True
                            if c2 == chain and not is_intra_chain:
                                contact_residues["intra-chain"].append(r1)
                                is_intra_chain = True 
                            elif c2 != chain and not is_inter_chain:
                                contact_residues["inter-chain"].append(r1)
                                is_inter_chain = True
                        elif r2.entity_type != gemmi.EntityType.Water and not is_res_other:
                            contact_residues["residue-other"].append(r1)    
                            is_res_other = True
                        if is_res_res and is_res_other and is_intra_chain and is_inter_chain:
                            break
                    if is_res_res and is_res_other and is_intra_chain and is_inter_chain:
                        break
        contact_residues = {k: list(set(v)) for k, v in contact_residues.items()}
        return contact_residues

    def custom_contact_extraction(self):
        return self._custom_contact_extraction(self.bio_st, self.cutoff)
    
    def chainwise_contacts(self, chain1: str, chain2: str=None):
        if chain2 is None:
            contacts = self.interfaces[chain1].keys()
            return contacts
        else:
            contact_pairs = [(r1, r2) for r1, c2 in self.interfaces[chain1].items() 
                                            if chain2 in c2 
                                            for r2 in c2[chain2]]
            return list(set(contact_pairs))

    @staticmethod
    def _custom_contact_extraction(structure: gemmi.Structure, cutoff: float):
        contacts = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        lysine_residues_nz = BindingInterfaceExtractor.lysine_residues_nz(structure)
        polymer_residues = BindingInterfaceExtractor.polymer_residues(structure)
        for chain in polymer_residues:
            for residue in polymer_residues[chain]:
                for atom in residue:
                    if atom.element.atomic_number > 1:
                        for chain2 in lysine_residues_nz:
                            if chain2 == chain:
                                continue
                            for lys_res, nz_atom in lysine_residues_nz[chain2]:
                                nz_pos = nz_atom.pos
                                nz_pos_arr = np.array([nz_pos.x, nz_pos.y, nz_pos.z])
                                atom_pos_arr = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                                dist = np.linalg.norm(atom_pos_arr - nz_pos_arr)
                                if dist < cutoff:
                                    contacts[chain2][lys_res][chain].append((residue, atom))
        return contacts

    @staticmethod
    def lysine_residues(structure: gemmi.Structure) -> list[gemmi.Residue]:
        lysine_residues = defaultdict(list)
        for chain in structure:
            for residue in chain:
                if residue.name in LYSINE_NAMES:
                    lysine_residues[chain.name].append(residue)
        return lysine_residues
    @staticmethod
    def polymer_residues(structure: gemmi.Structure) -> list[gemmi.Residue]:
        polymer_residues = defaultdict(list)
        for chain in structure:
            for residue in chain:
                if residue.entity_type == gemmi.EntityType.Polymer:
                    polymer_residues[chain.name].append(residue)
        return polymer_residues
    @staticmethod
    def lysine_residues_nz(structure: gemmi.Structure) -> dict[gemmi.Residue, gemmi.Atom]:
        lysine_residues = BindingInterfaceExtractor.lysine_residues(structure)
        lysine_residues_nz = defaultdict(list)
        for chain in lysine_residues:
            for residue in lysine_residues[chain]:
                for atom in residue:
                    if atom.name == 'NZ':
                        lysine_residues_nz[chain].append((residue, atom))
                        break
        return lysine_residues_nz
        
    # def print_interfaces(interfaces):
    #     """Pretty-print the interface results."""
    #     for (chain_a, chain_b), data in interfaces.items():
    #         orig_a, orig_b = data['original_chains']
    #         print(f"\n{'='*60}")
    #         print(f"Interface: chain {chain_a} (orig: {orig_a})  <->  "
    #             f"chain {chain_b} (orig: {orig_b})")
    #         print(f"  {len(data['chain_a_residues'])} interface residues in chain {chain_a}")
    #         print(f"  {len(data['chain_b_residues'])} interface residues in chain {chain_b}")
    #         print(f"  {len(data['contacts'])} contacts total")

    #         print(f"\n  Chain {chain_a} interface residues (name, seqnum, CA position):")
    #         for res_name, seqnum, pos in data['chain_a_residues']:
    #             print(f"    {res_name:4s} {seqnum:5d}   "
    #                 f"({pos.x:7.2f}, {pos.y:7.2f}, {pos.z:7.2f})")

    #         print(f"\n  Chain {chain_b} interface residues (name, seqnum, CA position):")
    #         for res_name, seqnum, pos in data['chain_b_residues']:
    #             print(f"    {res_name:4s} {seqnum:5d}   "
    #                 f"({pos.x:7.2f}, {pos.y:7.2f}, {pos.z:7.2f})")

    #         print(f"\n  Closest contacts:")
    #         for label_a, label_b, dist in data['contacts'][:10]:  # show top 10
    #             print(f"    {label_a:20s} <-> {label_b:20s}  {dist:.2f} A")


# ------------------------------------------------------------------
# Example usage
# ------------------------------------------------------------------
# if __name__ == "__main__":
#     import sys

#     pdb_path = sys.argv[1] if len(sys.argv) > 1 else "your_structure.cif"
#     cutoff   = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0

#     interfaces = find_interface_residues(pdb_path, cutoff=cutoff)
#     print_interfaces(interfaces)



  # def find_interfaces(self):
    #     """
    #     Find biologically relevant inter-chain contacts.
    #     Returns
    #     -------
    #     interfaces : dict
    #         Keys   -> (chain_name_a, chain_name_b) tuple, always sorted alphabetically
    #         Values -> {
    #             'chain_a_residues': list of (res_name, seq_num, gemmi.Position),
    #             'chain_b_residues': list of (res_name, seq_num, gemmi.Position),
    #             'contacts':         list of (res_a_label, res_b_label, dist_angstrom)
    #         }
    #     """
    #     # ------------------------------------------------------------------
    #     # 2. Check biological assemblies are present
    #     # ------------------------------------------------------------------
    #     if not self.structure.assemblies:
    #         raise ValueError(
    #             "No biological assembly annotations found in this file. "
    #             "Cannot distinguish biological contacts from crystal packing."
    #         )

    #     assembly = self.structure.assemblies[self.assembly_index]
    #     print(f"Using biological assembly: '{assembly.name}'")

    #     # ------------------------------------------------------------------
    #     # 3. Build the biological assembly explicitly
    #     #    This applies all rotation/translation operators stored in the file
    #     #    to reconstruct the full functional complex.
    #     #    AddNumber renames copied chains: A -> A1, A2 etc.
    #     # ------------------------------------------------------------------
    #     bio_st = gemmi.make_assembly(assembly, self.structure[0], gemmi.HowToNameCopiedChain.AddNumber)

    #     bioChain_to_author = self.build_bioChain_to_authorChain_mapping(bio_st)

    #     # ------------------------------------------------------------------
    #     # 4. Run NeighborSearch on the biological assembly
    #     #    We use bio_st.cell but if the biological assembly has no unit cell
    #     #    (common after make_assembly), gemmi uses the bounding box instead.
    #     # ------------------------------------------------------------------
    #     ns = gemmi.NeighborSearch(bio_st, self.structure.cell, self.cutoff).populate(include_h=False)

    #     # ------------------------------------------------------------------
    #     # 5. Run ContactSearch — ignore same-chain contacts
    #     # ------------------------------------------------------------------
    #     cs = gemmi.ContactSearch(self.cutoff)
    #     cs.ignore = gemmi.ContactSearch.Ignore.AdjacentResidues
    #     results = cs.find_contacts(ns)
    #     #print(f"Total inter-chain contacts found: {len(results)}")
    #     # ------------------------------------------------------------------
    #     # 6. Filter out symmetry copies (image_idx > 0)
    #     #    After make_assembly, all biologically relevant chains are explicit.
    #     #    Any image_idx > 0 hit is a crystal packing artifact on top of that.
    #     # ------------------------------------------------------------------
    #     bio_results = [r for r in results if r.image_idx == 0]
    #     print(f"Number of biologically relevant inter-chain contacts: {len(bio_results)}")

    #     # ------------------------------------------------------------------
    #     # 7. Collect interface residues per chain
    #     # ------------------------------------------------------------------

    #     iface = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    #     for r in bio_results:
    #         cra1 = r.partner1   # CRA = Chain, Residue, Atom
    #         cra2 = r.partner2

    #         chain1 = bioChain_to_author[cra1.chain.name]
    #         chain2 = bioChain_to_author[cra2.chain.name]
    #         res1   = cra1.residue
    #         res2   = cra2.residue

    #         if res1.entity_type != gemmi.EntityType.Polymer and res2.entity_type != gemmi.EntityType.Polymer:
    #             continue
    #         if res1.entity_type == gemmi.EntityType.Water or res2.entity_type == gemmi.EntityType.Water:
    #             continue
    #         if res1 == res2:
    #             continue
    #         iface[chain1][res1][chain2].append(res2)
    #         iface[chain2][res2][chain1].append(res1)
    #     return iface