from dsa.binding_interfaces.contacts.structure_contacts import StructureContactsStore, StructureContacts



class UniprotAlignedStructureContacts(StructureContacts):
    def __init__(self, uniprot_id: str, structure_id: str, structure_database: str = "pdb", distance_cutoff: float = 6.5, contact_filter_key: str = "lysine_only"):
        self.uniprot_id = uniprot_id
        structure_contacts_store = StructureContactsStore(database=structure_database, 
                                                            distance_cutoff=distance_cutoff, 
                                                            contact_filter_key=contact_filter_key)
        structure_contacts = structure_contacts_store.load(structure_id)
        self.__dict__.update(structure_contacts.__dict__)
        self.contacts, self.unp_residue1s, self.unp_residue2s = self.restrict_to_uniprot_id_and_get_aligned_contacts(uniprot_id)
        
    
    def restrict_to_uniprot_id_and_get_aligned_contacts(self, uniprot_id: str):
        filtered_contacts = []
        unp_residue1s = []
        unp_residue2s = []
        # self.structure.setup_aligner_if_none()
        # chains = self.structure.uniprot_aligner.chains_containing_uniprot_id(uniprot_id)
        # pdbe_chains = self.structure.uniprot_aligner.pdbe_mappings['pdb_chain_id'].unique().tolist()
        # if len(set(chains).intersection(set(pdbe_chains))) > 0:
        #     print("Found contacts in the same chain as the uniprot id")
        for c in self.contacts:
            unp_residue1 = c.residue1_in_uniprot_id(uniprot_id)
            unp_residue2 = c.residue2_in_uniprot_id(uniprot_id)
            if unp_residue1 is not None or unp_residue2 is not None:
                filtered_contacts.append(c)
                unp_residue1s.append(unp_residue1)
                unp_residue2s.append(unp_residue2)
        return filtered_contacts, unp_residue1s, unp_residue2s


    def get_residues_satisfying(self, contact_type: str|None = None, 
                                        bond_type: str|None = None, 
                                        method: str|None = None,
                                        distance_threshold: float|None = None) -> set[str]:
        """
        Returns the UniProt residues that satisfy the given contact type, bond type, method, and distance threshold.
        """
        unp_residues = set()
        for c, res1, res2 in zip(self.contacts, self.unp_residue1s, self.unp_residue2s):
            if not contact_type or c.satisfies_contact_type(contact_type): 
                if not bond_type or c.satisfies_bond_type(bond_type, method):
                    if not distance_threshold or c.is_distance_within(distance_threshold):
                        # if res1 and c.satisfies_bond_type1(bond_type):
                        #     unp_residues.add(res1)   
                        # if res2 and (not bond_type or c.satisfies_bond_type2(bond_type)):
                        #     unp_residues.add(res2)
                        if res1 and c.is_residue1_lysine():
                            unp_residues.add(res1)
                        if res2 and c.is_residue2_lysine():
                            unp_residues.add(res2)
        return unp_residues

    # def group_uniprot_residues_by_bond_types(self):
    #     residues_grouped = {}
    #     residues_grouped['inter-chain']={bp:set() for bp in self.BOND_TYPES}
    #     residues_grouped['intra-chain']={bp:set() for bp in self.BOND_TYPES}
    #     for c, (unp_residue1, unp_residue2) in zip(self.contacts, self.uniprot_aligned_contacts):
    #         inter_or_intra = "intra-chain" if c.is_intra_chain() else "inter-chain"
    #         if unp_residue1 is not None:
    #             for bp in c.residue1_bond_types():
    #                 residues_grouped[inter_or_intra][bp].add(unp_residue1)
    #         if unp_residue2 is not None:
    #             for bp in c.residue2_bond_types():
    #                 residues_grouped[inter_or_intra][bp].add(unp_residue2)
    #     return residues_grouped

    # def group_uniprot_residues_by_distances(self):
    #     residues_grouped = {}
    #     distance_thresholds = np.arange(3.0, self.distance_cutoff, 0.5)
    #     residues_grouped['inter-chain']={dt:set() for dt in distance_thresholds}
    #     residues_grouped['intra-chain']={dt:set() for dt in distance_thresholds}
    #     for c, (unp_residue1, unp_residue2) in zip(self.contacts, self.uniprot_aligned_contacts):
    #         inter_or_intra = "intra-chain" if c.is_intra_chain() else "inter-chain"
    #         if unp_residue1 is not None:
    #             for distance in distance_thresholds:
    #                 if c.dist <= distance:
    #                     residues_grouped[inter_or_intra][distance].add(unp_residue1)
    #         if unp_residue2 is not None:
    #             for distance in distance_thresholds:
    #                 if c.dist <= distance:
    #                     residues_grouped[inter_or_intra][distance].add(unp_residue2)
    #     return residues_grouped

