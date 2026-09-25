from collections import defaultdict
import gemmi

class StructureModel:
    def __init__(self, model: gemmi.Model):
        self.model = model
        #self.entity_to_uniprot = entity_to_uniprot
        self._chain_to_subchain, self._subchain_to_chain = self.map_chain_and_subchain()
        self._entity_to_chain = self.map_entity_to_chain()
        #self.chain_to_uniprot, self.uniprot_to_chain = self.map_chain_and_uniprot()
        # self.entity_to_subchain = self.map_entity_to_subchain()
        # self.subchain_to_chain = self.map_subchain_to_chain()
        # self.entity_to_chain = self.map_entity_to_chain()
        # self.chain_to_uniprot = self.map_chain_to_uniprot()

    @property
    def chain_ids(self) -> list[str]:
        return [chain.name for chain in self.model]
    
    @property
    def subchain_ids(self) -> list[str]:
        return list(self._subchain_to_chain.keys())

    def get_chain_reisdue_atom_objects(self, chain: str, 
                                             residue: tuple[str, str], 
                                             atom: tuple[str, str]) -> tuple[gemmi.Chain, gemmi.Residue, gemmi.Atom]:
        chain = self.model.find_chain(chain)
        residue_seq_id = gemmi.SeqId(residue[0])
        residue_name = residue[1]
        residue_group = chain[str(residue_seq_id)]
        for residue in residue_group:
            if residue.name == residue_name:
                break
        #residue = residue_group.find_residue(residue_name)
        atom = residue.find_atom(atom[0], altloc=atom[1])
        return chain, residue, atom

    def map_chain_and_subchain(self) -> dict[str, str]:
        chain_to_subchain = defaultdict(list)
        subchain_to_chain = defaultdict(str)
        for chain in self.model:
            for subchain in chain.subchains():
                chain_to_subchain[chain.name].append(subchain.subchain_id())
                subchain_to_chain[subchain.subchain_id()] = chain.name
        return chain_to_subchain, subchain_to_chain

    def map_entity_to_chain(self) -> dict[str, str]:
        entity_to_chain = defaultdict(set)
        for chain in self.model:
            for subchain in chain.subchains():
                for res in subchain:
                    if res.entity_id is not None and res.entity_id != "":
                        entity_to_chain[res.entity_id].add(chain.name)
                        break
        return entity_to_chain

    def get_chains_containing_entity(self, entity: str) -> list[str]:
        return self._entity_to_chain[entity]

    def get_subchains_in_chain(self, chain: str) -> list[str]:
        return self._chain_to_subchain[chain]

    def get_chain_containing_subchain(self, subchain: str) -> str:
        return self._subchain_to_chain[subchain]
    
    @property
    def entities(self) -> list[str]:
        return list(self._entity_to_chain.keys())
    

    def get_chain_sequence_as_dict(self, chain_id: str) -> str:
        sequence_dict = defaultdict(str)
        chain = self.model.find_chain(chain_id)
        for res in chain:
            sequence_dict[str(res.seqid)] = res.name
        return sequence_dict

    def get_gemmi_model(self) -> gemmi.Model:
        return self.model




    # def map_entity_to_subchain(self) -> dict[str, str]:
    #     entity_to_subchain = defaultdict(list)
    #     for chain in self.model:
    #         for subchain in chain.subchains():
    #             entity_to_subchain[entity].append(subchain)
    #     return entity_to_subchain
    
    # def asym_to_author_chain(self) -> dict[str, str]:
    #     mapping = {}
    #     for chain in self.structure[0]:
    #         # subchains() queries the underlying C++ structure directly
    #         for sub in chain.subchains():
    #             if sub and sub[0].subchain:
    #                 mapping[sub[0].subchain] = chain.name
    #     return mapping