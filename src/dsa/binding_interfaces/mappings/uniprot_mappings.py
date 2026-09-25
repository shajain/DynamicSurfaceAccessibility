import json
#from dsa.binding_interfaces.config import uniprot_to_PDB_ids_file, uniprot_to_alphafold_ids_file,  small_uniprot_to_structureIDs_file
from dsa.binding_interfaces.config import UNIPROT_STRUCTURE_MAPPINGS_DIR, small_uniprot_to_structureIDs_file
import matplotlib.pyplot as plt
from dsa.uniprot.config import ALL_UNIPROT_IDS
from collections import defaultdict
import numpy as np
import os



class UniprotStructureMappingsStore:
    MAPPERS = {"pdb": [], "alphafold": []}
    def __init__(self, structure_database: str = "pdb"):
        self.structure_db = structure_database
        if not os.path.exists(self.file_path()):
            self.create_empty_mappings_file()
        self.mapper = self.create_mapper_instance()
        

    def file_path(self):
        return UNIPROT_STRUCTURE_MAPPINGS_DIR / f"uniprot_to_{self.structure_db}_ids.json"

    def create_empty_mappings_file(self):
        if not os.path.exists(UNIPROT_STRUCTURE_MAPPINGS_DIR):
            os.makedirs(UNIPROT_STRUCTURE_MAPPINGS_DIR)
        with open(self.file_path(), "w") as f:
            json.dump({}, f)

    def load_mappings(self):
        with open(self.file_path(), "r") as f:
            return json.load(f)
        
    def create_mapper_instance(self):
        mappings = self.load_mappings()
        if self.structure_db == "pdb":
            return UniprotPDBMapper.factory(mappings=mappings)
        elif self.structure_db == "alphafold":
            return UniprotAlphaFoldMapper.factory(mappings=mappings)
        else:
            raise ValueError(f"Invalid structure database: {self.structure_db}")

    def save_mappings(self, mappings):
        with open(self.file_path(), "w") as f:
            json.dump(mappings, f)

    def update_mappings(self, mappings):
        existing_mappings = self.load_mappings()
        for uniprot in mappings.keys():
            if uniprot in existing_mappings:
                existing_mappings[uniprot].extend(mappings[uniprot])
                existing_mappings[uniprot] = list(set(existing_mappings[uniprot]))
            else:
                existing_mappings[uniprot] = mappings[uniprot]
        self.save_mappings(existing_mappings)
    
    def overwrite_mappings(self, mappings):
        existing_mappings = self.load_mappings()
        existing_mappings.update(mappings)
        self.save_mappings(existing_mappings)

    @classmethod
    def get_mapper(cls, structure_database: str = "pdb"):
        if not cls.MAPPERS[structure_database]:
            store = UniprotStructureMappingsStore(structure_database=structure_database)
            cls.MAPPERS[structure_database] = store.mapper
        return cls.MAPPERS[structure_database]



class UniprotStructureMapper:
    _instance: "UniprotPDBMapper | None" = None
    def __init__(self):
        # self.uniprot_to_strutureIDs_file = uniprot_to_strutureIDs_file
        self.mappings = defaultdict(list)
    
    def set_mappings(self, mappings: dict):
        self.mappings = mappings
    
    def _get_structures(self, uniprot_ids: list[str]):
        dct = {uniprot_id: self.mappings.get(uniprot_id, []) for uniprot_id in uniprot_ids}
        dct = {uniprot_id: [st.lower() for st in structure_ids] for uniprot_id, structure_ids in dct.items()}
        return dct
    
    def uniprot_to_structure_dict(self, uniprot_ids: list[str]=ALL_UNIPROT_IDS, remove_uniprot_without_structure: bool = False):
        dct = self._get_structures(uniprot_ids)
        if remove_uniprot_without_structure:
            dct = {uniprot_id: structure_ids for uniprot_id, structure_ids in dct.items() if len(structure_ids) > 0}
        return dct 

    def structures_for_uniprot_id(self, uniprot_id:str) -> list[str]:
        return self.mappings.get(uniprot_id, [])

    # @classmethod
    # def get_instance(cls, structure_db: str = "pdb"):
    #     if cls._instance is None:
    #         cls._instance = cls()
    #     return cls._instance

    def get_mapped_structures(self) -> list[str]:
        structures = [st.lower() for sts in list(self.mappings.values()) for st in sts]
        return list(set(structures))

    def uniprot_ids_with_structures(self) -> list[str]:
        return list(self.mappings.keys())

    @classmethod
    def factory(cls, mappings: dict):
        mapper = cls()
        mapper.set_mappings(mappings)
        return mapper

    # @classmethod
    # def saved_structures(cls, structure_db: str = "pdb") -> list[str]:
    #     instance = cls.get_instance()
    #     structures = [st for sts in list(instance.uniprot_to_strutureIDs.values()) for st in sts]
    #     return list(set(structures))

    @classmethod
    def representative_ids(cls, structure_db: str = "pdb", save: bool = False, from_saved_file: bool = False):
        if from_saved_file and small_uniprot_to_structureIDs_file[structure_db].exists():
            with open(small_uniprot_to_structureIDs_file[structure_db], "r") as f:
                return json.load(f)
        uniprot_to_structure_ids = cls.get_structures(structure_db=structure_db)
        count_wise_uniprot_ids = defaultdict(list)
        for uniprot_id, structure_ids in uniprot_to_structure_ids.items():
            count_wise_uniprot_ids[len(structure_ids)].append(uniprot_id)
        uniprot_ids = []
        for i in range(1,10):
            # Sample 10 index without replacement
            sampled_indices = np.random.choice(len(count_wise_uniprot_ids[i]), 10, replace=False)
            sampled_uniprot_ids = [count_wise_uniprot_ids[i][j] for j in sampled_indices]
            uniprot_ids.extend(sampled_uniprot_ids)
        representative_ids = {uniprot_id: uniprot_to_structure_ids[uniprot_id] for uniprot_id in uniprot_ids}
        if save:
            with open(small_uniprot_to_structureIDs_file[structure_db], "w") as f:
                json.dump(representative_ids, f)
        return representative_ids

    
    @classmethod
    def plot_structure_stats(cls, structure_db: str = "pdb", xlim: tuple[int, int] = (0, 100)):
        uniprot_to_structure_ids = cls.get_structures(structure_db=structure_db)
        counts = {uniprot_id: len(structure_ids) for uniprot_id, structure_ids in uniprot_to_structure_ids.items()}
        plt.hist(counts.values(), bins=range(xlim[0], xlim[1]+1))
        plt.xlim(xlim)
        plt.show()


class UniprotAlphaFoldMapper(UniprotStructureMapper):
    def __init__(self):
        super().__init__()
    def set_mappings(self, mappings: dict):
        if not mappings:
            mappings = {unp : [unp] for unp in ALL_UNIPROT_IDS}
        super().set_mappings(mappings)



class UniprotPDBMapper(UniprotStructureMapper):
    def __init__(self):
        super().__init__()

    def set_mappings(self, mappings: dict):
        mappings = {unp : [st.lower() for st in db_structs["pdb"]] for unp, db_structs in mappings.items()}
        super().set_mappings(mappings)

    
# def save_uniprot_to_strutureIDs(uniprot_mapping, source, mode='append'):
#     # {uniprot_id: {"pdbids": [...], "AlphaFold2": [...]}
#     # Create a dictionary with uniprot_id as the key. The value should also be a dictionary
#     # with "pdbids" as the key and the pdbids list as the value. It could have other keys corresponding
#     # to AlphaFold2, AlphaFold3, ESM3, etc. containg identifiers for those structures.
#     # mode = "append": Add ids if either the uniport_id doesnot exist or the corresponding value
#     # doesn't contain pdbids key, indicating that for this Uniprot ID the PDB ID was never searched. 
#     # mode = "overwrite": Even if the uniprot_id and pdbids key exists, it should overwrite the pdbids list.
#     with open(uniprot_to_strutureIDs_file, "r") as f:
#         uniprot_to_strutureIds = json.load(f)
#     for uniprot_id in uniprot_mapping.keys():
#         if mode == "append" and uniprot_id in uniprot_to_strutureIds and source in uniprot_to_strutureIds[uniprot_id]:
#             continue
#         ids = uniprot_mapping.get(uniprot_id, [])
#         if uniprot_id not in uniprot_to_strutureIds:
#             uniprot_to_strutureIds[uniprot_id] = {}
#         uniprot_to_strutureIds[uniprot_id][source] = ids
#     with open(uniprot_to_strutureIDs_file, "w") as f:
#         json.dump(uniprot_to_strutureIds, f)

# def remove_invalid_entries(file_path=uniprot_to_strutureIDs_file):
#     with open(file_path, "r") as f:
#         uniprot_to_strutureIds = json.load(f)
#     keys_to_remove = [k for k in uniprot_to_strutureIds.keys() if len(k.split(";")) != 1]
#     removed_keys = []
#     for key in keys_to_remove:
#         if key in uniprot_to_strutureIds:
#             del uniprot_to_strutureIds[key]
#             removed_keys.append(key)
#     with open(file_path, "w") as f:
#         json.dump(uniprot_to_strutureIds, f)
#     return removed_keys

# def get_mapped_uniprot_ids(source, file_path=uniprot_to_strutureIDs_file):
#     with open(file_path, "r") as f:
#         uniprot_to_strutureIds = json.load(f)
#     return [k for k in uniprot_to_strutureIds.keys() if source in uniprot_to_strutureIds[k]]


