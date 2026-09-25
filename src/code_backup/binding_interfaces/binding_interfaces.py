from dsa.binding_interfaces.config import ALL_UNIPROT_IDS
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappings
from collections import defaultdict
from dsa.binding_interfaces.process_structure.Structure import Structure
from typing import Literal
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence
from Bio.SeqUtils import seq1
from dsa.uniprot.config import example_uniprot_id
from dsa.binding_interfaces.config import BI_PROTEIN_DIR
import json


class BindingInterfaces:
    _instances: dict[tuple[str, float], "BindingInterfaces"] = {}
    INTERFACE_TYPES = ["all", "residue-residue", "residue-other", "intra-chain", "inter-chain"]
    OFFSETS = [-1, 0, 1, -2, 2, -3, 3]
    ATTRIBUTES_TO_SAVE = ["binding_interfaces", "unaligned_residues", "offsets", "is_aligned", "all_structures_processed"]

    def __init__(self, uniprot_ids: list[str] = ALL_UNIPROT_IDS, 
                       structure_db: str = "pdb", 
                       cutoff: float = 8.0, 
                       save: bool = True,
                       saved_binding_interfaces_only: bool = False):
        self.structure_db = structure_db
        self.uniprot_ids = uniprot_ids
        self.cutoff = cutoff
        self.uniprot_to_strutureIDs = UniprotStructureMappings.get_structures(self.uniprot_ids, self.structure_db, remove_uniprot_without_structure=True)
        self.all_structureIDs = self.get_all_structureIDs()
        self.file_path = BI_PROTEIN_DIR / ("_".join(["BI_protein", self.structure_db, str(self.cutoff)]) + ".json")
        self.binding_interfaces = defaultdict(lambda: defaultdict(list))
        self.unaligned_residues = defaultdict(lambda: defaultdict(list))
        self.offsets = defaultdict(lambda: defaultdict(int))
        self.is_aligned = defaultdict(lambda: defaultdict(bool))
        # It seems that there is some issue either with the cif file or the gemmi library that some residues that are on 
        # Polypeptide chains have nucleotides as residue name. See 7XCR entity 2 with chains B and J. It corresponds to 
        # protein A0A672GII6. Incorrect residues store the interface residues that are not protein residues. 
        # If you want to go through this example search for P62805 interfaces. A0A672GII6 is not in the uniprot_to_strutureIDs file.
        self.incorrect_residues = defaultdict(lambda: defaultdict(list))
        self.all_structures_processed = list()
        self.update_binding_interfaces_from_saved_file()
        new_structures = list(set(self.all_structureIDs) - set(self.all_structures_processed))
        if len(new_structures) > 0 and not saved_binding_interfaces_only:
            self.extract_binding_interfaces(new_structures)
            if save:
                self.save_binding_interfaces()
        



    def get_all_structureIDs(self):
        structures = [st for unp, sts in self.uniprot_to_strutureIDs.items() for st in sts]
        return list(set(structures))

    def processing_statistics(self):
        number_of_structures = len(self.all_structures_processed)
        number_of_structure_protein_pairs = sum([len(str_to_bi.keys()) for str_to_bi in self.binding_interfaces.values()])
        number_of_unaligned_structure_protein_pairs = sum([sum(is_aligned.values()) for is_aligned in self.is_aligned.values()])
        print(f"Total processed {number_of_structures} structures including saved ones")
        print(f"Total number of structure-protein pairs: {number_of_structure_protein_pairs}")
        print(f"Total number of unaligned structure-protein pairs: {number_of_unaligned_structure_protein_pairs}")
        return number_of_structures, number_of_structure_protein_pairs, number_of_unaligned_structure_protein_pairs

    def extract_binding_interfaces(self, structures):
        num_structures = 0
        for structure_id in structures:
            self.process_structure(structure_id)
            num_structures += 1
            if num_structures%500 == 0:
                print(f"Processed {num_structures} structures")
                self.processing_statistics()
        self.save_binding_interfaces()

    def process_structure(self, structure_id):
        contact_residues = Structure.extract_binding_interfaces(structure_id, database=self.structure_db, cutoff=self.cutoff, only_from_saved_file=True)
        if contact_residues is None:
            return
        self.all_structures_processed.append(structure_id)
        for uniprot_id, crs in contact_residues.items():
            #print(f"Uniprot:{uniprot_id} structure:{structure_id}")
            offset = None
            for ctype in BindingInterfaces.INTERFACE_TYPES:
                position_residue_pairs = crs[ctype]
                validation = self.validate_binding_interface_alignment(uniprot_id, position_residue_pairs, mode="lenient", offset=offset)
                if ctype == "all":
                    self.unaligned_residues[uniprot_id][structure_id].extend(validation["unaligned"])
                    self.incorrect_residues[uniprot_id][structure_id].extend(validation["incorrect"])
                    self.offsets[uniprot_id][structure_id] = validation["offset"]
                    self.is_aligned[uniprot_id][structure_id] = validation["is_aligned"]
                if validation["is_aligned"]:
                    self.binding_interfaces[uniprot_id][ctype].extend(validation["aligned"])
                    self.binding_interfaces[uniprot_id][ctype] = list(set(self.binding_interfaces[uniprot_id][ctype]))
                    offset = validation["offset"]
                else: 
                    break

    def _binding_interface(self, uniprot_id: str, ctype: Literal[*BindingInterfaces.INTERFACE_TYPES]="all") -> list[tuple[int, str]]:
        return self.binding_interfaces[uniprot_id][ctype]

    def update_binding_interfaces_from_saved_file(self):
        if not self.file_path.exists():
            return
        with open(self.file_path, "r") as f:
            data = json.load(f)
            BindingInterfaces.copy_to_object(self.__dict__, data)

    def save_binding_interfaces(self):
        self.file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.file_path, "w") as f:
            data = {key: self.__dict__[key] for key in BindingInterfaces.ATTRIBUTES_TO_SAVE}
            json.dump(data, f, indent=4)

    @staticmethod
    def copy_to_object(d1, d2):
        # this function is important to ensure that the saved dictionary variables when used to initialize the 
        # class variables the type of the class variables as defaultdict is preserved.
        if isinstance(d2, dict):
            for k, v in d2.items():
                if isinstance(v, dict):
                    d1[k] = BindingInterfaces.copy_to_object(d1[k], v)
                elif isinstance(v, list) and len(v)>0 and isinstance(v[0], list) and len(v[0]) ==2:
                    d1[k] = [(int(p), r) for p, r in v]
                else:
                    d1[k] = v
        else:
            d1 = d2
        return d1


    
    @staticmethod
    def validate_binding_interface_alignment(uniprot_id: str, position_residue_pairs: list[tuple[int, str]], 
                                                        mode: Literal["strict", "lenient"] = "lenient",
                                                        offset: int = None) -> bool:                                     
        if len(position_residue_pairs)==0:
            return {"is_aligned": True, "fraction_unaligned": 0.0, "unaligned": [], "aligned": [], "incorrect": [], "offset": None}
        offsets = BindingInterfaces.OFFSETS
        if offset is not None:
            offsets = [offset]

        best_offset = offsets[0]
        best_frac_unaligned = 1.0 
        for off in offsets:
            unaligned, aligned, incorrect = BindingInterfaces.align_against_sequence(uniprot_id, position_residue_pairs, offset=off)
            num_unaligned = len(unaligned)
            num_aligned = len(aligned)
            if num_unaligned + num_aligned == 0:
                fraction_unaligned = 1.0
            else:
                fraction_unaligned = num_unaligned/(num_unaligned + num_aligned)
            if mode == "strict":
                is_aligned = num_unaligned==0
            elif mode == "lenient":
                is_aligned = fraction_unaligned<=0.1
            if fraction_unaligned < best_frac_unaligned:
                best_frac_unaligned = fraction_unaligned
                best_offset = off
            if is_aligned:
                break
        if len(offsets)>1:
            output = BindingInterfaces.validate_binding_interface_alignment(uniprot_id, position_residue_pairs, mode="lenient", offset=best_offset)
        else:
            output = {
                "is_aligned": is_aligned,
                "fraction_unaligned": fraction_unaligned,
                "unaligned": unaligned,
                "aligned": aligned,
                "incorrect": incorrect,
                "offset": best_offset
            }
        return output
        
    @staticmethod
    def align_against_sequence(uniprot_id: str, position_residue_pairs: list[tuple[int, str]], offset: int = -1):
        #print(uniprot_id)
        sequence = UniprotToSequence.get_sequence(uniprot_id, only_saved_sequences=True)
        if sequence is None:
            return [], [], []
        unaligned = []
        aligned = []
        incorrect = []
        for (p, r) in position_residue_pairs:
            if seq1(r) == 'X':
                incorrect.append((p, r))
                continue
            if 0 <= p+offset < len(sequence):
                seq_res = sequence[p+offset]
                if seq1(r) == seq_res:
                    aligned.append((p, r))
                else:
                    unaligned.append((p, r, seq_res))
            else:
                unaligned.append((p, r, 'None'))
        # unaligned_residues = [(p, r, sequence[p+offset]) for p, r in position_residue_pairs 
        #                         if 0 <= p+offset < len(sequence) 
        #                         and seq1(r)!=sequence[p+offset]]
        # aligned_residues = [(p, r) for p, r in position_residue_pairs 
        #                     if 0 <= p+offset < len(sequence) 
        #                     and seq1(r)==sequence[p+offset]]
        return unaligned, aligned, incorrect

    @classmethod
    def binding_interface(cls, uniprot_id: str, structure_db: str, cutoff: float = 8.0, 
                            ctype: Literal[*BindingInterfaces.INTERFACE_TYPES]="all") -> list[tuple[int, str]]:
        if cls._instances.get((structure_db, cutoff)) is None:
            cls._instances[(structure_db, cutoff)] = cls(ALL_UNIPROT_IDS, structure_db=structure_db, cutoff=cutoff, saved_binding_interfaces_only=True)
        return cls._instances[(structure_db, cutoff)]._binding_interface(uniprot_id, ctype)

if __name__ == "__main__":
    BI = BindingInterfaces(ALL_UNIPROT_IDS, structure_db="pdb", cutoff=8.0, save=True)
    BindingInterfaces.binding_interface(example_uniprot_id, "pdb", 8.0, "all")