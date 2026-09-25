from dsa.uniprot.uniprot_to_sequence import UniProtToSequence
from Bio.SeqUtils import seq1, seq3
import numpy as np
from typing import Literal


def validate_against_sequence(uniprot_id: str, position_residue_pairs: list[tuple[int, str]], 
                                    mode: Literal["strict", "lenient"] = "strict") -> bool:
    mismatches = mismatches_against_sequence(uniprot_id, position_residue_pairs)
    num_mismatches = len(mismatches)
    fraction_mismatches = num_mismatches/len(position_residue_pairs)
    if mode == "strict":
        match = num_mismatches==0
    elif mode == "lenient":
        match = fraction_mismatches>0.1
    return match


def mismatches_against_sequence(uniprot_id: str, position_residue_pairs: list[tuple[int, str]]):
    sequence = UniProtToSequence.get_sequence(uniprot_id)
    mismatches = [(p, r, sequence[p-1]) for p, r in position_residue_pairs if seq1(r)!=sequence[p-1]]
    return mismatches

def extract_Lysines(cls, uniprot_id: str):
    sequence = UniProtToSequence.get_sequence(uniprot_id)
    lysines = [p for p, r in position_residue_pairs if r=='LYS']
    return lysines