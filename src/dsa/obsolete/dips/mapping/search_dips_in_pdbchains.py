from Bio import SeqIO, pairwise2
import requests

def search_dips_in_pdbchains(dips_seq, chain_sequences, threshold=0.95):
    """
    Align DIPS sequence against all alphabet chains.
    Returns best matching chain_id and identity score.
    """
    best_chain = None
    best_score = -1
    
    for chain_id, chain_seq in chain_sequences.items():
        alignments = pairwise2.align.globalms(
            dips_seq, chain_seq,
            2, -1, -0.5, -0.1,  # match, mismatch, gap_open, gap_extend
            score_only=True
        )
        # Normalize score by length
        max_possible = 2 * min(len(dips_seq), len(chain_seq))
        identity = alignments / max_possible if max_possible > 0 else 0
        
        if identity > best_score:
            best_score = identity
            best_chain = chain_id
            best_chain_seq = chain_seq
    match = True    
    if best_score < threshold:
        match = False
    return match, best_chain, best_chain_seq, best_score