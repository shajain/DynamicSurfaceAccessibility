from dsa.dips.mapping.search_dips_in_pdbchains import search_dips_in_pdbchains
from dsa.dips.mapping.get_sifts_data import pdb_sequence_to_uniprot



def _mapping_dict(side_key, dips_seq, matched_chain, matched_chain_seq, score, m=None):
    d = {
        "type": side_key,
        "dips_sequence": dips_seq,
        "matched_chain": matched_chain,
        "matched_chain_seq": matched_chain_seq,
        "match_score": score,
        "uniprot_start": m["uniprot_start"],
        "uniprot_end": m["uniprot_end"],
        "pdb_start": m["pdb_start"],
        "pdb_end": m["pdb_end"],
    }
    return d

def map_dips_to_uniprot(df_dips):
    """
    Input df_dips: columns pdb_id, seq_a, seq_b — DIPS sequences for each chain.

    Output: nested dict
        { uniprot_id: { pdb_id: { "seq_a": {...} | None, "seq_b": {...} | None } } }

    Under each pdb_id, two keys ``seq_a`` and ``seq_b`` each hold the full mapping dict
    for that side when it maps to this UniProt accession, or None if it does not. Each
    mapping dict has: dips_sequence, matched_chain, matched_chain_seq, match_score,
    uniprot_start, uniprot_end, pdb_start, pdb_end.

    The same PDB may appear under more than one UniProt key (different accessions). In
    typical DIPS rows, seq_a and seq_b are often the same sequence (e.g. homodimers), so
    frequently only one side is filled.

    If SIFTS lists several segments for the same chain, the first segment for that
    (uniprot_id, pdb_id, side) wins; later duplicate keys are skipped.
    """
    out = {}

    for _, row in df_dips.iterrows():
        # pdb_id is just one id not a list of ids
        pdb_id = row["pdb_id"]
        seq_a = row["seq_a"]
        seq_b = row["seq_b"]

        uniprot_mapping = pdb_sequence_to_uniprot(pdb_id, seq_a)
        # mapping_resp, seqres_resp = get_sifts_data(pdb_id)
        # chain_seqs = parse_chain_sequences(seqres_resp, pdb_id)
        # #uniprot mapping to populate output with the uniprot id key
        # uniprot_map = parse_uniprot_mapping(mapping_resp, pdb_id)

        # for (dips_seq, side_key) in [(seq_a, "seq_a"), (seq_b, "seq_b")]:
        #     if not isinstance(dips_seq, str) or not dips_seq:
        #         continue
        #     matched_chain, matched_chain_seq, score = search_dips_in_pdbchains(dips_seq, chain_seqs)
        #     if matched_chain:
        #         for m in uniprot_map[matched_chain]:
        #             uid = m["uniprot_id"]
        #             if uid not in out:
        #                 out[uid] = {}
        #             if pdb_id not in out[uid]:
        #                 out[uid][pdb_id] = {}
        #             mapping = _mapping_dict(side_key, dips_seq, matched_chain, matched_chain_seq, score, m)
        #             out[uid][pdb_id].append(mapping)
    return out
