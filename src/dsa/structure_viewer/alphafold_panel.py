"""
AlphaFold models shown next to a PDB structure, superposed on each PDB chain copy of their protein.

For every UniProt id of the PDB structure the AlphaFold model is loaded (as its own StructureData, with its contacts
and SASA). For every bio chain of that protein, the AlphaFold model is fitted onto the chain by least squares (Kabsch)
on the CA atoms of the residues present in both, matched by UniProt number and residue type. The page applies the
fit of the selected chain to the AlphaFold coordinates, so both panels share one frame and can be rotated together.
"""
import warnings
import numpy as np
from dsa.structure_viewer.residue_data import StructureData, load_alphafold

MAX_ALPHAFOLD_MODELS = 10     # AlphaFold models embedded in one page at most
MIN_FITTED_RESIDUES = 10      # fewer matched CA atoms than this: no fit, the panel stays empty


def ca_positions(data: StructureData) -> dict[int, np.ndarray]:
    """CA coordinates by residue index (residues are in model order)."""
    positions = {}
    index = 0
    for chain in data.model[0]:
        for residue in chain:
            ca = residue.find_atom("CA", "*")
            if ca is not None:
                positions[index] = np.array([ca.pos.x, ca.pos.y, ca.pos.z])
            index += 1
    return positions


def kabsch(moving: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Rotation R and translation t minimising |R @ moving + t - target| (points as rows), and the RMSD after the fit."""
    moving_center, target_center = moving.mean(axis=0), target.mean(axis=0)
    h = (moving - moving_center).T @ (target - target_center)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    translation = target_center - rotation @ moving_center
    fitted = moving @ rotation.T + translation
    rmsd = float(np.sqrt(((fitted - target) ** 2).sum(axis=1).mean()))
    return rotation, translation, rmsd


def load_alphafold_models(pdb: StructureData) -> dict[str, StructureData]:
    """AlphaFold models of the UniProt ids of the PDB structure (in order of appearance, at most MAX_ALPHAFOLD_MODELS)."""
    uniprot_ids = list(dict.fromkeys(pdb.residues["uniprot_id"].dropna()))
    if len(uniprot_ids) > MAX_ALPHAFOLD_MODELS:
        warnings.warn(f"{len(uniprot_ids)} UniProt ids; AlphaFold models are added for the first {MAX_ALPHAFOLD_MODELS} only")
        uniprot_ids = uniprot_ids[:MAX_ALPHAFOLD_MODELS]
    models = {}
    for uniprot_id in uniprot_ids:
        try:
            models[uniprot_id] = load_alphafold(uniprot_id)
        except Exception as e:
            warnings.warn(f"No AlphaFold model for {uniprot_id}: {e}")
    return models


def fit_chains(pdb: StructureData, models: dict[str, StructureData]) -> dict[str, dict]:
    """
    Per bio chain: the UniProt id of the chain (the most common one), and, if its AlphaFold model is available and at
    least MIN_FITTED_RESIDUES CA atoms match, the rotation / translation that put the AlphaFold model on the chain.
    """
    residues = pdb.residues
    pdb_ca = ca_positions(pdb)
    model_ca = {uniprot_id: ca_positions(model) for uniprot_id, model in models.items()}
    model_index = {uniprot_id: {(int(n), letter): i for i, n, letter, kind in
                                zip(model.residues.index, model.residues["resnum"], model.residues["resname1"], model.residues["kind"])
                                if kind == "protein"}
                   for uniprot_id, model in models.items()}
    fits = {}
    for bio_chain, chain_residues in residues.groupby("bio_chain", sort=False):
        mapped = chain_residues[chain_residues["uniprot_id"].notna() & chain_residues["unp_resnum"].notna()]
        if mapped.empty:
            continue
        uniprot_id = mapped["uniprot_id"].mode().iloc[0]
        fit = {"uniprot_id": uniprot_id, "rotation": None, "translation": None, "rmsd": None, "n": 0}
        if uniprot_id in models:
            pairs = []
            for i, unp_resnum, letter in zip(mapped.index, mapped["unp_resnum"], mapped["resname1"]):
                if mapped.at[i, "uniprot_id"] != uniprot_id:
                    continue
                j = model_index[uniprot_id].get((int(unp_resnum), letter))
                if j is not None and i in pdb_ca and j in model_ca[uniprot_id]:
                    pairs.append((model_ca[uniprot_id][j], pdb_ca[i]))
            fit["n"] = len(pairs)
            if len(pairs) >= MIN_FITTED_RESIDUES:
                moving, target = np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs])
                rotation, translation, rmsd = kabsch(moving, target)
                fit.update(rotation=rotation.round(6).tolist(), translation=translation.round(4).tolist(), rmsd=round(rmsd, 2))
        fits[bio_chain] = fit
    return fits
