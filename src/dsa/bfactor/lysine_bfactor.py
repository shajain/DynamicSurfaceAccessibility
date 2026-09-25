import warnings
from collections import defaultdict
import gemmi
import numpy as np
import pandas as pd
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure.structure import Structure
from dsa.misc.store import CompressedPickleStore
from dsa.bfactor.config import BFACTOR_STORE_DIR

PEPTIDE_POLYMER_TYPES = (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)


class LysineBFactor:
    """
    B-factors of the lysines of a PDB structure: the B-factor of the NZ atom and the average
    B-factor over the heavy atoms of each lysine residue.

    B-factors are read from the first model of the asymmetric unit (copies of a chain in the
    biological assembly share the same B-factors), after removing hydrogens and keeping only the
    first alternate conformation. There is one row per lysine per author chain.

    Raw B-factors depend on resolution and refinement, so they are also given as z-scores
    normalized over all protein heavy atoms of the same chain (*_z columns).
    """
    DATABASE = "pdb"

    def __init__(self, structure_id: str):
        self.structure_id = structure_id
        self.method = None
        self.resolution = None
        self.bfactor_df = self.compute()

    def compute(self) -> pd.DataFrame:
        structure = ClassFactory.create(Structure, self.DATABASE, self.structure_id)
        self.method = dict(structure.gemmi_structure.info).get("_exptl.method") or None
        self.resolution = structure.gemmi_structure.resolution or None
        residues = self.protein_residues(structure)
        if not residues:
            raise ValueError(f"No protein residues in {self.structure_id}")
        chain_stats = self.chain_bfactor_stats(residues)
        records = [self.make_record(chain, residue, chain_stats[chain])
                   for chain, residue in residues if residue.name == "LYS"]
        if not records:
            raise ValueError(f"No lysines in {self.structure_id}")
        bfactor_df = pd.DataFrame(records)
        self.add_uniprot(bfactor_df, structure, [(c, r) for c, r in residues if r.name == "LYS"])
        float_columns = ["nz_bfactor", "avg_bfactor", "nz_bfactor_z", "avg_bfactor_z"]
        bfactor_df[float_columns] = bfactor_df[float_columns].astype(float)
        bfactor_df["uniprot_id"] = bfactor_df["uniprot_id"].astype(object).where(bfactor_df["uniprot_id"].notna(), None)
        return bfactor_df

    @staticmethod
    def protein_residues(structure: Structure) -> list[tuple[str, gemmi.Residue]]:
        # residues of the asymmetric unit that belong to protein entities, as (author_chain, residue)
        peptide_entities = {e.name for e in structure.gemmi_structure.entities
                            if e.entity_type == gemmi.EntityType.Polymer and e.polymer_type in PEPTIDE_POLYMER_TYPES}
        model = structure.gemmi_structure[0].clone()
        model.remove_alternative_conformations()
        model.remove_hydrogens()
        return [(chain.name, residue.clone()) for chain in model for residue in chain
                if residue.entity_id in peptide_entities]

    @staticmethod
    def chain_bfactor_stats(residues: list[tuple[str, gemmi.Residue]]) -> dict[str, tuple[float, float]]:
        # (mean, std) of the B-factors of all protein heavy atoms of each chain
        chain_bfactors = defaultdict(list)
        for chain, residue in residues:
            chain_bfactors[chain].extend(atom.b_iso for atom in residue)
        return {chain: (float(np.mean(b)), float(np.std(b))) for chain, b in chain_bfactors.items()}

    def make_record(self, chain: str, residue: gemmi.Residue, chain_stats: tuple[float, float]) -> dict:
        nz = next((atom for atom in residue if atom.name == "NZ"), None)
        nz_bfactor = nz.b_iso if nz is not None else None
        avg_bfactor = float(np.mean([atom.b_iso for atom in residue])) if len(residue) > 0 else None
        return {
            'structure_id':  self.structure_id,
            'author_chain':  chain,
            'resnum':        residue.seqid.num,
            'icode':         residue.seqid.icode.strip(),
            'resname':       "K",
            'n_atoms':       len(residue),   # < 9 heavy atoms means the side chain is truncated
            'nz_bfactor':    round(nz_bfactor, 3) if nz_bfactor is not None else None,
            'avg_bfactor':   round(avg_bfactor, 3) if avg_bfactor is not None else None,
            'nz_bfactor_z':  z_score(nz_bfactor, chain_stats),
            'avg_bfactor_z': z_score(avg_bfactor, chain_stats),
            'uniprot_id':    None,   # filled in add_uniprot
            'unp_resnum':    None,
            'unp_resname':   None,
        }

    def add_uniprot(self, bfactor_df: pd.DataFrame, structure: Structure, lysines: list[tuple[str, gemmi.Residue]]):
        try:
            entity_to_uniprot = defaultdict(list)
            for uniprot_id, entities in structure.uniprot_to_entities.items():
                for entity in entities:
                    entity_to_uniprot[entity].append(uniprot_id)
            unp = [self.map_to_uniprot(structure, chain, residue, entity_to_uniprot[residue.entity_id])
                   for chain, residue in lysines]
            bfactor_df["uniprot_id"] = pd.Series([u[0] for u in unp], dtype=object)
            bfactor_df["unp_resnum"] = pd.array([u[1] for u in unp], dtype="Int64")
            bfactor_df["unp_resname"] = pd.Series([u[2] for u in unp], dtype=object)
        except Exception as e:
            # keep the B-factors even if the uniprot alignment fails; uniprot columns stay empty
            warnings.warn(f"UniProt mapping failed for {self.structure_id}: {e}")

    @staticmethod
    def map_to_uniprot(structure: Structure, author_chain: str, residue: gemmi.Residue,
                       uniprot_ids: list[str]) -> tuple[str, int, str]:
        for uniprot_id in uniprot_ids:
            unp_residue = structure.map_residue_to_uniprot(author_chain, residue, uniprot_id)
            if unp_residue is not None:
                return uniprot_id, unp_residue[0], unp_residue[1]
        return None, None, None

    def get_bfactor_for_residue(self, uniprot_id: str, residues: tuple[int, str]|list[tuple[int, str]]|list[int]|int,
                                type="nz_bfactor") -> list[list[float]]:
        """
        For each uniprot residue, the list of values from every chain that maps to it
        (empty list if the residue is not a lysine in the structure).
        """
        if not isinstance(residues, list):
            residues = [residues]
        df = self.bfactor_df[self.bfactor_df["uniprot_id"] == uniprot_id]
        values = []
        for residue in residues:
            resnum, resname = residue if isinstance(residue, tuple) else (residue, None)
            rows = df[df["unp_resnum"] == resnum]
            if resname is not None:
                rows = rows[rows["resname"] == resname]
            values.append(rows[type].tolist())
        return values

    @property
    def uniprot_ids(self) -> list[str]:
        return [u for u in self.bfactor_df["uniprot_id"].dropna().unique()]

    def summary(self) -> str:
        n_nz = self.bfactor_df["nz_bfactor"].notna().sum()
        return (f"{self.structure_id} ({self.method}, {self.resolution} A): "
                f"{self.bfactor_df['author_chain'].nunique()} chains, {len(self.bfactor_df)} lysines, "
                f"{n_nz} with NZ, uniprot ids {self.uniprot_ids}")

    def __str__(self):
        return f"LysineBFactor-{self.structure_id}"


class LysineBFactorStore(CompressedPickleStore):
    def __init__(self):
        super().__init__(directory=BFACTOR_STORE_DIR, extension="bfactor", object_type=LysineBFactor)


def z_score(value: float|None, chain_stats: tuple[float, float]) -> float|None:
    # None when the chain has no B-factor spread (e.g. NMR models with all B-factors set to 0)
    mean, std = chain_stats
    if value is None or std == 0:
        return None
    return round((value - mean) / std, 4)
