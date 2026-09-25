import gzip
import os
import tempfile
from pathlib import Path
import freesasa
import numpy as np
import pandas as pd
from Bio.SeqUtils import seq1
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure.structure_store import StructureFileStore
from dsa.misc.store import CompressedPickleStore
from dsa.sasa.config import SASA_MONOMER_STORE_DIR


class SASAMonomer:
    """
    Per-residue SASA of the AlphaFold monomer of a uniprot id, computed with FreeSASA
    (ProtOr classifier, which is required for relative SASA to be available).
    """
    DATABASE = "alphafold"

    def __init__(self, uniprot_id: str):
        self.uniprot_id = uniprot_id
        self.sasa_df = self.compute()

    def compute(self) -> pd.DataFrame:
        structure_file_store = ClassFactory.create(StructureFileStore, self.DATABASE)
        path_to_pdb = structure_file_store.get_file_path(self.uniprot_id)
        records = compute_monomer_sasa_records(self.uniprot_id, path_to_pdb)
        return self.index_sasa_df(pd.DataFrame(records))

    @staticmethod
    def index_sasa_df(sasa_df: pd.DataFrame) -> pd.DataFrame:
        sasa_df["resnum"] = sasa_df["resnum"].astype(int)
        sasa_df = sasa_df.set_index("resnum", drop=False)
        sasa_df.index.name = "resnum_idx"
        return sasa_df

    @classmethod
    def from_sasa_df(cls, uniprot_id: str, sasa_df: pd.DataFrame):
        # builds the object from an already computed sasa_df (e.g. the csv files written by compute_sasa.py)
        obj = cls.__new__(cls)
        obj.uniprot_id = uniprot_id
        obj.sasa_df = cls.index_sasa_df(sasa_df)
        return obj

    def get_sasa_for_residue(self, residues: tuple[int, str]|list[tuple[int, str]]|list[int]|int, type="asa"):
        resnums = []
        resnames = []
        if not isinstance(residues, list):
            residues = [residues]
        if len(residues) > 0:
            if isinstance(residues[0], tuple):
                resnums = [residue[0] for residue in residues]
                resnames = [residue[1] for residue in residues]
            elif isinstance(residues[0], int):
                resnums = residues
        missing = [r for r in resnums if r not in self.sasa_df.index]
        if missing:
            return [np.nan] * len(residues)
        if resnames:
            # check if resnames match the resnames in the sasa_df
            if not all(self.sasa_df.loc[resnums, "resname"] == resnames):
                return [np.nan] * len(residues)
        return self.sasa_df.loc[resnums, type].tolist()

    def get_residue_names(self, resnums):
        try:
            res_names = self.sasa_df.loc[resnums, "resname"].tolist()
        except KeyError:
            return ['X' for _ in resnums]
        return res_names

    def summary(self) -> str:
        return f"Number of residues: {len(self.sasa_df)}"

    def __str__(self):
        return f"SASAMonomer-{self.uniprot_id}"


class SASAMonomerStore(CompressedPickleStore):
    def __init__(self):
        super().__init__(directory=SASA_MONOMER_STORE_DIR, extension="sasa", object_type=SASAMonomer)


def compute_monomer_sasa_records(uniprot_id: str, path_to_pdb: Path) -> list[dict]:
    """
    Compute per-residue absolute and relative SASA using FreeSASA's
    built-in normalization (ProtOr radii, Wilke scale internally).

    Returns:
        list of dicts with keys:
            uniprot_id, chain, resnum, resname,
            asa, rel_asa, asa_mainchain, asa_sidechain,
            rel_asa_mainchain, rel_asa_sidechain, asa_nz
    """
    if not path_to_pdb.exists():
        raise FileNotFoundError(f"File {path_to_pdb} does not exist")
    uncompressed_path = None
    if path_to_pdb.suffix == '.gz':
        with gzip.open(path_to_pdb, 'rt') as f_in:
            pdb_content = f_in.read()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            tmp.write(pdb_content)
            uncompressed_path = tmp.name
    try:
        classifier = freesasa.Classifier.getStandardClassifier('protor')
        structure = freesasa.Structure(uncompressed_path or path_to_pdb.as_posix(), classifier)
        result = freesasa.calc(structure)
        residue_areas = result.residueAreas()
    finally:
        if uncompressed_path is not None:
            os.remove(uncompressed_path)

    records = []
    for chain_id, residues in residue_areas.items():
        for res_id, area in residues.items():
            records.append({
                'uniprot_id':        uniprot_id,
                'chain':             chain_id,
                'resnum':            res_id.strip(),
                'resname':           seq1(area.residueType),
                'asa':               round(area.total, 3),
                'rel_asa':           round(area.relativeTotal, 4) if area.relativeTotal is not None else None,
                'asa_mainchain':     round(area.mainChain, 3),
                'asa_sidechain':     round(area.sideChain, 3),
                'rel_asa_mainchain': round(area.relativeMainChain, 4) if area.relativeMainChain is not None else None,
                'rel_asa_sidechain': round(area.relativeSideChain, 4) if area.relativeSideChain is not None else None,
                'asa_nz':            None   # filled below
            })

    nz_lookup = {}
    for i in range(structure.nAtoms()):
        if (structure.atomName(i).strip() == 'NZ' and
            structure.residueName(i).strip() == 'LYS'):
            key = (structure.chainLabel(i), structure.residueNumber(i).strip())
            nz_atom_area = result.atomArea(i)
            nz_lookup[key] = round(nz_atom_area, 3) if nz_atom_area is not None else None

    for rec in records:
        key = (rec['chain'], rec['resnum'])
        if key in nz_lookup:
            rec['asa_nz'] = nz_lookup[key]

    return records
