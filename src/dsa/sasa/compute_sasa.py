import gzip
from subprocess import list2cmdline
from dsa.config import DATA_DIR
import freesasa
#from dsa.binding_interfaces.mappings.uniprot_to_alphaFold import get_alphafold_structure
from dsa.binding_interfaces.structure.structure import AlphaFoldStructure
from dsa.uniprot.config import ALL_UNIPROT_IDS
import pandas as pd
import tempfile
import os
from Bio.SeqUtils import seq1
import warnings
import numpy as np
SASA_DIR = DATA_DIR / "sasa"
SASA_MONOMER_DIR = SASA_DIR / "monomers"

class SASAMonomerStore:
    def __init__(self):
        self.sasa_monomer_dir = SASA_MONOMER_DIR

    def get_sasa_df_path(self, uniprot_id):
        return self.sasa_monomer_dir / f"{uniprot_id}.csv"

    def is_saved(self, uniprot_id):
        return self.get_sasa_df_path(uniprot_id).exists()

    def load(self, uniprot_id):
        path_to_sasa_df = self.get_sasa_df_path(uniprot_id)
        if self.is_saved(uniprot_id):
            sasa_df = pd.read_csv(path_to_sasa_df)
            return SASAMonomer.factory_method(uniprot_id, sasa_df)
        else:
            return None

    def saved_uniprot_ids(self):
        return [file.stem for file in self.sasa_monomer_dir.glob("*.csv")]

    def unsaved_uniprot_ids(self, uniprot_ids: list[str]):
        return [uniprot_id for uniprot_id in uniprot_ids if not self.is_saved(uniprot_id)]

    def compute_and_save(self, uniprot_id: str):
        try:
            sasa_monomer = SASAMonomer.factory_method(uniprot_id)
            sasa_df = sasa_monomer.sasa_df
            sasa_df.to_csv(self.get_sasa_df_path(uniprot_id), index=False)
        except Exception as e:
            print(f"Error computing and saving SASA for {uniprot_id}: {e}")
            return None



class SASAMonomer:
    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id
        # self.sasa_dir = SASA_DIR
        # self.sasa_monomer_dir = SASA_MONOMER_DIR
        self.sasa_df = None

    def compute_monomer_sasa(self):
        sasa_records = compute_monomer_sasa(self.uniprot_id)
        if sasa_records is None:
            return
        self.sasa_df = pd.DataFrame(sasa_records)
        self.index_sasa_df()

    def index_sasa_df(self):
        self.sasa_df["resnum"] = self.sasa_df["resnum"].astype(int)
        self.sasa_df = self.sasa_df.set_index("resnum", drop=False)
        self.sasa_df.index.name = "resnum_idx"

    def set_sasa_df(self, sasa_df):
        self.sasa_df = sasa_df

    def get_sasa_for_residue(self, residues: tuple[int, str]|list[tuple[int, str]]|list[int]|int, type="asa"):
        resnums = []
        resnames = []
        sasas = []
        # if isinstance(residues, tuple):
        #     resnums = residues[0]
        #     resnames = residues[1]
        if not isinstance(residues, list):
            residues = [residues]
        if len(residues) > 0:
            if isinstance(residues[0], tuple):
                resnums = [residue[0] for residue in residues]
                resnames = [residue[1] for residue in residues]
            elif isinstance(residues[0], int):
                resnums = residues
        # else:
        #     raise ValueError(f"Invalid residues: {residues}")
        # if type == "asa":
        #     return self.sasa_df.loc[resnums, "asa"].tolist()
        # elif type == "rel_asa":
        #     return self.sasa_df.loc[resnums, "rel_asa"].tolist()
        # elif type == "asa_nz":
        #     return self.sasa_df.loc[resnums, "asa_nz"].tolist()
        # else:
        #     raise ValueError(f"Invalid type: {type}")
        missing = [r for r in resnums if r not in self.sasa_df.index]
        if missing:
            # warnings.warn(f"Missing resnums {missing} for {self.uniprot_id}")
            return [np.nan] * len(residues)
        if resnames:
            # check if resnames match the resnames in the sasa_df
            if not all(self.sasa_df.loc[resnums, "resname"] == resnames):
                # warnings.warn(f"Resnames do not match the resnames in the sasa_df for residues {residues} in uniprot \
                #         {self.uniprot_id}")
                return [np.nan] * len(residues)
        sasas = self.sasa_df.loc[resnums, type].tolist()
        return sasas

    def get_residue_names(self, resnums):
        try:
            res_names = self.sasa_df.loc[resnums, "resname"].tolist()
        except KeyError:
            return ['X' for _ in resnums]
        return res_names
        #return self.sasa_df.loc[resnums, "resname"].tolist()

    @classmethod
    def factory_method(cls, uniprot_id, sasa_df=None):
        obj = cls(uniprot_id)
        if sasa_df is not None:
            obj.set_sasa_df(sasa_df)
            obj.index_sasa_df()
        else:
            obj.compute_monomer_sasa()
        return obj


def compute_monomer_sasa(uniprot_id):
    """
    Compute per-residue absolute and relative SASA using FreeSASA's
    built-in normalization (ProtOr radii, Wilke scale internally).

    Args:
        uniprot_id          : UniProt accession e.g. 'P12345'
        alphafold_interface : object with .get_pdb_string(uniprot_id) method

    Returns:
        list of dicts with keys:
            uniprot_id, chain, resnum, resname,
            asa, rel_asa, asa_mainchain, asa_sidechain, asa_polar, asa_apolar
    """
    #path_to_pdb = get_alphafold_structure(uniprot_id, format="pdb")
    path_to_pdb = AlphaFoldStructure(uniprot_id).get_file_path()
    failure = False
    # with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
    #     f.write(pdb_string)
    #     tmp_path = f.name
    uncompressed_path = None
    # this is a path object not a string
    if path_to_pdb.suffix == '.gz':
        with gzip.open(path_to_pdb, 'rt') as f_in:
            pdb_content = f_in.read()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            tmp.write(pdb_content)
            uncompressed_path = tmp.name
    try:
        # use ProtOr classifier — required for relative SASA to be available
        classifier = freesasa.Classifier.getStandardClassifier('protor')
        structure  = freesasa.Structure(uncompressed_path, classifier)
        result     = freesasa.calc(structure)
        residue_areas = result.residueAreas()   # rel_asa built-in
    except Exception as e:
        print(f"Error computing SASA for {uniprot_id}: {e}")
        failure = True
    finally:
        if uncompressed_path is not None:
            os.remove(uncompressed_path)
    if failure:
        return None
    records = []
    for chain_id, residues in residue_areas.items():
        for res_id, area in residues.items():
            records.append({
                'uniprot_id':    uniprot_id,
                'chain':         chain_id,
                'resnum':        res_id.strip(),
                'resname':       seq1(area.residueType),
                'asa':           round(area.total, 3),
                'rel_asa':       round(area.relativeTotal, 4) if area.relativeTotal is not None else None,
                'asa_mainchain': round(area.mainChain, 3),
                'asa_sidechain': round(area.sideChain, 3),
                'rel_asa_mainchain': round(area.relativeMainChain, 4) if area.relativeMainChain is not None else None,
                'rel_asa_sidechain': round(area.relativeSideChain, 4) if area.relativeSideChain is not None else None,
                'asa_nz':        None   # fill below in loop
            })
    nz_lookup = {}
    for i in range(structure.nAtoms()):
        if (structure.atomName(i).strip() == 'NZ' and
            structure.residueName(i).strip() == 'LYS'):
            key = (structure.chainLabel(i), structure.residueNumber(i).strip())
            nz_atom_area = result.atomArea(i)
            nz_lookup[key] = round(nz_atom_area, 3) if nz_atom_area is not None else None

    # fill NZ into records
    for rec in records:
        key = (rec['chain'], rec['resnum'])
        if key in nz_lookup:
            rec['asa_nz'] = nz_lookup[key]

    return records




# --- AlphaFold interface ---------------------------------------------------

# import requests

# class AlphaFold_DB:
#     def get_pdb_string(self, uniprot_id):
#         file
#         url = (f"https://alphafold.ebi.ac.uk/files/"
#                f"AF-{uniprot_id}-F1-model_v{self.version}.pdb")
#         r = requests.get(url)
#         if r.status_code == 200
#             return r.text
#         raise ValueError(f"No AlphaFold structure for {uniprot_id} "
#                          f"(HTTP {r.status_code})")


# --- usage ----------------------------------------------------------------

if __name__ == "__main__":
    import random
    refresh_interval = 10
    save = True

    sasa_monomer_store = SASAMonomerStore()

    # UniProt ids that still need processing: all ids minus those already saved.
    uniprot_ids_to_process = sasa_monomer_store.unsaved_uniprot_ids(ALL_UNIPROT_IDS)
    print(f"{len(uniprot_ids_to_process)} uniprot ids to process")
    saved_count = 0
    # A single lightweight extractor reused across iterations.
    while uniprot_ids_to_process:
        # Pick a random uniprot id so that parallel threads are unlikely to collide.
        uniprot_id = uniprot_ids_to_process.pop(random.randrange(len(uniprot_ids_to_process)))
        # The store re-checks if it is already saved (possibly by another thread) and only
        # extracts when needed; returns None if it was already saved.
        if save:
            sasa_monomer = sasa_monomer_store.compute_and_save(uniprot_id)
        else:
            sasa_monomer = SASAMonomer.factory_method(uniprot_id)
        if sasa_monomer is None:
            continue
        saved_count += 1
        if not save:
            print(sasa_monomer.sasa_df.head(100))
        # Periodically drop ids that were saved by other threads.
        if saved_count % refresh_interval == 0:
            uniprot_ids_to_process = sasa_monomer_store.unsaved_uniprot_ids(uniprot_ids_to_process)
            print(f"Refreshed queue: {len(uniprot_ids_to_process)} uniprot ids remaining")
