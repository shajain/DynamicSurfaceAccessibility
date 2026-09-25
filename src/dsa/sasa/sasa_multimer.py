import json
import os
import shutil
import string
import subprocess
import sys
import tempfile
import warnings
from collections import defaultdict
from pathlib import Path
import gemmi
import numpy as np
import pandas as pd
from Bio.SeqUtils import seq1
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure.structure import Structure
from dsa.misc.store import CompressedPickleStore
from dsa.sasa.config import SASA_DIR

SASA_MULTIMER_STORE_DIR = SASA_DIR / "multimers"

PEPTIDE_POLYMER_TYPES = (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)
# freesasa 2.1.2 keeps at most 5 characters of a residue label (sign + number + insertion code), so
# labels are 1..99999 without insertion code, then -999..9999 with insertion code A-Z, a-z
# (digits are not used as insertion codes: 1 + '0' would read as 10)
PLAIN_LABEL_MAX = 99999
INSERTION_CODES = string.ascii_uppercase + string.ascii_lowercase
LABEL_NUMBERS_WITH_INSERTION_CODE = range(-999, 10000)
MAX_RESIDUE_LABELS = PLAIN_LABEL_MAX + len(INSERTION_CODES) * len(LABEL_NUMBERS_WITH_INSERTION_CODE)


class SASAMultimer:
    """
    Per-residue SASA of the biological assembly of a PDB structure, computed with the FreeSASA
    command line tool (ProtOr radii) on the whole complex, so burial by other chains is included.

    Only protein residues are kept (modified amino acids included); nucleic acids, ligands,
    cofactors, glycans, waters, hydrogens and alternate conformations are removed.

    freesasa truncates chain labels to one character, so copies of a chain (A1, A2, ...) cannot be
    told apart by chain. Instead, every residue is given a unique residue label (number + insertion
    code, see residue_labels) in the file sent to freesasa, and the results are mapped back to the
    bio chain / author chain / uniprot residue through that label. There is one row per residue per
    chain copy.
    """
    DATABASE = "pdb"

    def __init__(self, structure_id: str):
        self.structure_id = structure_id
        self.sasa_df = self.compute()

    def compute(self) -> pd.DataFrame:
        structure = ClassFactory.create(Structure, self.DATABASE, self.structure_id)
        residues = self.protein_residues(structure)
        if not residues:
            raise ValueError(f"No protein residues in the biological assembly of {self.structure_id}")
        if len(residues) > MAX_RESIDUE_LABELS:
            raise ValueError(f"{len(residues)} protein residues in {self.structure_id}, "
                             f"more than the {MAX_RESIDUE_LABELS} unique labels freesasa can read")
        labels = residue_labels(len(residues))
        freesasa_residues = run_freesasa(self.write_labelled_model(structure, residues, labels))
        self.validate(residues, labels, freesasa_residues)
        records = [self.make_record(bio_chain, residue, freesasa_residues[label_key(label)])
                   for label, (bio_chain, residue) in zip(labels, residues)]
        sasa_df = pd.DataFrame(records)
        self.add_author_chains_and_uniprot(sasa_df, structure, residues)
        float_columns = ["asa", "rel_asa", "asa_mainchain", "asa_sidechain",
                         "rel_asa_mainchain", "rel_asa_sidechain", "asa_nz"]
        sasa_df[float_columns] = sasa_df[float_columns].astype(float)
        sasa_df["uniprot_id"] = sasa_df["uniprot_id"].astype(object).where(sasa_df["uniprot_id"].notna(), None)
        return sasa_df

    @staticmethod
    def protein_residues(structure: Structure) -> list[tuple[str, gemmi.Residue]]:
        # residues of the biological assembly that belong to protein entities, as (bio_chain, residue)
        peptide_entities = {e.name for e in structure.gemmi_structure.entities
                            if e.entity_type == gemmi.EntityType.Polymer and e.polymer_type in PEPTIDE_POLYMER_TYPES}
        model = structure.gemmi_bio_model.clone()
        model.remove_alternative_conformations()
        model.remove_hydrogens()
        return [(chain.name, residue.clone()) for chain in model for residue in chain
                if residue.entity_id in peptide_entities]

    @staticmethod
    def write_labelled_model(structure: Structure, residues: list[tuple[str, gemmi.Residue]],
                             labels: list[gemmi.SeqId]) -> gemmi.Structure:
        model = gemmi.Model("1")
        chains = {}
        for label, (bio_chain, residue) in zip(labels, residues):
            if bio_chain not in chains:
                model.add_chain(gemmi.Chain(bio_chain))
                chains[bio_chain] = model[len(model) - 1]
            labelled = residue.clone()
            labelled.seqid = label
            labelled.het_flag = 'A'   # written as ATOM so freesasa keeps modified amino acids (e.g. MSE)
            chains[bio_chain].add_residue(labelled)
        labelled_structure = gemmi.Structure()
        labelled_structure.add_model(model)
        labelled_structure.name = structure.gemmi_structure.name
        labelled_structure.cell = structure.gemmi_structure.cell
        labelled_structure.spacegroup_hm = structure.gemmi_structure.spacegroup_hm
        labelled_structure.setup_entities()
        return labelled_structure

    def validate(self, residues: list[tuple[str, gemmi.Residue]], labels: list[gemmi.SeqId],
                 freesasa_residues: dict[str, dict]):
        sent = {label_key(label): residue.name for label, (_, residue) in zip(labels, residues)}
        received = {label: res["name"].strip() for label, res in freesasa_residues.items()}
        if sent != received:
            missing = len(set(sent) - set(received))
            extra = len(set(received) - set(sent))
            raise ValueError(f"freesasa residues do not match the residues sent for {self.structure_id}: "
                             f"{missing} missing, {extra} unexpected, "
                             f"{sum(sent[k] != received[k] for k in set(sent) & set(received))} with a different name")

    def make_record(self, bio_chain: str, residue: gemmi.Residue, freesasa_residue: dict) -> dict:
        area = freesasa_residue["area"]
        relative_area = freesasa_residue.get("relative-area") or {}
        nz_area = next((atom["area"] for atom in freesasa_residue["atoms"] if atom["name"].strip() == "NZ"), None)
        return {
            'structure_id':      self.structure_id,
            'bio_chain':         bio_chain,
            'author_chain':      None,   # filled in add_author_chains_and_uniprot
            'resnum':            residue.seqid.num,
            'icode':             residue.seqid.icode.strip(),
            'resname':           seq1(residue.name),
            'resname3':          residue.name,
            'asa':               round(area["total"], 3),
            'rel_asa':           relative_fraction(relative_area.get("total")),
            'asa_mainchain':     round(area["main-chain"], 3),
            'asa_sidechain':     round(area["side-chain"], 3),
            'rel_asa_mainchain': relative_fraction(relative_area.get("main-chain")),
            'rel_asa_sidechain': relative_fraction(relative_area.get("side-chain")),
            'asa_nz':            round(nz_area, 3) if nz_area is not None and residue.name == "LYS" else None,
            'uniprot_id':        None,   # filled in add_author_chains_and_uniprot
            'unp_resnum':        None,
            'unp_resname':       None,
        }

    def add_author_chains_and_uniprot(self, sasa_df: pd.DataFrame, structure: Structure,
                                      residues: list[tuple[str, gemmi.Residue]]):
        bio_to_author = {bio_chain: structure.biochain_to_chain(bio_chain) for bio_chain in sasa_df["bio_chain"].unique()}
        sasa_df["author_chain"] = sasa_df["bio_chain"].map(bio_to_author)
        try:
            entity_to_uniprot = defaultdict(list)
            for uniprot_id, entities in structure.uniprot_to_entities.items():
                for entity in entities:
                    entity_to_uniprot[entity].append(uniprot_id)
            # copies of a chain share residue numbering, so map each (author chain, residue) once
            cache = {}
            unp = []
            for (bio_chain, residue) in residues:
                author_chain = bio_to_author[bio_chain]
                key = (author_chain, str(residue.seqid))
                if key not in cache:
                    cache[key] = self.map_to_uniprot(structure, author_chain, residue, entity_to_uniprot[residue.entity_id])
                unp.append(cache[key])
            sasa_df["uniprot_id"] = pd.Series([u[0] for u in unp], dtype=object)
            sasa_df["unp_resnum"] = pd.array([u[1] for u in unp], dtype="Int64")
            sasa_df["unp_resname"] = pd.Series([u[2] for u in unp], dtype=object)
        except Exception as e:
            # keep the SASA values even if the uniprot alignment fails; uniprot columns stay empty
            warnings.warn(f"UniProt mapping failed for {self.structure_id}: {e}")

    @staticmethod
    def map_to_uniprot(structure: Structure, author_chain: str, residue: gemmi.Residue,
                       uniprot_ids: list[str]) -> tuple[str, int, str]:
        for uniprot_id in uniprot_ids:
            unp_residue = structure.map_residue_to_uniprot(author_chain, residue, uniprot_id)
            if unp_residue is not None:
                return uniprot_id, unp_residue[0], unp_residue[1]
        return None, None, None

    def get_sasa_for_residue(self, uniprot_id: str, residues: tuple[int, str]|list[tuple[int, str]]|list[int]|int,
                             type="asa_nz") -> list[list[float]]:
        """
        For each uniprot residue, the list of values from every chain copy in the assembly
        that maps to it (empty list if the residue is not in the structure).
        """
        if not isinstance(residues, list):
            residues = [residues]
        df = self.sasa_df[self.sasa_df["uniprot_id"] == uniprot_id]
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
        return [u for u in self.sasa_df["uniprot_id"].dropna().unique()]

    def summary(self) -> str:
        n_nz = self.sasa_df["asa_nz"].notna().sum()
        return (f"{self.structure_id}: {self.sasa_df['bio_chain'].nunique()} chains, "
                f"{len(self.sasa_df)} residues, {n_nz} lysine NZ, uniprot ids {self.uniprot_ids}")

    def __str__(self):
        return f"SASAMultimer-{self.structure_id}"


class SASAMultimerStore(CompressedPickleStore):
    def __init__(self):
        super().__init__(directory=SASA_MULTIMER_STORE_DIR, extension="sasa", object_type=SASAMultimer)


def residue_labels(n: int) -> list[gemmi.SeqId]:
    # n unique residue labels that freesasa reads back unchanged (see PLAIN_LABEL_MAX)
    labels = [gemmi.SeqId(num, ' ') for num in range(1, min(n, PLAIN_LABEL_MAX) + 1)]
    for icode in INSERTION_CODES:
        if len(labels) >= n:
            break
        labels.extend(gemmi.SeqId(num, icode) for num in LABEL_NUMBERS_WITH_INSERTION_CODE)
    return labels[:n]


def label_key(label: gemmi.SeqId) -> str:
    # the residue "number" freesasa reports for this label
    return f"{label.num}{label.icode.strip()}"


def relative_fraction(percent: float|None) -> float|None:
    # freesasa json reports relative areas in percent; stored as a fraction to match SASAMonomer
    if percent is None or (isinstance(percent, float) and np.isnan(percent)):
        return None
    return round(percent / 100, 4)


def freesasa_executable() -> str:
    # prefer the freesasa installed next to the running python (conda env), then PATH
    candidate = Path(sys.executable).parent / "freesasa"
    if candidate.exists():
        return str(candidate)
    found = shutil.which("freesasa")
    if found is None:
        raise FileNotFoundError("freesasa command line tool not found; install with "
                                "`conda install -c conda-forge freesasa-c`")
    return found


def run_freesasa(structure: gemmi.Structure) -> dict[str, dict]:
    """
    Runs `freesasa --cif --format=json --depth=atom` on the structure and returns
    the freesasa residues keyed by residue number (the unique label).
    """
    groups = gemmi.MmcifOutputGroups(True)
    groups.auth_all = True   # freesasa 2.1.2 segfaults on mmCIF without auth_comp_id / auth_atom_id
    with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp:
        cif_path = tmp.name
    try:
        structure.make_mmcif_document(groups).write_file(cif_path)
        result = subprocess.run([freesasa_executable(), "--cif", "--radii=protor",
                                 "--format=json", "--depth=atom", cif_path],
                                capture_output=True, text=True)
    finally:
        os.remove(cif_path)
    if result.returncode != 0:
        raise RuntimeError(f"freesasa failed with return code {result.returncode}: {result.stderr[:500]}")
    chains = json.loads(result.stdout)["results"][0]["structure"][0]["chains"]
    freesasa_residues = {}
    for chain in chains:
        for residue in chain["residues"]:
            label = residue["number"].strip()
            if label in freesasa_residues:
                raise ValueError(f"Residue label {label} returned more than once by freesasa")
            freesasa_residues[label] = residue
    return freesasa_residues
