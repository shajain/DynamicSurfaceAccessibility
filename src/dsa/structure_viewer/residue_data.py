"""
Collects, for every residue of a structure shown in the viewer, the properties displayed on hover
and the lysine contacts / bonds highlighted on hover.

PDB structures are shown as their biological assembly (bio chains A1, B1, A2, ...), the same model the
multimer SASA and the contacts are computed on. AlphaFold structures are a single chain numbered by UniProt.

SASA columns, each with asa (A^2), rel_asa (fraction) and asa_nz (A^2, lysine NZ atom):
    multimer   PDB biological assembly, from SASAMultimerStore (burial by other chains included)
    chain      PDB chain on its own, computed here with the same freesasa pipeline (one run per author chain)
    alphafold  AlphaFold monomer of the residue's UniProt id, from SASAMonomerStore, matched by UniProt number
"""
from collections import defaultdict
from dataclasses import dataclass, field
import warnings
import gemmi
import pandas as pd
from Bio.SeqUtils import seq1
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure.structure import Structure
from dsa.binding_interfaces.contacts.structure_contacts import StructureContactsStore, StructureContacts
from dsa.binding_interfaces.contacts.helpers import is_lysine, ALL_LYSINE_NAMES
from dsa.sasa.sasa_multimer import SASAMultimer, SASAMultimerStore, residue_labels, label_key, run_freesasa
from dsa.sasa.sasa_monomer import SASAMonomerStore

BOND_TYPES = StructureContacts.BOND_TYPES    # salt_bridge, hydrogen_bond, covalent_isopeptide, cation_pi
CLOSE_CONTACT_DISTANCE = 4.5                 # contacts this close are shown as "contact <= 4.5 A" (the dist4.5 class)
SASA_COLUMNS = ["asa", "rel_asa", "asa_nz"]
ResidueKey = tuple[str, int, str]            # (bio_chain, seqid number, insertion code)


@dataclass
class StructureData:
    """Everything the viewer needs: the model to draw, one row per residue and the residue pairs in contact."""
    title: str
    database: str
    structure_id: str
    model: gemmi.Structure                 # what is drawn (waters, hydrogens and alternate conformations removed)
    residues: pd.DataFrame                 # one row per residue, in model order; index = residue index used by the page
    atom_residue: list[int]                # residue index of every atom, in model order (atom serial = position + 1)
    # residue pairs: res1, res2, types, dist, atom1, atom2 (serials), source ("pdb" / "alphafold"), scope ("inter" / "intra")
    bonds: list[dict] = field(default_factory=list)


def residue_key(chain: gemmi.Chain|str, residue: gemmi.Residue) -> ResidueKey:
    chain_name = chain if isinstance(chain, str) else chain.name
    return chain_name, residue.seqid.num, residue.seqid.icode.strip()


def load_pdb(structure_id: str, chains: list[str]|None = None) -> StructureData:
    contacts = StructureContactsStore("pdb").load(structure_id)
    structure = contacts.structure if contacts is not None else ClassFactory.create(Structure, "pdb", structure_id)
    if contacts is None:
        warnings.warn(f"No saved contacts for {structure_id}; bonds will not be shown")
    model, residues, atom_residue, serials = build_model(structure, chains)
    add_multimer_sasa(residues, structure_id)
    add_chain_sasa(residues, structure, structure_id)
    add_alphafold_sasa(residues)
    bonds = residue_bonds(contacts, residues, serials, source="pdb") if contacts is not None else []
    bonds += alphafold_bonds_on_pdb(residues, serials)
    return StructureData(f"PDB {structure_id}", "pdb", structure_id, model, residues, atom_residue, bonds)


def load_alphafold(uniprot_id: str) -> StructureData:
    store = StructureContactsStore("alphafold")
    contacts = store.load(uniprot_id) if store.get_file_path(uniprot_id).exists() else None
    structure = contacts.structure if contacts is not None else ClassFactory.create(Structure, "alphafold", uniprot_id)
    if contacts is None:
        warnings.warn(f"No saved contacts for AlphaFold {uniprot_id}; bonds will not be shown")
    model, residues, atom_residue, serials = build_model(structure, None)
    is_protein = residues["kind"] == "protein"
    # AlphaFold models are numbered by UniProt
    residues.loc[is_protein, "uniprot_id"] = uniprot_id
    residues.loc[is_protein, "unp_resnum"] = residues.loc[is_protein, "resnum"]
    residues.loc[is_protein, "unp_resname"] = residues.loc[is_protein, "resname1"]
    for prefix in ("multimer", "chain"):
        for column in SASA_COLUMNS:
            residues[f"{prefix}_{column}"] = None
    add_alphafold_sasa(residues)
    bonds = residue_bonds(contacts, residues, serials, source="alphafold") if contacts is not None else []
    return StructureData(f"AlphaFold {uniprot_id}", "alphafold", uniprot_id, model, residues, atom_residue, bonds)


def build_model(structure: Structure, chains: list[str]|None):
    """The model to draw, its residue table, the residue of every atom and the serial of every atom."""
    model = structure.gemmi_bio_model.clone()
    if chains is not None:
        missing = set(chains) - {chain.name for chain in model}
        if missing:
            raise ValueError(f"Chains {sorted(missing)} not in the biological assembly; "
                             f"available: {[chain.name for chain in model]}")
        for name in [chain.name for chain in model if chain.name not in chains]:
            model.remove_chain(name)
    drawn = gemmi.Structure()
    drawn.add_model(model)
    drawn.remove_waters()
    drawn.remove_hydrogens()
    drawn.remove_alternative_conformations()
    drawn.remove_empty_chains()

    peptide_entities = {e.name for e in structure.gemmi_structure.entities
                        if e.entity_type == gemmi.EntityType.Polymer
                        and e.polymer_type in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)}
    rows, atom_residue, serials = [], [], {}
    serial = 0
    for chain in drawn[0]:
        for residue in chain:
            index = len(rows)
            info = gemmi.find_tabulated_residue(residue.name)
            kind = ("protein" if residue.entity_id in peptide_entities else
                    "nucleic" if info is not None and info.is_nucleic_acid() else
                    "polymer" if residue.het_flag == "A" else "ligand")
            rows.append({"bio_chain": chain.name, "resnum": residue.seqid.num, "icode": residue.seqid.icode.strip(),
                         "resname3": residue.name, "resname1": one_letter(residue.name, info), "kind": kind,
                         # one-letter code K as in the contacts (LYS and modified lysines such as MLY, M3L, ALY, KCX),
                         # or a name in ALL_LYSINE_NAMES (LYZ and ALZ have no one-letter code in gemmi)
                         "is_lysine": residue.name in ALL_LYSINE_NAMES or (info is not None and is_lysine(residue))})
            for atom in residue:
                serial += 1
                atom.serial = serial
                atom_residue.append(index)
                serials[residue_key(chain, residue) + (atom.name,)] = serial
    drawn.setup_entities()   # after the residue kinds are read from the entity ids of the original structure
    residues = pd.DataFrame(rows)
    residues["author_chain"] = residues["bio_chain"].map(
        {bio_chain: author_chain(structure, bio_chain) for bio_chain in residues["bio_chain"].unique()})
    for column in ("uniprot_id", "unp_resnum", "unp_resname"):
        residues[column] = None
    return drawn, residues, atom_residue, serials


def one_letter(name: str, info) -> str:
    # gemmi's code (lower case for modified residues, e.g. M3L -> k) in upper case, else Biopython's, else X
    if info is not None and info.one_letter_code.strip():
        return info.one_letter_code.upper()
    return seq1(name, custom_map={}) or "X"


def author_chain(structure: Structure, bio_chain: str) -> str|None:
    try:
        return structure.biochain_to_chain(bio_chain)
    except Exception:
        return None


def add_multimer_sasa(residues: pd.DataFrame, structure_id: str):
    multimer = SASAMultimerStore().load(structure_id)
    if multimer is None:
        warnings.warn(f"No saved multimer SASA for {structure_id}; computing it")
        multimer = SASAMultimer(structure_id)
    sasa = multimer.sasa_df.rename(columns={c: f"multimer_{c}" for c in SASA_COLUMNS})
    keys = ["bio_chain", "resnum", "icode"]
    merged = residues[keys].merge(sasa[keys + ["uniprot_id", "unp_resnum", "unp_resname"] + [f"multimer_{c}" for c in SASA_COLUMNS]],
                                  on=keys, how="left")
    for column in ["uniprot_id", "unp_resnum", "unp_resname"] + [f"multimer_{c}" for c in SASA_COLUMNS]:
        residues[column] = merged[column].values


def add_chain_sasa(residues: pd.DataFrame, structure: Structure, structure_id: str):
    """SASA of each chain on its own. Copies of an author chain are identical alone, so each author chain is run once."""
    obj = SASAMultimer.__new__(SASAMultimer)   # only to reuse make_record; no SASA is computed here
    obj.structure_id = structure_id
    protein = SASAMultimer.protein_residues(structure)
    by_bio_chain = defaultdict(list)
    for bio_chain, residue in protein:
        by_bio_chain[bio_chain].append(residue)
    values = {}
    done_author_chains = {}
    for bio_chain, chain_residues in by_bio_chain.items():
        author = author_chain(structure, bio_chain)
        if author in done_author_chains:
            source = done_author_chains[author]
            for residue in chain_residues:
                key = residue_key(bio_chain, residue)
                values[key] = values.get(residue_key(source, residue))
            continue
        done_author_chains[author] = bio_chain
        pairs = [(bio_chain, residue) for residue in chain_residues]
        labels = residue_labels(len(pairs))
        freesasa_residues = run_freesasa(SASAMultimer.write_labelled_model(structure, pairs, labels))
        for label, (_, residue) in zip(labels, pairs):
            record = obj.make_record(bio_chain, residue, freesasa_residues[label_key(label)])
            values[residue_key(bio_chain, residue)] = {c: record[c] for c in SASA_COLUMNS}
    keys = list(zip(residues["bio_chain"], residues["resnum"], residues["icode"]))
    for column in SASA_COLUMNS:
        residues[f"chain_{column}"] = [(values.get(k) or {}).get(column) for k in keys]


def add_alphafold_sasa(residues: pd.DataFrame):
    """AlphaFold monomer SASA of the residue's UniProt id at the same UniProt number (None if the residue type differs)."""
    store = SASAMonomerStore()
    tables = {}
    for uniprot_id in residues["uniprot_id"].dropna().unique():
        if not store.get_file_path(uniprot_id).exists():
            continue
        monomer = store.load(uniprot_id)
        if monomer is not None:
            df = monomer.sasa_df
            tables[uniprot_id] = df[~df.index.duplicated()]
    out = {c: [] for c in SASA_COLUMNS}
    for uniprot_id, unp_resnum, unp_resname in zip(residues["uniprot_id"], residues["unp_resnum"], residues["unp_resname"]):
        df = tables.get(uniprot_id)
        row = None
        if df is not None and pd.notna(unp_resnum) and int(unp_resnum) in df.index:
            row = df.loc[int(unp_resnum)]
            if unp_resname is not None and row["resname"] != unp_resname:
                row = None
        for column in SASA_COLUMNS:
            out[column].append(None if row is None or pd.isna(row[column]) else float(row[column]))
    for column in SASA_COLUMNS:
        residues[f"alphafold_{column}"] = out[column]


def group_contacts(contacts: StructureContacts) -> dict[tuple[ResidueKey, ResidueKey], dict]:
    """
    The atom contacts grouped by residue pair (lysine NZ within the contact cutoff of the partner): the bond types
    satisfied, the closest atom pair (names) and its distance, and the one-letter codes of the two residues.
    """
    pairs = {}
    for contact in contacts.contacts:
        side1 = (residue_key(contact.chain1, contact.residue1), contact.atom1.name, contact.residue1.name)
        side2 = (residue_key(contact.chain2, contact.residue2), contact.atom2.name, contact.residue2.name)
        first, second = sorted([side1, side2], key=lambda side: side[0])
        entry = pairs.setdefault((first[0], second[0]), {"types": set(), "dist": None, "atom1": None, "atom2": None,
                                                         "name1": first[2], "name2": second[2]})
        entry["types"].update(bt for bt in BOND_TYPES if contact.satisfies_bond_type(bt))
        if entry["dist"] is None or contact.dist < entry["dist"]:
            entry.update(dist=round(contact.dist, 2), atom1=first[1], atom2=second[1])
    return pairs


def bond_types(types: set[str], dist: float) -> list[str]:
    types = sorted(types)
    if dist <= CLOSE_CONTACT_DISTANCE:
        types.append("contact_4.5")
    types.append("contact_6.5")
    return types


def residue_bonds(contacts: StructureContacts, residues: pd.DataFrame, serials: dict, source: str) -> list[dict]:
    """One entry per pair of drawn residues in contact in this structure; pairs with a residue not drawn are skipped."""
    index = {key: i for i, key in enumerate(zip(residues["bio_chain"], residues["resnum"], residues["icode"]))}
    bonds = []
    for (key1, key2), entry in group_contacts(contacts).items():
        if key1 not in index or key2 not in index:
            continue
        bonds.append({"res1": index[key1], "res2": index[key2], "types": bond_types(entry["types"], entry["dist"]),
                      "dist": entry["dist"], "atom1": serials.get(key1 + (entry["atom1"],)),
                      "atom2": serials.get(key2 + (entry["atom2"],)), "source": source,
                      "scope": "intra" if key1[0] == key2[0] else "inter"})
    return bonds


def alphafold_bonds_on_pdb(residues: pd.DataFrame, serials: dict) -> list[dict]:
    """
    The intra-chain contacts of the AlphaFold model of every UniProt id in the structure, placed on each PDB chain copy
    of that protein by UniProt number. A pair is placed only if both residues are in the chain copy and have the same
    one-letter code as in the AlphaFold model. The distance is the AlphaFold distance.
    """
    store = StructureContactsStore("alphafold")
    protein = residues[residues["uniprot_id"].notna() & residues["unp_resnum"].notna()]
    # (bio_chain, uniprot_id, unp_resnum) -> residue index
    position = {(c, u, int(n)): i for i, c, u, n in zip(protein.index, protein["bio_chain"], protein["uniprot_id"], protein["unp_resnum"])}
    bonds = []
    for uniprot_id in protein["uniprot_id"].unique():
        if not store.get_file_path(uniprot_id).exists():
            continue
        contacts = store.load(uniprot_id)
        if contacts is None:
            continue
        pairs = group_contacts(contacts)
        bio_chains = protein.loc[protein["uniprot_id"] == uniprot_id, "bio_chain"].unique()
        for (key1, key2), entry in pairs.items():
            letter1, letter2 = one_letter(entry["name1"], gemmi.find_tabulated_residue(entry["name1"])), \
                               one_letter(entry["name2"], gemmi.find_tabulated_residue(entry["name2"]))
            for bio_chain in bio_chains:
                i1, i2 = position.get((bio_chain, uniprot_id, key1[1])), position.get((bio_chain, uniprot_id, key2[1]))
                if i1 is None or i2 is None or residues.at[i1, "resname1"] != letter1 or residues.at[i2, "resname1"] != letter2:
                    continue
                pdb1 = (bio_chain, int(residues.at[i1, "resnum"]), residues.at[i1, "icode"])
                pdb2 = (bio_chain, int(residues.at[i2, "resnum"]), residues.at[i2, "icode"])
                bonds.append({"res1": i1, "res2": i2, "types": bond_types(entry["types"], entry["dist"]),
                              "dist": entry["dist"], "atom1": serials.get(pdb1 + (entry["atom1"],)),
                              "atom2": serials.get(pdb2 + (entry["atom2"],)), "source": "alphafold", "scope": "intra"})
    return bonds
