"""
Lookup table of structural properties of uniprot lysines.

Every property is a property of a single lysine, keyed by (uniprot_id, unp_resnum), so the table
is computed once for all lysines of interest and then looked up by MSDSAAnnotator for every
peptide and run. Each data source fills its own columns and leaves them empty (NA) when it has
no data for a protein, without affecting the other sources.
"""
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import pandas as pd
from dsa.binding_interfaces.contacts.protein_contacts import ProteinContactsStore
from dsa.binding_interfaces.contacts.bond_characterization import bond_characteristics  # registers the bond types
from dsa.binding_interfaces.contacts.bond_characterization.library import BondLibrary
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappingsStore
from dsa.sasa.sasa_monomer import SASAMonomerStore
from dsa.sasa.sasa_multimer import SASAMultimerStore
from dsa.bfactor.lysine_bfactor import LysineBFactorStore
from dsa.ms_dsa.evaluation1.config import LYSINE_PROPERTIES_FILE


NO_CONTACT_DISTANCE = 4.5
CONTACT_CLASSES = [f"dist{NO_CONTACT_DISTANCE:g}", *BondLibrary.bond_types_implemented(), "any_bond"]
CONTACT_TYPES = ["inter", "intra"]


class LysinePropertyTable:
    """
    One row per (uniprot_id, unp_resnum) with the columns:

    contacts (True if the lysine has the contact in any structure, NA if no contacts are saved):
        intra_chain_pdb, inter_chain_pdb, intra_chain_alphafold,
        {bond_type}_pdb, {bond_type}_alphafold
    AlphaFold monomer SASA (NA if the residue is missing or is not a lysine):
        asa_nz, rel_asa
    PDB biological assembly SASA, over every chain copy of every structure containing the lysine:
        asa_nz_multimer_min, asa_nz_multimer_mean, asa_nz_multimer_max, n_multimer
    PDB B-factor z-scores (normalized per chain), over every chain of every structure:
        nz_bfactor_z_min, nz_bfactor_z_mean, nz_bfactor_z_max, avg_bfactor_z_mean, n_bfactor
    contact classes (CONTACT_CLASSES: dist4.5, each bond type, any_bond), computed per chain copy:
        pdb_{inter|intra}_{class}                  True if some PDB chain copy of the lysine has the contact
        asa_nz_pdb_{inter|intra}_{class}_{min,mean,max}, n_pdb_{inter|intra}_{class}
                                                   multimer asa_nz over only the chain copies with the contact
        af_intra_{class}                           True if the AlphaFold monomer has the contact
        no_contact_pdb, asa_nz_pdb_no_contact_{min,mean,max}
                                                   no contact within NO_CONTACT_DISTANCE in any PDB chain copy;
                                                   multimer asa_nz over all its copies
        no_contact_af                              no contact within NO_CONTACT_DISTANCE in the AlphaFold monomer
        purely_intra_pdb, purely_intra_af          an intra-chain class and not pdb_inter_any_bond
    PDB and AlphaFold values are never combined into one number.

    Build with LysinePropertyTable.build({uniprot_id: [lysine positions]}), or load a saved table
    with LysinePropertyTable.load().
    """
    KEY = ["uniprot_id", "unp_resnum"]

    def __init__(self, table: pd.DataFrame):
        self.table = table

    @classmethod
    def build(cls, lysines: dict[str, list[int]], distance_cutoff: float = 6.5, contact_filter_key: str = "lysine_only",
              n_workers: int = 1, show_progress: bool = True) -> "LysinePropertyTable":
        # proteins are independent, so they are processed in parallel with n_workers > 1
        lysines = {uniprot_id: sorted(set(resnums)) for uniprot_id, resnums in lysines.items() if resnums}
        builder_args = (distance_cutoff, contact_filter_key)
        tables = {}
        if n_workers > 1:
            with Pool(n_workers, initializer=_init_worker_builder, initargs=builder_args) as pool:
                results = pool.imap_unordered(_worker_properties_of_protein, lysines.items())
                for i, (uniprot_id, table) in enumerate(results, start=1):
                    tables[uniprot_id] = table
                    report_progress(i, len(lysines), show_progress)
        else:
            builder = _LysinePropertyBuilder(*builder_args)
            for i, (uniprot_id, resnums) in enumerate(lysines.items(), start=1):
                tables[uniprot_id] = builder.properties_of_protein(uniprot_id, resnums)
                report_progress(i, len(lysines), show_progress)
        if not tables:
            return cls(pd.DataFrame(columns=cls.KEY))
        return cls(pd.concat([tables[uniprot_id] for uniprot_id in lysines], ignore_index=True))

    def extend(self, lysines: dict[str, list[int]], **build_kwargs) -> "LysinePropertyTable":
        # computes only the lysines that are not already in the table
        existing = set(zip(self.table["uniprot_id"], self.table["unp_resnum"]))
        missing = {uniprot_id: [r for r in resnums if (uniprot_id, r) not in existing]
                   for uniprot_id, resnums in lysines.items()}
        new = self.build(missing, **build_kwargs).table
        if new.empty:
            return self
        return LysinePropertyTable(pd.concat([self.table, new], ignore_index=True))

    def save(self, path: Path = LYSINE_PROPERTIES_FILE):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.table.to_parquet(path, index=False)

    @classmethod
    def load(cls, path: Path = LYSINE_PROPERTIES_FILE) -> "LysinePropertyTable":
        return cls(pd.read_parquet(path))

    @property
    def property_columns(self) -> list[str]:
        return [c for c in self.table.columns if c not in self.KEY]

    def __len__(self):
        return len(self.table)


class _LysinePropertyBuilder:
    """Computes the properties of the lysines of one protein; one method per data source."""
    STRUCTURE_DATABASE = "pdb"

    def __init__(self, distance_cutoff: float, contact_filter_key: str):
        self.contacts_stores = {database: ProteinContactsStore(database, distance_cutoff, contact_filter_key)
                                for database in ("pdb", "alphafold")}
        self.sasa_monomer_store = SASAMonomerStore()
        self.sasa_multimer_store = SASAMultimerStore()
        self.bfactor_store = LysineBFactorStore()
        self.pdb_mapper = UniprotStructureMappingsStore.get_mapper(self.STRUCTURE_DATABASE)
        self.saved_multimers = set(self.sasa_multimer_store.get_saved_keys())
        self.saved_bfactors = set(self.bfactor_store.get_saved_keys())

    def properties_of_protein(self, uniprot_id: str, resnums: list[int]) -> pd.DataFrame:
        # contacts and structures loaded for this protein are shared by the sources below
        self._contacts, self._structures = {}, {}
        table = pd.DataFrame({"uniprot_id": uniprot_id, "unp_resnum": pd.array(resnums, dtype="Int64")})
        for source in (self.contact_properties, self.monomer_sasa_properties,
                       self.multimer_sasa_properties, self.bfactor_properties, self.contact_class_properties):
            for column, values in source(uniprot_id, resnums).items():
                table[column] = values
        return table

    def contact_properties(self, uniprot_id: str, resnums: list[int]) -> dict[str, pd.array]:
        residues = [(resnum, "K") for resnum in resnums]
        flags = {}
        for database in self.contacts_stores:
            contacts = self.load_contacts(database, uniprot_id)
            queries = {f"{contact_type}_{database}": {"contact_type": contact_type}
                       for contact_type in (["intra_chain", "inter_chain"] if database == "pdb" else ["intra_chain"])}
            queries.update({f"{bond_type}_{database}": {"bond_type": bond_type, "method": BondLibrary.default_method}
                            for bond_type in BondLibrary.bond_types_implemented()})
            for column, query in queries.items():
                values = contacts.does_residue_satisfy(residues, **query) if contacts is not None else [None] * len(resnums)
                flags[column] = pd.array(values, dtype="boolean")
        return flags

    def monomer_sasa_properties(self, uniprot_id: str, resnums: list[int]) -> dict[str, list]:
        sasa_monomer = self.sasa_monomer_store.load(uniprot_id)
        if sasa_monomer is None:
            return {"asa_nz": [np.nan] * len(resnums), "rel_asa": [np.nan] * len(resnums)}
        lysines = sasa_monomer.sasa_df[sasa_monomer.sasa_df["resname"] == "K"]
        lysines = lysines[~lysines.index.duplicated()]
        return {column: lysines[column].reindex(resnums).astype(float).tolist() for column in ("asa_nz", "rel_asa")}

    def multimer_sasa_properties(self, uniprot_id: str, resnums: list[int]) -> dict[str, list]:
        values = self.values_over_structures(self.sasa_multimer_store, self.saved_multimers, "sasa_df",
                                             uniprot_id, ["asa_nz"])
        stats = values.groupby("unp_resnum")["asa_nz"].agg(["min", "mean", "max", "count"]).reindex(resnums)
        return {"asa_nz_multimer_min": stats["min"].tolist(),
                "asa_nz_multimer_mean": stats["mean"].tolist(),
                "asa_nz_multimer_max": stats["max"].tolist(),
                "n_multimer": stats["count"].fillna(0).astype(int).tolist()}

    def bfactor_properties(self, uniprot_id: str, resnums: list[int]) -> dict[str, list]:
        values = self.values_over_structures(self.bfactor_store, self.saved_bfactors, "bfactor_df",
                                             uniprot_id, ["nz_bfactor_z", "avg_bfactor_z"])
        grouped = values.groupby("unp_resnum")
        nz = grouped["nz_bfactor_z"].agg(["min", "mean", "max"]).reindex(resnums)
        return {"nz_bfactor_z_min": nz["min"].tolist(),
                "nz_bfactor_z_mean": nz["mean"].tolist(),
                "nz_bfactor_z_max": nz["max"].tolist(),
                "avg_bfactor_z_mean": grouped["avg_bfactor_z"].mean().reindex(resnums).tolist(),
                "n_bfactor": grouped["avg_bfactor_z"].count().reindex(resnums).fillna(0).astype(int).tolist()}

    def values_over_structures(self, store, saved_keys: set[str], df_attribute: str,
                               uniprot_id: str, columns: list[str]) -> pd.DataFrame:
        # rows of the lysines of uniprot_id from every saved structure mapped to it
        frames = []
        for structure_id in self.pdb_mapper.structures_for_uniprot_id(uniprot_id):
            if structure_id not in saved_keys:
                continue
            obj = self.load_structure_object(store, structure_id)
            if obj is None:
                continue
            df = getattr(obj, df_attribute)
            df = df[(df["uniprot_id"] == uniprot_id) & (df["resname"] == "K") & (df["unp_resname"] == "K")]
            frames.append(df[["unp_resnum"] + columns])
        if not frames:
            return pd.DataFrame({"unp_resnum": pd.Series(dtype="Int64"),
                                 **{c: pd.Series(dtype=float) for c in columns}})
        return pd.concat(frames, ignore_index=True)

    def load_contacts(self, database: str, uniprot_id: str):
        if database not in self._contacts:
            self._contacts[database] = self.contacts_stores[database].load(uniprot_id)
        return self._contacts[database]

    def load_structure_object(self, store, structure_id: str):
        key = (id(store), structure_id)
        if key not in self._structures:
            self._structures[key] = store.load(structure_id)
        return self._structures[key]

    def contact_class_properties(self, uniprot_id: str, resnums: list[int]) -> dict[str, list]:
        columns = {}
        pdb_contacts = self.load_contacts("pdb", uniprot_id)
        af_contacts = self.load_contacts("alphafold", uniprot_id)
        columns.update(self.pdb_contact_class_properties(uniprot_id, resnums, pdb_contacts))
        columns.update(self.af_contact_class_properties(uniprot_id, resnums, af_contacts))
        pdb_inter = pd.array(columns["pdb_inter_any_bond"], dtype="boolean")
        for source in ("pdb", "af"):
            intra = [f"{source}_intra_{c}" for c in CONTACT_CLASSES]
            any_intra = pd.DataFrame({c: pd.array(columns[c], dtype="boolean") for c in intra}).any(axis=1, skipna=True)
            known = pd.Series(pd.array(columns[intra[0]], dtype="boolean")).notna() & pd.Series(pdb_inter).notna()
            columns[f"purely_intra_{source}"] = pd.array((any_intra & ~pd.Series(pdb_inter).fillna(False))
                                                         .where(known, pd.NA), dtype="boolean")
        flags = ("pdb_inter_", "pdb_intra_", "af_intra_", "no_contact_")
        return {name: pd.array(values, dtype="boolean") if name.startswith(flags) else values
                for name, values in columns.items()}

    def pdb_contact_class_properties(self, uniprot_id: str, resnums: list[int], protein_contacts) -> dict[str, list]:
        # chain copies of the lysines, keyed (structure_id, bio_chain, resnum, icode), with their multimer asa_nz
        copies_asa = {}
        for structure_id in self.pdb_mapper.structures_for_uniprot_id(uniprot_id):
            if structure_id not in self.saved_multimers:
                continue
            sasa = self.load_structure_object(self.sasa_multimer_store, structure_id)
            if sasa is None:
                continue
            df = sasa.sasa_df
            df = df[(df["uniprot_id"] == uniprot_id) & (df["resname3"] == "LYS") & df["asa_nz"].notna()]
            for row in df.itertuples():
                copies_asa[(structure_id, row.bio_chain, int(row.resnum), row.icode or "")] = (int(row.unp_resnum), row.asa_nz)
        observed = {unp_resnum for unp_resnum, _ in copies_asa.values()}
        classified = self.lysine_copies_by_contact_class(protein_contacts) if protein_contacts is not None else None

        columns = {}
        for contact_type in CONTACT_TYPES:
            for contact_class in CONTACT_CLASSES:
                name = f"pdb_{contact_type}_{contact_class}"
                if classified is None:
                    columns[name] = [None] * len(resnums)
                    values = {}
                else:
                    copies = classified.get((contact_type, contact_class), set())
                    lysines = {unp_resnum for unp_resnum, _ in copies}
                    columns[name] = [r in lysines for r in resnums]
                    values = {}
                    for unp_resnum, copy in copies:
                        found = copies_asa.get(copy)
                        if found is not None:
                            values.setdefault(unp_resnum, []).append(found[1])
                self.add_statistics(columns, f"asa_nz_{name}", f"n_{name}", values, resnums)

        in_contact = set() if classified is None else {
            unp_resnum for contact_type in CONTACT_TYPES
            for unp_resnum, _ in classified.get((contact_type, CONTACT_CLASSES[0]), set())}
        no_contact = [None if classified is None or r not in observed else r not in in_contact for r in resnums]
        columns["no_contact_pdb"] = no_contact
        values = {}
        for unp_resnum, asa_nz in copies_asa.values():
            values.setdefault(unp_resnum, []).append(asa_nz)
        values = {r: v for r, v, flag in ((r, values.get(r), f) for r, f in zip(resnums, no_contact)) if flag and v}
        self.add_statistics(columns, "asa_nz_pdb_no_contact", "n_pdb_no_contact", values, resnums)
        return columns

    def af_contact_class_properties(self, uniprot_id: str, resnums: list[int], protein_contacts) -> dict[str, list]:
        columns = {}
        if protein_contacts is None:
            for contact_class in CONTACT_CLASSES:
                columns[f"af_intra_{contact_class}"] = [None] * len(resnums)
            columns["no_contact_af"] = [None] * len(resnums)
            return columns
        classified = self.lysine_copies_by_contact_class(protein_contacts)
        for contact_class in CONTACT_CLASSES:
            # the AlphaFold monomer has a single chain, so every contact is intra-chain
            lysines = {unp_resnum for unp_resnum, _ in classified.get(("intra", contact_class), set())}
            columns[f"af_intra_{contact_class}"] = [r in lysines for r in resnums]
        in_contact = {unp_resnum for unp_resnum, _ in classified.get(("intra", CONTACT_CLASSES[0]), set())}
        sasa_monomer = self.sasa_monomer_store.load(uniprot_id)
        observed = set()
        if sasa_monomer is not None:
            df = sasa_monomer.sasa_df
            observed = set(df.loc[(df["resname"] == "K") & df["asa_nz"].notna(), "resnum"].astype(int))
        columns["no_contact_af"] = [None if r not in observed else r not in in_contact for r in resnums]
        return columns

    @staticmethod
    def lysine_copies_by_contact_class(protein_contacts) -> dict[tuple[str, str], set[tuple[int, tuple]]]:
        """
        {(contact_type, contact_class): {(unp_resnum, (structure_id, bio_chain, resnum, icode))}} for the lysine
        sides of the contacts, a lysine side counting when it is a lysine mapped to the protein (as in
        UniprotAlignedStructureContacts.get_residues_satisfying).
        """
        classified = {}
        for structure_id, structure_contacts in protein_contacts.structure_contacts.items():
            for contact, unp1, unp2 in zip(structure_contacts.contacts, structure_contacts.unp_residue1s,
                                           structure_contacts.unp_residue2s):
                classes = contact_classes(contact)
                if not classes:
                    continue
                contact_type = "inter" if contact.is_inter_chain() else "intra"
                for chain, residue, unp, is_lysine in ((contact.chain1, contact.residue1, unp1, contact.is_residue1_lysine()),
                                                       (contact.chain2, contact.residue2, unp2, contact.is_residue2_lysine())):
                    if not (unp and is_lysine):
                        continue
                    copy = (structure_id, chain.name, residue.seqid.num, residue.seqid.icode.strip())
                    for contact_class in classes:
                        classified.setdefault((contact_type, contact_class), set()).add((int(unp[0]), copy))
        return classified

    @staticmethod
    def add_statistics(columns: dict, prefix: str, count_name: str, values: dict[int, list[float]], resnums: list[int]):
        columns[f"{prefix}_min"] = [min(values[r]) if r in values else np.nan for r in resnums]
        columns[f"{prefix}_mean"] = [float(np.mean(values[r])) if r in values else np.nan for r in resnums]
        columns[f"{prefix}_max"] = [max(values[r]) if r in values else np.nan for r in resnums]
        columns[count_name] = [len(values.get(r, [])) for r in resnums]


def contact_classes(contact) -> list[str]:
    # the classes of CONTACT_CLASSES that a contact satisfies
    classes = [CONTACT_CLASSES[0]] if contact.is_distance_within(NO_CONTACT_DISTANCE) else []
    bonds = [bond_type for bond_type in BondLibrary.bond_types_implemented()
             if contact.satisfies_bond_type(bond_type, BondLibrary.default_method)]
    return classes + bonds + (["any_bond"] if bonds else [])


_worker_builder = None


def _init_worker_builder(distance_cutoff: float, contact_filter_key: str):
    global _worker_builder
    _worker_builder = _LysinePropertyBuilder(distance_cutoff, contact_filter_key)


def _worker_properties_of_protein(item: tuple[str, list[int]]) -> tuple[str, pd.DataFrame]:
    uniprot_id, resnums = item
    return uniprot_id, _worker_builder.properties_of_protein(uniprot_id, resnums)


def report_progress(i: int, total: int, show_progress: bool):
    if show_progress and (i % 100 == 0 or i == total):
        print(f"Lysine properties: {i} of {total} uniprot ids")
