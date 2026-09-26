# Project background: DynamicSurfaceAccessibility (`dsa`)

Audience: anyone (human or Claude Code) picking up this repository on a new machine. Written 2026-09-26 from
the code as it stood then. For the session history see `docs/worklog/`; for terminology, design rationale and sanity
checks aimed at presenting the results see `src/dsa/presentation_notes.ipynb`.

## 1. Scientific goal

Lysines of proteins in solution are labelled with light / heavy dimethyl labels; the MS readout
**MS_DSA = light / (light + heavy)** per peptide (and run) reflects how accessible the peptide's lysine(s) were.
The project asks whether MS_DSA agrees with structure-based measures for the same lysines:
- solvent accessibility of the lysine NZ atom (`asa_nz`) in AlphaFold monomers and PDB biological assemblies,
- whether the lysine is at an inter-chain interface, in an intra-chain contact, or in no contact,
- bond types (salt bridge, hydrogen bond, cation-π, covalent isopeptide) and B-factors.

## 2. Data flow

```
data/ms_data/report.parquet  (DIA-NN report: runs, channels '0' light / '8' heavy, Ms1.Area, charges)
   │
   ├─ uniprot/            ALL_UNIPROT_IDS from the MS data; sequences cached in data/sequences.json
   ├─ binding_interfaces/mappings/   UniProt → PDB ids (RCSB search), data/uniprot_structure_mappings
   ├─ structures          data/structures/pdb/*.cif, data/structures/alphafold → /data/dbs/alphafold_db (pdb.gz)
   │                      PDB→UniProt residue alignment: PDBe (data/structures/alignment/pdbe) or SIFTS
   │
   ├─ contacts            StructureContacts per structure → ProteinContacts per protein (both PDB and AlphaFold)
   ├─ sasa/               SASAMonomer (AlphaFold) and SASAMultimer (PDB assemblies)
   ├─ bfactor/            LysineBFactor (PDB, per author chain)
   │
   ├─ ms_dsa/ms_dsa_v2.py MS_DSA + flags per (peptide, run)       → data/ms_data/report_ms_dsa_v2.csv
   └─ ms_dsa/evaluation1/ lysine table + per-run peptide tables     → data/ms_dsa/evaluation1/
          └─ notebooks eval6, eval7_monomer, eval7_multimer, eval8_venn (+ plots.py)
```

## 3. Code map (`src/dsa/`)

| path | role |
|---|---|
| `config.py` | `PROJECT_ROOT`, `DATA_DIR = src/dsa/data`, `RESULTS_DIR` |
| `misc/store.py` | `Store` base class (one file per key, `dump`/`load`/`get_saved_keys`); `PickleStore`, `CompressedPickleStore` (zstd), `JSONStore`, `SingleFileTabularStore` |
| `misc/compute.py` | `ComputeManager(keys, compute_function, store)`: processes random unsaved keys, saves non-None results, refreshes the key list every `update_interval`; safe to run many workers in parallel |
| `misc/downloader.py` | retrying HTTP `Downloader` base class |
| `uniprot/` | `ALL_UNIPROT_IDS` (from all MS parquet files), `UniprotToSequence` (UniProt REST, cached JSON) |
| `binding_interfaces/factory.py` | `ClassFactory.register/create(interface, key)` — picks `pdb` / `alphafold` implementations |
| `binding_interfaces/structure/` | `StructureFileStore` (PDB cif, AlphaFold pdb.gz v4/v6), `Structure` / `PDBStructure` / `AlphaFoldStructure` (gemmi; entities; first biological assembly via `make_assembly(AddNumber)` → bio chains `A1`, `A2`; `biochain_to_chain`; `map_residue_to_uniprot`), `StructureModel` |
| `binding_interfaces/structure_uniprot_alignment/` | `PDBSequenceAligner` (PDBe residue mappings, SIFTS fallback), `AlphaFoldSequenceAligner`; chain accepted if ≤ 10% mismatches (> 10 residues, < 20 mismatches); `alignment_store.py` (PDBe download) |
| `binding_interfaces/mappings/` | `UniprotStructureMappingsStore` / mapper: UniProt → PDB ids |
| `binding_interfaces/contacts/` | `ContactExtractor` (gemmi NeighborSearch, 6.5 Å, adjacent residues ignored, `image_idx == 0`), `LysineOnlyFilter` (lysine NZ as contact atom, both polymers, not adjacent), `Contact` (chains, residues, `dist`, `is_inter_chain`, `satisfies_bond_type`), `StructureContacts(+Store)`, `UniprotAlignedStructureContacts` (contacts of one protein with UniProt residues), `ProteinContacts(+Store)` (all structures of a protein; `grouped_residues`, `does_residue_satisfy`) |
| `binding_interfaces/contacts/bond_characterization/` | `BondLibrary`: salt bridge (NZ–carboxylate ≤ 4.0 Å), hydrogen bond (≤ 3.5 Å), cation-π (≤ 6.0 Å, ≤ 40°), covalent isopeptide (≤ 1.6 Å); method `assume_lysine` |
| `sasa/sasa_monomer.py` | `SASAMonomer(+Store)`: FreeSASA Python bindings on the AlphaFold model; per-residue `asa`, `rel_asa`, `asa_nz` → `data/sasa/monomers1` |
| `sasa/sasa_multimer.py` | `SASAMultimer(+Store)`: FreeSASA CLI on the protein-only biological assembly; unique residue labels preserve bio chains; one row per residue per chain copy with UniProt mapping → `data/sasa/multimers` |
| `sasa/generate_*`, `import_sasa_monomer_csvs.py` | batch jobs (ComputeManager); CSV → store import |
| `sasa/compute_sasa.py` | original CSV-based monomer SASA (superseded; `data/sasa/monomers`) |
| `bfactor/lysine_bfactor.py` | `LysineBFactor(+Store)`: NZ and residue-average B-factors per author chain, z-scored per chain → `data/bfactor/pdb` |
| `ms_dsa/ms_dsa.py` | original MS_DSA (all charges; `data/ms_data/ms_dsa.csv`) — still imported by `uniprot/config.py` |
| `ms_dsa/ms_dsa_v2.py` | MS_DSA and `MS_DSA_shared_charges` with flags (`channel_absent`, `no_shared_charges`, `zero_area`, …) |
| `ms_dsa/evaluation1/lysine_properties.py` | `LysinePropertyTable`: per-lysine contacts flags, monomer/multimer SASA, B-factors, contact-class columns (per chain copy) |
| `ms_dsa/evaluation1/ms_dsa_annotator.py` | `MSDSAAnnotator`: peptides → lysines (first exact match), peptide aggregation, per-run tables |
| `ms_dsa/evaluation1/build_annotations_v2.py` | builds/extends the lysine table and writes the run tables (`--rebuild` recomputes all lysines) |
| `ms_dsa/evaluation1/plots.py` | plotting helpers (contact-class groups, run-averaged histograms, pooled scatter, Venn) |
| `ms_dsa/Evaluation/binding_interface_eval.py` | original evaluation (per-run CSVs `ms_dsa_with_binding_interface_and_sasa_scores_{run}.csv`), superseded |
| `obsolete/`, `code_backup/`, `*_delete_later/` | retired code and notebooks |

## 4. Stores and data locations (all under `src/dsa/data/`, gitignored)

| data | location |
|---|---|
| MS report | `ms_data/report.parquet`; MS_DSA: `ms_data/ms_dsa.csv` (v1), `ms_data/report_ms_dsa_v2.csv` (v2) |
| PDB structures | `structures/pdb/*.cif` (~126k files) |
| AlphaFold structures | `structures/alphafold` → `/data/dbs/alphafold_db` |
| PDB→UniProt alignments | `structures/alignment/` (PDBe per structure; SIFTS `pdb_chain_uniprot.tsv`) |
| UniProt→structure mappings | `uniprot_structure_mappings/` |
| contacts | `binding_interfaces/structure_contacts/{pdb,alphafold}_65-LysineOnlyFilter`, `binding_interfaces/protein_contacts/{pdb,alphafold}_65-LysineOnlyFilter` |
| SASA | `sasa/monomers1` (AlphaFold, store), `sasa/multimers` (PDB, 37,683 structures), `sasa/monomers` (old CSVs) |
| B-factors | `bfactor/pdb` |
| evaluation | `ms_dsa/evaluation1/lysine_properties.parquet`, `peptide_lysines.parquet`, `runs/run_{0..4}.parquet` |

## 5. Pipeline order on a fresh machine (data is not in git)

`src/dsa/data/` is gitignored, so a new machine either mounts / copies the existing data (`/data/shajain/...`,
`/data/dbs/alphafold_db`) or recomputes it. The steps depend on each other and must run **in this order**. Commands
run from `src/` with the `dsa` env (`python -m dsa.<module>`); the `*.sh` scripts at the repo root start parallel
tmux workers. Several modules take their settings (e.g. `db = "pdb"` / `"alphafold"`) from their `__main__` block.

| # | step | needs | produces | how |
|---|---|---|---|---|
| 0 | **MS data** | the DIA-NN report (provided by the lab) | `data/ms_data/report.parquet` | copy in |
| 1 | **UniProt ids** | 0 | `ALL_UNIPROT_IDS` (built on import of `dsa.uniprot.config` from every MS parquet; also writes the v1 cache `ms_data/ms_dsa.csv`) | automatic |
| 2 | **UniProt sequences** | 1, internet (`rest.uniprot.org`) | `data/sequences.json` | automatic: `UniprotToSequence` fetches missing sequences when first used; or `python -m dsa.uniprot.uniprot_to_sequence` |
| 3 | **UniProt → PDB mapping + PDB files** | 1, internet (RCSB search, PDB download) | `data/uniprot_structure_mappings/`, `data/structures/pdb/*.cif` | `python -m dsa.binding_interfaces.mappings.uniprot_to_pdb` |
| 4 | **AlphaFold files** (mapping is the UniProt id itself: `AF-{id}-F1-model_v4/v6.pdb.gz`) | 1, internet (`alphafold.ebi.ac.uk`) or the local AlphaFold DB | `data/structures/alphafold/` (here a symlink to `/data/dbs/alphafold_db`) | `python -m dsa.binding_interfaces.mappings.uniprot_to_alphaFold` downloads the missing ones (note: it writes into the symlinked shared DB) |
| 5 | **PDB → UniProt residue alignments** | 3, internet (EBI SIFTS) | `data/structures/alignment/pdbe/` (per-structure SIFTS XML, the primary source) and `data/structures/alignment/pdb_chain_uniprot.tsv` (SIFTS flat file, fallback) | `python -m dsa.binding_interfaces.structure_uniprot_alignment.alignment_store`; the flat file has no downloader in the code — get `pdb_chain_uniprot.tsv` from the EBI SIFTS flat files |
| 6 | **Structure contacts** (per structure, PDB and AlphaFold) | 3, 4 | `data/binding_interfaces/structure_contacts/{pdb,alphafold}_65-LysineOnlyFilter` | `generate_structure_contacts.sh` → `dsa.binding_interfaces.contacts.generate_contacts` (set `db`) |
| 7 | **Protein contacts** (per protein: its structures' contacts mapped to UniProt residues) | 2, 3, 5, 6 | `data/binding_interfaces/protein_contacts/{pdb,alphafold}_65-LysineOnlyFilter` | `generate_protein_contacts.sh` → `dsa.binding_interfaces.contacts.generate_protein_contacts` (set `db`) |
| 8a | **Monomer SASA** | 4 | `data/sasa/monomers1` | `python -m dsa.sasa.generate_sasa_monomers` (or import old CSVs: `dsa.sasa.import_sasa_monomer_csvs`) |
| 8b | **Multimer SASA** | 3, 5 (UniProt mapping), FreeSASA CLI | `data/sasa/multimers` | `bash compute_multimer_sasa.sh [workers] [update_interval]` (~30 min, 24 workers) |
| 8c | **B-factors** | 3, 5 | `data/bfactor/pdb` | `bash compute_bfactor.sh` |
| 9 | **MS_DSA v2** | 0 | `data/ms_data/report_ms_dsa_v2.csv` | `python -m dsa.ms_dsa.ms_dsa_v2` |
| 10 | **Lysine table + run tables** | 2, 7, 8a–c, 9 | `data/ms_dsa/evaluation1/lysine_properties.parquet`, `peptide_lysines.parquet`, `runs/run_{0..4}.parquet` | `python -m dsa.ms_dsa.evaluation1.build_annotations_v2 --rebuild` (~13 min, 24 workers; without `--rebuild` only missing lysines are added) |
| 11 | **Notebooks** | 10 | figures | `jupyter nbconvert --to notebook --execute --inplace src/dsa/eval6.ipynb` (also `eval7_*`, `eval8_venn`, `presentation_notes`) |

Notes
- Steps 6–7 and 8a–c are the expensive ones (hours); step 7 is slow for proteins with many structures because
  unpickling contacts rebuilds each structure.
- Steps 8a–c are independent of each other and of 6–7, so they can run in parallel with them.
- `download_alignments.sh` refers to a module name (`pdbe_downloader`) that no longer exists; use step 5 instead.
- All batch steps are resumable: stores skip keys that are already saved.

## 6. Key design decisions

- **Biological assembly** (first annotated) for contacts and multimer SASA; protein residues only for SASA
  (modified amino acids kept; ligands, nucleic acids, water, hydrogens, altlocs removed).
- **FreeSASA chain names are truncated to one character**, so SASAMultimer never relies on them: residues get unique
  labels (≤ 5 characters, up to 671,947 labels) and are mapped back to bio chain / author chain / UniProt position.
- **Contact classes** (`dist4.5`, four bond types, `any_bond`) computed per chain copy; the inter-chain `asa_nz`
  of a lysine uses only the copies that have the inter-chain contact.
- **No contact** = no contact within 4.5 Å in any copy of that source (PDB or AlphaFold separately).
- **PDB and AlphaFold values are never combined into one number.**
- **MS_DSA per run**; shared-charge MS_DSA to avoid bias from charge states seen in only one channel; missing /
  zero channels are flagged, not silently dropped.
- Peptide aggregation: max-of-max, min-of-min, average of per-lysine means; flags true if any lysine (no contact:
  all lysines).

## 7. Known issues and caveats

- Two large proteins (Q71U36, P68371) have no saved PDB ProteinContacts; 9i1b has no multimer SASA (faulty
  assembly with duplicated chain copies).
- `Structure.uniprot_to_chains` and the SIFTS path in `PDBSequenceAligner` had bugs, fixed in commit a8a0ccf.
  SIFTS rows with an empty `PDB_BEG` are skipped (e.g. 6obj unmapped).
- `ComputeManager` retries a failed key only until its next refresh; with slow failing keys, stop workers once the
  saved count stops growing.
- `ms_dsa.py` (v1) caches to a fixed `ms_data/ms_dsa.csv` regardless of the input parquet; v2 names the cache after
  the input.
- Aligner warnings (`Invalid residue id: ('null', …)`) are noisy but harmless; worker logs can be large.
- Multimer `asa_nz` > 60 Å² (16 lysines) are model artifacts (stretched CE–NZ bonds, missing CE).
