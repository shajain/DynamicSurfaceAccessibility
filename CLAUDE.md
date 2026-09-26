# DynamicSurfaceAccessibility (package `dsa`)

Compares lysine labelling from mass spectrometry (MS_DSA, light/heavy dimethyl labels) with the structural
accessibility and contacts of the lysines in PDB and AlphaFold structures.

**Read first:** [docs/project_background.md](docs/project_background.md) (architecture, data flow, conventions).
**History:** [docs/worklog/](docs/worklog/) — dated notes of what was done and why (latest: 2026-09-26).
**Presentation notes:** `src/dsa/presentation_notes.ipynb` (terminology, design choices, sanity checks).

## Environment
- Python: conda env `dsa` in `/data/shajain/container_setup/miniforge3` → `/data/shajain/container_setup/miniforge3/envs/dsa/bin/python`.
  The system `python3` lacks the dependencies. `~/miniforge3` must be a symlink to that install (its scripts use
  `/home/shajain/miniforge3` as prefix); `/data/shajain/container_setup/setup_env.sh` creates it.
- The env has the FreeSASA command-line tool (`freesasa-c` 2.1.2) and the Python bindings.
- Run modules from `src/` as `python -m dsa.<module>`; the package is installed editable.
- **`src/dsa/data/` is gitignored** (≈ 2.5 GB + structures). AlphaFold structures: `data/structures/alphafold` is a
  symlink to `/data/dbs/alphafold_db`. On a new machine the data must be mounted/copied or regenerated **in order**:
  UniProt ids/sequences → UniProt→PDB mapping + PDB files, AlphaFold files → PDB→UniProt alignments → structure
  contacts → protein contacts → SASA (monomer, multimer) and B-factors → MS_DSA v2 → lysine/run tables → notebooks.
  Commands, inputs and outputs of every step: section 5 of `docs/project_background.md`.

## Conventions (follow them)
- **Store + compute framework:** results are saved through a `Store` subclass (`misc/store.py`, usually
  `CompressedPickleStore`), one file per key; batch jobs use `ComputeManager` (`misc/compute.py`) with parallel
  workers launched in tmux by the `compute_*.sh` scripts at the repo root.
- **Never combine AlphaFold and PDB values into one number**; plots may show them side by side.
- MS_DSA is always computed **per run**.
- Keep new work in new files unless asked to edit existing ones; the user often asks for plans before code.
- Add a dated file to `docs/worklog/` for substantial work.

## Current analysis entry points
- MS_DSA with flags: `ms_dsa/ms_dsa_v2.py` → `data/ms_data/report_ms_dsa_v2.csv`
- Per-lysine and per-peptide tables: `ms_dsa/evaluation1/` (`python -m dsa.ms_dsa.evaluation1.build_annotations_v2`
  [`--rebuild`]) → `data/ms_dsa/evaluation1/{lysine_properties,peptide_lysines}.parquet`, `runs/run_{0..4}.parquet`
- Plots: `ms_dsa/evaluation1/plots.py`; notebooks `src/dsa/eval6.ipynb`, `eval7_monomer.ipynb`,
  `eval7_multimer.ipynb`, `eval8_venn.ipynb` (older ones in `src/dsa/obsolete/`).
