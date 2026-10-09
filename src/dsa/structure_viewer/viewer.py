"""
Interactive 3D viewer for PDB biological assemblies and AlphaFold models (3Dmol.js in an HTML page).

    from dsa.structure_viewer.viewer import StructureViewer
    viewer = StructureViewer.from_pdb("4ah6")          # or StructureViewer.from_alphafold("Q6PI48")
    viewer.show()                                     # in a notebook
    viewer.save("4ah6.html")                          # standalone page, open in a browser

Hovering a residue shows its chain, UniProt residue, residue name and SASA (assembly, chain alone, AlphaFold),
and highlights the residues it is in contact with, coloured by bond type (salt bridge, hydrogen bond,
cation-pi, covalent isopeptide) or by distance (contact <= 4.5 A, <= 6.5 A). Clicking a residue pins it.

For a PDB structure a second panel shows the AlphaFold model of the selected chain's UniProt id, superposed on the
chain and rotated / zoomed together with it (active when exactly one chain is checked; see alphafold_panel.py).

3Dmol.js (vendor/3Dmol-min.js, BSD-3-Clause) is inlined into the page, so it works without internet access.
Its shaders light every surface only by max(dot(normal, light), 0) from one light, so surfaces turned away from the
light render black and dark patches move while rotating; the page gets an ambient term added (see three_dmol_source).
"""
import html
import json
import math
from pathlib import Path
import warnings
from dsa.structure_viewer.residue_data import StructureData, load_pdb, load_alphafold
from dsa.structure_viewer.alphafold_panel import load_alphafold_models, fit_chains

TEMPLATE = Path(__file__).with_name("viewer_template.html")
THREEDMOL = Path(__file__).parent / "vendor" / "3Dmol-min.js"
PALETTE = ["#4e79a7", "#f28e2b", "#59a14f", "#b07aa1", "#76b7b2", "#edc948", "#ff9da7", "#9c755f",
           "#86bcb6", "#d4a6c8", "#a0cbe8", "#ffbe7d", "#8cd17d", "#b6992d", "#f1ce63", "#79706e"]
LARGE_MODEL_ATOMS = 150_000
LIGHT_MIX = 0.45   # chains: palette colour mixed 45% towards white
DARK_MIX = 0.30    # lysines: the same colour mixed 30% towards black
AMBIENT_LIGHT = 0.4   # shading = AMBIENT_LIGHT + (1 - AMBIENT_LIGHT) * max(dot(normal, light), 0); 3Dmol's own is 0


class StructureViewer:
    def __init__(self, data: StructureData, color_by: str = "chain", alphafold_panel: bool = True):
        # color_by: "chain" (every bio chain its own colour) or "protein" (every UniProt id its own colour)
        # alphafold_panel: for a PDB structure, add the panel with the AlphaFold model of the selected chain
        if color_by not in ("chain", "protein"):
            raise ValueError("color_by must be 'chain' or 'protein'")
        self.data = data
        self.color_by = color_by
        self.alphafold_models, self.fits = {}, {}
        if alphafold_panel and data.database == "pdb":
            self.alphafold_models = load_alphafold_models(data)
            self.fits = fit_chains(data, self.alphafold_models)
        n_atoms = len(data.atom_residue)
        if n_atoms > LARGE_MODEL_ATOMS:
            warnings.warn(f"{n_atoms:,} atoms; the page may be slow. Pass chains=[...] to show part of the assembly.")

    @classmethod
    def from_pdb(cls, structure_id: str, chains: list[str]|None = None, **kwargs) -> "StructureViewer":
        """Biological assembly of a PDB structure; chains limits it to some bio chains (e.g. ["A1", "B1"])."""
        return cls(load_pdb(structure_id.lower(), chains), **kwargs)

    @classmethod
    def from_alphafold(cls, uniprot_id: str, **kwargs) -> "StructureViewer":
        return cls(load_alphafold(uniprot_id.upper()), **kwargs)

    @property
    def residues(self):
        """The residue table shown on hover (one row per residue), useful to look up values directly."""
        return self.data.residues

    @property
    def bonds(self):
        """The residue pairs in contact, with the residue table's index in res1 / res2."""
        import pandas as pd
        bonds = pd.DataFrame(self.data.bonds)
        if bonds.empty:
            return bonds
        label = lambda i: f"{self.residues.at[i, 'resname3']} {self.residues.at[i, 'resnum']}{self.residues.at[i, 'icode']} {self.residues.at[i, 'bio_chain']}"
        return bonds.assign(residue1=bonds["res1"].map(label), residue2=bonds["res2"].map(label))

    def chain_colors(self) -> dict[str, dict[str, str]]:
        """Per bio chain a light shade (the chain) and a dark shade of the same colour (its lysines)."""
        residues = self.data.residues
        if self.color_by == "protein":
            group = residues.groupby("bio_chain", sort=False)["uniprot_id"].agg(
                lambda s: s.dropna().iloc[0] if s.notna().any() else None)
            keys = list(dict.fromkeys(group.tolist()))
            colors = {key: PALETTE[i % len(PALETTE)] for i, key in enumerate(keys)}
            base = {chain: colors[key] if key is not None else "#9e9e9e" for chain, key in group.items()}
        else:
            chains = list(dict.fromkeys(residues["bio_chain"]))
            base = {chain: PALETTE[i % len(PALETTE)] for i, chain in enumerate(chains)}
        return {chain: {"light": mix(color, "#ffffff", LIGHT_MIX), "dark": mix(color, "#000000", DARK_MIX)}
                for chain, color in base.items()}

    def page_data(self) -> dict:
        page = {**structure_payload(self.data), "colors": self.chain_colors()}
        page["alphafold"] = {uniprot_id: structure_payload(model) for uniprot_id, model in self.alphafold_models.items()}
        page["fits"] = self.fits
        return page

    def to_html(self) -> str:
        data = json.dumps(self.page_data(), separators=(",", ":")).replace("</", "<\\/")
        page = TEMPLATE.read_text()
        page = page.replace("__TITLE__", html.escape(self.data.title)).replace("__DATA__", data)
        return page.replace("__3DMOL__", three_dmol_source())   # last, so the library text is not searched

    def save(self, path: str|Path) -> Path:
        path = Path(path)
        path.write_text(self.to_html())
        return path

    def show(self, height: int|None = None):
        """Displays the viewer in a notebook (inside an iframe, so several viewers can be shown at once)."""
        from IPython.display import HTML
        if height is None:
            height = 720 if self.alphafold_models else 650
        page = html.escape(self.to_html(), quote=True)
        with warnings.catch_warnings():
            # IPython suggests IFrame, which needs a URL; the page is passed inline with srcdoc instead
            warnings.filterwarnings("ignore", message="Consider using IPython.display.IFrame")
            return HTML(f'<iframe srcdoc="{page}" style="width:100%;height:{height}px;border:1px solid #d0d7de"></iframe>')

    def _repr_html_(self):
        return self.show().data


def structure_payload(data: StructureData) -> dict:
    """What the page needs of one structure: the model, the residue table (as columns) and the residue pairs."""
    residues = data.residues
    return {"title": data.title,
            "database": data.database,
            "cif": data.model.make_mmcif_document().as_string(),
            "residues": {column: [json_value(v) for v in residues[column].tolist()] for column in residues.columns},
            "atom_residue": data.atom_residue,
            "bonds": [{k: json_value(v) for k, v in b.items()} for b in data.bonds]}


def three_dmol_source(ambient: float = AMBIENT_LIGHT) -> str:
    """
    The vendored 3Dmol.js with an ambient term in its lighting: in every shader, max(dot(normal, light), 0) becomes
    ambient + (1 - ambient) * max(...). Surfaces facing the light keep their brightness; surfaces turned away keep
    `ambient` of their colour instead of going black. The vendored file itself is not changed.
    """
    source = THREEDMOL.read_text()
    patches = {"max( dotProduct, 0.0 )": 5, "max( -dotProduct, 0.0 )": 1}   # occurrences in 3Dmol.js 2.4.2
    for term, expected in patches.items():
        found = source.count(term)
        if found != expected:
            raise RuntimeError(f"3Dmol.js shader changed: {found} occurrences of '{term}', expected {expected}")
        source = source.replace(term, f"( {ambient:.3f} + {1 - ambient:.3f} * {term} )")
    return source


def mix(color: str, other: str, fraction: float) -> str:
    """color moved by fraction towards other (both #rrggbb)."""
    a = [int(color[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(other[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * fraction):02x}" for x, y in zip(a, b))


def json_value(value):
    # numpy / pandas scalars and missing values to plain JSON values
    if value is None:
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    try:
        import pandas as pd
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value
