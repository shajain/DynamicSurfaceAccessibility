"""
Plots of the annotated MS_DSA runs (MSDSAAnnotator.load_run), ported from eval2.ipynb to the v2 tables:

- plot_dsa_vs_sasa_with_binding_annotation : MS_DSA vs a SASA column, coloured by contact class
- binned_histogram                         : percent of peptides per bin, per contact class, with
                                             user supplied bin edges

Differences from eval2.ipynb, which read ms_dsa_with_binding_interface_and_sasa_scores_{run}.csv:
- rows are kept when their protein has both PDB and AlphaFold contacts and the peptide was found in the
  uniprot sequence (replaces the old "skipped" column)
- dsa_column selects MS_DSA (all charges) or MS_DSA_shared_charges
- exclude drops rows flagged in any of the given boolean columns, e.g. ["channel_absent", "zero_area"]
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CONTACT_COLUMNS = {
    ("inter_chain", "pdb"): "inter_chain_pdb",
    ("intra_chain", "pdb"): "intra_chain_pdb",
    ("intra_chain", "alphafold"): "intra_chain_alphafold",
}
GROUP_COLORS = {"Inter": "tab:red", "Intra AF": "tab:green", "Not Contact": "grey"}


def select_rows(ms_dsa: pd.DataFrame, exclude: list[str] = ()) -> pd.DataFrame:
    # rows with contact annotations (both sources) for a peptide located in its uniprot sequence
    keep = ms_dsa["peptide_mapped"] & ms_dsa["inter_chain_pdb"].notna() & ms_dsa["intra_chain_alphafold"].notna()
    for flag in exclude:
        keep &= ~ms_dsa[flag].astype(bool)
    return ms_dsa[keep]


def _mask(ms_dsa: pd.DataFrame, contact_type: str, source: str = "pdb", bond_type: str|None = "salt_bridge",
          strict: bool = False) -> pd.Series:
    mask = ms_dsa[CONTACT_COLUMNS[(contact_type, source)]].fillna(False).astype(bool)
    if strict and contact_type == "inter_chain":
        mask &= ~ms_dsa["intra_chain_pdb"].fillna(False).astype(bool)
    if bond_type is not None:
        mask &= ms_dsa[f"{bond_type}_{source}"].fillna(False).astype(bool)
    return mask


def _not_contact(ms_dsa: pd.DataFrame, bond_type: str|None = None, strict: bool = False) -> pd.Series:
    return (~_mask(ms_dsa, "inter_chain", "pdb", bond_type, strict)
            & ~_mask(ms_dsa, "intra_chain", "pdb", bond_type)
            & ~_mask(ms_dsa, "intra_chain", "alphafold", bond_type))


def _column_for(column: str|dict[str, str], group: str) -> str:
    # column is one column for every group, or {group: column} with the groups "Inter", "Intra AF", "Not Contact"
    return column[group] if isinstance(column, dict) else column


def _column_label(column: str|dict[str, str]) -> str:
    if not isinstance(column, dict):
        return column
    return "; ".join(f"{group}: {c}" for group, c in column.items())


def _exclude_note(exclude: list[str]) -> str:
    return f"\nexcluding {', '.join(exclude)}" if exclude else "\nall peptides"


def _frac_lt(df: pd.DataFrame, column: str, thresh: float = 0.9) -> float:
    values = df[column].dropna()
    return (values < thresh).mean() if len(values) else np.nan


def plot_dsa_vs_sasa_with_binding_annotation(
    ms_dsa: pd.DataFrame,
    contact_type: str = "inter_chain",
    bond_type: str|None = None,
    source: str = "pdb",
    sasa_column: str|dict[str, str] = "asa_nz_max",
    dsa_column: str = "MS_DSA",
    strict: bool = False,
    exclude: list[str] = (),
    run_index: int|None = None,
    ax: plt.Axes|None = None,
):
    ms_dsa = select_rows(ms_dsa, exclude)
    inter = ms_dsa[_mask(ms_dsa, "inter_chain", "pdb", bond_type, strict)]
    intra = ms_dsa[_mask(ms_dsa, "intra_chain", "alphafold", bond_type)]
    intra_pdb = ms_dsa[_mask(ms_dsa, "intra_chain", "pdb", bond_type)]
    not_contact = ms_dsa[_not_contact(ms_dsa, bond_type, strict)]
    highlight = ms_dsa[_mask(ms_dsa, contact_type, source, bond_type, strict)]

    for name, df in [("Inter PDB", inter), ("Intra PDB", intra_pdb), ("Intra AF", intra), ("Not contact", not_contact)]:
        print(f"{name}: {len(df)} rows  ({_frac_lt(df, dsa_column):.3f} with {dsa_column} < 0.9)")

    if ax is None:
        _, ax = plt.subplots()
    ax.scatter(not_contact[dsa_column], not_contact[_column_for(sasa_column, "Not Contact")], color="grey",
               label="Not Contact", s=20, alpha=0.5)
    ax.scatter(intra[dsa_column], intra[_column_for(sasa_column, "Intra AF")], color="green", label="Intra AF", s=10, alpha=0.5)
    ax.scatter(highlight[dsa_column], highlight[_column_for(sasa_column, "Inter")], color="red",
               label=f"{contact_type}/{source}", s=5, alpha=0.5)
    ax.set_xlabel(dsa_column)
    ax.set_ylabel("asa_nz" if isinstance(sasa_column, dict) else sasa_column)
    ax.set_title(f"{dsa_column} vs {_column_label(sasa_column)} on Run {run_index}{_exclude_note(exclude)}")
    ax.legend()
    return ax


def contact_groups(ms_dsa: pd.DataFrame, source: str = "pdb", bond_type: str|None = None,
                   strict: bool = False) -> dict[str, pd.DataFrame]:
    return {
        "Inter": ms_dsa[_mask(ms_dsa, "inter_chain", source, bond_type, strict)],
        "Intra AF": ms_dsa[_mask(ms_dsa, "intra_chain", "alphafold", bond_type)],
        "Not Contact": ms_dsa[_not_contact(ms_dsa, bond_type, strict)],
    }


def correlation_table(ms_dsa_by_run: dict[int, pd.DataFrame], dsa_column: str, value_column: str|dict[str, str],
                      source: str = "pdb", bond_type: str|None = None, strict: bool = False,
                      exclude: list[str] = (), method: str = "spearman") -> pd.DataFrame:
    """Correlation (and number of peptides) between dsa_column and value_column per run and contact class."""
    rows = []
    for run_index, ms_dsa in ms_dsa_by_run.items():
        for group, df in contact_groups(select_rows(ms_dsa, exclude), source, bond_type, strict).items():
            pair = df[[dsa_column, _column_for(value_column, group)]].dropna()
            rows.append({"run": run_index, "group": group, "n": len(pair),
                         method: pair.iloc[:, 0].corr(pair.iloc[:, 1], method=method) if len(pair) > 2 else np.nan})
    table = pd.DataFrame(rows).pivot(index="run", columns="group", values=[method, "n"])
    return table.reindex(columns=["Inter", "Intra AF", "Not Contact"], level=1).round(3)


def bin_labels(bins, decimals: int = 3) -> list[str]:
    edge = lambda x: f"{round(float(x), decimals):g}"   # e.g. 0.3666... -> 0.367
    labels = []
    for low, high in zip(bins[:-1], bins[1:]):
        if np.isinf(high):
            labels.append(f"≥{edge(low)}")
        elif np.isinf(low):
            labels.append(f"<{edge(high)}")
        else:
            labels.append(f"{edge(low)}–{edge(high)}")
    return labels


def binned_histogram(
    ms_dsa: pd.DataFrame,
    column: str|dict[str, str],
    bins,
    source: str = "pdb",
    bond_type: str|None = None,
    strict: bool = False,
    exclude: list[str] = (),
    require: str|None = None,
    run_index: int|None = None,
    ax: plt.Axes|None = None,
):
    """
    Percent of the peptides of each contact class falling in each bin of `column`
    (one column, or {group: column} to use a different column for "Inter", "Intra AF" and "Not Contact").
    bins are the bin edges, e.g. [0, 0.1, 0.9, 1] for MS_DSA or [0, 20, 40, np.inf] for asa_nz_max;
    as in np.histogram the last bin includes its right edge, and values outside the edges are not counted.
    require: only count peptides with a value in this column, e.g. the MS_DSA column of the other plots.
    """
    bins = np.asarray(bins, dtype=float)
    ms_dsa = select_rows(ms_dsa, exclude)
    if require is not None:
        ms_dsa = ms_dsa[ms_dsa[require].notna()]
    groups = contact_groups(ms_dsa, source, bond_type, strict)
    labels = bin_labels(bins)
    x = np.arange(len(labels))
    width = 0.8 / len(groups)
    if ax is None:
        _, ax = plt.subplots()
    for i, (name, df) in enumerate(groups.items()):
        values = df[_column_for(column, name)].dropna().values
        counts, _ = np.histogram(values, bins=bins)
        n_binned = counts.sum()
        pct = 100.0 * counts / n_binned if n_binned else np.zeros(len(counts))
        outside = len(values) - n_binned
        label = f"{name} ({n_binned}" + (f", {outside} outside bins)" if outside else ")")
        ax.bar(x + (i - (len(groups) - 1) / 2) * width, pct, width, label=label, color=GROUP_COLORS[name], alpha=0.7)
    ax.set_xticks(x, labels, rotation=30 if len(labels) > 3 else 0, ha="right" if len(labels) > 3 else "center")
    ax.set_xlabel("asa_nz" if isinstance(column, dict) else column)
    ax.set_ylabel("Percent of peptides")
    ax.set_title(f"{_column_label(column)} binned peptides on Run {run_index}{_exclude_note(exclude)}")
    ax.legend()
    return ax


def plot_runs(ms_dsa_by_run: dict[int, pd.DataFrame], plot, ncols: int = 3, panel_size: tuple[float, float] = (3.6, 2.9),
              suptitle: str|None = None, **kwargs) -> plt.Figure:
    """
    One figure with a panel per run, runs from left to right over as many rows as needed.
    plot is plot_dsa_vs_sasa_with_binding_annotation or binned_histogram; kwargs are passed to it.
    """
    nrows = int(np.ceil(len(ms_dsa_by_run) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_size[0] * ncols, panel_size[1] * nrows),
                             squeeze=False, constrained_layout=True)
    for i, (ax, (run_index, ms_dsa)) in enumerate(zip(axes.flat, ms_dsa_by_run.items())):
        plot(ms_dsa, run_index=run_index, ax=ax, **kwargs)
        legend = ax.get_legend()
        handles = getattr(legend, "legend_handles", None) or getattr(legend, "legendHandles", None) if legend else None
        # keep the order of a legend the plot made itself
        ax.legend(handles=handles, fontsize="x-small") if handles else ax.legend(fontsize="x-small")
        if i == 0 and suptitle is None:
            # the panel title without the run ("<plot> on Run <i>\n<filter>") becomes the figure title
            suptitle = ax.get_title().replace(f" on Run {run_index}", "")
        ax.set_title(f"Run {run_index}")
    for ax in list(axes.flat)[len(ms_dsa_by_run):]:
        ax.set_visible(False)
    fig.suptitle(suptitle)
    return fig


# ---------------------------------------------------------------------------------------------------------------
# Contact-class groups (lysine_properties.CONTACT_CLASSES). A group is a set of conditions on the peptide flags
# and the column plotted for it, so every group can take its value from its own source:
#     groups = {"Inter PDB": ([("pdb_inter_salt_bridge", True)], "asa_nz_pdb_inter_salt_bridge_max_max"), ...}
# ---------------------------------------------------------------------------------------------------------------
CLASS_GROUP_COLORS = {"Inter PDB": "tab:red", "Intra PDB": "tab:blue", "Intra AF": "tab:green",
                      "No contact PDB": "dimgrey", "No contact AF": "darkgrey"}
AF_STATISTIC_COLUMNS = {"max": "asa_nz_max", "min": "asa_nz_min", "avg": "asa_nz_avg"}
PDB_STATISTIC_SUFFIX = {"max": "max_max", "min": "min_min", "avg": "mean_avg"}


def class_groups(contact_class: str, statistic: str = "max", purely_intra: bool = False,
                 value_columns: dict[str, str]|None = None) -> dict[str, tuple[list[tuple[str, bool]], str]]:
    """
    Inter PDB, Intra PDB, Intra AF, No contact PDB and No contact AF for one contact class, with the asa_nz
    statistic ("max": max of max, "min": min of min, "avg": average of the per-lysine mean). PDB groups use the
    multimer asa_nz of the chain copies in the class, AlphaFold groups the monomer asa_nz.
    purely_intra: the intra groups exclude peptides with an inter-chain bond (pdb_inter_any_bond).
    value_columns: {group: column} to plot another column (e.g. a B-factor) instead of asa_nz.
    """
    suffix = PDB_STATISTIC_SUFFIX[statistic]
    not_inter = [("pdb_inter_any_bond", False)] if purely_intra else []
    groups = {
        "Inter PDB": ([(f"pdb_inter_{contact_class}", True)], f"asa_nz_pdb_inter_{contact_class}_{suffix}"),
        "Intra PDB": ([(f"pdb_intra_{contact_class}", True)] + not_inter, f"asa_nz_pdb_intra_{contact_class}_{suffix}"),
        "Intra AF": ([(f"af_intra_{contact_class}", True)] + not_inter, AF_STATISTIC_COLUMNS[statistic]),
        "No contact PDB": ([("no_contact_pdb", True)], f"asa_nz_pdb_no_contact_{suffix}"),
        "No contact AF": ([("no_contact_af", True)], AF_STATISTIC_COLUMNS[statistic]),
    }
    if value_columns:
        groups = {name: (conditions, value_columns.get(name, column)) for name, (conditions, column) in groups.items()}
    return groups


def select_group_rows(ms_dsa: pd.DataFrame, exclude: list[str] = (), require: str|None = None) -> pd.DataFrame:
    keep = ms_dsa["peptide_mapped"].copy()
    for flag in exclude:
        keep &= ~ms_dsa[flag].astype(bool)
    if require is not None:
        keep &= ms_dsa[require].notna()
    return ms_dsa[keep]


def group_rows(ms_dsa: pd.DataFrame, conditions: list[tuple[str, bool]]) -> pd.DataFrame:
    # NA flags satisfy no condition
    keep = pd.Series(True, index=ms_dsa.index)
    for column, value in conditions:
        flag = ms_dsa[column].astype("boolean")
        keep &= (flag == value).fillna(False)
    return ms_dsa[keep]


def group_scatter(ms_dsa: pd.DataFrame, groups: dict, dsa_column: str = "MS_DSA_shared_charges",
                  exclude: list[str] = (), title: str = "", ylabel: str = "asa_nz", run_index: int|None = None,
                  ax: plt.Axes|None = None):
    ms_dsa = select_group_rows(ms_dsa, exclude, require=dsa_column)
    if ax is None:
        _, ax = plt.subplots()
    sizes = {"No contact PDB": 16, "No contact AF": 16}
    # drawn in reverse order so that the first groups (Inter PDB) are on top; legend kept in the group order
    handles = {}
    for name, (conditions, column) in reversed(list(groups.items())):
        df = group_rows(ms_dsa, conditions)
        df = df[df[column].notna()]
        handles[name] = ax.scatter(df[dsa_column], df[column], s=sizes.get(name, 8), alpha=0.5,
                                   color=CLASS_GROUP_COLORS.get(name), label=f"{name} ({len(df)})")
    ax.legend(handles=[handles[name] for name in groups])
    ax.set_xlabel(dsa_column)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title} on Run {run_index}{_exclude_note(exclude)}")
    return ax


def group_histogram(ms_dsa: pd.DataFrame, groups: dict, bins, dsa_column: str = "MS_DSA_shared_charges",
                    exclude: list[str] = (), title: str = "", run_index: int|None = None, ax: plt.Axes|None = None):
    # percent of each group's peptides per bin of its value column (peptides with a value in dsa_column)
    bins = np.asarray(bins, dtype=float)
    ms_dsa = select_group_rows(ms_dsa, exclude, require=dsa_column)
    labels = bin_labels(bins)
    x = np.arange(len(labels))
    width = 0.8 / len(groups)
    if ax is None:
        _, ax = plt.subplots()
    for i, (name, (conditions, column)) in enumerate(groups.items()):
        values = group_rows(ms_dsa, conditions)[column].dropna().values
        counts, _ = np.histogram(values, bins=bins)
        n_binned = counts.sum()
        pct = 100.0 * counts / n_binned if n_binned else np.zeros(len(counts))
        outside = len(values) - n_binned
        label = f"{name} ({n_binned}" + (f", {outside} outside bins)" if outside else ")")
        ax.bar(x + (i - (len(groups) - 1) / 2) * width, pct, width, label=label,
               color=CLASS_GROUP_COLORS.get(name), alpha=0.8)
    ax.set_xticks(x, labels, rotation=30 if len(labels) > 3 else 0, ha="right" if len(labels) > 3 else "center")
    ax.set_ylabel("Percent of peptides")
    ax.set_title(f"{title} on Run {run_index}{_exclude_note(exclude)}")
    ax.legend()
    return ax


def group_sizes(ms_dsa_by_run: dict[int, pd.DataFrame], groups: dict, dsa_column: str = "MS_DSA_shared_charges",
                exclude: list[str] = ()) -> pd.DataFrame:
    # number of peptides per run and group that have both a dsa_column value and a group value
    rows = {}
    for run_index, ms_dsa in ms_dsa_by_run.items():
        selected = select_group_rows(ms_dsa, exclude, require=dsa_column)
        rows[run_index] = {name: int(group_rows(selected, conditions)[column].notna().sum())
                           for name, (conditions, column) in groups.items()}
    return pd.DataFrame(rows).T.rename_axis("run")


def group_correlations(ms_dsa_by_run: dict[int, pd.DataFrame], groups: dict, dsa_column: str = "MS_DSA_shared_charges",
                       exclude: list[str] = (), method: str = "spearman") -> pd.DataFrame:
    rows = {}
    for run_index, ms_dsa in ms_dsa_by_run.items():
        selected = select_group_rows(ms_dsa, exclude, require=dsa_column)
        row = {}
        for name, (conditions, column) in groups.items():
            pair = group_rows(selected, conditions)[[dsa_column, column]].dropna()
            row[name] = pair.iloc[:, 0].corr(pair.iloc[:, 1], method=method) if len(pair) > 2 else np.nan
        rows[run_index] = row
    return pd.DataFrame(rows).T.rename_axis("run").round(3)


# ---------------------------------------------------------------------------------------------------------------
# All runs combined: pooled scatter plots and histograms whose per-group percentages are averaged over the runs
# ---------------------------------------------------------------------------------------------------------------
SASA_SOURCE_GROUPS = {
    # Inter PDB always uses the multimer asa_nz; the intra and no-contact groups use one source
    "monomer": ["Inter PDB", "Intra AF", "No contact AF"],
    "multimer": ["Inter PDB", "Intra PDB", "No contact PDB"],
}
BFACTOR_STATISTIC_COLUMNS = {"max": "nz_bfactor_z_max_max", "min": "nz_bfactor_z_min_min", "avg": "avg_bfactor_z_mean_avg"}


def source_groups(contact_class: str, sasa_source: str, statistic: str = "max", purely_intra: bool = False,
                  value_column: str|None = None) -> dict:
    """class_groups restricted to one SASA source; value_column plots one column (e.g. MS_DSA or a B-factor) for all groups."""
    groups = class_groups(contact_class, statistic, purely_intra)
    groups = {name: groups[name] for name in SASA_SOURCE_GROUPS[sasa_source]}
    if value_column is not None:
        groups = {name: (conditions, value_column) for name, (conditions, _) in groups.items()}
    return groups


def pooled_scatter(ms_dsa_by_run: dict[int, pd.DataFrame], groups: dict, dsa_column: str = "MS_DSA_shared_charges",
                   exclude: list[str] = (), title: str = "", ylabel: str = "asa_nz", ax: plt.Axes|None = None):
    # every run on one plot, runs not distinguished
    pooled = pd.concat([select_group_rows(ms_dsa, exclude, require=dsa_column) for ms_dsa in ms_dsa_by_run.values()])
    return _pooled_scatter(pooled, groups, dsa_column, title, ylabel, ax)


def _pooled_scatter(pooled: pd.DataFrame, groups: dict, dsa_column: str, title: str, ylabel: str, ax: plt.Axes|None):
    if ax is None:
        _, ax = plt.subplots()
    handles = {}
    for name, (conditions, column) in reversed(list(groups.items())):
        df = group_rows(pooled, conditions)
        df = df[df[column].notna()]
        handles[name] = ax.scatter(df[dsa_column], df[column], s=14 if name.startswith("No contact") else 7, alpha=0.4,
                                   color=CLASS_GROUP_COLORS.get(name), label=f"{name} ({len(df)})")
    ax.legend(handles=[handles[name] for name in groups], fontsize="x-small")
    ax.set_xlabel(dsa_column)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return ax


def run_averaged_histogram(ms_dsa_by_run: dict[int, pd.DataFrame], groups: dict, bins, dsa_column: str = "MS_DSA_shared_charges",
                           exclude: list[str] = (), title: str = "", xlabel: str = "", show_spread: bool = True,
                           ax: plt.Axes|None = None):
    """
    For every run, the percent of each group's peptides per bin (normalized within the group), then the mean over
    the runs (runs where the group has no peptide are skipped); error bars are the standard deviation over runs.
    The legend gives the number of peptides summed over the runs.
    """
    bins = np.asarray(bins, dtype=float)
    labels = bin_labels(bins)
    x = np.arange(len(labels))
    width = 0.8 / len(groups)
    if ax is None:
        _, ax = plt.subplots()
    for i, (name, (conditions, column)) in enumerate(groups.items()):
        percents, total = [], 0
        for ms_dsa in ms_dsa_by_run.values():
            values = group_rows(select_group_rows(ms_dsa, exclude, require=dsa_column), conditions)[column].dropna().values
            counts, _ = np.histogram(values, bins=bins)
            if counts.sum():
                percents.append(100.0 * counts / counts.sum())
                total += counts.sum()
        mean = np.mean(percents, axis=0) if percents else np.zeros(len(labels))
        spread = np.std(percents, axis=0) if len(percents) > 1 and show_spread else None
        ax.bar(x + (i - (len(groups) - 1) / 2) * width, mean, width, yerr=spread, capsize=1.5,
               error_kw={"elinewidth": 0.6}, label=f"{name} ({total})", color=CLASS_GROUP_COLORS.get(name), alpha=0.8)
    ax.set_xticks(x, labels, rotation=30 if len(labels) > 3 else 0, ha="right" if len(labels) > 3 else "center")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Percent of peptides (mean over runs)")
    ax.set_title(title)
    ax.legend(fontsize="x-small")
    return ax


def pooled_correlations(ms_dsa_by_run: dict[int, pd.DataFrame], groups: dict, dsa_column: str = "MS_DSA_shared_charges",
                        exclude: list[str] = (), method: str = "spearman") -> pd.Series:
    pooled = pd.concat([select_group_rows(ms_dsa, exclude, require=dsa_column) for ms_dsa in ms_dsa_by_run.values()])
    result = {}
    for name, (conditions, column) in groups.items():
        pair = group_rows(pooled, conditions)[[dsa_column, column]].dropna()
        result[name] = round(pair.iloc[:, 0].corr(pair.iloc[:, 1], method=method), 3) if len(pair) > 2 else np.nan
    return pd.Series(result, name=method)


def contact_class_summary_figure(ms_dsa_by_run: dict[int, pd.DataFrame], contact_class: str, sasa_source: str,
                                 purely_intra: bool, exclude: list[str], sasa_bins, dsa_bins, bfactor_bins,
                                 dsa_column: str = "MS_DSA_shared_charges", scatter_statistic: str = "max",
                                 suptitle: str = "", panel_size: tuple[float, float] = (4.2, 3.2)) -> plt.Figure:
    """
    All runs combined, 3 x 3 panels:
        row 1: run-averaged asa_nz histograms for max, min and avg
        row 2: pooled scatter of dsa_column vs asa_nz (scatter_statistic) and the run-averaged dsa_column histogram
        row 3: run-averaged B-factor histograms (NZ max of max, NZ min of min, residue average)
    """
    fig, axes = plt.subplots(3, 3, figsize=(panel_size[0] * 3, panel_size[1] * 3), constrained_layout=True)
    common = dict(dsa_column=dsa_column, exclude=exclude)
    for ax, statistic in zip(axes[0], ("max", "min", "avg")):
        run_averaged_histogram(ms_dsa_by_run, source_groups(contact_class, sasa_source, statistic, purely_intra), sasa_bins,
                               title=f"asa_nz {statistic}", xlabel="asa_nz", ax=ax, **common)
    pooled_scatter(ms_dsa_by_run, source_groups(contact_class, sasa_source, scatter_statistic, purely_intra),
                   title=f"{dsa_column} vs asa_nz {scatter_statistic} (all runs)", ax=axes[1][0], **common)
    run_averaged_histogram(ms_dsa_by_run, source_groups(contact_class, sasa_source, purely_intra=purely_intra, value_column=dsa_column),
                           dsa_bins, title=dsa_column, xlabel=dsa_column, ax=axes[1][1], **common)
    axes[1][2].set_visible(False)
    for ax, (statistic, column) in zip(axes[2], BFACTOR_STATISTIC_COLUMNS.items()):
        run_averaged_histogram(ms_dsa_by_run, source_groups(contact_class, sasa_source, purely_intra=purely_intra, value_column=column),
                               bfactor_bins, title=f"B-factor {column}", xlabel="z-score", ax=ax, **common)
    fig.suptitle(suptitle)
    return fig


def contact_venn(ms_dsa_by_run: dict[int, pd.DataFrame], contact_class: str, sasa_source: str,
                 dsa_column: str = "MS_DSA_shared_charges", exclude: list[str] = (), unit: str = "peptide",
                 title: str = "", ax: plt.Axes|None = None):
    """
    3-set Venn diagram of the groups of source_groups(contact_class, sasa_source) (membership only, no value needed),
    pooled over the runs. unit="peptide": unique (Uniprot, peptide_seq) with a dsa_column value in some run
    (the contact flags do not depend on the run); unit="peptide-run": (Uniprot, peptide_seq, run) pairs.
    """
    from matplotlib_venn import venn3
    frames = []
    for run_index, ms_dsa in ms_dsa_by_run.items():
        frames.append(select_group_rows(ms_dsa, exclude, require=dsa_column).assign(_run=run_index))
    pooled = pd.concat(frames)
    key = ["Uniprot", "peptide_seq"] + (["_run"] if unit == "peptide-run" else [])
    groups = source_groups(contact_class, sasa_source)
    sets = {name: set(map(tuple, group_rows(pooled, conditions)[key].values)) for name, (conditions, _) in groups.items()}
    everything = set(map(tuple, pooled[key].values))
    in_none = len(everything - set().union(*sets.values()))
    if ax is None:
        _, ax = plt.subplots()
    venn = venn3(list(sets.values()), set_labels=[f"{name}\n({len(s)})" for name, s in sets.items()],
                 set_colors=[CLASS_GROUP_COLORS[name] for name in sets], alpha=0.5, ax=ax)
    for text in venn.subset_labels:
        if text is not None:
            text.set_fontsize(7)
    ax.set_title(f"{title}\n{len(everything)} {unit}s, {in_none} in none of the three")
    return sets
