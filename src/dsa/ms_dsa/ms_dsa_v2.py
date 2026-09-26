"""
MS_DSA per (peptide, run) from light/heavy dimethyl-labelled lysine peptides.

MS_DSA = light / (light + heavy), where light and heavy are summed Ms1.Area values, computed two ways:

- MS_DSA                 : all charge states of each channel are summed (same as ms_dsa.MS_DSA_Builder).
- MS_DSA_shared_charges  : only charge states measured in both channels are summed, so a channel seen
                           at extra charge states does not dominate. Not computed (NaN) when a channel
                           is absent or when no charge state is shared.

Everything is computed separately for each run. For MS_DSA a channel with no rows counts as area 0, so
MS_DSA is 0 or 1 when one channel is absent. Either value is NaN when both of its areas are 0 (0 / 0).
Such cases are flagged instead of removed:

- channel_absent / absent_channel       : no row at all for the light or the heavy channel.
- no_shared_charges                     : both channels have rows, but no charge state in common.
- zero_area / zero_area_channel         : both channels have rows, but a channel's summed area is 0
                                          (all charges).
- zero_area_shared / zero_area_channel_shared : a shared-charges area is 0 (only where
                                          MS_DSA_shared_charges is computed).
"""
import numpy as np
import pandas as pd
from pathlib import Path
from dsa.ms_dsa.config import MS_DATA_DIR, ms_data_relativePath

LIGHT_CHANNEL = '0'
HEAVY_CHANNEL = '8'
CHANNEL_NAMES = {LIGHT_CHANNEL: "light", HEAVY_CHANNEL: "heavy"}
GROUP_KEYS = ["Stripped.Sequence", "Run.Index"]
METADATA_COLUMNS = ["Run", "Protein.Names", "Genes", "Protein.Ids", "Protein.Group"]


class MS_DSA_Builder:
    def __init__(self, ms_data_relativePath: str|Path = ms_data_relativePath, rebuild: bool = False):
        self.ms_data_relativePath = Path(ms_data_relativePath)
        self.ms_path = MS_DATA_DIR / self.ms_data_relativePath
        assert self.ms_path.exists(), f"MS data file not found at {self.ms_path}"
        self.ms_dsa_path = self.ms_path.with_name(f"{self.ms_path.stem}_ms_dsa_v2.csv")
        self.ms = pd.read_parquet(self.ms_path)
        if self.ms_dsa_path.exists() and not rebuild:
            self.ms_dsa = self.load_ms_dsa()
        else:
            self.ms_dsa = self.build_ms_dsa()
            self.save_ms_dsa()

    def load_ms_dsa(self) -> pd.DataFrame:
        return pd.read_csv(self.ms_dsa_path)

    def save_ms_dsa(self):
        self.ms_dsa.to_csv(self.ms_dsa_path, index=False)

    def build_ms_dsa(self) -> pd.DataFrame:
        ms = self.keep_sequences_with_lysine(self.ms)
        ms = ms[ms["Channel"].isin(CHANNEL_NAMES)]
        metadata = self.group_metadata(ms)
        per_charge = self.areas_per_charge(ms)
        areas = self.channel_areas(per_charge)
        ms_dsa = metadata.join(areas)
        ms_dsa["MS_DSA"] = compute_dsa(ms_dsa["light_area"], ms_dsa["heavy_area"])
        ms_dsa["MS_DSA_shared_charges"] = compute_dsa(ms_dsa["light_area_shared"], ms_dsa["heavy_area_shared"])
        self.add_flags(ms_dsa)
        ms_dsa = ms_dsa.reset_index().rename(columns={"Stripped.Sequence": "peptide_seq"})
        return ms_dsa[self.output_columns()]

    @staticmethod
    def group_metadata(ms: pd.DataFrame) -> pd.DataFrame:
        grouped = ms.groupby(GROUP_KEYS)[METADATA_COLUMNS]
        inconsistent = grouped.nunique().gt(1).any()
        if inconsistent.any():
            raise ValueError(f"Columns {inconsistent[inconsistent].index.tolist()} differ within a (peptide, run) group")
        return grouped.first()

    @staticmethod
    def areas_per_charge(ms: pd.DataFrame) -> pd.DataFrame:
        # one row per (peptide, run, charge), columns light / heavy; NaN when the channel has no row.
        # Ms1.Area is float32 in the report; sum in float64
        areas = ms["Ms1.Area"].astype("float64")
        per_charge = (areas.groupby([ms[k] for k in GROUP_KEYS + ["Precursor.Charge", "Channel"]]).sum()
                        .unstack("Channel")
                        .reindex(columns=list(CHANNEL_NAMES))
                        .rename(columns=CHANNEL_NAMES))
        return per_charge

    @staticmethod
    def channel_areas(per_charge: pd.DataFrame) -> pd.DataFrame:
        by_group = per_charge.groupby(level=GROUP_KEYS)
        shared = per_charge.notna().all(axis=1)
        shared_areas = per_charge[shared].groupby(level=GROUP_KEYS).sum()
        areas = pd.DataFrame({
            "n_light_charges":  by_group["light"].count(),
            "n_heavy_charges":  by_group["heavy"].count(),
            "n_shared_charges": shared.groupby(level=GROUP_KEYS).sum(),
            "light_area":       by_group["light"].sum(min_count=1),   # NaN when the channel is absent
            "heavy_area":       by_group["heavy"].sum(min_count=1),
        })
        for channel in ["light", "heavy"]:
            # NaN when no charge state is shared (which includes every group with a channel absent)
            areas[f"{channel}_area_shared"] = shared_areas[channel].reindex(areas.index)
        return areas

    @staticmethod
    def add_flags(ms_dsa: pd.DataFrame):
        light_absent = ms_dsa["n_light_charges"] == 0
        heavy_absent = ms_dsa["n_heavy_charges"] == 0
        ms_dsa["channel_absent"] = light_absent | heavy_absent
        ms_dsa["absent_channel"] = np.select([light_absent, heavy_absent], ["light", "heavy"], default=None)
        ms_dsa["no_shared_charges"] = ~ms_dsa["channel_absent"] & (ms_dsa["n_shared_charges"] == 0)
        ms_dsa["multiple_charges"] = (ms_dsa["n_light_charges"] > 1) | (ms_dsa["n_heavy_charges"] > 1)
        computed = {"": ~ms_dsa["channel_absent"], "_shared": ms_dsa["n_shared_charges"] > 0}
        for suffix, flagged in computed.items():
            light_zero = flagged & (ms_dsa[f"light_area{suffix}"] == 0)
            heavy_zero = flagged & (ms_dsa[f"heavy_area{suffix}"] == 0)
            ms_dsa[f"zero_area{suffix}"] = light_zero | heavy_zero
            ms_dsa[f"zero_area_channel{suffix}"] = np.select(
                [light_zero & heavy_zero, light_zero, heavy_zero], ["both", "light", "heavy"], default=None)

    @staticmethod
    def output_columns() -> list[str]:
        return ["peptide_seq", "Run", "Run.Index", "MS_DSA", "MS_DSA_shared_charges",
                "Protein.Names", "Genes", "Protein.Ids", "Protein.Group",
                "channel_absent", "absent_channel", "no_shared_charges", "zero_area", "zero_area_channel",
                "zero_area_shared", "zero_area_channel_shared", "multiple_charges",
                "n_light_charges", "n_heavy_charges", "n_shared_charges",
                "light_area", "heavy_area", "light_area_shared", "heavy_area_shared"]

    def get_ms_dsa_df(self) -> pd.DataFrame:
        return self.ms_dsa

    def get_uniprot_ids(self) -> list[str]:
        return list({uniprot_id for ids in self.ms_dsa["Protein.Ids"] for uniprot_id in ids.split(";")})

    @staticmethod
    def keep_sequences_with_lysine(ms: pd.DataFrame) -> pd.DataFrame:
        return ms[ms["Stripped.Sequence"].str.strip().str.upper().str.contains("K")]

    @staticmethod
    def get_all_uniprot_ids(ms_data_dir: Path = MS_DATA_DIR) -> list[str]:
        uniprot_ids = set()
        for file in ms_data_dir.glob("**/*.parquet"):
            uniprot_ids.update(MS_DSA_Builder(file.relative_to(ms_data_dir)).get_uniprot_ids())
        return list(uniprot_ids)


def compute_dsa(light_area: pd.Series, heavy_area: pd.Series) -> pd.Series:
    # an absent channel (NaN area) counts as 0; 0 / 0 gives NaN
    light, heavy = light_area.fillna(0), heavy_area.fillna(0)
    total = light + heavy
    return (light / total).where(total > 0)


def summarize(ms_dsa: pd.DataFrame) -> str:
    n = len(ms_dsa)
    lines = [f"(peptide, run) groups: {n:,}"]
    lines.append(f"  no_shared_charges: {ms_dsa['no_shared_charges'].sum():,}")
    for column in ["absent_channel", "zero_area_channel", "zero_area_channel_shared"]:
        counts = ms_dsa[column].value_counts().to_dict()
        lines.append(f"  {column}: {counts} (any: {ms_dsa[column].notna().sum():,}, "
                     f"{ms_dsa[column].notna().mean():.1%})")
    for column in ["MS_DSA", "MS_DSA_shared_charges"]:
        dsa = ms_dsa[column]
        lines.append(f"  {column}: 0 -> {(dsa == 0).sum():,}, 1 -> {(dsa == 1).sum():,}, "
                     f"between -> {((dsa > 0) & (dsa < 1)).sum():,}, NaN -> {dsa.isna().sum():,}")
    return "\n".join(lines)


if __name__ == "__main__":
    builder = MS_DSA_Builder(ms_data_relativePath=ms_data_relativePath, rebuild=True)
    print(f"Saved {builder.ms_dsa_path}")
    print(summarize(builder.ms_dsa))
