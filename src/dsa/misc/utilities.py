from collections import defaultdict
import pandas as pd


def defaultdict_to_dict(obj):
    if isinstance(obj, defaultdict):
        return {k: defaultdict_to_dict(v) for k, v in obj.items()}
    return obj




def nested_dict_to_table(obj, value_fn=len, value_name="count", column_level=None):
    """
    obj: nested dict; leaves are lists/sets/tuples
    value_fn: applied to each leaf (default: len)
    column_level: path index to use as columns (e.g. 0 -> inter-chain / intra-chain).
                  None keeps every key level as a MultiIndex row label.
    """
    rows = []

    def walk(x, path):
        if isinstance(x, dict):
            for k, v in x.items():
                walk(v, path + (k,))
        else:
            rows.append({"path": path, value_name: value_fn(x)})

    walk(obj, ())
    if not rows:
        return pd.DataFrame(columns=[value_name])

    depth = max(len(r["path"]) for r in rows)
    records = []
    for r in rows:
        rec = {f"level_{i}": (r["path"][i] if i < len(r["path"]) else None)
               for i in range(depth)}
        rec[value_name] = r[value_name]
        records.append(rec)

    df = pd.DataFrame.from_records(records)
    level_cols = [f"level_{i}" for i in range(depth)]

    if column_level is None:
        return df.set_index(level_cols)

    col = level_cols[column_level]
    idx = [c for c in level_cols if c != col]
    return df.pivot_table(index=idx, columns=col, values=value_name, aggfunc="sum")