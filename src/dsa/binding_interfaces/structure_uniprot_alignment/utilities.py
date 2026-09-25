import re
import gemmi

def extract_num_str(s) -> str | None:
    m = re.search(r'-?\d+', str(s))
    return m.group() if m else None

def single_letter_resname(res_name: str) -> str:
    if len(res_name) == 3:
        info = gemmi.find_tabulated_residue(res_name)
        return info.one_letter_code.upper() if info is not None else 'X'
    elif len(res_name) == 1:
        return res_name.upper()
    else:
        return 'X'

def int_if_str(s) -> int | None:
    if isinstance(s, str):
        return int(s)
    elif isinstance(s, int):
        return s
    else:
        return None

def str2int(s) -> int | None:
    if isinstance(s, int):
        return s
    elif isinstance(s, str):
        num_str = extract_num_str(s)
        if num_str is not None:
            return int(num_str)
        else:
            return None
    ValueError(f"Invalid input: {s}")