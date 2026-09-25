import numbers
import gemmi
import pandas as pd

def single_letter_residue_code(residue: gemmi.Residue|str) -> str:
    if isinstance(residue, gemmi.Residue):
        residue_name = residue.name
    elif isinstance(residue, str):
        residue_name = residue
    else:
        raise ValueError(f"Invalid residue type: {type(residue)}")
    res_info = gemmi.find_tabulated_residue(residue_name)
    return res_info.one_letter_code.upper()


def is_actually_a_number(val):
    # This catches pd.NA, None, AND np.nan all at once!
    if pd.isna(val): 
        return False
        
    return isinstance(val, numbers.Number)