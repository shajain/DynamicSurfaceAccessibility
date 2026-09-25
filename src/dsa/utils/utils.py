from collections import defaultdict

def to_defaultdict_list(d):
    """Recursively convert nested dicts back to defaultdicts of nested lists and dicts"""
    if isinstance(d, dict):
        return defaultdict(dict2defaultdict(type(d.values()[0])), {k: to_defaultdict_list(v)})
            elif isinstance(v, list):
                return defaultdict(list, {k: to_defaultdict_list(v)})
    elif isinstance(d, list):
        return [to_defaultdict_list(i) for i in d]
    return d
    
def dict2defaultdict(typee):
    return defaultdict if typee == dict else typee


def to_nested_defaultdict(d, depth=3):
    """
    Converts nested dict to defaultdict structure:
    depth 3: str -> defaultdict(str -> defaultdict(str -> list))
    """
    if not isinstance(d, dict):
        return d  # base case — list or primitive value
    
    if depth == 1:
        # innermost level — values are lists
        return defaultdict(list, {k: v for k, v in d.items()})
    else:
        # outer levels — values are dicts
        return defaultdict(dict, {
            k: to_nested_defaultdict(v, depth - 1) 
            for k, v in d.items()
        })