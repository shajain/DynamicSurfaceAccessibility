import gemmi
class BondLibrary:
    _REGISTRY = {}  # Global dictionary sitting in memory
    default_method = "assume_lysine"

    @classmethod
    def register(cls, bond_type: str, method="assume_lysine"):
        # STEP A: This outer function runs when Python reads the @ line.
        # It creates a closure that remembers `source_type` and `interface`.
        def decorator(check_bond):
            # STEP B: This inner function receives the actual class object.
            # It mutates the dictionary inside StructureFactory.
            cls._REGISTRY.setdefault(bond_type.lower(), {})[method] = check_bond
            # STEP C: It returns the subclass unmodified so the rest 
            # of your code can use PDBStructure normally.
            return check_bond
        return decorator

    @classmethod
    def check_bond(cls, residue1: gemmi.Residue, residue2: gemmi.Residue, bond_type: str, method: str=None):
        method = method or cls.default_method
        check_bond = cls._REGISTRY[bond_type.lower()][method]
        return check_bond(residue1, residue2)

    @classmethod
    def bond_types_implemented(cls):
        return list(cls._REGISTRY.keys())

    @classmethod
    def methods_implemented(cls, bond_type: str):
        return list(cls._REGISTRY[bond_type.lower()].keys())
    