class ClassFactory:
    _REGISTRY = {}  # Global dictionary sitting in memory

    @classmethod
    def register(cls, key: str, interface: type):
        # STEP A: This outer function runs when Python reads the @ line.
        # It creates a closure that remembers `source_type` and `interface`.
        def decorator(subclass):
            # STEP B: This inner function receives the actual class object.
            # It mutates the dictionary inside StructureFactory.
            cls._REGISTRY.setdefault(key.lower(), {})[interface] = subclass
            # STEP C: It returns the subclass unmodified so the rest 
            # of your code can use PDBStructure normally.
            return subclass
        return decorator

    @classmethod
    def create(cls, interface: type, key: str, *args, **kwargs):
        concrete_cls = cls._REGISTRY[key.lower()][interface]
        return concrete_cls(*args, **kwargs)