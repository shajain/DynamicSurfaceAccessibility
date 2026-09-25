class FileStrings:
    
    
    @staticmethod
    def stringify_cutoff(distance_cutoff: float) -> str:
        return f"{int(distance_cutoff * 10)}"
    
    @staticmethod
    def unstringify_cutoff(cutoff_string: str) -> float:
        return float(cutoff_string) / 10

    