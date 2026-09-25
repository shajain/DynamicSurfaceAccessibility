from pandas.core.dtypes.dtypes import str_type
from dsa.uniprot.config import ALL_UNIPROT_IDS
from dsa.uniprot.config import sequence_file
import urllib.request
import json
import pandas as pd



    
class UniprotToSequence:
    _instance: "UniprotToSequence | None" = None

    def __init__(self, uniprot_ids: list[str] = ALL_UNIPROT_IDS, uniprot_to_sequence_file: str = sequence_file):
        self.uniprot_ids = uniprot_ids
        self.sequence_file = uniprot_to_sequence_file
        self.sequences = self.load()
        self.update()

    @classmethod
    def _get_instance(cls, uniprot_to_sequence_file: str = sequence_file) -> "UniProtToSequence":
        if cls._instance is None:
            cls._instance = cls(uniprot_ids=ALL_UNIPROT_IDS, uniprot_to_sequence_file=uniprot_to_sequence_file)
        return cls._instance


    def update(self):
        any_new_sequence = False
        for uniprot_id in self.uniprot_ids:
            if uniprot_id not in self.sequences:
                sequence = self.fetch_from_api(uniprot_id)
                self.sequences[uniprot_id] = sequence
                any_new_sequence = True
        if any_new_sequence:
            self.save()

    def fetch_from_api(self, uniprot_id: str):
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"  
        with urllib.request.urlopen(url) as r:
            fasta = r.read().decode()
        lines = fasta.strip().split('\n')
        seq = ''.join(lines[1:])  # skip header line
        return seq

    def load(self):
        if not self.sequence_file.exists():
            print(f"Sequence file not found at {self.sequence_file}")
            print(f"Creating an empty sequence file at {self.sequence_file}")
            with open(self.sequence_file, "w") as f:
                json.dump({}, f)
            return  {}
        with open(self.sequence_file, "r") as f:
            return json.load(f)

    def save(self, uniprot_id: str=None, sequence: str=None):
        if uniprot_id is not None and sequence is not None:
            self.sequences[uniprot_id] = sequence
        with open(self.sequence_file, "w") as f:
            json.dump(self.sequences, f)


    def _get_sequence(self, uniprot_id: str, only_saved_sequences: bool = False):
        if uniprot_id not in self.sequences and not only_saved_sequences:
            sequence = self.fetch_from_api(uniprot_id)
            self.save(uniprot_id, sequence)
        return self.sequences.get(uniprot_id, None)
    
    def stored_sequences(self, format = dict | pd.DataFrame):
        if format == dict:
            return self.sequences
        elif format== pd.DataFrame:
            return pd.DataFrame(self.sequences.items(), columns=["uniprot_id", "sequence"])
        else:
            raise ValueError(f"Invalid format: {format}")


    @classmethod
    def get_sequence(cls, uniprot_id: list[str]|str, only_saved_sequences: bool = False) -> str | list[str]:
        if isinstance(uniprot_id, str):
            return UniprotToSequence._get_instance()._get_sequence(uniprot_id, only_saved_sequences=only_saved_sequences)
        elif isinstance(uniprot_id, list):
            return [UniprotToSequence._get_instance()._get_sequence(id, only_saved_sequences=only_saved_sequences) for id in uniprot_id]
        return None



if __name__ == "__main__":
    uniprot_to_sequence = UniprotToSequence()
    print(uniprot_to_sequence.stored_sequences(format=dict))
    print(uniprot_to_sequence.stored_sequences(format=pd.DataFrame))


