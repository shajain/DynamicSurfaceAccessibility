import pickle
from pathlib import Path
from abc import ABC, abstractmethod
import warnings
import gemmi
import zstandard as zstd
import json
import pandas as pd

class Store(ABC):
    def __init__(self, directory: Path|str, extension: str="", object_type: type = None):
        self._extensions = extension.split(".")
        if isinstance(directory, str):
            directory = Path(directory)
        self.directory = directory
        self.directory.mkdir(parents=True, exist_ok=True)
        self.extension = extension
        self.type = object_type

    def key_to_file_name(self, key: str) -> str:
        return f"{key}.{self.extension}"

    def file_name_to_key(self, file_name: str) -> str:
        return file_name.split(".")[0]

    def get_file_path(self, key: str=None) -> Path:
        return self.directory / self.key_to_file_name(key)
    
    def dump(self, key: str, data: any):
        file_path = self.get_file_path(key)
        if self.type is not None and type(data) == self.type:
            data = self.serialize(data)
        self.dump_file(file_path, data)
    
    @abstractmethod
    def dump_file(self, file_path: Path, data: any):
        pass

    @abstractmethod
    def load_file(self, file_path: Path) -> any:
        pass

    def deserialize(self, data: any) -> any:
        return data

    def serialize(self, data: any) -> any:
        return data

    # def dump_from_object(self, key: str, data: any):
    #     data = self.serialize(data)
    #     self.dump(key, data)
    
    # def load_as_object(self, key: str) -> any:
    #     data = self.load(key)
    #     return self.deserialize(data)

    def load(self, key: str):
        file_path = self.get_file_path(key)
        #assert file_path.exists(), f"File {file_path} does not exist."
        data = None
        try:
            data = self.load_file(file_path)
        except Exception as e:
            warnings.warn(f"Error loading file {file_path}: {e}")
        if data is not None and self.type and type(data) != self.type:
            data = self.deserialize(data)
        return data

    def update_saved_keys(self):
        self.saved_keys = self.get_saved_keys()

    def get_saved_keys(self, from_keys: list[str] = None) -> list[str]:
        if from_keys is None:
            saved_keys = [self.file_name_to_key(file.name) for file in self.directory.glob(f"*.{self.extension}")]
        else:
            saved_keys = [self.file_name_to_key(file.name) for file in self.directory.glob(f"*.{self.extension}")]
            saved_keys = list(set(saved_keys) & set(from_keys))
        return list(set(saved_keys))

    def update_extensions(self, sub_extension: str):
        if sub_extension not in self._extensions:
            self._extensions.append(sub_extension)
            _extensions = [e for e in self._extensions if e.strip()]
            self.extension = ".".join(_extensions)

    
class PickleStore(Store):
    def __init__(self, directory: Path|str, extension: str|list = "", object_type: type = None):
        super().__init__(directory, extension=extension, object_type=object_type)
        self.update_extensions("pkl")

    def dump_file(self, file_path: Path, data: any):
        with open(file_path, "wb") as f:
            pickle.dump(data, f)

    def load_file(self, file_path: Path) -> any:
        with open(file_path, "rb") as f:
            return pickle.load(f)

    def serialize(self, data: any) -> any:
        return data
    
    def deserialize(self, data: any) -> any:
        return data

class CompressedPickleStore(PickleStore):
    def __init__(self, directory: Path|str, extension: str = "", object_type: type = None):
        super().__init__(directory, extension=extension, object_type=object_type)
        self.update_extensions("zst")
        self.compressor = zstd.ZstdCompressor(level=3)
        self.decompressor = zstd.ZstdDecompressor()

    def dump_file(self, file_path: Path, data: any):
        with open(file_path, "wb") as f:
            with self.compressor.stream_writer(f) as writer:
                pickle.dump(data, writer, protocol=pickle.HIGHEST_PROTOCOL)
            
    def load_file(self, file_path: Path) -> any:
        with open(file_path, "rb") as f:
            with self.decompressor.stream_reader(f) as reader:
                return pickle.load(reader)


class JSONStore(Store): 
    def __init__(self, directory: Path|str, extension: str = "", object_type: type = None):
        super().__init__(directory, extension=extension, object_type=object_type)

    def dump_file(self, file_path: Path, data: any):
        with open(file_path, "w") as f:
            json.dump(data, f)

    def load_file(self, file_path: Path) -> any:
        with open(file_path, "r") as f:
            return json.load(f)

    def serialize(self, data: any) -> any:
        return data.to_dict()

    def deserialize(self, data: any) -> any:
        return self.type.from_dict(data)


class SingleFileStore(Store):    
    def __init__(self, directory: Path|str, file_name: str, object_type: type = None):
        # get extension from file_name
        extension = ".".join(file_name.split(".")[1:])
        self.file_name = file_name
        super().__init__(directory, extension, object_type=object_type)
        self.data = self.load_file(self.get_file_path())

    def key_to_file_name(self, key: str=None) -> str:
        NotImplementedError("key_to_file_name is not implemented for SingleFileStore")

    def file_name_to_key(self, file_name: str) -> str:
        NotImplementedError("file_name_to_key is not implemented for SingleFileStore")

    def get_file_path(self) -> Path:
        return self.directory / self.file_name

    

class SingleFileJsonStore(SingleFileStore):
    def __init__(self, directory: Path|str, file_name: str, object_type: type = None):
        super().__init__(directory, file_name, object_type)

    def dump_file(self, file_path: Path, data: any):
        json.dump(data, file_path)

    def load_file(self, file_path: Path) -> any:
        data = json.load(file_path)
        return data

    def deserialize(self, data: any) -> any:
        return self.type.from_dict(data)

    def serialize(self, data: any) -> any:
        return data.to_dict()


class SingleFileTabularStore(SingleFileStore):
    def __init__(self, directory: Path|str, file_name: str, dtype: dict[str, type]|type = str):
        self.file_format, self.separator = self._extract_format_and_separator(file_name)
        object_type = pd.DataFrame
        self.dtype = dtype
        super().__init__(directory, file_name, object_type)

    def _extract_format_and_separator(self, file_name: str) -> str:
        if "tsv" in file_name.split("."):
            file_format = "tsv"
            separator = '\t'
        elif "csv" in file_name.split("."):
            file_format = "csv"
            separator = ","
        else:
            raise ValueError(f"Invalid file name: {file_name}")
        return file_format, separator
    
    def dump_file(self, file_path: Path, data: any):
        data.to_csv(file_path, sep=self.separator, dtype=self.dtype)

    def load_file(self, file_path: Path) -> any:
        data = pd.read_csv(file_path, sep=self.separator, dtype=self.dtype, comment="#")
        return data