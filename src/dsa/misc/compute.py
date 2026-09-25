import random
from typing import Callable
from dsa.misc.store import Store

class ComputeManager:
    def __init__(self, keys_to_process: list[str], 
                       compute_function: Callable, 
                       store: Store, 
                       update_interval: int = 25,
                       show_progress: bool = True):
        self.keys_to_process = keys_to_process
        self.compute_function = compute_function
        self.store = store
        self.successful_keys = []
        self.unsuccessful_keys = []
        self.saved_keys = self.get_saved_keys()
        self.unprocessed_keys = self.get_unprocessed_keys()
        self.update_interval = update_interval
        self.show_progress = show_progress


    def compute(self):
        iter_count = 0
        if self.show_progress:
            print(self.progress_message())
        while self.unprocessed_keys_count > 0:
            key = self.get_random_unprocessed_key()
            data = self.compute_function(key)
            iter_count += 1
            if iter_count % self.update_interval == 0:
                self.update_keys()
                if self.show_progress:
                    print(self.progress_message(key, data))
            if data is None:
                self.unsuccessful_keys.append(key)
                continue
            self.store.dump(key, data)
            self.successful_keys.append(key)
        
    def progress_message(self, key: str=None, data: any=None) -> str:
        msg = (f"Processed {self.processed_keys_count} of {self.keys_to_process_count} keys")
        msg += f"\n successful: {self.successful_keys_count}, unsuccessful: {self.unsuccessful_keys_count}"
        msg += f"\n current key: {key}" if key is not None else ""
        if hasattr(data, "summary"):
            msg += f"\n data: {data.summary()}"
        else:
            msg += f"\n data: {data}"
        return msg


    def get_random_unprocessed_key(self) -> str:
        return random.choice(self.unprocessed_keys)

    def update_keys(self):
        self.saved_keys = self.get_saved_keys()
        self.unprocessed_keys = self.get_unprocessed_keys()


    def get_saved_keys(self) -> list[str]:
        return self.store.get_saved_keys(from_keys=self.keys_to_process)

    def get_unprocessed_keys(self) -> list[str]:
        unprocessed_keys = [key for key in self.keys_to_process if key not in self.saved_keys]
        unprocessed_keys = unprocessed_keys + self.unsuccessful_keys
        return unprocessed_keys


    @property
    def keys_to_process_count(self) -> int:
        return len(self.keys_to_process)

    @property
    def processed_keys_count(self) -> int:
        return len(self.processed_keys)

    @property
    def unprocessed_keys_count(self) -> int:
        return len(self.unprocessed_keys)

    @property
    def saved_keys_count(self) -> int:
        return len(self.saved_keys)

    @property
    def successful_keys_count(self) -> int:
        return len(self.successful_keys)

    @property
    def unsuccessful_keys_count(self) -> int:
        return len(self.unsuccessful_keys)

    @property
    def processed_keys(self) -> list[str]:
        return self.successful_keys + self.unsuccessful_keys
