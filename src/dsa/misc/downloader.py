import requests


import requests
from requests.adapters import HTTPAdapter, Retry
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional


class Downloader(ABC):
    """
    Base class for downloading objects identified by an accession number
    from a remote service, over a persistent, retry-enabled HTTP session.

    Subclasses must implement `build_url()` to map an accession to a
    full download URL for that specific service.
    """

    def __init__(
        self,
        base_url: str,
        max_retries: int = 3,
        backoff_factor: float = 0.5,
        timeout: float = 10.0,
        headers: Optional[dict] = None,
    ):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.session = self._build_session(max_retries, backoff_factor, headers)

    def _build_session(
        self,
        max_retries: int,
        backoff_factor: float,
        headers: Optional[dict],
    ) -> requests.Session:
        session = requests.Session()
        retries = Retry(
            total=max_retries,
            backoff_factor=backoff_factor,
            status_forcelist=[500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retries)
        session.mount("https://", adapter)
        session.mount("http://", adapter)

        if headers:
            session.headers.update(headers)

        return session

    @abstractmethod
    def build_url(self, accession: str) -> str:
        """
        Construct the exact URL to download the object of interest
        for a given accession. Must be implemented by subclasses.
        """
        raise NotImplementedError

    def download(self, accession: str) -> Optional[bytes]:
        """
        Download the raw content for a given accession.
        Returns None (and warns) on failure rather than raising,
        so batch/loop callers don't need to wrap every call in try/except.
        """
        url = self.build_url(accession)
        response = None
        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"Error downloading {accession} from {url}: {e}")
            return None
        data = self.process_response(response, accession)
        return data

    def process_response(self, response, accession: str):
        """
        Process the response content.
        """
        return response.content
            

    # def download_to_file(self, accession: str, output_path: Path) -> bool:
    #     """
    #     Download and save raw content directly to disk.
    #     Returns True on success, False on failure.
    #     """
    #     content = self.download(accession)
    #     if content is None:
    #         return False

    #     output_path = Path(output_path)
    #     output_path.parent.mkdir(parents=True, exist_ok=True)
    #     output_path.write_bytes(content)
    #     return True

    def close(self):
        """Close the underlying session's connection pool."""
        self.session.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()    