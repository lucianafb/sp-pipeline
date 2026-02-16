"""Data source modules for SP-Pipeline."""

from abc import ABC, abstractmethod
from typing import Any

from sp_pipeline.models import SPRecord


class BaseSource(ABC):
    """Abstract base class for all data sources."""

    @property
    @abstractmethod
    def name(self) -> str:
        ...

    @abstractmethod
    def query(self, params: dict[str, Any]) -> list[SPRecord]:
        ...

    @abstractmethod
    def is_available(self) -> bool:
        ...
