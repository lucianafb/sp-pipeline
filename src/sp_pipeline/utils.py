"""Utility functions for SP-Pipeline."""

import hashlib
import json
import logging
import time
from pathlib import Path
from datetime import datetime, timedelta
from typing import Any, Optional

from tqdm import tqdm


def setup_logging(level: str = "INFO", log_file: Optional[str] = None) -> logging.Logger:
    """Configure logging for the pipeline.

    Args:
        level: Logging level string.
        log_file: Optional path to log file.

    Returns:
        Configured logger.
    """
    logger = logging.getLogger("sp_pipeline")
    logger.setLevel(getattr(logging, level.upper()))

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger.addHandler(console)

    # File handler
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


class RateLimiter:
    """Simple rate limiter for API calls."""

    def __init__(self, calls_per_second: int = 3):
        self.min_interval = 1.0 / calls_per_second
        self.last_call = 0.0

    def wait(self):
        """Wait if necessary to respect rate limit."""
        elapsed = time.time() - self.last_call
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_call = time.time()


class QueryCache:
    """Simple file-based cache for API query results."""

    def __init__(self, cache_dir: str, ttl_days: int = 30, enabled: bool = True):
        self.cache_dir = Path(cache_dir)
        self.ttl = timedelta(days=ttl_days)
        self.enabled = enabled
        if self.enabled:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _key_to_path(self, key: str) -> Path:
        """Convert a cache key to a file path."""
        hashed = hashlib.sha256(key.encode()).hexdigest()[:16]
        return self.cache_dir / f"{hashed}.json"

    def get(self, key: str) -> Optional[Any]:
        """Retrieve a cached result.

        Args:
            key: Cache key (typically the query parameters serialized).

        Returns:
            Cached data or None if miss/expired.
        """
        if not self.enabled:
            return None

        path = self._key_to_path(key)
        if not path.exists():
            return None

        try:
            with open(path, "r") as f:
                cached = json.load(f)

            # Check expiration
            cached_time = datetime.fromisoformat(cached["timestamp"])
            if datetime.now() - cached_time > self.ttl:
                path.unlink()
                return None

            return cached["data"]
        except (json.JSONDecodeError, KeyError):
            path.unlink(missing_ok=True)
            return None

    def set(self, key: str, data: Any):
        """Store a result in cache.

        Args:
            key: Cache key.
            data: Data to cache (must be JSON-serializable).
        """
        if not self.enabled:
            return

        path = self._key_to_path(key)
        cache_entry = {
            "timestamp": datetime.now().isoformat(),
            "key": key,
            "data": data,
        }
        with open(path, "w") as f:
            json.dump(cache_entry, f)

    def clear(self):
        """Clear all cached data."""
        if self.cache_dir.exists():
            for f in self.cache_dir.glob("*.json"):
                f.unlink()


def make_cache_key(**kwargs) -> str:
    """Create a deterministic cache key from query parameters."""
    sorted_items = sorted(kwargs.items())
    return json.dumps(sorted_items, sort_keys=True)


def progress_bar(iterable, desc: str = "", total: Optional[int] = None):
    """Wrap an iterable with a progress bar."""
    return tqdm(iterable, desc=desc, total=total, unit="records")


def extract_cleavage_motif(full_sequence: str, cleavage_pos: int, window: int = 3) -> str:
    """Extract the cleavage site motif from a sequence.

    Extracts residues around the cleavage site in the format:
    [-3][-2][-1]-[+1][+2][+3]

    Args:
        full_sequence: Complete protein sequence.
        cleavage_pos: Position of cleavage (1-based, last residue of SP).
        window: Number of residues on each side.

    Returns:
        Motif string, e.g., "VFA-AP".
    """
    if not full_sequence or cleavage_pos <= 0:
        return ""

    idx = cleavage_pos - 1  # Convert to 0-based
    before = full_sequence[max(0, idx - window + 1) : idx + 1]
    after = full_sequence[idx + 1 : idx + 1 + window]

    return f"{before}-{after}" if after else before
