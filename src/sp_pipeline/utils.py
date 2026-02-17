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


# ---------------------------------------------------------------------------
# Biophysical feature calculation (Kyte-Doolittle scale)
# ---------------------------------------------------------------------------

_KYTE_DOOLITTLE = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}

_CHARGE_PH7 = {
    "K": +1.0, "R": +1.0, "H": +0.1,
    "D": -1.0, "E": -1.0,
}


def _hydrophobicity_profile(sequence: str, window: int = 3) -> list[float]:
    seq = sequence.upper()
    if len(seq) < window:
        return [sum(_KYTE_DOOLITTLE.get(aa, 0.0) for aa in seq) / max(len(seq), 1)]
    return [
        sum(_KYTE_DOOLITTLE.get(seq[i + j], 0.0) for j in range(window)) / window
        for i in range(len(seq) - window + 1)
    ]


def compute_sp_features(sp_sequence: str) -> dict:
    """Compute biophysical features for a signal peptide sequence.

    Returns dict with hydrophobicity_mean, net_charge_ph7,
    n_region (n-terminal charged), h_region (hydrophobic core), c_region (cleavage end).
    """
    if not sp_sequence:
        return {"hydrophobicity_mean": None, "net_charge_ph7": None,
                "n_region": None, "h_region": None, "c_region": None}

    seq = sp_sequence.upper()
    n = len(seq)
    hydro_mean = round(sum(_KYTE_DOOLITTLE.get(aa, 0.0) for aa in seq) / n, 3)
    charge = round(sum(_CHARGE_PH7.get(aa, 0.0) for aa in seq), 2)

    profile = _hydrophobicity_profile(seq, window=3)

    h_start = min(3, n // 3)
    for i, val in enumerate(profile):
        if val > 1.0:
            h_start = i
            break

    h_end = h_start
    for i in range(len(profile) - 1, h_start - 1, -1):
        if profile[i] > 1.0:
            h_end = min(i + 2, n - 2)
            break

    c_start = min(h_end + 1, n - 1)

    return {
        "hydrophobicity_mean": hydro_mean,
        "net_charge_ph7": charge,
        "n_region": seq[:max(1, h_start)] or None,
        "h_region": seq[h_start : h_end + 1] or None,
        "c_region": seq[c_start:] or None,
    }


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
