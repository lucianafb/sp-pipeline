"""Deduplication and idempotency for SP-Pipeline.

Ensures that:
1. Records from multiple sources don't create duplicates
2. Re-running the pipeline appends without duplicating
"""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from sp_pipeline.models import SPRecord, CSV_COLUMNS

logger = logging.getLogger("sp_pipeline.dedup")


def deduplicate_records(records: list[SPRecord]) -> list[SPRecord]:
    """Remove duplicate records based on unique_key.

    When duplicates exist, preference is given to:
    1. Experimental evidence over predicted
    2. UniProt over GenBank/NCBI

    Args:
        records: List of SPRecord instances.

    Returns:
        Deduplicated list.
    """
    seen: dict[str, SPRecord] = {}
    evidence_priority = {
        "experimental": 4,
        "curated_annotation": 3,
        "predicted_signalp": 2,
        "predicted_phobius": 2,
        "predicted_automatic": 1,
    }
    source_priority = {
        "UniProt": 3,
        "GenBank": 2,
        "NCBI": 1,
    }

    for record in records:
        key = record.unique_key
        if key not in seen:
            seen[key] = record
        else:
            existing = seen[key]
            # Compare priority: prefer experimental, then UniProt
            existing_score = (
                evidence_priority.get(existing.evidence_method, 0)
                + source_priority.get(existing.source_db, 0)
            )
            new_score = (
                evidence_priority.get(record.evidence_method, 0)
                + source_priority.get(record.source_db, 0)
            )
            if new_score > existing_score:
                seen[key] = record

    deduped = list(seen.values())
    removed = len(records) - len(deduped)
    if removed > 0:
        logger.info(f"Deduplication: removed {removed} duplicates, {len(deduped)} unique records")

    return deduped


def merge_with_existing(
    new_records: list[SPRecord],
    existing_csv: Optional[str] = None,
) -> pd.DataFrame:
    """Merge new records with an existing CSV file (idempotent append).

    Args:
        new_records: New SPRecord instances.
        existing_csv: Path to existing CSV file to merge with.

    Returns:
        Merged DataFrame ready for export.
    """
    new_df = pd.DataFrame([r.to_dict() for r in new_records])

    if existing_csv and Path(existing_csv).exists():
        logger.info(f"Merging with existing file: {existing_csv}")
        existing_df = pd.read_csv(existing_csv)

        # Create unique keys for both
        def make_key(row):
            return f"{row.get('source_db', '')}|{row.get('entry_id', '')}|{row.get('sp_start', '')}|{row.get('sp_end', '')}"

        existing_df["_key"] = existing_df.apply(make_key, axis=1)
        new_df["_key"] = new_df.apply(make_key, axis=1)

        # Find truly new records
        existing_keys = set(existing_df["_key"])
        new_only = new_df[~new_df["_key"].isin(existing_keys)]

        logger.info(f"Existing: {len(existing_df)}, New unique: {len(new_only)}")

        # Combine
        merged = pd.concat([existing_df, new_only], ignore_index=True)
        merged = merged.drop(columns=["_key"])
    else:
        merged = new_df

    # Ensure column order
    for col in CSV_COLUMNS:
        if col not in merged.columns:
            merged[col] = None
    merged = merged[CSV_COLUMNS]

    return merged
