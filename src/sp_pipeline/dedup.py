"""Deduplication and idempotency for SP-Pipeline.

Ensures that:
1. Records with identical SP sequences are collapsed into one representative record
   (all_accessions and duplicate_count track the merged entries)
2. Re-running the pipeline appends without duplicating
"""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from sp_pipeline.models import SPRecord, CSV_COLUMNS

logger = logging.getLogger("sp_pipeline.dedup")

# Evidence priority: higher = better (kept when collapsing duplicates)
_EVIDENCE_PRIORITY = {
    "experimental": 4,
    "curated_annotation": 3,
    "predicted_signalp": 2,
    "predicted_phobius": 2,
    "predicted_automatic": 1,
}

_SOURCE_PRIORITY = {
    "UniProt": 3,
    "GenBank": 2,
    "NCBI": 1,
}


def deduplicate_records(records: list[SPRecord]) -> list[SPRecord]:
    """Collapse records that share the exact same SP sequence and type.

    Groups by (sp_sequence, sp_type). Within each group:
    - Keeps the record with the best evidence (experimental > curated > predicted)
    - Merges all entry_ids into all_accessions (pipe-separated, sorted)
    - Sets duplicate_count to the number of entries in the group

    This correctly handles:
    - DENV-1/2/3/4 prM signal → 1 record with duplicate_count=4
    - Multiple alphavirus strains with identical E3 → 1 representative
    - Same protein from UniProt + NCBI → collapsed by sequence

    Args:
        records: List of SPRecord instances.

    Returns:
        Deduplicated list with all_accessions and duplicate_count set.
    """
    if not records:
        return records

    # Group by (sp_sequence, sp_type)
    groups: dict[tuple, list[SPRecord]] = {}
    for record in records:
        key = (record.sp_sequence, record.sp_type)
        groups.setdefault(key, []).append(record)

    deduped = []
    for (seq, sp_type), group in groups.items():
        # Sort by evidence + source priority (descending: best first)
        def score(r: SPRecord) -> int:
            return (
                _EVIDENCE_PRIORITY.get(r.evidence_method, 0)
                + _SOURCE_PRIORITY.get(r.source_db, 0)
            )

        group.sort(key=score, reverse=True)
        best = group[0]

        # Merge all accessions (unique, sorted)
        all_ids = sorted(set(r.entry_id for r in group))
        best.all_accessions = "|".join(all_ids)
        best.duplicate_count = len(all_ids)

        deduped.append(best)

    removed = len(records) - len(deduped)
    if removed > 0:
        logger.info(
            f"Deduplication: removed {removed} duplicates, {len(deduped)} unique records"
        )

    return deduped


def merge_with_existing(
    new_records: list[SPRecord],
    existing_csv: Optional[str] = None,
) -> pd.DataFrame:
    """Merge new records with an existing CSV file (idempotent append).

    Uses sp_sequence + sp_type as the unique key — same logic as dedup.

    Args:
        new_records: New SPRecord instances (already deduplicated).
        existing_csv: Path to existing CSV file to merge with.

    Returns:
        Merged DataFrame ready for export.
    """
    new_df = pd.DataFrame([r.to_dict() for r in new_records])

    if existing_csv and Path(existing_csv).exists():
        logger.info(f"Merging with existing file: {existing_csv}")
        existing_df = pd.read_csv(existing_csv)

        # Key: sp_sequence + sp_type
        def make_key(row):
            return f"{row.get('sp_sequence', '')}|{row.get('sp_type', '')}"

        existing_df["_key"] = existing_df.apply(make_key, axis=1)
        new_df["_key"] = new_df.apply(make_key, axis=1)

        existing_keys = set(existing_df["_key"])
        new_only = new_df[~new_df["_key"].isin(existing_keys)]

        logger.info(f"Existing: {len(existing_df)}, New unique: {len(new_only)}")

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
