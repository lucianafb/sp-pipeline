"""Export modules for SP-Pipeline.

Handles writing results to CSV, TSV, and optional FASTA formats.
"""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from sp_pipeline.models import SPRecord, CSV_COLUMNS

logger = logging.getLogger("sp_pipeline.exporters")


def export_csv(
    df: pd.DataFrame,
    output_path: str,
    delimiter: str = ",",
) -> str:
    """Export DataFrame to CSV/TSV file.

    Args:
        df: DataFrame with SP records.
        output_path: Path to output file.
        delimiter: Delimiter character.

    Returns:
        Path to the written file.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(path, index=False, sep=delimiter)
    logger.info(f"Exported {len(df)} records to {path}")

    return str(path)


def export_fasta(
    records: list[SPRecord],
    output_path: str,
    sequence_type: str = "sp",
) -> str:
    """Export sequences to FASTA format.

    Args:
        records: List of SPRecord instances.
        output_path: Path to output FASTA file.
        sequence_type: "sp" for signal peptide only, "full" for full protein.

    Returns:
        Path to the written file.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        for record in records:
            if sequence_type == "sp":
                seq = record.sp_sequence
                header = (
                    f">{record.entry_id}|{record.source_db}|SP|"
                    f"{record.sp_start}-{record.sp_end}|"
                    f"{record.organism}|{record.sp_type}"
                )
            else:
                seq = record.full_sequence
                header = (
                    f">{record.entry_id}|{record.source_db}|FULL|"
                    f"{record.organism}|{record.protein_name}"
                )

            if seq:
                f.write(f"{header}\n")
                # Wrap at 80 characters
                for i in range(0, len(seq), 80):
                    f.write(f"{seq[i:i+80]}\n")

    logger.info(f"Exported {len(records)} sequences to {path}")
    return str(path)


def generate_summary(df: pd.DataFrame) -> str:
    """Generate a text summary of the query results.

    Args:
        df: DataFrame with SP records.

    Returns:
        Summary string.
    """
    lines = [
        "=" * 60,
        "SP-Pipeline Query Summary",
        "=" * 60,
        f"Total records: {len(df)}",
        "",
    ]

    # By query group
    if "query_group" in df.columns:
        lines.append("Records by query group:")
        for group, count in df["query_group"].value_counts().items():
            lines.append(f"  {group}: {count}")
        lines.append("")

    # By source
    if "source_db" in df.columns:
        lines.append("Records by source:")
        for source, count in df["source_db"].value_counts().items():
            lines.append(f"  {source}: {count}")
        lines.append("")

    # By SP type
    if "sp_type" in df.columns:
        lines.append("Records by SP type:")
        for sp_type, count in df["sp_type"].value_counts().items():
            lines.append(f"  {sp_type}: {count}")
        lines.append("")

    # By evidence
    if "evidence_method" in df.columns:
        lines.append("Records by evidence method:")
        for method, count in df["evidence_method"].value_counts().items():
            lines.append(f"  {method}: {count}")
        lines.append("")

    # By organism (top 10)
    if "organism" in df.columns:
        lines.append("Top 10 organisms:")
        for org, count in df["organism"].value_counts().head(10).items():
            lines.append(f"  {org}: {count}")
        lines.append("")

    # SP length stats
    if "sp_length" in df.columns and df["sp_length"].notna().any():
        lengths = df["sp_length"].dropna()
        lines.extend([
            "Signal peptide length statistics:",
            f"  Mean: {lengths.mean():.1f}",
            f"  Median: {lengths.median():.1f}",
            f"  Min: {lengths.min():.0f}",
            f"  Max: {lengths.max():.0f}",
            f"  Std: {lengths.std():.1f}",
        ])

    lines.append("=" * 60)
    return "\n".join(lines)
