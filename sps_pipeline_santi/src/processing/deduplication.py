"""
Deduplicación de péptidos señal por secuencia exacta.

Agrupa registros con la misma secuencia de SP y conserva el de mejor evidencia,
fusionando metadata (accessions, organisms, etc.).
"""

import logging
import pandas as pd

from ..processing.evidence import get_evidence_priority

logger = logging.getLogger(__name__)


def deduplicate_by_sequence(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplica un DataFrame de SPs por secuencia exacta (sp_sequence_aa).

    Para cada grupo de secuencias idénticas:
    - Conserva el registro con mejor categoría de evidencia
    - Fusiona los accessions en 'all_accessions'
    - Cuenta duplicados en 'duplicate_count'

    Args:
        df: DataFrame con columnas estándar del pipeline.

    Returns:
        DataFrame deduplicado.
    """
    if df.empty:
        return df

    original_count = len(df)
    logger.info(f"Deduplicando {original_count} registros por secuencia exacta...")

    # Agregar columna de prioridad para ordenar
    df = df.copy()
    df["_evidence_priority"] = df["evidence_category"].apply(get_evidence_priority)

    # Agrupar por secuencia exacta + sp_type (no mezclar Type1 con Type2)
    grouped = df.groupby(["sp_sequence_aa", "sp_type"])

    deduped_records = []

    for (seq, sp_type), group in grouped:
        # Ordenar por prioridad de evidencia (menor = mejor)
        group_sorted = group.sort_values("_evidence_priority")

        # Tomar el mejor registro como representante
        best = group_sorted.iloc[0].to_dict()

        # Fusionar accessions
        all_accessions = sorted(group["accession"].unique())
        best["all_accessions"] = "|".join(all_accessions)
        best["duplicate_count"] = len(all_accessions)

        deduped_records.append(best)

    result = pd.DataFrame(deduped_records)

    # Limpiar columna temporal
    if "_evidence_priority" in result.columns:
        result = result.drop(columns=["_evidence_priority"])

    deduped_count = len(result)
    removed = original_count - deduped_count
    logger.info(f"Deduplicación completa: {original_count} → {deduped_count} ({removed} duplicados eliminados)")

    return result.reset_index(drop=True)
