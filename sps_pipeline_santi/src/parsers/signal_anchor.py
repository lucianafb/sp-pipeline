"""
Parser de Signal Anchors Type 2 (No clivables).

Extrae señales de anclaje N-terminales desde features de UniProt.
Estrategia dual: busca features "Signal" no clivables y "Transmembrane" cerca del N-terminal.
"""

import logging
from typing import Optional

from ..processing.evidence import classify_evidence
from ..processing.features import compute_features

logger = logging.getLogger(__name__)

# Umbral: TM region que empieza antes de esta posición es candidata a signal anchor
SIGNAL_ANCHOR_MAX_START = 60


def parse_signal_anchors(entry: dict, query_group: str = "") -> list[dict]:
    """
    Extrae signal anchors (Type 2) de una entrada de UniProt.

    Busca:
    1. Features "Signal" (algunas entradas anotan signal anchors como Signal)
    2. Features "Transmembrane" que empiecen cerca del N-terminal (< 60 aa)

    Args:
        entry: Entrada JSON de UniProt.
        query_group: Nombre del preset/grupo de query.

    Returns:
        Lista de dicts con datos del signal anchor.
    """
    results = []

    accession = entry.get("primaryAccession", "N/A")
    genes_list = entry.get("genes", [])
    gene_name = "N/A"
    if genes_list:
        gene_name = genes_list[0].get("geneName", {}).get("value", "N/A")

    organism_info = entry.get("organism", {})
    organism = organism_info.get("scientificName", "N/A")
    taxonomy_id = organism_info.get("taxonId", "N/A")

    sequence_info = entry.get("sequence", {})
    full_sequence = sequence_info.get("value", "")

    entry_type = entry.get("entryType", "")
    is_reviewed = "Swiss-Prot" in entry_type
    source = "Swiss-Prot" if is_reviewed else "TrEMBL"

    features = entry.get("features", [])

    # Track si ya encontramos un signal anchor para esta entrada
    found_anchor = False

    for ft in features:
        ft_type = ft.get("type", "")

        location = ft.get("location", {})
        start_info = location.get("start", {})
        end_info = location.get("end", {})

        start = start_info.get("value")
        end = end_info.get("value")

        if not isinstance(start, int) or not isinstance(end, int):
            continue

        # Estrategia 1: Feature "Signal" en entradas con keyword signal-anchor
        # (UniProt ocasionalmente los anota así)
        if ft_type == "Signal" and not found_anchor:
            sa_seq = full_sequence[start - 1 : end]
            if not sa_seq:
                continue

            evidence_category, evidence_codes = classify_evidence(ft)
            sp_features = compute_features(sa_seq)

            results.append(_build_record(
                accession, gene_name, organism, taxonomy_id, query_group,
                sa_seq, start, end, evidence_category, evidence_codes,
                sp_features, full_sequence, source, is_reviewed,
            ))
            found_anchor = True

        # Estrategia 2: Transmembrane cerca del N-terminal
        elif ft_type == "Transmembrane" and start <= SIGNAL_ANCHOR_MAX_START and not found_anchor:
            sa_seq = full_sequence[start - 1 : end]
            if not sa_seq:
                continue

            evidence_category, evidence_codes = classify_evidence(ft)
            sp_features = compute_features(sa_seq)

            results.append(_build_record(
                accession, gene_name, organism, taxonomy_id, query_group,
                sa_seq, start, end, evidence_category, evidence_codes,
                sp_features, full_sequence, source, is_reviewed,
            ))
            found_anchor = True

    return results


def _build_record(
    accession, gene, organism, taxonomy_id, query_group,
    sa_seq, start, end, evidence_category, evidence_codes,
    sp_features, full_sequence, source, is_reviewed,
) -> dict:
    """Construye un registro estándar de signal anchor."""
    return {
        "accession": accession,
        "gene": gene,
        "organism": organism,
        "taxonomy_id": taxonomy_id,
        "query_group": query_group,
        "sp_type": "SIGNAL_ANCHOR_TYPE2",
        "sp_subtype": "",
        "evidence_category": evidence_category,
        "evidence_codes": "|".join(evidence_codes),
        "sp_sequence_aa": sa_seq,
        "sp_length": len(sa_seq),
        "start_pos": start,
        "end_pos": end,
        "cleavage_site_motif": "N/A (Non-cleavable)",
        "hydrophobicity_mean": sp_features["hydrophobicity_mean"],
        "net_charge_ph7": sp_features["net_charge_ph7"],
        "n_region": sp_features["n_region"],
        "h_region": sp_features["h_region"],
        "c_region": sp_features["c_region"],
        "full_sequence": full_sequence,
        "duplicate_count": 1,
        "all_accessions": accession,
        "source": source,
        "reviewed": is_reviewed,
    }
