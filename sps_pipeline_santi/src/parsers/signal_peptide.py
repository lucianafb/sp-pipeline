"""
Parser de Péptidos Señal Type 1 (Clivables).

Extrae signal peptides clásicos N-terminales desde features de UniProt.
"""

import logging
from typing import Optional

from ..processing.evidence import classify_evidence
from ..processing.features import compute_features

logger = logging.getLogger(__name__)


def extract_cleavage_site(full_seq: str, end_pos: int, context: int = 5) -> str:
    """
    Extrae motivo del sitio de corte alrededor de end_pos.

    UniProt usa indexación 1-based. end_pos = último AA del SP.
    El corte ocurre DESPUÉS de end_pos.

    Formato: P5-P4-P3-P2-P1|P1'-P2'-P3'-P4'-P5'
    """
    # Convertir a 0-based: end_pos en 1-based → cut_index en 0-based
    cut_index = end_pos  # En 0-based, el corte está en posición end_pos

    if cut_index - context < 0 or cut_index + context > len(full_seq):
        # Ajustar contexto disponible
        pre_start = max(0, cut_index - context)
        post_end = min(len(full_seq), cut_index + context)
        pre_cut = full_seq[pre_start:cut_index]
        post_cut = full_seq[cut_index:post_end]
    else:
        pre_cut = full_seq[cut_index - context : cut_index]
        post_cut = full_seq[cut_index : cut_index + context]

    return f"{pre_cut}|{post_cut}"


def parse_signal_peptides(entry: dict, query_group: str = "") -> list[dict]:
    """
    Extrae péptidos señal Type 1 de una entrada de UniProt.

    Args:
        entry: Entrada JSON de UniProt.
        query_group: Nombre del preset/grupo de query.

    Returns:
        Lista de dicts con datos del SP.
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

    # Determinar si es Swiss-Prot o TrEMBL
    entry_type = entry.get("entryType", "")
    is_reviewed = "Swiss-Prot" in entry_type
    source = "Swiss-Prot" if is_reviewed else "TrEMBL"

    features = entry.get("features", [])

    for ft in features:
        ft_type = ft.get("type", "")

        if ft_type != "Signal":
            continue

        location = ft.get("location", {})
        start_info = location.get("start", {})
        end_info = location.get("end", {})

        start = start_info.get("value")
        end = end_info.get("value")

        # Validar posiciones numéricas
        if not isinstance(start, int) or not isinstance(end, int):
            logger.debug(f"Skipping {accession}: posiciones no numéricas ({start}, {end})")
            continue

        # Extraer secuencia del SP (UniProt 1-based)
        sp_seq = full_sequence[start - 1 : end]

        if not sp_seq:
            continue

        # Evidencia
        evidence_category, evidence_codes = classify_evidence(ft)

        # Sitio de corte
        cleavage_motif = extract_cleavage_site(full_sequence, end)

        # Features computados
        sp_features = compute_features(sp_seq)

        results.append({
            "accession": accession,
            "gene": gene_name,
            "organism": organism,
            "taxonomy_id": taxonomy_id,
            "query_group": query_group,
            "sp_type": "SIGNAL_PEPTIDE_TYPE1",
            "sp_subtype": "",
            "evidence_category": evidence_category,
            "evidence_codes": "|".join(evidence_codes),
            "sp_sequence_aa": sp_seq,
            "sp_length": len(sp_seq),
            "start_pos": start,
            "end_pos": end,
            "cleavage_site_motif": cleavage_motif,
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
        })

    return results
