"""
Parser de señales virales para Flavivirus y Alphavirus.

Maneja las particularidades de las señales internas de poliproteínas virales.

IMPORTANTE: Las entradas virales en UniProt NO tienen features tipo "Signal".
Las señales se identifican a partir de:
  - Features "Chain" (límites de las proteínas procesadas)
  - Features "Transmembrane" (dominios TM que actúan como señales internas)
  - Las secuencias señal se infieren de las regiones entre el fin de un chain
    y el inicio del siguiente.

Flavivirus:
  Poliproteína: C → prM → E → NS1 → ...
  - anchC (C-term de Capsid) → señal para prM
  - prM (C-term TM) → señal para E
  - E (C-term TM2) → señal para NS1

Alphavirus:
  Polipéptido estructural: C → E3 → E2 → 6K → E1
  - E3 = señal para E2 (Sec61-dependiente)
  - 6K = señal para E1 (Sec61-dependiente)
"""

import logging
from typing import Optional

from ..processing.evidence import classify_evidence
from ..processing.features import compute_features

logger = logging.getLogger(__name__)


def _build_record(
    accession, gene, organism, taxonomy_id, query_group,
    sp_seq, start, end, sp_type, sp_subtype,
    evidence_category, evidence_codes, cleavage_motif,
    full_sequence, source, is_reviewed,
) -> dict:
    """Construye un registro estándar."""
    sp_features = compute_features(sp_seq)
    return {
        "accession": accession,
        "gene": gene,
        "organism": organism,
        "taxonomy_id": taxonomy_id,
        "query_group": query_group,
        "sp_type": sp_type,
        "sp_subtype": sp_subtype,
        "evidence_category": evidence_category,
        "evidence_codes": "|".join(evidence_codes) if evidence_codes else "",
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
    }


def _extract_cleavage(full_seq: str, end_pos: int, context: int = 5) -> str:
    """Extrae motivo de sitio de corte."""
    cut_index = end_pos
    pre_start = max(0, cut_index - context)
    post_end = min(len(full_seq), cut_index + context)
    pre_cut = full_seq[pre_start:cut_index]
    post_cut = full_seq[cut_index:post_end]
    if not pre_cut or not post_cut:
        return "N/A"
    return f"{pre_cut}|{post_cut}"


def _get_entry_metadata(entry: dict) -> dict:
    """Extrae metadata común de una entrada de UniProt."""
    accession = entry.get("primaryAccession", "N/A")
    genes_list = entry.get("genes", [])
    gene_name = "N/A"
    if genes_list:
        gene_name = genes_list[0].get("geneName", {}).get("value", "N/A")

    organism_info = entry.get("organism", {})
    organism = organism_info.get("scientificName", "N/A")
    taxonomy_id = organism_info.get("taxonId", "N/A")

    full_sequence = entry.get("sequence", {}).get("value", "")

    entry_type = entry.get("entryType", "")
    is_reviewed = "Swiss-Prot" in entry_type
    source = "Swiss-Prot" if is_reviewed else "TrEMBL"

    return {
        "accession": accession,
        "gene": gene_name,
        "organism": organism,
        "taxonomy_id": taxonomy_id,
        "full_sequence": full_sequence,
        "source": source,
        "is_reviewed": is_reviewed,
    }


def _get_chains(features: list) -> dict:
    """Extrae y organiza Chain features."""
    chains = {}
    for ft in features:
        if ft.get("type") == "Chain":
            desc = ft.get("description", "")
            loc = ft.get("location", {})
            s = loc.get("start", {}).get("value")
            e = loc.get("end", {}).get("value")
            if isinstance(s, int) and isinstance(e, int):
                chains[desc] = {"start": s, "end": e, "feature": ft}
    return chains


def _get_transmembrane_regions(features: list) -> list:
    """Extrae features Transmembrane."""
    tms = []
    for ft in features:
        if ft.get("type") == "Transmembrane":
            loc = ft.get("location", {})
            s = loc.get("start", {}).get("value")
            e = loc.get("end", {}).get("value")
            if isinstance(s, int) and isinstance(e, int):
                tms.append({"start": s, "end": e, "feature": ft})
    return sorted(tms, key=lambda x: x["start"])


# =============================================================================
# FLAVIVIRUS
# =============================================================================

# Nombres de chains de interés en flavivirus y sus señales
FLAVI_SIGNAL_MAP = {
    # chain_name_pattern: (signal_role, target_protein)
    "Capsid": ("anchC", "prM"),
    "Protein prM": ("prM C-term", "E"),
    "Envelope": ("E C-term", "NS1"),
    # Versiones alternativas
    "Capsid protein C": ("anchC", "prM"),
    "Genome polyprotein": (None, None),  # Es la poliproteína completa
}


def _match_chain_name(chain_name: str, patterns: dict) -> Optional[tuple]:
    """Busca coincidencia parcial de nombre de chain."""
    chain_lower = chain_name.lower()
    for pattern, value in patterns.items():
        if pattern.lower() in chain_lower:
            return value
    return None


def parse_flavivirus_signals(entry: dict, query_group: str = "") -> list[dict]:
    """
    Extrae señales de poliproteína de Flavivirus.

    Estrategia:
    1. Buscar features "Signal" (raro pero posible)
    2. Usar boundaries de Chain features para identificar señales internas
    3. Extraer TM regions en los C-terminales de proteínas estructurales
    """
    results = []
    meta = _get_entry_metadata(entry)
    features = entry.get("features", [])
    full_seq = meta["full_sequence"]

    if not full_seq:
        return results

    chains = _get_chains(features)
    tms = _get_transmembrane_regions(features)

    # Estrategia 1: Buscar features "Signal" explícitos (raro)
    for ft in features:
        if ft.get("type") == "Signal":
            loc = ft.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            if isinstance(start, int) and isinstance(end, int):
                sp_seq = full_seq[start - 1 : end]
                if sp_seq:
                    ev_cat, ev_codes = classify_evidence(ft)
                    results.append(_build_record(
                        meta["accession"], meta["gene"], meta["organism"],
                        meta["taxonomy_id"], query_group,
                        sp_seq, start, end,
                        "SIGNAL_PEPTIDE_TYPE1" if start <= 5 else "VIRAL_INTERNAL_SIGNAL",
                        ft.get("description", "Signal"),
                        ev_cat, ev_codes, _extract_cleavage(full_seq, end),
                        full_seq, meta["source"], meta["is_reviewed"],
                    ))

    # Estrategia 2: Inferir señales desde TM regions en C-terminales de chains
    # Las TM regions que están al final de una proteína estructural funcionan
    # como señales para la siguiente proteína
    for chain_name, chain_info in chains.items():
        # Ignorar la poliproteína completa
        if "polyprotein" in chain_name.lower() or "genome" in chain_name.lower():
            continue

        chain_start = chain_info["start"]
        chain_end = chain_info["end"]

        # Buscar TM regions dentro de este chain (especialmente cerca del C-term)
        for tm in tms:
            # TM dentro del chain y cerca del C-terminal (últimos 60 aa del chain)
            if tm["start"] >= chain_start and tm["end"] <= chain_end + 5:
                if abs(tm["end"] - chain_end) < 30:  # Cerca del C-term
                    tm_seq = full_seq[tm["start"] - 1 : tm["end"]]
                    if not tm_seq or len(tm_seq) < 10:
                        continue

                    # Determinar qué proteína transloca esta TM
                    sp_subtype = f"{chain_name} C-term TM → next protein signal"

                    # Buscar la match en nuestro mapa
                    for pattern, (role, target) in FLAVI_SIGNAL_MAP.items():
                        if role and pattern.lower() in chain_name.lower():
                            sp_subtype = f"{role}→{target} signal"
                            break

                    ev_cat, ev_codes = classify_evidence(tm["feature"])

                    # Verificar que no sea duplicado
                    is_dup = any(r["start_pos"] == tm["start"] and r["end_pos"] == tm["end"] for r in results)
                    if not is_dup:
                        results.append(_build_record(
                            meta["accession"], meta["gene"], meta["organism"],
                            meta["taxonomy_id"], query_group,
                            tm_seq, tm["start"], tm["end"],
                            "VIRAL_INTERNAL_SIGNAL", sp_subtype,
                            ev_cat, ev_codes,
                            _extract_cleavage(full_seq, tm["end"]),
                            full_seq, meta["source"], meta["is_reviewed"],
                        ))

    # Estrategia 3: Si no encontramos nada via chains pero hay TMs,
    # tomar la primera TM region cerca del N-terminal como señal
    if not results and tms:
        first_tm = tms[0]
        if first_tm["start"] < 150:  # Razonablemente cerca del N-terminal
            tm_seq = full_seq[first_tm["start"] - 1 : first_tm["end"]]
            if tm_seq and len(tm_seq) >= 10:
                ev_cat, ev_codes = classify_evidence(first_tm["feature"])
                results.append(_build_record(
                    meta["accession"], meta["gene"], meta["organism"],
                    meta["taxonomy_id"], query_group,
                    tm_seq, first_tm["start"], first_tm["end"],
                    "VIRAL_INTERNAL_SIGNAL", "N-proximal TM (putative signal)",
                    ev_cat, ev_codes,
                    _extract_cleavage(full_seq, first_tm["end"]),
                    full_seq, meta["source"], meta["is_reviewed"],
                ))

    return results


# =============================================================================
# ALPHAVIRUS
# =============================================================================

# Chains to EXCLUDE (precursors, mature proteins that are not signals)
ALPHA_CHAINS_EXCLUDE = [
    "precursor",
    "spike",
    "capsid",
    "envelope glycoprotein e1",
    "glycoprotein e2",
    "non-structural",
    "polyprotein",
    "nsp",
]

# Chains de interés positivo en alphavirus — E3 and 6K are signal peptides
ALPHA_CHAINS_POSITIVE = {
    "assembly protein e3": "E3→p62/E2 signal (Sec61-dependent)",
    "e3": "E3→p62/E2 signal (Sec61-dependent)",
    "6k protein": "6K→E1 signal (Sec61-dependent)",
    "6kda protein": "6K→E1 signal (Sec61-dependent)",
    "6k": "6K→E1 signal (Sec61-dependent)",
}


def parse_alphavirus_signals(entry: dict, query_group: str = "") -> list[dict]:
    """
    Extrae señales E3 y 6K de Alphavirus.

    Estrategia:
    1. Buscar features "Signal" explícitos (raro)
    2. Extraer chains E3 y 6K directamente — estas proteínas SON las señales
    3. Si no hay chains anotados, buscar TM regions
    """
    results = []
    meta = _get_entry_metadata(entry)
    features = entry.get("features", [])
    full_seq = meta["full_sequence"]

    if not full_seq:
        return results

    chains = _get_chains(features)
    tms = _get_transmembrane_regions(features)

    # Estrategia 1: Buscar features "Signal" explícitos
    for ft in features:
        if ft.get("type") == "Signal":
            loc = ft.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            if isinstance(start, int) and isinstance(end, int):
                sp_seq = full_seq[start - 1 : end]
                if sp_seq:
                    ev_cat, ev_codes = classify_evidence(ft)
                    results.append(_build_record(
                        meta["accession"], meta["gene"], meta["organism"],
                        meta["taxonomy_id"], query_group,
                        sp_seq, start, end,
                        "VIRAL_INTERNAL_SIGNAL", ft.get("description", "Signal"),
                        ev_cat, ev_codes, _extract_cleavage(full_seq, end),
                        full_seq, meta["source"], meta["is_reviewed"],
                    ))

    # Estrategia 2: Extraer chains E3 y 6K
    # En alphavirus, E3 y 6K SON las señales (no solo contienen señales)
    for chain_name, chain_info in chains.items():
        chain_lower = chain_name.lower()

        # First: check if this chain should be excluded
        if any(excl in chain_lower for excl in ALPHA_CHAINS_EXCLUDE):
            continue

        # Then: check if it matches a positive pattern
        role = None
        for pattern, signal_role in ALPHA_CHAINS_POSITIVE.items():
            if pattern in chain_lower:
                role = signal_role
                break

        if role is None:
            continue

        chain_start = chain_info["start"]
        chain_end = chain_info["end"]
        chain_seq = full_seq[chain_start - 1 : chain_end]

        if not chain_seq:
            continue

        # Verificar que no sea duplicado de un Signal feature ya encontrado
        is_dup = any(
            r["start_pos"] == chain_start and r["end_pos"] == chain_end
            for r in results
        )
        if is_dup:
            continue

        ev_cat, ev_codes = classify_evidence(chain_info["feature"])

        results.append(_build_record(
            meta["accession"], meta["gene"], meta["organism"],
            meta["taxonomy_id"], query_group,
            chain_seq, chain_start, chain_end,
            "VIRAL_INTERNAL_SIGNAL", role,
            ev_cat, ev_codes,
            _extract_cleavage(full_seq, chain_end),
            full_seq, meta["source"], meta["is_reviewed"],
        ))

    # Estrategia 3: Fallback — TM regions si no encontramos chains
    if not results and tms:
        for tm in tms:
            if tm["start"] < 60:
                tm_seq = full_seq[tm["start"] - 1 : tm["end"]]
                if tm_seq and len(tm_seq) >= 10:
                    ev_cat, ev_codes = classify_evidence(tm["feature"])
                    results.append(_build_record(
                        meta["accession"], meta["gene"], meta["organism"],
                        meta["taxonomy_id"], query_group,
                        tm_seq, tm["start"], tm["end"],
                        "VIRAL_INTERNAL_SIGNAL", "TM-based signal (putative)",
                        ev_cat, ev_codes,
                        "N/A (TM-based)",
                        full_seq, meta["source"], meta["is_reviewed"],
                    ))
                    break  # Solo tomar la primera TM near N-term

    return results
