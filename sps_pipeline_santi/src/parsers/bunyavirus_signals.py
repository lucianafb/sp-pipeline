"""
Parser de señales para Orthobunyavirus y Orthohantavirus.

Ambas familias comparten la arquitectura del precursor GPC (Glycoprotein Precursor /
Envelopment Polyprotein):

    NH₂─[ SP ]─[ Gn ]─────────────[ Gc ]─COOH

donde:
  - SP: Péptido señal clásico N-terminal (Type 1, clivable, Sec61-dependiente)
  - Gn: Glycoproteína N. Anclada a membrana por TM.
  - Gc: Glycoproteína C. Proteína de fusión, también con TM.

La señal principal que buscamos es el SP N-terminal del GPC (~15-20 aa).
Opcionalmente, el boundary TM entre Gn y Gc puede actuar como señal interna.

Estrategias de extracción:
  1. Feature "Signal" explícito (disponible en entries reviewed de LACV, ANDV, SNV, BCCV)
  2. Inferencia desde Chain boundaries (si hay Chain "Glycoprotein N", el SP es pos 1..inicio_Gn-1)
  3. Fallback TM: primera TM region cerca del N-terminal
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


def _get_protein_name(entry: dict) -> str:
    """Extrae el nombre de la proteína."""
    pd = entry.get("proteinDescription", {})
    rn = pd.get("recommendedName", {})
    if rn:
        return rn.get("fullName", {}).get("value", "")
    sn = pd.get("submissionNames", [])
    if sn:
        return sn[0].get("fullName", {}).get("value", "")
    return ""


def _is_gpc_entry(entry: dict) -> bool:
    """Determina si la entrada es un precursor GPC / envelopment polyprotein."""
    name = _get_protein_name(entry).lower()
    keywords = [
        "envelopment polyprotein",
        "glycoprotein precursor",
        "polyprotein m",
        "m polyprotein",
        "glycoprotein",
        "envelope glycoprotein",
    ]
    return any(kw in name for kw in keywords)


def _get_features_by_type(features: list, ft_type: str) -> list:
    """Extrae features de un tipo específico."""
    result = []
    for ft in features:
        if ft.get("type") == ft_type:
            loc = ft.get("location", {})
            s = loc.get("start", {}).get("value")
            e = loc.get("end", {}).get("value")
            if isinstance(s, int) and isinstance(e, int):
                result.append({
                    "start": s,
                    "end": e,
                    "description": ft.get("description", ""),
                    "feature": ft,
                })
    return sorted(result, key=lambda x: x["start"])


def _find_gn_gc_chains(chains: list) -> tuple:
    """
    Busca chains correspondientes a Glycoprotein N (Gn) y Glycoprotein C (Gc).
    Retorna (gn_chain, gc_chain) o (None, None).
    """
    gn = None
    gc = None
    for ch in chains:
        desc = ch["description"].lower()
        # Gn patterns
        if any(p in desc for p in ["glycoprotein n", "glycoprotein gn", "protein gn"]):
            gn = ch
        # Gc patterns
        elif any(p in desc for p in ["glycoprotein c", "glycoprotein gc", "protein gc"]):
            gc = ch
    return gn, gc


# =============================================================================
# PARSER PRINCIPAL
# =============================================================================

def parse_bunyavirus_signals(entry: dict, query_group: str = "") -> list[dict]:
    """
    Extrae señales del precursor GPC de Orthobunyavirus y Orthohantavirus.

    Solo procesa entries que sean GPC/envelopment polyprotein.

    Estrategia:
    1. Feature "Signal" explícito → SP Type 1 clásico
    2. Chains Gn/Gc → inferir SP como pos 1..(inicio_Gn - 1)
    3. Fallback TM → primera TM near N-terminal como señal putativa
    """
    results = []
    meta = _get_entry_metadata(entry)
    features = entry.get("features", [])
    full_seq = meta["full_sequence"]

    if not full_seq:
        return results

    # Solo procesar entries que parezcan ser GPC/glycoproteínas
    prot_name = _get_protein_name(entry)
    if not _is_gpc_entry(entry):
        # No es una glycoproteína / GPC — skip
        # (esto filtra nucleoproteínas, RNA polimerasas, etc.)
        logger.debug(f"  Skip {meta['accession']}: '{prot_name}' no es GPC")
        return results

    signals = _get_features_by_type(features, "Signal")
    chains = _get_features_by_type(features, "Chain")
    tms = _get_features_by_type(features, "Transmembrane")

    # -------------------------------------------------------------------------
    # Estrategia 1: Feature "Signal" explícito
    # -------------------------------------------------------------------------
    for sig in signals:
        sp_seq = full_seq[sig["start"] - 1 : sig["end"]]
        if not sp_seq or len(sp_seq) < 5:
            continue

        ev_cat, ev_codes = classify_evidence(sig["feature"])
        results.append(_build_record(
            meta["accession"], meta["gene"], meta["organism"],
            meta["taxonomy_id"], query_group,
            sp_seq, sig["start"], sig["end"],
            "SIGNAL_PEPTIDE_TYPE1",
            f"GPC signal peptide ({prot_name[:40]})",
            ev_cat, ev_codes,
            _extract_cleavage(full_seq, sig["end"]),
            full_seq, meta["source"], meta["is_reviewed"],
        ))

    # -------------------------------------------------------------------------
    # Estrategia 2: Inferir SP desde Chain boundaries
    # Si hay Chain para Gn que empiece en pos > 1, el SP es 1..(inicio_Gn - 1)
    # -------------------------------------------------------------------------
    if not results and chains:
        gn, gc = _find_gn_gc_chains(chains)

        if gn and gn["start"] > 1:
            # El SP es la secuencia antes de Gn
            sp_end = gn["start"] - 1
            sp_seq = full_seq[0:sp_end]

            if sp_seq and 5 <= len(sp_seq) <= 50:
                ev_cat, ev_codes = classify_evidence(gn["feature"])
                results.append(_build_record(
                    meta["accession"], meta["gene"], meta["organism"],
                    meta["taxonomy_id"], query_group,
                    sp_seq, 1, sp_end,
                    "SIGNAL_PEPTIDE_TYPE1",
                    f"GPC signal peptide (inferred from Gn start)",
                    ev_cat, ev_codes,
                    _extract_cleavage(full_seq, sp_end),
                    full_seq, meta["source"], meta["is_reviewed"],
                ))

        # Señal interna Gn→Gc (TM boundary)
        if gn and gc:
            # Buscar TM region entre el final de Gn y el inicio de Gc
            gn_end = gn["end"]
            gc_start = gc["start"]
            for tm in tms:
                # TM que está en la frontera Gn/Gc (cerca del final de Gn)
                if abs(tm["end"] - gn_end) < 30 or abs(tm["start"] - gc_start) < 30:
                    tm_seq = full_seq[tm["start"] - 1 : tm["end"]]
                    if tm_seq and len(tm_seq) >= 10:
                        # Verificar que no es duplicado
                        is_dup = any(
                            r["start_pos"] == tm["start"] and r["end_pos"] == tm["end"]
                            for r in results
                        )
                        if not is_dup:
                            ev_cat, ev_codes = classify_evidence(tm["feature"])
                            results.append(_build_record(
                                meta["accession"], meta["gene"], meta["organism"],
                                meta["taxonomy_id"], query_group,
                                tm_seq, tm["start"], tm["end"],
                                "VIRAL_INTERNAL_SIGNAL",
                                f"Gn→Gc boundary TM signal",
                                ev_cat, ev_codes,
                                _extract_cleavage(full_seq, tm["end"]),
                                full_seq, meta["source"], meta["is_reviewed"],
                            ))
                            break  # Solo uno

    # -------------------------------------------------------------------------
    # Estrategia 3: Fallback — TM near N-terminal
    # Para TrEMBL entries sin Signal ni Chain annotations
    # -------------------------------------------------------------------------
    if not results and tms:
        # Buscar primera TM acerca del N-terminal
        first_tm = tms[0]
        if first_tm["start"] < 80:
            # Tomar secuencia antes de la primera TM como SP putativo
            sp_end = first_tm["start"] - 1
            if sp_end > 5:
                sp_seq = full_seq[0:sp_end]
                if sp_seq and 5 <= len(sp_seq) <= 50:
                    ev_cat, ev_codes = classify_evidence(first_tm["feature"])
                    results.append(_build_record(
                        meta["accession"], meta["gene"], meta["organism"],
                        meta["taxonomy_id"], query_group,
                        sp_seq, 1, sp_end,
                        "SIGNAL_PEPTIDE_TYPE1",
                        f"GPC signal peptide (putative, from TM position)",
                        ev_cat, ev_codes,
                        _extract_cleavage(full_seq, sp_end),
                        full_seq, meta["source"], meta["is_reviewed"],
                    ))

            # También la TM misma podría ser significativa
            tm_seq = full_seq[first_tm["start"] - 1 : first_tm["end"]]
            if tm_seq and len(tm_seq) >= 10:
                is_dup = any(
                    r["start_pos"] == first_tm["start"] and r["end_pos"] == first_tm["end"]
                    for r in results
                )
                if not is_dup:
                    ev_cat, ev_codes = classify_evidence(first_tm["feature"])
                    results.append(_build_record(
                        meta["accession"], meta["gene"], meta["organism"],
                        meta["taxonomy_id"], query_group,
                        tm_seq, first_tm["start"], first_tm["end"],
                        "VIRAL_INTERNAL_SIGNAL",
                        f"N-proximal TM (putative GPC signal)",
                        ev_cat, ev_codes,
                        "N/A (TM-based)",
                        full_seq, meta["source"], meta["is_reviewed"],
                    ))

    if results:
        logger.info(
            f"  {meta['accession']}: {prot_name[:40]} → "
            f"{len(results)} señal(es) extraída(s)"
        )

    return results
