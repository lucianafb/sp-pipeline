"""
Clasificación de evidencia basada en códigos ECO de UniProt.

Categoriza features como EXPERIMENTAL, BY_SIMILARITY o PREDICTED
según los códigos de evidencia asociados.
"""

from typing import Optional


# Mapeo de códigos ECO a categorías
# Ref: https://www.uniprot.org/help/evidences
ECO_CATEGORIES = {
    # Evidencia experimental directa
    "ECO:0000269": "EXPERIMENTAL",    # Experimental - manual assertion
    "ECO:0000303": "EXPERIMENTAL",    # Non-traceable author statement (published)
    "ECO:0000305": "EXPERIMENTAL",    # Curator inference

    # Inferencia por similitud
    "ECO:0000250": "BY_SIMILARITY",   # Sequence similarity
    "ECO:0000255": "BY_SIMILARITY",   # Sequence model (PROSITE, Pfam rules)

    # Predicción automática
    "ECO:0000256": "PREDICTED",       # Automatic annotation (SAM, etc.)
    "ECO:0007829": "PREDICTED",       # Automatic assertion (evidence used in manual)
    "ECO:0000259": "PREDICTED",       # Automatic annotation by PROSITE
    "ECO:0000312": "PREDICTED",       # Imported from another database
}

# Orden de prioridad para deduplicación (menor = mejor)
EVIDENCE_PRIORITY = {
    "EXPERIMENTAL": 0,
    "BY_SIMILARITY": 1,
    "PREDICTED": 2,
    "PREDICTED (No Evidence)": 3,
}


def classify_evidence(feature: dict) -> tuple[str, list[str]]:
    """
    Clasifica la evidencia de un feature de UniProt.

    Args:
        feature: Dict de un feature de UniProt (JSON).

    Returns:
        Tuple de (categoría, lista_de_codigos_ECO).
        Categoría: "EXPERIMENTAL", "BY_SIMILARITY", "PREDICTED", o "PREDICTED (No Evidence)".
    """
    evidences = feature.get("evidences", [])

    if not evidences:
        return "PREDICTED (No Evidence)", []

    eco_codes = []
    best_category = "PREDICTED"

    for ev in evidences:
        code = ev.get("evidenceCode", "")
        if code:
            eco_codes.append(code)

        category = ECO_CATEGORIES.get(code, "PREDICTED")
        # Quedarse con la mejor categoría encontrada
        if EVIDENCE_PRIORITY.get(category, 99) < EVIDENCE_PRIORITY.get(best_category, 99):
            best_category = category

    return best_category, eco_codes


def get_evidence_priority(category: str) -> int:
    """Retorna prioridad numérica de una categoría (menor = mejor)."""
    return EVIDENCE_PRIORITY.get(category, 99)
