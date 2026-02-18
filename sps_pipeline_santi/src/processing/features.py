"""
Features computados para péptidos señal.

Calcula hidrofobicidad (Kyte-Doolittle), carga neta a pH 7,
y límites aproximados de n-region, h-region y c-region.
"""

# Escala Kyte-Doolittle de hidrofobicidad
KYTE_DOOLITTLE = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}

# Carga de aminoácidos a pH ~7
CHARGE_PH7 = {
    "K": +1.0,  # Lys pKa ~10.5
    "R": +1.0,  # Arg pKa ~12.5
    "H": +0.1,  # His pKa ~6.0 (parcialmente protonado a pH 7)
    "D": -1.0,  # Asp pKa ~3.7
    "E": -1.0,  # Glu pKa ~4.1
}


def hydrophobicity_mean(sequence: str) -> float:
    """Calcula la hidrofobicidad media Kyte-Doolittle de una secuencia."""
    if not sequence:
        return 0.0
    values = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in sequence.upper()]
    return round(sum(values) / len(values), 3)


def net_charge_ph7(sequence: str) -> float:
    """Calcula la carga neta estimada a pH 7."""
    if not sequence:
        return 0.0
    charge = sum(CHARGE_PH7.get(aa, 0.0) for aa in sequence.upper())
    return round(charge, 2)


def hydrophobicity_profile(sequence: str, window: int = 5) -> list[float]:
    """
    Calcula perfil de hidrofobicidad con ventana deslizante.

    Returns:
        Lista de valores de hidrofobicidad media por ventana.
    """
    seq = sequence.upper()
    if len(seq) < window:
        return [hydrophobicity_mean(seq)]

    profile = []
    for i in range(len(seq) - window + 1):
        window_seq = seq[i : i + window]
        values = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in window_seq]
        profile.append(sum(values) / len(values))

    return profile


def estimate_sp_regions(sequence: str) -> dict:
    """
    Estima los límites aproximados de las 3 regiones del péptido señal:
    - n-region: N-terminal positivamente cargada
    - h-region: Núcleo hidrofóbico
    - c-region: C-terminal con sitio de corte

    Usa heurística basada en perfil de hidrofobicidad y carga.

    Returns:
        Dict con 'n_region', 'h_region', 'c_region', cada uno con 'start' y 'end' (1-based).
    """
    seq = sequence.upper()
    length = len(seq)

    if length < 5:
        return {
            "n_region": {"start": 1, "end": length, "sequence": seq},
            "h_region": {"start": 0, "end": 0, "sequence": ""},
            "c_region": {"start": 0, "end": 0, "sequence": ""},
        }

    # Calcular perfil de hidrofobicidad
    profile = hydrophobicity_profile(seq, window=3)

    # Encontrar el inicio de la h-region: primer tramo donde hidrofobicidad > 1.0
    h_start = 0
    for i, val in enumerate(profile):
        if val > 1.0:
            h_start = i
            break
    else:
        # Si no hay región muy hidrofóbica, estimar
        h_start = min(3, length // 3)

    # Encontrar el fin de la h-region: último tramo donde hidrofobicidad > 1.0
    h_end = h_start
    for i in range(len(profile) - 1, h_start - 1, -1):
        if profile[i] > 1.0:
            h_end = i + 2  # +window-1 para compensar offset
            break

    # Asegurar que h_end no sobrepase la secuencia
    h_end = min(h_end, length - 2)

    # c-region: desde fin de h-region hasta el final
    c_start = h_end + 1
    if c_start >= length:
        c_start = length - 1

    # n-region: desde el inicio hasta el inicio de h-region
    n_end = h_start

    return {
        "n_region": {
            "start": 1,
            "end": max(1, n_end),
            "sequence": seq[:max(1, n_end)],
        },
        "h_region": {
            "start": max(1, n_end + 1),
            "end": h_end + 1,
            "sequence": seq[n_end : h_end + 1] if h_end >= n_end else "",
        },
        "c_region": {
            "start": h_end + 2,
            "end": length,
            "sequence": seq[h_end + 1 :] if h_end + 1 < length else "",
        },
    }


def compute_features(sp_sequence: str) -> dict:
    """
    Calcula todos los features para un péptido señal dado.

    Returns:
        Dict con hydrophobicity_mean, net_charge_ph7, y regiones estimadas.
    """
    regions = estimate_sp_regions(sp_sequence)

    return {
        "hydrophobicity_mean": hydrophobicity_mean(sp_sequence),
        "net_charge_ph7": net_charge_ph7(sp_sequence),
        "sp_length": len(sp_sequence),
        "n_region": regions["n_region"]["sequence"],
        "h_region": regions["h_region"]["sequence"],
        "c_region": regions["c_region"]["sequence"],
    }
