"""
Script de diagnóstico para encontrar las queries correctas de UniProt
para los presets virales: flavivirus, alphavirus_E3, influenza_NA.

Uso:
    python scripts/test_viral_queries.py
"""

import json
import requests

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,id,protein_name,organism_name,organism_id,ft_signal,ft_transmem"


def count_results(query: str, reviewed_only: bool = False) -> tuple[int, list[dict]]:
    """Fetch first page of results and return total + sample entries."""
    if reviewed_only:
        query = f"({query}) AND (reviewed:true)"
    url = f"{UNIPROT_API}?query={query}&fields={FIELDS}&format=json&size=5"
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        data = r.json()
        # Total count is in the X-Total-Results header
        total = int(r.headers.get("X-Total-Results", 0))
        return total, data.get("results", [])
    except Exception as e:
        print(f"  ERROR: {e}")
        return 0, []


def show_sample(entries: list[dict], n: int = 3):
    for entry in entries[:n]:
        acc = entry.get("primaryAccession", "")
        org = entry.get("organism", {}).get("scientificName", "")
        protein_desc = entry.get("proteinDescription", {})
        rec_name = protein_desc.get("recommendedName", {})
        name = rec_name.get("fullName", {}).get("value", "")
        if not name:
            subs = protein_desc.get("submissionNames", [])
            name = subs[0].get("fullName", {}).get("value", "") if subs else "?"
        has_signal = bool(entry.get("features") and
                          any(f.get("type") == "Signal" for f in entry.get("features", [])))
        has_transmem = bool(entry.get("features") and
                            any(f.get("type") == "Transmembrane" for f in entry.get("features", [])))
        feats = []
        if has_signal:
            feats.append("ft_signal")
        if has_transmem:
            feats.append("ft_transmem")
        print(f"    [{acc}] {name[:60]} | {org[:40]} | {', '.join(feats) or 'no features'}")


def test_section(title: str):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


# ──────────────────────────────────────────────────────────────────────────────
# FLAVIVIRUS
# ──────────────────────────────────────────────────────────────────────────────
test_section("FLAVIVIRUS — Taxonomy IDs")

flavi_candidates = [
    ("11050", "Flaviviridae (family)"),
    ("11051", "Orthoflavivirus (genus, former Flavivirus)"),
    ("11052", "Dengue virus"),
    ("12637", "Dengue virus 1"),
    ("11053", "Yellow fever virus"),
    ("11082", "West Nile virus"),
    ("64320", "Zika virus"),
    ("11103", "Hepatitis C virus (not flavivirus proper but same family)"),
]

for taxid, label in flavi_candidates:
    q = f"(taxonomy_id:{taxid}) AND (ft_signal:*)"
    total, sample = count_results(q)
    print(f"  taxid:{taxid} ({label}): {total} results with ft_signal")
    if total > 0:
        # Also check with reviewed only
        total_rev, _ = count_results(q, reviewed_only=False)
        show_sample(sample, 2)

test_section("FLAVIVIRUS — Protein name searches")

flavi_protein_queries = [
    ("Genome polyprotein", "ft_signal:* AND protein_name:\"Genome polyprotein\""),
    ("Polyprotein any taxid 11051", "(taxonomy_id:11051) AND (protein_name:polyprotein)"),
    ("Signal in polyprotein taxid 11051", "(taxonomy_id:11051) AND (ft_signal:*) AND (protein_name:polyprotein)"),
    ("Reviewed polyprotein 11051", "(taxonomy_id:11051) AND (protein_name:\"Genome polyprotein\") AND (reviewed:true)"),
    ("Dengue polyprotein", "(taxonomy_id:12637) AND (protein_name:\"Genome polyprotein\")"),
    ("Dengue ft_signal", "(taxonomy_id:12637) AND (ft_signal:*)"),
    ("Zika ft_signal", "(taxonomy_id:64320) AND (ft_signal:*)"),
    ("WNV ft_signal", "(taxonomy_id:11082) AND (ft_signal:*)"),
    ("YFV ft_signal", "(taxonomy_id:11053) AND (ft_signal:*)"),
]

for label, query in flavi_protein_queries:
    total, sample = count_results(query)
    print(f"\n  [{label}]: {total} results")
    print(f"  Query: {query}")
    if total > 0:
        show_sample(sample, 2)

# ──────────────────────────────────────────────────────────────────────────────
# ALPHAVIRUS
# ──────────────────────────────────────────────────────────────────────────────
test_section("ALPHAVIRUS — Taxonomy IDs")

alpha_candidates = [
    ("11019", "Alphavirus (genus)"),
    ("11034", "Chikungunya virus"),
    ("11036", "Sindbis virus"),
    ("11029", "Semliki Forest virus"),
    ("11049", "Venezuelan equine encephalitis virus"),
]

for taxid, label in alpha_candidates:
    q = f"(taxonomy_id:{taxid}) AND (ft_signal:*)"
    total, sample = count_results(q)
    print(f"  taxid:{taxid} ({label}): {total} results with ft_signal")
    if total > 0:
        show_sample(sample, 2)

test_section("ALPHAVIRUS — Protein name searches")

alpha_queries = [
    ("Structural polyprotein ft_signal 11019", "(taxonomy_id:11019) AND (protein_name:\"Structural polyprotein\") AND (ft_signal:*)"),
    ("Structural polyprotein 11019", "(taxonomy_id:11019) AND (protein_name:\"Structural polyprotein\")"),
    ("Polyprotein ft_signal 11019", "(taxonomy_id:11019) AND (protein_name:polyprotein) AND (ft_signal:*)"),
    ("CHIKV structural polyprotein", "(taxonomy_id:11034) AND (protein_name:polyprotein)"),
    ("CHIKV ft_signal", "(taxonomy_id:11034) AND (ft_signal:*)"),
    ("Sindbis ft_signal", "(taxonomy_id:11036) AND (ft_signal:*)"),
    ("SFV ft_signal", "(taxonomy_id:11029) AND (ft_signal:*)"),
    ("E3 protein 11019", "(taxonomy_id:11019) AND (protein_name:E3)"),
    ("pE2 protein 11019", "(taxonomy_id:11019) AND (protein_name:pE2)"),
]

for label, query in alpha_queries:
    total, sample = count_results(query)
    print(f"\n  [{label}]: {total} results")
    print(f"  Query: {query}")
    if total > 0:
        show_sample(sample, 2)

# ──────────────────────────────────────────────────────────────────────────────
# INFLUENZA NEURAMINIDASE (Type II — no cleavable SP)
# ──────────────────────────────────────────────────────────────────────────────
test_section("INFLUENZA NEURAMINIDASE — ft_transmem searches")

na_queries = [
    ("NA all influenza ft_transmem", "(taxonomy_id:11308) AND (protein_name:neuraminidase) AND (ft_transmem:*)"),
    ("NA influenza A ft_transmem", "(taxonomy_id:11320) AND (protein_name:neuraminidase) AND (ft_transmem:*)"),
    ("NA influenza B ft_transmem", "(taxonomy_id:11520) AND (protein_name:neuraminidase) AND (ft_transmem:*)"),
    ("NA reviewed ft_transmem", "(taxonomy_id:11320 OR taxonomy_id:11520) AND (protein_name:neuraminidase) AND (ft_transmem:*) AND (reviewed:true)"),
    ("NA influenza A all", "(taxonomy_id:11320) AND (protein_name:neuraminidase)"),
    ("NA influenza A ft_signal (should be 0)", "(taxonomy_id:11320) AND (protein_name:neuraminidase) AND (ft_signal:*)"),
]

for label, query in na_queries:
    total, sample = count_results(query)
    print(f"\n  [{label}]: {total} results")
    print(f"  Query: {query}")
    if total > 0:
        show_sample(sample, 2)

test_section("SUMMARY — Recommended queries")
print("""
Based on the above results, update _build_query() in uniprot.py:

FLAVIVIRUS:
  Use the taxonomy IDs that return results above.
  Search for 'Genome polyprotein' entries which contain Signal features.

ALPHAVIRUS:
  Use (taxonomy_id:11019) AND (protein_name:"Structural polyprotein")
  or specific CHIKV/SINV/SFV taxids with ft_signal.

INFLUENZA NA:
  Use (taxonomy_id:11320 OR ...) AND (protein_name:neuraminidase) AND (ft_transmem:*)
  NOT ft_signal (NA has no cleavable signal peptide).
""")
