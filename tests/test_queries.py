"""Quick test script to debug UniProt queries for viral signal peptides.
Run this directly: python test_queries.py
"""

import requests
import json

BASE = "https://rest.uniprot.org/uniprotkb/search"

queries = {
    # Test 1: Flavivirus - search by family Flaviviridae
    "flavivirus_family_11050": "(taxonomy_id:11050) AND (ft_signal:*) AND (reviewed:true)",
    
    # Test 2: Flavivirus - search by polyprotein name
    "flavivirus_polyprotein": '(taxonomy_id:11050) AND (protein_name:"Genome polyprotein") AND (reviewed:true)',
    
    # Test 3: Dengue specifically
    "dengue_signal": "(taxonomy_id:12637) AND (ft_signal:*)",
    
    # Test 4: Dengue polyprotein
    "dengue_polyprotein": '(taxonomy_id:12637) AND (protein_name:"Genome polyprotein")',
    
    # Test 5: Broad flavivirus - just polyprotein with signal
    "flavi_broad": '(taxonomy_id:11050) AND (protein_name:"polyprotein")',
    
    # Test 6: Flavivirus genus 11051
    "flavivirus_genus_11051": "(taxonomy_id:11051) AND (ft_signal:*)",
    
    # Test 7: Alphavirus structural polyprotein
    "alphavirus_structural": "(taxonomy_id:11019) AND (ft_signal:*) AND (reviewed:true)",
    
    # Test 8: Alphavirus broad
    "alphavirus_broad": '(taxonomy_id:11019) AND (protein_name:"polyprotein") AND (reviewed:true)',
    
    # Test 9: Influenza NA transmembrane
    "influenza_NA_transmem": "(taxonomy_id:11320 OR taxonomy_id:11520) AND (protein_name:neuraminidase) AND (ft_transmem:*)",
    
    # Test 10: Influenza NA broad
    "influenza_NA_broad": "(taxonomy_id:11320 OR taxonomy_id:11520) AND (protein_name:neuraminidase) AND (reviewed:true)",
}

print("Testing UniProt queries...\n")
print(f"{'Query name':<35} {'Results':>8}  Query string")
print("-" * 120)

for name, query in queries.items():
    try:
        url = f"{BASE}?query={query}&format=json&size=1"
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            data = r.json()
            # Try to get total from headers or response
            results = data.get("results", [])
            # Count is not always in JSON, check headers
            total = len(results)
            # A better way: check the x-total-results header
            total_header = r.headers.get("x-total-results", "?")
            print(f"{name:<35} {total_header:>8}  {query[:80]}")
        else:
            print(f"{name:<35} {'ERROR':>8}  HTTP {r.status_code}")
    except Exception as e:
        print(f"{name:<35} {'ERROR':>8}  {e}")

print("\nDone!")
