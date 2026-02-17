#!/usr/bin/env python3
"""Complete SP-Pipeline verification and testing guide.

This script verifies that SP-Pipeline is working correctly and shows
examples of how to use and customize queries.

Usage:
  python verify_pipeline.py                    # Show all tests & guides
  python verify_pipeline.py output/results.csv # Verify a specific CSV
  python verify_pipeline.py --quick            # Quick test of all presets
  python verify_pipeline.py --examples         # Show custom query examples
"""

import sys
import subprocess
import argparse
import pandas as pd
from pathlib import Path
from datetime import datetime


def run_command(cmd: list[str], timeout: int = 120) -> tuple[int, str, str]:
    """Run shell command and return (returncode, stdout, stderr)."""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return 1, "", f"Command timed out after {timeout}s"
    except Exception as e:
        return 1, "", str(e)


def test_preset(preset_name: str, output_dir: str = "output") -> dict:
    """Test a preset and return results."""
    output_file = Path(output_dir) / f"test_{preset_name}.csv"
    Path(output_dir).mkdir(exist_ok=True)

    print(f"\n  Testing {preset_name}...", end=" ", flush=True)

    cmd = ["sp-pipeline", "query", "--preset", preset_name, "-o", str(output_file)]
    returncode, stdout, stderr = run_command(cmd)

    result = {
        "preset": preset_name,
        "status": "✓" if returncode == 0 else "✗",
        "record_count": 0,
    }

    if returncode == 0 and output_file.exists():
        try:
            df = pd.read_csv(output_file)
            result["record_count"] = len(df)
            print(f"{len(df)} records")
        except Exception as e:
            print(f"Error reading CSV: {e}")
            result["status"] = "✗"
    else:
        print("Failed")

    return result


def verify_csv_structure(csv_path: str) -> dict:
    """Verify the structure of a CSV output file."""
    results = {"errors": [], "warnings": [], "info": []}

    try:
        df = pd.read_csv(csv_path)
        results["info"].append(f"✓ CSV loaded: {len(df)} records")
    except Exception as e:
        results["errors"].append(f"✗ Error loading CSV: {e}")
        return results

    # Check expected columns
    expected_columns = [
        'entry_id', 'source_db', 'accession', 'organism', 'taxonomy_id',
        'protein_name', 'gene_name', 'full_sequence', 'sp_sequence',
        'sp_start', 'sp_end', 'sp_length', 'sp_type', 'cleavage_site_motif',
        'evidence_method', 'viral_subtype', 'viral_host', 'query_group',
        'date_retrieved'
    ]

    missing_cols = set(expected_columns) - set(df.columns)
    if missing_cols:
        results["errors"].append(f"✗ Missing columns: {missing_cols}")
    else:
        results["info"].append(f"✓ All expected columns present")

    # Check for empty query_group (critical bug indicator)
    empty_query_groups = df['query_group'].isna().sum()
    if empty_query_groups > 0:
        results["errors"].append(f"✗ {empty_query_groups} records with empty query_group")
    else:
        results["info"].append("✓ All records have query_group assigned")

    # Verify data consistency
    if 'sp_start' in df.columns and 'sp_end' in df.columns and 'sp_length' in df.columns:
        df['calculated_length'] = df['sp_end'] - df['sp_start'] + 1
        mismatches = (df['sp_length'] != df['calculated_length']).sum()
        if mismatches > 0:
            results["warnings"].append(f"⚠ {mismatches} records with inconsistent sp_length")
        else:
            results["info"].append("✓ sp_length consistent with sp_start/sp_end")

    # Check sp_sequence length
    if 'sp_sequence' in df.columns and 'sp_length' in df.columns:
        df['seq_len'] = df['sp_sequence'].str.len()
        seq_mismatches = (df['sp_length'] != df['seq_len']).sum()
        if seq_mismatches > 0:
            results["warnings"].append(f"⚠ {seq_mismatches} records with len(sp_sequence) != sp_length")

    # Statistics
    results["info"].append("\n📊 Statistics:")
    results["info"].append(f"  Total records: {len(df)}")

    if 'query_group' in df.columns:
        groups = df['query_group'].value_counts()
        results["info"].append(f"  Query groups: {', '.join([f'{k}: {v}' for k, v in groups.items()])}")

    if 'sp_type' in df.columns:
        types = df['sp_type'].value_counts()
        results["info"].append(f"  SP types: {', '.join([f'{k}: {v}' for k, v in types.items()])}")

    if 'sp_length' in df.columns:
        results["info"].append(f"  SP length:")
        results["info"].append(f"    Mean: {df['sp_length'].mean():.1f} aa")
        results["info"].append(f"    Median: {df['sp_length'].median():.1f} aa")
        results["info"].append(f"    Range: {df['sp_length'].min()}-{df['sp_length'].max()} aa")

    return results


def show_custom_queries():
    """Show custom query examples for users."""
    print("\n" + "="*70)
    print("CUSTOM QUERY EXAMPLES - How to customize UniProt searches")
    print("="*70)

    examples = [
        {
            "title": "Curated Alphavirus E3 signal peptides (base query)",
            "query": "(taxonomy_id:11019) AND (protein_name:\"Structural polyprotein\") AND (reviewed:true)",
            "explanation": [
                "taxonomy_id:11019 = Alphavirus genus",
                "protein_name:\"Structural polyprotein\" = structural polyprotein only",
                "reviewed:true = curated entries only",
            ],
        },
        {
            "title": "Chikungunya virus (CHIKV) E3 specifically",
            "query": "(taxonomy_id:11034) AND (protein_name:\"Structural polyprotein\") AND (reviewed:true)",
            "explanation": [
                "taxonomy_id:11034 = Chikungunya virus (specific species)",
                "Can be used at https://www.uniprot.org/ for manual search",
            ],
        },
        {
            "title": "Dengue virus Flavivirus",
            "query": "(taxonomy_id:12637) AND (protein_name:\"Genome polyprotein\") AND (reviewed:true)",
            "explanation": [
                "taxonomy_id:12637 = Dengue virus (all 4 serotypes)",
                "Parser automatically extracts SPs from internal cleavage sites",
            ],
        },
        {
            "title": "All Alphavirus variants (including unreviewed)",
            "query": "(taxonomy_id:11019) AND (protein_name:\"Structural polyprotein\")",
            "explanation": [
                "Removed reviewed:true to include computational predictions",
                "Warning: may return 1000+ entries",
            ],
        },
    ]

    for i, ex in enumerate(examples, 1):
        print(f"\n[{i}] {ex['title']}")
        print(f"    Query: {ex['query']}")
        print(f"    Explanation:")
        for line in ex['explanation']:
            print(f"      - {line}")

    print("\n" + "="*70)
    print("ALPHAVIRUS TAXONOMY IDs")
    print("="*70)

    alpha_taxids = [
        ("11019", "Alphavirus (entire genus)"),
        ("11034", "Chikungunya virus (CHIKV)"),
        ("11036", "Sindbis virus (SINV)"),
        ("11029", "Semliki Forest virus (SFV)"),
        ("11049", "Venezuelan equine encephalitis virus (VEEV)"),
    ]

    for taxid, name in alpha_taxids:
        print(f"  {taxid}: {name}")

    print("\n" + "="*70)
    print("FLAVIVIRUS TAXONOMY IDs")
    print("="*70)

    flavi_taxids = [
        ("12637", "Dengue virus"),
        ("64320", "Zika virus"),
        ("11082", "West Nile virus"),
        ("11053", "Yellow fever virus"),
        ("11080", "Japanese encephalitis virus"),
        ("11090", "Tick-borne encephalitis virus"),
    ]

    for taxid, name in flavi_taxids:
        print(f"  {taxid}: {name}")


def show_customization_guide():
    """Show how to customize presets and queries."""
    print("\n" + "="*70)
    print("CUSTOMIZATION GUIDE")
    print("="*70)

    print("\n[1] MODIFY EVIDENCE LEVEL")
    print("    Default: curated (reviewed:true)")
    print("    Options: --evidence experimental | curated_annotation | all")
    print("    Example:")
    print("      sp-pipeline query --preset alphavirus_E3 --evidence all \\")
    print("        -o output/alphavirus_all.csv")

    print("\n[2] MODIFY SEARCH MODE (viral presets)")
    print("    Default: representative (best curated entries)")
    print("    Options: --mode exhaustive (all available)")
    print("    Example:")
    print("      sp-pipeline query --preset flavivirus --mode exhaustive \\")
    print("        -o output/flavivirus_all.csv")

    print("\n[3] EXPORT FASTA")
    print("    Export signal peptide sequences in FASTA format")
    print("    Example:")
    print("      sp-pipeline query --preset human_type1 \\")
    print("        -o output/human_type1.csv --fasta output/human_type1.fasta")

    print("\n[4] APPEND TO EXISTING CSV")
    print("    Add new results without duplicating")
    print("    Example:")
    print("      sp-pipeline query --preset alphavirus_E3 \\")
    print("        -o output/results.csv --append")

    print("\n[5] CLEAR CACHE")
    print("    Force fresh queries from API")
    print("    Command: sp-pipeline clear-cache")

    print("\n[6] VERIFY SOURCES")
    print("    Check which data sources are available")
    print("    Commands:")
    print("      sp-pipeline check     # Check all sources")
    print("      sp-pipeline presets   # List available presets\n")


def show_quick_test():
    """Run quick tests of all presets."""
    print("\n" + "="*70)
    print("QUICK PRESET TEST")
    print("="*70)

    presets = [
        "human_type1", "human_type2",
        "alphavirus_E3", "flavivirus", "influenza_HA", "influenza_NA"
    ]

    results = []
    for preset in presets:
        result = test_preset(preset)
        results.append(result)

    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70 + "\n")

    print(f"{'Preset':<20} {'Status':<10} {'Records':<10}")
    print("-" * 40)
    for r in results:
        print(f"{r['preset']:<20} {r['status']:<10} {r['record_count']:<10}")

    # Show expected ranges
    print("\n" + "="*70)
    print("EXPECTED RESULTS (v0.2.0+)")
    print("="*70 + "\n")

    expected = {
        "human_type1": "500+",
        "human_type2": "400+",
        "alphavirus_E3": "50+",
        "flavivirus": "50+",
        "influenza_HA": "200+",
        "influenza_NA": "100+",
    }

    for r in results:
        exp = expected.get(r['preset'], "?")
        status = "✓" if r['record_count'] > 0 else "✗"
        print(f"{status} {r['preset']:<20} expected {exp:<10} got {r['record_count']}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="SP-Pipeline verification and testing guide"
    )
    parser.add_argument(
        "csv_file",
        nargs="?",
        help="CSV file to verify (optional)"
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Quick test of all presets"
    )
    parser.add_argument(
        "--examples",
        action="store_true",
        help="Show custom query examples only"
    )

    args = parser.parse_args()

    # Show header
    if not args.examples:
        print("""
╔═══════════════════════════════════════════════════════════════════════╗
║                   SP-PIPELINE VERIFICATION SUITE                      ║
║                          Version 0.2.0+                               ║
╚═══════════════════════════════════════════════════════════════════════╝
        """)

    # Quick test mode
    if args.quick:
        show_quick_test()
        show_customization_guide()
        return

    # Custom queries example mode
    if args.examples:
        show_custom_queries()
        show_customization_guide()
        return

    # Verify CSV mode
    if args.csv_file:
        csv_path = args.csv_file
        if not Path(csv_path).exists():
            print(f"❌ Error: File '{csv_path}' not found")
            sys.exit(1)

        results = verify_csv_structure(csv_path)

        # Print results
        print("\n" + "="*70)
        print("VERIFICATION REPORT")
        print("="*70 + "\n")

        if results["errors"]:
            print("❌ CRITICAL ERRORS:")
            for error in results["errors"]:
                print(f"  {error}")
            print()

        if results["warnings"]:
            print("⚠️  WARNINGS:")
            for warning in results["warnings"]:
                print(f"  {warning}")
            print()

        if results["info"]:
            print("ℹ️  INFORMATION:")
            for info in results["info"]:
                print(f"  {info}")
            print()

        # Summary
        print("="*70)
        if results["errors"]:
            print("❌ VERIFICATION FAILED - Critical errors found")
            sys.exit(1)
        elif results["warnings"]:
            print("⚠️  VERIFICATION COMPLETED WITH WARNINGS")
        else:
            print("✅ VERIFICATION SUCCESSFUL - Pipeline working correctly")
        return

    # Default: show all guides
    show_custom_queries()
    show_customization_guide()
    print("\n" + "="*70)
    print("TO RUN A QUICK TEST")
    print("="*70)
    print("\n  python verify_pipeline.py --quick\n")


if __name__ == "__main__":
    main()
