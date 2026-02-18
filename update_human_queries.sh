#!/bin/bash
# Update human signal peptide queries
# This script regenerates all human SP outputs with fresh data

set -e  # Exit on error

echo "=================================================="
echo "SP-Pipeline - Update Human Queries"
echo "=================================================="
echo ""

# Backup old files
echo "📦 Backing up old files..."
mkdir -p output/old_queries
if [ -f "output/results.csv" ]; then
    mv output/results.csv output/old_queries/results_$(date +%Y%m%d).csv
    echo "  ✓ Backed up results.csv"
fi
if [ -f "output/human_type1_exp.csv" ]; then
    mv output/human_type1_exp.csv output/old_queries/human_type1_exp_$(date +%Y%m%d).csv
    echo "  ✓ Backed up human_type1_exp.csv"
fi
if [ -f "output/human_type2.csv" ]; then
    mv output/human_type2.csv output/old_queries/human_type2_$(date +%Y%m%d).csv
    echo "  ✓ Backed up human_type2.csv"
fi
echo ""

# ============================================
# HUMAN TYPE 1 - ALL EVIDENCE
# ============================================
echo "=================================================="
echo "1. HUMAN TYPE 1 (All Evidence)"
echo "=================================================="
echo "Parameters:"
echo "  - Preset: human_type1"
echo "  - Evidence: all"
echo "  - Output: output/human_type1_all.csv"
echo ""

sp-pipeline query \
  --preset human_type1 \
  --evidence all \
  -o output/human_type1_all.csv

echo "✓ Human type 1 (all) completed"
echo ""

# ============================================
# HUMAN TYPE 1 - EXPERIMENTAL ONLY
# ============================================
echo "=================================================="
echo "2. HUMAN TYPE 1 (Experimental Only)"
echo "=================================================="
echo "Parameters:"
echo "  - Preset: human_type1"
echo "  - Evidence: experimental"
echo "  - Output: output/human_type1_exp.csv"
echo ""

sp-pipeline query \
  --preset human_type1 \
  --evidence experimental \
  -o output/human_type1_exp.csv

echo "✓ Human type 1 (experimental) completed"
echo ""

# ============================================
# HUMAN TYPE 1 - CURATED ONLY
# ============================================
echo "=================================================="
echo "3. HUMAN TYPE 1 (Curated Only)"
echo "=================================================="
echo "Parameters:"
echo "  - Preset: human_type1"
echo "  - Evidence: curated_annotation"
echo "  - Output: output/human_type1_curated.csv"
echo ""

sp-pipeline query \
  --preset human_type1 \
  --evidence curated_annotation \
  -o output/human_type1_curated.csv

echo "✓ Human type 1 (curated) completed"
echo ""

# ============================================
# HUMAN TYPE 2 - ALL EVIDENCE
# ============================================
echo "=================================================="
echo "4. HUMAN TYPE 2 (All Evidence)"
echo "=================================================="
echo "Parameters:"
echo "  - Preset: human_type2"
echo "  - Evidence: all"
echo "  - Output: output/human_type2.csv"
echo ""

sp-pipeline query \
  --preset human_type2 \
  --evidence all \
  -o output/human_type2.csv

echo "✓ Human type 2 completed"
echo ""

# ============================================
# COMBINED - TYPE 1 + TYPE 2
# ============================================
echo "=================================================="
echo "5. COMBINED (Type 1 + Type 2, Curated Only)"
echo "=================================================="
echo "Parameters:"
echo "  - Presets: human_type1 + human_type2"
echo "  - Evidence: curated_annotation"
echo "  - Output: output/human_combined_curated.csv"
echo ""

sp-pipeline query \
  --preset human_type1 \
  --preset human_type2 \
  --evidence curated_annotation \
  -o output/human_combined_curated.csv

echo "✓ Combined (curated) completed"
echo ""

# ============================================
# SUMMARY
# ============================================
echo "=================================================="
echo "✅ ALL HUMAN QUERIES COMPLETED"
echo "=================================================="
echo ""
echo "Output files generated:"
echo ""

# Show statistics
python3 -c "
import pandas as pd
import os

files = [
    ('human_type1_all.csv', 'Human type 1 (all evidence)'),
    ('human_type1_exp.csv', 'Human type 1 (experimental)'),
    ('human_type1_curated.csv', 'Human type 1 (curated)'),
    ('human_type2.csv', 'Human type 2 (all evidence)'),
    ('human_combined_curated.csv', 'Combined type 1+2 (curated)')
]

print('File                          | Records | Description')
print('-' * 75)

for fname, desc in files:
    path = f'output/{fname}'
    if os.path.exists(path):
        df = pd.read_csv(path)
        print(f'{fname:30} | {len(df):7} | {desc}')
    else:
        print(f'{fname:30} | NOT FOUND')

print()
print('To verify any file:')
print('  python3 verify_pipeline.py output/<filename>')
" 2>/dev/null || echo "  (pandas not available for statistics)"

echo ""
echo "Old files backed up to: output/old_queries/"
