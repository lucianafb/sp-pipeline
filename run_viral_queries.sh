#!/bin/bash
# Viral queries for SP-Pipeline
# Requires: export SP_PIPELINE_NCBI_EMAIL="your@email.com"

set -e  # Exit on error

echo "=================================================="
echo "SP-Pipeline - Viral Queries"
echo "=================================================="
echo ""

# Check that NCBI is configured
if [ -z "$SP_PIPELINE_NCBI_EMAIL" ]; then
    echo "❌ ERROR: NCBI not configured"
    echo ""
    echo "Please run first:"
    echo "  export SP_PIPELINE_NCBI_EMAIL=\"your@email.com\""
    echo ""
    exit 1
fi

echo "✓ NCBI configured: $SP_PIPELINE_NCBI_EMAIL"
echo ""

# Check available sources
echo "Checking data sources..."
sp-pipeline check
echo ""

# ============================================
# ALPHAVIRUS
# ============================================
echo "=================================================="
echo "1. ALPHAVIRUS"
echo "=================================================="
echo "Parameters:"
echo "  - Evidence: all"
echo "  - Mode: exhaustive (all available)"
echo "  - FASTA: no"
echo "  - Output: output/viral_queries.csv"
echo ""

sp-pipeline query \
  --preset alphavirus_E3 \
  --evidence all \
  --mode exhaustive \
  -o output/viral_queries.csv

echo "✓ Alphavirus completed"
echo ""

# ============================================
# FLAVIVIRUS
# ============================================
echo "=================================================="
echo "2. FLAVIVIRUS"
echo "=================================================="
echo "Parameters:"
echo "  - Evidence: all"
echo "  - Mode: exhaustive (all available)"
echo "  - FASTA: no"
echo "  - Output: output/viral_queries.csv (append)"
echo ""

sp-pipeline query \
  --preset flavivirus \
  --evidence all \
  --mode exhaustive \
  --append \
  -o output/viral_queries.csv

echo "✓ Flavivirus completed"
echo ""

# ============================================
# INFLUENZA HA
# ============================================
echo "=================================================="
echo "3. INFLUENZA HA (Hemagglutinin)"
echo "=================================================="
echo "Parameters:"
echo "  - Evidence: all"
echo "  - Mode: exhaustive (all available)"
echo "  - FASTA: no"
echo "  - Output: output/viral_queries.csv (append)"
echo ""

sp-pipeline query \
  --preset influenza_HA \
  --evidence all \
  --mode exhaustive \
  --append \
  -o output/viral_queries.csv

echo "✓ Influenza HA completed"
echo ""

# ============================================
# INFLUENZA NA
# ============================================
echo "=================================================="
echo "4. INFLUENZA NA (Neuraminidase)"
echo "=================================================="
echo "Parameters:"
echo "  - Evidence: all"
echo "  - Mode: representative (curated subset)"
echo "  - FASTA: no"
echo "  - Output: output/viral_queries.csv (append)"
echo ""

sp-pipeline query \
  --preset influenza_NA \
  --evidence all \
  --mode representative \
  --append \
  -o output/viral_queries.csv

echo "✓ Influenza NA completed"
echo ""

# ============================================
# SUMMARY
# ============================================
echo "=================================================="
echo "✅ ALL QUERIES COMPLETED"
echo "=================================================="
echo ""
echo "Output file: output/viral_queries.csv"
echo ""

# Verify results
if [ -f "output/viral_queries.csv" ]; then
    TOTAL_LINES=$(wc -l < output/viral_queries.csv)
    TOTAL_RECORDS=$((TOTAL_LINES - 1))  # Subtract header
    echo "Total records: $TOTAL_RECORDS"
    echo ""
    
    # Show statistics by query_group
    echo "Records by preset:"
    python3 -c "
import pandas as pd
df = pd.read_csv('output/viral_queries.csv')
print(df['query_group'].value_counts().to_string())
print()
print('Records by source:')
print(df['source_db'].value_counts().to_string())
print()
print('Records by SP type:')
print(df['sp_type'].value_counts().to_string())
" 2>/dev/null || echo "  (pandas not available for detailed statistics)"
    
    echo ""
    echo "To verify results:"
    echo "  python3 verify_pipeline.py output/viral_queries.csv"
else
    echo "❌ Error: Output file was not generated"
fi
