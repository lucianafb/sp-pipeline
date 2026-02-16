# SP-Pipeline

Automated pipeline for querying signal peptide data from curated biological databases.

## Overview

SP-Pipeline retrieves, standardizes, and exports signal peptide (SP) annotations from multiple databases into a unified CSV format. It supports both curated/experimental annotations and computational predictions, with built-in presets for common queries.

## Features

- **Multiple data sources**: UniProt/SwissProt, NCBI/GenBank, NCBI Influenza Virus Resource
- **Preset queries**: Ready-to-use queries for human SPs (type 1 & 2), influenza HA/NA, alphavirus E3, and flavivirus SPs
- **Configurable evidence filtering**: Toggle between experimental, predicted, or all annotations
- **Exhaustive vs. representative modes**: For viral queries, choose between all available sequences or representative subsets
- **De novo predictions** (optional): Run SignalP and Phobius on retrieved sequences
- **Idempotent exports**: Re-run the pipeline without duplicating records
- **Cache**: Local caching to avoid redundant API calls
- **FASTA export**: Optional export of SP sequences in FASTA format

## Installation

```bash
# Clone the repository
git clone https://github.com/your-org/sp-pipeline.git
cd sp-pipeline

# Install in development mode
pip install -e .

# Or install directly
pip install .
```

### Requirements

- Python >= 3.9
- Internet access (for API queries)

## Quick Start

```bash
# List available presets
sp-pipeline presets

# Query human type 1 signal peptides
sp-pipeline query --preset human_type1 --output output/human_sp_type1.csv

# Query human type 1 and type 2 together
sp-pipeline query --preset human_type1 --preset human_type2 --output output/human_sp_all.csv

# Query influenza HA signal peptides (representative mode)
sp-pipeline query --preset influenza_HA --mode representative --output output/influenza_ha.csv

# Query influenza HA (exhaustive - all available sequences)
sp-pipeline query --preset influenza_HA --mode exhaustive --output output/influenza_ha_all.csv

# Query with experimental evidence only
sp-pipeline query --preset human_type1 --evidence experimental --output output/human_exp.csv

# Query flavivirus SPs including internal signals
sp-pipeline query --preset flavivirus --output output/flavivirus_sp.csv

# Append new results to existing CSV
sp-pipeline query --preset alphavirus_E3 --output output/results.csv --append

# Also export FASTA
sp-pipeline query --preset human_type1 --output output/human.csv --fasta

# Custom query
sp-pipeline custom --organism "Mus musculus" --sp-type type1_cleaved --output output/mouse.csv
```

## Available Presets

| Preset | Description |
|--------|-------------|
| `human_type1` | Human type 1 (cleaved) signal peptides |
| `human_type2` | Human type 2 signal anchors (uncleaved) |
| `influenza_HA` | Influenza hemagglutinin SPs (all subtypes, A/B/C/D) |
| `influenza_NA` | Influenza neuraminidase SPs (all subtypes) |
| `alphavirus_E3` | Alphavirus E3/pE2 signal peptides |
| `flavivirus` | Flavivirus SPs (prM + internal signals) |

## Output Format

The pipeline produces a CSV with the following columns:

| Column | Description |
|--------|-------------|
| `entry_id` | Database entry ID (e.g., P01234) |
| `source_db` | Source database (UniProt, GenBank, NCBI) |
| `accession` | Secondary accession |
| `organism` | Organism scientific name |
| `taxonomy_id` | NCBI Taxonomy ID |
| `protein_name` | Protein name |
| `gene_name` | Gene name |
| `full_sequence` | Complete protein sequence |
| `sp_sequence` | Signal peptide sequence |
| `sp_start` | SP start position (1-based) |
| `sp_end` | SP end position (1-based, inclusive) |
| `sp_length` | Signal peptide length |
| `sp_type` | Classification: `type1_cleaved`, `type2_signal_anchor`, `internal_signal` |
| `cleavage_site_motif` | Cleavage site motif (e.g., AFA-QP) |
| `evidence_method` | Evidence: `experimental`, `curated_annotation`, `predicted_signalp`, etc. |
| `viral_subtype` | Viral subtype (e.g., H1N1, DENV-2, CHIKV) |
| `viral_host` | Viral host organism |
| `query_group` | Preset name that retrieved this record |
| `date_retrieved` | Date of data retrieval |

## Configuration

### NCBI Access

To use NCBI as a data source, you need to provide an email address:

```bash
# Option 1: Environment variable
export SP_PIPELINE_NCBI_EMAIL="your.email@institution.edu"

# Option 2: Config file
# Edit config/default_config.yaml or create ~/.config/sp-pipeline/config.yaml
```

For higher rate limits (10 req/s instead of 3), get an API key from [NCBI](https://www.ncbi.nlm.nih.gov/account/settings/):

```bash
export SP_PIPELINE_NCBI_API_KEY="your_api_key_here"
```

### Custom Configuration

Create `~/.config/sp-pipeline/config.yaml` to override defaults:

```yaml
ncbi:
  email: "your.email@institution.edu"
  api_key: "your_key"

cache:
  enabled: true
  ttl_days: 60

query:
  default_evidence: "experimental"
  default_mode: "representative"

output:
  include_fasta: true
```

Or pass a config file directly:

```bash
sp-pipeline query --preset human_type1 --config my_config.yaml
```

## CLI Commands

| Command | Description |
|---------|-------------|
| `sp-pipeline query` | Run preset queries |
| `sp-pipeline custom` | Run a custom query |
| `sp-pipeline presets` | List available presets |
| `sp-pipeline check` | Check data source availability |
| `sp-pipeline clear-cache` | Clear the local query cache |

## Project Structure

```
sp-pipeline/
├── config/
│   └── default_config.yaml       # Default configuration
├── src/sp_pipeline/
│   ├── cli.py                     # Command-line interface
│   ├── config.py                  # Configuration management
│   ├── models.py                  # Data models (SPRecord)
│   ├── pipeline.py                # Main orchestrator
│   ├── dedup.py                   # Deduplication
│   ├── exporters.py               # CSV/FASTA export
│   ├── utils.py                   # Utilities (cache, rate limiting)
│   ├── sources/                   # Data source modules
│   │   ├── uniprot.py             # UniProt API
│   │   └── ncbi.py                # NCBI Entrez API
│   └── predictors/                # Prediction wrappers (Phase 3)
│       ├── signalp.py
│       └── phobius.py
├── notebooks/
│   └── sp_pipeline_colab.ipynb    # Google Colab interface
├── tests/
├── output/                        # Generated results
└── cache/                         # Local query cache
```

## Development Roadmap

- [x] **Phase 1**: Core pipeline + UniProt source + human presets
- [x] **Phase 2**: NCBI source + viral presets (influenza, alphavirus, flavivirus)
- [ ] **Phase 3**: SignalP/Phobius prediction wrappers
- [ ] **Phase 4**: Google Colab notebook, extended tests, documentation

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Install in dev mode: `pip install -e ".[dev]"`
4. Run tests: `pytest`
5. Submit a pull request

## License

MIT License

## Citation

If you use SP-Pipeline in your research, please cite:

```
[Citation pending]
```
