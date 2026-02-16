"""Data models for the SP-Pipeline."""

from dataclasses import dataclass, field, asdict
from datetime import date
from typing import Optional


@dataclass
class SPRecord:
    """Standardized signal peptide record.

    This is the canonical data structure used throughout the pipeline.
    Every source module must convert its raw data into SPRecord instances.
    """

    # Identifiers
    entry_id: str                          # e.g., P01234, AAB12345
    source_db: str                         # UniProt, GenBank, NCBI, RefSeq
    accession: Optional[str] = None        # Secondary accession if available

    # Organism info
    organism: str = ""                     # e.g., "Homo sapiens"
    taxonomy_id: Optional[int] = None      # NCBI Taxonomy ID

    # Protein info
    protein_name: str = ""
    gene_name: Optional[str] = None

    # Sequences
    full_sequence: str = ""                # Full protein sequence
    sp_sequence: str = ""                  # Signal peptide sequence only

    # SP coordinates (1-based, inclusive)
    sp_start: int = 0
    sp_end: int = 0
    sp_length: int = 0

    # SP classification
    sp_type: str = ""                      # type1_cleaved, type2_signal_anchor, internal_signal
    cleavage_site_motif: Optional[str] = None  # e.g., "AXA", "VFA-AP"

    # Evidence
    evidence_method: str = ""              # experimental, predicted_signalp, predicted_phobius, curated_annotation

    # Viral-specific fields
    viral_subtype: Optional[str] = None    # e.g., H1, N2, DENV-2, CHIKV
    viral_host: Optional[str] = None       # e.g., human, avian, swine

    # Metadata
    query_group: str = ""                  # Preset name that retrieved this record
    date_retrieved: str = field(default_factory=lambda: date.today().isoformat())

    def __post_init__(self):
        """Compute derived fields."""
        if self.sp_sequence and not self.sp_length:
            self.sp_length = len(self.sp_sequence)
        if self.full_sequence and self.sp_end and not self.sp_sequence:
            self.sp_sequence = self.full_sequence[self.sp_start - 1 : self.sp_end]
            self.sp_length = len(self.sp_sequence)

    def to_dict(self) -> dict:
        """Convert to dictionary for DataFrame export."""
        return asdict(self)

    @property
    def unique_key(self) -> str:
        """Generate a unique key for deduplication.

        Uses entry_id + source + sp_start + sp_end to handle cases where
        one protein has multiple signal peptides (e.g., flavivirus internal signals).
        """
        return f"{self.source_db}|{self.entry_id}|{self.sp_start}|{self.sp_end}"


# CSV column order for export
CSV_COLUMNS = [
    "entry_id",
    "source_db",
    "accession",
    "organism",
    "taxonomy_id",
    "protein_name",
    "gene_name",
    "full_sequence",
    "sp_sequence",
    "sp_start",
    "sp_end",
    "sp_length",
    "sp_type",
    "cleavage_site_motif",
    "evidence_method",
    "viral_subtype",
    "viral_host",
    "query_group",
    "date_retrieved",
]
