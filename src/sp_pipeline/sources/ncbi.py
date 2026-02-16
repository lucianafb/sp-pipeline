"""NCBI/GenBank data source for SP-Pipeline.

Queries NCBI databases via Entrez API for signal peptide annotations.
Particularly useful for viral sequences with better coverage than UniProt.

Entrez docs: https://www.ncbi.nlm.nih.gov/books/NBK25497/
"""

import logging
import re
import time
from typing import Any, Optional

from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature

from sp_pipeline.models import SPRecord
from sp_pipeline.sources import BaseSource
from sp_pipeline.utils import RateLimiter, QueryCache, make_cache_key, extract_cleavage_motif

logger = logging.getLogger("sp_pipeline.sources.ncbi")

# NCBI taxonomy IDs
TAXONOMY_IDS = {
    "influenza_a": "11320",
    "influenza_b": "11520",
    "influenza_c": "11552",
    "influenza_d": "1513237",
    "influenza_all": "11308",    # Orthomyxoviridae can be too broad; use specific
    "alphavirus": "11019",
    "flavivirus": "11050",
    "human": "9606",
}

BATCH_SIZE = 200  # Entrez fetch batch size


class NCBISource(BaseSource):
    """NCBI/GenBank data source for signal peptide queries."""

    def __init__(
        self,
        email: str,
        api_key: Optional[str] = None,
        cache: Optional[QueryCache] = None,
    ):
        if not email:
            raise ValueError(
                "NCBI Entrez requires an email. Set it in config.yaml or "
                "SP_PIPELINE_NCBI_EMAIL environment variable."
            )
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        self._cache = cache
        rate = 10 if api_key else 3
        self._rate_limiter = RateLimiter(calls_per_second=rate)

    @property
    def name(self) -> str:
        return "NCBI"

    def is_available(self) -> bool:
        """Check if NCBI Entrez is reachable."""
        try:
            self._rate_limiter.wait()
            handle = Entrez.einfo()
            handle.read()
            handle.close()
            return True
        except Exception:
            return False

    def query(self, params: dict[str, Any]) -> list[SPRecord]:
        """Query NCBI for signal peptide records.

        Args:
            params: Query parameters.

        Returns:
            List of SPRecord instances.
        """
        query_group = params.get("query_group", "custom")
        search_query = self._build_query(params)
        logger.info(f"NCBI search query: {search_query}")

        # Check cache
        cache_key = make_cache_key(source="ncbi", query=search_query)
        if self._cache:
            cached = self._cache.get(cache_key)
            if cached is not None:
                logger.info(f"Cache hit: {len(cached)} records")
                return [SPRecord(**r) for r in cached]

        # Search for IDs
        ids = self._search(search_query, params)
        logger.info(f"Found {len(ids)} NCBI IDs")

        if not ids:
            return []

        # Fetch records in batches
        records = []
        for i in range(0, len(ids), BATCH_SIZE):
            batch_ids = ids[i : i + BATCH_SIZE]
            batch_records = self._fetch_batch(batch_ids, params)
            records.extend(batch_records)
            logger.debug(f"Fetched batch {i // BATCH_SIZE + 1}: {len(batch_records)} records")

        logger.info(f"Parsed {len(records)} SP records from NCBI")

        # Set query group
        for r in records:
            r.query_group = query_group

        # Cache results
        if self._cache:
            self._cache.set(cache_key, [r.to_dict() for r in records])

        return records

    def _build_query(self, params: dict[str, Any]) -> str:
        """Build an NCBI Entrez search query.

        Args:
            params: Query parameters.

        Returns:
            Entrez query string.
        """
        parts = []

        # Organism / taxonomy filter
        organism_group = params.get("organism_group", "")
        if organism_group:
            if organism_group == "influenza":
                # Search across all influenza types
                parts.append(
                    "(txid11320[Organism] OR txid11520[Organism] OR "
                    "txid11552[Organism] OR txid1513237[Organism])"
                )
            elif organism_group in TAXONOMY_IDS:
                tax_id = TAXONOMY_IDS[organism_group]
                parts.append(f"txid{tax_id}[Organism]")
        elif params.get("taxonomy_id"):
            parts.append(f"txid{params['taxonomy_id']}[Organism]")
        elif params.get("organism"):
            parts.append(f'"{params["organism"]}"[Organism]')

        # Protein name filter
        protein = params.get("protein", "")
        if protein:
            protein_map = {
                "hemagglutinin": "hemagglutinin[Protein Name]",
                "neuraminidase": "neuraminidase[Protein Name]",
                "E3_pE2": "(E3[Protein Name] OR p62[Protein Name] OR pE2[Protein Name])",
            }
            if protein in protein_map:
                parts.append(protein_map[protein])
            elif protein != "all":
                parts.append(f"{protein}[Protein Name]")

        # Signal peptide annotation
        parts.append("signal peptide[Feature key]")

        # Limit to protein database
        # (handled by database selection in esearch)

        query = " AND ".join(parts) if parts else "signal peptide[Feature key]"
        return query

    def _search(self, query: str, params: dict) -> list[str]:
        """Search NCBI and return a list of IDs.

        Args:
            query: Entrez search query.
            params: Original parameters (for max_results).

        Returns:
            List of NCBI IDs.
        """
        max_results = params.get("max_results", 5000)

        try:
            self._rate_limiter.wait()
            handle = Entrez.esearch(
                db="protein",
                term=query,
                retmax=max_results,
                usehistory="y",
            )
            results = Entrez.read(handle)
            handle.close()

            count = int(results.get("Count", 0))
            ids = results.get("IdList", [])
            logger.info(f"NCBI search returned {count} total hits, retrieved {len(ids)} IDs")

            return ids

        except Exception as e:
            logger.error(f"NCBI search error: {e}")
            return []

    def _fetch_batch(self, ids: list[str], params: dict) -> list[SPRecord]:
        """Fetch and parse a batch of NCBI records.

        Args:
            ids: List of NCBI protein IDs.
            params: Query parameters.

        Returns:
            List of SPRecord instances.
        """
        records = []

        try:
            self._rate_limiter.wait()
            handle = Entrez.efetch(
                db="protein",
                id=",".join(ids),
                rettype="gb",
                retmode="text",
            )

            for record in SeqIO.parse(handle, "genbank"):
                parsed = self._parse_genbank_record(record, params)
                records.extend(parsed)

            handle.close()

        except Exception as e:
            logger.error(f"NCBI fetch error for batch: {e}")

        return records

    def _parse_genbank_record(self, record, params: dict) -> list[SPRecord]:
        """Parse a BioPython SeqRecord into SPRecord(s).

        Args:
            record: BioPython SeqRecord from GenBank.
            params: Query parameters.

        Returns:
            List of SPRecord instances.
        """
        results = []
        full_sequence = str(record.seq)

        # Extract organism and taxonomy
        organism = record.annotations.get("organism", "")
        taxonomy = record.annotations.get("taxonomy", [])
        tax_id = None
        for feat in record.features:
            if feat.type == "source":
                db_xrefs = feat.qualifiers.get("db_xref", [])
                for xref in db_xrefs:
                    if xref.startswith("taxon:"):
                        tax_id = int(xref.split(":")[1])

        # Protein name
        protein_name = record.description
        gene_name = None
        for feat in record.features:
            if feat.type == "CDS" or feat.type == "Protein":
                gene_name = feat.qualifiers.get("gene", [None])[0]
                pn = feat.qualifiers.get("product", [None])[0]
                if pn:
                    protein_name = pn

        # Find signal peptide features
        for feat in record.features:
            if feat.type == "sig_peptide":
                start = int(feat.location.start) + 1  # Convert to 1-based
                end = int(feat.location.end)

                sp_seq = full_sequence[start - 1 : end]
                motif = extract_cleavage_motif(full_sequence, end)

                sp_record = SPRecord(
                    entry_id=record.id,
                    source_db="GenBank",
                    accession=record.name,
                    organism=organism,
                    taxonomy_id=tax_id,
                    protein_name=protein_name,
                    gene_name=gene_name,
                    full_sequence=full_sequence,
                    sp_sequence=sp_seq,
                    sp_start=start,
                    sp_end=end,
                    sp_type="type1_cleaved",
                    cleavage_site_motif=motif,
                    evidence_method="curated_annotation",
                    viral_subtype=self._extract_viral_subtype(organism, params),
                )
                results.append(sp_record)

        return results

    def _extract_viral_subtype(self, organism: str, params: dict) -> Optional[str]:
        """Extract viral subtype from organism name.

        Args:
            organism: Organism name string.
            params: Query parameters.

        Returns:
            Viral subtype or None.
        """
        organism_group = params.get("organism_group", "")
        if not organism_group:
            return None

        if organism_group == "influenza":
            match = re.search(r"\(H(\d+)N(\d+)\)", organism)
            if match:
                return f"H{match.group(1)}N{match.group(2)}"
            if "Influenza B" in organism:
                return "B"
            elif "Influenza C" in organism:
                return "C"
            elif "Influenza D" in organism:
                return "D"

        elif organism_group == "flavivirus":
            name = organism.lower()
            mapping = {
                "dengue": "DENV",
                "zika": "ZIKV",
                "west nile": "WNV",
                "yellow fever": "YFV",
                "japanese encephalitis": "JEV",
                "tick-borne": "TBEV",
            }
            for key, val in mapping.items():
                if key in name:
                    # Try to get serotype for dengue
                    if key == "dengue":
                        match = re.search(r"type (\d)", name)
                        return f"DENV-{match.group(1)}" if match else "DENV"
                    return val

        elif organism_group == "alphavirus":
            name = organism.lower()
            mapping = {
                "chikungunya": "CHIKV",
                "sindbis": "SINV",
                "semliki": "SFV",
                "venezuelan": "VEEV",
                "eastern": "EEEV",
                "western": "WEEV",
                "ross river": "RRV",
                "o'nyong": "ONNV",
            }
            for key, val in mapping.items():
                if key in name:
                    return val

        return None
