"""UniProt data source for SP-Pipeline.

Queries UniProt REST API for signal peptide annotations.
Supports both SwissProt (reviewed/curated) and TrEMBL (unreviewed).

UniProt REST API docs: https://rest.uniprot.org/
"""

import logging
import time
from typing import Any, Optional

import requests

from sp_pipeline.models import SPRecord
from sp_pipeline.sources import BaseSource
from sp_pipeline.utils import RateLimiter, QueryCache, make_cache_key, extract_cleavage_motif

logger = logging.getLogger("sp_pipeline.sources.uniprot")

UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"
BATCH_SIZE = 500  # UniProt API pagination size


class UniProtSource(BaseSource):
    """UniProt/SwissProt data source for signal peptide queries."""

    def __init__(self, cache: Optional[QueryCache] = None):
        self._cache = cache
        self._rate_limiter = RateLimiter(calls_per_second=5)

    @property
    def name(self) -> str:
        return "UniProt"

    def is_available(self) -> bool:
        """Check if UniProt API is reachable."""
        try:
            r = requests.get(f"{UNIPROT_API_BASE}/search?query=*&size=1", timeout=10)
            return r.status_code == 200
        except requests.RequestException:
            return False

    def query(self, params: dict[str, Any]) -> list[SPRecord]:
        """Query UniProt for signal peptide records.

        Args:
            params: Query parameters including:
                - organism: str (e.g., "Homo sapiens")
                - taxonomy_id: int (e.g., 9606)
                - sp_type: str (type1_cleaved, type2_signal_anchor, all)
                - evidence: str (experimental, predicted, curated_annotation, all)
                - organism_group: str (for viral queries: influenza, alphavirus, flavivirus)
                - protein: str (protein name filter)
                - query_group: str (preset name)

        Returns:
            List of SPRecord instances.
        """
        query_group = params.get("query_group", "custom")
        uniprot_query = self._build_query(params)
        logger.info(f"UniProt query: {uniprot_query}")

        # Check cache
        cache_key = make_cache_key(source="uniprot", query=uniprot_query)
        if self._cache:
            cached = self._cache.get(cache_key)
            if cached is not None:
                logger.info(f"Cache hit: {len(cached)} records")
                return [SPRecord(**r) for r in cached]

        # Fetch from API
        raw_results = self._fetch_all(uniprot_query)
        logger.info(f"Fetched {len(raw_results)} raw entries from UniProt")

        # Parse into SPRecords
        records = []
        for entry in raw_results:
            parsed = self._parse_entry(entry, params)
            records.extend(parsed)

        logger.info(f"Parsed {len(records)} SP records")

        # Cache results
        if self._cache:
            self._cache.set(cache_key, [r.to_dict() for r in records])

        # Set query group
        # Set query group
        for r in records:
            r.query_group = query_group

        # Filter by evidence
        evidence = params.get("evidence", "all")
        records = _filter_by_evidence(records, evidence)
        logger.info(f"After evidence filter ({evidence}): {len(records)} records")

        return records

    def _build_query(self, params: dict[str, Any]) -> str:
        """Build a UniProt query string from parameters.

        UniProt query syntax:
        https://www.uniprot.org/help/query-fields
        """
        parts = []

        # Signal peptide annotation must exist
        sp_type = params.get("sp_type", "all")
        if sp_type == "type1_cleaved":
            parts.append("(ft_signal:*)")
        elif sp_type == "type2_signal_anchor":
            parts.append("(ft_transit:*)")  # TRANSIT or use annotation
            # Signal anchors are annotated differently - use TRANSMEM near N-terminus
            # or SIGNAL with "Signal-anchor" in description
            parts.append("(cc_subcellular_location:\"Signal-anchor\")")
            # Remove the transit part, replace with combined approach
            parts.clear()
            parts.append(
                "((ft_signal:*) AND (cc_subcellular_location:\"signal-anchor\"))"
                " OR (ft_transmem:* AND annotation_score:5)"
            )
        elif sp_type in ("all", "internal_signal"):
            parts.append("(ft_signal:*)")

        # Organism filter
        if "taxonomy_id" in params and params["taxonomy_id"]:
            parts.append(f"(taxonomy_id:{params['taxonomy_id']})")
        elif "organism" in params and params["organism"]:
            parts.append(f"(organism_name:\"{params['organism']}\")")

        # Organism group (for viral queries)
        organism_group = params.get("organism_group", "")
        if organism_group == "influenza":
            parts.append("(taxonomy_id:11308)")  # Influenza A, B, C, D
            protein = params.get("protein", "")
            if protein == "hemagglutinin":
                parts.append("(protein_name:hemagglutinin)")
            elif protein == "neuraminidase":
                parts.append("(protein_name:neuraminidase)")
        elif organism_group == "alphavirus":
            parts.append("(taxonomy_id:11019)")  # Alphavirus genus
            protein = params.get("protein", "")
            if protein in ("E3_pE2", "E3", "pE2"):
                parts.append(
                    "((protein_name:\"Protein E3\") OR (protein_name:\"p62\") "
                    "OR (protein_name:\"PE2\"))"
                )
        elif organism_group == "flavivirus":
            parts.append("(taxonomy_id:11050)")  # Flavivirus genus

        # Evidence filter
        evidence = params.get("evidence", "all")
        if evidence == "experimental":
            parts.append("(reviewed:true)")  # SwissProt only as proxy for experimental
        elif evidence == "curated_annotation":
            parts.append("(reviewed:true)")

        query = " AND ".join(parts) if parts else "ft_signal:*"
        return query

    def _fetch_all(self, query: str) -> list[dict]:
        """Fetch all results from UniProt API with pagination.

        Args:
            query: UniProt query string.

        Returns:
            List of raw entry dictionaries.
        """
        results = []
        url = (
            f"{UNIPROT_API_BASE}/search?"
            f"query={query}&"
            f"fields=accession,id,protein_name,gene_names,organism_name,"
            f"organism_id,sequence,ft_signal,ft_transmem,cc_subcellular_location&"
            f"format=json&size={BATCH_SIZE}"
        )

        while url:
            self._rate_limiter.wait()
            try:
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                data = response.json()
                results.extend(data.get("results", []))
                logger.debug(f"Fetched batch: {len(data.get('results', []))} entries")

                # Pagination: check for 'Link' header
                url = None
                link_header = response.headers.get("Link", "")
                if 'rel="next"' in link_header:
                    # Extract next URL
                    url = link_header.split(";")[0].strip("<>")

            except requests.RequestException as e:
                logger.error(f"UniProt API error: {e}")
                break

        return results

    def _parse_entry(self, entry: dict, params: dict) -> list[SPRecord]:
        """Parse a UniProt JSON entry into SPRecord(s).

        A single entry can produce multiple records if it has
        multiple signal peptides or signal anchors.

        Args:
            entry: Raw UniProt JSON entry.
            params: Original query parameters.

        Returns:
            List of SPRecord instances.
        """
        records = []

        # Basic info
        entry_id = entry.get("primaryAccession", "")
        organism = entry.get("organism", {}).get("scientificName", "")
        tax_id = entry.get("organism", {}).get("taxonId")

        # Protein name
        protein_desc = entry.get("proteinDescription", {})
        rec_name = protein_desc.get("recommendedName", {})
        protein_name = rec_name.get("fullName", {}).get("value", "")
        if not protein_name:
            sub_names = protein_desc.get("submissionNames", [])
            if sub_names:
                protein_name = sub_names[0].get("fullName", {}).get("value", "")

        # Gene name
        gene_names = entry.get("genes", [])
        gene_name = gene_names[0].get("geneName", {}).get("value", "") if gene_names else ""

        # Sequence
        seq_data = entry.get("sequence", {})
        full_sequence = seq_data.get("value", "")

        # Parse features for signal peptides
        features = entry.get("features", [])
        for feat in features:
            feat_type = feat.get("type", "")

            if feat_type == "Signal":
                # Type 1 cleaved signal peptide
                location = feat.get("location", {})
                start = location.get("start", {}).get("value", 0)
                end = location.get("end", {}).get("value", 0)

                if start and end and full_sequence:
                    sp_seq = full_sequence[start - 1 : end]
                    motif = extract_cleavage_motif(full_sequence, end)

                    # Determine evidence
                    evidences = feat.get("evidences", [])
                    evidence_method = self._determine_evidence(evidences)

                    record = SPRecord(
                        entry_id=entry_id,
                        source_db="UniProt",
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
                        evidence_method=evidence_method,
                        viral_subtype=self._extract_viral_subtype(entry, params),
                    )
                    records.append(record)

            elif feat_type == "Transmembrane" and params.get("sp_type") in (
                "type2_signal_anchor",
                "all",
            ):
                # Potential type 2 signal anchor: transmembrane near N-terminus
                location = feat.get("location", {})
                start = location.get("start", {}).get("value", 0)

                if start and start <= 60:  # Near N-terminus heuristic
                    end = location.get("end", {}).get("value", 0)
                    if end and full_sequence:
                        sp_seq = full_sequence[0:end]  # Signal anchor = N-term to end of TM
                        evidence_method = self._determine_evidence(
                            feat.get("evidences", [])
                        )

                        record = SPRecord(
                            entry_id=entry_id,
                            source_db="UniProt",
                            organism=organism,
                            taxonomy_id=tax_id,
                            protein_name=protein_name,
                            gene_name=gene_name,
                            full_sequence=full_sequence,
                            sp_sequence=sp_seq,
                            sp_start=1,
                            sp_end=end,
                            sp_type="type2_signal_anchor",
                            cleavage_site_motif=None,
                            evidence_method=evidence_method,
                            viral_subtype=self._extract_viral_subtype(entry, params),
                        )
                        records.append(record)

        return records

    def _determine_evidence(self, evidences: list[dict]) -> str:
        """Determine the evidence level from UniProt evidence codes.

        Args:
            evidences: List of evidence dictionaries from UniProt.

        Returns:
            Evidence method string.
        """
        if not evidences:
            return "curated_annotation"

        for ev in evidences:
            code = ev.get("evidenceCode", "")
            # ECO:0000269 = experimental evidence
            if code == "ECO:0000269":
                return "experimental"
            # ECO:0000305 = curator inference
            elif code == "ECO:0000305":
                return "curated_annotation"
            # ECO:0000256 = sequence analysis (automatic)
            elif code == "ECO:0000256":
                return "predicted_automatic"

        return "curated_annotation"

    def _extract_viral_subtype(self, entry: dict, params: dict) -> Optional[str]:
        """Extract viral subtype from organism name or features.

        Args:
            entry: Raw UniProt entry.
            params: Query parameters.

        Returns:
            Viral subtype string or None.
        """
        organism_group = params.get("organism_group", "")
        if not organism_group:
            return None

        organism = entry.get("organism", {}).get("scientificName", "")

        if organism_group == "influenza":
            # Try to extract H/N subtype from organism name
            # e.g., "Influenza A virus (A/Puerto Rico/8/1934(H1N1))"
            import re

            match = re.search(r"\(H(\d+)N(\d+)\)", organism)
            if match:
                return f"H{match.group(1)}N{match.group(2)}"
            # Check for type B, C, D
            if "Influenza B" in organism:
                return "B"
            elif "Influenza C" in organism:
                return "C"
            elif "Influenza D" in organism:
                return "D"

        elif organism_group == "flavivirus":
            # Extract virus name abbreviation
            name = organism.lower()
            if "dengue" in name:
                import re
                match = re.search(r"type (\d)", name)
                return f"DENV-{match.group(1)}" if match else "DENV"
            elif "zika" in name:
                return "ZIKV"
            elif "west nile" in name:
                return "WNV"
            elif "yellow fever" in name:
                return "YFV"
            elif "japanese encephalitis" in name:
                return "JEV"
            elif "tick-borne" in name:
                return "TBEV"

        elif organism_group == "alphavirus":
            name = organism.lower()
            if "chikungunya" in name:
                return "CHIKV"
            elif "sindbis" in name:
                return "SINV"
            elif "semliki" in name:
                return "SFV"
            elif "venezuelan" in name:
                return "VEEV"
            elif "eastern" in name:
                return "EEEV"
            elif "western" in name:
                return "WEEV"
            elif "ross river" in name:
                return "RRV"
            elif "o'nyong" in name:
                return "ONNV"

        return None


def _filter_by_evidence(records, evidence):
    """Post-query filter by evidence method."""
    if evidence == "all":
        return records
    elif evidence == "experimental":
        return [r for r in records if r.evidence_method == "experimental"]
    elif evidence == "predicted":
        return [r for r in records if "predicted" in r.evidence_method]
    elif evidence == "curated_annotation":
        return [r for r in records if r.evidence_method == "curated_annotation"]
    return records
