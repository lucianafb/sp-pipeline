"""UniProt data source for SP-Pipeline.

Queries UniProt REST API for signal peptide annotations.
Supports both SwissProt (reviewed/curated) and TrEMBL (unreviewed).

UniProt REST API docs: https://rest.uniprot.org/
"""

import logging
import re
import time
from typing import Any, Optional

import requests

from sp_pipeline.models import SPRecord
from sp_pipeline.sources import BaseSource
from sp_pipeline.utils import RateLimiter, QueryCache, make_cache_key, extract_cleavage_motif, compute_sp_features

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
        """Query UniProt for signal peptide records."""
        query_group = params.get("query_group", "custom")
        uniprot_query = self._build_query(params)
        logger.info(f"UniProt query: {uniprot_query}")

        # Check cache
        cache_key = make_cache_key(source="uniprot", query=uniprot_query)
        if self._cache:
            cached = self._cache.get(cache_key)
            if cached is not None:
                logger.info(f"Cache hit: {len(cached)} records")
                records = [SPRecord(**r) for r in cached]
                for r in records:
                    r.query_group = query_group
                return records

        # Determine which fields to request based on query type
        fields = self._get_fields(params)

        # Fetch from API
        raw_results = self._fetch_all(uniprot_query, fields)
        logger.info(f"Fetched {len(raw_results)} raw entries from UniProt")

        # Parse into SPRecords
        records = []
        for entry in raw_results:
            parsed = self._parse_entry(entry, params)
            records.extend(parsed)

        logger.info(f"Parsed {len(records)} SP records")

        # Cache results (before setting query_group so cache is group-agnostic)
        if self._cache:
            self._cache.set(cache_key, [r.to_dict() for r in records])

        # Set query group
        for r in records:
            r.query_group = query_group

        return records

    def _get_fields(self, params: dict[str, Any]) -> str:
        """Get the API fields to request based on query type."""
        base_fields = (
            "accession,id,protein_name,gene_names,organism_name,"
            "organism_id,sequence,ft_signal,ft_transmem,cc_subcellular_location"
        )
        # For viral polyprotein queries, also request chain/propeptide/site features
        organism_group = params.get("organism_group", "")
        if organism_group in ("flavivirus", "alphavirus", "bunyavirus"):
            base_fields += ",ft_chain,ft_peptide,ft_topo_dom,ft_site"
        return base_fields

    def _build_query(self, params: dict[str, Any]) -> str:
        """Build a UniProt query string from parameters."""
        parts = []

        sp_type = params.get("sp_type", "all")
        organism_group = params.get("organism_group", "")

        # --- HUMAN QUERIES ---
        if not organism_group:
            if sp_type == "type1_cleaved":
                parts.append("(ft_signal:*)")
            elif sp_type == "type2_signal_anchor":
                # Signal anchors: transmembrane near N-terminus
                # Query for proteins with TRANSMEM annotation that are type II
                # (N-terminus outside, single-pass)
                parts.append("(ft_transmem:*)")
                parts.append('(cc_subcellular_location:"Single-pass type II")')
            elif sp_type in ("all", "internal_signal"):
                parts.append("(ft_signal:*)")

            # Organism filter
            if params.get("taxonomy_id"):
                parts.append(f"(taxonomy_id:{params['taxonomy_id']})")
            elif params.get("organism"):
                parts.append(f'(organism_name:"{params["organism"]}")')

            # Evidence filter for query (reviewed = curated)
            evidence = params.get("evidence", "all")
            if evidence in ("experimental", "curated_annotation"):
                parts.append("(reviewed:true)")

        # --- INFLUENZA QUERIES ---
        elif organism_group == "influenza":
            protein = params.get("protein", "")

            if protein == "hemagglutinin":
                # HA has a classical N-terminal signal peptide (ft_signal)
                parts.append("(ft_signal:*)")
                parts.append(
                    "(taxonomy_id:11320 OR taxonomy_id:11520 OR "
                    "taxonomy_id:11552 OR taxonomy_id:1513237)"
                )
                parts.append("(protein_name:hemagglutinin)")

            elif protein == "neuraminidase":
                # NA is a type II membrane protein — it has NO cleavable signal peptide.
                # The N-terminal transmembrane anchor (aa ~7-27) serves as signal anchor.
                # Query with ft_transmem, NOT ft_signal (which returns 0 results for NA).
                # Use reviewed:true to limit to ~144 curated entries (unreviewed = 90k+).
                parts.append("(ft_transmem:*)")
                parts.append(
                    "(taxonomy_id:11320 OR taxonomy_id:11520 OR "
                    "taxonomy_id:11552 OR taxonomy_id:1513237)"
                )
                parts.append("(protein_name:neuraminidase)")
                parts.append("(reviewed:true)")

        # --- ALPHAVIRUS QUERIES ---
        elif organism_group == "alphavirus":
            # Alphavirus structural polyprotein: C - E3 - E2 - 6K - E1
            # E3 (Chain feature) is the signal peptide for E2.
            # 6K (Chain feature) is the signal peptide for E1.
            # There are NO ft_signal features on these entries in UniProt —
            # signals are encoded as Chain+Site(signal peptidase) features.
            # Query for structural polyprotein entries with Chain annotations.
            parts.append("(taxonomy_id:11019)")
            parts.append('(protein_name:"Structural polyprotein")')
            parts.append("(reviewed:true)")

        # --- BUNYAVIRUS QUERIES ---
        elif organism_group == "bunyavirus":
            # Bunyavirus GPC: NH₂-[SP]-[Gn]-[Gc]-COOH
            # The N-terminal SP is a classic cleaved Type 1 signal peptide.
            # Genus-level taxonomy IDs for the two main groups:
            #   Orthohantavirus: 1980415
            #   Orthobunyavirus: 11572
            # Specific viruses can be passed via taxonomy_id param.
            if params.get("taxonomy_id"):
                parts.append(f"(taxonomy_id:{params['taxonomy_id']})")
            else:
                parts.append("(taxonomy_id:1980415 OR taxonomy_id:11572)")
            # GPC / M-segment glycoprotein precursor entries
            parts.append(
                "(protein_name:glycoprotein OR protein_name:\"envelope protein\" "
                "OR protein_name:\"M segment\" OR protein_name:\"polyprotein M\")"
            )
            parts.append("(reviewed:true)")

        # --- FLAVIVIRUS QUERIES ---
        elif organism_group == "flavivirus":
            # Flavivirus signal peptides are internal signals within the Genome polyprotein.
            # They are NOT annotated as ft_signal — instead they appear as Propeptide
            # and Transmembrane features at chain boundaries in the polyprotein.
            # Genus taxonomy_id:11051 returns 0 results in UniProt REST API;
            # must use specific virus taxonomy IDs.
            FLAVI_TAXIDS = (
                "taxonomy_id:12637 OR "   # Dengue virus (all serotypes)
                "taxonomy_id:64320 OR "   # Zika virus
                "taxonomy_id:11082 OR "   # West Nile virus
                "taxonomy_id:11053 OR "   # Yellow fever virus
                "taxonomy_id:11080 OR "   # Japanese encephalitis virus
                "taxonomy_id:11090"       # Tick-borne encephalitis virus
            )
            parts.append(f"({FLAVI_TAXIDS})")
            parts.append('(protein_name:"Genome polyprotein")')
            parts.append("(reviewed:true)")

        query = " AND ".join(parts) if parts else "ft_signal:*"
        return query

    def _fetch_all(self, query: str, fields: Optional[str] = None) -> list[dict]:
        """Fetch all results from UniProt API with pagination, retry and progress bar."""
        from tqdm import tqdm

        results = []

        if not fields:
            fields = (
                "accession,id,protein_name,gene_names,organism_name,"
                "organism_id,sequence,ft_signal,ft_transmem,cc_subcellular_location"
            )

        url = (
            f"{UNIPROT_API_BASE}/search?"
            f"query={query}&"
            f"fields={fields}&"
            f"format=json&size={BATCH_SIZE}"
        )

        MAX_RETRIES = 3
        BACKOFF = 2  # seconds base for exponential backoff

        total = None
        pbar = None

        while url:
            self._rate_limiter.wait()

            # Retry with exponential backoff
            for attempt in range(MAX_RETRIES):
                try:
                    response = requests.get(url, timeout=(10, 60))

                    if response.status_code == 200:
                        break
                    if response.status_code == 429:
                        wait = int(response.headers.get("Retry-After", BACKOFF ** (attempt + 1)))
                        logger.warning(f"Rate limited (429). Waiting {wait}s...")
                        time.sleep(wait)
                        continue
                    if response.status_code >= 500:
                        wait = BACKOFF ** (attempt + 1)
                        logger.warning(f"Server error ({response.status_code}). Retry in {wait}s...")
                        time.sleep(wait)
                        continue
                    response.raise_for_status()

                except requests.exceptions.Timeout:
                    time.sleep(BACKOFF ** (attempt + 1))
                except requests.exceptions.ConnectionError as e:
                    logger.warning(f"Connection error: {e}. Retry...")
                    time.sleep(BACKOFF ** (attempt + 1))
            else:
                logger.error(f"UniProt API failed after {MAX_RETRIES} attempts")
                break

            data = response.json()
            batch = data.get("results", [])
            results.extend(batch)

            # Init progress bar on first response
            if total is None:
                total = int(response.headers.get("x-total-results", 0))
                if total > 0:
                    pbar = tqdm(total=total, desc="Downloading", unit=" entries")

            if pbar:
                pbar.update(len(batch))

            logger.debug(f"Fetched batch: {len(batch)} entries")

            # Pagination
            url = None
            link_header = response.headers.get("Link", "")
            if 'rel="next"' in link_header:
                url = link_header.split(";")[0].strip("<>")

        if pbar:
            pbar.close()

        return results

    def _parse_entry(self, entry: dict, params: dict) -> list[SPRecord]:
        """Parse a UniProt JSON entry into SPRecord(s)."""
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

        # Parse features
        features = entry.get("features", [])
        organism_group = params.get("organism_group", "")
        sp_type_filter = params.get("sp_type", "all")

        for feat in features:
            feat_type = feat.get("type", "")

            # --- Signal peptide (type 1, cleaved) ---
            if feat_type == "Signal":
                location = feat.get("location", {})
                start = location.get("start", {}).get("value", 0)
                end = location.get("end", {}).get("value", 0)

                if start and end and full_sequence:
                    sp_seq = full_sequence[start - 1 : end]
                    motif = extract_cleavage_motif(full_sequence, end)
                    evidences = feat.get("evidences", [])
                    evidence_method = self._determine_evidence(evidences)
                    bio = compute_sp_features(sp_seq)

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
                        **bio,
                    )
                    records.append(record)

            # --- Transmembrane (type 2, signal anchor or NA-like) ---
            elif feat_type == "Transmembrane":
                # Viral polyproteins are handled by specialized parsers —
                # skip generic TM processing to avoid creating duplicate records.
                if organism_group in ("flavivirus", "alphavirus", "bunyavirus"):
                    continue

                # For human type 2 query OR influenza NA
                if sp_type_filter in ("type2_signal_anchor", "all") or (
                    organism_group == "influenza"
                    and params.get("protein") == "neuraminidase"
                ):
                    location = feat.get("location", {})
                    start = location.get("start", {}).get("value", 0)
                    end = location.get("end", {}).get("value", 0)

                    # For signal anchors: transmembrane must be near N-terminus
                    if start and start <= 60 and end and full_sequence:
                        sp_seq = full_sequence[0:end]
                        evidence_method = self._determine_evidence(
                            feat.get("evidences", [])
                        )
                        bio = compute_sp_features(sp_seq)

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
                            **bio,
                        )
                        records.append(record)

        # --- Special handling for flavivirus polyproteins ---
        if organism_group == "flavivirus" and "polyprotein" in protein_name.lower():
            polyprotein_records = self._parse_flavivirus_polyprotein(
                entry, full_sequence, params
            )
            records.extend(polyprotein_records)

        # --- Special handling for alphavirus polyproteins ---
        if organism_group == "alphavirus" and (
            "polyprotein" in protein_name.lower()
            or "structural" in protein_name.lower()
        ):
            polyprotein_records = self._parse_alphavirus_polyprotein(
                entry, full_sequence, params
            )
            records.extend(polyprotein_records)

        # --- Special handling for bunyavirus GPC ---
        if organism_group == "bunyavirus":
            gpc_records = self._parse_bunyavirus_gpc(entry, full_sequence, params)
            records.extend(gpc_records)

        return records

    def _parse_flavivirus_polyprotein(
        self, entry: dict, full_sequence: str, params: dict
    ) -> list[SPRecord]:
        """Parse signal peptides from flavivirus genome polyprotein.

        Flavivirus polyproteins have NO ft_signal features in UniProt.
        Signal peptides are identified via Site features annotated as
        "Cleavage; by host signal peptidase":
          - First site (e.g. DENV aa 14-15): N-terminal SP = residues 1 to site.start-1
            (the ER anchor / capsid anchor that functions as SP for prM)
          - Subsequent sites (e.g. DENV aa 180-181, 675-676): internal SP = the
            transmembrane region that ends immediately before the cleavage site
            (TM of M acts as SP for E; TM of E acts as SP for NS1)
        """
        records = []
        entry_id = entry.get("primaryAccession", "")
        organism = entry.get("organism", {}).get("scientificName", "")
        tax_id = entry.get("organism", {}).get("taxonId")
        protein_name = "Genome polyprotein"
        viral_subtype = self._extract_viral_subtype(entry, params)
        features = entry.get("features", [])

        # Collect signal peptidase cleavage sites
        sp_sites = []
        for feat in features:
            if feat.get("type") != "Site":
                continue
            desc = feat.get("description", "").lower()
            if "signal peptidase" not in desc:
                continue
            loc = feat.get("location", {})
            cut = loc.get("start", {}).get("value", 0)  # position AFTER cleavage
            if cut:
                sp_sites.append({"cut": cut, "evidences": feat.get("evidences", [])})

        if not sp_sites:
            return records

        sp_sites.sort(key=lambda s: s["cut"])

        # Collect TM regions for internal signal lookup
        tm_regions = []
        for feat in features:
            if feat.get("type") != "Transmembrane":
                continue
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value", 0)
            end = loc.get("end", {}).get("value", 0)
            if start and end:
                tm_regions.append(
                    {"start": start, "end": end, "evidences": feat.get("evidences", [])}
                )

        # Process each signal peptidase site
        for i, site in enumerate(sp_sites):
            cut = site["cut"]
            if not full_sequence:
                continue

            if i == 0:
                # First site: N-terminal SP (ER anchor of capsid → cleaved to release prM)
                # Biologically, this anchor is always short (< 30 aa in flaviviruses).
                # If the first annotated site is far from the N-terminus, skip — it means
                # the C→prM anchor site is not annotated in this entry.
                sp_start = 1
                sp_end = cut - 1
                if sp_end < sp_start or sp_end > 30:
                    continue
                sp_seq = full_sequence[sp_start - 1 : sp_end]
                motif = extract_cleavage_motif(full_sequence, sp_end)
                bio = compute_sp_features(sp_seq)
                records.append(SPRecord(
                    entry_id=entry_id,
                    source_db="UniProt",
                    organism=organism,
                    taxonomy_id=tax_id,
                    protein_name=protein_name,
                    gene_name="",
                    full_sequence=full_sequence,
                    sp_sequence=sp_seq,
                    sp_start=sp_start,
                    sp_end=sp_end,
                    sp_type="type1_cleaved",
                    cleavage_site_motif=motif,
                    evidence_method=self._determine_evidence(site["evidences"]),
                    viral_subtype=viral_subtype,
                    **bio,
                ))
            else:
                # Subsequent sites: find the TM region that ends just before this cut
                # (the TM of the previous protein acts as signal for the next)
                matching_tm = None
                for tm in tm_regions:
                    if abs(tm["end"] - (cut - 1)) <= 5:
                        matching_tm = tm
                        break

                if matching_tm:
                    sp_seq = full_sequence[matching_tm["start"] - 1 : matching_tm["end"]]
                    bio = compute_sp_features(sp_seq)
                    records.append(SPRecord(
                        entry_id=entry_id,
                        source_db="UniProt",
                        organism=organism,
                        taxonomy_id=tax_id,
                        protein_name=protein_name,
                        gene_name="",
                        full_sequence=full_sequence,
                        sp_sequence=sp_seq,
                        sp_start=matching_tm["start"],
                        sp_end=matching_tm["end"],
                        sp_type="internal_signal",
                        cleavage_site_motif=None,
                        evidence_method=self._determine_evidence(
                            matching_tm["evidences"]
                        ),
                        viral_subtype=viral_subtype,
                        **bio,
                    ))

        return records

    def _parse_alphavirus_polyprotein(
        self, entry: dict, full_sequence: str, params: dict
    ) -> list[SPRecord]:
        """Parse signal peptides from alphavirus structural polyprotein.

        Structure: C - E3 - E2 - 6K - E1
        E3 (Chain feature "Assembly protein E3") is the signal peptide for E2.
        6K (Chain feature "6K protein") is the signal peptide for E1.

        In UniProt, these are annotated as Chain features with Site(signal peptidase)
        cleavage sites — NOT as ft_signal features.
        """
        records = []
        entry_id = entry.get("primaryAccession", "")
        organism = entry.get("organism", {}).get("scientificName", "")
        tax_id = entry.get("organism", {}).get("taxonId")
        protein_name = (
            entry.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "Structural polyprotein")
        )
        viral_subtype = self._extract_viral_subtype(entry, params)
        features = entry.get("features", [])

        # Target chains: E3 and 6K both serve as signal peptides (cleaved by signal peptidase).
        # Exclude "Precursor of protein E3/E2" (which spans both E3+E2).
        SP_CHAIN_KEYWORDS = ("assembly protein e3", "6k protein")

        for feat in features:
            if feat.get("type") != "Chain":
                continue
            desc = feat.get("description", "").lower()
            if "precursor" in desc:
                continue
            if not any(kw in desc for kw in SP_CHAIN_KEYWORDS):
                continue

            location = feat.get("location", {})
            start = location.get("start", {}).get("value", 0)
            end = location.get("end", {}).get("value", 0)

            if start and end and full_sequence:
                sp_seq = full_sequence[start - 1 : end]
                motif = extract_cleavage_motif(full_sequence, end)
                bio = compute_sp_features(sp_seq)
                record = SPRecord(
                    entry_id=entry_id,
                    source_db="UniProt",
                    organism=organism,
                    taxonomy_id=tax_id,
                    protein_name=f"{protein_name} [{feat.get('description', '')}]",
                    gene_name="",
                    full_sequence=full_sequence,
                    sp_sequence=sp_seq,
                    sp_start=start,
                    sp_end=end,
                    sp_type="type1_cleaved",
                    cleavage_site_motif=motif,
                    evidence_method=self._determine_evidence(
                        feat.get("evidences", [])
                    ),
                    viral_subtype=viral_subtype,
                    **bio,
                )
                records.append(record)

        return records

    def _parse_bunyavirus_gpc(
        self, entry: dict, full_sequence: str, params: dict
    ) -> list[SPRecord]:
        """Parse signal peptides from bunyavirus GPC (Glycoprotein Precursor).

        Architecture: NH₂-[SP]-[Gn]-[Gc]-COOH
        The N-terminal SP is a classic Type 1 cleaved signal peptide.

        Extraction strategies (in order):
        1. Explicit Signal feature annotation
        2. Infer SP from Chain boundaries: if Chain "Glycoprotein N" starts at pos > 1,
           SP = residues 1..(Gn_start - 1)
        3. Fallback: first Transmembrane region near N-terminus as putative SP
        """
        records = []
        entry_id = entry.get("primaryAccession", "")
        organism = entry.get("organism", {}).get("scientificName", "")
        tax_id = entry.get("organism", {}).get("taxonId")
        protein_name = (
            entry.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "Glycoprotein precursor")
        )
        viral_subtype = self._extract_viral_subtype(entry, params)
        features = entry.get("features", [])

        def _build_record(sp_seq, start, end, sp_type, ev_list):
            motif = extract_cleavage_motif(full_sequence, end)
            bio = compute_sp_features(sp_seq)
            return SPRecord(
                entry_id=entry_id,
                source_db="UniProt",
                organism=organism,
                taxonomy_id=tax_id,
                protein_name=protein_name,
                gene_name="",
                full_sequence=full_sequence,
                sp_sequence=sp_seq,
                sp_start=start,
                sp_end=end,
                sp_type=sp_type,
                cleavage_site_motif=motif,
                evidence_method=self._determine_evidence(ev_list),
                viral_subtype=viral_subtype,
                **bio,
            )

        # Collect features by type
        signals = []
        chains = []
        tms = []
        for feat in features:
            ft_type = feat.get("type", "")
            loc = feat.get("location", {})
            s = loc.get("start", {}).get("value")
            e = loc.get("end", {}).get("value")
            if not (isinstance(s, int) and isinstance(e, int)):
                continue
            item = {"start": s, "end": e, "desc": feat.get("description", ""), "evidences": feat.get("evidences", [])}
            if ft_type == "Signal":
                signals.append(item)
            elif ft_type == "Chain":
                chains.append(item)
            elif ft_type == "Transmembrane":
                tms.append(item)

        chains.sort(key=lambda x: x["start"])
        tms.sort(key=lambda x: x["start"])

        # Strategy 1: Explicit Signal feature
        for sig in signals:
            sp_seq = full_sequence[sig["start"] - 1 : sig["end"]]
            if sp_seq and len(sp_seq) >= 5:
                records.append(_build_record(sp_seq, sig["start"], sig["end"], "type1_cleaved", sig["evidences"]))

        if records:
            return records

        # Strategy 2: Infer SP from Gn chain start
        gn_chain = None
        for ch in chains:
            desc = ch["desc"].lower()
            if any(p in desc for p in ("glycoprotein n", "glycoprotein gn", "protein gn")):
                gn_chain = ch
                break

        if gn_chain and gn_chain["start"] > 1:
            sp_end = gn_chain["start"] - 1
            sp_seq = full_sequence[0:sp_end]
            if sp_seq and 5 <= len(sp_seq) <= 50:
                records.append(_build_record(sp_seq, 1, sp_end, "type1_cleaved", gn_chain["evidences"]))

        if records:
            return records

        # Strategy 3: Fallback — first TM near N-terminus
        for tm in tms:
            if tm["start"] < 80:
                sp_end = tm["start"] - 1
                if sp_end >= 5:
                    sp_seq = full_sequence[0:sp_end]
                    if sp_seq and 5 <= len(sp_seq) <= 50:
                        records.append(_build_record(sp_seq, 1, sp_end, "type1_cleaved", tm["evidences"]))
                break  # Only use the first TM

        return records

    def _determine_evidence(self, evidences: list[dict]) -> str:
        """Determine the evidence level from UniProt evidence codes."""
        if not evidences:
            return "curated_annotation"

        for ev in evidences:
            code = ev.get("evidenceCode", "")
            if code == "ECO:0000269":
                return "experimental"
            elif code == "ECO:0000305":
                return "curated_annotation"
            elif code == "ECO:0000256":
                return "predicted_automatic"

        return "curated_annotation"

    def _extract_viral_subtype(self, entry: dict, params: dict) -> Optional[str]:
        """Extract viral subtype from organism name."""
        organism_group = params.get("organism_group", "")
        if not organism_group:
            return None

        organism = entry.get("organism", {}).get("scientificName", "")

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
            if "dengue" in name:
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
            elif "st. louis" in name or "saint louis" in name:
                return "SLEV"
            elif "murray valley" in name:
                return "MVEV"
            elif "powassan" in name:
                return "POWV"

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
                "mayaro": "MAYV",
                "una": "UNAV",
                "barmah": "BFV",
            }
            for key, val in mapping.items():
                if key in name:
                    return val

        elif organism_group == "bunyavirus":
            name = organism.lower()
            # Orthohantaviruses
            hanta_mapping = {
                "andes": "ANDV",
                "sin nombre": "SNV",
                "hantaan": "HTNV",
                "seoul": "SEOV",
                "puumala": "PUUV",
                "dobrava": "DOBV",
                "bayou": "BAYV",
                "black creek canal": "BCCV",
                "new york": "NYV",
            }
            # Orthobunyaviruses
            bunya_mapping = {
                "la crosse": "LACV",
                "bunyamwera": "BUNV",
                "schmallenberg": "SBV",
                "california encephalitis": "CEV",
                "tahyna": "TAHV",
                "oropouche": "OROV",
            }
            for key, val in {**hanta_mapping, **bunya_mapping}.items():
                if key in name:
                    return val

        return None


