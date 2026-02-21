"""Main pipeline orchestrator for SP-Pipeline.

Coordinates data sources, deduplication, and export.
"""

import logging
from typing import Any, Optional

from sp_pipeline.config import load_config, get_preset, list_presets
from sp_pipeline.dedup import deduplicate_records, merge_with_existing
from sp_pipeline.exporters import export_csv, export_fasta, generate_summary
from sp_pipeline.models import SPRecord
from sp_pipeline.predictors.signalp import SignalPPredictor
from sp_pipeline.sources.uniprot import UniProtSource
from sp_pipeline.sources.ncbi import NCBISource
from sp_pipeline.utils import setup_logging, QueryCache

logger = logging.getLogger("sp_pipeline")


class SPPipeline:
    """Main pipeline class that orchestrates the full query workflow."""

    def __init__(self, config_path: Optional[str] = None):
        """Initialize the pipeline.

        Args:
            config_path: Optional path to custom config file.
        """
        self.config = load_config(config_path)

        # Setup logging
        log_config = self.config.get("logging", {})
        setup_logging(
            level=log_config.get("level", "INFO"),
            log_file=log_config.get("file") or None,
        )

        # Setup cache
        cache_config = self.config.get("cache", {})
        self.cache = QueryCache(
            cache_dir=cache_config.get("directory", "/tmp/sp_pipeline_cache"),
            ttl_days=cache_config.get("ttl_days", 30),
            enabled=cache_config.get("enabled", True),
        )

        # Initialize sources
        self._sources = {}
        self._init_sources()

        # Initialize predictors
        self._signalp: SignalPPredictor | None = None
        self._init_predictors()

    def _init_sources(self):
        """Initialize data source modules."""
        # UniProt (always available)
        self._sources["uniprot"] = UniProtSource(cache=self.cache)

        # NCBI (requires email)
        ncbi_config = self.config.get("ncbi", {})
        email = ncbi_config.get("email", "")
        if email:
            self._sources["ncbi"] = NCBISource(
                email=email,
                api_key=ncbi_config.get("api_key") or None,
                cache=self.cache,
            )
        else:
            logger.warning(
                "NCBI source not configured. Set ncbi.email in config or "
                "SP_PIPELINE_NCBI_EMAIL env var to enable NCBI queries."
            )

    def _init_predictors(self):
        """Initialize prediction modules (SignalP, etc.)."""
        predictor_cfg = self.config.get("predictors", {}).get("signalp", {})
        token = predictor_cfg.get("api_token", "")
        if token and predictor_cfg.get("enabled", False):
            self._signalp = SignalPPredictor(
                api_token=token,
                organism=predictor_cfg.get("organism", "eukarya"),
                batch_size=predictor_cfg.get("batch_size", 500),
            )
            logger.info("SignalP predictor initialized (organism=%s)", predictor_cfg.get("organism", "eukarya"))
        else:
            self._signalp = None

    def run(
        self,
        presets: list[str],
        evidence: Optional[str] = None,
        mode: Optional[str] = None,
        include_predictions: bool = False,
        output: str = "output/results.csv",
        append: bool = False,
        include_fasta: bool = False,
    ) -> str:
        """Run the pipeline with one or more presets.

        Args:
            presets: List of preset names to run.
            evidence: Override evidence filter (experimental/predicted/all).
            mode: Override mode (exhaustive/representative).
            include_predictions: Run SignalP/Phobius de novo.
            output: Output CSV path.
            append: Append to existing CSV.
            include_fasta: Also export FASTA file.

        Returns:
            Path to the output CSV file.
        """
        all_records: list[SPRecord] = []

        for preset_name in presets:
            logger.info(f"Running preset: {preset_name}")
            try:
                preset_config = get_preset(self.config, preset_name)
            except ValueError as e:
                logger.error(str(e))
                continue

            # Apply overrides
            params = preset_config.copy()
            params["query_group"] = preset_name
            if evidence:
                params["evidence"] = evidence
            if mode:
                params["mode"] = mode

            # Query each configured source
            sources_to_use = params.get("sources", ["uniprot"])
            for source_name in sources_to_use:
                source = self._sources.get(source_name)
                if source is None:
                    logger.warning(f"Source '{source_name}' not available, skipping")
                    continue

                if not source.is_available():
                    logger.warning(f"Source '{source_name}' is not reachable, skipping")
                    continue

                logger.info(f"Querying {source.name}...")
                try:
                    records = source.query(params)
                    all_records.extend(records)
                    logger.info(f"Got {len(records)} records from {source.name}")
                except Exception as e:
                    logger.error(f"Error querying {source.name}: {e}")

        if not all_records:
            logger.warning("No records found for any preset/source combination")
            return ""

        # Deduplicate (before predictions — gives clean input for SignalP)
        logger.info(f"Total records before dedup: {len(all_records)}")
        all_records = deduplicate_records(all_records)
        logger.info(f"Total records after dedup: {len(all_records)}")

        # Run SignalP predictions if requested
        if include_predictions and self._signalp and self._signalp.is_available():
            logger.info("Running SignalP predictions on %d sequences...", len(all_records))
            seqs_to_predict = [
                {
                    "entry_id": r.entry_id,
                    "full_sequence": r.full_sequence,
                    "organism": r.organism,
                    "query_group": r.query_group,
                    "source_db": r.source_db,
                    "protein_name": r.protein_name,
                    "gene_name": r.gene_name,
                    "taxonomy_id": r.taxonomy_id,
                }
                for r in all_records
                if r.full_sequence
            ]
            predicted = self._signalp.predict(seqs_to_predict)
            logger.info(f"SignalP predicted {len(predicted)} signal peptides")
            all_records.extend(predicted)
            all_records = deduplicate_records(all_records)
        elif include_predictions:
            logger.warning(
                "SignalP predictions requested but not configured. "
                "Set SP_PIPELINE_SIGNALP_TOKEN env var and "
                "predictors.signalp.enabled=true in config."
            )

        # Filter by evidence
        if evidence and evidence != "all":
            before = len(all_records)
            if evidence == "experimental":
                all_records = [r for r in all_records if r.evidence_method == "experimental"]
            elif evidence == "predicted":
                all_records = [r for r in all_records if "predicted" in r.evidence_method]
            elif evidence == "curated_annotation":
                all_records = [r for r in all_records if r.evidence_method == "curated_annotation"]
            logger.info(f"Evidence filter ({evidence}): {before} -> {len(all_records)} records")

        # Merge with existing if appending
        existing_csv = output if append else None
        df = merge_with_existing(all_records, existing_csv)

        # Export CSV
        output_path = export_csv(
            df,
            output,
            delimiter=self.config.get("output", {}).get("delimiter", ","),
        )

        # Export FASTA if requested
        if include_fasta or self.config.get("output", {}).get("include_fasta", False):
            fasta_path = output.rsplit(".", 1)[0] + "_sp.fasta"
            export_fasta(all_records, fasta_path, sequence_type="sp")

        # Print summary
        summary = generate_summary(df)
        logger.info("\n" + summary)

        return output_path

    def run_custom(
        self,
        organism: Optional[str] = None,
        taxonomy_id: Optional[int] = None,
        sp_type: str = "all",
        evidence: str = "all",
        sources: Optional[list[str]] = None,
        output: str = "output/custom_results.csv",
        **kwargs,
    ) -> str:
        """Run a custom query without using presets.

        Args:
            organism: Organism name filter.
            taxonomy_id: NCBI taxonomy ID filter.
            sp_type: SP type filter.
            evidence: Evidence filter.
            sources: List of sources to query.
            output: Output CSV path.
            **kwargs: Additional query parameters.

        Returns:
            Path to the output CSV file.
        """
        params = {
            "organism": organism,
            "taxonomy_id": taxonomy_id,
            "sp_type": sp_type,
            "evidence": evidence,
            "query_group": "custom",
            **kwargs,
        }

        sources_to_use = sources or ["uniprot"]
        all_records: list[SPRecord] = []

        for source_name in sources_to_use:
            source = self._sources.get(source_name)
            if source is None:
                logger.warning(f"Source '{source_name}' not available")
                continue

            try:
                records = source.query(params)
                all_records.extend(records)
            except Exception as e:
                logger.error(f"Error querying {source.name}: {e}")

        if not all_records:
            logger.warning("No records found")
            return ""

        all_records = deduplicate_records(all_records)
        df = merge_with_existing(all_records)

        output_path = export_csv(df, output)
        summary = generate_summary(df)
        logger.info("\n" + summary)

        return output_path

    def list_available_presets(self) -> list[dict[str, str]]:
        """List all available presets.

        Returns:
            List of preset info dicts.
        """
        return list_presets(self.config)

    def check_sources(self) -> dict[str, bool]:
        """Check availability of all configured sources and predictors.

        Returns:
            Dict mapping source/predictor name to availability status.
        """
        status = {}
        for name, source in self._sources.items():
            try:
                status[name] = source.is_available()
            except Exception:
                status[name] = False

        if self._signalp:
            status["signalp"] = self._signalp.is_available()

        return status
