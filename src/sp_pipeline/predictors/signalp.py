"""SignalP 6.0 predictor via BioLib cloud API.

Uses pybiolib to submit sequences to DTU/SignalP-6 on BioLib,
poll for results, and parse the prediction summary.
"""

import logging
import os
import re
import tempfile
from pathlib import Path
from typing import Optional

from sp_pipeline.models import SPRecord
from sp_pipeline.utils import compute_sp_features, extract_cleavage_motif

logger = logging.getLogger("sp_pipeline.predictors.signalp")

_SIGNALP_SP_TYPE_MAP = {
    "SP(Sec/SPI)": "type1_cleaved",
    "LIPO(Sec/SPII)": "type1_cleaved",
    "TAT(Tat/SPI)": "type1_cleaved",
    "TATLIPO(Tat/SPII)": "type1_cleaved",
    "PILIN(Sec/SPIII)": "type1_cleaved",
}

_ORGANISM_MAP = {
    "eukarya": "euk",
    "euk": "euk",
    "other": "other",
    "gram+": "other",
    "gram-": "other",
    "archaea": "other",
}

MIN_SEQUENCE_LENGTH = 10
MAX_SEQUENCE_LENGTH = 10_000


class SignalPPredictor:
    """Predict signal peptides using SignalP 6.0 via BioLib."""

    BIOLIB_APP = "DTU/SignalP-6"

    def __init__(
        self,
        api_token: str,
        organism: str = "eukarya",
        batch_size: int = 500,
    ):
        self.api_token = api_token
        self.organism = _ORGANISM_MAP.get(organism.lower(), "euk")
        self.batch_size = min(batch_size, 1000)

    def is_available(self) -> bool:
        """Check if the predictor is configured with a valid token."""
        if not self.api_token:
            return False
        try:
            import biolib  # noqa: F401
            return True
        except ImportError:
            logger.warning("pybiolib not installed. Run: pip install pybiolib")
            return False

    def predict(self, sequences: list[dict]) -> list[SPRecord]:
        """Run SignalP predictions on a list of sequences.

        Args:
            sequences: List of dicts with at least 'entry_id' and 'full_sequence'.
                       Optional keys: 'organism', 'query_group', 'source_db',
                       'protein_name', 'gene_name', 'taxonomy_id'.

        Returns:
            SPRecord list — only entries where SignalP finds a signal peptide,
            with evidence_method='predicted_signalp'.
        """
        valid_seqs = [
            s for s in sequences
            if s.get("full_sequence")
            and MIN_SEQUENCE_LENGTH <= len(s["full_sequence"]) <= MAX_SEQUENCE_LENGTH
        ]

        if not valid_seqs:
            logger.warning("No valid sequences to predict (need %d–%d aa)", MIN_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH)
            return []

        logger.info(
            "Submitting %d sequences to SignalP 6.0 (organism=%s, batches of %d)",
            len(valid_seqs), self.organism, self.batch_size,
        )

        seq_lookup = {s["entry_id"]: s for s in valid_seqs}
        all_predictions: list[SPRecord] = []

        for batch_idx in range(0, len(valid_seqs), self.batch_size):
            batch = valid_seqs[batch_idx : batch_idx + self.batch_size]
            batch_num = batch_idx // self.batch_size + 1
            total_batches = (len(valid_seqs) + self.batch_size - 1) // self.batch_size
            logger.info("Batch %d/%d (%d sequences)", batch_num, total_batches, len(batch))

            try:
                raw_output = self._run_batch(batch)
                predictions = self._parse_results(raw_output, seq_lookup)
                all_predictions.extend(predictions)
            except Exception as e:
                logger.error("SignalP batch %d failed: %s", batch_num, e)

        return all_predictions

    def _run_batch(self, sequences: list[dict]) -> str:
        """Submit a batch to BioLib and return the prediction summary text."""
        import biolib

        old_token = os.environ.get("BIOLIB_TOKEN")
        try:
            os.environ["BIOLIB_TOKEN"] = self.api_token

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False,
            ) as fasta_file:
                fasta_path = fasta_file.name
                for seq in sequences:
                    fasta_file.write(f">{seq['entry_id']}\n")
                    full_seq = seq["full_sequence"]
                    for i in range(0, len(full_seq), 80):
                        fasta_file.write(full_seq[i : i + 80] + "\n")

            app = biolib.load(self.BIOLIB_APP)
            args = (
                f"--fastafile {fasta_path} "
                f"--organism {self.organism} "
                f"--format txt "
                f"--mode fast"
            )
            logger.debug("BioLib args: %s", args)
            job = app.cli(args=args)

            output_text = self._extract_prediction_text(job)
            return output_text
        finally:
            if old_token is not None:
                os.environ["BIOLIB_TOKEN"] = old_token
            elif "BIOLIB_TOKEN" in os.environ:
                del os.environ["BIOLIB_TOKEN"]
            if "fasta_path" in locals():
                Path(fasta_path).unlink(missing_ok=True)

    def _extract_prediction_text(self, job) -> str:
        """Extract prediction summary text from a BioLib job result.

        Tries output files first, then falls back to stdout.
        """
        try:
            output_files = job.list_output_files()
            logger.debug("BioLib output files: %s", [f.name if hasattr(f, 'name') else str(f) for f in output_files])

            for candidate in output_files:
                fname = candidate.name if hasattr(candidate, "name") else str(candidate)
                if "prediction_results" in fname and fname.endswith(".txt"):
                    file_obj = job.get_output_file(fname)
                    if hasattr(file_obj, "get_file_handle"):
                        return file_obj.get_file_handle().read().decode("utf-8", errors="replace")
                    if isinstance(file_obj, bytes):
                        return file_obj.decode("utf-8", errors="replace")
                    return str(file_obj)
        except Exception as e:
            logger.debug("Could not read output files, falling back to stdout: %s", e)

        stdout = job.get_stdout()
        if isinstance(stdout, bytes):
            return stdout.decode("utf-8", errors="replace")
        return str(stdout)

    def _parse_results(
        self,
        output_text: str,
        seq_lookup: dict[str, dict],
    ) -> list[SPRecord]:
        """Parse SignalP prediction summary into SPRecord instances.

        Only returns records where a signal peptide was predicted (not OTHER).
        """
        records: list[SPRecord] = []
        if not output_text.strip():
            logger.warning("Empty SignalP output")
            return records

        for line in output_text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parsed = self._parse_prediction_line(line)
            if parsed is None:
                continue

            entry_id = parsed["entry_id"]
            prediction = parsed["prediction"]

            if prediction == "OTHER":
                continue

            seq_data = seq_lookup.get(entry_id)
            if seq_data is None:
                logger.debug("Sequence '%s' not found in lookup — skipping", entry_id)
                continue

            full_sequence = seq_data.get("full_sequence", "")
            sp_type = self._map_sp_type(prediction)
            cs_pos = parsed.get("cs_pos")

            sp_start = 1
            sp_end = cs_pos if cs_pos else 0
            sp_sequence = full_sequence[sp_start - 1 : sp_end] if (full_sequence and sp_end > 0) else ""

            features = compute_sp_features(sp_sequence)
            cleavage_motif = (
                extract_cleavage_motif(full_sequence, sp_end) if (full_sequence and sp_end > 0) else None
            )

            record = SPRecord(
                entry_id=entry_id,
                source_db=seq_data.get("source_db", "predicted"),
                organism=seq_data.get("organism", ""),
                taxonomy_id=seq_data.get("taxonomy_id"),
                protein_name=seq_data.get("protein_name", ""),
                gene_name=seq_data.get("gene_name"),
                full_sequence=full_sequence,
                sp_sequence=sp_sequence,
                sp_start=sp_start,
                sp_end=sp_end,
                sp_type=sp_type,
                cleavage_site_motif=cleavage_motif,
                evidence_method="predicted_signalp",
                query_group=seq_data.get("query_group", ""),
                hydrophobicity_mean=features.get("hydrophobicity_mean"),
                net_charge_ph7=features.get("net_charge_ph7"),
                n_region=features.get("n_region"),
                h_region=features.get("h_region"),
                c_region=features.get("c_region"),
            )
            records.append(record)

        logger.info("SignalP predicted %d signal peptides from output", len(records))
        return records

    def _parse_prediction_line(self, line: str) -> Optional[dict]:
        """Parse a single tab-separated prediction line.

        Expected format (Eukarya):
            ID  Prediction  OTHER  SP(Sec/SPI)  CS Position
        With CS Position like: "CS pos: 19-20. Pr: 0.9780"

        Returns dict with entry_id, prediction, probabilities, cs_pos, cs_prob.
        """
        parts = line.split("\t")
        if len(parts) < 3:
            return None

        entry_id = parts[0].strip()
        prediction = parts[1].strip()

        cs_pos = None
        cs_prob = None
        cs_field = parts[-1].strip() if len(parts) >= 4 else ""
        cs_match = re.search(r"CS pos:\s*(\d+)-(\d+)\.\s*Pr:\s*([\d.]+)", cs_field)
        if cs_match:
            cs_pos = int(cs_match.group(1))
            cs_prob = float(cs_match.group(3))

        return {
            "entry_id": entry_id,
            "prediction": prediction,
            "cs_pos": cs_pos,
            "cs_prob": cs_prob,
        }

    @staticmethod
    def _map_sp_type(prediction: str) -> str:
        """Map SignalP prediction label to internal sp_type."""
        if prediction in ("SP", "SP(Sec/SPI)", "TAT", "TAT(Tat/SPI)"):
            return "type1_cleaved"
        if prediction in ("LIPO", "LIPO(Sec/SPII)", "TATLIPO", "TATLIPO(Tat/SPII)"):
            return "type1_cleaved"
        if prediction in ("PILIN", "PILIN(Sec/SPIII)"):
            return "type1_cleaved"
        return "type1_cleaved"
