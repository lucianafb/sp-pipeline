"""Basic tests for SP-Pipeline."""

import pytest
from sp_pipeline.models import SPRecord, CSV_COLUMNS
from sp_pipeline.dedup import deduplicate_records
from sp_pipeline.utils import extract_cleavage_motif


class TestSPRecord:
    def test_create_basic_record(self):
        record = SPRecord(
            entry_id="P01234",
            source_db="UniProt",
            organism="Homo sapiens",
            protein_name="Test protein",
            full_sequence="MKLSLVAAMLLLLSAARAEEEED",
            sp_start=1,
            sp_end=15,
        )
        assert record.entry_id == "P01234"
        assert record.sp_sequence == "MKLSLVAAMLLLLSA"
        assert record.sp_length == 15

    def test_unique_key(self):
        record = SPRecord(
            entry_id="P01234",
            source_db="UniProt",
            sp_start=1,
            sp_end=20,
        )
        assert record.unique_key == "UniProt|P01234|1|20"

    def test_to_dict(self):
        record = SPRecord(entry_id="P01234", source_db="UniProt")
        d = record.to_dict()
        assert isinstance(d, dict)
        assert d["entry_id"] == "P01234"

    def test_csv_columns_complete(self):
        record = SPRecord(entry_id="P01234", source_db="UniProt")
        d = record.to_dict()
        for col in CSV_COLUMNS:
            assert col in d, f"Missing column: {col}"


class TestDeduplication:
    def test_removes_exact_duplicates(self):
        records = [
            SPRecord(entry_id="P01234", source_db="UniProt", sp_start=1, sp_end=20),
            SPRecord(entry_id="P01234", source_db="UniProt", sp_start=1, sp_end=20),
        ]
        result = deduplicate_records(records)
        assert len(result) == 1

    def test_keeps_different_records(self):
        records = [
            SPRecord(entry_id="P01234", source_db="UniProt", sp_start=1, sp_end=20),
            SPRecord(entry_id="P05678", source_db="UniProt", sp_start=1, sp_end=25),
        ]
        result = deduplicate_records(records)
        assert len(result) == 2

    def test_prefers_experimental_evidence(self):
        records = [
            SPRecord(
                entry_id="P01234", source_db="GenBank", sp_start=1, sp_end=20,
                evidence_method="curated_annotation",
            ),
            SPRecord(
                entry_id="P01234", source_db="UniProt", sp_start=1, sp_end=20,
                evidence_method="experimental",
            ),
        ]
        # Same unique key because source_db differs but entry_id+positions same?
        # Actually unique_key includes source_db, so these are different.
        # Let's test with same source:
        records2 = [
            SPRecord(
                entry_id="P01234", source_db="UniProt", sp_start=1, sp_end=20,
                evidence_method="curated_annotation",
            ),
            SPRecord(
                entry_id="P01234", source_db="UniProt", sp_start=1, sp_end=20,
                evidence_method="experimental",
            ),
        ]
        result = deduplicate_records(records2)
        assert len(result) == 1
        assert result[0].evidence_method == "experimental"


class TestUtils:
    def test_cleavage_motif(self):
        seq = "MKLSLVAAMLLLLSAARAEEEED"
        motif = extract_cleavage_motif(seq, 15, window=3)
        assert "-" in motif
        assert len(motif) > 0

    def test_cleavage_motif_empty(self):
        assert extract_cleavage_motif("", 0) == ""
        assert extract_cleavage_motif("", 5) == ""
