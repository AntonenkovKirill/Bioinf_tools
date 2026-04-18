from __future__ import annotations

import sys 
import gzip
import logging
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append(str(Path(__file__).resolve().parents[1]))

from main import (
    _normalize_bounds,
    filter_fastq,
    setup_logger,
)


def make_fastq_record(
    seq: str,
    record_id: str,
    phred_quality: list[int],
) -> SeqRecord:
    """
    Create a SeqRecord suitable for FASTQ writing.
    """
    record = SeqRecord(Seq(seq), id=record_id, description="")
    record.letter_annotations["phred_quality"] = phred_quality
    return record


def write_fastq(path: Path, records: list[SeqRecord]) -> None:
    """
    Write records to plain FASTQ.
    """
    with path.open("w", encoding="utf-8") as handle:
        SeqIO.write(records, handle, "fastq")


def write_fastq_gz(path: Path, records: list[SeqRecord]) -> None:
    """
    Write records to gzipped FASTQ.
    """
    with gzip.open(path, "wt") as handle:
        SeqIO.write(records, handle, "fastq")


def read_fastq_ids(path: Path) -> list[str]:
    """
    Read FASTQ file and return record IDs.
    Supports plain and gzipped FASTQ.
    """
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            return [record.id for record in SeqIO.parse(handle, "fastq")]

    with path.open("r", encoding="utf-8") as handle:
        return [record.id for record in SeqIO.parse(handle, "fastq")]


class TestFilterFastqCore:
    def test_filter_fastq_keeps_only_records_passing_all_filters(
        self,
        tmp_path: Path,
    ) -> None:
        """
        Test that only records satisfying GC, length and quality filters are kept.
        """
        input_path = tmp_path / "input.fastq"
        output_path = tmp_path / "result.fastq"

        records = [
            make_fastq_record("GCGC", "keep_me", [40, 40, 40, 40]),
            make_fastq_record("ATAT", "low_gc", [40, 40, 40, 40]),
            make_fastq_record("GCGCGC", "too_long", [40] * 6),
            make_fastq_record("GGCC", "low_quality", [10, 10, 10, 10]),
        ]
        write_fastq(input_path, records)

        result = filter_fastq(
            input_fastq=str(input_path),
            output_fastq=str(output_path),
            gc_bounds=(0.5, 1.0),
            length_bounds=(4, 4),
            quality_threshold=30,
        )

        filtered_file = tmp_path / "filtered" / "result.fastq"
        kept_ids = read_fastq_ids(filtered_file)

        assert filtered_file.exists()
        assert kept_ids == ["keep_me"]
        assert result["total_reads"] == 4
        assert result["kept_reads"] == 1
        assert result["filtered_out"] == 3

    def test_filter_fastq_creates_output_inside_filtered_directory(
        self,
        tmp_path: Path,
    ) -> None:
        """
        Test that output file is created in filtered/ subdirectory.
        """
        input_path = tmp_path / "input.fastq"
        output_path = tmp_path / "my_output.fastq"

        records = [
            make_fastq_record("GCGC", "r1", [35, 35, 35, 35]),
        ]
        write_fastq(input_path, records)

        result = filter_fastq(
            input_fastq=str(input_path),
            output_fastq=str(output_path),
        )

        expected_output = tmp_path / "filtered" / "my_output.fastq"

        assert expected_output.exists()
        assert Path(result["output_fastq"]) == expected_output

    def test_filter_fastq_supports_gzipped_input_and_output(
        self,
        tmp_path: Path,
    ) -> None:
        """
        Test support for .fastq.gz input and output.
        """
        input_path = tmp_path / "input.fastq.gz"
        output_path = tmp_path / "result.fastq.gz"

        records = [
            make_fastq_record("GCGC", "r1", [40, 40, 40, 40]),
            make_fastq_record("ATAT", "r2", [40, 40, 40, 40]),
        ]
        write_fastq_gz(input_path, records)

        result = filter_fastq(
            input_fastq=str(input_path),
            output_fastq=str(output_path),
            gc_bounds=(0.75, 1.0),
            length_bounds=(4, 4),
            quality_threshold=30,
        )

        filtered_file = tmp_path / "filtered" / "result.fastq.gz"
        kept_ids = read_fastq_ids(filtered_file)

        assert filtered_file.exists()
        assert kept_ids == ["r1"]
        assert result["kept_reads"] == 1

    def test_filter_fastq_returns_zero_fraction_for_empty_input(
        self,
        tmp_path: Path,
    ) -> None:
        """
        Test that empty input FASTQ is handled correctly.
        """
        input_path = tmp_path / "empty.fastq"
        output_path = tmp_path / "result.fastq"

        write_fastq(input_path, [])

        result = filter_fastq(
            input_fastq=str(input_path),
            output_fastq=str(output_path),
        )

        filtered_file = tmp_path / "filtered" / "result.fastq"

        assert filtered_file.exists()
        assert result["total_reads"] == 0
        assert result["kept_reads"] == 0
        assert result["filtered_out"] == 0
        assert result["kept_fraction"] == 0.0
        assert read_fastq_ids(filtered_file) == []

def test_filter_fastq_filters_by_quality_threshold(
    tmp_path: Path,
) -> None:
    """
    Test that records below quality threshold are filtered out.
    """
    input_path = tmp_path / "input.fastq"
    output_path = tmp_path / "result.fastq"

    records = [
        make_fastq_record("GCGC", "high_quality", [40, 40, 40, 40]),
        make_fastq_record("GCGC", "low_quality", [10, 10, 10, 10]),
    ]
    write_fastq(input_path, records)

    result = filter_fastq(
        input_fastq=str(input_path),
        output_fastq=str(output_path),
        quality_threshold=30,
    )

    filtered_file = tmp_path / "filtered" / "result.fastq"
    kept_ids = read_fastq_ids(filtered_file)

    assert kept_ids == ["high_quality"]
    assert result["kept_reads"] == 1
    assert result["filtered_out"] == 1


class TestErrorsAndValidation:
    def test_filter_fastq_raises_file_not_found_for_missing_input(
        self,
        tmp_path: Path,
    ) -> None:
        """
        Test that missing input file raises FileNotFoundError.
        """
        missing_input = tmp_path / "missing.fastq"
        output_path = tmp_path / "result.fastq"

        with pytest.raises(FileNotFoundError):
            filter_fastq(
                input_fastq=str(missing_input),
                output_fastq=str(output_path),
            )

    def test_filter_fastq_raises_value_error_for_invalid_gc_bounds(
        self,
        tmp_path: Path,
    ) -> None:
        """
        Test that invalid GC bounds raise ValueError.
        """
        input_path = tmp_path / "input.fastq"
        output_path = tmp_path / "result.fastq"

        records = [
            make_fastq_record("GCGC", "r1", [40, 40, 40, 40]),
        ]
        write_fastq(input_path, records)

        with pytest.raises(ValueError, match="GC lower bound"):
            filter_fastq(
                input_fastq=str(input_path),
                output_fastq=str(output_path),
                gc_bounds=(0.9, 0.1),
            )

    def test_normalize_bounds_accepts_number_and_returns_default_low(
        self,
    ) -> None:
        """
        Test helper normalization for single numeric bound.
        """
        result = _normalize_bounds(10, 0.0, 100.0)
        assert result == (0.0, 10.0)


class TestLogging:
    def test_setup_logger_and_filter_fastq_write_messages_to_log_file(
        self,
        tmp_path: Path,
    ) -> None:
        """
        Test that logger writes INFO messages to log file during filtering.
        """
        input_path = tmp_path / "input.fastq"
        output_path = tmp_path / "result.fastq"
        log_path = tmp_path / "logs" / "fastq_filter.log"

        records = [
            make_fastq_record("GCGC", "r1", [40, 40, 40, 40]),
        ]
        write_fastq(input_path, records)

        logger = setup_logger(str(log_path))
        assert isinstance(logger, logging.Logger)

        filter_fastq(
            input_fastq=str(input_path),
            output_fastq=str(output_path),
            logger=logger,
        )

        assert log_path.exists()

        log_text = log_path.read_text(encoding="utf-8")
        assert "INFO" in log_text
        assert "Starting FASTQ filtering" in log_text
        assert "Filtering completed successfully" in log_text