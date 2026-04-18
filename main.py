from __future__ import annotations

import argparse
import gzip
import logging
import sys
from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import ClassVar, Dict, Iterator, Mapping, Union, overload

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

SliceOrInt = Union[slice, int]


@dataclass(frozen=True, slots=True)
class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.

    Provides:
    - len(seq_obj)
    - indexing and slicing: seq_obj[i], seq_obj[i:j:k]
    - iteration support
    - string representation
    - alphabet validation interface
    """

    _sequence: str

    def __post_init__(self) -> None:
        if not isinstance(self._sequence, str):
            raise TypeError("Sequence must be a string.")
        if len(self._sequence) == 0:
            raise ValueError("Sequence must be non-empty.")
        if not self.check_alphabet():
            raise ValueError(
                f"Invalid alphabet for {self.__class__.__name__}: "
                f"{self._sequence!r}"
            )

    def __len__(self) -> int:
        return len(self._sequence)

    @overload
    def __getitem__(self, key: int) -> str:
        ...

    @overload
    def __getitem__(self, key: slice) -> "BiologicalSequence":
        ...

    def __getitem__(self, key: SliceOrInt) -> Union[str, "BiologicalSequence"]:
        if isinstance(key, slice):
            return self.__class__(self._sequence[key])
        return self._sequence[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._sequence)

    def __str__(self) -> str:
        return f"{self.__class__.__name__}('{self._sequence}')"

    def __repr__(self) -> str:
        return str(self)

    @abstractmethod
    def check_alphabet(self) -> bool:
        """Validate sequence symbols against its alphabet."""
        raise NotImplementedError


class NucleicAcidSequence(BiologicalSequence, ABC):
    """
    Base class for nucleic acids (DNA/RNA).

    Implements:
    - check_alphabet()
    - complement()
    - reverse()
    - reverse_complement()

    Polymorphism is achieved via subclass class variables:
    - _alphabet
    - _complement_map
    """

    _alphabet: ClassVar[frozenset[str]]
    _complement_map: ClassVar[Mapping[str, str]]

    def __post_init__(self) -> None:
        object.__setattr__(self, "_sequence", self._sequence.upper())
        super().__post_init__()

    def check_alphabet(self) -> bool:
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError(
                "NucleicAcidSequence is abstract; "
                "instantiate DNASequence or RNASequence."
            )
        return set(self._sequence).issubset(self._alphabet)

    def complement(self) -> "NucleicAcidSequence":
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError(
                "NucleicAcidSequence is abstract; "
                "instantiate DNASequence or RNASequence."
            )
        complemented = "".join(self._complement_map[base] for base in self._sequence)
        return self.__class__(complemented)

    def reverse(self) -> "NucleicAcidSequence":
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError(
                "NucleicAcidSequence is abstract; "
                "instantiate DNASequence or RNASequence."
            )
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    """DNA sequence."""

    _alphabet: ClassVar[frozenset[str]] = frozenset({"A", "T", "G", "C"})
    _complement_map: ClassVar[Mapping[str, str]] = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
    }

    def transcribe(self) -> "RNASequence":
        return RNASequence(self._sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    """RNA sequence."""

    _alphabet: ClassVar[frozenset[str]] = frozenset({"A", "U", "G", "C"})
    _complement_map: ClassVar[Mapping[str, str]] = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
    }


class AminoAcidSequence(BiologicalSequence):
    """Protein sequence."""

    _alphabet: ClassVar[frozenset[str]] = frozenset(
        {
            "A", "C", "D", "E", "F",
            "G", "H", "I", "K", "L",
            "M", "N", "P", "Q", "R",
            "S", "T", "V", "W", "Y",
        }
    )

    def __post_init__(self) -> None:
        object.__setattr__(self, "_sequence", self._sequence.upper())
        super().__post_init__()

    def check_alphabet(self) -> bool:
        return set(self._sequence).issubset(self._alphabet)

    def aa_composition(self) -> Dict[str, int]:
        """Return amino acid composition as a frequency dictionary."""
        return dict(Counter(self._sequence))


def _open_maybe_gzip(path: Path, mode: str):
    """
    Open plain text or gzipped file depending on suffix.

    Parameters
    ----------
    path : Path
        Input or output file path.
    mode : str
        File mode, expected text mode such as 'rt' or 'wt'.

    Returns
    -------
    file object
    """
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode, encoding="utf-8")


def _normalize_bounds(
    bounds,
    default_low: float,
    default_high: float,
) -> tuple[float, float]:
    """
    Normalize filtering bounds.

    Supported formats:
    - None -> (default_low, default_high)
    - single int/float -> (default_low, value)
    - tuple/list of two numeric values -> (low, high)
    """
    if bounds is None:
        return default_low, default_high

    if isinstance(bounds, (int, float)):
        return default_low, float(bounds)

    if isinstance(bounds, (tuple, list)) and len(bounds) == 2:
        low, high = float(bounds[0]), float(bounds[1])
        return low, high

    raise ValueError(
        "Bounds must be None, a single number, or a tuple/list of two numbers."
    )


def _validate_filter_arguments(
    gc_bounds: tuple[float, float],
    length_bounds: tuple[int, int],
    quality_threshold: float,
) -> None:
    """
    Validate filtering parameter ranges.
    """
    gc_low, gc_high = gc_bounds
    len_low, len_high = length_bounds

    if gc_low > gc_high:
        raise ValueError("GC lower bound cannot be greater than GC upper bound.")

    if not (0.0 <= gc_low <= 1.0 and 0.0 <= gc_high <= 1.0):
        raise ValueError("GC bounds must be within the range [0.0, 1.0].")

    if len_low < 0 or len_high < 0:
        raise ValueError("Length bounds must be non-negative integers.")

    if len_low > len_high:
        raise ValueError(
            "Length lower bound cannot be greater than length upper bound."
        )

    if quality_threshold < 0:
        raise ValueError("Quality threshold must be non-negative.")


def setup_logger(log_file: str) -> logging.Logger:
    """
    Configure file logger.

    Parameters
    ----------
    log_file : str
        Path to the log file.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger("fastq_filter_logger")
    logger.setLevel(logging.INFO)

    if logger.handlers:
        logger.handlers.clear()

    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s"
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0.0, 1.0),
    length_bounds=(0, 2**32),
    quality_threshold: float = 0.0,
    logger: logging.Logger | None = None,
) -> dict:
    """
    Filter reads from FASTQ by:
    - GC fraction (inclusive, 0..1)
    - length (inclusive)
    - mean Phred quality (inclusive)

    Supports:
    - .fastq
    - .fastq.gz

    Output is written into:
    <output_parent>/filtered/<output_filename>

    Parameters
    ----------
    input_fastq : str
        Path to input FASTQ/FASTQ.GZ file.
    output_fastq : str
        Output filename/path. Final file is written into `filtered/` directory.
    gc_bounds : tuple[float, float] | float
        GC lower and upper bounds, inclusive.
    length_bounds : tuple[int, int] | int
        Read length lower and upper bounds, inclusive.
    quality_threshold : float
        Minimum mean Phred quality threshold.
    logger : logging.Logger | None
        Optional logger for file logging.

    Returns
    -------
    dict
        Summary statistics dictionary.
    """
    in_path = Path(input_fastq)
    out_path_raw = Path(output_fastq)

    if not in_path.exists():
        message = f"Input FASTQ file does not exist: {in_path}"
        if logger is not None:
            logger.error(message)
        raise FileNotFoundError(message)

    gc_low, gc_high = _normalize_bounds(gc_bounds, 0.0, 1.0)
    len_low, len_high = _normalize_bounds(length_bounds, 0.0, float(2**32))

    _validate_filter_arguments(
        gc_bounds=(gc_low, gc_high),
        length_bounds=(int(len_low), int(len_high)),
        quality_threshold=quality_threshold,
    )

    out_dir = out_path_raw.parent / "filtered"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / out_path_raw.name

    if logger is not None:
        logger.info(
            "Starting FASTQ filtering: input=%s, output=%s, gc_bounds=%s, "
            "length_bounds=%s, quality_threshold=%s",
            in_path,
            out_path,
            (gc_low, gc_high),
            (int(len_low), int(len_high)),
            quality_threshold,
        )

    total_reads = 0
    kept_reads = 0

    def passes_filters(record: SeqRecord) -> bool:
        seq_len = len(record.seq)
        if not (len_low <= seq_len <= len_high):
            return False

        record_gc = gc_fraction(record.seq)
        if not (gc_low <= record_gc <= gc_high):
            return False

        phred_scores = record.letter_annotations.get("phred_quality")
        if not phred_scores:
            return False

        if mean(phred_scores) < quality_threshold:
            return False

        return True

    with _open_maybe_gzip(in_path, "rt") as fin, _open_maybe_gzip(out_path, "wt") as fout:
        buffer: list[SeqRecord] = []

        for record in SeqIO.parse(fin, "fastq"):
            total_reads += 1

            if passes_filters(record):
                kept_reads += 1
                buffer.append(record)

                if len(buffer) >= 1000:
                    SeqIO.write(buffer, fout, "fastq")
                    buffer.clear()

        if buffer:
            SeqIO.write(buffer, fout, "fastq")

    result = {
        "input_fastq": str(in_path),
        "output_fastq": str(out_path),
        "total_reads": total_reads,
        "kept_reads": kept_reads,
        "filtered_out": total_reads - kept_reads,
        "kept_fraction": kept_reads / total_reads if total_reads else 0.0,
        "gc_bounds": (gc_low, gc_high),
        "length_bounds": (int(len_low), int(len_high)),
        "quality_threshold": quality_threshold,
    }

    if logger is not None:
        logger.info(
            "Filtering completed successfully: total_reads=%d, kept_reads=%d, "
            "filtered_out=%d, output=%s",
            result["total_reads"],
            result["kept_reads"],
            result["filtered_out"],
            result["output_fastq"],
        )

    return result


def build_parser() -> argparse.ArgumentParser:
    """
    Build CLI argument parser for FASTQ filtering.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Filter FASTQ/FASTQ.GZ reads by GC content, read length, "
            "and mean Phred quality."
        )
    )

    parser.add_argument(
        "-i",
        "--input-fastq",
        required=True,
        help="Path to input FASTQ or FASTQ.GZ file.",
    )
    parser.add_argument(
        "-o",
        "--output-fastq",
        required=True,
        help=(
            "Output FASTQ filename or path. "
            "The filtered file will be saved into the 'filtered/' directory."
        ),
    )
    parser.add_argument(
        "--gc-min",
        type=float,
        default=0.0,
        help="Minimum GC fraction (default: 0.0).",
    )
    parser.add_argument(
        "--gc-max",
        type=float,
        default=1.0,
        help="Maximum GC fraction (default: 1.0).",
    )
    parser.add_argument(
        "--length-min",
        type=int,
        default=0,
        help="Minimum read length (default: 0).",
    )
    parser.add_argument(
        "--length-max",
        type=int,
        default=2**32,
        help=f"Maximum read length (default: {2**32}).",
    )
    parser.add_argument(
        "-q",
        "--quality-threshold",
        type=float,
        default=0.0,
        help="Minimum mean Phred quality (default: 0.0).",
    )
    parser.add_argument(
        "--log-file",
        default="fastq_filter.log",
        help="Path to log file (default: fastq_filter.log).",
    )

    return parser


def main() -> None:
    """
    CLI entry point for FASTQ filtering.
    """
    parser = build_parser()
    args = parser.parse_args()

    logger = setup_logger(args.log_file)

    try:
        result = filter_fastq(
            input_fastq=args.input_fastq,
            output_fastq=args.output_fastq,
            gc_bounds=(args.gc_min, args.gc_max),
            length_bounds=(args.length_min, args.length_max),
            quality_threshold=args.quality_threshold,
            logger=logger,
        )

        print("FASTQ filtering completed successfully.")
        print(f"Input file: {result['input_fastq']}")
        print(f"Output file: {result['output_fastq']}")
        print(f"Total reads: {result['total_reads']}")
        print(f"Kept reads: {result['kept_reads']}")
        print(f"Filtered out: {result['filtered_out']}")
        print(f"Kept fraction: {result['kept_fraction']:.4f}")

    except Exception as error:
        logger.error("FASTQ filtering failed: %s", error)
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()