from __future__ import annotations
from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass
from typing import ClassVar, Dict, Iterator, Mapping, Union, overload
from pathlib import Path
import gzip
from statistics import mean
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
    - pretty printing: str(seq_obj)
    - alphabet validation: check_alphabet()
    """

    _sequence: str

    def __post_init__(self) -> None:
        if not isinstance(self._sequence, str):
            raise TypeError("Sequence must be a string.")
        if len(self._sequence) == 0:
            raise ValueError("Sequence must be non-empty.")
        if not self.check_alphabet():
            raise ValueError(f"Invalid alphabet for {self.__class__.__name__}: {self._sequence!r}")

    def __len__(self) -> int:
        return len(self._sequence)

    @overload
    def __getitem__(self, key: int) -> str: ...

    @overload
    def __getitem__(self, key: slice) -> "BiologicalSequence": ...

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

    Polymorphism is achieved via subclass class-variables:
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
            raise NotImplementedError("NucleicAcidSequence is abstract; instantiate DNASequence or RNASequence.")
        return set(self._sequence).issubset(self._alphabet)

    def complement(self) -> "NucleicAcidSequence":
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError("NucleicAcidSequence is abstract; instantiate DNASequence or RNASequence.")
        comp = "".join(self._complement_map[base] for base in self._sequence)
        return self.__class__(comp)

    def reverse(self) -> "NucleicAcidSequence":
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError("NucleicAcidSequence is abstract; instantiate DNASequence or RNASequence.")
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
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
    _alphabet: ClassVar[frozenset[str]] = frozenset({"A", "U", "G", "C"})
    _complement_map: ClassVar[Mapping[str, str]] = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
    }


class AminoAcidSequence(BiologicalSequence):
    """
    Protein sequence (amino acids).
    """

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
        """
        Example meaningful method:
        returns amino acid composition (counts).
        """
        return dict(Counter(self._sequence))


def _open_maybe_gzip(path: Path, mode: str):
    """
    Open plain text or gzipped file depending on suffix.
    mode: 'rt' / 'wt' (text modes).
    """
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode)


def _normalize_bounds(bounds, default_low: float, default_high: float) -> tuple[float, float]:
    """
    Convert bounds that may be:
      - a single number: interpreted as (default_low, number)
      - a (low, high) tuple/list
    """
    if bounds is None:
        return default_low, default_high
    if isinstance(bounds, (int, float)):
        return default_low, float(bounds)
    if isinstance(bounds, (tuple, list)) and len(bounds) == 2:
        return float(bounds[0]), float(bounds[1])
    raise ValueError("Bounds must be a number, a 2-tuple/list, or None.")


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0.0, 1.0),
    length_bounds=(0, 2**32),
    quality_threshold: float = 0.0,
) -> dict:
    """
    Filter reads from FASTQ by:
      - GC fraction (0..1), inclusive
      - length, inclusive
      - mean Phred quality, inclusive

    Uses Biopython (SeqIO/SeqRecord/SeqUtils).
    Writes filtered reads to output_fastq.
    Additionally creates '<output_dir>/filtered/<basename>' like your previous behavior.

    Returns summary dict with counts and output path.
    """
    in_path = Path(input_fastq)
    out_path_raw = Path(output_fastq)

    # keep your previous behavior: save into .../filtered/<basename>
    out_dir = out_path_raw.parent / "filtered"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / out_path_raw.name

    gc_low, gc_high = _normalize_bounds(gc_bounds, 0.0, 1.0)
    len_low, len_high = _normalize_bounds(length_bounds, 0.0, float(2**32))

    total = 0
    kept = 0

    def passes_filters(rec: SeqRecord) -> bool:
        seq_len = len(rec.seq)
        if not (len_low <= seq_len <= len_high):
            return False

        # GC fraction in [0..1]
        rec_gc = gc_fraction(rec.seq)
        if not (gc_low <= rec_gc <= gc_high):
            return False

        # Mean Phred quality
        phred = rec.letter_annotations.get("phred_quality")
        if not phred:
            return False
        if mean(phred) < quality_threshold:
            return False

        return True

    with _open_maybe_gzip(in_path, "rt") as fin, _open_maybe_gzip(out_path, "wt") as fout:
        reader = SeqIO.parse(fin, "fastq")
        writer = SeqIO.write

        buffer: list[SeqRecord] = []
        for rec in reader:
            total += 1
            if passes_filters(rec):
                kept += 1
                buffer.append(rec)

                # Write in chunks to avoid large memory usage
                if len(buffer) >= 1000:
                    writer(buffer, fout, "fastq")
                    buffer.clear()

        if buffer:
            writer(buffer, fout, "fastq")

    return {
        "input_fastq": str(in_path),
        "output_fastq": str(out_path),
        "total_reads": total,
        "kept_reads": kept,
        "filtered_out": total - kept,
        "kept_fraction": (kept / total) if total else 0.0,
        "gc_bounds": (gc_low, gc_high),
        "length_bounds": (int(len_low), int(len_high)),
        "quality_threshold": quality_threshold,
    }
