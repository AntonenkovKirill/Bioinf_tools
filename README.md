---
title: "Report Project 1"
author: "Group #18"
date: "2026-03-01"
output:
  html_document:
    df_print: default
    highlight: zenburn
    toc: TRUE
    toc_depth: 3
mainfont: NanumGothic
fontsize: 12pt
editor_options:
  markdown:
    wrap: 72
---

# Bioinf_tools (HW16 – OOP Refactoring)

This project provides tools for working with biological sequences
(DNA, RNA, proteins) and for filtering FASTQ files.

In HW16 the project was fully refactored:

- implemented OOP architecture for biological sequences,
- removed legacy procedural functions,
- rewritten FASTQ filtering using Biopython,
- removed old module-based structure (all logic now lives in `main.py`).

---

# 1. OOP Biological Sequences

## BiologicalSequence (abstract base class)

Provides:

- `len(sequence)`
- indexing and slicing (`seq[i]`, `seq[i:j]`)
- iteration support
- string representation
- alphabet validation interface (`check_alphabet()`)

This class defines the common interface for all biological sequences.

---

## NucleicAcidSequence (abstract)

Base class for DNA and RNA.

Implements:

- `complement()`
- `reverse()`
- `reverse_complement()`

Polymorphism is achieved using class-level attributes:

- `_alphabet`
- `_complement_map`

No conditional logic (e.g. `if DNA`) is used inside these methods.

Direct instantiation of `NucleicAcidSequence` is prohibited.

---

## DNASequence

Inherits from `NucleicAcidSequence`.

Additional method:

- `transcribe()` → returns `RNASequence`

---

## RNASequence

Inherits from `NucleicAcidSequence`.

Uses RNA-specific alphabet and complement rules.

---

## AminoAcidSequence

Represents protein sequences.

Additional method:

- `aa_composition()` — returns amino acid frequency dictionary.

---

# 2. FASTQ Filtering (Biopython Implementation)

Function:

```python
filter_fastq(
    input_fastq,
    output_fastq,
    gc_bounds=(0.0, 1.0),
    length_bounds=(0, 2**32),
    quality_threshold=0.0
)
```

Filtering criteria:

- GC fraction (0–1)

- read length

- mean Phred quality score

Implementation details:

- Bio.SeqIO for reading/writing FASTQ

- SeqRecord

- Bio.SeqUtils.gc_fraction

- record.letter_annotations["phred_quality"]

- supports .fastq and .fastq.gz

The function:

- writes output into a filtered/ subdirectory

- returns a summary dictionary with statistics.

---

## Examle Usage

```
from main import DNASequence, filter_fastq

dna = DNASequence("ATGC")
print(dna.reverse_complement())
print(dna.transcribe())

summary = filter_fastq(
    "input.fastq.gz",
    "output.fastq",
    gc_bounds=(0.35, 0.65),
    length_bounds=(50, 300),
    quality_threshold=30
)

print(summary)
```
