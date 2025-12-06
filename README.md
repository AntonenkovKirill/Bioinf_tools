---
title: "Report Project 1"
author: "Group #18"
date: "2025-10-01"
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

This project contains a set of Python-scripts and modules for analysis and
filtering of genetic sequences, with an emphasis on working with
DNA/RNA and fastq data.

## Functions

**run_dna_rna_tools** — universal procedure for working with
DNA/RNA sequences.

**filter_fastq** — filter the rows in a fastq-file by different
biological characteristics.

**convert_multiline_fasta_to_oneline** — reads input
fasta-file in which the sequence (DNA/RNA/protein) can be
is broken up into several lines. Then saves to a new fasta-file in
where each sequence fits into one string.

**parse_blast_output** — processes the BLAST output as a txt file,
for each query in the Sequences producing significant
alignments selects Description from the first line of matches, sorts
names alphabetically and saves in output_file with a single column.

The project is structured modularly: all supporting functions are
separate files for reuse and ease of testing.

### Files structure

**main.py** — main script, contains the **filter_fastq_module.py** and **run_dna_rna_module.py** functions and examples of their
call.

**filter_fastq_module.py** — module for filter of
fastq-sequences.

**run_dna_rna_module.py** — module c processing functions
DNA/RNA sequences.

**bio_files_processor.py** — a script containing functions
*convert_multiline_fasta_to_oneline* и *parse_blast_output.*

## run_dna_rna_tools function

Allows you to work with multiple DNA or RNA sequences at once,
by applying one of the available processing procedures.

### Main procedures

**is_nucleic_acid** — returns a bool result: is
sequence a valid (DNA only or RNA only, cannot be mixed
T and U).

**transcribe** — transcribes DNA to RNA (replaces T/t with U/u).

**reverse** — reverse the sequence.

**complement** — returns the complementary sequence (saving
case of symbols).

**reverse_complement** — returns the complementary inverse
sequence.

Arguments: Any number of string sequences is passed,
The last argument is the name of the procedure. If you submit one
sequence, returns the string; if multiple - list.

### Usage

```         
python 

run_dna_rna_tools('ATG', 'reverse') \# 'GTA'
run_dna_rna_tools('ATGC', 'AGTC', 'reverse') \# ['CGTA', 'CTGA']
```

## Function filter_fastq

Filters the rows contained in the fastq file as specified
biological criteria: percentage of GC, read length and mean quality
(phred33).

Arguments:

**input_fastq** — fastq-file sent to input.

**output_fastq** — a fastq-file containing only those reeds that
satisfy conditions expressed as arguments *gc_bounds,
length_bounds и quality_threshold*.

**seqs** — dictionary of the type {name: (sequence string, string
quality)}

**gc_bounds** — range for GC percent (default (0, 100)). Can
to pass a single number - it is considered that it is maximum, minimum = 0.

**length_bounds** — length range (default (0, 2 * *32)).

**quality_threshold** — the average quality threshold (default is 0).

The reads must meet all filter conditions:

-   Consist of standard nucleotide characters only.

-   GC-composition is in range.

-    Length is within range.

-   verage quality is not lower than the threshold.

Result: returns the same file format as input, but with
by filtered ridges.

### Usage

```         
python 

filter_fastq(
input_fastq, 
output_fastq, 
gc_bounds=(40, 60), 
length_bounds=(100, 300),
quality_threshold=30
) 
```

## Function convert_multiline_fasta_to_oneline

Translation of genetic sequences recorded in the fasta-file in multi-line format to single-line format.

**input_fasta** — fasta-file sent to input.

**output_fasta** — Final fasta-file.

Result: returns fasta-file with same genetic
sequences, but written in a straight line.

## Docstrings and documentation

All modules and scripts contain detailed annotations and documentation for
each function describing input and output, as well as the logic. It is recommended to use the source code of the modules to get
additional information on the implementation of each procedure.

