import os
from Modules.run_dna_rna_module import reverse
from Modules.run_dna_rna_module import is_nucleic_acid
from Modules.run_dna_rna_module import transcribe
from Modules.run_dna_rna_module import complement
from Modules.run_dna_rna_module import reverse_complement
from Modules.filter_fastq_module import calc_gc_content
from Modules.filter_fastq_module import check_gc_bounds
from Modules.filter_fastq_module import check_length_bounds
from Modules.filter_fastq_module import mean_phred33
from Modules.filter_fastq_module import check_quality
from Modules.filter_fastq_module import is_standard_sequence
from Modules.filter_fastq_module import read_fastq
from Modules.filter_fastq_module import write_fastq

def run_dna_rna_tools(*args):
    """
    Check, transcription, reverse or complementary conversion genetic sequencing
    Arguments:
        sequences / str
        procedure / str

    Return a changed genetic sequencing or bool result of "is_nucleic_acid" function
    """
    if len(args) < 2:
        raise ValueError("Необходимо вписать хотя бы одну последовательность и процедуру")
    procedure = args[-1]
    sequences = args[:-1]

    for seq in sequences:
        if not is_nucleic_acid(seq):
            raise ValueError("Некорректная последовательность (смешанные T и U или посторонние символы)")

    results = []
    for seq in sequences:
        if procedure == 'is_nucleic_acid':
            results.append(is_nucleic_acid(seq))
        elif procedure == 'transcribe':
            results.append(transcribe(seq))
        elif procedure == 'reverse':
            results.append(reverse(seq))
        elif procedure == 'complement':
            results.append(complement(seq))
        elif procedure == 'reverse_complement':
            results.append(reverse_complement(seq))
        else:
            raise ValueError("Неизвестная процедура")

    if len(results) == 1:
        return results[0]
    return results


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2 ** 32),
    quality_threshold=0
) -> dict:
    """
    Filtration of fastq-reads from input_fastq по GC, length, standard composition and quality.
    Arguments:
        input_fastq - fastq-file
        output_fastq -

    Return result of filtration of fastq-reads based on variables gc_bounds, length_bounds и quality_threshold
    """
    seqs = read_fastq(input_fastq)
    filtered = {}
    for name, (seq, qual) in seqs.items():
        if not is_standard_sequence(seq):
            continue
        if (
            check_gc_bounds(seq, gc_bounds)
            and check_length_bounds(seq, length_bounds)
            and check_quality(qual, quality_threshold)
        ):
            filtered[name] = (seq, qual)
    write_fastq(filtered, output_fastq)
    return filtered


import os

def calc_gc_content(seq: str) -> float:
    """
    Counting GC in reads
    Arguments:
        seq: str

    Return float.
    """
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    valid_bases = {'A', 'T', 'G', 'C'}
    base_count = sum(seq.count(b) for b in valid_bases)
    if base_count == 0:
        return 0.0
    return gc_count / base_count * 100


def check_gc_bounds(seq: str, gc_bounds) -> bool:
    """
    Checking compliance of GC-composition with gc_bounds (percentage, inclusive)
    Arguments:
        seq: str
        gc_bounds: int / float

    Return bool.
    """
    gc = calc_gc_content(seq)
    if isinstance(gc_bounds, (int, float)):
        return 0 <= gc <= gc_bounds
    elif isinstance(gc_bounds, (tuple, list)) and len(gc_bounds) == 2:
        left, right = gc_bounds
        return left <= gc <= right
    else:
        raise ValueError("gc_bounds должен быть числом или кортежем/списком из двух чисел")


def check_length_bounds(seq: str, length_bounds) -> bool:
    """
    Checking of sequencing length
    Arguments:
        seq: str
        length_bounds: int / float

    Return bool.
    """
    length = len(seq)
    if isinstance(length_bounds, (int, float)):
        return 0 <= length <= length_bounds
    elif isinstance(length_bounds, (tuple, list)) and len(length_bounds) == 2:
        left, right = length_bounds
        return left <= length <= right
    else:
        raise ValueError("length_bounds должен быть числом или кортежем/списком из двух чисел")


def mean_phred33(qual_str: str) -> float:
    """
    Average quality on the phred33 scale
    Arguments:
        qual_str: str

    Return float.
    """
    if not qual_str:
        return 0.0
    qualities = [ord(char) - 33 for char in qual_str]
    return sum(qualities) / len(qualities)


def check_quality(qual_str: str, quality_threshold) -> bool:
    """
    Checking of average quality
    Arguments:
        qual_str: str
        quality_threshold: int / float

    Return bool.
    """
    return mean_phred33(qual_str) >= quality_threshold


def is_standard_sequence(seq: str) -> bool:
    """
    Checks whether the sequence consists only of standard characters
    Arguments:
        seq: str

    Return bool.
    """
    valid_bases = set('ATGCatgc')
    return all(base in valid_bases for base in seq)


def read_fastq(path: str) -> dict:
    """
    Read FASTQ-file and transform it into dict {name: (sequence, quality)}.
    """
    result = {}
    with open(path, "r") as fin:
        while True:
            name_line = fin.readline()
            if not name_line:
                break
            seq_line = fin.readline().strip()
            plus_line = fin.readline()
            qual_line = fin.readline().strip()
            name = name_line.strip()[1:]
            result[name] = (seq_line, qual_line)
    return result


def write_fastq(seqs: dict, output_fastq: str):
    """
    Write dict of reads into FASTQ-file.
    Saving result into subdirectory 'filtered'. If it doesn't exist - create it.
    """
    out_dir = os.path.join(os.path.dirname(output_fastq), 'filtered')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, os.path.basename(output_fastq))
    with open(out_path, "w") as fout:
        for name, (seq, qual) in seqs.items():
            fout.write(f"@{name}\n{seq}\n+\n{qual}\n")

def is_nucleic_acid(seq: str) -> bool:
    seq_upper = seq.upper()
    dna_bases = {'A', 'T', 'G', 'C'}
    rna_bases = {'A', 'U', 'G', 'C'}

    set_seq = set(seq_upper)
    if set_seq.issubset(dna_bases):
        return True
    if set_seq.issubset(rna_bases):
        return True
    return False


def transcribe(seq: str) -> str:
    def transcribe_char(c):
        if c == 'T':
            return 'U'
        if c == 't':
            return 'u'
        return c

    return ''.join(transcribe_char(c) for c in seq)


def reverse(seq: str) -> str:
    return seq[::-1]


def complement(seq: str) -> str:
    complement_map = {
        'A': 'T', 'a': 't',
        'T': 'A', 't': 'a',
        'U': 'A', 'u': 'a',
        'G': 'C', 'g': 'c',
        'C': 'G', 'c': 'g'
    }
    return ''.join(complement_map.get(base, base) for base in seq)


def reverse_complement(seq: str) -> str:
    return reverse(complement(seq))
