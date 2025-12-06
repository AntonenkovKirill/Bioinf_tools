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

