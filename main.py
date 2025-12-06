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

