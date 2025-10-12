from run_dna_rna_module import reverse
from run_dna_rna_module import is_nucleic_acid
from run_dna_rna_module import transcribe
from run_dna_rna_module import complement
from run_dna_rna_module import reverse_complement
from filter_fastq_module import calc_gc_content
from filter_fastq_module import check_gc_bounds
from filter_fastq_module import check_length_bounds
from filter_fastq_module import mean_phred33
from filter_fastq_module import check_quality
from filter_fastq_module import is_standard_sequence


def run_dna_rna_tools(*args):
    """
    Проверка, транскрипция, разворот или комплементарный перевод генетических последовательностей
    Аргументы:
    sequences / str
    procedure / str

    Возвращает измененную генетическую последовательность или булевый результат выполнения процедуры is_nucleic_acid
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
    seqs: dict,
    gc_bounds=(0, 100),
    length_bounds=(0, 2 ** 32),
    quality_threshold=0
) -> dict:
    """
    Фильтрация fastq-ридов по GC, длине, стандартному составу и качеству.
    На вход: seqs — словарь {name: (sequence, quality)}

    Возвращает результат фильтрации переменной seqs на основании аргументов gc_bounds, length_bounds и quality_threshold
    """
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
    return filtered
