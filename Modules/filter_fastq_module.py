def calc_gc_content(seq: str) -> float:
    """
    Расчет процента GC для ридов
    Аргумент:
    seq: str

    Возвращает float.
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
    Проверка соответствия GC-состава границам gc_bounds (в процентах, включительно)
    Аргументы:
    seq: str
    gc_bounds: int / float

    Возвращает bool.
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
    Проверка длины последовательности
    Аргументы:
    seq: str
    length_bounds: int / float

    Возвращает bool.
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
    Среднее качество по шкале phred33
    Аргументы:
    qual_str: str

    Возвращает float.
    """
    if not qual_str:
        return 0.0
    qualities = [ord(char) - 33 for char in qual_str]
    return sum(qualities) / len(qualities)


def check_quality(qual_str: str, quality_threshold) -> bool:
    """
    Проверка среднего качества
    Аргументы:
    qual_str: str
    quality_threshold: int / float

    Возвращает bool.
    """
    return mean_phred33(qual_str) >= quality_threshold


def is_standard_sequence(seq: str) -> bool:
    """
    Проверяет, состоит ли последовательность только из стандартных символов.
    Аргументы:
    seq: str

    Возвращает bool.
    """
    valid_bases = set('ATGCatgc')
    return all(base in valid_bases for base in seq)