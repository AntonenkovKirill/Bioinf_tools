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