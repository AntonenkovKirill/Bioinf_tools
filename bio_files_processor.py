def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str) -> None:
    """
    Read multiline fasta-file and save new oneline fasta-file
    Arguments:
    input_fasta {str} - input fasta-file
    output_fasta {str} - output fasta-file

    Result: oneline fasta-file saved in output_fasta
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        header = None
        seq_lines = []
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    outfile.write(header + '\n')
                    outfile.write(''.join(seq_lines) + '\n')
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            outfile.write(header + '\n')
            outfile.write(''.join(seq_lines) + '\n')
