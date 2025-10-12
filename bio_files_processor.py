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

def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Read BLAST-output from input_file,
    for every QUERY in "Sequences producing significant alignments:"
    choose Description from the first string,
    sort by alphabet and save into output_file with one column.

    Arguments:
        input_file {str}: input fasta-file
        output_file {str}: output fasta-file
    """
    descriptions = []
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    in_significant_section = False
    first_description_taken = False

    for line in lines:
        line = line.strip()
        if line.startswith("QUERY"):
            first_description_taken = False
            in_significant_section = False
        elif line == "Sequences producing significant alignments:":
            in_significant_section = True
        elif in_significant_section:
            if line == "":
                in_significant_section = False
            elif not first_description_taken:
                descriptions.append(line)
                first_description_taken = True

    descriptions.sort()

    with open(output_file, 'w') as outfile:
        for desc in descriptions:
            outfile.write(desc + '\n')
