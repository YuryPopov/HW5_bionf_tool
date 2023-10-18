def convert_multiline_fasta_to_oneline(
        input_fasta: str, output_fasta: str) -> None:
    if not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'

    with open(input_fasta, 'r') as in_f, open(output_fasta, 'w') as out_f:
        lines = in_f.readlines()
        sequences = []
        current_seq = ''

        for line in lines:
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                current_seq = line
            else:
                current_seq += line

        if current_seq:
            sequences.append(current_seq)

        for sequence in sequences:
            out_f.write(sequence)

    print(f'File {output_fasta} successfully created')
