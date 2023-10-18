import os
from scripts.dna_tool import dna_rna_tool
from scripts.protein_tool import protein_tool
from scripts.fastq_tool import fastq_tool


def transform_file_to_dict(input_path: str) -> dict:
    result_dict = dict()

    with open(input_path, 'r') as file:
        lines = file.readlines()
        i = 0
        while i < len(lines):
            header = lines[i].strip()
            sequence = lines[i + 1].strip()
            quality = lines[i + 3].strip()
            result_dict[header] = (sequence, quality)
            i += 4

    return result_dict
    

def save_to_file(result, output_filename) -> None:
    output_directory = 'fastq_filtration_result'
    os.makedirs(output_directory, exist_ok=True)
    output_path = os.path.join(output_directory, output_filename)
    output_path = './fastq_filtration_result/' + output_filename + '.fastq'
    with open(output_path, 'w') as output_file:
        for header, (sequence, quality) in result.items():
            output_file.write(header + '\n')
            output_file.write(sequence + '\n')
            output_file.write('+\n')
            output_file.write(quality + '\n')


def bioinforma(sequence_data: list) -> None:
    """
    Function takes list of following arguments:
    'data' - contains path to source fastq file or data
    'type' - type of sequence,
    'operation' - it describes procedure
    'output_filename' - for saving new fastq file
    Function processes protein, DNA, RNA or
    FASTQ file according to functions from
    other modules (see Readme)
    """

    seq_type = sequence_data[1].lower()
    operators = sequence_data[2]
    data = sequence_data[0]
    if seq_type == 'fastq':
        data = transform_file_to_dict(data)
        if sequence_data[:-1] == operators:
            output_filename = os.path.basename(data)[:6]
        else:
            output_filename = sequence_data[:-1]

    if seq_type == 'protein':
        return protein_tool(
            seqs=data,
            option=operators)
    elif seq_type == 'fastq':
        return save_to_file(fastq_tool(
            seqs=data,
            gc_bounds=operators['gc_bounds'],
            length_bounds=operators['length_bounds'],
            quality_threshold=operators['quality_threshold']), output_filename)
    elif seq_type == 'na':
        return dna_rna_tool(data, operators)
    else:
        return 'Enter valid sequence type'
