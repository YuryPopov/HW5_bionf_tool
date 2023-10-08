from scripts.dna_tool import dna_rna_tool
from scripts.protein_tool import protein_tool
from scripts.fastq_tool import fastq_tool


def bioinforma(sequence_data: list):
    """
    Function takes list with following structure:
    [[sequences] or {file}, 'type', 'operation' or {values for FASTQ}]
    Function processes protein, DNA, RNA or
    FASTQ file according to functions from
    other modules (see Readme)
    """
    seq_type = sequence_data[1].lower()
    operators = sequence_data[2]
    data = sequence_data[0]

    if seq_type == 'protein':
        return protein_tool(
            seqs=data,
            option=operators)
    elif seq_type == 'fastq':
        return fastq_tool(
            seqs=data,
            gc_bounds=operators['gc_bounds'],
            length_bounds=operators['length_bounds'],
            quality_threshold=operators['quality_threshold'])
    elif seq_type == 'na':
        return dna_rna_tool(data, operators)
    else:
        return 'Enter valid sequence type'
