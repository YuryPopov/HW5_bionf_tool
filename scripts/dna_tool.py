import logging
import sys

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(stream=sys.stdout)
handler.setFormatter(logging.Formatter
                     (fmt='[%(asctime)s: %(levelname)s] %(message)s'))
logger.addHandler(handler)

DNA_BASES = ['a', 'g', 'c', 't', 'A', 'G', 'C', 'T']
RNA_BASES = ['a', 'g', 'c', 'u', 'A', 'G', 'C', 'U']
COMMAND_LIST = [
    'transcribe',
    'reverse',
    'complement',
    'reverse_complement',
    'protein'
]
WARNING = 'Enter valid sequence'

DNA_TO_RNA = {
    'a': 'a',
    'A': 'A',
    'c': 'C',
    'C': 'C',
    't': 'u',
    'T': 'U',
    'g': 'g',
    'G': 'G'
}

RNA_TO_DNA = {
    'a': 'a',
    'A': 'A',
    'c': 'C',
    'C': 'C',
    'u': 't',
    'U': 'T',
    'g': 'g',
    'G': 'G'
}

dna_complement = {
    'a': 't',
    'A': 'T',
    'c': 'g',
    'C': 'G',
    't': 'a',
    'T': 'A',
    'g': 'c',
    'G': 'C'
}

aminoacids = {
    'GCU': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',
    'CGU': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',
    'AGA': 'R',
    'AGG': 'R',
    'AAU': 'N',
    'AAC': 'N',
    'GAU': 'D',
    'GAC': 'D',
    'UGU': 'C',
    'UGC': 'C',
    'CAA': 'Q',
    'CAG': 'Q',
    'GAA': 'E',
    'GAG': 'E',
    'GGU': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G',
    'CAU': 'H',
    'CAC': 'H',
    'AUU': 'I',
    'AUC': 'I',
    'AUA': 'I',
    'UUA': 'L',
    'UUG': 'L',
    'CUU': 'L',
    'CUC': 'L',
    'CUA': 'L',
    'CUG': 'L',
    'AAA': 'K',
    'AAG': 'K',
    'AUG': 'M',
    'UUU': 'F',
    'UUC': 'F',
    'CCU': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',
    'UCU': 'S',
    'UCC': 'S',
    'UCA': 'S',
    'UCG': 'S',
    'AGU': 'S',
    'AGC': 'S',
    'ACU': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',
    'UGG': 'W',
    'UAU': 'Y',
    'UAC': 'Y',
    'GUU': 'V',
    'GUC': 'V',
    'GUA': 'V',
    'GUG': 'V',
    'UAA': '.',
    'UAG': '.',
    'UGA': '.'
    }


def check_sequences(sequence: str) -> str:
    '''
    Function strip spaces before and after sequence,
    and check, that sequence is valid (contains only
    bases' symbols and match with DNA or RNA bases)
    '''
    check = set(list(sequence))
    if check.issubset(DNA_BASES) or check.issubset(RNA_BASES):
        return True
    else:
        logger.warning((f'Ooooops. This sequence {sequence} is incorrect!'))


def transcribe(sequence: str) -> str:
    bases = set(list(sequence))
    rna_sequence = ''
    if bases.issubset(DNA_BASES):
        for base in sequence:
            rna_sequence += DNA_TO_RNA[base]
    else:
        logger.warning((f'Ooooops. This sequence {sequence} is incorrect!'))

    if rna_sequence:
        return rna_sequence


def reverse(sequence: str) -> str:
    return sequence[::-1]


def complement(sequence: str) -> str:
    complement_sequence = ''
    for base in sequence:
        complement_sequence += dna_complement[base]
    return complement_sequence


def reverse_complement(sequence: str) -> str:
    complement_seq = complement(sequence)
    return complement_seq[::-1]


def convert_to_protein(sequence: str) -> str:
    '''
    This function convert RNA sequence to protein.
    Function creates list of strings (codon_chain) with codons.
    Finally, function searches for a codon from the chain among the dictionary
    values and returns the corresponding key.
    '''
    bases = set(list(sequence))
    if not bases.issubset(RNA_BASES):
        logger.warning((f'Ooooops. This sequence {sequence} is incorrect!'))
    else:
        codon_chain = []
        codon = ''
        for base in sequence:
            codon += base
            if len(codon) == 3:
                codon_chain.append(codon)
                codon = ''

        protein_sequence = []
        for codon in codon_chain:
            if codon in aminoacids.keys():
                protein_sequence.append(aminoacids[codon])

        end_protein_sequence = ''.join(protein_sequence)
        if end_protein_sequence:
            return end_protein_sequence


def reverse_transcription(sequence: str) -> str:
    bases = set(list(sequence))
    dna_sequence = ''
    if bases.issubset(RNA_BASES):
        for base in sequence:
            dna_sequence += RNA_TO_DNA[base]
        return dna_sequence
    else:
        logger.warning((f'Ooooops. This sequence {sequence} is incorrect!'))


def dna_rna_tool(*args):
    sequences, command = args
    results = []
    passed_sequences = []
    for seq in sequences:
        checked_seq = check_sequences(seq.strip())
        if checked_seq:
            passed_sequences.append(seq.strip())

    d_of_functions = {
                    'transcribe': transcribe,
                    'reverse': reverse,
                    'reverse_complement': reverse_complement,
                    'complement': complement,
                    'reverse_transcription': reverse_transcription,
                    'protein': convert_to_protein}
    for seq in passed_sequences:
        result = d_of_functions[command.lower()](seq)
        if result:
            results.append(result)

    return results
