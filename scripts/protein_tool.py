import logging
import sys

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(stream=sys.stdout)
handler.setFormatter(logging.Formatter
                     (fmt='[%(asctime)s: %(levelname)s] %(message)s'))
logger.addHandler(handler)

AMINO_ACIDS_NAMES = {'A': 'Ala',
                     'R': 'Arg',
                     'N': 'Asn',
                     'D': 'Asp',
                     'V': 'Val',
                     'H': 'His',
                     'G': 'Gly',
                     'Q': 'Gln',
                     'E': 'Glu',
                     'I': 'Ile',
                     'L': 'Leu',
                     'K': 'Lys',
                     'M': 'Met',
                     'P': 'Pro',
                     'S': 'Ser',
                     'Y': 'Tyr',
                     'T': 'Thr',
                     'W': 'Trp',
                     'F': 'Phe',
                     'C': 'Cys'}

GRAVY_AA_VALUES = {'L': 3.8,
                   'K': -3.9,
                   'M': 1.9,
                   'F': 2.8,
                   'P': -1.6,
                   'S': -0.8,
                   'T': -0.7,
                   'W': -0.9,
                   'Y': -1.3,
                   'V': 4.2,
                   'A': 1.8,
                   'R': -4.5,
                   'N': -3.5,
                   'D': -3.5,
                   'C': 2.5,
                   'Q': -3.5,
                   'E': -3.5,
                   'G': -0.4,
                   'H': -3.2,
                   'I': 4.5}

VALID_SYMBOLS = set(AMINO_ACIDS_NAMES)


def calc_gravy(seq: str) -> float:
    """
    Calculate GRAVY (grand average of hydropathy) value
    of given amino acids sequence
    """
    gravy_aa_sum = 0
    for amino_ac in seq:
        gravy_aa_sum += GRAVY_AA_VALUES[amino_ac]
    return round(gravy_aa_sum / len(seq), 3)


def calc_total_charge(charged_amino_ac_numbers_list: list,
                      ph_value: float) -> float:
    """
    Calculate the approximate total charge of some amino acid sequence
    for given pH value
    based only on a list of the number of key charged amino acids.
    """
    n_terminal_charge = 1 / (1 + 10 ** (ph_value - 8.2))
    c_terminal_charge = -1 / (1 + 10 ** (3.65 - ph_value))
    cys_charge = -charged_amino_ac_numbers_list[0] / (1 + 10 ** (8.18 - ph_value))
    asp_charge = -charged_amino_ac_numbers_list[1] / (1 + 10 ** (3.9 - ph_value))
    glu_charge = -charged_amino_ac_numbers_list[2] / (1 + 10 ** (4.07 - ph_value))
    tyr_charge = -charged_amino_ac_numbers_list[3] / (1 + 10 ** (10.46 - ph_value))
    his_charge = charged_amino_ac_numbers_list[4] / (1 + 10 ** (ph_value - 6.04))
    lys_charge = charged_amino_ac_numbers_list[5] / (1 + 10 ** (ph_value - 10.54))
    arg_charge = charged_amino_ac_numbers_list[6] / (1 + 10 ** (ph_value - 12.48))
    total_charge = (n_terminal_charge +
                    c_terminal_charge +
                    cys_charge +
                    asp_charge +
                    glu_charge +
                    tyr_charge +
                    his_charge +
                    lys_charge +
                    arg_charge)
    return total_charge


def calc_iso_point(seq: str):
    """
    Calculate approximate isoelectric point of given amino acids sequence
    """
    charged_amino_ac_numbers = {
        "C": 0, "D": 0, "E": 0, "Y": 0, "H": 0, "K": 0, "R": 0
    }
    for amino_ac in seq:
        if amino_ac in charged_amino_ac_numbers:
            charged_amino_ac_numbers[amino_ac] += 1
    total_charge_tmp = 1
    ph_iso_point = -0.1
    while total_charge_tmp > 0:
        ph_iso_point += 0.1
        total_charge_tmp = calc_total_charge(
            charged_amino_ac_numbers,
            ph_iso_point)
    return round(ph_iso_point, 1)


def transform_to_three_letters(seq: str) -> str:
    """
    Transform 1-letter aminoacid symbols in
    sequence to 3-letter symbols separated by
    hyphens.
    """
    three_letter_seq = ''
    for amino_acid in seq:
        three_letter_seq += AMINO_ACIDS_NAMES[amino_acid] + '-'
    return three_letter_seq[:-1]


def sequence_length(seq: str) -> int:
    """
    Function counts number of aminoacids in
    given sequence
    """
    return len(seq)


def calc_protein_mass(seq: str) -> int:
    """
    Calculate protein molecular weight using the average
    molecular weight of amino acid - 110 Da
    """
    return sequence_length(seq) * 110


def get_min_max_protein(sequences: list[str], operator='both'):
    """
    Function take list of sequences and returns heaviest and lightest
    protein sequences. In case of some heaviest/
    lightest proteins with same masses - returns them.
    Argument 'operator' can be defined by user as 'min' or 'max'.
    In this case function returns only lightest or heaviest protein
    respectively.
    """
    sorted_max = sorted(sequences, key=len, reverse=True)
    sorted_min = sorted_max[::-1]
    list_max = []
    list_min = []
    if len(sequences) == 1:
        return f'''{sequences} - {calc_protein_mass(sequences[0])}
        \n{sequences} - {calc_protein_mass(sequences[0])}'''
    else:
        for index, value in enumerate(sorted_max[:-1]):
            if len(value) == len(sorted_max[index + 1]):
                list_max.append(value)
            else:
                list_max.append(value)
                break
    for index, value in enumerate(sorted_min[:-1]):
        if len(value) == len(sorted_min[index + 1]):
            list_min.append(value)
        else:
            list_min.append(value)
            break

    if operator == 'max':
        return list_max
    elif operator == 'min':
        return list_min
    else:
        return list_max, list_min


def check_sequences(seq: str):
    """
    Throw log notification if at least one sequence
    contains non valid symbols
    """
    letters = set(list(seq.upper()))
    if isinstance(seq, str) and letters.issubset(VALID_SYMBOLS):
        return seq
    else:
        logger.warning((f'Ooooops. This sequence {seq} is incorrect!'))


# Didn't place at the beginning because the functions are defined above
FUNC_STR_INPUT = {
    'gravy': calc_gravy,
    'iso': calc_iso_point,
    'rename': transform_to_three_letters,
    'lengths': sequence_length,
    'molw': calc_protein_mass}

FUNC_LIST_INPUT = {
    'heavy': 'max',
    'light': 'min',
    'heavy_light': 'both'}


def protein_tool(seqs: list[str], option: str):
    """
    Perform some simple operations on amino acids sequences.
    """
    passed_sequences = []
    for seq in seqs:
        checked_seq = check_sequences(seq)
        if checked_seq:
            passed_sequences.append(checked_seq)

    if option in FUNC_STR_INPUT.keys():
        results = []
        for seq in passed_sequences:
            result_tmp = FUNC_STR_INPUT[option](seq.upper())
            results.append(result_tmp)
        return results
    elif option in FUNC_LIST_INPUT.keys():
        return get_min_max_protein(
            passed_sequences,
            operator=FUNC_LIST_INPUT[option])
    else:
        raise ValueError("Enter valid operation")
