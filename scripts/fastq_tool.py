def count_gc_content(sequence: str) -> float:
    counter = 0
    for base in sequence:
        if base == 'G' or base == 'C':
            counter += 1
    gc_content = counter / len(sequence) * 100
    return round(gc_content, 2)


def check_quality(seq_quality: str, threshold: float) -> bool:
    quality_ascii = []
    for char in seq_quality:
        quality_ascii.append(ord(char))

    quality_perc = list(map(lambda x: int(x) - 33, quality_ascii))
    mean_quality = sum(quality_perc) / len(quality_perc)

    return True if mean_quality > threshold else False


def fastq_tool(
        seqs: dict,
        gc_bounds=(0, 100),
        length_bounds=(0, 2**32),
        quality_threshold=0) -> dict:
    '''
    Function takes dict of FASTQ data:
    {"name of seq": ("sequence", "sequence_quality")}, and some
    conditions for filtering:
    - boundaries of GC-content(default value (0, 100)),
    - boundaries of length (default value (0, 2**32),
    - quality treshold (default value 0%).
    Function calculates needed indicators, and filters dict
    of sequences according to these indicators.
    Function returns filtered dictionary, which consists of only
    sequences, which meet the conditions.
    '''
    result_dict = {}

    for name, data in seqs.items():
        gc_content = count_gc_content(data[0])
        seq_length = len(data[0])
        quality = check_quality(
            seq_quality=data[1],
            threshold=quality_threshold)

        if (gc_bounds[0] < gc_content < gc_bounds[1]
                and length_bounds[0] < seq_length < length_bounds[1]
                and quality):
            result_dict[name] = data

    return result_dict
