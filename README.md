# BioInforma
BioInforma is tool for working with sequences and FASTQ files.

## How it works:
Open script __main.py__ and use function _main.py_

Function takes list which should contain following objects:

* __data__: it is _list_ of strings of RNA/DNA/protein sequences or _dictionary_ with information from FASTQ-file (see below).
* __type of data__: _string_ 'NA', 'PROTEIN', 'FASTQ'.
* __operations__: _string_ with commands for working with DNA, RNA or protein or _dictionary_ with data for FASTQ filtering (allowed commands and filtering data see below)

Function returns list of sequences or filtered dictionary in case of FASTQ operations.


## DNA/RNA dependent procedures:
1. Transcribe - transform DNA to RNA (ATG to AUG) (___Doesn't work with RNA sequences___)
2. Reverse - return reversed sequence (ATG to GTA)
3. Complement - return complement sequence (TTa to AAt)
4. Reverse complement - return reversed complement sequence (GAT to ATC)
5. Reverse transcription - transform RNA to DNA (AUG to ATG) (___Doesn't work with DNA sequences___)
6. Convert to protein - transform RNA to protein chain coded with one-letter aminoacid names (if you enter DNA sequence, it firstly transcribe to RNA using _Transcribe_ function)
__Pay attention__
* If your sequence doesn't contain Timin and Uracil, it is considered RNA.
* If number of bases in sequence isn't multiple of three, last uncompleted bases are thrown away (e.g. AUGUU transform to M).
* Stop-codons are coded as '.' (dot) in protein sequence. 
* You can enter one or several sequences, but in both cases sequences should be passed as lists.
* Sequences are _case sensitive_.

### These commands are used to call DNA/RNA procedures 
 __('transcribe', 'reverse', 'complement', 'reverse_complement', 'reverse_transcription' or 'protein')__


## Protein sequences dependent procedures:
* 'lengths' - return list with numbers of AA in each sequence(s)
* 'molw' - return list of protein molecular weight (use the average molecular weight of AA, 110 Da)
* 'iso' - return list of approximate isoelectric point of given amino acids sequence
* 'gravy' - return list of GRAVY (grand average of hydropathy) values
* 'rename' - return list of sequences in 3-letter AA code (AA separated by hyphens)
* 'heavy' - return the sequence(s) with maximum molecular weight
* 'light' - return the sequence(s) with minimum molecular weight
* 'heavy_light' - return two lists with sequence(s) with maximum and minimum molecular weight

You can input sequence(s) as list with any strings of sequences. 
__Pay attention__ that your sequence(s) should contain 1-letter symbols (case does not matters) of 20 common amino acids ('U' for selenocysteine and 'O' for pyrrolysine doesn't allowed).

## Working with FASTQ dictionaries:
FASTQ data transfer to function as a dict with following structure:
{'name': ('sequence', 'quality')}
Also, next values for filtering should be send to function as a _dict_:
{'gc_bounds': (min, max),
'length_bounds':(min, max),
'quality_threshold': number}

Function filters FASTQ data according these values:
* GC-bound - lower and upper bound of GC content of sequence (_int_ or _float_) in range 0 - 100%.
* Length_bound - length of sequence (_int_ or _float_) in range 0 - 2^32 bases.
* quality_threshold - threshhold of quality (lower bound) (_int_ or _float_) in range 0 - 100%.

## Inapropriate sequences
If some sequence in list you've sent to function is inappropriate (e.g. DNA contains 'U' letter), this sequence doesn't pass check and will not processed.
This sequence will be showed in terminal as ERROR

See example below.
```python
python3
>>> from main import bioinforma
>>> print(bioinforma([['AGATTT', 'CCCUUU'], 'na', 'transcribe']))
### [2023-10-07 23:39:21,678: WARNING] Ooooops. This sequence CCCUUU is incorrect!
### ['AGAUUU']

## BioInforma using examples 
__(you schould be in folder with main.py script)__

_Nucleic acids_
```python
python3
>>> from main import bioinforma
>>> print(bioinforma([['AGATTT', 'CCCGGG'], 'na', 'transcribe']))
### ['AGAUUU', 'CCCGGG']
>>> print(bioinforma([['AAUUCCC'], 'NA', 'reverse']))
### ['CCCUUAA']

```

_Proteins_
```python
python3
>>> from main import bioinforma
>>> print(bioinforma([['ACGTWWA', 'ILATTWP'], 'protein', 'iso']))
### [5.8, 6.0]
>>> print(bioinforma([['ilattwp'], 'protein', 'gravy']))
### [0.886]
>>> print(bioinforma([['ACGTwwa'], 'protein', 'rename']))
### ['Ala-Cys-Gly-Thr-Trp-Trp-Ala']

```

_FASTQ_
```python
python3
>>> from main import bioinforma
>>> EXAMPLE_FASTQ = {
    '@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D')
    }  

>>> operations = {'gc_bounds':(0, 100), 'length_bounds':(0, 2**32),'quality_threshold': 32}

>>> print(bioinforma([EXAMPLE_FASTQ, 'fastq', operations]))
### {'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD')}

```

## For discussion:
github.com/YuryPopov
