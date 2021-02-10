#!/usr/bin/env python3

#Python Modules
import numpy as np
from subprocess import run, PIPE
import os
import sys
from collections import defaultdict
import csv

#Modules from this package


#Modules that need to be installed
try:
    from Bio import SeqIO
except ImportError:
    print('Error: Biopython not found, can be installed with\npip3 install biopython', file=sys.stderr)
    sys.exit(1)

try:
    import pysam
except ImportError:
    print('Error: Pysam not found, can be installed with\npip3 install pysam', file=sys.stderr)
    sys.exit(1)

#Set-up arguments...
from argparse import ArgumentParser, SUPPRESS

parser = ArgumentParser(description='Align an reference RNA file to read sequences.\n Output will be a set of CSV files containing information about the length of the reads, number of reads aligned to a reference sequence and the length of overhangs of the alignment. \n Reference RNA file will be automatically indexed', add_help=False)

optional = parser.add_argument_group('Optional Arguments')
required = parser.add_argument_group('Required Arguments')
flags = parser.add_argument_group('Flags')

#Add back help...
optional.add_argument(
    '-h',
    '--help',
    action='help',
    default=SUPPRESS,
    help='show this help message and exit'
)

required.add_argument('-r', '--reference', help='Path to the reference seqeunces', required=True)
required.add_argument('-q', '--reads', help='Path to the read sequences', required=True)
optional.add_argument('-n', '--name',  help='Prefix for the output files')
optional.add_argument('-d', '--directory', default = os.curdir, help='Directory to store the output files')
optional.add_argument('-m', '--min_score', default=-1, type=int, help='Minimum score to accept, default is the shortest read length')
flags.add_argument('-e', '--remove_exact', action='store_true', help='Remove exact read matches to the reference sequence')
flags.add_argument('-u', '--make_unique', action='store_true', help='Make FASTA headers unique in reference and reads i.e. >Read_1 >Read_2')

parser._action_groups.append(optional)
parser._action_groups.append(flags)

args = parser.parse_args()


