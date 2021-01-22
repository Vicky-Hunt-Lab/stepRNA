#!/usr/bin/env python3

from Bio import SeqIO
import numpy as np
from subprocess import run, PIPE
import os

from argparse import ArgumentParser


parser = ArgumentParser(description='Align 26G RNAs to passenger sequences')

parser.add_argument('-r', '--reference', help='Path to the 26G indexed reference basename')
parser.add_argument('-p', '--p_reads', help='Path to the passenger read sequences')
parser.add_argument('-m', '--min_score', default=True, help='Minimum score to accept, default is the shortest read length')

args = parser.parse_args()

# Script to take a 26G reference index and passenger sequence reads (15-28 nt) and align them 


def shortest_seq(read_file, file_type='fasta'):
    '''Find the shortest length sequene from a fasta or fastq file'''
    with open(read_file, 'rU') as handle:
        parser = SeqIO.parse(handle, 'fasta')
        minimum = np.inf
        for record in parser:
            if len(record.seq) < minimum:
                minimum = len(record.seq)
    return minimum

ref = args.reference
reads = args.p_reads
min_score = args.min_score

if min_score:
    min_score = 3 * shortest_seq(reads, file_type = 'fasta')
else:
    min_score = 3 * min_score



def replace_ext(path, extension):
    ext_path = os.path.basename(path) 
    ext_path = os.path.splitext(ext_path)[0]
    return ext_path + extension

sam_file = replace_ext(reads, '.sam')

command = ['bowtie2', '-x', ref, '-U', reads, '-f', '-N', '0', '-L', '10', '--no-1mm-upfront', '--local', '--ma', '3', '--mp', '28,28', '--rfg', '28,28', '--rdg', '28,28', '--score-min', 'L,{},0'.format(min_score), '-S', sam_file]
print(command)

bowtie = run(command)


