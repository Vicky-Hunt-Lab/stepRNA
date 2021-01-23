#!/usr/bin/env python3

from Bio import SeqIO
import os
from subprocess import run

from argparse import ArgumentParser

parser = ArgumentParser(description='Convert fasta headers from non-unique headers to read[x]')

parser.add_argument('-f', '--fasta', help='Fasta files to make unique headers')
parser.add_argument('-n', '--name', default = 'Read', help='Basename for header')

args = parser.parse_args()

fasta_in = args.fasta
name = args.name

fasta_out = os.path.splitext(fasta_in)[0] + '_corrected' + os.path.splitext(fasta_in)[1]

with open(fasta_in) as fin, open(fasta_out, 'w') as fout:
    records = SeqIO.parse(fin, 'fasta')
    for x, record in enumerate(records):
        record.id = name + '_{}'.format(x + 1) 
        record.description = ''
        SeqIO.write(record, fout, 'fasta')

run(['mv', fasta_out, fasta_in])
