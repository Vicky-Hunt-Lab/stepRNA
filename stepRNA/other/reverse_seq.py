#!/usr/bin/env python3

from Bio import SeqIO
import os

from argparse import ArgumentParser

parser = ArgumentParser(description='Reverse passenger sequences to set correct direction for alignment')

parser.add_argument('-r', '--reads', help='Path to the reads to reverse')

args = parser.parse_args()

fasta_in = args.reads

fasta_out = os.path.splitext(fasta_in)[0] + '_reversed' + os.path.splitext(fasta_in)[1]

with open(fasta_in) as fin, open(fasta_out, 'w') as fout:
    records = SeqIO.parse(fin, 'fasta')
    for record in records:
        record.seq = record.seq[::-1]
        SeqIO.write(record, fout, 'fasta')

