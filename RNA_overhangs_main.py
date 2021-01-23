#!/usr/bin/env python3

from subprocess import run
import os
import glob

# Add arguments
# Reference, reads, convert_headers?, min_length (of matches)

from argparse import ArgumentParser

parser = ArgumentParser(description='Script to run the whole pipeline')

parser.add_argument('-r', '--reference', help='Path to the 26G reference basename')
parser.add_argument('-f', '--reads', help='Path to the Passenger reads (FASTA file)')
parser.add_argument('-c', '--conv_headers', nargs='?', const=True, help='Convert headers to unique names?')
parser.add_argument('-m', '--min_length', default=True, help='Set a minimum length for the exact alignments, default is the shorted read length in the Passenger reads')

args = parser.parse_args()

reference = args.reference
reads = args.reads
min_length = args.min_length

# Create unique headers if needed...
if args.conv_headers:
    for fin in [reads, reference]:
        run(['uniq_fasta_head.py', '-f', fin])

# Reverse the direction of the passenger sequences... (check with Becky)

def replace_ext(path, extension):
    ext_path = os.path.basename(path) 
    ext_path = os.path.splitext(ext_path)[0]
    return ext_path + extension

run(['reverse_seq.py','-r', reads])
reads = replace_ext(reads, '_reversed.fa')

# Search for indexed files, if not present then generate them...
def bowtie_indexed(reference):
    ref_name = os.path.splitext(reference)[0]
    directory = reference.split('/')[:-1]
    if directory == []:
        directory = './'
    search = glob.glob('{}*.bt2'.format(ref_name.split('/')[-1]))
    count = 0
    for name in search:
        if name in os.listdir(directory):
            count += 1
    if count == 6:
        return True

# Index if no basename
if bowtie_indexed(reference):
    ref_name = os.path.splitext(reference)[0]
else:
    ref_name = os.path.splitext(reference)[0]
    run(['bowtie2-build', reference, ref_name]) 

# Do alignment
command = ['DICER_RNA_Alignment.py', '-r', ref_name, '-p', reads]

if min_length:
    pass
else:
    command = command + ['-m', min_length]

run(command)

