#!/usr/bin/env python3

#Python default modules
import numpy as np
from subprocess import run, PIPE
import os
import sys

#Modules that need to be installed
try:
    from Bio import SeqIO
except ImportError:
    print('Error: Biopython not found, can be installed with\npip3 install biopython', file=sys.stderr)
    exit()

try:
    import pysam
except ImportError:
    print('Error: Pysam not found, can be installed with\npip3 install pysam', file=sys.stderr)
    exit()

#Set-up arguments...
from argparse import ArgumentParser

parser = ArgumentParser(description='Align an reference RNA file to passenger sequences\n Reference RNA file will be automatically referenced and unique headers will be made if required')

parser.add_argument('-r', '--reference', help='Path to the 26G indexed reference basename')
parser.add_argument('-p', '--p_reads', help='Path to the passenger read sequences')
parser.add_argument('-m', '--min_score', default=-1, type=int, help='Minimum score to accept, default is the shortest read length')
parser.add_argument('-d', '--remove_exact', action='store_true', help='Remove exact read matches to the reference sequence')
parser.add_argument('-u', '--make_unique', action='store_true', help='Make FASTA headers unique in reference and reads i.e. >Read_1 >Read_2')

args = parser.parse_args()

#Functions to use...
def shortest_seq(read_file, file_type='fasta'):
    '''Find the shortest length sequene from a fasta or fastq file'''
    with open(read_file, 'rU') as handle:
        parser = SeqIO.parse(handle, 'fasta')
        minimum = np.inf
        for record in parser:
            if len(record.seq) < minimum:
                minimum = len(record.seq)
    return minimum

def replace_ext(path, extension):
    '''Replace the file extension with one of choice'''
    ext_path = os.path.basename(path) 
    ext_path = os.path.splitext(ext_path)[0]
    return ext_path + extension

def sam_to_bam(sam_file):
    '''Convert the sam file output to a sorted bam'''
    sorted_bam = replace_ext(sam_file, '.sorted.bam')
    #Initialise files
    infile = pysam.AlignmentFile(sam_file, 'r')
    outfile = pysam.AlignmentFile(sorted_bam, 'wb', template=infile)
    infile.close()
    outfile.close()
    #Convert and sort
    try:
        pysam.view('-bS', '-o', sorted_bam, sam_file, catch_stdout=False)
        pysam.sort('-o', sorted_bam, sorted_bam)
        pysam.index(sorted_bam)
        os.remove(sam_file)
    except:
        print('Something went wrong when converting the SAM file to BAM file. Exiting...')
        sys.exit(1)
    return sorted_bam

#check_unique() is very slow...
def check_unqiue(seq_file, filetype='fasta'):
    '''Check a file for unique FASTA/Q IDs.
    Returns False if all IDs are unique'''
    lst = []
    for record in SeqIO.parse(seq_file, filetype):
        print(record.id)
        if record.id in lst:
            return True
        else:
            lst.append(record.id)
    return False

def make_unique(seq_file, filetype='fasta', name='Read', keep_ori=False):
    '''Make unqiue file names for each file in a FASTA/Q file'''
    temp_file = replace_ext(seq_file, '_corrected.fasta')
    with open(temp_file, 'w') as temp:
        for x, record in enumerate(SeqIO.parse(seq_file, filetype)):
            record.id = name + '_{}'.format(x + 1) 
            record.description = ''
            SeqIO.write(record, temp, filetype)
    if keep_ori:
       return temp_file 
    else:
        run(['mv', temp_file, seq_file])
        return seq_file

def rm_ref_matches(refs, reads, ref_type='fasta', read_type='fasta'):
    '''Search through the reference genome and reads and remove any exact matches
    Saves a new file named BASE_rmref.fasta'''
    rm_ref = replace_ext(reads, '_rmref.fasta')
    print('Writing to {}'.format(rm_ref))
    fout = open(rm_ref, 'w')
    ref_records = SeqIO.parse(refs, ref_type)
    read_records = SeqIO.parse(reads, read_type)
    for read_record in read_records:
        keep = True
        count = 0
        for ref_record in ref_records:
            count += 1
            if read_record.seq == ref_record.seq:
                keep = False
                break
        if keep:
            SeqIO.write(read_record, fout, read_type)
        ref_records = SeqIO.parse(refs, read_type)
    return rm_ref

#Command to run...
# Parse arguments...
ref = args.reference
reads = args.p_reads
min_score = args.min_score

#Remove exact matches to reference if set...
if args.remove_exact:
    reads = rm_ref_matches(ref, reads)

#Make unique headers if set...
if args.make_unique:
    reads = make_unique(reads)
    ref = make_unique(reference)

#Build a reference (suppress verbosity)...
command = 'bowtie2-build {} {}'.format(reads, ref)
command = command.split()
run(command, stderr=PIPE)

# Set min_score if not present, else set as match bonus * min_score...
if min_score != -1:
    min_score = 3 * min_score
else:
    min_score = 3 * shortest_seq(reads, file_type = 'fasta')

# Run bowtie command...
sam_file = replace_ext(reads, '.sam')
command = ['bowtie2', '-x', ref, '-U', reads, '-f', '-N', '0', '-L', '10', '--no-1mm-upfront', '--local', '--ma', '3', '--mp', '28,28', '--score-min', 'L,{},0'.format(min_score), '-S', sam_file]
bowtie = run(command)
#Convert sam to bam...
sorted_bam = sam_to_bam(sam_file)

#Count overhangs...

#Put those that are counted into a new BAM file...

#Put overhangs infomation into a csv and print to terminal...

#Print histogram of overhangs to terminal...

#Summarise type of referernce:read. For example:
# 9 possibilities...
# ------------
#   -------     Is one 'type'
