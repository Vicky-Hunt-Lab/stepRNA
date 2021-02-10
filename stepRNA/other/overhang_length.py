#!/usr/bin/env/ python3

from subprocess import run, PIPE
import os
import pysam

from argparse import ArgumentParser

parser = ArgumentParser(description='Extract CIGAR information about the length of overhangs')

parser.add_argument('-b', '--bamfile', help='Path to the sorted BAM file')

args = parser.parse_args()

samfile = args.samfile

# Important columns in sam file are:
# Qname [0], Rnam[2], Pos[3], Cigar[5], Seq[9]
# Qname is the query name, might be helpful if they were differenet names
# Rname is the reference name, helpful as above
# Pos is the leftmost starting position:
    # If a  read goes overlaps the 5' end then this will be 1, otherwise it shows:
    # 5'______________________ 3'
    #       3'______________ 5'
# Cigar shows the substituations (overhang in this case) and matches along with lengths. E.g. 20M1S is a rightmost 1nt overhang; 3S20M is a leftmost 3nt overhang (makes POS redundant??)
# Seq is the aligned portion, should be equal to sum of Cigar string

align_file = pysam.AlignmentFile(bamfile, 'rb')
for line in align_file:
    if line.cigarstring != None:
        test = line.qend - line.qstart - len(line.seq)
        if test < 0:
            print('##########')
            print(test)
            print(line.cigarstring)
        print(test)
        if line.qstart == 3:
            print(line.seq, 'is {}bp long'.format(len(line.seq)))
            print(line.cigarstring, '=', line.cigartuples)
            print(line.cigar)
            print(line.qstart)
        

def left_match(line):
    if line.reference_end == line.query_length - line.query_alignment_start:
        
align_file = pysam.AlignmentFile(bamfile, 'rb')
for line in align_file:
    if line.cigarstring != None:
        if line.reference_start == line.reference_alignment_start:

# Need to sort out fasta headers first to de-duplicate them...
