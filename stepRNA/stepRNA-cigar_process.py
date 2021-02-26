#!/usr/bin/env python3

from collections import defaultdict
import json
import os

from stepRNA.processing import MakeBam
from stepRNA.commands import left_overhang, right_overhang
from stepRNA.general import check_dir
from stepRNA.output import refs_counts

import pysam


def main(sorted_bam, filepath, write_json=False):
    '''Process CIGAR strings.

    sorted_bam [STR] - path to a sorted BAM file
    
    Outputs: JSON dictionaries and BAM files for each overhang length'''
    samfile = pysam.AlignmentFile(sorted_bam, 'rb')
    right_dic = defaultdict(lambda:0)
    left_dic = defaultdict(lambda:0)
    type_dic = defaultdict(lambda:0)
    read_len_dic = defaultdict(lambda:0)
    refs_read_dic = defaultdict(lambda:0)
    for name in samfile.references:
        refs_read_dic[name] = 0
    MakeBam_dic = defaultdict(lambda:None) 
    all_passed = MakeBam(samfile)
    for line in samfile:
            if line.cigarstring != None:
                if ('D' or 'I') not in line.cigarstring:
                    ref_pos = line.get_reference_positions(full_length = True)
                    try:
                        right, right_type = right_overhang(samfile, line, ref_pos)
                        left, left_type = left_overhang(samfile, line, ref_pos)
                        #Add to MakeBam
                        def add_to_MakeBam(dic, length, additional, record):
                            length =  additional + '_' + str(length) 
                            if dic[length + '_overhang'] == None:
                                dic[length + '_overhang'] = MakeBam(samfile)
                            dic.get(length + '_overhang').add_record(record)
                        add_to_MakeBam(MakeBam_dic, right, right_type, line)
                        add_to_MakeBam(MakeBam_dic, left, left_type, line)
                        all_passed.add_record(line)
                        # Create dictionaries to sort information...
                        right_dic[right] += 1 # right overhang count
                        left_dic[left] += 1 # left overhang count
                        type_dic[left_type + '_' + right_type] += 1 # type of overhang count
                        read_len_dic[line.query_length] += 1 # read length count
                        refs_read_dic[line.reference_name] += 1 # number of reads algining to reference
                    except Exception:
                        continue
    #directory = os.path.dirname(filepath)
    #outdir = os.path.join(directory, 'AlignmentFiles')
    outdir = filepath + '_AlignmentFiles'
    check_dir(outdir)
    for key in MakeBam_dic:
        outfile = os.path.join(outdir, '{}_{}.bam'.format(os.path.basename(filepath), key))
        MakeBam_dic[key].save_to_file(outfile)
    all_passed.save_to_file(os.path.join(outdir, '{}_passed.bam'.format(os.path.basename(filepath))))
    if write_json:
        for dic in right_dic, left_dic, type_dic, read_len_dic, refs_read_dic:
            #Needs to be sorted out!!!
            def write_to_json(dic, filepath, suffix):
                with open(prefix + suffix, 'w') as fout:
                        json.dump(dic, fout)
    return right_dic, left_dic, type_dic, read_len_dic, refs_read_dic

if __name__ == "__main__":
    from argparse import ArgumentParser, SUPPRESS

    parser = ArgumentParser(description='Process CIGAR strings to calculate the overhang lengths', add_help=False)

    optional = parser.add_argument_group('optional arguments')
    required = parser.add_argument_group('required arguments')
    flags = parser.add_argument_group('flags')
    #Add back help...
    optional.add_argument(
        '-h',
        '--help',
        action='help',
        default=SUPPRESS,
        help='show this help message and exit'
    )

    #Add required arguments...
    required.add_argument('--bamfile', '-b', help='Path to a sorted BAMfile', required=True)
    optional.add_argument('--prefix', '-p', help='Prefix to add to the file. Default is file basename')

    flags.add_argument('--write_json', '-j', action='store_true' , help='Write counts dictionaries to JSON files')

    args = parser.parse_args()

    sorted_bam = args.bamfile
    write_json = args.write_json
    if args.prefix is None:
        prefix = os.path.splitext(sorted_bam)[0]
    else:
        prefix = args.prefix
    right_dic, left_dic, type_dic, read_len_dic, refs_read_dic = main(sorted_bam, prefix, write_json)
