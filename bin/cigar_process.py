from collections import defaultdict

from stepRNA.processing import MakeBam
from stepRNA.commands import left_overhang, right_overhang

import pysam

from argparse import ArgumentParser, SUPPRESS

parser = ArgumentParser(description='Process CIGAR strings to calculate the overhang lengths', add_help=False)

optional = parser.add_argument_group('optional arguments')
required = parser.add_argument_group('required arguments')
#Add back help...
optional.add_argument(
    '-h',
    '--help',
    action='help',
    default=SUPPRESS,
    help='show this help message and exit'
)

#Add required arguments...
required.add_argument('--bamfile', '-b', help='Path to a sorted BAMfile')

#optional.add_argument('--filetype', '-t', default='fasta', choices = ['fasta', 'fastq'], help='Filetype in')

args = parser.parse_args()

sorted_bam = args.bamfile

def main(sorted_bam):
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
    for line in samfile:
            if line.cigarstring != None:
                if ('D' or 'I') not in line.cigarstring:
                    ref_pos = line.get_reference_positions(full_length = True)
                    try:
                        right, right_type = right_overhang(samfile, line, ref_pos)
                        left, left_type = left_overhang(samfile, line, ref_pos)
                        #Add to MakeBam
                        def add_to_MakeBam(MakeBam_dic, right, additional):
                            right = str(right) + '_' + additional 
                            if MakeBam_dic[right + '_overhang'] == None:
                                MakeBam_dic[right + '_overhang'] = MakeBam(samfile)
                            MakeBam_dic.get(right + '_overhang').add_record(line)
                        add_to_MakeBam(MakeBam_dic, right, 'right')
                        add_to_MakeBam(MakeBam_dic, left, 'left')
                        # Create dictionaries to sort information...
                        right_dic[right] += 1 # right overhang count
                        left_dic[left] += 1 # left overhang count
                        type_dic[left_type + '_' + right_type] += 1 # type of overhang count
                        read_len_dic[line.query_length] += 1 # read length count
                        refs_read_dic[line.reference_name] += 1 # number of reads algining to reference
                        #write_to_bam(line, left_type, right_type, prefix=prefix) # separate reads to 'bam' files
                    except Exception:
                        continue
    for key in MakeBam_dic:
        print('test_data/' + key + '.bam')
        MakeBam_dic[key].save_to_file('test_data/' + key + '.bam')
        
if __name__ == "__main__":
    main(sorted_bam)
