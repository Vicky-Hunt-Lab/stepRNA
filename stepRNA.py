#!/usr/bin/env python3

#Python default modules
import numpy as np
from subprocess import run, PIPE
import os
import sys
from collections import defaultdict
import csv

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
from argparse import ArgumentParser

parser = ArgumentParser(description='Align an reference RNA file to passenger sequences\n Reference RNA file will be automatically referenced and unique headers will be made if required')

parser.add_argument('-r', '--reference', help='Path to the 26G indexed reference basename')
parser.add_argument('-q', '--reads', help='Path to the passenger read sequences')
parser.add_argument('-n', '--name',  help='Prefix for the output files')
parser.add_argument('-d', '--directory', default = os.curdir, help='Directory to store the output files')
parser.add_argument('-m', '--min_score', default=-1, type=int, help='Minimum score to accept, default is the shortest read length')
parser.add_argument('-e', '--remove_exact', action='store_true', help='Remove exact read matches to the reference sequence')
parser.add_argument('-u', '--make_unique', action='store_true', help='Make FASTA headers unique in reference and reads i.e. >Read_1 >Read_2')

args = parser.parse_args()

#Functions to use...
def check_dir(path):
    '''Check whether the directory specified is present. Create one if not.'''
    if os.path.isdir(path):
        print('Files being outputted to: {}'.format(os.path.abspath(path)))
    else:
        try:
            os.mkdir(os.path.relpath(path))
            print('Created directory: {}'.format(os.path.abspath(path)))
        except FileNotFoundError:
            print('Something went wrong when making the directory\n \
Are you sure it is a valid path entered?')
            sys.exit(1)
    return os.path.abspath(path)
            
def mini_maxi(read_file, file_type='fasta'):
    '''Find the shortest length sequene from a fasta or fastq file'''
    with open(read_file, 'rU') as handle:
        parser = SeqIO.parse(handle, 'fasta')
        minimum = np.inf
        maximum = 0
        for record in parser:
            length = len(record.seq)
            if length < minimum:
                minimum = len(record.seq)
            if length > maximum:
                maximum = len(record.seq)
    return minimum, maximum

def replace_ext(path, extension):
    '''Replace the file extension with one of choice'''
    ext_path = os.path.splitext(path)[0]
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

def left_overhang(sorted_bam, line, ref_positions):
    '''Get the length of the left overhang; 0 = no overhang, -ve = reference overhang, +ive = query overhang'''
    ref_positions = line.get_reference_positions(full_length=True) 
    if ref_positions[0] == 0:
        return 0, 'LE'
    elif ref_positions[0] == None:
        if line.reference_start != 0:
            raise Exception
        else:
            return line.query_alignment_start, 'LQ'
    else:
        #If reference_position[0] > 0
        return -ref_positions[0], 'LR'

def right_overhang(sorted_bam, line, ref_positions):
    '''Get the length of the right overhang; 0 = no overhang, -ve = reference overhang, +ive = query overhang'''
    ref_length = sorted_bam.lengths[sorted_bam.get_tid(line.reference_name)]
    if ref_positions[-1] == ref_length - 1:
        return 0, 'RE'
    elif ref_positions[-1] == None:
        if line.reference_end != ref_length:
           raise Exception 
        else:
            return ref_length - line.query_alignment_end, 'RQ'
    else:
        #If reference_position[-1] < ref_length
        return ref_positions[-1] - (ref_length - 1), 'RR'

def write_to_bam(line, left_type, right_type, prefix):
    '''Take a pysam.AlignmentFile record and write it to a file'''
    with open(prefix + '_' + left_type + '_' + right_type + '.bam', 'a+') as bam_out:
        bam_out.write(line.to_string() + '\n')

def make_csv(dics, csv_name, headers, show=True):
    '''Make a csv continaing the overhang information
    dics=[right_dic, left_dic] format'''
    keys = set()
    for dic in dics:
        for key in dic:
            keys.add(key)
    keys = list(keys)
    keys.sort()
    with open(csv_name, 'w') as csv_out:
        writer = csv.writer(csv_out, delimiter=',')
        writer.writerow(headers)
        if show:
            print('\t'.join(headers))
        for key in keys:
            if dics[0][key] == None:
                dics[0][key] = 0
            if dics[1][key] == None:
                dics[1][key] = 0
            writer.writerow([key, dics[1].get(key), dics[0].get(key)])
            if show:
                print('{}\t{}\t{}'.format(key, dics[1].get(key), dics[0].get(key)))

def make_type_csv(dic, csv_name, headers, show=True, sort=False):
    if sort:
        keys = set()
        for key in dic:
            keys.add(key)
        keys = list(keys)
        keys.sort()
    else:
        keys=dic.keys()
    with open(csv_name, 'w') as csv_out:
        writer = csv.writer(csv_out, delimiter = ',')
        writer.writerow(headers)
        if show:
            print('\t'.join(headers))
        for key in keys:
            writer.writerow([key, dic[key]])
            if show:
                print('{}\t{}'.format(key, dic[key]))

def print_hist(density_list, keys):
    for item in range(len(density_list)):
        try:
            length = int(density_list[item]) * '.'
            print('{}\t| {}'.format(keys[item], length))
        except ValueError:
            print('{}\t{}'.format(keys[item], ''))

#Command to run...
# Parse arguments...
ref = args.reference
reads = args.reads
min_score = args.min_score
outdir = check_dir(args.directory)
if args.name is None:
    filename = reads 
else:
    filename = args.name

#Join together output directory and filename to make a prefix...
prefix = os.path.join(outdir, filename)

#Remove exact matches to reference if set...
if args.remove_exact:
    reads = rm_ref_matches(ref, reads)

#Make unique headers if set...
if args.make_unique:
    reads = make_unique(reads)
    ref = make_unique(reference)

#Build a reference (suppress verbosity)...
ref_base = os.path.splitext(ref)[0]
command = 'bowtie2-build {} {}'.format(ref, ref_base)
command = command.split()
run(command, stdout=PIPE)

# Set min_score if not present, else set as match bonus * min_score...
minimum, maximum = mini_maxi(reads, file_type = 'fasta')

if min_score != -1:
    min_score = 3 * min_score
else:
    min_score = 3 * minimum 


# Run bowtie command...
sam_file = replace_ext(prefix, '.sam')
command = ['bowtie2', '-x', ref_base, '-U', reads, '-f', '-N', '0', '-L', '10', '--no-1mm-upfront', '--nofw','--local', '--ma', '3', '--mp', '{},{}'.format(maximum, maximum), '--score-min', 'L,{},0'.format(min_score), '-S', sam_file]

bowtie = run(command)
#Convert sam to bam...
sorted_bam = sam_to_bam(sam_file)

#Count overhangs...
bam_in = pysam.AlignmentFile(sorted_bam, 'rb')
right_dic = defaultdict(lambda:0)
left_dic = defaultdict(lambda:0)
type_dic = defaultdict(lambda:0)
read_len_dic = defaultdict(lambda:0)
refs_read_dic = defaultdict(lambda:0)
for name in bam_in.references:
    refs_read_dic[name] = 0

for line in bam_in:
    if line.cigarstring != None:
        if ('D' or 'I') not in line.cigarstring:
            ref_pos = line.get_reference_positions(full_length = True)
            try:
                right, right_type = right_overhang(bam_in, line, ref_pos)
                left, left_type = left_overhang(bam_in, line, ref_pos)
                # Create dictionaries to sort information...
                right_dic[right] += 1 # right overhang count
                left_dic[left] += 1 # left overhang count
                type_dic[left_type + '_' + right_type] += 1 # type of overhang count
                read_len_dic[line.query_length] += 1 # read length count
                refs_read_dic[line.reference_name] += 1 # number of reads algining to reference
                write_to_bam(line, left_type, right_type, prefix=prefix) # separate reads to 'bam' files
            except Exception:
                continue

#Put overhangs infomation into a csv and print to terminal...
print('\n## Overhang counts ##')
make_csv([right_dic, left_dic], prefix + '_overhang.csv', ['OH','Left','Right'])
print('\n## Overhang types ##')
make_type_csv(type_dic, prefix + '_type.csv', ['OH_type', 'count'])
print('\n## Read lengths ##')
make_type_csv(read_len_dic, prefix + '_read_len.csv', ['Read_length', 'count'], sort=True)
make_type_csv(refs_read_dic, prefix + '_ref_hits.csv', ['Reference', 'count'], show=False)
print()
with open(prefix + '_overhang_summary.csv') as summary:
    left_dens = []
    right_dens = []
    left_tot = 0
    right_tot = 0
    keys = []
    csv_reader = csv.reader(summary, delimiter=',')
    head = next(csv_reader)
    for line in csv_reader:
        keys.append(line[0])
        left_dens.append(int(line[1]))
        left_tot += int(line[1])
        right_dens.append(int(line[2]))
        right_tot += int(line[2])


for key in range(len(keys)):
    left_dens[key] = 100 * left_dens[key] / left_tot
    right_dens[key] = 100 * right_dens[key] / right_tot
    
#Print histogram of overhangs to terminal...
print('Left Handside Overhang Histogram')
print_hist(left_dens, keys)

print('Right Handside Overhang Histogram')
print_hist(right_dens, keys)
