# File processing functions

# Python core modules
import os
import sys
from subprocess import run, PIPE

#Package modules
from stepRNA.general import replace_ext

try:
    import pysam
except ImportError:
    print('Error: Pysam not found, can be installed with\npip3 install pysam', file=sys.stderr)
    sys.exit(1)

def sam_to_bam(sam_file):
    '''Convert the sam file output to a sorted bam
    
    sam_file [STR] - path to the SAM alignment file
    
    This function will convert the SAM file to a sorted BAM and index it in the directory where the SAM file is.

    Returns: A sorted BAM file string. Replaces the SAM file with sorted BAM file'''
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

def make_unique(seq_file, filetype='fasta', name='Read', keep_ori=False):
    '''Make unqiue file names for each file in a FASTA/Q file.
    
    seq_file [STR] - path to the Fasta/Fastq file
    filetype [STR] - specify whether "fasta" or fastq"; default: "fasta"
    name [STR] - prefix for the header. Default: "Read"
    keep_ori [True/False] - keep the original file
    
    Will output the new file into the input file directory.

    Returns: New filepath string'''
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
    
    refs [STR] - path to the reference seqeunces
    reads [STR] - path to the read seqeunces
    ref/read_type [STR] - specify whether "fasta" or "fastq" for the reference or read files

    Saves a new file named READ_BASE_rmref.fasta
    
    Returns: New filepath string'''
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
