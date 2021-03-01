# Python core modules
from collections import defaultdict
import os
import sys
from subprocess import run, PIPE

#Installed modules
import pysam
from Bio import SeqIO

#Package modules
from stepRNA.general import replace_ext

def sam_to_bam(sam_file, logger):
    '''Convert the sam file output to a sorted bam
    
    sam_file [STR] - path to the SAM alignment file
    
    This function will convert the SAM file to a sorted BAM and index it in the directory where the SAM file is.

    Returns: A sorted BAM file string. Replaces the SAM file with sorted BAM file'''
    logger.write('Converting from SAM FILE to SORTED BAM FILE')
    sorted_bam = replace_ext(sam_file, '.sorted.bam')
    infile = pysam.AlignmentFile(sam_file, 'r')
    outfile = pysam.AlignmentFile(sorted_bam, 'wb', template=infile)
    infile.close()
    outfile.close()
    try:
        pysam.view('-bS', '-o', sorted_bam, sam_file, catch_stdout=False)
        pysam.sort('-o', sorted_bam, sorted_bam)
        pysam.index(sorted_bam)
        os.remove(sam_file)
        logger.write('SAM to BAM conversion complete')
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
        print('Unique headers made')
        return temp_file 
    else:
        run(['mv', temp_file, seq_file])
        print('Unique headers made')
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
    with open(rm_ref, 'w') as fout:
        ref_records = SeqIO.parse(refs, ref_type)
        read_records = SeqIO.parse(reads, read_type)
        ref_seqs = [] 
        for ref_record in ref_records:
            ref_seqs.append(ref_record.seq)
        ref_seqs = sorted(ref_seqs)
        read_seqs = []
        for read_record in read_records:
            read_seqs.append([read_record.seq, read_record])
        read_seqs = sorted(read_seqs, key = lambda x: x[0])
        writetofasta = []
        while len(ref_seqs) != 0:
            ref_seq = ref_seqs.pop(0)
            output = []
            index = 0
            index_counter = 0
            for read_seq in read_seqs:
                index_counter += 1
                if read_seq[0].tomutable() == ref_seq.tomutable():
                    index = index_counter
                else:
                    output.append(read_seq[1])
                    if index == 0:
                        continue
                    else:
                        read_seqs = read_seqs[index:]
                        SeqIO.write(output, fout, read_type)
                        break
        SeqIO.write(output, fout, read_type)
        return rm_ref 

class MakeBam():
    '''Initialise with an open pysam AlignmentFile and remove the SQ information.
    
    Add pysam.AlignmentFile records with add_record
    Use save_to_file to save to filename'''
    def __init__(self, samfile):
        self.header_dic = {}
        for key in samfile.header.keys():
            if key == 'SQ':
                self.header_dic[key] = []
            else:
                self.header_dic[key] = samfile.header.get(key)
        self.name_lst = []
        self.length_lst = []
        self.records = []
        self.ind = 0
    def add_record(self, line):
        '''Add a SAM file record to the header dictionary and file'''
        if line.reference_name in self.name_lst:
            self.records.append([line, self.name_lst.index(line.reference_name)])
        else: 
            #Retains index information
            self.records.append([line, self.ind])
            self.name_lst.append(line.reference_name)
            self.length_lst.append(line.header.get_reference_length(line.reference_name))
            self.ind += 1

    def save_to_file(self, filename, filetype = 'bam'):
        '''Save the information to a file: SAM or BAM. Default = bam'''
        #Forms the new BAM header
        for name, length in zip(self.name_lst, self.length_lst):
            self.header_dic['SQ'].append({'SN': name, 'LN' : length})
        #Potential for future different file type creation
        if filetype == 'bam':
            with pysam.AlignmentFile(filename, 'wb', header=self.header_dic) as outfile:
                for line, ind in self.records:
                    line.reference_id = ind 
                    outfile.write(line)
