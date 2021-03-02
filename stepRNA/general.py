#General functions to use

#Python default modules
import numpy as np  
import os.path
import sys

try:
    from Bio import SeqIO
except ImportError:
    print('Error: Biopython not found, can be installed with\npip3 install biopython', file=sys.stderr)
    sys.exit(1)


def check_dir(path, show=True):
    '''Check whether the directory specified is present. Create one if not.'''
    if os.path.isdir(path):
        if show:
            print('Files being outputted to: {}'.format(os.path.abspath(path)))
    else:
        try:
            os.mkdir(os.path.relpath(path))
            if show:
                print('Created directory: {}'.format(os.path.abspath(path)))
        except FileNotFoundError:
            if show:
                print('Something went wrong when making the directory\n \
Are you sure it is a valid path entered?')
            sys.exit(1)
    return os.path.abspath(path)
            
def mini_maxi(read_file, file_type='fasta'):
    '''Find the shortest length sequence from a fasta or fastq file
    read_file [STR] - path to the Fasta/Fastq file
    file_type [STR] - specify filetype ["fasta"/"fastq"]; default is fasta'''
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
    '''Replace the extention of the string for a path with one of choice.
    path [STR] - path to the filename to replace
    extension [STR] - extenstion to add'''
    ext_path = os.path.splitext(path)[0]
    return ext_path + extension

import sys

def progressbar(it, prefix="", size=60, file=sys.stdout):
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()        
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()
