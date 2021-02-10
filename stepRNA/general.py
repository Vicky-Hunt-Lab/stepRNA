#General functions to use

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


