from os import path

from Bio import SeqIO

def detect_extension(filepath):
    extension = path.splitext(filepath)[-1].strip('.')
    if extension is 'fa':
        return 'fasta'
    if extension is 'fq':
        return 'fastq'
    else:
        return extension

def open_parser(seqfile):
    filetype = detect_extension(seqfile)
    if not ((filetype == 'fasta') or (filetype == 'fastq')):
        raise ValueError('File extesion is not fasta or fastq, please change to the appropritate one using $mv FILENAME.EXT FILENAME.fasta')

    return SeqIO.parse(seqfile, filetype), filetype

def extract_seqs(seq_dict, record):
    seq = record.seq.back_transcribe()
    try:
        seq_dict[seq].append(record)
    except KeyError:
        seq_dict[seq] = [record]

def seq_matches(seq_dict1, seq_dict2, filetype, outfile):
    with open(outfile, 'w') as fout:
        for seq in seq_dict2.keys():
            if seq in seq_dict1:
                for r in seq_dict1[seq]:
                    SeqIO.write(r, fout, filetype)
#                for r in seq_dict2[seq]:
#                    SeqIO.write(r, fout, filetype)

def seq_nonmatches(seq_dict1, seq_dict2, filetype, outfile):
    with open(outfile, 'w') as fout:
        for seq in seq_dict1.keys():
            if seq not in seq_dict2:
                for r in seq_dict1[seq]:
                    SeqIO.write(r, fout, filetype)
#                for r in seq_dict2[seq]:
#                    SeqIO.write(r, fout, filetype)

def main(seqfile1, seqfile2, outfile, match = True):
    seqparser1, filetype = open_parser(seqfile1)
    seq_dict1 = {}
    for record in seqparser1:
        extract_seqs(seq_dict1, record)

    seqparser2, _ = open_parser(seqfile2)
    seq_dict2 = {}
    for record in seqparser2:
        extract_seqs(seq_dict2, record)

    if match:
        seq_matches(seq_dict1, seq_dict2, filetype, outfile) 
    else:
        seq_nonmatches(seq_dict1, seq_dict2, filetype, outfile) 

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description = 'Filter out non matches between two input files')

    parser.add_argument('sequence_file1', help = 'A sequence file input')
    parser.add_argument('sequence_file2', help = 'A sequence file input to compare to')
    parser.add_argument('outfile', help = 'Output filepath')
    parser.add_argument('--nonmatches', action='store_false', help = 'Extract non-matches')

    args = parser.parse_args()

    main(args.sequence_file1, args.sequence_file2, args.outfile, args.nonmatches)
