#!/usr/bin/env python3

import os
from subprocess import run, PIPE

from stepRNA.processing import sam_to_bam
from stepRNA.general import mini_maxi, replace_ext, check_dir
from stepRNA.output import Logger


def main(ref_base, reads, prefix, min_score, logger):
    '''Run the main Bowtie2 command'''
    minimum, maximum = mini_maxi(reads, file_type = 'fasta')
    if min_score != -1:
        min_score = 3 * min_score
    else:
        min_score = 3 * minimum
    # Run bowtie command...
    sam_file = replace_ext(prefix, '.sam')
    command = ['bowtie2', '-x', ref_base, '-U', reads, '-f', '-N', '0', '-L', '10', '--no-1mm-upfront', '--nofw','--local', '--ma', '3', '--mp', '{},{}'.format(maximum, maximum), '--score-min', 'L,{},0'.format(min_score), '-S', sam_file]
    bowtie = run(command, stderr=PIPE)
    logger.write('Alignment statistics (from Bowtie2):\n')
    logger.write(bowtie.stderr.decode('utf-8'))
    logger.write('Alignment complete\n')
    #Convert sam to bam...
    sorted_bam = sam_to_bam(sam_file, logger)
    return sorted_bam

if __name__ == "__main__":
    description = 'Run a bowtie command with customisable features and convert to a sorted BAM file'

    from argparse import ArgumentParser, SUPPRESS

    parser = ArgumentParser(description=description, add_help=False)

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
    required.add_argument('--ref_base', '-r', help='Bowtie-build reference basename', required=True)
    required.add_argument('--reads', '-q', help='Reads to align to the reference', required=True)
    #Add optional arugments...
    optional.add_argument('-m', '--min_score', default=-1, type=int, help='Minimum score to accept, default is the shortest read length')
    optional.add_argument('-n', '--name',  help='Prefix for the output files')
    optional.add_argument('-d', '--directory', default = os.curdir, help='Directory to store the output files')

    args = parser.parse_args()

    ref_base = args.ref_base
    reads = args.reads
    min_score = args.min_score
    outdir = check_dir(args.directory)
    if args.name is None:
        filename = os.path.splitext(reads)[0]
    else:
        filename = args.name
    #Join together output directory and filename to make a prefix...
    prefix = os.path.join(outdir, filename)
    logger = Logger(prefix + '_output.log')
    #Run main script
    main(ref_base, reads, prefix, min_score, logger)
