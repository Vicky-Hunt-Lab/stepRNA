#Modules
from stepRNA.processing import make_unique


def main(fin, filetype='filetype', name='Read', keep_ori=False):
    return make_unique(fin, filetype, name, keep_ori)

if __name__ == "__main__":
    from argparse import ArgumentParser, SUPPRESS

    parser = ArgumentParser(description='Make unique FASTA/Q headers in a file', add_help=False)

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
    required.add_argument('--fin', '-f', help='FASTA file')
    #Add optional arugments...
    optional.add_argument('--filetype', '-t', default='fasta', choices = ['fasta', 'fastq'], help='Filetype in')
    optional.add_argument('--name', '-n', default='Read', type=str, help='Header prefix. Default: "Read"')
    flags.add_argument('--keep_ori', '-o', action='store_true', help='Keep the origional file or replace. Default is to replace')

    #parser._action_groups.append(optional)
    args = parser.parse_args()
    main(fin, filetype, name, keep_ori)
