from subprocess import run, PIPE
import os


def main(ref):
    ref_base = os.path.splitext(ref)[0]
    command = ['bowtie2-build', ref, ref_base]
    run(command, stdout=PIPE)
    return ref_base

if __name__ == "__main__":
    description = 'Index a reference sequence for ready for bowtie alignment'

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
    required.add_argument('--ref', '-r', help='Bowtie-build reference basename')

    #parser._action_groups.append(optional)
    args = parser.parse_args()

    ref = args.ref
    main(ref)
