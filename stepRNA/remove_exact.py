#Modules:
from stepRNA.processing import rm_ref_matches


def main(ref, reads):
    """Remove the exact matches to the reference seqeunces from the read seqeuences
    
    REF [STR] - path to the reference seqeunces (FASTA file)
    READS [STR] - path to the read sequences (FASTA file)
    
    Returns: Path to the read sequences"""
    return rm_ref_matches(ref, reads)

if __name__ == "__main__":
    #ArgumentParser
    from argparse import ArgumentParser, SUPPRESS

    parser = ArgumentParser(description='Remove the exact matches to the reference seqeunces from the read seqeuences', add_help=False)

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
    required.add_argument('-r', '--reference', help='Path to the reference seqeunces')
    required.add_argument('-q', '--reads', help='Path to the read sequences')

    parser._action_groups.append(optional)

    args = parser.parse_args()

    ref = args.reference
    reads = args.reads
    main(ref, reads)
