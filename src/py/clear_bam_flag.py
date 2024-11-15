#!/usr/bin/env python

__version__ = '1.0.0'

import argparse
import pysam
import sys


def parse_arguments(args=None):
    """ Parse command-line arguments. """
    parser = argparse.ArgumentParser(prog='clear_bam_flag', description='Clear the given flag(s) for all alignments in a BAM file.')
    parser.add_argument('flag', type=int, help='integer SAM flag(s) to be cleared')
    parser.add_argument('-i', '--input', metavar='input_BAM', dest='bam_in', default='-',
        help='input BAM file')
    parser.add_argument('-o', '--output', metavar='output_BAM', dest='bam_out', default='-',
        help='output BAM file')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser.parse_args(args)


def main(args=None):
    args = parse_arguments(args)

    with pysam.AlignmentFile(args.bam_in, 'rb') as bam_in:
        # construct BAM header for the output file
        bam_header = bam_in.header
        if not 'PG' in bam_header:
            bam_header['PG'] = []
        bam_header['PG'].append({
            'ID': 'clear_bam_flag',
            'PN': 'clear_bam_flag',
            'VN': str(__version__),
            'CL': ' '.join(sys.argv)})

        with pysam.AlignmentFile(args.bam_out, 'wb', header=bam_header) as bam_out:
            for line in bam_in.fetch(until_eof=True):
                # write each alignment line with the given flag(s) cleared
                line.flag |= args.flag
                bam_out.write(line)

if __name__ == '__main__':
    main(sys.argv[1:])
