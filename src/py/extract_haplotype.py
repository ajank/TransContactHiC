#!/usr/bin/env python

__version__ = '1.0.0'

import argparse
import gzip
import pysam
import sys


def parse_arguments(args=None):
    """ Parse command-line arguments. """
    parser = argparse.ArgumentParser(prog='extract_haplotype', description='Extract same-haplotype Hi-C reads from a BAM file. Useful for obtaining haplotype-specific contact maps. Alignments will be annotated in the CO field according to the underlying SNVs, and all alignments for a sequencing read will get an aggregated annotation in the RG field.')
    parser.add_argument('-v', '--variants',  metavar='TSV.gz', dest='tsv', required=True,
        help='tab-separated .tsv.gz file containing SNVs, with header row and the folowing columns: chrom, 1-based position, haplotype 1, haplotype 2')
    parser.add_argument('-i', '--input', metavar='input_BAM', dest='bam_in', required=True,
        help='input BAM file')
    parser.add_argument('-o', '--output', metavar='output_BAM', dest='bam_out', required=True,
        help='output BAM file')
    parser.add_argument('-r', '--requested-haplotype', dest='haplotype', required=False, help='requested haplotype')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    return parser.parse_args(args)


def read_tsv(file):
    """ Read SNV annotations in a custom tab-separated format. """
    haplotypes = []

    for line in gzip.open(file, mode='rt'):
        line = line.strip('\n')
        if not line or line.startswith('#'):
            continue

        line = line.split()
        if not haplotypes:
            haplotypes = line[2:]
            continue

        nucleotides = line[2:]
        if any(nucleotides.count(nucl) > 1 for nucl in nucleotides):
            warnings.warn('Warning: excluding SNVs that do not determine a haplotype!')
        snv_dict = {nucl: hapl for nucl, hapl in zip(nucleotides, haplotypes) if nucleotides.count(nucl) == 1}

        # line[0] is chrom, line[1] is 1-based pos
        yield (line[0], int(line[1]) - 1), snv_dict


def annotate_haplotype(alignment_buffer, snvs, requested_haplotype = None):
    """ Annotate haplotypes for one sequencing read, according to the provided SNVs. """
    haplotypes = set()

    # find the underlying SNVs, save the annotation in the CO field
    for line in alignment_buffer:
        comment = []

        for query_pos, reference_pos in line.get_aligned_pairs(matches_only=True):
            snv = snvs.get((line.reference_name, reference_pos))
            if snv:
                haplotype = snv.get(line.query_sequence[query_pos])
                if haplotype:
                    haplotypes.add(haplotype)
                    # annotation: haplotype, reference position, query position, line 1/2
                    comment.append(haplotype + '|' + str(reference_pos) + '|' + str(query_pos) + '|'
                                   + '|'.join(snv.keys()))

        if comment:
            line.set_tag('CO', ','.join(comment), replace=False)

    # determine the haplotype annotation for all alignments in the sequencing read
    if len(haplotypes) == 0:
        read_group = 'unknown'
    elif len(haplotypes) == 1:
        read_group = list(haplotypes)[0]
    else:
        read_group = 'ambiguous'

    # if a specific haplotype is requested, discard all the alignments if it is not matched
    if requested_haplotype:
        if read_group != requested_haplotype:
            alignment_buffer.clear()

    # save the haplotype annotation in the RG field for all alignments
    for line in alignment_buffer:
        line.set_tag('RG', read_group)


def write_buffer(alignment_buffer, output_file):
    """ Write the alignment lines for one sequencing read to the output BAM file. """
    for line in alignment_buffer:
        output_file.write(line)


def main(args=None):
    args = parse_arguments(args)

    print('Parsing SNVs...')
    snvs = dict(read_tsv(args.tsv))

    print('Parsing BAM...')
    with pysam.AlignmentFile(args.bam_in, 'rb') as bam_in:
        # construct BAM header for the output file
        bam_header = bam_in.header
        if not 'PG' in bam_header:
            bam_header['PG'] = []
        bam_header['PG'].append({
            'ID': 'extract_haplotype',
            'PN': 'extract_haplotype',
            'VN': str(__version__),
            'CL': ' '.join(sys.argv)})

        with pysam.AlignmentFile(args.bam_out, 'wb', header=bam_header) as bam_out:
            buffer = []

            for line in bam_in.fetch(until_eof=True):
                # all the alignment lines for one sequencing read are in the buffer, we need to flush them first
                if buffer and buffer[0].query_name != line.query_name:
                    annotate_haplotype(buffer, snvs, args.haplotype)
                    write_buffer(buffer, bam_out)
                    buffer.clear()

                buffer.append(line)

            # flush the alignment lines for the last sequencing read
            if buffer:
                annotate_haplotype(buffer, snvs)
                write_buffer(buffer, bam_out)


if __name__ == '__main__':
    main(sys.argv[1:])
