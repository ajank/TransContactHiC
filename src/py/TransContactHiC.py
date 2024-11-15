#!/usr/bin/env python

__version__ = '1.0.0'

import argparse
import gzip
import pysam
import sys
import warnings


"""
Overview of TransContactHiC:

The data on chromatin interactions are structured as follows: BAM alignment <
pairtools alignment (one or two BAM alignments) < query sequence (both read1
and read2).

The input files are: .bam file (multiple BAM alignments for each query
sequence), .pairs.gz file (multiple pairs of pairtools alignments for each
query sequence) and .tsv.gz file containing Single Nucleotide Variants.

The pairtools alignments are first matched with BAM alignments. Haplotypes are
then assigned to each of the pairtools alignments based on the read sequence
stored in BAM alignments and the provided SNVs.

The main output is .tsv.gz file containing one row for each suitable alignment
pair from the input .pairs.gz file, along with haplotype information. A .bam
file with the annotated BAM alignments can also be saved.
"""


def parse_region(region):
    """
    Parse a string specifying region and possibly haplotype of the
    interactions, e.g. chr1:10001-20000 or chr1:10001-20000:hapl1
    or ::hapl2.

    Returns: tuple (chrom, start, end, haplotype) where the unspecified
    elements equal None.
    """

    region_split = region.split(':')
    chrom = region_split[0] if region_split[0] else None
    haplotype = region_split[2] if len(region_split) > 2 and region_split[2] else None

    interval = region_split[1] if len(region_split) > 1 and region_split[1] else ''
    interval_split = interval.split('-')
    start = int(interval_split[0]) if interval_split[0] else None
    end = int(interval_split[1]) if len(interval_split) > 1 and interval_split[1] else None

    return chrom, start, end, haplotype


def region_matched(region, chrom, pos, haplotype):
    """
    Test if a given genomic position and haplotype falls within a
    given 'region' specified as an argument to --region.

    Returns: bool.
    """

    region_chrom, region_start, region_end, region_haplotype = region
    return ((not region_chrom or region_chrom == chrom) and
        (not region_start or region_start <= pos) and
        (not region_end or pos <= region_end) and
        (not region_haplotype or region_haplotype == haplotype))


def parse_arguments(args=None):
    """
    Parse the command-line arguments from 'args'.

    Returns: argparse.Namespace object.
    """

    parser = argparse.ArgumentParser(prog='TransContactHiC', description='Extract interactions from a Hi-C .bam file and the corresponding .pairs.gz file. Optionally annotate the haplotype of each alignment, according to the underlying SNVs.')
    parser.add_argument('-v', '--variants', metavar='variants.tsv.gz', required=True,
        help='tab-separated .tsv.gz file containing SNVs, with header row and the following columns: chrom, 1-based position, haplotype 1, haplotype 2')
    parser.add_argument('-b', '--input-bam', metavar='input.bam', required=True,
        help='input .bam file')
    parser.add_argument('-p', '--input-pairs', metavar='input.pairs.gz', required=True,
        help='input .pairs.gz file')
    parser.add_argument('-o', '--output-tsv', metavar='output.tsv.gz', required=True,
        help='output .tsv.gz file')
    parser.add_argument('--output-bam', metavar='output.bam', default=None,
        help='output .bam file to save the extracted interactions')
    parser.add_argument('--min-mapq', type=int, default=1, help='the minimal MAPQ score to consider an alignment; the same argument should given to pairtools while generating the .pairs.gz file. Default: 1')
    parser.add_argument('--region', action='append', default=[],
        help='region and possibly haplotype in which interactions shall take place, e.g. chr1:10001-20000 or chr1:10001-20000:hapl1 or ::hapl2; can be specified one or two times')
    parser.add_argument('--fully-phased', action='store_true', help='keep only fully phased interactions')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args(args)
    args.region = [parse_region(region) for region in args.region]

    return args


input_pairs_column_names = dict()
input_pairs_next_line = None

def read_buffer_pairtools(input_pairs):
    """
    Read all the alignment pairs for a single query sequence from the
    pairtools .pairs.gz file.

    Required columns in the .pairs.gz file: readID (must be the very first
    column), chrom1, pos51, pos31, strand1, chrom2, pos52, pos32, strand2.

    Recommended options for pairtools: pairtools parse --walks-policy all \
      --add-columns pos5,pos3 --drop-sam --no-flip. The .pairs.gz file shall
        not be sorted.

    Returns: tuple (read_id, list of split lines from the .pairs.gz file).
    """

    global input_pairs_column_names, input_pairs_next_line
    read_id = None
    buffer = []

    if not input_pairs_next_line:
        input_pairs_next_line = input_pairs.readline()

    while True:
        if not input_pairs_next_line:
            # end of input file, return the buffer collected so far
            break
        if input_pairs_next_line.startswith('#columns: '):
            # column headers
            col_names = input_pairs_next_line.strip('\n').split(' ')[1:]
            # ensure that readID is the very first column
            assert col_names[0] == 'readID'
            # keep the subsequent column names in a dictionary
            input_pairs_column_names = { col_name: col_id for col_id, col_name in enumerate(col_names[1:]) }
            # ensure that all the columns we need are provided
            for col in ['chrom1', 'pos51', 'pos31', 'strand1', 'chrom2', 'pos52', 'pos32', 'strand2']:
                assert col in input_pairs_column_names
        elif input_pairs_next_line[0] != '#':
            line = input_pairs_next_line.strip('\n').split('\t')
            if not read_id:
                # first alignment for this query sequence
                read_id = line[0]
            elif read_id != line[0]:
                # beyond last alignment for this query sequence
                # keep the alignment for later, and return the buffer collected so far
                break
            buffer.append(line[1:])
        # fetch the next alignment
        input_pairs_next_line = input_pairs.readline()

    return read_id, buffer


input_bam_next_line = None

def read_buffer_bam(input_bam_iterator, read_id):
    """
    Read all the alignments for the query sequence named 'read_id' from
    the .bam file.

    Returns: a list of BAM alignments.
    """

    global input_bam_next_line
    buffer = []

    while True:
        # fetch the next alignment
        if not input_bam_next_line:
            try:
                input_bam_next_line = next(input_bam_iterator)
            except StopIteration:
                # end of input file, return the buffer collected so far
                break

        # skip alignments until we reach the ones for the given query sequence
        if read_id == input_bam_next_line.query_name:
            # add alignment for the given query sequence to the buffer
            buffer.append(input_bam_next_line)
            input_bam_next_line = None
        elif buffer:
            # beyond last alignment for the given query sequence
            # keep the alignment for later, and return the buffer collected so far
            break

    return buffer


def AlignedSegment_qstart(segment):
    """
    Start index of the aligned query portion of the sequence (0-based, inclusive).

    Similar to query_alignment_start, but including the leading soft-clipped bases.
    """

    qstart = 0
    if segment.cigartuples:
        for op, count in segment.cigartuples:
            if op == pysam.CMATCH:
                break
            qstart += count
    return qstart


def AlignedSegment_qend(segment):
    """
    End index of the aligned query portion of the sequence (0-based, exclusive).

    Similar to query_alignment_end, but including the tailing soft-clipped bases.
    """

    return AlignedSegment_qstart(segment) + segment.query_alignment_length


def parse_buffer_pairtools(buffer_pairtools):
    """
    Parse the list of split lines from the .pairs.gz file (as returned by the
    function read_buffer_pairtools()) into two lists: pairtools_alignments and
    alignment_pairs.

    The list pairtools_alignments will contain sublists [chrom, start, end,
    abam_list, haplotype, phased] which define a
    single pairtools alignment. The last three items (initially [], None,
    None) are intended for: a list of BAM alignment lines (AlignedSegment),
    haplotype, and a flag indicating whether the alignment is fully phased.

    The list alignment_pairs will be a list of tuples (index1 in
    pairtools_alignments, index2 in pairtools_alignments, strand1, strand2), where
    the strand information is retained from pairtools alignment pairs and can
    differ from the strand in BAM alignments.

    Returns: tuple (pairtools_alignments, alignment_pairs).
    """

    pairtools_alignments = []
    alignment_pairs = []

    for pairtools_line in buffer_pairtools:
        def val(col_name):
            return pairtools_line[input_pairs_column_names[col_name]]
        assert val('pos1') == val('pos51')
        assert val('pos2') == val('pos52')

        strand1 = val('strand1')
        start1 = int(val('pos51')) if strand1 == '+' else int(val('pos31'))
        end1 = int(val('pos31')) if strand1 == '+' else int(val('pos51'))
        assert start1 <= end1
        strand2 = val('strand2')
        start2 = int(val('pos52')) if strand2 == '+' else int(val('pos32'))
        end2 = int(val('pos32')) if strand2 == '+' else int(val('pos52'))
        assert start2 <= end2

        pa = [val('chrom1'), start1, end1, [], None, None]
        if not pairtools_alignments:
            # first alignment in the first pairtools pair
            pairtools_alignments.append(pa)
        else:
            # usually the first alignment in a pairtools pair will be the same
            # as the second alignment in the previous pair; if is not the case
            # (e.g. unmapped longer gap in a read), add the first one as well
            if pa != pairtools_alignments[-1]:
                pairtools_alignments.append(pa)

        pa = [val('chrom2'), start2, end2, [], None, None]
        pairtools_alignments.append(pa)

        if val('chrom1') != '!' and val('chrom2') != '!':
            alignment_pairs.append((len(pairtools_alignments) - 2, len(pairtools_alignments) - 1,
                val('strand1'), val('strand2')))

    return pairtools_alignments, alignment_pairs


def split_buffer_bam(buffer_bam, min_mapq):
    """
    Split the contents of buffer_bam into two lists, one for BAM alignments
    originating from read 1, the other one for BAM alignments originating
    from read 2. Each list will contain sublists [qstart, qend, abam, index
    in pairtools_alignments], where qstart and qend are positions in read1 or
    read2, abam is the BAM alignment line (AlignedSegment) and the index in
    pairtools_alignments is initially set to None.

    Returns: a tuple (buffer_bam_read1, buffer_bam_read2) composed of the two
    lists.
    """

    buffer_bam_read1 = []
    buffer_bam_read2 = []
    read1_length = None
    read2_length = None
    for abam in buffer_bam:
        if not abam.is_supplementary:
            if abam.is_read1:
                assert read1_length is None
                read1_length = abam.query_length
            else:
                assert read2_length is None
                read2_length = abam.query_length
    assert read1_length
    assert read2_length

    for abam in buffer_bam:
        if abam.mapq < min_mapq:
            continue
        read_length = read1_length if abam.is_read1 else read2_length
        qstart = read_length - AlignedSegment_qend(abam) if abam.is_reverse else AlignedSegment_qstart(abam)
        qend = read_length - AlignedSegment_qstart(abam) if abam.is_reverse else AlignedSegment_qend(abam)
        if abam.is_read1:
            buffer_bam_read1.append([qstart, qend, abam, None])
        else:
            buffer_bam_read2.append([qstart, qend, abam, None])

    buffer_bam_read1.sort()
    buffer_bam_read2.sort()

    return buffer_bam_read1, buffer_bam_read2


def match_buffers_from_pairtools_and_bam(buffer_pairtools, buffer_bam, min_mapq):
    """
    For the query sequence, match a list of split lines from the .pairs.gz
    file with a list of BAM alignments as detailed below. Exclude the BAM
    alignments with mapping quality (MAPQ) lower than min_mapq.

    Returns: tuple (pairtools_alignments, alignment_pairs) as split by the
    function parse_buffer_pairtools(), with BAM alignment lines added to the
    matching pairtools alignments in pairtools_alignments.
    """

    pairtools_alignments, alignment_pairs = parse_buffer_pairtools(buffer_pairtools)
    buffer_bam_read1, buffer_bam_read2 = split_buffer_bam(buffer_bam, min_mapq)

    def debug_print():
        """
        Print the contents of the data structures for debugging purposes. Not
        used in production code.
        """

        print(buffer_bam[0].query_name)
        print("pairtools alignments:")
        for i, pa in enumerate(pairtools_alignments):
            print(i, pa)
        print("pairtools alignment pairs:")
        print(alignment_pairs)
        print("BAM alignments (followed by the index of matched pairtools alignment):")
        for qstart, qend, abam, pa_index in buffer_bam_read1 + list(reversed(buffer_bam_read2)):
            print('R1' if abam.is_read1 else 'R2', qstart + 1, qend, abam.reference_name,
                abam.reference_start + 1, abam.reference_end, '-' if abam.is_reverse else '+',
                abam.mapping_quality, abam.cigarstring, pa_index)
        print()

    def skip_unmapped_pairtools_alignments(pa_index, increment = 1):
        """
        Shift the pairtools alignment index 'pa_index' in the direction specified
        by 'increment' (-1 or +1) to skip any unmapped pairtools alignments.

        Returns: the shifted index. If the pairtools_alignments is exhausted
        during the procedure, returns -1 or len(pairtools_alignments).
        """

        if 0 <= pa_index and pa_index < len(pairtools_alignments):
            pa = pairtools_alignments[pa_index]
            chrom = pa[0]
            while chrom == '!':
                pa_index += increment
                if pa_index < 0 or len(pairtools_alignments) <= pa_index:
                    break
                pa = pairtools_alignments[pa_index]
                chrom = pa[0]
        return pa_index

    def match_pairtools_alignments_with_bam(buffer_bam_singleread, pa_index, increment = 1):
        """
        Match the pairtools alignments with BAM alignments
        in 'buffer_bam_singleread', starting from the pairtools alignment
        index 'pa_index' and continuing in the direction specified
        by 'increment' (-1 or +1). This is done by iterating over BAM
        alignments from 'buffer_bam_singleread, shifting 'pa_index' only to skip
        unmapped pairtools alignments or when a match is achieved.

        Returns: the shifted index. If the pairtools_alignments is exhausted
        during the procedure, returns -1 or len(pairtools_alignments).
        """

        # iterate over BAM alignments
        for i, abam_list in enumerate(buffer_bam_singleread):
            pa_index = skip_unmapped_pairtools_alignments(pa_index, increment)

            if pa_index < 0 or pa_index >= len(pairtools_alignments):
                # leave some BAM alignments unmatched after exhausting all the pairtools alignments
                break

            pa = pairtools_alignments[pa_index]
            chrom = pa[0]
            start = pa[1]
            end = pa[2]
            abam = abam_list[2]
            if chrom == abam.reference_name and \
                ((not abam.is_reverse and start == abam.reference_start + 1) or \
                (abam.is_reverse and abam.reference_end == end)):
                    # pairtools alignment matched with BAM alignment
                    buffer_bam_singleread[i][3] = pa_index
                    pa[3].append(abam)
                    pa_index += increment
            else:
                # skip the unmatched BAM alignment
                next

        pa_index = skip_unmapped_pairtools_alignments(pa_index, increment)
        return pa_index

    # Match the pairtools alignments separately with BAM alignments from
    # read1 (starting from the beginning of pairtools_alignments), and from
    # read2 (starting from the end of pairtools_alignments).
    pa_index_r1 = match_pairtools_alignments_with_bam(buffer_bam_read1, 0)
    pa_index_r2 = match_pairtools_alignments_with_bam(buffer_bam_read2, len(pairtools_alignments) - 1, -1)

    # After the procedure, all the mapped pairtools alignments in
    # pairtools_alignments should be matched to at least one BAM alignment:
    # on read1 (pairtools_alignments indices [0, pa_index_r1 - 1]) or on read2
    # (pairtools_alignments indices [pa_index_r2 + 1, len(pairtools_alignments) - 1]).
    # If this is not the case, print a warning.
    if pa_index_r1 < pa_index_r2 + 1:
        warnings.warn(f"failed to resolve pairtools alignments for read {buffer_bam[0].query_name}")
        # debug_print()

    return pairtools_alignments, alignment_pairs


def read_tsv(file):
    """
    Read SNV annotations in a custom tab-separated format.

    Input file should be a tab-separated text file containing SNV annotations,
    with header row and the following columns: chrom, 1-based position,
    haplotype 1, haplotype 2. Haplotype names are inferred from the header.

    Returns: a generator of tuples ((chrom, 1-based position), snv_dict),
    where snv_dict is a dictionary having SNV-distinguishing nucleotides as
    keys, and haplotypes as values.
    """

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
            warnings.warn("excluding SNVs that do not determine a haplotype")
        snv_dict = {nucl: hapl for nucl, hapl in zip(nucleotides, haplotypes) if nucleotides.count(nucl) == 1}

        # line[0] is chrom, line[1] is 1-based pos
        yield (line[0], int(line[1]) - 1), snv_dict


def annotate_haplotype(pairtools_alignments, snvs):
    """
    Annotate haplotypes for each pairtools alignment in pairtools_alignments,
    according to the provided SNVs.
    """

    for i, pa in enumerate(pairtools_alignments):
        haplotypes = set()

        # find the underlying SNVs
        for abam in pa[3]:
            comment = []

            for query_pos, reference_pos in abam.get_aligned_pairs(matches_only=True):
                snv = snvs.get((abam.reference_name, reference_pos))
                if snv:
                    haplotype = snv.get(abam.query_sequence[query_pos])
                    if haplotype:
                        haplotypes.add(haplotype)
                        # annotation: haplotype, reference position, query position, line 1/2
                        comment.append(haplotype + '|' + str(reference_pos) + '|' + str(query_pos) + '|'
                                       + '|'.join(snv.keys()))

            if comment:
                # save the information on the underlying SNV in the XS field
                abam.set_tag('XS', ','.join(comment), replace=False)

        # determine the haplotype annotation for all BAM alignments in the pairtools alignment
        if len(haplotypes) == 0:
            haplotype = 'unknown'
            phased = False
        elif len(haplotypes) == 1:
            haplotype = list(haplotypes)[0]
            phased = True
        else:
            haplotype = 'ambiguous'
            phased = False

        # for all BAM alignments in the pairtools alignment:
        for abam in pa[3]:
            # save the pairtools alignment index in the XG field
            abam.set_tag('XG', i)
            # save the haplotype annotation in the XH field
            abam.set_tag('XH', haplotype)

        pa[4] = haplotype
        pa[5] = phased


def write_buffer(pairtools_alignments, alignment_pairs, buffer_bam, output_tsv, output_bam = None, fully_phased = False, regions = []):
    """
    Write the alignment information for the query sequence to the
    output .tsv.gz file. If an output .bam file is provided, also write the
    annotated BAM alignments there.

    The output .tsv.gz file will have one row for each suitable line from
    the .pairs.gz file, with the following columns:
        chrom1: the chromosome of the alignment on side 1
        pos1: the 1-based genomic position of the outer-most (5’) mapped bp on side 1
        strand1: the strand of the alignment on side 1
        haplotype1: the haplotype of the alignment on side 1
        chrom2: the chromosome of the alignment on side 2
        pos2: the 1-based genomic position of the outer-most (5’) mapped bp on side 2
        strand2: the strand of the alignment on side 2
        haplotype2: the haplotype of the alignment on side 1
    Note: strand1 and strand2 are taken from alignment_pairs, other values from pairtools_alignments.

    The output .bam file will contain all the alignments of the query
    sequences that contributed to the .pairs.gz file. The following extra
    tags will be added:
        XS [string]: SNVs used to separate haplotypes
        XG [integer]: the index of pairtools alignment (0, 1, ...)
        XH [string]: haplotype (annotated for each pairtools alignment)
    """

    wrote_tsv = False

    for i, j, strand1, strand2 in alignment_pairs:
        alignments_fully_phased = pairtools_alignments[i][5] and pairtools_alignments[j][5]

        # write an alignment pair if we do not require fully phased alignments,
        # or if this pair of alignments is indeed fully phased
        if not fully_phased or alignments_fully_phased:

            chrom1 = pairtools_alignments[i][0]
            pos1 = pairtools_alignments[i][1 if strand1 == '+' else 2]
            haplotype1 = pairtools_alignments[i][4]

            chrom2 = pairtools_alignments[j][0]
            pos2 = pairtools_alignments[j][1 if strand2 == '+' else 2]
            haplotype2 = pairtools_alignments[j][4]

            # write an alignment pair if we do not filter for genomic region/haplotype,
            # or if this pair of alignments matches the filters
            if (not regions or
                (region_matched(regions[0], chrom1, pos1, haplotype1) and
                    (len(regions) < 2 or region_matched(regions[1], chrom2, pos2, haplotype2))) or
                (region_matched(regions[0], chrom2, pos2, haplotype2) and
                    (len(regions) < 2 or region_matched(regions[1], chrom1, pos1, haplotype1)))):

                output_tsv.write(chrom1 + '\t' + str(pos1) + '\t' + strand1 + '\t' + haplotype1 + '\t' +
                    chrom2 + '\t' + str(pos2) + '\t' + strand2 + '\t' + haplotype2 + '\n')

                wrote_tsv = True

    if output_bam and wrote_tsv:
        for abam in buffer_bam:
            output_bam.write(abam)


def main(args=None):
    """
    Main TransContactHiC function.
    """

    args = parse_arguments(args)

    print("Parsing SNVs...")
    snvs = dict(read_tsv(args.variants))

    print("Parsing .pairs.gz and BAM...")
    with gzip.open(args.input_pairs, mode='rt') as input_pairs:
        with pysam.AlignmentFile(args.input_bam, 'rb') as input_bam:
            input_bam_iterator = input_bam.fetch(until_eof=True)
            with gzip.open(args.output_tsv, 'wt') as output_tsv:
                output_tsv.write('chrom1\tpos1\tstrand1\thaplotype1\tchrom2\tpos2\tstrand2\thaplotype2\n')

                # construct BAM header for the output file
                bam_header = input_bam.header
                if not 'PG' in bam_header:
                    bam_header['PG'] = []
                bam_header['PG'].append({
                    'ID': 'TransContactHiC',
                    'PN': 'TransContactHiC',
                    'VN': str(__version__),
                    'CL': ' '.join(sys.argv)})

                output_bam = pysam.AlignmentFile(args.output_bam, 'wb', header=bam_header) if args.output_bam else None

                while True:
                    read_id, buffer_pairtools = read_buffer_pairtools(input_pairs)
                    if not read_id:
                        break
                    buffer_bam = read_buffer_bam(input_bam_iterator, read_id)
                    pairtools_alignments, alignment_pairs = match_buffers_from_pairtools_and_bam(buffer_pairtools, buffer_bam, args.min_mapq)
                    annotate_haplotype(pairtools_alignments, snvs)
                    write_buffer(pairtools_alignments, alignment_pairs, buffer_bam, output_tsv, output_bam,
                        fully_phased = args.fully_phased, regions = args.region)

                if args.output_bam:
                    output_bam.close()

if __name__ == '__main__':
    main(sys.argv[1:])
