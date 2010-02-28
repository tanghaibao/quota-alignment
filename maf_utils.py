#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This script uses the great bx-python tool <http://bitbucket.org/james_taylor/bx-python/>
to parse the .maf file in order to convert the blocks into .qa format.

"""

import os
import sys
from bx.align import maf
from bx import interval_index_file


# Stolen from Brad Chapman's blog
# http://bcbio.wordpress.com/2009/07/26/sorting-genomic-alignments-using-python/
def build_index(in_file, index_file):
    indexes = interval_index_file.Indexes()

    in_handle = open(in_file)
    reader = maf.Reader(in_handle)
    while 1:
        pos = reader.file.tell()
        rec = reader.next()
        if rec is None: break
        for c in rec.components:
            indexes.add(c.src, c.forward_strand_start,
                    c.forward_strand_end, pos, max=c.src_size )

    index_handle = open(index_file, "w")
    print >>sys.stderr, "build %s for fast record retrieval" % index_file
    indexes.write(index_handle)


def get_alignments(maf_file):
    
    base, ext = os.path.splitext(maf_file)

    # build index file for fast retrieval
    index_file = maf_file + ".index"
    if not os.path.exists(index_file):
        build_index(maf_file, index_file)
    index = maf.Indexed(maf_file, index_file)

    fp = file(maf_file)
    reader = maf.Reader(fp)

    alignments = []
    while 1:
        rec = reader.next()
        if rec is None: break
        intervals = []
        for c in rec.components:
            chr, left, right, strand, weight = c.src, c.forward_strand_start, \
                    c.forward_strand_end, c.strand, rec.score
            intervals.append((chr, left, right, strand, weight))
        alignments.append(intervals)

    fp.close()

    return alignments 



if __name__ == '__main__':
    
    from optparse import OptionParser

    usage = "%prog [options] infile"
    parser = OptionParser(usage=usage, version="%prog 1.0")

    (options, args) = parser.parse_args()
    try:
        maf_file = args[0]
    except:
        print >>sys.stderr, "please send input file name"
        sys.exit(parser.print_help())

    alignments = get_alignments(maf_file)
    for alignment in alignments:
        print alignment

