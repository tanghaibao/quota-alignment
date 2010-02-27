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
    with open(in_file) as in_handle:
        reader = maf.Reader(in_handle)
        while 1:
            pos = reader.file.tell()
            rec = reader.next()
            if rec is None:
                break
            for c in rec.components:
                indexes.add(c.src, c.forward_strand_start,
                        c.forward_strand_end, pos, max=c.src_size )

    with open(index_file, "w") as index_handle:
        print >>sys.stderr, "build %s for fast retrieval" % index_file
        indexes.write(index_handle)


def main(options, args):
    
    in_file = args[0]
    base, ext = os.path.splitext(in_file)

    # build index file
    index_file = in_file + ".index"
    if not os.path.exists(index_file):
        build_index(in_file, index_file)
    index = maf.Indexed(in_file, index_file)

    fp = file(in_file)
    reader = maf.Reader(fp)

    j = 0
    rec_info = []
    while 1:
        pos = reader.file.tell()
        rec_info.append((j/2, pos))   # position of alignment j in file
        rec = reader.next()
        if rec is None:
            break
        for c in rec.components:
            chromosome, left, right, weight = c.src, c.forward_strand_start, \
                    c.forward_strand_end, rec.score
            j += 1

    fp.close()



if __name__ == '__main__':
    
    from optparse import OptionParser

    usage = "usage: %prog [options] infile"
    parser = OptionParser(usage=usage, version="%prog 1.0")

    (options, args) = parser.parse_args()
    if len(args)!=1:
        parser.error("incorrect number of arguments")
        parser.print_help()
        sys.exit()

    main(options, args)

