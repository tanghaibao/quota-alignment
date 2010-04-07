#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog qa_file --qbed query.bed --sbed subject.bed

convert qa_file back to the original gene names
"""

import os.path as op
import itertools
import sys

from bed_utils import Bed, RawLine

def qa_to_pairs(qa_file, qbed, sbed):

    for line in open(qa_file):
        if line[0] == "#":
            print line,
            continue
        s = RawLine(line)
        query = qbed[s.pos_a].accn
        subject = sbed[s.pos_b].accn
        print "\t".join((query, subject))


if __name__ == "__main__":

    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed",
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed",
            help="path to sbed")

    (options, args) = parser.parse_args()

    if not (len(args) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    qbed = Bed(options.qbed)
    sbed = Bed(options.sbed)

    qa_file = args[0]

    qa_to_pairs(qa_file, qbed, sbed)

