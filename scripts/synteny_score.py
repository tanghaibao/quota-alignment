#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog rice_sorghum.blastp.filtered --qbed=rice.nolocaldups.bed --sbed=sorghum.nolocaldups.bed

Given a blast, we find the syntenic regions for every single gene. The algorithm works by expanding the query gene to a window centered on the gene. A single linkage algorithm follows that outputs the synteny block. 

The result looks like the following:
Os01g0698300    Sb03g032090     S
Os01g0698500    Sb03g032140     G

The pairs (A, B) -- A is query, and then B is the syntenic region found
G is "Gray gene", which means it does not have match to the region (fractionated or inserted). In this case, a right flanker is used to represent the region.
S is "Syntelog", which means it has a match to the region. In this case, the match itself is used to represent the region.

"""

import sys
import numpy as np
import itertools
import os.path as op
from bisect import bisect_left

sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from grouper import Grouper
from bed_utils import Bed, BlastLine, get_order


simple_bed = lambda bed: [(b.seqid, i) for (i, b) in enumerate(bed.beds)]

def transposed(data):
    x, y = zip(*data)
    return zip(y, x)


def find_synteny_region(query, data, window, cutoff):
    # get all synteny blocks for a query, algorithm is single linkage
    # anchors are a window centered on query 
    # two categories of syntenic regions depending on what query is:
    # (Syntelog): syntenic region is denoted by the syntelog
    # (Gray gene): syntenic region is marked by the closest flanker 
    
    regions = []
    ysorted = sorted(data, key=lambda x:x[1])
    g = Grouper()

    a, b = itertools.tee(ysorted)
    next(b, None)
    for ia, ib in itertools.izip(a, b):
        if ib[1]-ia[1] < window: g.join(ia, ib)

    for group in sorted(g):
        if len(group) < cutoff: continue
        group.sort()
        pos = bisect_left(group, (query, 0))
        if pos == len(group): 
            syn_region = (group[-1], "G")
        else:
            flanker = group[pos]
            syn_region = (flanker, "S") if flanker[0]==query else (flanker, "G")
        regions.append(syn_region)

    return regions 


def batch_query(qbed, sbed, all_data, window, cutoff, transpose=False):
    # process all genes present in the bed file 
    if transpose: 
        all_data = transposed(all_data)
        qbed, sbed = sbed, qbed

    all_data.sort()
    for seqid, ranks in itertools.groupby(simple_bed(qbed), key=lambda x: x[0]):
        ranks = (x[1] for x in ranks)
        for r in ranks:
            rmin_pos = bisect_left(all_data, (r-window, 0))
            rmax_pos = bisect_left(all_data, (r+window+1, 0))
            data = all_data[rmin_pos:rmax_pos]
            regions = find_synteny_region(r, data, window, cutoff)
            for pivot, label in regions:
                print "%s\t%s\t%s" % (qbed[r].accn, sbed[pivot[1]].accn, label)


def main(blast_file, options):
    qbed_file, sbed_file = options.qbed, options.sbed
    
    window = options.window
    cutoff = options.cutoff

    print >>sys.stderr, "read annotation files %s and %s" % (qbed_file, sbed_file)
    qbed = Bed(qbed_file)
    sbed = Bed(sbed_file)

    qorder = get_order(qbed)
    sorder = get_order(sbed)

    fp = file(blast_file)
    print >>sys.stderr, "read BLAST file %s (total %d lines)" % \
            (blast_file, sum(1 for line in fp))
    fp.seek(0)
    blasts = sorted([BlastLine(line) for line in fp], \
            key=lambda b: b.score, reverse=True)

    all_data = []
    for b in blasts:
        query, subject = b.query, b.subject
        if query not in qorder or subject not in sorder: continue
        qi, q = qorder[query]
        si, s = sorder[subject]
        all_data.append((qi, si))

    batch_query(qbed, sbed, all_data, window, cutoff, transpose=False)
    batch_query(qbed, sbed, all_data, window, cutoff, transpose=True)


if __name__ == '__main__':
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed",
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed",
            help="path to sbed")

    params_group = optparse.OptionGroup(parser, "BLAST filters")
    params_group.add_option("--window", dest="window", type="int", default=30,
            help="synteny window size [default: %default]")
    params_group.add_option("--cutoff", dest="cutoff", type="int", default=4, 
            help="the minimum number of anchors to call synteny [default: %default]")

    parser.add_option_group(params_group)

    (options, blast_files) = parser.parse_args()

    if not (len(blast_files) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    main(blast_files[0], options)

