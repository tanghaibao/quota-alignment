#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog rice_sorghum.blastp.filtered --qbed=rice.nolocaldups.bed --sbed=sorghum.nolocaldups.bed

Given a blast, we find the syntenic regions for every single gene. The algorithm works by expanding the query gene to a window centered on the gene. A single linkage algorithm follows that outputs the synteny block. 

The result looks like the following:
Os01g0698300    Sb03g032090     S    7     +
Os01g0698500    Sb03g032140     G    11    +

The pairs (A, B) -- A is query, and then B is the syntenic region found
G is "Gray gene", which means it does not have match to the region (fractionated or inserted). In this case, a right flanker is used to represent the region.
S is "Syntelog", which means it has a match to the region. In this case, the match itself is used to represent the region.
The number in the 4th column is the synteny score. For the same query, it is ordered with decreasing synteny score. The last column means orientation. "+" is same direction.
"""

import sys
import itertools
import os.path as op
from bisect import bisect_left

sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from grouper import Grouper
from bed_utils import Bed, BlastLine
from lis import longest_increasing_subsequence, longest_decreasing_subsequence


def transposed(data):
    x, y = zip(*data)
    return zip(y, x)


def find_synteny_region(query, data, window, cutoff, colinear=False):
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
        
        group.sort()
        orientation = "+"

        # run a mini-dagchainer here, take the direction that gives us most anchors 
        if colinear:
            y_indexed_group = [(y, i) for i, (x, y) in enumerate(group)]
            lis = longest_increasing_subsequence(y_indexed_group)
            lds = longest_decreasing_subsequence(y_indexed_group)
            
            lis_len, lds_len = len(lis), len(lds)
            if lis_len >= lds_len: 
                score = lis_len
                group = [group[i] for y, i in lis]
            else:
                score = lds_len
                group = [group[i] for y, i in lds]
                orientation = "-"
        else:
            xpos, ypos = zip(*group)
            # get the number of unique positions
            score = min(len(set(xpos)), len(set(ypos)))

        pos = bisect_left(group, (query, 0))
        flanker = group[-1] if pos==len(group) else group[pos]
        gray_gene = "G"
        if flanker[0]==query: 
            gray_gene = "S"
            score += 1 # extra bonus for finding syntelog

        if score >= cutoff: 
            syn_region = [flanker, gray_gene, score, orientation]
            regions.append(syn_region)

    return sorted(regions, key=lambda x: -x[2]) # decreasing synteny score


def batch_query(qbed, sbed, all_data, window, cutoff, colinear=False, transpose=False):
    # process all genes present in the bed file 
    if transpose: 
        all_data = transposed(all_data)
        qbed, sbed = sbed, qbed

    all_data.sort()
    for seqid, ranks in itertools.groupby(qbed.get_simple_bed(), key=lambda x: x[0]):
        ranks = [x[1] for x in ranks]
        for r in ranks:
            rmin = max(r-window, ranks[0])
            rmax = min(r+window+1, ranks[-1])
            rmin_pos = bisect_left(all_data, (rmin, 0))
            rmax_pos = bisect_left(all_data, (rmax, 0))
            data = all_data[rmin_pos:rmax_pos]
            regions = find_synteny_region(r, data, window, cutoff, colinear=colinear)
            if not regions:
                print "%s\tna\tna\tna" % (qbed[r].accn)
            for pivot, label, score, orientation in regions:
                print "%s\t%s\t%s\t%d\t%s" % (qbed[r].accn, sbed[pivot[1]].accn, \
                        label, score, orientation)


def main(blast_file, options):
    qbed_file, sbed_file = options.qbed, options.sbed
    
    window = options.window
    cutoff = options.cutoff
    colinear = options.colinear

    print >>sys.stderr, "read annotation files %s and %s" % (qbed_file, sbed_file)
    qbed = Bed(qbed_file)
    sbed = Bed(sbed_file)

    qorder = qbed.get_order()
    sorder = sbed.get_order()

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

    batch_query(qbed, sbed, all_data, window, cutoff, colinear=colinear, transpose=False)
    batch_query(qbed, sbed, all_data, window, cutoff, colinear=colinear, transpose=True)


if __name__ == '__main__':
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed",
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed",
            help="path to sbed")

    params_group = optparse.OptionGroup(parser, "Synteny parameters")
    params_group.add_option("--window", dest="window", type="int", default=20,
            help="synteny window size [default: %default]")
    params_group.add_option("--cutoff", dest="cutoff", type="int", default=4, 
            help="the minimum number of anchors to call synteny [default: %default]")
    params_group.add_option("--nocolinear", dest="colinear", action="store_false",
            default=True, help="don't expect collinearity? [default: collinear regions]")

    parser.add_option_group(params_group)

    (options, blast_files) = parser.parse_args()

    if not (len(blast_files) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    main(blast_files[0], options)

