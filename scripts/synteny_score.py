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
import sqlite3
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
                group = [group[i] for (y, i) in lis]
            else:
                score = lds_len
                group = [group[i] for (y, i) in lds]
                orientation = "-"
        else:
            xpos, ypos = zip(*group)
            score = min(len(set(xpos)), len(set(ypos)))

        group.sort()
        pos = bisect_left(group, (query, 0))
        left_flanker = group[0] if pos==0 else group[pos-1]
        right_flanker = group[-1] if pos==len(group) else group[pos] 

        # pick the closest flanker
        if abs(query - left_flanker[0]) < abs(query - right_flanker[0]):
            flanker = left_flanker
        else:
            flanker = right_flanker

        qflanker, syntelog = flanker
        if qflanker==query: 
            gray = "S"
            score += 1 # extra bonus for finding syntelog
        else:
            gray = "G"

        if score < cutoff: continue

        # y-boundary of the block
        left, right = group[0][1], group[-1][1]
        # this characterizes a syntenic region (left, right). syntelog is -1 if it's a gray gene
        syn_region = (syntelog, left, right, gray, orientation, score)
        regions.append(syn_region)

    return sorted(regions, key=lambda x: -x[-1]) # decreasing synteny score


def batch_query(qbed, sbed, all_data, options, c=None, transpose=False):

    cutoff = int(options.cutoff * options.window)
    window = options.window / 2
    sqlite = options.sqlite
    colinear = (options.scoring=="collinear")

    # process all genes present in the bed file 
    if transpose: 
        all_data = transposed(all_data)
        qbed, sbed = sbed, qbed

    all_data.sort()
    simple_bed = lambda x: (sbed[x].seqid, sbed[x].start)
    
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
                print "%s\t%s" % (qbed[r].accn, "\t".join(["na"]*5))
            for syntelog, left, right, gray, orientation, score in regions:
                query = qbed[r].accn

                left_chr, left_pos = simple_bed(left)
                right_chr, right_pos = simple_bed(right)

                anchor = sbed[syntelog].accn
                anchor_chr, anchor_pos = simple_bed(syntelog)
                # below is useful for generating the syntenic region in the coge url
                left_dist = abs(anchor_pos - left_pos) if anchor_chr==left_chr else 0
                right_dist = abs(anchor_pos - right_pos) if anchor_chr==right_chr else 0
                flank_dist = (max(left_dist, right_dist) / 10000 + 1) * 10000

                data = (query, anchor, gray, score, flank_dist, orientation)
                print "\t".join(map(str, data))
                if sqlite:
                    c.execute("insert into synteny values (?,?,?,?,?,?)", data)


def main(blast_file, options):
    qbed_file, sbed_file = options.qbed, options.sbed
    sqlite = options.sqlite
    
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

    c = None
    if options.sqlite:
        conn = sqlite3.connect(options.sqlite)
        c = conn.cursor()
        c.execute("drop table if exists synteny")
        c.execute("create table synteny (query integer, anchor text, gray varchar(1), score integer, dr integer, orientation varchar(1))")

    batch_query(qbed, sbed, all_data, options, c=c, transpose=False)
    batch_query(qbed, sbed, all_data, options, c=c, transpose=True)

    if sqlite:
        c.execute("create index q on synteny (query)")
        conn.commit()
        c.close()


if __name__ == '__main__':
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed",
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed",
            help="path to sbed")
    parser.add_option("--sqlite", dest="sqlite", default=None,
            help="write sqlite database")

    params_group = optparse.OptionGroup(parser, "Synteny parameters")
    params_group.add_option("--window", dest="window", type="int", default=40,
            help="synteny window size [default: %default]")
    params_group.add_option("--cutoff", dest="cutoff", type="float", default=.1, 
            help="the minimum number of anchors to call synteny [default: %default]")
    supported_scoring = ("collinear", "density")
    params_group.add_option("--scoring", dest="scoring", choices=supported_scoring, default="collinear",
            help="scoring scheme, must be one of " + str(supported_scoring) +" [default: %default]")

    parser.add_option_group(params_group)

    (options, blast_files) = parser.parse_args()

    if not (len(blast_files) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    main(blast_files[0], options)

