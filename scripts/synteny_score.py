#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog rice_sorghum.blastp.filtered --qbed=rice.nolocaldups.bed --sbed=sorghum.nolocaldups.bed

Given a blast, we find the syntenic regions for every single gene. The algorithm
works by expanding the query gene to a window centered on the gene. A single
linkage algorithm follows that outputs the synteny block.

The result looks like the following:
Os01g0698300    Sb03g032090     S    7     +
Os01g0698500    Sb03g032140     G    11    +

The pairs (A, B) -- A is query, and then B is the syntenic region found
G is "Gray gene", which means it does not have match to the region (fractionated
or inserted). In this case, a right flanker is used to represent the region.
S is "Syntelog", which means it has a match to the region. In this case, the
match itself is used to represent the region.

The number in the 4th column is the synteny score. For the same query, it is
ordered with decreasing synteny score. The last column means orientation. "+" is
the same direction.

Finally, accuracy can be improved if the raw BLAST is available for validation.
The dup hits were removed in filtered BLAST, but should still be in the raw
BLAST. We can have a 2nd pass at the raw BLAST, flipping any false 'proxy' or
'syntelog'. In 2nd pass, we are not allowed to change the synteny score, only
validating the hits. The validation is enabled by appending the following
settings:

    --lift rice_sorghum.blast --qbed=rice.bed --sbed=sorghum.bed
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


def get_flanker(group, query):
    """
    >>> get_flanker([(370, 15184), (372, 15178), (373, 15176), (400, 15193)],  385)
    ((373, 15176), (400, 15193), True)

    >>> get_flanker([(124, 13639), (137, 13625)], 138)
    ((137, 13625), (137, 13625), False)
    """
    group.sort()
    #print >>sys.stderr, group, query
    pos = bisect_left(group, (query, 0))
    left_flanker = group[0] if pos==0 else group[pos-1]
    right_flanker = group[-1] if pos == len(group) else group[pos]
    # pick the closest flanker
    if abs(query - left_flanker[0]) < abs(query - right_flanker[0]):
        flanker, other = left_flanker, right_flanker
    else:
        flanker, other = right_flanker, left_flanker

    flanked = not (pos==0 or pos == len(group) or flanker == query)
    #print >>sys.stderr, flanker, flanked, "\n\n"

    return flanker, other, flanked


def find_synteny_region(query, sbed, data, window, cutoff, colinear=False):
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
        pos1, pos2 = ia[1], ib[1]
        if pos2 - pos1 < window and sbed[pos1].seqid==sbed[pos2].seqid:
            g.join(ia, ib)

    #print list(g)

    for group in sorted(g):
        (qflanker, syntelog), (far_flanker, far_syntelog), flanked = get_flanker(group, query)

        # run a mini-dagchainer here, take the direction that gives us most anchors
        orientation = "+"
        if colinear:
            y_indexed_group = [(y, i) for i, (x, y) in enumerate(group)]
            lis = longest_increasing_subsequence(y_indexed_group)
            lds = longest_decreasing_subsequence(y_indexed_group)

            if len(lis) >= len(lds):
                track = lis
            else:
                track = lds
                orientation = "-"

            group = [group[i] for (y, i) in track]

        xpos, ypos = zip(*group)
        score = min(len(set(xpos)), len(set(ypos)))

        if qflanker==query:
            gray = "S"
        else:
            gray = "G" if not flanked else "F"
            score -= 1 # slight penalty for not finding syntelog

        if score < cutoff: continue

        # y-boundary of the block
        left, right = group[0][1], group[-1][1]
        # this characterizes a syntenic region (left, right). syntelog is -1 if it's a gray gene
        syn_region = (syntelog, far_syntelog, left, right, gray, orientation, score)
        regions.append(syn_region)

    return sorted(regions, key=lambda x: -x[-1]) # decreasing synteny score


def batch_query(qbed, sbed, all_data, options, c=None, transpose=False):

    cutoff = int(options.cutoff * options.window)
    window = options.window / 2
    sqlite = options.sqlite
    colinear = (options.scoring=="collinear")
    qnote, snote = options.qnote, options.snote

    # process all genes present in the bed file
    if transpose:
        all_data = transposed(all_data)
        qbed, sbed = sbed, qbed
        qnote, snote = snote, qnote

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
            regions = find_synteny_region(r, sbed, data, window, cutoff, colinear=colinear)
            #if not regions: print "%s\t%s" % (qbed[r].accn, "\t".join(["na"]*5))
            for syntelog, far_syntelog, left, right, gray, orientation, score in regions:
                query = qbed[r].accn

                left_chr, left_pos = simple_bed(left)
                right_chr, right_pos = simple_bed(right)

                anchor = sbed[syntelog].accn
                anchor_chr, anchor_pos = simple_bed(syntelog)
                # below is useful for generating the syntenic region in the coge url
                left_dist = abs(anchor_pos - left_pos) if anchor_chr==left_chr else 0
                right_dist = abs(anchor_pos - right_pos) if anchor_chr==right_chr else 0
                flank_dist = (max(left_dist, right_dist) / 10000 + 1) * 10000

                far_syntelog = sbed[far_syntelog].accn

                left_pos, right_pos = sorted((left_pos, right_pos))
                data = [query, anchor, gray, score, flank_dist, orientation, far_syntelog]
                if sqlite:
                    c.execute("insert into synteny values (?,?,?,?,?,?,?,?)", data[:6]+[qnote,snote])
                else:
                    print "\t".join(map(str, data))


def main(blast_file, options):
    qbed_file, sbed_file = options.qbed, options.sbed
    sqlite = options.sqlite

    print >>sys.stderr, "read annotation files %s and %s" % (qbed_file, sbed_file)

    try:
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

    except IndexError:
        print >>sys.stderr, "No results were found in the query or source bed file"

    c = None
    if options.sqlite:
        conn = sqlite3.connect(options.sqlite)
        c = conn.cursor()
        c.execute("drop table if exists synteny")
        c.execute("create table synteny (query text, anchor text, gray varchar(1), score integer, dr integer, "
                "orientation varchar(1), qnote text, snote text)")

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

    coge_group = optparse.OptionGroup(parser, "CoGe-specific options")
    coge_group.add_option("--sqlite", dest="sqlite", default=None,
            help="write sqlite database")
    coge_group.add_option("--qnote", dest="qnote", default="null",
            help="query dataset group id")
    coge_group.add_option("--snote", dest="snote", default="null",
            help="subject dataset group id")

    params_group = optparse.OptionGroup(parser, "Synteny parameters")
    params_group.add_option("--window", dest="window", type="int", default=40,
            help="synteny window size [default: %default]")
    params_group.add_option("--cutoff", dest="cutoff", type="float", default=.1,
            help="the minimum number of anchors to call synteny [default: %default]")
    supported_scoring = ("collinear", "density")
    params_group.add_option("--scoring", dest="scoring", choices=supported_scoring, default="collinear",
            help="scoring scheme, must be one of " + str(supported_scoring) +" [default: %default]")

    params_group = optparse.OptionGroup(parser, "Accuracy improvement")
    params_group.add_option("--lift",
            help="Raw BLAST file for validation [default: disabled]")
    params_group.add_option("--qbedlift", help="Path to qbed for validation")
    params_group.add_option("--sbedlift", help="Path to sbed for validation")

    parser.add_option_group(coge_group)
    parser.add_option_group(params_group)

    (options, blast_files) = parser.parse_args()

    if not (len(blast_files) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    main(blast_files[0], options)

