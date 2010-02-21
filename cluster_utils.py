#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
I/O related functions to perform on the cluster (or dag) file

a dag file looks like:
## alignment a3068_scaffold_1 vs. b8_1 Alignment #1  score = 102635.0 (num aligned pairs: 2053): 
a3068_scaffold_1        scaffold_1||548785||550552||scaffold_100153.1||-1||CDS||30767256||87.37 140     140     b8_1    1||427548||427811||AT1G02210||-1||CDS||20183105||87.37     172     172     1.000000e-250   50

a plain cluster file looks like:
# cluster score 1263.482
2       2840    3       3880    1.647
2       2859    3       3867    2.560

"""

import math
import sys
from itertools import groupby
from optparse import OptionParser


# copied from brentp's dag_chainer.py
CONSTANT_MATCH_SCORE=None
MAX_MATCH_SCORE=50.0

def scoringF(evalue, constant_match=CONSTANT_MATCH_SCORE, max_match=MAX_MATCH_SCORE):
    if not constant_match is None:
        return constant_match
    if evalue == 0.0: return max_match
    matchScore = 10 * -math.log10(evalue);
    matchScore = int(matchScore +.5) / 10
    return max_match if matchScore > max_match else matchScore


def read_clusters(filename, precision=1, dag_fmt=False):

    fp = file(filename)

    # clusters contain only bounds, point_clusters contain all points
    clusters = [] 
    
    total_lines = sum(1 for row in fp)
    fmt = "dag" if dag_fmt else "cluster"
    print >>sys.stderr, "total lines in %s file (%d)" % (fmt, total_lines)
    fp.seek(0)
    row = fp.readline()
    j = 1
    while row:
        if row.strip() == "": break
        row = fp.readline()
        cluster = []
        while row and row[0] != "#":
            atoms = row.rstrip().split("\t")
            if row.strip()== "": break
            if dag_fmt:
                ca, a, cb, b, evalue = atoms[0], atoms[2], atoms[4], \
                                       atoms[6], atoms[8]
                score = int(scoringF(float(evalue)))
            else: # handle my own cluster fmt
                ca, a, cb, b, score = atoms
            score = int(float(score) * precision) 
            a, b = int(a), int(b)
            gene1, gene2 = (ca, a), (cb, b)
            cluster.append((gene1, gene2, score))
            row = fp.readline()

        if len(cluster) == 0: continue
        clusters.append(cluster)

    return clusters


def write_clusters(filehandle, clusters):
    for cluster in clusters:
        cluster_score = sum(x[-1] for x in cluster)
        filehandle.write("# cluster score %d \n" % (cluster_score)) 
        for gene1, gene2, score in cluster:
            # gene is (name, posn)
            filehandle.write("%s\t%d\t" % gene1 )
            filehandle.write("%s\t%d\t" % gene2 )
            filehandle.write("%d\n" % score )
    print >>sys.stderr, "wrote clusters to %s" % filehandle.name


def make_range(clusters):
    # convert to interval ends from a list of anchors
    eclusters = [] 
    for cluster in clusters:
        xlist, ylist, scores = zip(*cluster)
        score = sum(scores)

        xchr, xmin = min(xlist) 
        xchr, xmax = max(xlist)
        ychr, ymin = min(ylist) 
        ychr, ymax = max(ylist)

        eclusters.append(((xchr, xmin, xmax), (ychr, ymin, ymax), score))

    return eclusters


def make_projection(clusters):
    # let the x-projection of the blocks 1..n
    # output the y-projection sequence for downstream permutation analysis
    # both lists are 2-level as we have integer sequences for multiple chromosomes
    
    clusters.sort()
    x_projection, y_projection = [], []

    for i, cluster in enumerate(clusters):
        block_id = i+1
        cluster.sort()
        xlist, ylist, scores = zip(*cluster)

        xchr, xfirst = xlist[0]
        x_projection.append((xchr, xfirst, block_id))

        ychr, yfirst = ylist[0]
        ychr, ylast = ylist[-1]
        if yfirst < ylast: 
            sign = 1
        else:
            yfirst, ylast = ylast, yfirst
            sign = -1
        y_projection.append((ychr, yfirst, sign * block_id))

    y_projection.sort()

    return x_projection, y_projection 


def print_intseq(projection, filehandle):
    # from a list of (chr, pos, signed_id) => a nested list of multichromosome
    # signed integers
    g = groupby(projection, lambda x: x[0])
    intseq = [list(x[2] for x in blocks) for chr, blocks in g]
    for s in intseq:
        print >>filehandle, " ".join(str(x) for x in s) + "$"


def print_grimm(clusters, filehandle=sys.stdout):
    # GRIMM-style output, for more info, see http://nbcr.sdsc.edu/GRIMM/grimm.cgi
    x_projection, y_projection = make_projection(clusters)

    print >>filehandle, ">genome x"
    print_intseq(x_projection, filehandle)
    print >>filehandle, ">genome y"
    print_intseq(y_projection, filehandle)


def interval_union(intervals):
    # return total size of intervals, expect interval as (chr, left, right)
    intervals.sort()

    total_len = 0
    cur_chr, cur_left, cur_right = intervals[0] # left-most interval
    for interval in intervals:
        # open a new interval if left > cur_right or chr != cur_chr
        if interval[1] > cur_right or interval[0] != cur_chr:
            total_len += cur_right - cur_left + 1
            cur_chr, cur_left, cur_right = interval
        else:
            # update cur_right
            cur_right = max(interval[2], cur_right)

    # the last one
    total_len += cur_right - cur_left + 1

    return total_len


def calc_coverage(clusters, self=False):
    # calculates the length that's covered, for coverage statistics
    eclusters = make_range(clusters)

    intervals_x = [x[0] for x in eclusters]
    intervals_y = [x[1] for x in eclusters]
    
    if self:
        total_len_x = interval_union(intervals_x+intervals_y)
        total_len_y = total_len_x
    else:
        total_len_x = interval_union(intervals_x)
        total_len_y = interval_union(intervals_y)

    return total_len_x, total_len_y



if __name__ == '__main__':

    usage = "Convert between .cluster (or .dag to .cluster) \n" \
            "%prog [options] input output \n" \
            ".. if output not given, will write to stdout"
    parser = OptionParser(usage)

    parser.add_option("-d", "--dag", dest="dag_fmt",
            action="store_true", default=False,
            help="dag formatted input [default: %default]")
    parser.add_option("-p", "--precision", dest="precision",
            action="store", type="int", default=1,
            help="convert scores into int(score*precision) " \
                "since MIP algorithm only deals with integer weights "\
                "[default: %default]")
    parser.add_option("-c", "--calc_coverage", dest="calc_coverage",
            action="store_true", default=False,
            help="calculate the total length these clusters occupy "\
                "[default: %default]")
    parser.add_option("-r", "--print_grimm", dest="print_grimm",
            action="store_true", default=False,
            help="print two integer sequences for permutation GRIMM analysis "\
                 "[default: %default]")

    (options, args) = parser.parse_args()

    try:
        dag_file = args[0]
        if len(args) == 2:
            cluster_file = args[1]
            fw = file(cluster_file, "w")
        else:
            fw = sys.stdout
    except:
        sys.exit(parser.print_help())

    # file format conversion
    clusters = read_clusters(dag_file, options.precision, options.dag_fmt)

    if options.print_grimm:
        print_grimm(clusters)
        sys.exit(0)

    write_clusters(fw, clusters)

    if options.calc_coverage:
        total_len_x, total_len_y = calc_coverage(clusters)
        print >>sys.stderr, "Total length on x-axis:", total_len_x 
        print >>sys.stderr, "Total length on y-axis:", total_len_y

