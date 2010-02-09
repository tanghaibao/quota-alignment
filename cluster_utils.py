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


if __name__ == '__main__':

    usage = "Convert between .cluster (or .dag to .cluster) \n" \
            "%prog [options] input output "
    parser = OptionParser(usage)

    parser.add_option("-d", "--dag", dest="dag_fmt",
            action="store_true", default=False,
            help="dag formatted input [default: %default]")
    parser.add_option("-p", "--precision", dest="precision",
            action="store", type="int", default=1,
            help="convert scores into int(score*precision) " \
                "since MIP algorithm only deals with integer weights "\
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
    write_clusters(fw, clusters)

