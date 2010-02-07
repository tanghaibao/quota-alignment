#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
this python program can do two things:
1. merge dags by recursively merging dag bounds
2. build conflict graph where edges represent 1d-`overlap' between blocks
3. feed the data into the linear programming solver.

a dag file looks like:
## alignment a3068_scaffold_1 vs. b8_1 Alignment #1  score = 102635.0 (num aligned pairs: 2053): 
a3068_scaffold_1        scaffold_1||548785||550552||scaffold_100153.1||-1||CDS||30767256||87.37 140     140     b8_1    1||427548||427811||AT1G02210||-1||CDS||20183105||87.37     172     172     1.000000e-250   50

a plain cluster file looks like:
# mergedcluster score 1263.482
2       2840    3       3880    1.647
2       2859    3       3867    2.560

"""

import sys
import pprint
import cStringIO
import operator
from mystruct import Grouper
from subprocess import Popen
from optparse import OptionParser
from dag_chainer import scoringF


def range_mergeable(a, b):
    # 1-d version of box_mergeable
    a_chr, a_min, a_max = a
    b_chr, b_min, b_max = b
    # must be on the same chromosome
    if a_chr!=b_chr: return False 
    
    # make sure it is end-to-end merge, and within distance cutoff
    return (a_min <= b_max + Nmax) and \
           (b_min <= a_max + Nmax)


def box_mergeable(boxa, boxb):

    boxa_xrange, boxa_yrange, _ = boxa
    boxb_xrange, boxb_yrange, _ = boxb

    return range_mergeable(boxa_xrange, boxb_xrange) and \
           range_mergeable(boxa_yrange, boxb_yrange)


def range_relaxed_overlap(a, b):
    # 1-d version of box_relaxed_overlap
    a_chr, a_min, a_max = a
    b_chr, b_min, b_max = b
    # must be on the same chromosome
    if a_chr!=b_chr: return False 
    
    # Kmax is the allowable overlap level, it is usually Nmax
    # however some very small segments might slip through
    # therefore we also control for the segment size
    Kmax = min(a_max-a_min, b_max-b_min, Nmax)
    # to handle ranges that slightly overlap
    return (a_min <= b_max - Kmax) and \
           (b_min <= a_max - Kmax)


def box_relaxed_overlap(boxa, boxb):

    boxa_xrange, boxa_yrange, _ = boxa
    boxb_xrange, boxb_yrange, _ = boxb

    return range_relaxed_overlap(boxa_xrange, boxb_xrange) or \
           range_relaxed_overlap(boxa_yrange, boxb_yrange)


def make_range(clusters):
    # convert to interval ends from a list of anchors
    eclusters = [] 
    for cluster in clusters:
        xlist = [a[0] for a in cluster]
        ylist = [a[1] for a in cluster]
        score = sum(a[-1] for a in cluster)
        
        xchr, xmin = min(xlist) 
        xchr, xmax = max(xlist)
        ychr, ymin = min(ylist) 
        ychr, ymax = max(ylist)

        eclusters.append(((xchr, xmin, xmax), (ychr, ymin, ymax), score))

    return eclusters


#___merge clusters to combine inverted blocks______________________________

def merge_clusters(chain, clusters):

    # there are breakpoints induced by inversions and by translocations
    # inversion-breakpoints are excessive breakpoints that I want to remove
    
    chain_num = len(chain)
    eclusters = make_range(clusters)
    #pprint.pprint(eclusters)

    mergeables = Grouper() # disjoint sets of clusters that can be merged
    # check all pairwise combinations
    for i in xrange(chain_num):
        ci = chain[i]
        mergeables.join(ci)
        for j in xrange(i+1, chain_num):
            cj = chain[j]
            if box_mergeable(eclusters[ci], eclusters[cj]):
                mergeables.join(ci, cj)

    to_merge = {} 
    for mergeable in mergeables:
        for m in mergeable:
            to_merge[m] = min(mergeables[m])

    merged_chain = []
    for c in chain:
        if to_merge[c]==c: # i.e. parent of mergeables
            merged_chain.append(c)

    # refresh clusters list, merge chains
    for k, v in to_merge.iteritems():
        if to_merge[k]!=k: # i.e. not map to self
            clusters[v].extend(clusters[k])

    # maintain the x-sort
    [cluster.sort() for cluster in clusters]

    # nothing is merged
    updated = (len(merged_chain) != chain_num)
    return merged_chain, updated


def recursive_merge_clusters(chain, clusters):

    # as some rearrangment patterns are recursive, the extension of blocks
    # will take several iterations
    while 1: 
        chain, updated = merge_clusters(chain, clusters)
        print >>sys.stderr, "merging..."
        if not updated: break

    return chain, clusters


#___formulate linear programming instance____________________________________

def construct_graph(clusters):

    # check pairwise cluster comparison, if they overlap then mark edge as `conflict'

    eclusters = make_range(clusters)
    # (1-based index, synteny_score)
    nodes = [(i+1, c[-1]) for i, c in enumerate(eclusters)]
    edges = []
    nnodes = len(nodes)
    for i in xrange(nnodes):
        for j in xrange(i+1, nnodes):
            if box_relaxed_overlap(eclusters[i], eclusters[j]): 
                edges.append((i+1, j+1))

    return nodes, edges


def format_lp(nodes, edges):

    """
    \* Problem: lp_test *\

    Maximize
     obj: x(1) + 2 x(2) + 3 x(3) + 4 x(4)

    Subject To
     r_1: x(1) + x(2) <= 1

    End
    
    """
    lp_handle = cStringIO.StringIO()
    lp_handle.write("\* Problem: synteny *\ \n\n")
    
    lp_handle.write("Maximize\n obj: ")
    for i, score in nodes:
        lp_handle.write("+ %d x(%d) " % (score, i))
    lp_handle.write("\n\n")
    
    lp_handle.write("Subject To\n")
    for i, edge in enumerate(edges):
        a, b = edge
        lp_handle.write(" r_%d: x(%d) + x(%d) <= 1 \n" %(i, a, b))
    lp_handle.write("\n")

    lp_handle.write("Binary\n")
    for i, score in nodes:
        lp_handle.write(" x(%d) \n" %i )
    lp_handle.write("\n")
    
    lp_handle.write("End\n")

    lp_data = lp_handle.getvalue()
    lp_handle.close()

    return lp_data


def run_lp_solver(lp_data):

    print >>sys.stderr, "Write problem spec to work/lp_data"

    lpfile = "work/data.lp" # problem instance
    outfile = "work/data.out" # verbose output
    listfile = "work/data.list" # simple output

    fw = file(lpfile, "w")
    fw.write(lp_data)
    fw.close()

    try:
        proc = Popen("glpsol --cuts --fpump --lp work/data.lp -o %s -w %s" % (outfile, listfile), shell=True)
    except OSError as detail:
        print >>sys.stderr, "Error:", detail
        print >>sys.stderr, "You need to install program 'glpsol' on your path"
        print >>sys.stderr, "[http://www.gnu.org/software/glpk/]"
        sys.exit(1)

    proc.communicate()

    return listfile


def parse_lp_output(listfile):
    
    filtered_list = []

    fp = file(listfile)
    header = fp.readline()
    columns, rows = header.split()
    rows = int(rows)
    data = fp.readlines()
    # the info are contained in the last several lines
    filtered_list = [int(x) for x in data[-rows:]]
    filtered_list = [i for i, x in enumerate(filtered_list) if x==1]

    return filtered_list


def solve_lp(clusters):
    
    nodes, edges = construct_graph(clusters)

    lp_data = format_lp(nodes, edges)
    listfile = run_lp_solver(lp_data)
    filtered_list = parse_lp_output(listfile)
    
    # non-overlapping set on both axis
    filtered_clusters = [clusters[x] for x in filtered_list]

    return filtered_clusters


#___I/O routines________________________________________________________

def read_dag(filename, dag_fmt=True):

    """
## alignment a3068_scaffold_1 vs. b8_1 Alignment #1  score = 102635.0 (num aligned pairs: 2053): 
a3068_scaffold_1        scaffold_1||548785||550552||scaffold_100153.1||-1||CDS||30767256||87.37 140     140     b8_1    1||427548||427811||AT1G02210||-1||CDS||20183105||87.37     172     172     1.000000e-250   50
    """

    fp = file(filename)

    # clusters contain only bounds, point_clusters contain all points
    clusters = [] 
    
    total_lines = sum(1 for row in fp)
    print >>sys.stderr, "total lines in dag file (%d)" % total_lines
    fp.seek(0)
    row = fp.readline()
    j = 1
    while row:
        if row.strip()=="": break
        #cluster_score = int(float(row.rsplit("(", 1)[0].split()[-1]))
        row = fp.readline()
        cluster = []
        while row and row[0]!="#":
            atoms = row.strip().split("\t")
            if row.strip()=="": break
            if dag_fmt:
                ca, _, a, _, cb, _, b, _, evalue, _ = atoms
                score = int(scoringF(float(evalue)))
            else: # handle my own cluster fmt
                ca, a, cb, b, score = atoms
                score = int(float(score) * 1000) # often synteny_score 
            a, b = int(a), int(b)
            gene1, gene2 = (ca, a), (cb, b)
            cluster.append((gene1, gene2, score))
            row = fp.readline()

        clusters.append(cluster)

    return clusters


def write_chain(fw_chain, chain, clusters):

    # output the anchor points in each chain
    for c in chain:
        lines = clusters[c]
        print >>fw_chain, "# mergedcluster score %.3f" % sum(x[-1] for x in lines)
        for line in lines:
            gene1, gene2, synteny_score = line
            chr1, pos1 = gene1
            chr2, pos2 = gene2
            print >>fw_chain, \
                    "%s\t%d\t%s\t%d\t%.3f" % (chr1, pos1, chr2, pos2, synteny_score)


def write_clusters(filehandle, clusters):
    for cluster in clusters:
        cluster_score = sum(x[-1] for x in cluster)
        filehandle.write("# cluster score %.3f \n" % (cluster_score)) 
        for gene1, gene2, synteny_score in cluster:
            filehandle.write("%s\t%d\t" % gene1 )
            filehandle.write("%s\t%d\t" % gene2 )
            filehandle.write("%.3f\n" % synteny_score )


#________________________________________________________________________

if __name__ == '__main__':

    usage = "Quota synteny alignment \n" \
            "%prog [options] cluster_file "
    parser = OptionParser(usage)
    parser.add_option("-d", "--dag", dest="dag_fmt",
            action="store_true", default=False,
            help="the format of the input file is dag formatted? "\
                    "[default: %default]")
    parser.add_option("-m", "--merge", dest="merge",
            action="store_true", default=False,
            help="merge blocks first that are explained by local inversions, "\
                    "merged clusters are stored in cluster_file.merged "\
                    "[default: %default]")
    parser.add_option("-n", "--Nmax", dest="Nmax", 
            type="int", default=20,
            help="distance cutoff to determine whether two blocks are overlapping "\
                    "[default: %default gene steps] ")
    parser.add_option("-q", "--quota", dest="quota", 
            type="string", default="1:1",
            help="screen blocks to constrain mapping (sometimes useful for orthology), "\
                    "put in the format like (#subgenomes expected for genome x):"\
                    "(#subgenomes expected for genome y) "\
                    "[default: %default]")

    (options, args) = parser.parse_args()

    try:
        cluster_file = args[0]
    except:
        sys.exit(parser.print_help())

    clusters = read_cluster(cluster_file)

    Nmax = options.Nmax

    if options.merge: 
            
        chain = range(len(clusters))
        chain, clusters = recursive_merge_clusters(chain, clusters)

        merged_cluster_file = cluster_file + ".merged"
        fw = file(merged_cluster_file, "w")
        write_chain(fw, chain, clusters)
        clusters = [clusters[c] for c in chain]

    clusters = solve_lp(clusters)

    filtered_cluster_file = cluster_file + ".filtered"
    fw = file(filtered_cluster_file, "w")
    write_clusters(fw, sorted(clusters))


