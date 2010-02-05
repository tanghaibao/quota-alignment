#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
this python program can do two things:
1. merge dags by recursively merging dag bounds
2. formatted input and feed into the max independent set solver.

dag file looks like
## alignment a3068_scaffold_1 vs. b8_1 Alignment #1  score = 102635.0 (num aligned pairs: 2053): a3068_scaffold_1        
scaffold_1||548785||550552||scaffold_100153.1||-1||CDS||30767256||87.37 140     140     b8_1    1||427548||427811||AT1G02210||-1||CDS||20183105||87.37     172     172     1.000000e-250   50

"""

import sys
import pprint
import cStringIO
from mystruct import Grouper
from subprocess import Popen, PIPE
from optparse import OptionParser


def read_dag(filename):

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
        synteny_score = int(float(row.rsplit("(", 1)[0].split()[-1]))
        row = fp.readline()
        cluster = []
        while row and row[0]!="#":
            atoms = row.strip().split("\t")
            if row.strip()=="": break
            ca, _, a, _, cb, _, b, _, _ = atoms
            a, b = int(a), int(b)
            gene1, gene2 = (ca, a), (cb, b)
            cluster.append((gene1, gene2, 0))
            row = fp.readline()

        clusters.append(cluster)

    return clusters


def distance(gene1, gene2):
    chr1, pos1 = gene1
    chr2, pos2 = gene2
    if chr1 != chr2: return Nmax + 1 # this ensures un-chainable
    return abs(pos1 - pos2)


def distance_x(cluster_i, cluster_j):
    # x-sorted
    min_i_x, max_i_x = min(cluster_i)[0], max(cluster_i)[0]
    min_j_x, max_j_x = min(cluster_j)[0], max(cluster_j)[0]

    del_x1 = distance(min_i_x, max_j_x)
    del_x2 = distance(max_i_x, min_j_x)
    del_x3 = distance(min_i_x, min_j_x)
    del_x4 = distance(max_i_x, max_j_x)
    
    return min(del_x1, del_x2, del_x3, del_x4)

def distance_y(cluster_i, cluster_j):
    # y-sorted
    min_i_y = min(cluster_i, key=operator.itemgetter(1))[1]
    max_i_y = max(cluster_i, key=operator.itemgetter(1))[1]
    min_j_y = min(cluster_j, key=operator.itemgetter(1))[1]
    max_j_y = max(cluster_j, key=operator.itemgetter(1))[1]

    del_y1 = distance(min_i_y, max_j_y)
    del_y2 = distance(max_i_y, min_j_y)
    del_y3 = distance(min_i_y, min_j_y)
    del_y4 = distance(max_i_y, max_j_y)

    return min(del_y1, del_y2, del_y3, del_y4)


def merge_clusters(chain, clusters):

    # there are, in general, two kinds of breakpoints
    # those that are induced by inversions, and those by translocations
    # inversion-breakpoints are excessive breakpoints that I want to remove
    
    chain_num = len(chain)
    mergeables = Grouper() # disjoint sets of clusters that can be merged
    for j in xrange(chain_num):
        cj = chain[j]
        mergeables.join(cj, cj)
        for i in xrange(j-1, -1, -1):
            ci = chain[i]
            del_x = distance_x(clusters[ci], clusters[cj])
            if del_x > Nmax: continue 

            del_y = distance_y(clusters[ci], clusters[cj])
            if del_x + del_y > Nmax: continue
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


def range_overlap(rangei, rangej):

    amin, amax = rangei
    bmin, bmax = rangej

    # chromosomes must match
    if amin[0]!=bmin[0]: return False
    Kmax = min(amax[1]-amin[1], bmax[1]-bmin[1], Nmax/2)
    # to cope with ends that slightly overlap
    # but the overlap needs to be less than half of either segment
    return (amin[1] <= bmax[1] - Kmax) \
            and (bmin[1] <= amax[1] - Kmax)


def box_overlap(nodei, nodej):

    for rangei in nodei[:2]:
        for rangej in nodej[:2]:
            if range_overlap(rangei, rangej):
                return True

    return False


def construct_graph(clusters):

    """
    check pairwise cluster comparison, if they overlap then print edge
    (i.e. edges represent `conflict' clusters)

    """

    eclusters = [] # need ((min_x, max_x), (min_y, max_y), score) format
    for cluster in clusters:
        xlist = [a[0] for a in cluster]
        ylist = [a[1] for a in cluster]
        score = sum(a[-1] for a in cluster)
        eclusters.append(((min(xlist), max(xlist)), \
                (min(ylist), max(ylist)), score))

    # (1-based index, synteny_score)
    nodes = [(i+1, c[-1]) for i, c in enumerate(eclusters)]
    edges = []
    nnodes = len(nodes)
    for i in xrange(nnodes):
        nodei = eclusters[i]
        for j in xrange(i+1, nnodes):
            nodej = eclusters[j]
            if box_overlap(nodei, nodej): 
                edges.append((i+1, j+1))

    return nodes, edges


def format_lp(nodes, edges):

    """
    \* Problem: lp_test *\

    Maximize
     obj: x(1) + 2 x(2) + 3 x(3) + 4 x(4)

    Subject To
     r_1: x(1) + x(2) <= 1

    Bounds
     0 <= x(1) <= 1
     0 <= x(2) <= 1

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

    fw = file("work/lp_data", "w")
    fw.write(lp_data)
    fw.close()

    outfile = "work/lp_data.out"
    listfile = "work/lp_data.list"
    try:
        proc = Popen("glpsol --lp work/lp_data -o %s -w %s" % (outfile, listfile), shell=True)
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


def max_ind_set(clusters):
    
    nodes, edges = construct_graph(clusters)

    lp_data = format_lp(nodes, edges)
    listfile = run_lp_solver(lp_data)
    filtered_list = parse_lp_output(listfile)
    
    # non-overlapping set on both axis
    filtered_clusters = [clusters[x] for x in filtered_list]

    return filtered_clusters


### Output rountines

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


def write_bounds(filehandle, clusters):
    pass


if __name__ == '__main__':

    usage = "usage: %prog [options] dag_file bounds_file"
    parser = OptionParser(usage)
    parser.add_option("-r", "--relax_overlap", dest="Nmax",
            action="store", type="int", default=40,
            help="define overlap between clusters that within less than certain distance")
    parser.add_option("-m", "--merge_bounds", dest="merge",
            action="store_true", default=False,
            help="merge clusters that are explained by local inversions")
    parser.add_option("-s", "--screen_bounds", dest="screen",
            action="store_true", default=False,
            help="screen blocks to get one-on-one mapping (best orthology)")

    (options, args) = parser.parse_args()

    try:
        dag_file = args[0]
        bounds_file = args[1]
    except:
        sys.exit(parser.print_help())

    clusters = read_dag(dag_file)

    Nmax = options.Nmax

    if not options.merge: sys.exit(0)
        
    chain = range(len(clusters))
    chain, clusters = recursive_merge_clusters(chain, clusters)

    pprint.pprint(clusters)

    if not options.screen: sys.exit(0)

    clusters = max_ind_set(clusters)

    fw = file(bounds_file, "w")
    write_clusters(sorted(clusters), fw)


