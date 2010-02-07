#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
this python program can do two things:
1. merge dags by recursively merging dag bounds
2. build conflict graph where edges represent 1d-`overlap' between blocks
3. feed the data into the linear programming solver.

"""

import sys
import pprint
import cStringIO
from optparse import OptionParser

from grouper import Grouper
from cluster_utils import read_clusters, write_clusters
from clq_solvers import BKSolver
from lp_solvers import GLPKSolver, SCIPSolver


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


#___merge clusters to combine inverted blocks (optional)___________________________

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


#___populate conflict graph before sending into MIP solver__________________________

def construct_conflict_graph(clusters):

    # check pairwise cluster comparison, if they overlap then mark edge as `conflict'

    eclusters = make_range(clusters)
    # (1-based index, cluster score)
    nodes = [(i+1, c[-1]) for i, c in enumerate(eclusters)]
    # represents the contraints over x-axis and y-axis
    edges_x, edges_y = [], []

    nnodes = len(nodes)
    for i in xrange(nnodes):
        for j in xrange(i+1, nnodes):
            if range_relaxed_overlap(eclusters[i][0], eclusters[j][0]): 
                edges_x.append((i, j))
            if range_relaxed_overlap(eclusters[i][1], eclusters[j][1]): 
                edges_y.append((i, j))

    return nodes, edges_x, edges_y


def format_clq(nodes, edges, q):
    
    """
    Example:

    5 0
    0 1
    0 2

    """
    clq_handle = cStringIO.StringIO()
    clq_handle.write("%d %d\n" % (len(nodes), q))
    for i, j in edges:
        clq_handle.write("%d %d\n" % (i, j))

    clq_data = clq_handle.getvalue()
    clq_handle.close()

    return clq_data


#___formulate mixed integer programming instance____________________________________

def format_lp(nodes, constraints_x, qa, constraints_y, qb):

    """
    Example:

    Maximize
     4 x1 + 2 x2 + 3 x3 + x4
    Subject To
     x1 + x2 <= 1
    End
    
    """
    lp_handle = cStringIO.StringIO()

    lp_handle.write("Maximize\n ")
    for i, score in nodes:
        lp_handle.write("+ %d x%d " % (score, i))
    lp_handle.write("\n")
    
    lp_handle.write("Subject To\n")
    for c in constraints_x:
        additions = " + ".join("x%d" % (x+1) for x in c)
        lp_handle.write(" %s <= %d\n" % (additions, qa))
    for c in constraints_y:
        additions = " + ".join("x%d" % (x+1) for x in c)
        lp_handle.write(" %s <= %d\n" % (additions, qb))

    lp_handle.write("Binary\n")
    for i, score in nodes:
        lp_handle.write(" x%d\n" %i )
    
    lp_handle.write("End\n")

    lp_data = lp_handle.getvalue()
    lp_handle.close()

    return lp_data


def solve_lp(clusters, quota, solver="SCIP"):
    
    qb, qa = quota # flip it
    nodes, edges_x, edges_y = construct_conflict_graph(clusters)
    clq_data_x = format_clq(nodes, edges_x, qa)
    constraints_x = BKSolver(clq_data_x).results

    clq_data_y = format_clq(nodes, edges_y, qb)
    constraints_y = BKSolver(clq_data_y).results

    lp_data = format_lp(nodes, constraints_x, qa, constraints_y, qb)
    if solver=="SCIP":
        filtered_list = SCIPSolver(lp_data).results
        if not filtered_list:
            print >>sys.stderr, "SCIP fails... trying GLPK"
            filtered_list = GLPKSolver(lp_data).results
            
    elif solver=="GLPK":
        filtered_list = GLPKSolver(lp_data).results
        if not filtered_list:
            print >>sys.stderr, "GLPK fails... tryinig SCIP"
            filtered_list = SCIPSolver(lp_data).results
    
    # non-overlapping set on both axis
    filtered_clusters = [clusters[x] for x in filtered_list]

    return filtered_clusters



if __name__ == '__main__':

    usage = "Quota synteny alignment \n" \
            "%prog [options] cluster_file "
    parser = OptionParser(usage)

    parser.add_option("-m", "--merge", dest="merge",
            action="store_true", default=False,
            help="merge blocks first that are explained by local inversions, "\
                    "merged clusters are stored in cluster_file.merged "\
                    "[default: %default]")
    parser.add_option("-n", "--Nmax", dest="Nmax", 
            type="int", default=20,
            help="distance cutoff to determine whether two blocks are overlapping "\
                    "[default: %default genes] ")
    parser.add_option("-q", "--quota", dest="quota", 
            type="string", default="1:1",
            help="screen blocks to constrain mapping (useful for orthology), "\
                    "put in the format like (#subgenomes expected for genome X):"\
                    "(#subgenomes expected for genome Y) "\
                    "[default: %default]")
    parser.add_option("-s", "--solver", dest="solver",
            type="string", default="SCIP",
            help="use MIP solver, only SCIP or GLPK are currently implemented "\
                    "[default: %default]")

    (options, args) = parser.parse_args()

    try:
        cluster_file = args[0]
    except:
        sys.exit(parser.print_help())

    # option sanity check
    solver_options = ("SCIP", "GLPK")
    assert options.solver in solver_options, \
            "solver must be one of %s" % solver_options

    try:
        qa, qb = options.quota.split(":")
        qa, qb = int(qa), int(qb)
    except:
        print >>sys.stderr, "quota string should be the form x:x (like 2:4, 1:3, etc.)"
        sys.exit(1)

    assert qa <= 12 and qb <= 12, \
            "quota %s set too loose, make quota less than 12 each" % options.quota
    quota = (qa, qb) 

    clusters = read_clusters(cluster_file)

    Nmax = options.Nmax

    if options.merge: 
        chain = range(len(clusters))
        chain, clusters = recursive_merge_clusters(chain, clusters)

        merged_cluster_file = cluster_file + ".merged"
        fw = file(merged_cluster_file, "w")
        clusters = [clusters[c] for c in chain]
        write_clusters(fw, clusters)

    clusters = solve_lp(clusters, quota, solver=options.solver)

    filtered_cluster_file = cluster_file + ".filtered"
    fw = file(filtered_cluster_file, "w")
    write_clusters(fw, sorted(clusters))


