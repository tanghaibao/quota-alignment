#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
this python program does the following
1. merge dags by recursively merging dag bounds
2. build conflict graph where edges represent 1d-`overlap' between blocks
3. feed the data into the linear programming solver.

"""

import os
import sys
import pprint
import cStringIO
from optparse import OptionParser

from grouper import Grouper
from cluster_utils import read_clusters, write_clusters, \
        make_range, calc_coverage
from clq_solvers import BKSolver
from lp_solvers import GLPKSolver, SCIPSolver


def range_mergeable(a, b):
    # 1-d version of box_mergeable
    a_chr, a_min, a_max = a
    b_chr, b_min, b_max = b
    # must be on the same chromosome
    if a_chr!=b_chr: return False 
    
    # make sure it is end-to-end merge, and within distance cutoff
    return (a_min <= b_max + Kmax) and \
           (b_min <= a_max + Kmax)


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
    
    # some very small segments might slip through
    # therefore we also control for the segment size
    Mmax = min(a_max-a_min, b_max-b_min, Nmax)
    # to handle ranges that slightly overlap
    return (a_min <= b_max - Mmax) and \
           (b_min <= a_max - Mmax)


def box_relaxed_overlap(boxa, boxb):

    boxa_xrange, boxa_yrange, _ = boxa
    boxb_xrange, boxb_yrange, _ = boxb

    return range_relaxed_overlap(boxa_xrange, boxb_xrange) or \
           range_relaxed_overlap(boxa_yrange, boxb_yrange)


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
            clusters[v] = list(set(clusters[v])|set(clusters[k]))

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

    # non-self
    if not (constraints_x is constraints_y):
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


def solve_lp(clusters, quota, work_dir="work", self_match=False, solver="SCIP", verbose=False):
    
    qb, qa = quota # flip it
    nodes, edges_x, edges_y = construct_conflict_graph(clusters)

    if self_match:
        clq_data = format_clq(nodes, edges_x+edges_y, qa)
        constraints = BKSolver(clq_data, work_dir=work_dir).results

        lp_data = format_lp(nodes, constraints, qa, constraints, qb)

    else:
        clq_data_x = format_clq(nodes, edges_x, qa)
        constraints_x = BKSolver(clq_data_x, work_dir=work_dir).results

        clq_data_y = format_clq(nodes, edges_y, qb)
        constraints_y = BKSolver(clq_data_y, work_dir=work_dir).results

        lp_data = format_lp(nodes, constraints_x, qa, constraints_y, qb)

    if solver=="SCIP":
        filtered_list = SCIPSolver(lp_data, work_dir=work_dir, verbose=verbose).results
        if not filtered_list:
            print >>sys.stderr, "SCIP fails... trying GLPK"
            filtered_list = GLPKSolver(lp_data, work_dir=work_dir, verbose=verbose).results
            
    elif solver=="GLPK":
        filtered_list = GLPKSolver(lp_data, work_dir=work_dir, verbose=verbose).results
        if not filtered_list:
            print >>sys.stderr, "GLPK fails... trying SCIP"
            filtered_list = SCIPSolver(lp_data, work_dir=work_dir, verbose=verbose).results
    
    # non-overlapping set on both axis
    filtered_clusters = [clusters[x] for x in filtered_list]

    return filtered_clusters



if __name__ == '__main__':

    usage = "Quota synteny alignment \n" \
            "%prog [options] cluster_file "
    parser = OptionParser(usage)

    parser.add_option("-m", "--merge", dest="merge",
            action="store_true", default=False,
            help="`block merging` procedure -- merge blocks that are close to "\
                    "each other, merged clusters are stored in cluster_file.merged "\
                    "[default: %default]")
    parser.add_option("-k", "--Kmax", dest="Kmax", 
            type="int", default=0,
            help="merge blocks that are close to each other within distance cutoff"\
                    "(cutoff for `block merging`) "\
                    "[default: %default genes] ")
    parser.add_option("-q", "--quota", dest="quota", 
            type="string", default="1:1",
            help="`quota mapping` procedure -- screen blocks to constrain mapping"\
                    " (useful for orthology), "\
                    "put in the format like (#subgenomes expected for genome X):"\
                    "(#subgenomes expected for genome Y) "\
                    "[default: %default]")
    parser.add_option("-n", "--Nmax", dest="Nmax", 
            type="int", default=40,
            help="distance cutoff to tolerate two blocks that are "\
                    "slightly overlapping (cutoff for `quota mapping`) "\
                    "[default: %default genes] ")
    parser.add_option("-f", "--self_match", dest="self_match",
            action="store_true", default=False,
            help="you might turn this on when you use this to screen paralogous blocks, "\
                    "if you have already reduced mirrored blocks into non-redundant set")
    parser.add_option("-s", "--solver", dest="solver",
            type="string", default="SCIP",
            help="use MIP solver, only SCIP or GLPK are currently implemented "\
                    "[default: %default]")
    parser.add_option("-v", dest="verbose", action="store_true", default=False,
                      help="show verbose solver output")

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

    if options.self_match and qa!=qb:
        raise Exception, "when comparing genome to itself, quota must be the same number (like 1:1, 2:2) "\
            "you have %s" % options.quota
    if qa > 12 or qb > 12:
        raise Exception, "quota %s set too loose, make quota less than 12 each" % options.quota
    quota = (qa, qb) 

    clusters = read_clusters(cluster_file)
    for cluster in clusters:
        assert len(cluster) > 0

    total_len_x, total_len_y = calc_coverage(clusters, options.self_match)

    Kmax = options.Kmax
    Nmax = options.Nmax

    # below runs `block merging`
    if options.merge: 
        chain = range(len(clusters))
        chain, clusters = recursive_merge_clusters(chain, clusters)

        merged_cluster_file = cluster_file + ".merged"
        fw = file(merged_cluster_file, "w")
        clusters = [clusters[c] for c in chain]
        write_clusters(fw, clusters)

    # below runs `quota mapping`
    if "-q" not in sys.argv and "--quota" not in sys.argv:
        sys.exit(0)

    op = os.path
    work_dir = op.join(op.dirname(op.abspath(cluster_file)), "work")
    print >>sys.stderr, "write intermediate files to", work_dir
    clusters = solve_lp(clusters, quota, work_dir=work_dir, \
            self_match=options.self_match, \
            solver=options.solver, verbose=options.verbose)

    filtered_cluster_file = cluster_file + ".filtered"
    fw = file(filtered_cluster_file, "w")
    write_clusters(fw, sorted(clusters))

    filtered_len_x, filtered_len_y = calc_coverage(clusters, options.self_match)
    if options.self_match:
        print >>sys.stderr, "coverage: %.1f%% (self-match)" % \
                (filtered_len_x*100./total_len_x)
    else:
        print >>sys.stderr, "genome x coverage: %.1f%%" % \
                (filtered_len_x*100./total_len_x)
        print >>sys.stderr, "genome y coverage: %.1f%%" % \
                (filtered_len_y*100./total_len_y)


