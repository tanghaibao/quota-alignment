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
import itertools

from cluster_utils import read_clusters, write_clusters, \
        make_range, calc_coverage
from box_utils import get_1D_overlap, get_2D_overlap
from lp_solvers import GLPKSolver, SCIPSolver


#___merge clusters to combine inverted blocks (optional)___________________________

def merge_clusters(chain, clusters):

    # due to the problem of chaining, some chains might overlap each other
    # these need to be removed

    chain_num = len(chain)
    eclusters = make_range(clusters, extend=Kmax)
    #pprint.pprint(eclusters)

    mergeables = get_2D_overlap(chain, eclusters)

    merged_chain = []
    for mergeable in mergeables:
        merged_mother = min(mergeable)
        g = (clusters[x] for x in mergeable)
        merged_cluster = itertools.chain(*g)

        clusters[merged_mother] = list(set(merged_cluster))
        merged_chain.append(merged_mother)

    # maintain the x-sort
    [cluster.sort() for cluster in clusters]

    # nothing is merged
    updated = (len(merged_chain) != chain_num)
    return merged_chain, updated


def recursive_merge_clusters(chain, clusters):

    # as some rearrangment patterns are recursive, the extension of blocks
    # will take several iterations
    while 1: 
        print >>sys.stderr, "merging... (%d)" % len(chain)
        chain, updated = merge_clusters(chain, clusters)
        if not updated: break

    return chain, clusters


#___populate conflict graph before sending into MIP solver__________________________

def construct_conflict_graph(clusters):

    # check pairwise cluster comparison, if they overlap then mark edge as `conflict'

    eclusters = make_range(clusters, extend=-Nmax)
    # (1-based index, cluster score)
    nodes = [(i+1, c[-1]) for i, c in enumerate(eclusters)]

    eclusters_x, eclusters_y, scores = zip(*eclusters)

    # represents the contraints over x-axis and y-axis
    constraints_x = get_1D_overlap(eclusters_x)
    constraints_y = get_1D_overlap(eclusters_y)

    return nodes, constraints_x, constraints_y 


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
    nodes, constraints_x, constraints_y = construct_conflict_graph(clusters)

    if self_match:
        constraints = constraints_x | constraints_y
        lp_data = format_lp(nodes, constraints, qa, constraints, qb)
    else:
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

    from optparse import OptionParser, OptionGroup

    usage = "Quota synteny alignment \n" \
            "%prog [options] cluster_file "
    parser = OptionParser(usage)

    merge_group = OptionGroup(parser, "Merge function")
    merge_group.add_option("--merge", dest="merge",
            action="store_true", default=False,
            help="`block merging` procedure -- merge blocks that are close to "\
                    "each other, merged clusters are stored in cluster_file.merged "\
                    "[default: %default]")
    merge_group.add_option("--Km", dest="Kmax", 
            type="int", default=0,
            help="merge blocks that are close to each other within distance cutoff"\
                    "(cutoff for `block merging`) "\
                    "[default: %default genes] ")
    parser.add_option_group(merge_group)

    quota_group = OptionGroup(parser, "Quota mapping function")
    quota_group.add_option("--quota", dest="quota", 
            type="string", default="1:1",
            help="`quota mapping` procedure -- screen blocks to constrain mapping"\
                    " (useful for orthology), "\
                    "put in the format like (#subgenomes expected for genome X):"\
                    "(#subgenomes expected for genome Y) "\
                    "[default: %default]")
    quota_group.add_option("--Dm", dest="Nmax", 
            type="int", default=40,
            help="distance cutoff to tolerate two blocks that are "\
                    "slightly overlapping (cutoff for `quota mapping`) "\
                    "[default: %default units (gene dist or bp dist, depending on the input)]")
    parser.add_option_group(quota_group)

    other_group = OptionGroup(parser, "Other options")
    other_group.add_option("--self_match", dest="self_match",
            action="store_true", default=False,
            help="you might turn this on when you use this to screen paralogous blocks, "\
                 "especially if you have reduced mirrored blocks into non-redundant set")
    other_group.add_option("--solver", dest="solver",
            type="string", default="SCIP",
            help="use MIP solver, only SCIP or GLPK are currently implemented "\
                 "[default: %default]")
    other_group.add_option("--verbose", dest="verbose", action="store_true", default=False,
                      help="show verbose solver output")
    parser.add_option_group(other_group)

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
        raise Exception, "when comparing genome to itself, quota must be the same number " \
                "(like 1:1, 2:2) you have %s" % options.quota
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
    if "--quota" not in sys.argv:
        sys.exit(0)

    op = os.path
    work_dir = op.join(op.dirname(op.abspath(cluster_file)), "work")

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


