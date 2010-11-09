#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This python program does the following:
1. merge 2D-overlapping blocks 
2. build constraints that represent 1D-overlap among blocks
3. feed the data into the linear programming solver
"""

import os
import sys
import cStringIO
import itertools

from cluster_utils import read_clusters, write_clusters, \
        make_range, calc_coverage
from box_utils import get_1D_overlap, get_2D_overlap
from lp_solvers import GLPKSolver, SCIPSolver


def merge_clusters(chain, clusters, Dmax=0, min_size=0):
    """
    Due to the problem of chaining, some chains might overlap each other
    these need to be merged 
    """
    eclusters = make_range(clusters, extend=Dmax)

    mergeables = get_2D_overlap(chain, eclusters)

    merged_chain = []
    for mergeable in mergeables:
        merged_mother = min(mergeable)
        g = (clusters[x] for x in mergeable)
        merged_cluster = itertools.chain(*g)
        merged_cluster = list(set(merged_cluster))

        clusters[merged_mother] = merged_cluster
        if len(merged_cluster) >= min_size:
            merged_chain.append(merged_mother)

    # maintain the x-sort
    [cluster.sort() for cluster in clusters]

    print >>sys.stderr, "merging... (%d->%d)" % (len(chain), len(merged_chain))

    return merged_chain


def get_constraints(clusters, quota=(1,1), Nmax=0):
    """
    Check pairwise cluster comparison, if they overlap then mark edge as conflict
    """
    qa, qb = quota
    eclusters = make_range(clusters, extend=-Nmax)
    # (1-based index, cluster score)
    nodes = [(i+1, c[-1]) for i, c in enumerate(eclusters)]

    eclusters_x, eclusters_y, scores = zip(*eclusters)

    # represents the contraints over x-axis and y-axis
    constraints_x = get_1D_overlap(eclusters_x, qa)
    constraints_y = get_1D_overlap(eclusters_y, qb)

    return nodes, constraints_x, constraints_y 


def format_lp(nodes, constraints_x, qa, constraints_y, qb):
    """
    Maximize
     4 x1 + 2 x2 + 3 x3 + x4
    Subject To
     x1 + x2 <= 1
    End
    """
    lp_handle = cStringIO.StringIO()

    lp_handle.write("Maximize\n ")
    records = 0
    for i, score in nodes:
        lp_handle.write("+ %d x%d " % (score, i))
        # SCIP does not like really long string per row
        records += 1
        if records%10==0: lp_handle.write("\n")
    lp_handle.write("\n")
    
    num_of_constraints = 0
    lp_handle.write("Subject To\n")
    for c in constraints_x:
        additions = " + ".join("x%d" % (x+1) for x in c)
        lp_handle.write(" %s <= %d\n" % (additions, qa))
    num_of_constraints += len(constraints_x)

    # non-self
    if not (constraints_x is constraints_y):
        for c in constraints_y:
            additions = " + ".join("x%d" % (x+1) for x in c)
            lp_handle.write(" %s <= %d\n" % (additions, qb))
        num_of_constraints += len(constraints_y)

    print >>sys.stderr, "number of variables (%d), number of constraints (%d)" % \
            (len(nodes), num_of_constraints)

    lp_handle.write("Binary\n")
    for i, score in nodes:
        lp_handle.write(" x%d\n" %i )
    
    lp_handle.write("End\n")

    lp_data = lp_handle.getvalue()
    lp_handle.close()

    return lp_data


def solve_lp(clusters, quota, work_dir="work", Nmax=0, 
        self_match=False, solver="SCIP", verbose=False):
    """
    Solve the formatted LP instance
    """
    qb, qa = quota # flip it
    nodes, constraints_x, constraints_y = get_constraints(clusters, (qa, qb), Nmax=Nmax)

    if self_match:
        constraints_x = constraints_y = constraints_x | constraints_y

    lp_data = format_lp(nodes, constraints_x, qa, constraints_y, qb)

    if solver=="SCIP":
        filtered_list = SCIPSolver(lp_data, work_dir, verbose=verbose).results
        if not filtered_list:
            print >>sys.stderr, "SCIP fails... trying GLPK"
            filtered_list = GLPKSolver(lp_data, work_dir, verbose=verbose).results
            
    elif solver=="GLPK":
        filtered_list = GLPKSolver(lp_data, work_dir, verbose=verbose).results
        if not filtered_list:
            print >>sys.stderr, "GLPK fails... trying SCIP"
            filtered_list = SCIPSolver(lp_data, work_dir, verbose=verbose).results
    
    # non-overlapping set on both axis
    filtered_clusters = [clusters[x] for x in filtered_list]

    return filtered_clusters



if __name__ == '__main__':

    from optparse import OptionParser, OptionGroup

    usage = "Quota synteny alignment \n" \
            "%prog [options] qa_file "
    parser = OptionParser(usage)

    merge_group = OptionGroup(parser, "Merge function")
    merge_group.add_option("--merge", dest="merge",
            action="store_true", default=False,
            help="`block merging` procedure -- merge blocks that are close to "\
                    "each other, merged clusters are stored in qa_file.merged "\
                    "[default: %default]")
    merge_group.add_option("--Dm", dest="Dmax", 
            type="int", default=0,
            help="merge blocks that are close to each other within distance cutoff "\
                    "(cutoff for `block merging`) "\
                    "[default: %default units (gene or bp dist)] ")
    merge_group.add_option("--min_size", dest="min_size",
            type="int", default=1,
            help="keep blocks that contain more than certain number of anchors "\
                    "[default: %default anchor points] ")
    parser.add_option_group(merge_group)

    quota_group = OptionGroup(parser, "Quota mapping function")
    quota_group.add_option("--quota", dest="quota", 
            type="string", default=None,
            help="`quota mapping` procedure -- screen blocks to constrain mapping"\
                    " (useful for orthology), "\
                    "put in the format like (#subgenomes expected for genome X):"\
                    "(#subgenomes expected for genome Y) "\
                    "[default: %default]")
    quota_group.add_option("--Nm", dest="Nmax", 
            type="int", default=40,
            help="distance cutoff to tolerate two blocks that are "\
                    "slightly overlapping (cutoff for `quota mapping`) "\
                    "[default: %default units (gene or bp dist)]")
    parser.add_option_group(quota_group)

    supported_solvers = ("SCIP", "GLPK")
    other_group = OptionGroup(parser, "Other options")
    other_group.add_option("--format", dest="format", default="qa",
            help="one of ('qa', 'raw'). if 'raw' each line is treated as a cluster and should"
                           " be used with --merge .\n[default: %default]")
    other_group.add_option("--self", dest="self_match",
            action="store_true", default=False,
            help="you might turn this on when screening paralogous blocks, "\
                 "esp. if you have reduced mirrored blocks into non-redundant set")
    other_group.add_option("--solver", dest="solver",
            default="SCIP", choices=supported_solvers,
            help="use MIP solver, must be one of %s " % (supported_solvers,) +\
                 "[default: %default]")
    other_group.add_option("--verbose", dest="verbose", action="store_true", 
            default=False, help="show verbose solver output")
    parser.add_option_group(other_group)

    (options, args) = parser.parse_args()

    try:
        qa_file = args[0]
    except:
        sys.exit(parser.print_help())

    # sanity check for the quota
    if options.quota:
        try:
            qa, qb = options.quota.split(":")
            qa, qb = int(qa), int(qb)
        except:
            print >>sys.stderr, "quota string should be the form x:x (2:4, 1:3, etc.)"
            sys.exit(1)

        if options.self_match and qa!=qb:
            raise Exception, "when comparing genome to itself, " \
                    "quota must be the same number " \
                    "(like 1:1, 2:2) you have %s" % options.quota
        if qa > 12 or qb > 12:
            raise Exception, "quota %s too loose, make it <=12 each" % options.quota
        quota = (qa, qb) 

    self_match = options.self_match

    clusters = read_clusters(qa_file, fmt=options.format)
    for cluster in clusters:
        assert len(cluster) > 0

    # below runs `block merging`
    if options.merge: 
        chain = range(len(clusters))
        chain = merge_clusters(chain, clusters, Dmax=options.Dmax, min_size=options.min_size)

        merged_qa_file = qa_file + ".merged"
        fw = file(merged_qa_file, "w")
        clusters = [clusters[c] for c in chain]
        write_clusters(fw, clusters)

    total_len_x, total_len_y = calc_coverage(clusters, self_match=self_match)

    if not options.quota:
        sys.exit(0)

    # below runs `quota mapping`
    op = os.path
    work_dir = op.join(op.dirname(op.abspath(qa_file)), "work")

    clusters = solve_lp(clusters, quota, work_dir=work_dir, \
            Nmax=options.Nmax, self_match=self_match, \
            solver=options.solver, verbose=options.verbose)

    filtered_qa_file = qa_file + ".filtered"
    fw = file(filtered_qa_file, "w")

    write_clusters(fw, sorted(clusters))

    filtered_len_x, filtered_len_y = calc_coverage(clusters, self_match=self_match)
    if self_match:
        print >>sys.stderr, "coverage: %.1f%% (self-match)" % \
                (filtered_len_x*100./total_len_x)
    else:
        print >>sys.stderr, "genome X coverage: %.1f%%" % \
                (filtered_len_x*100./total_len_x)
        print >>sys.stderr, "genome Y coverage: %.1f%%" % \
                (filtered_len_y*100./total_len_y)


