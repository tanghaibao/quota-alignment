#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This script implements algorithm for finding intersecting rectangles, both on the 2d dotplot,
and on the projection onto axis

"""

from grouper import Grouper

def range_overlap(a, b):
    # 1-d version of box_mergeable
    a_chr, a_min, a_max = a
    b_chr, b_min, b_max = b
    # must be on the same chromosome
    if a_chr!=b_chr: return False 
    
    # make sure it is end-to-end merge, and within distance cutoff
    return (a_min <= b_max) and \
           (b_min <= a_max)


def box_overlap(boxa, boxb):

    boxa_xrange, boxa_yrange, _ = boxa
    boxb_xrange, boxb_yrange, _ = boxb

    return range_overlap(boxa_xrange, boxb_xrange) and \
           range_overlap(boxa_yrange, boxb_yrange)


def get_1D_overlap(eclusters):
    """
    Naive implementation for 1D overlapping
    returns pairs of ids that are in conflict

    """

    overlap_pairs = []

    nnodes = len(eclusters) 
    for i in xrange(nnodes):
        for j in xrange(i+1, nnodes):
            if range_overlap(eclusters[i], eclusters[j]): 
                overlap_pairs.append((i, j))

    return overlap_pairs


def get_2D_overlap(chain, eclusters):
    """
    Naive implementation for 2d overlapping
    check all pairwise combinations, O(n^2) complexity
    returns Grouper() object

    """

    chain_num = len(chain)

    mergeables = Grouper() # disjoint sets of clusters that can be merged
    # check all pairwise combinations
    for i in xrange(chain_num):
        ci = chain[i]
        mergeables.join(ci)
        for j in xrange(i+1, chain_num):
            cj = chain[j]
            if box_overlap(eclusters[ci], eclusters[cj]):
                mergeables.join(ci, cj)

    return mergeables


# Faster version
def get_2D_overlap_fast(chain, eclusters, Kmax):
    """
    Implements a sweep line algorithm when n is large:
    assume block has x_ends, and y_ends for the bounds

    1. sort x_ends, and take a sweep line to scan the x_ends
    2. if this is left end, test y-axis intersection of current block with `active` set;
       also put this block in the `active` set
    3. if this is right end, remove block from the `active` set

    """

    mergeables = Grouper()

    return mergeables



if __name__ == '__main__':
    pass
