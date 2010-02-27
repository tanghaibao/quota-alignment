#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This script implements a high performance algorithm for finding 2d intersecting rectangles, assume a list of boxes with (x1, x2, y1, y2) showing the extent:

1. take the x-ends and sort, then use a sweep line to scan the x-ends
2. if it is left end, test the y-axis intersection with an `active` set, also update `active`
3. if it is right end, remove the block from the `active` set

"""

from grouper import Grouper

def range_mergeable(a, b):
    # 1-d version of box_mergeable
    a_chr, a_min, a_max = a
    b_chr, b_min, b_max = b
    # must be on the same chromosome
    if a_chr!=b_chr: return False 
    
    # make sure it is end-to-end merge, and within distance cutoff
    return (a_min <= b_max) and \
           (b_min <= a_max)


def box_mergeable(boxa, boxb):

    boxa_xrange, boxa_yrange, _ = boxa
    boxb_xrange, boxb_yrange, _ = boxb

    return range_mergeable(boxa_xrange, boxb_xrange) and \
           range_mergeable(boxa_yrange, boxb_yrange)


def get_2Doverlap(chain, eclusters):
    """
    Naive implementation -- check all pairwise combinations, O(n^2) complexity

    """

    chain_num = len(chain)

    mergeables = Grouper() # disjoint sets of clusters that can be merged
    # check all pairwise combinations
    for i in xrange(chain_num):
        ci = chain[i]
        mergeables.join(ci)
        for j in xrange(i+1, chain_num):
            cj = chain[j]
            if box_mergeable(eclusters[ci], eclusters[cj]):
                mergeables.join(ci, cj)

    return mergeables


# Faster version
def get_2Doverlap_fast(chain, eclusters, Kmax):
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
