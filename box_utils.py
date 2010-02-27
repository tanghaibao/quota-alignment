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


def get_1D_overlap(eclusters, depth=1):
    """
    Naive implementation for 1D overlapping
    returns cliques of ids that are in conflict

    """

    overlap_set = set() 
    active = set()

    ends = []
    for i, (chr, left, right) in enumerate(eclusters):
        ends.append((chr, left, 0, i))
        ends.append((chr, right, 1, i))
    ends.sort()

    chr_last = ""
    for chr, pos, left_right, i in ends:
        if chr != chr_last: active.clear()
        if left_right==0: active.add(i) 
        else: active.remove(i)    

        if len(active) > depth:
            overlap_set.add(tuple(sorted(active)))

        chr_last = chr

    return list(overlap_set)


def get_2D_overlap(chain, eclusters):
    """
    Implements a sweep line algorithm when n is large:
    assume block has x_ends, and y_ends for the bounds

    1. sort x_ends, and take a sweep line to scan the x_ends
    2. if this is left end, test y-axis intersection of current block with `active` set;
       also put this block in the `active` set
    3. if this is right end, remove block from the `active` set

    """

    mergeables = Grouper()
    active = set()

    x_ends = []
    for i, (range_x, range_y, score) in enumerate(eclusters):
        chr, left, right = range_x
        x_ends.append((chr, left, 0, i))
        x_ends.append((chr, right, 1, i))
    x_ends.sort()

    chr_last = ""
    for chr, pos, left_right, i in x_ends:
        if chr != chr_last: active.clear()
        if left_right==0: 
            active.add(i) 
            for x in active:
                # check y-overlap
                if range_overlap(eclusters[x][1], eclusters[i][1]):
                    mergeables.join(x, i)
        else: # right end
            active.remove(i) 

        chr_last = chr

    return mergeables


if __name__ == '__main__':
    pass
