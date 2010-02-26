#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from grouper import Grouper

def range_mergeable(a, b, Kmax):
    # 1-d version of box_mergeable
    a_chr, a_min, a_max = a
    b_chr, b_min, b_max = b
    # must be on the same chromosome
    if a_chr!=b_chr: return False 
    
    # make sure it is end-to-end merge, and within distance cutoff
    return (a_min <= b_max + Kmax) and \
           (b_min <= a_max + Kmax)


def box_mergeable(boxa, boxb, Kmax):

    boxa_xrange, boxa_yrange, _ = boxa
    boxb_xrange, boxb_yrange, _ = boxb

    return range_mergeable(boxa_xrange, boxb_xrange, Kmax) and \
           range_mergeable(boxa_yrange, boxb_yrange, Kmax)

def get_2doverlap(chain, eclusters, Kmax):

    chain_num = len(chain)

    mergeables = Grouper() # disjoint sets of clusters that can be merged
    # check all pairwise combinations
    for i in xrange(chain_num):
        ci = chain[i]
        mergeables.join(ci)
        for j in xrange(i+1, chain_num):
            cj = chain[j]
            if box_mergeable(eclusters[ci], eclusters[cj], Kmax):
                mergeables.join(ci, cj)

    return mergeables


if __name__ == '__main__':
    pass
