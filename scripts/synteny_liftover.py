#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog anchor_file blast_file [options] 

typical use for this program is given an anchor list (for example, syntenic genes), choose from the second blast_file for which the pairs are close to the anchors.

Anchor list should have the following format, each row defines a pair:

    geneA geneB
    geneC geneD

Use KD-tree for processing distance query.

"""

import sys
from bed_utils import Bed
import numpy as np
import collections
from scipy.spatial import cKDTree
from blast_to_raw import BlastLine


def main(anchor_file, blast_file, options):

    qbed_file, sbed_file = options.qbed, options.sbed
    # is this a self-self blast?
    is_self = (qbed_file == sbed_file)
    if is_self:
        print >>sys.stderr, "... looks like a self-self BLAST to me"

    print >>sys.stderr, "read annotation files %s and %s" % (qbed_file, sbed_file)
    qbed = Bed(qbed_file)
    sbed = Bed(sbed_file)
    qorder = qbed.get_order()
    sorder = sbed.get_order()
    _ = lambda x: x.rsplit(".", 1)[0]

    fp = file(blast_file)
    print >>sys.stderr, "read BLAST file %s (total %d lines)" % \
            (blast_file, sum(1 for line in fp))
    fp.seek(0)
    blasts = sorted([BlastLine(line) for line in fp], \
            key=lambda b: b.score, reverse=True)
    filtered_blasts = []
    seen = set()
    for b in blasts:
        query, subject = _(b.query), _(b.subject)
        if query not in qorder or subject not in sorder: continue
        qi, q = qorder[query]
        si, s = sorder[subject]

        if is_self and qi > si:
            # remove redundant a<->b to one side when doing self-self BLAST
            query, subject = subject, query
            qi, si = si, qi
            q, s = s, q

        key = query, subject
        if key in seen: continue
        seen.add(key)
        b.query, b.subject = key

        b.qi, b.si = qi, si
        b.qseqid, b.sseqid = q.seqid, s.seqid
        
        filtered_blasts.append(b)


    all_anchors = collections.defaultdict(list)
    fp = file(anchor_file)
    for row in fp:
        if row[0]=='#': continue
        a, b = row.split()
        if a not in qorder or b not in sorder: continue
        qi, q = qorder[a]
        si, s = sorder[b]
        all_anchors[(q.seqid, s.seqid)].append((qi, si))

    # grouping the hits based on chromosome pair for sending in find_nearby
    all_hits = collections.defaultdict(list)
    for b in filtered_blasts:
        all_hits[(b.qseqid, b.sseqid)].append((b.qi, b.si))

    # select hits that are close to the anchor list
    j = 0
    fw = sys.stdout
    for chr_pair in sorted(all_hits.keys()):
        hits = np.array(all_hits[chr_pair])
        anchors = np.array(all_anchors[chr_pair])

        print >>sys.stderr, chr_pair, len(anchors)
        if len(anchors)==0: continue
        tree = cKDTree(anchors, leafsize=16)
        #print tree.data
        dists, idxs = tree.query(hits, p=1, distance_upper_bound=options.dist)
        #print [(d, idx) for (d, idx) in zip(dists, idxs) if idx!=tree.n]

        for i, (dd, idx) in enumerate(zip(dists, idxs)):
            if dd==0: continue # same anchors
            if idx!=tree.n:
                qi, si = hits[i]
                query, subject = qbed[qi]["accn"], sbed[si]["accn"]
                print >>fw, "\t".join((query, subject, "lifted"))
                j+=1
    
    print >>sys.stderr, j, "new pairs found"



if __name__ == '__main__':
    
    import optparse

    p = optparse.OptionParser(__doc__)

    p.add_option("--qbed", dest="qbed", help="path to qbed")
    p.add_option("--sbed", dest="sbed", help="path to sbed")

    p.add_option("--dist", dest="dist",
            default=10, type="int", 
            help="the extent of flanking regions to search [default: %default]") 

    options, files = p.parse_args()

    if not (len(files) == 2 and options.qbed and options.sbed):
        sys.exit(p.print_help())

    anchor_file, blast_file = files

    main(anchor_file, blast_file, options)

