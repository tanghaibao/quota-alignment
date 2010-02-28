#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This script uses bx-python tool <http://bitbucket.org/james_taylor/bx-python/>
to parse the .maf file in order to convert the blocks into .qa format.
"""

import os
import sys
from bx.align import maf
from bx import interval_index_file

from cluster_utils import read_clusters


def alignment_to_cluster(alignment):
    """
    From the pairwise alignment into a cluster format, providing two (fake) anchors
    to mark the block boundary
    """
    region_a, region_b = alignment
    chr_a, start_a, stop_a, strand_a, score_a = region_a
    chr_b, start_b, stop_b, strand_b, score_b = region_b
    if strand_a!=strand_b: 
        start_b, stop_b = stop_b, start_b

    cluster = []
    cluster.append(((chr_a, start_a), (chr_b, start_b), score_a))
    cluster.append(((chr_a, stop_a), (chr_b, stop_b), 0))

    return cluster


def get_clusters(maf_file):
    """
    Called in cluster_utils, from maf_file to get clusters to convert in .qa format
    """
    base, ext = os.path.splitext(maf_file)

    fp = file(maf_file)
    reader = maf.Reader(fp)

    clusters = []
    for rec in reader:
        alignment = []
        for c in rec.components:
            chr, left, right, strand, weight = c.src, c.forward_strand_start, \
                    c.forward_strand_end, c.strand, rec.score
            alignment.append((chr, left, right, strand, weight))

        clusters.append(alignment_to_cluster(alignment))

    fp.close()

    return clusters 


def screen_maf(qa_file, maf_file):
    """
    Screen the .maf file based on the cluster info in the qa_file
    """
    clusters = read_clusters(qa_file)
    filtered_maf = maf_file + ".filtered"

    screened_alignments = set() 
    for cluster in clusters:
        for anchor in cluster:
            score = anchor[-1]
            if score!=0:
                screened_alignments.add(anchor)

    fp = file(maf_file)
    reader = maf.Reader(fp)

    fw = file(filtered_maf, "w")
    writer = maf.Writer(fw)

    for rec in reader:
        alignment = []
        for c in rec.components:
            chr, left, right, strand, score = c.src, c.forward_strand_start, \
                    c.forward_strand_end, c.strand, rec.score
            alignment.append((chr, left, right, strand, score))

        cluster = alignment_to_cluster(alignment)
        if cluster[0] in screened_alignments:
            writer.write(rec)

    fp.close()

    print >>sys.stderr, "write (%d) alignments to '%s'" % \
            (len(screened_alignments), filtered_maf)


if __name__ == '__main__':
    
    from optparse import OptionParser

    usage = "Use information in qa_file to screen maf_file\n" \
            "%prog [options] qa_file maf_file"
    parser = OptionParser(usage)

    (options, args) = parser.parse_args()
    try:
        qa_file, maf_file = args
    except:
        sys.exit(parser.print_help())

    screen_maf(qa_file, maf_file)

