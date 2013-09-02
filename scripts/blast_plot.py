#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog blast_file --qbed query.bed --sbed subject.bed

visualize the blast_file in a dotplot
"""

import os.path as op
import itertools
import sys
import collections
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

from bed_utils import Bed, BlastLine
from qa_plot import get_breaks, get_len
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from grouper import Grouper


def score(cluster):
    x, y = zip(*cluster)
    return min(len(set(x)), len(set(y)))


def single_linkage(points, xdist, ydist, N):

    # This is the core single linkage algorithm
    # this behaves in O(n) complexity: we iterate through the pairs, for each pair
    # we look back on the adjacent pairs to find links

    clusters = Grouper()
    n = len(points)
    points.sort()
    for i in xrange(n):
        for j in xrange(i-1, -1, -1):
            # x-axis distance
            del_x = points[i][0]-points[j][0]
            if del_x > xdist: break
            # y-axis distance
            del_y = points[i][1]-points[j][1]
            if abs(del_y) > ydist: continue
            # otherwise join
            clusters.join(points[i], points[j])
    clusters = [cluster for cluster in list(clusters) if score(cluster)>=N]
    return clusters


def batch_linkage(points, qbed, sbed, xdist=20, ydist=20, N=6):

    # this runs single_linkage() per chromosome pair
    chr_pair_points = collections.defaultdict(list)
    for qi, si in points:
        q, s = qbed[qi], sbed[si]
        chr_pair_points[(q.seqid, s.seqid)].append((qi, si))

    clusters = []
    for points in chr_pair_points.values():
        clusters.extend(single_linkage(points, xdist, ydist, N))

    return clusters


def draw_box(clusters, ax, color="b"):

    for cluster in clusters:
        xrect, yrect = zip(*cluster)
        xmin, xmax, ymin, ymax = min(xrect), max(xrect), \
                                min(yrect), max(yrect)
        ax.add_patch(Rectangle((xmin, ymin), xmax-xmin, ymax-ymin,\
                                ec=color, fc='y', alpha=.5))
        #ax.plot(xrect, yrect, 'r.', ms=3)


def dotplot(blast_file, qbed, sbed, image_name, is_self=False,
            synteny=False, bpscale=False):

    blast_fh = file(blast_file)
    blasts = [BlastLine(line) for line in blast_fh]
    seen = set()

    qorder = qbed.get_order()
    sorder = sbed.get_order()

    data = []
    for b in blasts:
        query, subject = b.query, b.subject
        if bpscale:
            qi = (b.qstart + b.qstop) / 2
            si = (b.sstart + b.sstop) / 2
        else:
            if query not in qorder or subject not in sorder: continue
            key = query, subject
            if key in seen: continue
            seen.add(key)

            qi, q = qorder[query]
            si, s = sorder[subject]

        data.append((qi, si))

    fig = plt.figure(1,(8,8))
    root = fig.add_axes([0,0,1,1]) # the whole canvas
    ax = fig.add_axes([.1,.1,.8,.8]) # the dot plot

    x, y = zip(*data)
    ax.scatter(x, y, c='k', s=1, lw=0, alpha=.9)
    print >> sys.stderr, "A total of {0} dots to output.".format(len(data))

    if synteny:
        clusters = batch_linkage(data, qbed, sbed)
        draw_box(clusters, ax)

    xsize = sum([s for (seqid, s) in get_len(qbed, bpscale=bpscale)])
    ysize = sum([s for (seqid, s) in get_len(sbed, bpscale=bpscale)])
    xlim = (0, xsize)
    ylim = (0, ysize)

    xchr_labels, ychr_labels = [], []
    ignore = True # tag to mark whether to plot chromosome name (skip small ones)
    ignore_size = 100
    # plot the chromosome breaks
    for (seqid, beg, end) in get_breaks(qbed, bpscale=bpscale):
        ignore = abs(end-beg) < ignore_size
        xchr_labels.append((seqid, (beg + end)/2, ignore))
        ax.plot([beg, beg], ylim, "g-", alpha=.5)

    for (seqid, beg, end) in get_breaks(sbed, bpscale=bpscale):
        ignore = abs(end-beg) < ignore_size
        ychr_labels.append((seqid, (beg + end)/2, ignore))
        ax.plot(xlim, [beg, beg], "g-", alpha=.5)

    # plot the chromosome labels
    for label, pos, ignore in xchr_labels:
        pos = .1 + pos * .8/xlim[1]
        if not ignore:
            root.text(pos, .93, r"%s" % label, color="b",
                size=9, alpha=.5, rotation=45)

    for label, pos, ignore in ychr_labels:
        pos = .1 + pos * .8/ylim[1]
        if not ignore:
            root.text(.91, pos, r"%s" % label, color="b",
                size=9, alpha=.5, ha="left", va="center")

    # create a diagonal to separate mirror image for self comparison
    if is_self:
        ax.plot(xlim, ylim, 'm-', alpha=.5, lw=2)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # i always like the latex font
    _ = lambda x: r"$\mathsf{%s}$" % x.replace("_", " ").replace(" ", r"\ ")
    to_ax_label = lambda fname: _(op.basename(fname).split(".")[0])

    # add genome names
    ax.set_xlabel(to_ax_label(qbed.filename))
    ax.set_ylabel(to_ax_label(sbed.filename))

    # beautify the numeric axis
    [tick.set_visible(False) for tick in ax.get_xticklines() + ax.get_yticklines()]
    formatter = ticker.FuncFormatter(lambda x, pos: r"$%d$" % x)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    plt.setp(ax.get_xticklabels() + ax.get_yticklabels(), color='gray', size=10)

    root.set_axis_off()
    print >>sys.stderr, "print image to `%s`" % image_name
    plt.savefig(image_name, dpi=1000)


if __name__ == "__main__":

    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed",
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed",
            help="path to sbed")
    parser.add_option("--synteny", dest="synteny",
            default=False, action="store_true",
            help="run a fast synteny scan and display synteny blocks")
    parser.add_option("--bpscale", default=False, action="store_true",
            help="Use actual bp position in bed file [default: %default]")
    parser.add_option("--format", dest="format", default="png",
            help="generate image of format (png, pdf, ps, eps, svg, etc.)"
            "[default: %default]")

    (options, args) = parser.parse_args()

    if not (len(args) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    is_self = False
    if options.qbed==options.sbed:
        print >>sys.stderr, "Looks like this is self-self comparison"
        is_self = True

    qbed = Bed(options.qbed)
    sbed = Bed(options.sbed)
    synteny = options.synteny

    blast_file = args[0]

    image_name = op.splitext(blast_file)[0] + "." + options.format
    dotplot(blast_file, qbed, sbed, image_name, is_self=is_self,
            synteny=synteny, bpscale=options.bpscale)
