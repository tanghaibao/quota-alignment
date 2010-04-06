#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog blast_file --qbed query.bed --sbed subject.bed

visualize the blast_file in a dotplot
"""

import os.path as op
import itertools
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from bed_utils import Bed, BlastLine


def get_breaks(bed):
    # get chromosome break positions
    simple_bed = bed.get_simple_bed()
    for seqid, ranks in itertools.groupby(simple_bed, key=lambda x:x[0]):
        ranks = list(ranks)
        # chromosome, extent of the chromosome
        yield seqid, ranks[0][1], ranks[-1][1]


def dotplot(blast_file, qbed, sbed, image_name):

    blast_fh = file(blast_file)
    blasts = [BlastLine(line) for line in blast_fh]
    seen = set()

    qorder = qbed.get_order()
    sorder = sbed.get_order()

    data = []
    for b in blasts:
        query, subject = b.query, b.subject
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
    ax.scatter(x, y, c='k', s=.05, lw=0, alpha=.9)

    xlim = (0, len(qbed))
    ylim = (0, len(sbed))

    xchr_labels, ychr_labels = [], []
    # plot the chromosome breaks
    for (seqid, beg, end) in get_breaks(qbed):
        xchr_labels.append((seqid, (beg + end)/2))
        ax.plot([beg, beg], ylim, "g-", alpha=.5)

    for (seqid, beg, end) in get_breaks(sbed):
        ychr_labels.append((seqid, (beg + end)/2))
        ax.plot(xlim, [beg, beg], "g-", alpha=.5)

    # plot the chromosome labels
    for label, pos in xchr_labels:
        pos = .1 + pos * .8/xlim[1]
        root.text(pos, .91, r"%s" % label, color="b",
            size=9, alpha=.5, rotation=45)

    for label, pos in ychr_labels:
        pos = .1 + pos * .8/ylim[1]
        root.text(.91, pos, r"%s" % label, color="b",
            size=9, alpha=.5, ha="left", va="center")

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # i always like the latex font
    _ = lambda x: r"$\rm{%s}$" % x.replace("_", " ").replace(" ", r"\ ")
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
    print >>sys.stderr, "print image to %s" % image_name
    plt.savefig(image_name, dpi=1000)


if __name__ == "__main__":

    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed",
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed",
            help="path to sbed")

    (options, args) = parser.parse_args()

    if not (len(args) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    qbed = Bed(options.qbed)
    sbed = Bed(options.sbed)

    blast_file = args[0]

    image_name = op.splitext(blast_file)[0] + ".png"
    dotplot(blast_file, qbed, sbed, image_name)

