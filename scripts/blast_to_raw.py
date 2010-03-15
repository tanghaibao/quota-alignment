#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog blast_file --qbed query.bed --sbed subject.bed

accepts .bed format: <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>
or .flat format: <http://github.com/brentp/flatfeature/>
and a blast file.

local dup filter:
if the input is query.bed and subject.bed, the script files query.localdups and subject.localdups are created containing the parent|offspring dups, as inferred by subjects hitting the same query or queries hitting the same subject.

top N filter:
just returns the best N hits given a query or subject

cscore filter:
see supplementary info for sea anemone genome paper <http://www.sciencemag.org/cgi/content/abstract/317/5834/86>, formula below
cscore(A,B) = score(A,B) / max(best score for A, best score for B)

Finally a .raw file (which is the input for the quota-align pipeline <http://github.com/tanghaibao/quota-alignment/>) is created
"""

import sys
import os.path as op
import collections
import itertools

from math import log10
from bed_utils import Bed, BlastLine
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from grouper import Grouper


# helper functions in the BLAST filtering to get rid alternative splicings
def gene_name(st):
    # this is ugly, but different annotation groups are inconsistent
    # with how the alternative splicings are named;
    # mostly it can be done by removing the suffix
    # except for papaya (evm...) and maize (somewhat complicated)
    if st.startswith("ev"):
        return st
    if st.startswith("Os"):
        return st.rsplit("-",1)[0]
    return st.rsplit(".", 1)[0]


# get the gene order given a Bed or Flat object
get_order = lambda bed: dict((f['accn'], (i, f)) for (i, f) in enumerate(bed))


def main(blast_file, options):

    is_flat_fmt = options.qbed.endswith(".flat")
    qbed_file, sbed_file = options.qbed, options.sbed

    # is this a self-self blast?
    is_self = (qbed_file == sbed_file)
    if is_self:
        print >>sys.stderr, "... looks like a self-self BLAST to me"
    
    tandem_Nmax = options.tandem_Nmax
    top_N = options.top_N
    cscore = options.cscore

    print >>sys.stderr, "read annotation files %s and %s" % (qbed_file, sbed_file)
    if is_flat_fmt:
        from flatfeature import Flat
        qbed = Flat(qbed_file)
        sbed = Flat(sbed_file)
    else:
        qbed = Bed(qbed_file)
        sbed = Bed(sbed_file)

    qorder = get_order(qbed) 
    sorder = get_order(sbed) 

    fp = file(blast_file)
    print >>sys.stderr, "read BLAST file %s (total %d lines)" % \
            (blast_file, sum(1 for line in fp))
    fp.seek(0)
    blasts = sorted([BlastLine(line) for line in fp], \
            key=lambda b: b.score, reverse=True)

    filtered_blasts = []
    seen = set() 
    for b in blasts:
        query, subject = b.query, b.subject
        query, subject = gene_name(query), gene_name(subject)
        if query not in qorder or subject not in sorder: continue
        qi, q = qorder[query]
        si, s = sorder[subject]
        
        # remove redundancy a<->b when doing self-self BLAST
        if is_self and qi > si:
            query, subject = subject, query
            qi, si = si, qi
            q, s = s, q

        key = query, subject
        if key in seen: continue
        seen.add(key)
        b.query, b.subject = key

        b.qi, b.si = qi, si
        b.qseqid, b.sseqid = q['seqid'], s['seqid']
        
        filtered_blasts.append(b)

    if tandem_Nmax:
        print >>sys.stderr, "running the local dups filter (tandem_Nmax=%d)..." % tandem_Nmax

        qtandems = tandem_grouper(qbed, filtered_blasts,
                flip=True, tandem_Nmax=tandem_Nmax)
        standems = tandem_grouper(sbed, filtered_blasts, 
                flip=False, tandem_Nmax=tandem_Nmax)

        qdups_fh = open(op.splitext(qbed_file)[0] + ".localdups", "w")

        if is_self:
            for s in standems: qtandems.join(*s)
            qdups_to_mother = write_localdups(qdups_fh, qtandems, qbed)
            sdups_to_mother = qdups_to_mother
        else:
            qdups_to_mother = write_localdups(qdups_fh, qtandems, qbed)
            sdups_fh = open(op.splitext(sbed_file)[0] + ".localdups", "w")
            sdups_to_mother = write_localdups(sdups_fh, standems, sbed)

        # write out new .bed after tandem removal
        write_new_bed(qbed, qdups_to_mother, is_flat_fmt=is_flat_fmt)
        if not is_self:
            write_new_bed(sbed, sdups_to_mother, is_flat_fmt=is_flat_fmt)
        
        before_filter = len(filtered_blasts)
        filtered_blasts = list(filter_tandem(filtered_blasts, \
                qdups_to_mother, sdups_to_mother))
        print >>sys.stderr, "after filter (%d->%d)..." % \
                (before_filter, len(filtered_blasts))

        qnew_name = "%s.nolocaldups%s" % op.splitext(qbed.filename)
        snew_name = "%s.nolocaldups%s" % op.splitext(sbed.filename)

        if is_flat_fmt:
            from flatfeature import Flat
            qbed_new = Flat(qnew_name)
            sbed_new = Flat(snew_name)
        else:
            qbed_new = Bed(qnew_name)
            sbed_new = Bed(snew_name)

        qorder = get_order(qbed_new) 
        sorder = get_order(sbed_new) 

    if top_N:
        before_filter = len(filtered_blasts)
        print >>sys.stderr, "running the top-N filter (top_N=%d)..." % top_N
        filtered_blasts = list(filter_top_n(filtered_blasts, top_N=top_N))
        print >>sys.stderr, "after filter (%d->%d)..." % (before_filter, len(filtered_blasts))

    if cscore:
        before_filter = len(filtered_blasts)
        print >>sys.stderr, "running the cscore filter (cscore>=%.2f)..." % cscore
        filtered_blasts = list(filter_cscore(filtered_blasts, cscore=cscore))
        print >>sys.stderr, "after filter (%d->%d)..." % (before_filter, len(filtered_blasts))

    # this is the final output we will write to after BLAST filters
    raw_name = "%s.raw" % op.splitext(blast_file)[0]
    raw_fh = open(raw_name, "w")

    write_raw(qorder, sorder, filtered_blasts, raw_fh)
    #write_new_blast(filtered_blasts) 


def write_localdups(dups_fh, tandems, bed):

    print >>sys.stderr, "write local dups to file", dups_fh.name

    tandem_groups = []
    for group in tandems:
        rows = [bed[i] for i in group]
        # within the tandem groups, genes are sorted with decreasing size
        rows.sort(key=lambda a: abs(a['end'] - a['start']), reverse=True)
        tandem_groups.append([row['accn'] for row in rows])

    dups_to_mother = {}
    for accns in sorted(tandem_groups):
        print >>dups_fh, "\t".join(accns)
        for dup in accns[1:]:
            dups_to_mother[dup] = accns[0]

    return dups_to_mother


def write_new_bed(bed, children, is_flat_fmt=False):
    # generate local dup removed annotation files
    out_name = "%s.nolocaldups%s" % op.splitext(bed.filename)
    print >>sys.stderr, "write tandem-filtered bed file %s" % out_name
    fh = open(out_name, "w")
    for i, row in enumerate(bed):
        if row['accn'] in children: continue
        if is_flat_fmt:
            if i == 0: print >>fh, "\t".join(Flat.names)
            print >>fh, Flat.row_string(row)
        else:
            print >>fh, row
    fh.close()


def write_raw(qorder, sorder, filtered_blasts, raw_fh):

    print >>sys.stderr, "write raw file %s" % raw_fh.name
    for b in filtered_blasts: 
        qi, q = qorder[b.query]
        si, s = sorder[b.subject]
        qseqid, sseqid = q['seqid'], s['seqid']

        score = 50 if b.evalue == 0 else min(int(-log10(b.evalue)), 50)
        print >>raw_fh, "\t".join(map(str, (qseqid, qi, sseqid, si, score)))


def write_new_blast(filtered_blasts, fh=sys.stdout):
    for b in filtered_blasts:
        print >>fh, b

# ---------------- All BLAST filters ----------------

def filter_cscore(blast_list, cscore=.5):

    best_score = {}
    for b in blast_list:
        if b.query not in best_score or b.score > best_score[b.query]:
            best_score[b.query] = b.score
        if b.subject not in best_score or b.score > best_score[b.subject]:
            best_score[b.subject] = b.score

    for b in blast_list:
        cur_cscore = b.score / max(best_score[b.query], best_score[b.subject])
        if cur_cscore > cscore:
            yield b


def filter_top_n(blast_list, top_N=10):

    q_hits = collections.defaultdict(int)
    s_hits = collections.defaultdict(int)
    top_N *= 2

    for b in blast_list:
        q_hits[b.query] += 1
        s_hits[b.subject] += 1
        if q_hits[b.query] + s_hits[b.subject] > top_N: continue
        yield b
    

def filter_tandem(blast_list, qdups_to_mother, sdups_to_mother):
    
    mother_blast = []
    for b in blast_list:
        if b.query in qdups_to_mother: b.query = qdups_to_mother[b.query]
        if b.subject in sdups_to_mother: b.subject = sdups_to_mother[b.subject]
        mother_blast.append(b)
    
    mother_blast.sort(key=lambda b: b.score, reverse=True)
    seen = set() 
    for b in mother_blast:
        if b.query==b.subject: continue
        key = b.query, b.subject
        if key in seen: continue
        seen.add(key)
        yield b


def tandem_grouper(bed, blast_list, tandem_Nmax=10, flip=True):
    if not flip:
        simple_blast = [(b.query, (b.sseqid, b.si)) for b in blast_list if b.evalue < 1e-10] 
    else:
        simple_blast = [(b.subject, (b.qseqid, b.qi)) for b in blast_list if b.evalue < 1e-10] 

    simple_blast.sort()

    standems = Grouper()
    for name, hits in itertools.groupby(simple_blast, key=lambda x:x[0]):
        # these are already sorted.
        hits = [x[1] for x in hits]
        for ia, a in enumerate(hits[:-1]):
            b = hits[ia + 1]
            # on the same chromosome and rank difference no larger than tandem_Nmax
            if b[1] - a[1] <= tandem_Nmax and b[0] == a[0]: 
                standems.join(a[1], b[1])

    return standems


if __name__ == "__main__":
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed", 
            help="path to qbed or qflat")
    parser.add_option("--sbed", dest="sbed", 
            help="path to sbed or sflat")

    filter_group = optparse.OptionGroup(parser, "BLAST filters")
    filter_group.add_option("--tandem_Nmax", dest="tandem_Nmax", type="int", default=None, 
            help="merge tandem genes within distance [default: %default]")
    filter_group.add_option("--top_N", dest="top_N", type="int", default=None,
            help="retain only top-N hit from BLAST [default: %default]")
    filter_group.add_option("--cscore", type="float", default=None,
            help="retain hits that have good bitscore [default: %default]")
    
    parser.add_option_group(filter_group)

    (options, blast_files) = parser.parse_args()

    if not (len(blast_files) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    main(blast_files[0], options)

