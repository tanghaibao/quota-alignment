#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog blast_file -q query.bed -s subject.bed

accepts .bed format: <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>
or .flat format: <http://github.com/brentp/flatfeature/>
and a blast file.

if the input is query.bed and subject.bed, the script files query.localdups and subject.localdups are created containing the parent|offspring dups, as inferred by subjects hitting the same query or queries hitting the same subject

and new blast is printed to stdout, and a .raw file (which is the input for the quota-align pipeline <http://github.com/tanghaibao/quota-alignment/>) is created
"""

import sys
import os.path as op
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from grouper import Grouper
import collections
from math import log10

# Helper functions in the BLAST filtering to get rid alternative splicings
def gene_name(st):
    # this is ugly, but different groups are inconsistent
    # with how the alternative splicings are named;
    # mostly it can be done by removing the suffix
    # except for papaya (evm...) and maize (somewhat complicated)
    if st.startswith("ev"):
        return st
    if st.startswith("Os"):
        return st.rsplit("-",1)[0]
    return st.rsplit(".", 1)[0]


class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si')

    def __init__(self, sline):
        args = sline.split("\t")
        self.query = args[0]
        self.subject = args[1]
        self.pctid = float(args[2])
        self.hitlen = int(args[3])
        self.nmismatch = int(args[4])
        self.ngaps = int(args[5])
        self.qstart = int(args[6])
        self.qstop = int(args[7])
        self.sstart = int(args[8])
        self.sstop = int(args[9])
        self.evalue = float(args[10])
        self.score = float(args[11])

    def __repr__(self):
        return "BLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        return "\t".join(map(str, (getattr(self, attr) \
                for attr in BlastLine.__slots__)))


class BedRow(object):
    __slots__ = ('seqid', 'start', 'end', 'accn', 'stuff')
    def __init__(self, seqid, start, end, accn, stuff):
       self.seqid = seqid 
       self.start = int(start)
       self.end = int(end)
       self.accn = accn 
       self.stuff = stuff

    def __str__(self):
        to_write = [self.seqid, self.start, self.end, self.accn]
        if self.stuff:
            to_write.extend(self.stuff)
        return "\t".join(map(str, to_write)) 

    def __getitem__(self, key): 
        return getattr(self, key)


class Bed(list):
    
    def __init__(self, filename):
        self.filename = filename
        beds = []
        for line in open(filename):
            if line[0] == "#": continue
            if line.startswith('track'): continue
            line = line.strip().split("\t")
            chr, start, end, name = line[:4]
            start, end = int(start), int(end)
            br = BedRow(chr, start, end, name, line[4:])
            beds.append(br)
        self.seqids = sorted(set(b.seqid for b in beds))
        self.beds = sorted(beds, key=lambda a: (a.seqid, a.start)) 

    def __getitem__(self, i):
        return self.beds[i]

    def __len__(self):
        return len(self.beds)

    def __iter__(self):
        for b in self.beds:
            yield b


def main(qbed_file, sbed_file, blast_file, Nmax=10, is_flat_fmt=True):

    print >>sys.stderr, "read annotation files %s and %s" % (qbed_file, sbed_file)
    if is_flat_fmt:
        from flatfeature import Flat
        qbed = Flat(qbed_file)
        sbed = Flat(sbed_file)
    else:
        qbed = Bed(qbed_file)
        sbed = Bed(sbed_file)

    qorder = dict((f['accn'], (i, f)) for (i, f) in enumerate(qbed))
    sorder = dict((f['accn'], (i, f)) for (i, f) in enumerate(sbed))

    print >>sys.stderr, "read BLAST file %s" % blast_file
    seen = set() 
    blasts = sorted([BlastLine(line) for line in open(blast_file)], 
                    key=lambda b: b.score, reverse=True)

    filtered_blasts = []
    for b in blasts:
        query, subject = gene_name(b.query), gene_name(b.subject)
        if query not in qorder or subject not in sorder: continue
        qi, q = qorder[query]
        si, s = sorder[subject]
        if (qi, si) in seen: continue
        seen.add((qi, si))
        b.qi = qi
        b.si = si
        b.qseqid, b.sseqid = q['seqid'], s['seqid']
        filtered_blasts.append(b)

    qmothers = {}
    qchildren = set()
    qdups_fh = open(op.splitext(qbed_file)[0] + ".localdups", "w")
    for accns in sorted(tandem_grouper(sbed, qbed, filtered_blasts, flip=True)):
        print >>qdups_fh, "|".join(accns)
        for dup in accns[1:]:
            qmothers[dup] = accns[0]
        qchildren.update(accns[1:])

    smothers = {}
    schildren = set()
    sdups_fh = open(op.splitext(sbed_file)[0] + ".localdups", "w")
    for accns in sorted(tandem_grouper(qbed, sbed, filtered_blasts, flip=False, Nmax=Nmax)):
        print >>sdups_fh, "|".join(accns)
        for dup in accns[1:]:
            smothers[dup] = accns[0]
        schildren.update(accns[1:])

    print >>sys.stderr, "write local dups to files %s and %s" % \
            (qdups_fh.name, sdups_fh.name)

    # generate filtered annotation files
    for bed, children in [(qbed, qchildren), (sbed, schildren)]:
        out_name = "%s.filtered%s" % op.splitext(bed.filename)
        fh = open(out_name, "w")
        for i, row in enumerate(bed):
            if row['accn'] in children: continue
            if is_flat_fmt:
                if i == 0: print >>fh, "\t".join(Flat.names)
                print >>fh, Flat.row_string(row)
            else:
                print >>fh, row
        fh.close()

    write_new_files(qbed, sbed, filtered_blasts, qmothers, smothers, 
            blast_file, is_flat_fmt)


def write_new_files(qbed, sbed, filtered_blasts, qmothers, smothers, 
        blast_file, is_flat_fmt):
    qnew_name = op.basename("%s.filtered%s" % op.splitext(qbed.filename))
    snew_name = op.basename("%s.filtered%s" % op.splitext(sbed.filename))

    raw_name = op.basename("%s.filtered.raw" % op.splitext(blast_file)[0])
    raw_fh = open(raw_name, "w")

    if is_flat_fmt:
        from flatfeature import Flat
        qbed_new = Flat(qnew_name)
        sbed_new = Flat(snew_name)
    else:
        qbed_new = Bed(qnew_name)
        sbed_new = Bed(snew_name)

    qorder = dict((f['accn'], (i, f)) for (i, f) in enumerate(qbed_new))
    sorder = dict((f['accn'], (i, f)) for (i, f) in enumerate(sbed_new))

    print >>sys.stderr, "write raw file %s" % raw_fh.name
    for b in filter_to_mothers_and_write_blast(filtered_blasts, qmothers, smothers):
        query, subject = gene_name(b.query), gene_name(b.subject)
        if query not in qorder or subject not in sorder: continue
        qi, q = qorder[query]
        si, s = sorder[subject]
        qseqid, sseqid = q['seqid'], s['seqid']

        score = 50 if b.evalue == 0 else min(int(-log10(b.evalue)), 50)
        print >>raw_fh, "\t".join(map(str, (qseqid, qi, sseqid, si, score)))


def filter_to_mothers_and_write_blast(blast_list, qmothers, smothers):
    
    mother_blast = []
    for b in blast_list:
        if b.query in qmothers: b.query = qmothers[b.query]
        if b.subject in smothers: b.subject = smothers[b.subject]
        mother_blast.append(b)
    
    mother_blast.sort(key=lambda b: b.evalue)
    seen = {}
    for b in mother_blast:
        key = b.query, b.subject
        if key in seen: continue
        seen[key] = True
        print b
        yield b


# for accns in sorted(tandem_grouper(sbed, qbed, filtered_blasts, flip=True)):
def tandem_grouper(abed, bbed, blast_list, Nmax=10, flip=False):
    for seqid in abed.seqids:
        ai_to_bi = collections.defaultdict(list)
        for b in blast_list:
            if flip:
                if b.qseqid == seqid:
                    ai_to_bi[b.si].append(b.qi)
            else:
                if b.sseqid == seqid:
                    ai_to_bi[b.qi].append(b.si)

        standems = Grouper()
        for qi in ai_to_bi:
            ordered = sorted(ai_to_bi[qi])
            for ia, a in enumerate(ordered[:-1]):
                for b in ordered[ia + 1:]:
                    assert b > a, (b, a)
                    if b - a > Nmax: break
                    standems.join(a, b)

        for group in standems:
            rows = [bbed[i] for i in group]
            rows.sort(key=lambda a: (a['end'] - a['start']), reverse=True)
            yield [row['accn'] for row in rows]


if __name__ == "__main__":
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed", 
            help="path to qbed or qflat")
    parser.add_option("--sbed", dest="sbed", 
            help="path to sbed or sflat")
    parser.add_option("--Nmax", dest="Nmax", type="int", default=10, 
            help="merge tandem genes within distance "\
                 "[default: %default genes]")

    (options, blast_files) = parser.parse_args()

    if not (len(blast_files) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    is_flat_fmt = options.qbed.endswith(".flat")

    main(options.qbed, options.sbed, blast_files[0], Nmax=options.Nmax,
            is_flat_fmt=is_flat_fmt)

