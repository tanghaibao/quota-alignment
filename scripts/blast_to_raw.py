"""
accepts .bed format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

or .flat format: http://github.com/brentp/flatfeature/

and a blast file.
if the input is query.qbed and subject.qbed,
files query.localdups and subject.localdups are created containing the parent|offspring dups
dups inferred by subjects hitting the same query or queries hitting the same subject

and new blast is printed to stdout, and a .raw file is created
"""

import sys
import os.path as op
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from grouper import Grouper
import collections
from math import log10

class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si')

    def __init__(self, sline):
        args = sline.split("\t")
        self.query  =args[0]
        self.subject  = args[1]
        self.pctid =float(args[2])
        self.hitlen =int(args[3])
        self.nmismatch =int(args[4])
        self.ngaps =int(args[5])
        self.qstart =int(args[6])
        self.qstop =int(args[7])
        self.sstart =int(args[8])
        self.sstop =int(args[9])
        self.evalue =float(args[10])
        self.score =float(args[11])

    def __repr__(self):
        return "BLine('%s' to '%s', eval=%.3f, score=%.1f)" % (self.query, self.subject, self.evalue, self.score)


    def to_blast_line(self):
        return "\t".join(map(str, (getattr(self, attr) for attr in BlastLine.__slots__)))

class BedRow(object):
    __slots__ = ('seqid', 'start', 'end', 'name', 'stuff')
    def __init__(self, seqid, start, end, name, stuff):
       self.seqid = seqid
       self.start = int(start)
       self.end = int(end)
       self.name = name
       self.stuff = stuff

    def to_bed_line(self):
        return "\t".join(map(str, [self.seqid, self.start, self.end, self.name] + self.stuff))

class Bed(list):
    
    def __init__(self, filename):
        self.filename = filename
        beds = []
        for line in open(filename):
            if line[0] == "#": continue
            if line.startwith('track'): continue
            line = line.split("\t")
            chr, start, end, name = line[:4]
            br = BedRow(chr, start, end, name, line[4:])
            beds.append(br)
        self.seqids = sorted(set(b.seqid for b in beds))
        self.beds = sorted(beds, key=lambda a: (a.chr, a.start)) 

    def __getitem__(self, i):
        return self.beds[i]
    def __len__(self):
        return len(self.beds)


def main(qbed_file, sbed_file, blast_file):
    if qbed_file.endswith(".flat"):
        from flatfeature import Flat
        qbed = Flat(qbed_file)
        sbed = Flat(sbed_file)
        qorder = dict((f['accn'], (i, f)) for (i, f) in enumerate(qbed))
        sorder = dict((f['accn'], (i, f)) for (i, f) in enumerate(sbed))
    else:
        qbed = Bed(qbed_file)
        sbed = Bed(sbed_file)
        qorder = dict((f.name, (i, f)) for (i, f) in enumerate(qbed))
        sorder = dict((f.name, (i, f)) for (i, f) in enumerate(sbed))


    seen = {}
    blasts = sorted([BlastLine(line) for line in open(blast_file)], 
                    key=lambda b: b.evalue)
    filtered_blasts = []
    for b in blasts:
        qi, q = qorder[b.query]
        si, s = sorder[b.subject]
        if (qi, si) in seen: continue
        seen[qi, si] = 1
        b.qi = qi
        b.si = si
        if qbed_file.endswith(".flat"):
            b.qseqid, b.sseqid = q['seqid'], s['seqid']
        else:
            b.qseqid, b.sseqid = q.seqid, s.seqid
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
    for accns in sorted(tandem_grouper(qbed, sbed, filtered_blasts, flip=False)):
        print >>sdups_fh, "|".join(accns)
        for dup in accns[1:]:
            smothers[dup] = accns[0]
        schildren.update(accns[1:])


    for bed, children in [(qbed, qchildren), (sbed, schildren)]:
        out_name = "%s.filtered%s" % op.splitext(bed.filename)
        fh = open(out_name, "w")
        for i, row in enumerate(bed):
            if bed.filename.endswith(".flat"):
                if i == 0: print >>fh, "\t".join(Flat.names)
                if row['accn'] in children: continue
                print >>fh, Flat.row_string(row)
            else:
                if row.name in children: continue
                print >>fh, row.to_bed_line()
        fh.close()
    write_new_files(qbed, sbed, filtered_blasts, qmothers, smothers, blast_file)


def write_new_files(qbed, sbed, filtered_blasts, qmothers, smothers, blast_file):
    qnew_name = "%s.filtered%s" % op.splitext(qbed.filename)
    snew_name = "%s.filtered%s" % op.splitext(sbed.filename)

    raw_name = "%s.filtered.raw" % op.splitext(blast_file)[0]
    raw_fh = open(raw_name, "w")
    if qbed.filename.endswith(".flat"):
        from flatfeature import Flat
        qbed_new = Flat(qnew_name)
        sbed_new = Flat(snew_name)
        qorder = dict((f['accn'], (i, f)) for (i, f) in enumerate(qbed_new))
        sorder = dict((f['accn'], (i, f)) for (i, f) in enumerate(sbed_new))
    else:
        qbed_new = Bed(qnew_name)
        sbed_new = Bed(snew_name)
        qorder = dict((f.name, (i, f)) for (i, f) in enumerate(qbed_new))
        sorder = dict((f.name, (i, f)) for (i, f) in enumerate(sbed_new))

    for b in filter_to_mothers_and_write_blast(filtered_blasts, qmothers, smothers):
        # write raw file.
        qi, q = qorder[b.query]
        si, s = sorder[b.subject]
        if qbed.filename.endswith(".flat"):
            qseqid, sseqid = q['seqid'], s['seqid']
        else:
            qseqid, sseqid = q.seqid, s.seqid
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
        print b.to_blast_line()
        yield b


#    for accns in sorted(tandem_grouper(sbed, qbed, filtered_blasts, flip=True)):
def tandem_grouper(abed, bbed, blast_list, Tandem_Nmax=10, flip=False):
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
                    if b - a > Tandem_Nmax: break
                    standems.join(a, b)

        for group in standems:
            rows = [bbed[i] for i in group]
            rows.sort(key=lambda a: (a['end'] - a['start']), reverse=True)
            yield [row['accn'] for row in rows]


if __name__ == "__main__":
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("-q", dest="qbed", help="path to qbed or qflat")
    parser.add_option("-s", dest="sbed", help="path to sbed or sflat")

    (options, blast_files) = parser.parse_args()

    if not (len(blast_files) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    main(options.qbed, options.sbed, blast_files[0])
