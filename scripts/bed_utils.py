"""
Classes to handle the .bed file and .raw file
"""
# get the gene order given a Bed object
get_order = lambda bed: dict((f['accn'], (i, f)) for (i, f) in enumerate(bed))

class BedLine(object):
    # the Bed format supports more columns. we only need
    # the first 4, but keep the information in 'stuff'.
    __slots__ = ("seqid", "start", "end", "accn", "stuff")

    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.seqid = args[0]
        self.start = int(args[1])
        self.end = int(args[2])
        self.accn = args[3]
        self.stuff = args[4:] if len(args) > 4 else None

    def __str__(self):
        s = "\t".join(map(str, [getattr(self, attr) \
                    for attr in BedLine.__slots__[:-1]]))
        if self.stuff:
            s += "\t" + "\t".join(self.stuff)
        return s

    def __getitem__(self, key):
        return getattr(self, key)


class Bed(list):

    def __init__(self, filename):
        self.filename = filename
        beds = []
        for line in open(filename):
            if line[0] == "#": continue
            if line.startswith('track'): continue
            beds.append(BedLine(line))

        self.seqids = sorted(set(b.seqid for b in beds))
        self.beds = sorted(beds, key=lambda a: (a.seqid, a.start, a.accn))

    def __getitem__(self, i):
        return self.beds[i]

    def __len__(self):
        return len(self.beds)

    def __iter__(self):
        for b in self.beds:
            yield b

    def get_order(self):
        return dict((f.accn, (i, f)) for (i, f) in enumerate(self))

class RawLine(object):
    __slots__ = ("seqid_a", "pos_a", "seqid_b", "pos_b", "score")

    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.seqid_a = args[0]
        self.pos_a = int(args[1])
        self.seqid_b = args[2]
        self.pos_b = int(args[3])
        self.score = int(args[4])

    def __str__(self):
        return "\t".join(map(str, [getattr(self, attr) \
                for attr in RawLine.__slots__]))

    def __getitem__(self, key):
        return getattr(self, key)


class Raw(list):

    def __init__(self, filename):
        self.filename = filename
        raws = []
        for line in open(filename):
            if line[0] == "#": continue
            raws.append(RawLine(line))
        self.raws = raws

    def __len__(self):
        return len(self.raws)

    def __iter__(self):
        for b in self.raws:
            yield b


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
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        return "\t".join(map(str, [getattr(self, attr) \
                for attr in BlastLine.__slots__][:-4]))

