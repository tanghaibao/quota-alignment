class BedLine(object):
    __slots__ = ("seqid", "start", "end", "accn")

    def __init__(self, sline): 
        args = sline.strip().split("\t")
        self.seqid = args[0] 
        self.start = int(args[1])
        self.end = int(args[2])
        self.accn = args[3] 

    def __str__(self):
        return "\t".join(map(str, [getattr(self, attr) \
                for attr in BedLine.__slots__]))

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
        self.beds = sorted(beds, key=lambda a: (a.seqid, a.start)) 

    def __getitem__(self, i):
        return self.beds[i]

    def __len__(self):
        return len(self.beds)

    def __iter__(self):
        for b in self.beds:
            yield b
	
    def get_order(self):
        return dict((f['accn'], (i, f)) for (i, f) in enumerate(self.beds))
		

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
                for attr in BedLine.__slots__]))

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