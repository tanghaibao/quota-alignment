#!/usr/bin/env python
# -*- coding: UTF-8 -*-


"""
%prog gff_file

convert the gff format into the bed format 
"""

import os.path as op
import sys

try:
    from BCBio.GFF import GFFParser
except:
    print >>sys.stderr, "BCBio.GFF module not found, try download and install from`" \
            "<http://github.com/chapmanb/bcbb/tree/master/gff/BCBio/>"
    sys.exit(1)


def gff_to_bed(gff_file, bed_fh=sys.stdout, cds=True):

    parser = GFFParser()
    seqids = parser.parse(gff_file, None)

    for seqid in seqids:
        for feat in seqid.features:
            subf = feat.sub_features
            if feat.type in ("chromosome", "protein"): continue
            is_cds = any(f.type=="mRNA" or f.type=="CDS" for f in subf) and\
                    feat.type=="gene"
            if cds == is_cds:
                print >>bed_fh, "\t".join(str(x) for x in (seqid.id, feat.location.start, \
                        feat.location.end, feat.id, feat.type))


if __name__ == "__main__":

    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--noncoding", dest="cds", action="store_false", 
            default=True, help="extract coding features?")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        sys.exit(parser.print_help())

    gff_file = args[0]

    gff_to_bed(gff_file, cds=options.cds)

